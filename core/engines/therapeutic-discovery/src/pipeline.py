"""
Drug Discovery Pipeline Orchestrator.

Implements the 10-stage pipeline:
0. Initialize - Load config and setup
1. Normalize Target - Validate and enrich target hypothesis
2. Structure Discovery - Find/fetch protein structures
3. Structure Prep - Prepare structures for docking
4. Molecule Generation - Generate candidate molecules
5. Chemistry QC - Filter molecules by drug-likeness
6. Conformers - Generate 3D conformers
7. Docking - Dock molecules to protein
8. Ranking - Score and rank candidates
9. Reporting - Generate final report

Based on phase-5-6.pdf specification.
"""
import os
import json
import time
import uuid
from pathlib import Path
from typing import List, Dict, Any, Optional, Callable
from datetime import datetime
from loguru import logger

from .models import (
    TargetHypothesis,
    StructureInfo,
    StructureManifest,
    GeneratedMolecule,
    MoleculeProperties,
    DockingResult,
    RankedCandidate,
    PipelineRun,
    PipelineConfig,
    PediatricSafetyFlag,
    PediatricSafetyAssessment,
)
from .nim_clients import create_nim_clients, NIMServiceManager
from .checkpoint import CheckpointManager


def pediatric_safety_assessment(molecule_data: dict) -> PediatricSafetyAssessment:
    """Assess drug candidate safety for pediatric patients.

    Pediatric patients have unique pharmacokinetic considerations:
    - Immature hepatic metabolism (CYP3A4 reaches adult levels by age 1-2)
    - Developing blood-brain barrier (more permeable in young children)
    - Higher body water percentage (affects hydrophilic drug distribution)
    - Cardiac sensitivity (lower threshold for QT prolongation)
    - Renal maturation (GFR reaches adult levels by age 2)

    This function returns a clinical-shaped verdict, so it refuses to produce one it cannot
    support. The only `critical` rule -- and therefore the only thing that can set
    `pediatric_safe=False` -- is hERG cardiac liability. When no hERG prediction is supplied,
    the verdict is `None` ("not assessed"), never `True`. Before 2026-11-04 the caller fed this
    function `herg_ic50` from a hand-written expression, `max(0.5, 30 - logP*4 - (MW-350)*0.02)`,
    cited as "Ref: Aronov (2005) J Med Chem 48:1341-1348" -- a paper that publishes no such
    equation -- and only when `logP > 3.5 and MW > 350`. Everything else got
    `herg_ic50 = inf` by default, raised no cardiac flag, and was reported
    `pediatric_safe: true, recommendation: "Suitable for pediatric development"`. A candidate
    with dofetilide's physicochemistry (MW 441, logP 2.3, hERG IC50 ~10 nM) passed that way.
    """
    safety_flags = []
    unassessed = []

    # Molecular weight check for BBB crossing
    mw = molecule_data.get("molecular_weight", 0)
    if mw > 500:
        safety_flags.append(PediatricSafetyFlag(
            flag="HIGH_MW_BBB_RISK",
            detail=f"MW {mw} may limit CNS penetration for pediatric brain tumors",
            severity="medium",
        ))

    # Lipophilicity check
    logp = molecule_data.get("logP", 0)
    if logp > 5:
        safety_flags.append(PediatricSafetyFlag(
            flag="HIGH_LIPOPHILICITY",
            detail=f"LogP {logp} may cause hepatotoxicity in immature liver",
            severity="high",
        ))

    # hERG cardiac liability. Two accepted inputs, both of which must come from a model:
    #   * `herg_inhibition_probability` — ADMETPredictor's hERG endpoint, P(inhibitor) in [0,1]
    #   * `herg_ic50` — a measured or predicted IC50 in uM, from a named source
    # Absence is not safety: it makes the verdict unassessable.
    herg_prob = molecule_data.get("herg_inhibition_probability")
    herg_ic50 = molecule_data.get("herg_ic50")
    cardiac_assessed = False

    if isinstance(herg_prob, (int, float)):
        cardiac_assessed = True
        if herg_prob >= 0.5:
            safety_flags.append(PediatricSafetyFlag(
                flag="HERG_CARDIAC_RISK",
                detail=(f"predicted hERG inhibition probability {herg_prob:.2f} — pediatric "
                        "patients are more sensitive to QT prolongation"),
                severity="critical",
            ))
    elif isinstance(herg_ic50, (int, float)):
        cardiac_assessed = True
        if herg_ic50 < 10:
            safety_flags.append(PediatricSafetyFlag(
                flag="HERG_CARDIAC_RISK",
                detail="hERG IC50 < 10uM — pediatric patients more sensitive to QT prolongation",
                severity="critical",
            ))
    else:
        unassessed.append(
            "no cardiac liability model was run: neither a hERG inhibition probability nor a "
            "hERG IC50 was supplied, and none is estimated from physicochemical properties"
        )

    # Teratogenicity check (for adolescent patients)
    if molecule_data.get("teratogenic_risk", False):
        safety_flags.append(PediatricSafetyFlag(
            flag="TERATOGENICITY",
            detail="Teratogenic risk — requires pregnancy prevention program for adolescent patients",
            severity="high",
        ))

    # High TPSA may limit oral bioavailability in pediatric formulations
    tpsa = molecule_data.get("tpsa", 0)
    if tpsa > 140:
        safety_flags.append(PediatricSafetyFlag(
            flag="HIGH_TPSA_ORAL_RISK",
            detail=f"TPSA {tpsa} may limit oral bioavailability — pediatric formulation challenges",
            severity="medium",
        ))

    # Many rotatable bonds — potential for variable pediatric absorption
    rotatable = molecule_data.get("rotatable_bonds", 0)
    if rotatable > 10:
        safety_flags.append(PediatricSafetyFlag(
            flag="HIGH_FLEXIBILITY",
            detail=f"{rotatable} rotatable bonds — variable GI absorption in pediatric patients",
            severity="low",
        ))

    has_critical = len([f for f in safety_flags if f.severity == "critical"]) > 0

    if has_critical:
        verdict = False
        recommendation = "Critical pediatric safety flag — not suitable without further work"
    elif not cardiac_assessed:
        # No critical rule could be evaluated, so there is no verdict to give.
        verdict = None
        recommendation = (
            "NOT ASSESSED — " + "; ".join(unassessed) +
            ". Run ADMETPredictor().predict(smiles) and pass its hERG endpoint as "
            "herg_inhibition_probability to obtain a verdict."
        )
    elif safety_flags:
        verdict = True
        recommendation = "Requires additional pediatric safety evaluation"
    else:
        verdict = True
        recommendation = "No pediatric safety flag raised by the rules applied"

    return PediatricSafetyAssessment(
        pediatric_safe=verdict,
        safety_flags=safety_flags,
        recommendation=recommendation,
        cardiac_liability_assessed=cardiac_assessed,
        unassessed_reasons=unassessed,
    )


class DrugDiscoveryPipeline:
    """
    Main pipeline orchestrator for drug discovery.
    """

    STAGES = [
        "Initialize",
        "Normalize Target",
        "Structure Discovery",
        "Structure Prep",
        "Molecule Generation",
        "Chemistry QC",
        "Conformers",
        "Docking",
        "Ranking",
        "Reporting",
    ]

    def __init__(
        self,
        config: PipelineConfig,
        nim_manager: NIMServiceManager = None,
        progress_callback: Callable[[int, str], None] = None,
    ):
        self.config = config
        self.nim = nim_manager or create_nim_clients()
        self.progress_callback = progress_callback

        # Setup output directory
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Pipeline state
        self.run_id = str(uuid.uuid4())[:8]
        self.run = PipelineRun(
            run_id=self.run_id,
            target_gene=self.config.target_gene,
            started_at=datetime.now(),
            config=self.config.dict(),
        )

        # Lazy ADMET predictor for the pediatric cardiac-liability rule.
        # None = not loaded yet; False = loading failed, do not retry.
        self._admet = None

        # Checkpoint manager
        self.checkpoint_manager = CheckpointManager(self.output_dir)

        # Stage timing
        self.stage_timings: Dict[int, float] = {}

        # Intermediate results
        self.target: Optional[TargetHypothesis] = None
        self.structures: Optional[StructureManifest] = None
        self.molecules: List[GeneratedMolecule] = []
        self.docking_results: List[DockingResult] = []
        self.ranked_candidates: List[RankedCandidate] = []

    def _report_progress(self, stage: int, message: str):
        """Report progress to callback if available."""
        if self.progress_callback:
            self.progress_callback(stage, message)
        logger.info(f"[Stage {stage}] {message}")

    def _save_checkpoint(self, stage: int):
        """Save checkpoint after stage completes."""
        try:
            self.checkpoint_manager.save(self, stage)
        except Exception as e:
            logger.warning(f"Checkpoint save failed for stage {stage}: {e}")

    def run_pipeline(
        self,
        target: TargetHypothesis = None,
        resume_from_run_id: str = None,
    ) -> PipelineRun:
        """
        Run the complete pipeline.

        Args:
            target: Target hypothesis (uses config.target_gene if not provided)
            resume_from_run_id: If set, resume from the latest checkpoint of this run

        Returns:
            PipelineRun with results
        """
        resume_stage = -1
        if resume_from_run_id:
            latest = self.checkpoint_manager.find_latest_checkpoint(resume_from_run_id)
            if latest is not None:
                checkpoint = self.checkpoint_manager.load(resume_from_run_id, latest)
                if checkpoint:
                    self.checkpoint_manager.restore_pipeline_state(self, checkpoint)
                    self.run_id = resume_from_run_id
                    resume_stage = latest
                    self._report_progress(
                        latest,
                        f"Resuming from checkpoint (stage {latest}: {self.STAGES[latest]})",
                    )
            else:
                logger.warning(
                    f"No checkpoint found for run {resume_from_run_id}, starting fresh"
                )

        if resume_stage < 0:
            self.run = PipelineRun(
                run_id=self.run_id,
                target_gene=self.config.target_gene,
                started_at=datetime.now(),
                config=self.config.dict(),
            )

        try:
            stages = [
                (0, lambda: self.stage_0_initialize(target)),
                (1, self.stage_1_normalize_target),
                (2, self.stage_2_structure_discovery),
                (3, self.stage_3_structure_prep),
                (4, self.stage_4_molecule_generation),
                (5, self.stage_5_chemistry_qc),
                (6, self.stage_6_conformers),
                (7, self.stage_7_docking),
                (8, self.stage_8_ranking),
                (9, self.stage_9_reporting),
            ]

            for stage_num, stage_fn in stages:
                if resume_stage < stage_num:
                    t0 = time.time()
                    stage_fn()
                    self.stage_timings[stage_num] = round(time.time() - t0, 3)

            self.run.stage_timings = self.stage_timings
            self.run.status = "completed"
            self.run.completed_at = datetime.now()

        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            self.run.status = "failed"
            self.run.error_message = str(e)
            raise

        return self.run

    def stage_0_initialize(self, target: TargetHypothesis = None):
        """Stage 0: Initialize pipeline."""
        self._report_progress(0, "Initializing pipeline")

        # Create target if not provided
        if target is None:
            self.target = TargetHypothesis(
                gene=self.config.target_gene,
                rationale="Pipeline target",
            )
        else:
            self.target = target

        # Check NIM services
        available = self.nim.get_available_services()
        self._report_progress(0, f"Available NIM services: {available}")

        self.run.current_stage = 0
        self.run.stages_completed.append(0)
        self._save_checkpoint(0)

    def stage_1_normalize_target(self):
        """Stage 1: Normalize and validate target."""
        self._report_progress(1, f"Normalizing target: {self.target.gene}")

        # Validate gene symbol
        gene = self.target.gene.upper()
        self.target.gene = gene

        # Log if UniProt ID is missing (not required for pipeline to run)
        if not self.target.uniprot_id:
            logger.warning(f"No UniProt ID available for {gene}. Continuing without it.")

        self.run.current_stage = 1
        self.run.stages_completed.append(1)
        self._save_checkpoint(1)

    def stage_2_structure_discovery(self):
        """Stage 2: Discover protein structures."""
        self._report_progress(2, f"Discovering structures for {self.target.gene}")

        structures = []

        # Use provided PDB IDs
        for pdb_id in self.target.pdb_ids:
            structures.append(StructureInfo(
                pdb_id=pdb_id,
                method="Unknown",
            ))

        # Load from local structure cache or auto-fetch from RCSB
        cache_file = Path(self.config.structure_cache) / f"{self.target.gene.lower()}_structures.json"
        if cache_file.exists():
            with open(cache_file) as f:
                cached = json.load(f)
                for item in cached:
                    if not any(s.pdb_id == item.get("structure_id", "").replace("PDB:", "") for s in structures):
                        structures.append(StructureInfo(
                            pdb_id=item.get("structure_id", "").replace("PDB:", ""),
                            method=item.get("method", "Unknown"),
                            resolution=item.get("resolution"),
                            ligand_smiles=item.get("inhibitor_smiles"),
                        ))
        elif self.target.pdb_ids:
            # No cache — auto-fetch from RCSB and cache for future runs
            try:
                from .cryoem_evidence import CryoEMEvidenceManager
                em_manager = CryoEMEvidenceManager(
                    structures_dir=Path(self.config.structure_cache)
                )
                fetched = em_manager.auto_load_structures(
                    self.target.gene, self.target.pdb_ids
                )
                for entry in fetched:
                    pdb_id = entry.structure_id.replace("PDB:", "").strip()
                    if not any(s.pdb_id == pdb_id for s in structures):
                        structures.append(StructureInfo(
                            pdb_id=pdb_id,
                            method=entry.method,
                            resolution=entry.resolution,
                            ligand_smiles=entry.inhibitor_smiles,
                        ))
            except Exception as e:
                logger.warning(f"Auto-fetch from RCSB failed: {e}")

        self.structures = StructureManifest(
            target_gene=self.target.gene,
            structures=structures,
        )

        # Set primary structure
        best = self.structures.get_best_structure()
        if best:
            self.structures.primary_structure = best.pdb_id

        self._report_progress(2, f"Found {len(structures)} structures")
        self.run.current_stage = 2
        self.run.stages_completed.append(2)
        self._save_checkpoint(2)

    def _structure_cache_dirs(self):
        """Every directory that may hold a cached PDB file, in search order.

        `StructureViewer` -- the component that actually downloads structures -- caches them
        under `<engine>/data/structures/pdb_cache/`, which this pipeline never consulted, and
        the project's own PDB cache lives outside the repo entirely. `PDB_CACHE_DIR` makes that
        reachable without editing code.
        """
        dirs = []
        env = os.environ.get("PDB_CACHE_DIR")
        if env:
            dirs.append(Path(env))
        cache = Path(self.config.structure_cache)
        dirs += [cache / "pdb_cache", cache]
        dirs.append(Path(__file__).resolve().parent.parent / "data" / "structures" / "pdb_cache")
        return dirs

    def _resolve_structure_file(self, pdb_id: str):
        """Path to coordinates for `pdb_id`, or None. Never invents a file.

        Looks in every known cache directory, then (unless `HCLS_NO_PDB_FETCH` is set) fetches
        the entry from RCSB and caches it, because docking needs coordinates and an accession is
        not coordinates.
        """
        pdb_id = pdb_id.upper().replace("PDB:", "").strip()
        names = (f"{pdb_id}.pdb", f"{pdb_id.lower()}.pdb", f"pdb{pdb_id.lower()}.ent",
                 f"{pdb_id.lower()}.cif")
        for root in self._structure_cache_dirs():
            for name in names:
                candidate = root / name
                if candidate.is_file():
                    return candidate

        if os.environ.get("HCLS_NO_PDB_FETCH"):
            return None

        try:
            import requests
            response = requests.get(
                f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=30)
            response.raise_for_status()
            target_dir = Path(self.config.structure_cache) / "pdb_cache"
            target_dir.mkdir(parents=True, exist_ok=True)
            out = target_dir / f"{pdb_id.lower()}.pdb"
            out.write_text(response.text)
            logger.info(f"Fetched {pdb_id} from RCSB into {out}")
            return out
        except Exception as e:
            logger.warning(f"Could not resolve {pdb_id} to coordinates: {type(e).__name__}: {e}")
            return None

    def stage_3_structure_prep(self):
        """Stage 3: locate coordinates. No preparation step is implemented.

        The body used to be four lines of comment ("In production, would: 1. Download PDB files
        2. Remove water, ions 3. Add hydrogens 4. Identify binding site") followed by
        `struct.prepared = True` for every structure -- so `report_<run>.json` recorded
        `"prepared": true` for entries that had never been downloaded, let alone stripped of
        water or protonated, and a consumer read that as docking-ready input. (The same bug was
        fixed in hcls-orchestrator/run_pipeline.py; its honesty test only greps that file, so
        this copy survived.) This stage now does the one thing it can do honestly -- resolve
        each accession to coordinates on disk -- and reports `prepared: False`.
        """
        self._report_progress(3, "Locating structure coordinates (no preparation implemented)")

        located = 0
        for struct in self.structures.structures:
            if not struct.structure_file:
                path = self._resolve_structure_file(struct.pdb_id)
                if path:
                    struct.structure_file = str(path)
            if struct.structure_file:
                located += 1
            struct.prepared = False
            struct.preparation = "not performed — no preparation step is implemented"

        self._report_progress(
            3,
            f"{located}/{len(self.structures.structures)} structures have coordinates on disk; "
            "no water/ion removal, protonation or binding-site detection was performed",
        )
        self.run.current_stage = 3
        self.run.stages_completed.append(3)
        self._save_checkpoint(3)

    def stage_4_molecule_generation(self):
        """Stage 4: Generate candidate molecules."""
        self._report_progress(4, f"Generating {self.config.num_molecules} molecules")

        # Get seed SMILES
        seed = self.config.reference_compound_smiles

        # Try to get from structure with co-crystallized ligand
        if not seed and self.structures:
            for struct in self.structures.structures:
                if struct.ligand_smiles:
                    seed = struct.ligand_smiles
                    break

        if not seed:
            seed = "c1ccccc1"  # Benzene as last resort

        # Generate molecules using MolMIM
        generated = self.nim.molmim.generate(
            seed_smiles=seed,
            num_molecules=self.config.num_molecules,
        )

        # Convert to GeneratedMolecule objects
        for i, mol_data in enumerate(generated):
            self.molecules.append(GeneratedMolecule(
                id=f"{self.run_id}-mol-{i:04d}",
                smiles=mol_data.get("smiles", ""),
                source_seed=seed,
                generation_method=mol_data.get("method", "MolMIM"),
                target_gene=self.target.gene,
                generation_score=mol_data.get("score", 0.5),
            ))

        self.run.molecules_generated = len(self.molecules)
        self._report_progress(4, f"Generated {len(self.molecules)} molecules")
        self.run.current_stage = 4
        self.run.stages_completed.append(4)
        self._save_checkpoint(4)

    def stage_5_chemistry_qc(self):
        """Stage 5: Filter molecules by drug-likeness."""
        self._report_progress(5, "Running chemistry QC filters")

        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, QED
            HAS_RDKIT = True
        except ImportError:
            HAS_RDKIT = False
            logger.warning("RDKit not available, skipping detailed QC")

        passed = []
        for mol in self.molecules:
            if not HAS_RDKIT:
                mol.passed_qc = True
                passed.append(mol)
                continue

            rdmol = Chem.MolFromSmiles(mol.smiles)
            if rdmol is None:
                continue

            # Calculate properties
            mw = Descriptors.MolWt(rdmol)
            logp = Descriptors.MolLogP(rdmol)
            hbd = Lipinski.NumHDonors(rdmol)
            hba = Lipinski.NumHAcceptors(rdmol)
            tpsa = Descriptors.TPSA(rdmol)
            rotatable = Lipinski.NumRotatableBonds(rdmol)

            # Count Lipinski violations
            violations = 0
            if mw > 500:
                violations += 1
            if logp > 5:
                violations += 1
            if hbd > 5:
                violations += 1
            if hba > 10:
                violations += 1

            mol.properties = MoleculeProperties(
                molecular_weight=round(mw, 2),
                logP=round(logp, 2),
                hbd=hbd,
                hba=hba,
                tpsa=round(tpsa, 2),
                rotatable_bonds=rotatable,
                lipinski_violations=violations,
                qed=round(QED.qed(rdmol), 3),
            )

            # Apply filters
            if (mw <= self.config.max_mw and
                logp <= self.config.max_logp and
                violations <= self.config.max_lipinski_violations):
                mol.passed_qc = True
                passed.append(mol)

        self.molecules = passed
        self._report_progress(5, f"{len(passed)} molecules passed QC")
        self.run.current_stage = 5
        self.run.stages_completed.append(5)
        self._save_checkpoint(5)

    def stage_6_conformers(self):
        """Stage 6: Generate 3D conformers."""
        self._report_progress(6, "Generating 3D conformers")

        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            HAS_RDKIT = True
        except ImportError:
            HAS_RDKIT = False

        conformer_dir = self.output_dir / "conformers"
        conformer_dir.mkdir(exist_ok=True)

        for mol in self.molecules:
            if not HAS_RDKIT:
                continue

            rdmol = Chem.MolFromSmiles(mol.smiles)
            if rdmol is None:
                continue

            # Add hydrogens and generate 3D
            rdmol = Chem.AddHs(rdmol)
            result = AllChem.EmbedMolecule(rdmol, randomSeed=42)
            if result == 0:
                AllChem.MMFFOptimizeMolecule(rdmol)

                # Save conformer
                conf_file = conformer_dir / f"{mol.id}.sdf"
                writer = Chem.SDWriter(str(conf_file))
                writer.write(rdmol)
                writer.close()

                mol.conformer_file = str(conf_file)

        self._report_progress(6, f"Generated conformers for {len([m for m in self.molecules if m.conformer_file])} molecules")
        self.run.current_stage = 6
        self.run.stages_completed.append(6)
        self._save_checkpoint(6)

    def stage_7_docking(self):
        """Stage 7: Dock molecules to protein."""
        self._report_progress(7, "Running molecular docking")

        if not self.structures or not self.structures.structures:
            logger.warning("No structures available for docking")
            self.run.current_stage = 7
            self.run.stages_completed.append(7)
            return

        # Dock against COORDINATES. `primary_structure` is an accession ("5FTK"); passing it
        # straight to dock() meant the literal four characters were POSTed as the protein (and
        # uploaded to NVCF as chemical/x-pdb in cloud mode), so the stage docked against nothing
        # and stage 8 still weighted the result at 0.4. Resolve it, or refuse the stage.
        primary_id = self.structures.primary_structure or self.structures.structures[0].pdb_id
        receptor = next(
            (s for s in self.structures.structures
             if s.pdb_id == primary_id and s.structure_file),
            None,
        ) or next((s for s in self.structures.structures if s.structure_file), None)

        if receptor is None:
            path = self._resolve_structure_file(primary_id)
            if path is None:
                raise RuntimeError(
                    f"Docking requires structure coordinates and none could be resolved for "
                    f"{primary_id}. Searched "
                    f"{[str(d) for d in self._structure_cache_dirs()]} and the RCSB download "
                    "(set PDB_CACHE_DIR to point at an existing cache, or unset "
                    "HCLS_NO_PDB_FETCH to allow the fetch). Refusing to dock against a PDB "
                    "accession string."
                )
            receptor = next(s for s in self.structures.structures if s.pdb_id == primary_id)
            receptor.structure_file = str(path)

        primary = receptor.pdb_id
        receptor_coordinates = Path(receptor.structure_file).read_text()

        # Dock each molecule
        dock_failures = 0
        for mol in self.molecules:
            try:
                poses = self.nim.diffdock.dock(
                    protein_pdb=receptor_coordinates,
                    ligand_smiles=mol.smiles,
                    num_poses=self.config.num_poses,
                )
            except RuntimeError as e:
                dock_failures += 1
                if dock_failures >= 3:
                    raise RuntimeError(
                        f"DiffDock service appears down — {dock_failures} consecutive failures. "
                        f"Last error: {e}"
                    )
                logger.warning(f"Docking failed for {mol.id}: {e}")
                continue

            dock_failures = 0  # Reset on success
            if poses:
                best_pose = poses[0]
                self.docking_results.append(DockingResult(
                    molecule_id=mol.id,
                    structure_id=primary,
                    # Carry the method's own units through instead of forcing everything into
                    # `docking_score`: DiffDock reports a pose confidence, not an affinity.
                    score_type=best_pose.get("score_type", "binding_energy"),
                    docking_score=best_pose.get("docking_score"),
                    confidence_score=best_pose.get("confidence_score"),
                    confidence=best_pose.get("confidence", 0) or 0,
                    hydrogen_bonds=best_pose.get("hydrogen_bonds", 0),
                    contacts=best_pose.get("contacts", []),
                    method="DiffDock",
                ))

        self.run.molecules_docked = len(self.docking_results)
        self._report_progress(7, f"Docked {len(self.docking_results)} molecules")
        self.run.current_stage = 7
        self.run.stages_completed.append(7)
        self._save_checkpoint(7)

    def _predict_herg(self, smiles: str):
        """P(hERG inhibition) from the engine's trained ADMET suite, or None.

        The engine ships a real ADMET predictor with a trained `hERG` endpoint
        (src/admet_predictor.py); the pediatric assessment used a hand-written expression
        instead. Returning None when the model is unavailable is deliberate: the assessment then
        reports "not assessed" rather than defaulting to safe.
        """
        if not smiles:
            return None
        if self._admet is False:                     # a previous call already failed
            return None
        try:
            if self._admet is None:
                from .admet_predictor import ADMETPredictor
                self._admet = ADMETPredictor()
            result = self._admet.predict(smiles)
        except Exception as e:                       # model absent, weights missing, bad SMILES
            logger.warning(
                f"No hERG prediction available ({type(e).__name__}: {e}); pediatric cardiac "
                "liability will be reported as not assessed."
            )
            self._admet = False
            return None
        for key, value in (result.predictions or {}).items():
            if "herg" in key.lower():
                try:
                    return float(value)
                except (TypeError, ValueError):
                    return None
        return None

    def stage_8_ranking(self):
        """Stage 8: Score and rank candidates."""
        self._report_progress(8, "Ranking candidates")

        # Build molecule lookup
        mol_lookup = {m.id: m for m in self.molecules}
        dock_lookup = {d.molecule_id: d for d in self.docking_results}

        candidates = []
        for mol in self.molecules:
            dock_result = dock_lookup.get(mol.id)

            # Calculate composite score
            gen_score = mol.generation_score
            qed_score = mol.properties.qed if mol.properties else 0.5

            # Normalize the docking term ACCORDING TO ITS UNITS. Until 2026-11-03 every
            # result was normalized as if it were a kcal/mol affinity (`-score / 12`), which
            # silently inverted the sign of a DiffDock confidence (higher = better) and gave
            # it a unit it never had. See DockingResult.score_type.
            dock_score = dock_result.docking_score if dock_result else None
            if dock_result is None:
                dock_normalized = 0.0  # No docking data = worst score
            elif dock_result.score_type == "diffdock_confidence":
                # `confidence` is already P(pose within 2 A RMSD) in [0, 1].
                dock_normalized = max(0.0, min(1.0, dock_result.confidence))
            elif dock_score is not None:
                # Binding energy, kcal/mol: -12 (excellent) → 1.0, 0 (no binding) → 0.0
                dock_normalized = max(0.0, min(1.0, -dock_score / 12.0))
            else:
                dock_normalized = 0.0  # scored by no quantity we can rank on

            composite = (
                self.config.generation_weight * gen_score +
                self.config.docking_weight * dock_normalized +
                self.config.qed_weight * qed_score
            )

            candidates.append(RankedCandidate(
                rank=1,  # Temporary value, updated after sorting
                molecule_id=mol.id,
                smiles=mol.smiles,
                name=mol.name,
                target_gene=mol.target_gene,
                score_type=dock_result.score_type if dock_result else "binding_energy",
                docking_score=dock_score,
                confidence_score=dock_result.confidence_score if dock_result else None,
                generation_score=gen_score,
                qed_score=qed_score,
                composite_score=round(composite, 4),
                properties=mol.properties,
                binding_residues=dock_result.contacts if dock_result else [],
            ))

        # Run pediatric safety assessment on each candidate
        for c in candidates:
            mol_data = {}
            if c.properties:
                mol_data["molecular_weight"] = c.properties.molecular_weight
                mol_data["logP"] = c.properties.logP
                mol_data["tpsa"] = c.properties.tpsa
                mol_data["rotatable_bonds"] = c.properties.rotatable_bonds
                # The hERG input used to be invented here:
                #     if logp > 3.5 and mw > 350:
                #         estimated_herg = max(0.5, 30 - (logp*4) - (mw-350)*0.02)
                # cited as "Ref: Aronov (2005) J Med Chem 48:1341-1348". Aronov 2005 publishes
                # no such equation, the expression is not a hERG model, and the gate meant most
                # candidates got no value at all and were then declared pediatric-safe by
                # default. Ask the real ADMET model; if it is unavailable, pass nothing and let
                # the assessment return "not assessed".
                herg = self._predict_herg(c.smiles)
                if herg is not None:
                    mol_data["herg_inhibition_probability"] = herg
                # Screen for known teratogenic structural alerts via RDKit
                # Thalidomide-like glutarimides, retinoids, warfarin scaffolds
                try:
                    from rdkit import Chem
                    smiles = c.smiles
                    if smiles:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            teratogenic_smarts = [
                                "[#7]1C(=O)CCC1=O",       # Glutarimide (thalidomide class)
                                "C=C(C)c1ccc(cc1)C(=O)O", # Retinoid acid scaffold
                            ]
                            for smarts in teratogenic_smarts:
                                pat = Chem.MolFromSmarts(smarts)
                                if pat and mol.HasSubstructMatch(pat):
                                    mol_data["teratogenic_risk"] = True
                                    break
                except Exception:
                    pass
            c.pediatric_safety = pediatric_safety_assessment(mol_data)
            # `pediatric_safe` is three-valued now: `not safe` would have swept "not assessed"
            # in with "critical flag raised", which are different claims.
            if c.pediatric_safety.pediatric_safe is False:
                c.alerts.append("PEDIATRIC_CRITICAL_FLAG")
            elif c.pediatric_safety.pediatric_safe is None:
                c.alerts.append("PEDIATRIC_NOT_ASSESSED")
            elif c.pediatric_safety.safety_flags:
                c.alerts.append("PEDIATRIC_REVIEW_NEEDED")

        # Sort by composite score (higher is better)
        candidates.sort(key=lambda x: x.composite_score, reverse=True)

        # Assign ranks
        for i, c in enumerate(candidates):
            c.rank = i + 1

        self.ranked_candidates = candidates[:self.config.top_n_candidates]
        self.run.top_candidates = [c.molecule_id for c in self.ranked_candidates]

        self._report_progress(8, f"Top {len(self.ranked_candidates)} candidates ranked")
        self.run.current_stage = 8
        self.run.stages_completed.append(8)
        self._save_checkpoint(8)

    def stage_9_reporting(self):
        """Stage 9: Generate final report."""
        self._report_progress(9, "Generating report")

        report = {
            "run_id": self.run_id,
            "target": self.target.dict() if self.target else None,
            "config": self.config.dict(),
            "summary": {
                "molecules_generated": len(self.molecules),
                "molecules_docked": len(self.docking_results),
                "top_n": self.config.top_n_candidates,
                # `pediatric_flagged_count` used to be "not pediatric_safe", which counted
                # every candidate whose cardiac liability was never evaluated as flagged.
                # Three states, three counts.
                "pediatric_safe_count": len([
                    c for c in self.ranked_candidates
                    if c.pediatric_safety and c.pediatric_safety.pediatric_safe is True
                ]),
                "pediatric_flagged_count": len([
                    c for c in self.ranked_candidates
                    if c.pediatric_safety and c.pediatric_safety.pediatric_safe is False
                ]),
                "pediatric_not_assessed_count": len([
                    c for c in self.ranked_candidates
                    if c.pediatric_safety and c.pediatric_safety.pediatric_safe is None
                ]),
            },
            "top_candidates": [c.dict() for c in self.ranked_candidates],
            "structures_used": [s.dict() for s in (self.structures.structures if self.structures else [])],
            "stage_timings": self.stage_timings,
            "completed_at": datetime.now().isoformat(),
        }

        # Save report
        report_file = self.output_dir / f"report_{self.run_id}.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2, default=str)

        self._report_progress(9, f"Report saved to {report_file}")
        self.run.current_stage = 9
        self.run.stages_completed.append(9)
        self._save_checkpoint(9)

        return report


def run_vcp_demo_pipeline(
    output_dir: str = "outputs",
    resume_from: str = None,
    structure_cache: str = "data/structures",
) -> PipelineRun:
    """
    Run the VCP FTD demo pipeline.

    Convenience function for the demo narrative.

    Args:
        output_dir: Directory for pipeline outputs
        resume_from: Run ID to resume from (uses latest checkpoint)
        structure_cache: Where to look for (and cache) PDB coordinates. Stage 7 refuses to dock
            without them, so this is injectable rather than hard-coded.
    """
    config = PipelineConfig(
        target_gene="VCP",
        reference_compound_smiles="CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5",
        num_molecules=20,
        output_dir=output_dir,
        structure_cache=structure_cache,
    )

    target = TargetHypothesis(
        gene="VCP",
        protein="Valosin-containing protein (p97)",
        uniprot_id="P55072",
        rationale="VCP mutations cause FTD-ALS. CB-5083 is a known VCP inhibitor.",
        therapeutic_area="Neurodegeneration",
        mechanism="AAA+ ATPase inhibition",
        confidence="high",
        priority=5,
        pdb_ids=["5FTK", "8OOI", "9DIL"],
        status="validated",
    )

    pipeline = DrugDiscoveryPipeline(config)
    return pipeline.run_pipeline(target, resume_from_run_id=resume_from)
