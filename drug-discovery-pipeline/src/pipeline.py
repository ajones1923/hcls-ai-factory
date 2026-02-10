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
)
from .nim_clients import create_nim_clients, NIMServiceManager
from .checkpoint import CheckpointManager


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

        # Checkpoint manager
        self.checkpoint_manager = CheckpointManager(self.output_dir)

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
            if resume_stage < 0:
                self.stage_0_initialize(target)

            if resume_stage < 1:
                self.stage_1_normalize_target()

            if resume_stage < 2:
                self.stage_2_structure_discovery()

            if resume_stage < 3:
                self.stage_3_structure_prep()

            if resume_stage < 4:
                self.stage_4_molecule_generation()

            if resume_stage < 5:
                self.stage_5_chemistry_qc()

            if resume_stage < 6:
                self.stage_6_conformers()

            if resume_stage < 7:
                self.stage_7_docking()

            if resume_stage < 8:
                self.stage_8_ranking()

            if resume_stage < 9:
                self.stage_9_reporting()

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

    def stage_3_structure_prep(self):
        """Stage 3: Prepare structures for docking."""
        self._report_progress(3, "Preparing structures for docking")

        # In production, would:
        # 1. Download PDB files
        # 2. Remove water, ions
        # 3. Add hydrogens
        # 4. Identify binding site

        for struct in self.structures.structures:
            struct.prepared = True

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

        # Use primary structure
        primary = self.structures.primary_structure or self.structures.structures[0].pdb_id

        # Dock each molecule
        dock_failures = 0
        for mol in self.molecules:
            try:
                poses = self.nim.diffdock.dock(
                    protein_pdb=primary,
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
                    docking_score=best_pose.get("docking_score", 0),
                    confidence=best_pose.get("confidence", 0),
                    hydrogen_bonds=best_pose.get("hydrogen_bonds", 0),
                    contacts=best_pose.get("contacts", []),
                    method="DiffDock",
                ))

        self.run.molecules_docked = len(self.docking_results)
        self._report_progress(7, f"Docked {len(self.docking_results)} molecules")
        self.run.current_stage = 7
        self.run.stages_completed.append(7)
        self._save_checkpoint(7)

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
            dock_score = dock_result.docking_score if dock_result else 0
            qed_score = mol.properties.qed if mol.properties else 0.5

            # Normalize docking score (lower/more negative is better binding)
            # Range: -12 kcal/mol (excellent) → 1.0, 0 kcal/mol (no binding) → 0.0
            if dock_result is not None:
                dock_normalized = max(0.0, min(1.0, -dock_score / 12.0))
            else:
                dock_normalized = 0.0  # No docking data = worst score

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
                docking_score=dock_score,
                generation_score=gen_score,
                qed_score=qed_score,
                composite_score=round(composite, 4),
                properties=mol.properties,
                binding_residues=dock_result.contacts if dock_result else [],
            ))

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
            },
            "top_candidates": [c.dict() for c in self.ranked_candidates],
            "structures_used": [s.dict() for s in (self.structures.structures if self.structures else [])],
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
) -> PipelineRun:
    """
    Run the VCP FTD demo pipeline.

    Convenience function for the demo narrative.

    Args:
        output_dir: Directory for pipeline outputs
        resume_from: Run ID to resume from (uses latest checkpoint)
    """
    config = PipelineConfig(
        target_gene="VCP",
        reference_compound_smiles="CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5",
        num_molecules=20,
        output_dir=output_dir,
        structure_cache="data/structures",
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
