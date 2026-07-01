"""
Unified clinical report generation for the HCLS AI Factory precision medicine platform.

Generates comprehensive clinical PDF reports combining outputs from all three
pipeline stages (Genomics, RAG/Chat, Drug Discovery) plus CAR-T Intelligence
into a 5-section unified report:

    1. Genomic Profile     -- Patient variant summary, pathogenic variants, annotations
    2. Target Analysis     -- Drug target rationale, protein structures, druggability
    3. CAR-T Evaluation    -- Cell therapy options, trials, safety, manufacturing
    4. Drug Candidates     -- Ranked molecules from MolMIM/DiffDock/RDKit
    5. Treatment Recommendation -- Meta-agent synthesized clinical recommendation

Output formats: Markdown, JSON, PDF (reportlab with NVIDIA branding).

Section 11.3.3 / Section 10.3 Gap 6 of the HCLS AI Factory architecture.
"""

from __future__ import annotations

import json
import logging
import os
import re
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Brand colors and constants
# ---------------------------------------------------------------------------

NVIDIA_GREEN = "#76B900"
DARK_GRAY = "#333333"
LIGHT_GRAY = "#F5F5F5"
MEDIUM_GRAY = "#666666"
WHITE = "#FFFFFF"
TABLE_HEADER_BG = "#2d2d44"
TABLE_ALT_ROW = "#F8F8FC"
PATHOGENIC_RED = "#C41E3A"
STRUCTURE_BLUE = "#0066CC"

REPORT_TITLE = "HCLS AI Factory -- Precision Medicine Pipeline"
REPORT_SUBTITLE = "Unified Clinical Report"
REPORT_FOOTER = "HCLS AI Factory -- Precision Medicine Pipeline"

# PDF page dimensions (letter)
PAGE_WIDTH_INCH = 8.5
PAGE_HEIGHT_INCH = 11.0
MARGIN_INCH = 0.75

# Maximum items to include in tables
MAX_VARIANTS_IN_REPORT = 30
MAX_DRUG_CANDIDATES_IN_REPORT = 20
MAX_STRUCTURES_IN_REPORT = 10
MAX_TRIALS_IN_REPORT = 15


# ---------------------------------------------------------------------------
# Data classes for report input
# ---------------------------------------------------------------------------

@dataclass
class VariantRecord:
    """A genomic variant for inclusion in the report."""

    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""
    gene: str = ""
    consequence: str = ""
    impact: str = ""
    clinical_significance: str = ""
    rsid: str = ""
    am_pathogenicity: float | None = None
    am_class: str = ""
    quality: float = 0.0
    genotype: str = ""
    text_summary: str = ""
    disease_associations: str = ""

    @property
    def variant_id(self) -> str:
        return f"{self.chrom}_{self.pos}_{self.ref}_{self.alt}"

    @property
    def is_pathogenic(self) -> bool:
        sig = (self.clinical_significance or "").lower()
        return "pathogenic" in sig and "likely_benign" not in sig and "benign" not in sig

    @property
    def is_high_impact(self) -> bool:
        return (self.impact or "").upper() == "HIGH"

    def to_dict(self) -> dict[str, Any]:
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "gene": self.gene,
            "consequence": self.consequence,
            "impact": self.impact,
            "clinical_significance": self.clinical_significance,
            "rsid": self.rsid,
            "am_pathogenicity": self.am_pathogenicity,
            "am_class": self.am_class,
            "disease_associations": self.disease_associations,
        }


@dataclass
class TargetRecord:
    """A drug target for inclusion in the report."""

    gene: str = ""
    protein: str = ""
    uniprot_id: str = ""
    rationale: str = ""
    therapeutic_area: str = ""
    mechanism: str = ""
    confidence: str = "medium"
    priority: int = 0
    druggable: bool = False
    drug_status: str = ""
    known_drugs: list[str] = field(default_factory=list)
    diseases: list[str] = field(default_factory=list)
    pdb_ids: list[str] = field(default_factory=list)
    variant_count: int = 0
    function: str = ""
    pathway: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene": self.gene,
            "protein": self.protein,
            "uniprot_id": self.uniprot_id,
            "rationale": self.rationale,
            "therapeutic_area": self.therapeutic_area,
            "mechanism": self.mechanism,
            "confidence": self.confidence,
            "priority": self.priority,
            "druggable": self.druggable,
            "drug_status": self.drug_status,
            "known_drugs": self.known_drugs,
            "diseases": self.diseases,
            "pdb_ids": self.pdb_ids,
            "variant_count": self.variant_count,
        }


@dataclass
class StructureRecord:
    """A protein structure for inclusion in the report."""

    pdb_id: str = ""
    method: str = ""
    resolution: float | None = None
    chain: str = "A"
    description: str = ""
    binding_site_residues: list[str] = field(default_factory=list)
    ligand_id: str = ""
    ligand_smiles: str = ""
    is_primary: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "pdb_id": self.pdb_id,
            "method": self.method,
            "resolution": self.resolution,
            "chain": self.chain,
            "description": self.description,
            "binding_site_residues": self.binding_site_residues,
            "ligand_id": self.ligand_id,
            "is_primary": self.is_primary,
        }


@dataclass
class DrugCandidateRecord:
    """A ranked drug candidate for inclusion in the report."""

    rank: int = 0
    molecule_id: str = ""
    smiles: str = ""
    name: str = ""
    target_gene: str = ""
    docking_score: float = 0.0
    generation_score: float = 0.0
    qed_score: float | None = None
    composite_score: float = 0.0
    molecular_weight: float = 0.0
    logp: float = 0.0
    hbd: int = 0
    hba: int = 0
    tpsa: float = 0.0
    rotatable_bonds: int = 0
    lipinski_violations: int = 0
    best_pose_file: str = ""
    binding_residues: list[str] = field(default_factory=list)
    alerts: list[str] = field(default_factory=list)
    passes_filters: bool = True

    def to_dict(self) -> dict[str, Any]:
        return {
            "rank": self.rank,
            "molecule_id": self.molecule_id,
            "smiles": self.smiles,
            "name": self.name,
            "target_gene": self.target_gene,
            "docking_score": self.docking_score,
            "qed_score": self.qed_score,
            "composite_score": self.composite_score,
            "molecular_weight": self.molecular_weight,
            "logp": self.logp,
            "lipinski_violations": self.lipinski_violations,
            "passes_filters": self.passes_filters,
        }


@dataclass
class CARTEvaluationRecord:
    """CAR-T therapy evaluation data for inclusion in the report."""

    target_antigen: str = ""
    applicable: bool = False
    rationale: str = ""
    constructs: list[dict[str, Any]] = field(default_factory=list)
    clinical_trials: list[dict[str, Any]] = field(default_factory=list)
    safety_profile: dict[str, Any] = field(default_factory=dict)
    manufacturing: dict[str, Any] = field(default_factory=dict)
    approved_products: list[str] = field(default_factory=list)
    biomarkers: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "target_antigen": self.target_antigen,
            "applicable": self.applicable,
            "rationale": self.rationale,
            "constructs": self.constructs,
            "clinical_trials": self.clinical_trials,
            "safety_profile": self.safety_profile,
            "manufacturing": self.manufacturing,
            "approved_products": self.approved_products,
            "biomarkers": self.biomarkers,
        }


@dataclass
class ReportData:
    """Complete data for generating a unified clinical report."""

    # Metadata
    patient_id: str = ""
    target_gene: str = ""
    report_id: str = field(default_factory=lambda: str(uuid.uuid4())[:12])
    generated_at: str = field(default_factory=lambda: datetime.now().isoformat())
    pipeline_version: str = "0.1.0"

    # Section 1: Genomic Profile
    variants: list[VariantRecord] = field(default_factory=list)
    total_variants_analyzed: int = 0
    pathogenic_count: int = 0
    high_impact_count: int = 0

    # Section 2: Target Analysis
    targets: list[TargetRecord] = field(default_factory=list)
    structures: list[StructureRecord] = field(default_factory=list)
    primary_structure_pdb: str = ""

    # Section 3: CAR-T Evaluation
    cart_evaluation: CARTEvaluationRecord = field(
        default_factory=CARTEvaluationRecord
    )

    # Section 4: Drug Candidates
    drug_candidates: list[DrugCandidateRecord] = field(default_factory=list)
    reference_compound: str = ""
    reference_smiles: str = ""
    generation_method: str = "MolMIM"
    docking_method: str = "DiffDock"

    # Section 5: Treatment Recommendation
    recommendation: str = ""
    confidence: float = 0.0
    sources: list[str] = field(default_factory=list)
    follow_up_questions: list[str] = field(default_factory=list)

    # Pipeline metrics
    pipeline_metrics: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "metadata": {
                "patient_id": self.patient_id,
                "target_gene": self.target_gene,
                "report_id": self.report_id,
                "generated_at": self.generated_at,
                "pipeline_version": self.pipeline_version,
            },
            "genomic_profile": {
                "variants": [v.to_dict() for v in self.variants],
                "total_variants_analyzed": self.total_variants_analyzed,
                "pathogenic_count": self.pathogenic_count,
                "high_impact_count": self.high_impact_count,
            },
            "target_analysis": {
                "targets": [t.to_dict() for t in self.targets],
                "structures": [s.to_dict() for s in self.structures],
                "primary_structure_pdb": self.primary_structure_pdb,
            },
            "cart_evaluation": self.cart_evaluation.to_dict(),
            "drug_candidates": {
                "candidates": [d.to_dict() for d in self.drug_candidates],
                "reference_compound": self.reference_compound,
                "reference_smiles": self.reference_smiles,
                "generation_method": self.generation_method,
                "docking_method": self.docking_method,
            },
            "treatment_recommendation": {
                "recommendation": self.recommendation,
                "confidence": self.confidence,
                "sources": self.sources,
                "follow_up_questions": self.follow_up_questions,
            },
            "pipeline_metrics": self.pipeline_metrics,
        }


# ---------------------------------------------------------------------------
# Data collection helpers
# ---------------------------------------------------------------------------

class _DataCollector:
    """
    Internal helper that collects report data from Milvus, pipeline outputs,
    the knowledge base, and the meta-agent.
    """

    def __init__(
        self,
        milvus_host: str = "localhost",
        milvus_port: int = 19530,
        embedding_model: str = "BAAI/bge-small-en-v1.5",
        output_base_dir: str | Path = "outputs",
        knowledge: dict[str, Any] | None = None,
    ):
        self.milvus_host = milvus_host
        self.milvus_port = milvus_port
        self.embedding_model = embedding_model
        self.output_base_dir = Path(output_base_dir)
        self._knowledge = knowledge

    def _load_knowledge(self) -> dict[str, Any]:
        """Load the knowledge connections database."""
        if self._knowledge is not None:
            return self._knowledge

        try:
            import sys

            rag_src = Path(__file__).resolve().parent.parent.parent / "core/engines/precision-intelligence" / "src"
            if rag_src.exists():
                sys.path.insert(0, str(rag_src.parent))
                from src.knowledge import KNOWLEDGE_CONNECTIONS  # type: ignore

                self._knowledge = KNOWLEDGE_CONNECTIONS
                sys.path.pop(0)
            else:
                self._knowledge = {}
        except Exception:
            self._knowledge = {}

        return self._knowledge

    def collect_genomic_profile(
        self,
        target_gene: str,
        patient_id: str = "",
    ) -> tuple[list[VariantRecord], int, int, int]:
        """
        Query Milvus genomic_evidence for patient variants.

        Returns: (variants, total_analyzed, pathogenic_count, high_impact_count)
        """
        variants: list[VariantRecord] = []
        total = 0
        pathogenic = 0
        high_impact = 0

        try:
            from pymilvus import Collection, connections, utility

            connections.connect(
                alias="report_genomic",
                host=self.milvus_host,
                port=self.milvus_port,
            )

            if not utility.has_collection("genomic_evidence"):
                logger.warning("genomic_evidence collection not found")
                return variants, total, pathogenic, high_impact

            coll = Collection("genomic_evidence")
            coll.load()
            total = coll.num_entities

            # Query for target gene variants
            safe_gene = re.sub(r'[^A-Za-z0-9_\-]', '', target_gene)
            output_fields = [
                "chrom", "pos", "ref", "alt", "qual", "gene",
                "consequence", "impact", "genotype", "text_summary",
                "clinical_significance", "rsid", "am_pathogenicity", "am_class",
            ]

            results = coll.query(
                expr=f'gene == "{safe_gene}"',
                output_fields=output_fields,
                limit=MAX_VARIANTS_IN_REPORT,
            )

            for r in results:
                am_score = r.get("am_pathogenicity")
                if am_score == -1.0:
                    am_score = None

                v = VariantRecord(
                    chrom=r.get("chrom", ""),
                    pos=r.get("pos", 0),
                    ref=r.get("ref", ""),
                    alt=r.get("alt", ""),
                    gene=r.get("gene", ""),
                    consequence=r.get("consequence", ""),
                    impact=r.get("impact", ""),
                    clinical_significance=r.get("clinical_significance", ""),
                    rsid=r.get("rsid", ""),
                    am_pathogenicity=am_score,
                    am_class=r.get("am_class", ""),
                    quality=r.get("qual", 0.0),
                    genotype=r.get("genotype", ""),
                    text_summary=r.get("text_summary", ""),
                )
                variants.append(v)
                if v.is_pathogenic:
                    pathogenic += 1
                if v.is_high_impact:
                    high_impact += 1

            # Also query for pathogenic variants across all genes
            try:
                pathogenic_results = coll.query(
                    expr='impact == "HIGH"',
                    output_fields=["gene"],
                    limit=1000,
                )
                high_impact = len(pathogenic_results)
            except Exception:
                pass

            logger.info(
                "Collected %d variants for %s (%d pathogenic, %d high-impact of %d total)",
                len(variants), target_gene, pathogenic, high_impact, total,
            )

        except ImportError:
            logger.warning("pymilvus not available -- skipping genomic profile collection")
        except Exception as exc:
            logger.error("Failed to collect genomic profile: %s", exc)

        return variants, total, pathogenic, high_impact

    def collect_target_analysis(
        self,
        target_gene: str,
    ) -> tuple[list[TargetRecord], list[StructureRecord]]:
        """
        Collect target analysis data from knowledge base and pipeline outputs.

        Returns: (targets, structures)
        """
        targets: list[TargetRecord] = []
        structures: list[StructureRecord] = []

        knowledge = self._load_knowledge()
        safe_gene = re.sub(r'[^A-Za-z0-9_\-]', '', target_gene).upper()

        gene_info = knowledge.get(safe_gene, {})
        if gene_info:
            target = TargetRecord(
                gene=safe_gene,
                protein=gene_info.get("protein", ""),
                rationale=f"{gene_info.get('function', '')}. Pathway: {gene_info.get('pathway', '')}.",
                therapeutic_area=", ".join(gene_info.get("diseases", [])[:2]),
                druggable=gene_info.get("druggable", False),
                drug_status=gene_info.get("drug_status", ""),
                known_drugs=gene_info.get("drugs", []) if isinstance(gene_info.get("drugs"), list) else [gene_info.get("drugs", "")],
                diseases=gene_info.get("diseases", []),
                pdb_ids=gene_info.get("pdb_ids", []),
                function=gene_info.get("function", ""),
                pathway=gene_info.get("pathway", ""),
            )
            targets.append(target)

            # Create structure records from PDB IDs
            for i, pdb_id in enumerate(gene_info.get("pdb_ids", [])[:MAX_STRUCTURES_IN_REPORT]):
                struct = StructureRecord(
                    pdb_id=pdb_id,
                    is_primary=(i == 0),
                )
                structures.append(struct)

        # Check for pipeline structure manifests
        if self.output_base_dir.exists():
            for manifest_file in self.output_base_dir.rglob("structure_manifest.json"):
                try:
                    with open(manifest_file) as f:
                        manifest = json.load(f)
                    if manifest.get("target_gene", "").upper() == safe_gene:
                        for struct_data in manifest.get("structures", []):
                            struct = StructureRecord(
                                pdb_id=struct_data.get("pdb_id", ""),
                                method=struct_data.get("method", ""),
                                resolution=struct_data.get("resolution"),
                                chain=struct_data.get("chain", "A"),
                                description=struct_data.get("description", ""),
                                binding_site_residues=struct_data.get("binding_site_residues", []),
                                ligand_id=struct_data.get("ligand_id", ""),
                                ligand_smiles=struct_data.get("ligand_smiles", ""),
                            )
                            # Avoid duplicates
                            existing_ids = {s.pdb_id for s in structures}
                            if struct.pdb_id not in existing_ids:
                                structures.append(struct)
                except Exception:
                    pass

        logger.info(
            "Collected %d targets and %d structures for %s",
            len(targets), len(structures), target_gene,
        )
        return targets, structures

    def collect_cart_evaluation(
        self,
        target_gene: str,
    ) -> CARTEvaluationRecord:
        """
        Collect CAR-T evaluation data.

        CAR-T is primarily applicable for hematologic/oncologic targets with
        surface antigen expression. For non-applicable targets, returns a
        record with applicable=False.
        """
        # Known CAR-T target antigens
        cart_antigens = {
            "CD19", "BCMA", "CD22", "CD20", "CD33", "CD123", "HER2",
            "GD2", "EGFR", "MSLN", "PSMA", "MUC1", "EpCAM",
        }

        safe_gene = re.sub(r'[^A-Za-z0-9_\-]', '', target_gene).upper()

        if safe_gene not in cart_antigens:
            return CARTEvaluationRecord(
                target_antigen=safe_gene,
                applicable=False,
                rationale=(
                    f"{safe_gene} is not a cell surface antigen amenable to CAR-T targeting. "
                    f"Small molecule and other therapeutic modalities are more appropriate."
                ),
            )

        # Build CAR-T evaluation from knowledge
        evaluation = CARTEvaluationRecord(
            target_antigen=safe_gene,
            applicable=True,
            rationale=f"{safe_gene} is a validated CAR-T target antigen.",
        )

        # Try to get structured knowledge
        try:
            from hcls_common.meta_agent import CARTDataSource

            cart_source = CARTDataSource(
                milvus_host=self.milvus_host,
                milvus_port=self.milvus_port,
            )
            knowledge = cart_source.get_knowledge(target_antigen=safe_gene)
            data = knowledge.get("data", {})

            if isinstance(data, dict):
                targets_info = data.get("targets", {})
                evaluation.approved_products = targets_info.get("approved_products", [])
                evaluation.safety_profile = data.get("toxicities", {})
                evaluation.biomarkers = data.get("biomarkers", {})
                evaluation.manufacturing = data.get("manufacturing", {})
        except Exception as exc:
            logger.debug("Could not load CAR-T knowledge: %s", exc)

        # Try to search CAR-T Milvus collections
        try:
            from hcls_common.meta_agent import CARTDataSource

            cart_source = CARTDataSource(
                milvus_host=self.milvus_host,
                milvus_port=self.milvus_port,
            )
            evidence = cart_source.search_evidence(
                query=f"{safe_gene} CAR-T cell therapy",
                collections=["cart_clinical_trials"],
                target_antigen=safe_gene,
                top_k=5,
            )
            if "cart_clinical_trials" in evidence:
                evaluation.clinical_trials = evidence["cart_clinical_trials"]
        except Exception:
            pass

        return evaluation

    def collect_drug_candidates(
        self,
        target_gene: str,
    ) -> tuple[list[DrugCandidateRecord], str, str]:
        """
        Collect drug candidates from pipeline output.

        Returns: (candidates, reference_compound, reference_smiles)
        """
        candidates: list[DrugCandidateRecord] = []
        ref_compound = ""
        ref_smiles = ""

        safe_gene = re.sub(r'[^A-Za-z0-9_\-]', '', target_gene).upper()

        if not self.output_base_dir.exists():
            return candidates, ref_compound, ref_smiles

        # Search for ranked candidate files
        for json_file in sorted(
            self.output_base_dir.rglob("ranked_candidates.json"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        ):
            try:
                with open(json_file) as f:
                    data = json.load(f)

                cand_list = (
                    data if isinstance(data, list)
                    else data.get("candidates", data.get("ranked_candidates", []))
                )

                for cand_data in cand_list[:MAX_DRUG_CANDIDATES_IN_REPORT]:
                    props = cand_data.get("properties", {})
                    candidate = DrugCandidateRecord(
                        rank=cand_data.get("rank", 0),
                        molecule_id=cand_data.get("molecule_id", ""),
                        smiles=cand_data.get("smiles", ""),
                        name=cand_data.get("name", ""),
                        target_gene=cand_data.get("target_gene", safe_gene),
                        docking_score=cand_data.get("docking_score", 0.0),
                        generation_score=cand_data.get("generation_score", 0.0),
                        qed_score=cand_data.get("qed_score"),
                        composite_score=cand_data.get("composite_score", 0.0),
                        molecular_weight=props.get("molecular_weight", cand_data.get("molecular_weight", 0.0)),
                        logp=props.get("logP", cand_data.get("logp", 0.0)),
                        hbd=props.get("hbd", cand_data.get("hbd", 0)),
                        hba=props.get("hba", cand_data.get("hba", 0)),
                        tpsa=props.get("tpsa", cand_data.get("tpsa", 0.0)),
                        rotatable_bonds=props.get("rotatable_bonds", cand_data.get("rotatable_bonds", 0)),
                        lipinski_violations=props.get("lipinski_violations", cand_data.get("lipinski_violations", 0)),
                        best_pose_file=cand_data.get("best_pose_file", ""),
                        binding_residues=cand_data.get("binding_residues", []),
                        alerts=cand_data.get("alerts", []),
                        passes_filters=cand_data.get("passes_filters", True),
                    )
                    candidates.append(candidate)

                # Get reference compound info
                config = data.get("config", {}) if isinstance(data, dict) else {}
                ref_smiles = config.get("reference_compound_smiles", "")
                ref_compound = config.get("reference_compound", "")

                if candidates:
                    break  # Use most recent file

            except Exception as exc:
                logger.warning("Failed to read candidates from %s: %s", json_file, exc)

        logger.info("Collected %d drug candidates for %s", len(candidates), target_gene)
        return candidates, ref_compound, ref_smiles

    def collect_recommendation(
        self,
        target_gene: str,
        patient_id: str = "",
        variants: list[VariantRecord] | None = None,
        targets: list[TargetRecord] | None = None,
        candidates: list[DrugCandidateRecord] | None = None,
        cart_eval: CARTEvaluationRecord | None = None,
    ) -> tuple[str, float, list[str]]:
        """
        Generate treatment recommendation using Claude (meta-agent synthesis).

        Returns: (recommendation_text, confidence, follow_up_questions)
        """
        # Build context from collected data
        context_parts: list[str] = []

        if variants:
            pathogenic = [v for v in variants if v.is_pathogenic]
            context_parts.append(
                f"**Genomic findings:** {len(variants)} variants in {target_gene}, "
                f"{len(pathogenic)} pathogenic."
            )
            for v in pathogenic[:5]:
                context_parts.append(
                    f"  - {v.rsid or v.variant_id}: {v.consequence}, {v.clinical_significance}"
                )

        if targets:
            t = targets[0]
            context_parts.append(
                f"\n**Target:** {t.gene} ({t.protein}), druggable={t.druggable}, "
                f"status={t.drug_status}"
            )
            if t.diseases:
                context_parts.append(f"  Diseases: {', '.join(t.diseases[:3])}")

        if candidates:
            context_parts.append(
                f"\n**Drug candidates:** {len(candidates)} ranked molecules."
            )
            for c in candidates[:3]:
                context_parts.append(
                    f"  - Rank {c.rank}: docking={c.docking_score:.2f}, "
                    f"QED={c.qed_score or 0:.3f}, MW={c.molecular_weight:.1f}"
                )

        if cart_eval and cart_eval.applicable:
            context_parts.append(
                f"\n**CAR-T:** Applicable (target: {cart_eval.target_antigen}). "
                f"Approved products: {', '.join(cart_eval.approved_products) or 'None'}."
            )
        elif cart_eval:
            context_parts.append(
                f"\n**CAR-T:** Not applicable for {target_gene}. {cart_eval.rationale}"
            )

        context = "\n".join(context_parts)

        # Try to get LLM synthesis
        try:
            api_key = os.getenv("ANTHROPIC_API_KEY")
            if not api_key:
                return self._static_recommendation(
                    target_gene, variants, targets, candidates, cart_eval
                ), 0.5, []

            import anthropic

            client = anthropic.Anthropic(api_key=api_key, timeout=60)

            synthesis_prompt = f"""\
Based on the following clinical evidence for patient {patient_id or 'N/A'} with \
target gene {target_gene}, provide a concise treatment recommendation in markdown format.

{context}

Include:
1. Summary of key findings (2-3 sentences)
2. Recommended therapeutic approach (prioritized list)
3. Key considerations and risks
4. Suggested next steps
5. Confidence assessment (low/medium/high with rationale)

Format the response in clear markdown sections.
"""

            response = client.messages.create(
                model="claude-sonnet-4-20250514",
                max_tokens=2048,
                temperature=0.3,
                system=(
                    "You are a precision medicine advisor for the HCLS AI Factory. "
                    "Provide evidence-based treatment recommendations based on genomic, "
                    "structural, and drug discovery data."
                ),
                messages=[{"role": "user", "content": synthesis_prompt}],
            )

            recommendation = response.content[0].text
            confidence = 0.7  # Moderate confidence for AI-generated recommendation
            follow_ups = [
                "What are the ADMET properties of the top candidates?",
                "Are there ongoing clinical trials for this target?",
                "What is the predicted CNS penetration of these molecules?",
            ]

            return recommendation, confidence, follow_ups

        except ImportError:
            logger.warning("anthropic not available -- using static recommendation")
        except Exception as exc:
            logger.error("LLM recommendation failed: %s", exc)

        return self._static_recommendation(
            target_gene, variants, targets, candidates, cart_eval
        ), 0.5, []

    def _static_recommendation(
        self,
        target_gene: str,
        variants: list[VariantRecord] | None,
        targets: list[TargetRecord] | None,
        candidates: list[DrugCandidateRecord] | None,
        cart_eval: CARTEvaluationRecord | None,
    ) -> str:
        """Generate a static recommendation when LLM is not available."""
        parts = [f"## Treatment Recommendation for {target_gene}\n"]

        if variants:
            pathogenic = [v for v in variants if v.is_pathogenic]
            parts.append(
                f"**Genomic evidence** identified {len(pathogenic)} pathogenic variant(s) "
                f"in {target_gene} from {len(variants)} total variants analyzed."
            )

        if targets:
            t = targets[0]
            if t.druggable:
                parts.append(
                    f"\n**{target_gene} is a druggable target** ({t.drug_status}). "
                    f"Known drugs: {', '.join(t.known_drugs) if t.known_drugs else 'None approved'}."
                )
            else:
                parts.append(
                    f"\n**{target_gene} druggability assessment:** Currently limited. "
                    f"Alternative modalities may be required."
                )

        if candidates:
            best = candidates[0]
            parts.append(
                f"\n**Lead candidate:** Rank {best.rank} with docking score "
                f"{best.docking_score:.2f} kcal/mol, QED {best.qed_score or 0:.3f}."
            )

        if cart_eval:
            if cart_eval.applicable:
                parts.append(
                    f"\n**CAR-T therapy** is a viable option for {cart_eval.target_antigen}."
                )
            else:
                parts.append(
                    f"\n**CAR-T therapy** is not applicable for this target."
                )

        parts.append(
            "\n\n*Note: This recommendation was generated without LLM synthesis. "
            "Configure ANTHROPIC_API_KEY for AI-powered treatment recommendations.*"
        )

        return "\n".join(parts)


# ---------------------------------------------------------------------------
# ClinicalReportGenerator
# ---------------------------------------------------------------------------

class ClinicalReportGenerator:
    """
    Generates unified clinical reports combining outputs from all pipeline stages.

    Supports three output formats:
    - Markdown: Full markdown with tables, suitable for Streamlit/web display.
    - JSON: Structured JSON for programmatic consumption.
    - PDF: Professional PDF with NVIDIA branding via reportlab.

    Usage::

        generator = ClinicalReportGenerator()
        report_data = generator.generate_report(
            patient_id="HG002",
            target_gene="VCP",
        )
        # Save as PDF
        generator.save_pdf(report_data, "outputs/clinical_report.pdf")
        # Get markdown
        md = generator.to_markdown(report_data)
        # Get JSON
        js = generator.to_json(report_data)
    """

    def __init__(
        self,
        milvus_host: str = "localhost",
        milvus_port: int = 19530,
        embedding_model: str = "BAAI/bge-small-en-v1.5",
        output_base_dir: str | Path = "outputs",
        knowledge: dict[str, Any] | None = None,
    ):
        """
        Initialize the report generator.

        Args:
            milvus_host: Milvus vector database host.
            milvus_port: Milvus vector database port.
            embedding_model: Sentence transformer model for queries.
            output_base_dir: Base directory for pipeline outputs.
            knowledge: Optional pre-loaded knowledge connections dict.
        """
        self._collector = _DataCollector(
            milvus_host=milvus_host,
            milvus_port=milvus_port,
            embedding_model=embedding_model,
            output_base_dir=output_base_dir,
            knowledge=knowledge,
        )
        self.output_base_dir = Path(output_base_dir)

        logger.info("ClinicalReportGenerator initialized")

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def generate_report(
        self,
        patient_id: str,
        target_gene: str,
        include_cart: bool = True,
        include_recommendation: bool = True,
        variants: list[VariantRecord] | None = None,
        targets: list[TargetRecord] | None = None,
        structures: list[StructureRecord] | None = None,
        drug_candidates: list[DrugCandidateRecord] | None = None,
        cart_evaluation: CARTEvaluationRecord | None = None,
        recommendation: str | None = None,
    ) -> ReportData:
        """
        Assemble all 5 report sections.

        If data objects (variants, targets, etc.) are provided directly, they
        are used as-is. Otherwise, data is collected from Milvus, pipeline
        outputs, and the knowledge base.

        Args:
            patient_id: Patient identifier.
            target_gene: Primary target gene symbol.
            include_cart: Whether to include CAR-T evaluation section.
            include_recommendation: Whether to generate treatment recommendation.
            variants: Pre-collected variant records.
            targets: Pre-collected target records.
            structures: Pre-collected structure records.
            drug_candidates: Pre-collected drug candidate records.
            cart_evaluation: Pre-collected CAR-T evaluation.
            recommendation: Pre-written recommendation text.

        Returns:
            Populated ReportData object.
        """
        start_time = time.monotonic()
        report = ReportData(
            patient_id=patient_id,
            target_gene=target_gene,
        )

        # Section 1: Genomic Profile
        if variants is not None:
            report.variants = variants
            report.pathogenic_count = sum(1 for v in variants if v.is_pathogenic)
            report.high_impact_count = sum(1 for v in variants if v.is_high_impact)
        else:
            (
                report.variants,
                report.total_variants_analyzed,
                report.pathogenic_count,
                report.high_impact_count,
            ) = self._collector.collect_genomic_profile(target_gene, patient_id)

        logger.info("Section 1 (Genomic Profile): %d variants", len(report.variants))

        # Section 2: Target Analysis
        if targets is not None:
            report.targets = targets
        else:
            collected_targets, collected_structures = (
                self._collector.collect_target_analysis(target_gene)
            )
            report.targets = collected_targets
            if structures is None:
                report.structures = collected_structures

        if structures is not None:
            report.structures = structures

        if report.structures:
            primary = next(
                (s for s in report.structures if s.is_primary),
                report.structures[0],
            )
            report.primary_structure_pdb = primary.pdb_id

        logger.info(
            "Section 2 (Target Analysis): %d targets, %d structures",
            len(report.targets),
            len(report.structures),
        )

        # Section 3: CAR-T Evaluation
        if include_cart:
            if cart_evaluation is not None:
                report.cart_evaluation = cart_evaluation
            else:
                report.cart_evaluation = self._collector.collect_cart_evaluation(
                    target_gene
                )

            logger.info(
                "Section 3 (CAR-T): applicable=%s", report.cart_evaluation.applicable
            )

        # Section 4: Drug Candidates
        if drug_candidates is not None:
            report.drug_candidates = drug_candidates
        else:
            (
                report.drug_candidates,
                report.reference_compound,
                report.reference_smiles,
            ) = self._collector.collect_drug_candidates(target_gene)

        logger.info(
            "Section 4 (Drug Candidates): %d candidates", len(report.drug_candidates)
        )

        # Section 5: Treatment Recommendation
        if recommendation is not None:
            report.recommendation = recommendation
        elif include_recommendation:
            (
                report.recommendation,
                report.confidence,
                report.follow_up_questions,
            ) = self._collector.collect_recommendation(
                target_gene=target_gene,
                patient_id=patient_id,
                variants=report.variants,
                targets=report.targets,
                candidates=report.drug_candidates,
                cart_eval=report.cart_evaluation if include_cart else None,
            )

        elapsed_ms = (time.monotonic() - start_time) * 1000
        report.pipeline_metrics["report_generation_ms"] = round(elapsed_ms, 1)

        logger.info(
            "Report generated in %.1f ms (report_id=%s)",
            elapsed_ms,
            report.report_id,
        )
        return report

    # ------------------------------------------------------------------
    # Output format: Markdown
    # ------------------------------------------------------------------

    def to_markdown(self, report: ReportData) -> str:
        """
        Generate a full Markdown report with tables.

        Returns: Markdown string.
        """
        parts: list[str] = []

        # Title
        parts.append(f"# {REPORT_TITLE}")
        parts.append(f"## {REPORT_SUBTITLE}")
        parts.append("")
        parts.append(f"**Patient ID:** {report.patient_id}")
        parts.append(f"**Target Gene:** {report.target_gene}")
        parts.append(f"**Report ID:** {report.report_id}")
        parts.append(f"**Generated:** {report.generated_at}")
        parts.append("")
        parts.append("---")
        parts.append("")

        # Table of Contents
        parts.append("## Table of Contents")
        parts.append("1. [Genomic Profile](#1-genomic-profile)")
        parts.append("2. [Target Analysis](#2-target-analysis)")
        parts.append("3. [CAR-T Evaluation](#3-car-t-evaluation)")
        parts.append("4. [Small Molecule Candidates](#4-small-molecule-candidates)")
        parts.append("5. [Treatment Recommendation](#5-treatment-recommendation)")
        parts.append("")
        parts.append("---")
        parts.append("")

        # Section 1: Genomic Profile
        parts.append(self._section_genomic_profile_md(report))
        parts.append("")

        # Section 2: Target Analysis
        parts.append(self._section_target_analysis_md(report))
        parts.append("")

        # Section 3: CAR-T Evaluation
        parts.append(self._section_cart_evaluation_md(report))
        parts.append("")

        # Section 4: Drug Candidates
        parts.append(self._section_drug_candidates_md(report))
        parts.append("")

        # Section 5: Treatment Recommendation
        parts.append(self._section_treatment_recommendation_md(report))
        parts.append("")

        # Footer
        parts.append("---")
        parts.append(f"*{REPORT_FOOTER}*")
        parts.append(
            f"*Platform: NVIDIA DGX Spark | Parabricks 4.6 | BioNeMo NIM | Claude*"
        )

        return "\n".join(parts)

    def _section_genomic_profile_md(self, report: ReportData) -> str:
        """Markdown for Section 1: Genomic Profile."""
        parts: list[str] = []
        parts.append("## 1. Genomic Profile")
        parts.append("")

        if report.total_variants_analyzed > 0:
            parts.append(
                f"**Total variants analyzed:** {report.total_variants_analyzed:,}"
            )
        parts.append(f"**Pathogenic variants:** {report.pathogenic_count}")
        parts.append(f"**High-impact variants:** {report.high_impact_count}")
        parts.append("")

        if report.variants:
            # Variant table
            parts.append("### Variants in Target Gene")
            parts.append("")
            parts.append(
                "| Gene | Position | Ref/Alt | rsID | Consequence | Impact | ClinVar | AlphaMissense |"
            )
            parts.append(
                "|------|----------|---------|------|-------------|--------|---------|---------------|"
            )
            for v in report.variants[:MAX_VARIANTS_IN_REPORT]:
                am_str = (
                    f"{v.am_pathogenicity:.2f} ({v.am_class})"
                    if v.am_pathogenicity is not None
                    else "N/A"
                )
                parts.append(
                    f"| {v.gene} | {v.chrom}:{v.pos} | {v.ref}/{v.alt} | "
                    f"{v.rsid or 'N/A'} | {v.consequence} | {v.impact} | "
                    f"{v.clinical_significance} | {am_str} |"
                )
            parts.append("")

            # Pathogenicity summary
            pathogenic = [v for v in report.variants if v.is_pathogenic]
            if pathogenic:
                parts.append("### Pathogenic Variants")
                parts.append("")
                for v in pathogenic:
                    parts.append(
                        f"- **{v.rsid or v.variant_id}** ({v.gene}): "
                        f"{v.consequence}, {v.clinical_significance}"
                    )
                    if v.disease_associations:
                        parts.append(f"  - Diseases: {v.disease_associations}")
                parts.append("")
        else:
            parts.append(
                "*No variant data available. Connect to Milvus or provide variant records.*"
            )
            parts.append("")

        return "\n".join(parts)

    def _section_target_analysis_md(self, report: ReportData) -> str:
        """Markdown for Section 2: Target Analysis."""
        parts: list[str] = []
        parts.append("## 2. Target Analysis")
        parts.append("")

        if report.targets:
            for target in report.targets:
                parts.append(f"### {target.gene} ({target.protein})")
                parts.append("")

                # Target properties table
                parts.append("| Property | Value |")
                parts.append("|----------|-------|")
                parts.append(f"| Gene | {target.gene} |")
                parts.append(f"| Protein | {target.protein} |")
                if target.uniprot_id:
                    parts.append(f"| UniProt | {target.uniprot_id} |")
                if target.function:
                    parts.append(f"| Function | {target.function} |")
                if target.pathway:
                    parts.append(f"| Pathway | {target.pathway} |")
                parts.append(f"| Druggable | {'Yes' if target.druggable else 'No'} |")
                parts.append(f"| Drug Status | {target.drug_status} |")
                if target.known_drugs:
                    parts.append(
                        f"| Known Drugs | {', '.join(target.known_drugs)} |"
                    )
                parts.append("")

                # Disease associations
                if target.diseases:
                    parts.append("**Disease Associations:**")
                    for disease in target.diseases:
                        parts.append(f"- {disease}")
                    parts.append("")

                # Rationale
                if target.rationale:
                    parts.append(f"**Rationale:** {target.rationale}")
                    parts.append("")

        # Structures
        if report.structures:
            parts.append("### Protein Structures")
            parts.append("")
            parts.append("| PDB ID | Method | Resolution | Description | Primary |")
            parts.append("|--------|--------|------------|-------------|---------|")
            for s in report.structures[:MAX_STRUCTURES_IN_REPORT]:
                res_str = f"{s.resolution:.1f} A" if s.resolution else "N/A"
                primary_str = "Yes" if s.is_primary else ""
                parts.append(
                    f"| {s.pdb_id} | {s.method} | {res_str} | "
                    f"{s.description} | {primary_str} |"
                )
            parts.append("")

            if report.primary_structure_pdb:
                parts.append(
                    f"**Primary docking structure:** {report.primary_structure_pdb}"
                )
                parts.append("")

        if not report.targets and not report.structures:
            parts.append("*No target analysis data available.*")
            parts.append("")

        return "\n".join(parts)

    def _section_cart_evaluation_md(self, report: ReportData) -> str:
        """Markdown for Section 3: CAR-T Evaluation."""
        parts: list[str] = []
        parts.append("## 3. CAR-T Evaluation")
        parts.append("")

        cart = report.cart_evaluation
        if not cart.target_antigen:
            parts.append("*CAR-T evaluation not performed.*")
            parts.append("")
            return "\n".join(parts)

        if not cart.applicable:
            parts.append(
                f"**Status:** Not applicable for {cart.target_antigen}"
            )
            parts.append("")
            parts.append(f"**Rationale:** {cart.rationale}")
            parts.append("")
            return "\n".join(parts)

        parts.append(f"**Target Antigen:** {cart.target_antigen}")
        parts.append(f"**Status:** CAR-T therapy is applicable")
        parts.append("")

        if cart.approved_products:
            parts.append("### Approved CAR-T Products")
            for product in cart.approved_products:
                parts.append(f"- {product}")
            parts.append("")

        if cart.safety_profile:
            parts.append("### Safety Profile")
            parts.append("")
            parts.append("| Toxicity | Details |")
            parts.append("|----------|---------|")
            for key, value in cart.safety_profile.items():
                display_key = key.upper().replace("_", " ")
                parts.append(f"| {display_key} | {value} |")
            parts.append("")

        if cart.biomarkers:
            parts.append("### Biomarkers")
            for category, markers in cart.biomarkers.items():
                if isinstance(markers, list):
                    parts.append(f"**{category.title()}:** {', '.join(markers)}")
                else:
                    parts.append(f"**{category.title()}:** {markers}")
            parts.append("")

        if cart.manufacturing:
            parts.append("### Manufacturing Considerations")
            for key, value in cart.manufacturing.items():
                display_key = key.replace("_", " ").title()
                if isinstance(value, list):
                    parts.append(f"**{display_key}:**")
                    for item in value:
                        parts.append(f"- {item}")
                else:
                    parts.append(f"**{display_key}:** {value}")
            parts.append("")

        if cart.clinical_trials:
            parts.append("### Relevant Clinical Trials")
            parts.append("")
            for trial in cart.clinical_trials[:MAX_TRIALS_IN_REPORT]:
                if isinstance(trial, dict):
                    trial_id = trial.get("id", trial.get("trial_id", "N/A"))
                    parts.append(f"- **{trial_id}**: {json.dumps(trial, default=str)[:200]}")
            parts.append("")

        return "\n".join(parts)

    def _section_drug_candidates_md(self, report: ReportData) -> str:
        """Markdown for Section 4: Drug Candidates."""
        parts: list[str] = []
        parts.append("## 4. Small Molecule Candidates")
        parts.append("")

        if report.reference_compound:
            parts.append(f"**Reference Compound:** {report.reference_compound}")
        if report.reference_smiles:
            parts.append(f"**Reference SMILES:** `{report.reference_smiles}`")
        parts.append(f"**Generation Method:** {report.generation_method}")
        parts.append(f"**Docking Method:** {report.docking_method}")
        parts.append("")

        if report.drug_candidates:
            # Summary table
            parts.append("### Ranked Candidates")
            parts.append("")
            parts.append(
                "| Rank | Molecule | Docking (kcal/mol) | QED | MW (Da) | LogP | "
                "Lipinski Viol. | Composite |"
            )
            parts.append(
                "|------|----------|-------------------|-----|---------|------|"
                "----------------|-----------|"
            )
            for c in report.drug_candidates[:MAX_DRUG_CANDIDATES_IN_REPORT]:
                qed_str = f"{c.qed_score:.3f}" if c.qed_score is not None else "N/A"
                parts.append(
                    f"| {c.rank} | {c.molecule_id or c.name or 'N/A'} | "
                    f"{c.docking_score:.2f} | {qed_str} | {c.molecular_weight:.1f} | "
                    f"{c.logp:.2f} | {c.lipinski_violations} | {c.composite_score:.3f} |"
                )
            parts.append("")

            # Detailed candidate entries
            parts.append("### Candidate Details")
            parts.append("")
            for c in report.drug_candidates[:10]:
                parts.append(f"#### Candidate #{c.rank}")
                parts.append(f"**SMILES:** `{c.smiles}`")
                parts.append("")
                parts.append(
                    f"- Docking score: {c.docking_score:.2f} kcal/mol"
                )
                parts.append(
                    f"- QED: {c.qed_score:.3f}" if c.qed_score is not None else "- QED: N/A"
                )
                parts.append(f"- Molecular weight: {c.molecular_weight:.1f} Da")
                parts.append(f"- LogP: {c.logp:.2f}")
                parts.append(f"- H-bond donors/acceptors: {c.hbd}/{c.hba}")
                parts.append(f"- TPSA: {c.tpsa:.1f}")
                parts.append(f"- Rotatable bonds: {c.rotatable_bonds}")
                parts.append(
                    f"- Lipinski violations: {c.lipinski_violations}"
                )
                if c.alerts:
                    parts.append(f"- Alerts: {', '.join(c.alerts)}")
                parts.append("")
        else:
            parts.append(
                "*No drug candidates available. Run the drug discovery pipeline or provide candidate records.*"
            )
            parts.append("")

        return "\n".join(parts)

    def _section_treatment_recommendation_md(self, report: ReportData) -> str:
        """Markdown for Section 5: Treatment Recommendation."""
        parts: list[str] = []
        parts.append("## 5. Treatment Recommendation")
        parts.append("")

        if report.recommendation:
            parts.append(report.recommendation)
            parts.append("")
        else:
            parts.append("*No treatment recommendation generated.*")
            parts.append("")

        if report.confidence:
            parts.append(
                f"**Confidence:** {report.confidence:.0%}"
            )
            parts.append("")

        if report.sources:
            parts.append("**Data Sources:**")
            for src in report.sources:
                parts.append(f"- {src}")
            parts.append("")

        if report.follow_up_questions:
            parts.append("### Suggested Follow-Up Questions")
            for q in report.follow_up_questions:
                parts.append(f"- {q}")
            parts.append("")

        return "\n".join(parts)

    # ------------------------------------------------------------------
    # Output format: JSON
    # ------------------------------------------------------------------

    def to_json(self, report: ReportData, indent: int = 2) -> str:
        """
        Generate a structured JSON report.

        Returns: JSON string.
        """
        return json.dumps(report.to_dict(), indent=indent, default=str)

    def save_json(self, report: ReportData, output_path: str | Path) -> Path:
        """Save report as JSON file."""
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            f.write(self.to_json(report))
        logger.info("JSON report saved: %s", path)
        return path

    # ------------------------------------------------------------------
    # Output format: PDF (reportlab)
    # ------------------------------------------------------------------

    def to_pdf(self, report: ReportData) -> bytes:
        """
        Generate a PDF report with NVIDIA branding.

        Returns: PDF bytes.
        """
        from reportlab.lib import colors as rl_colors
        from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import inch
        from reportlab.platypus import (
            HRFlowable,
            PageBreak,
            Paragraph,
            SimpleDocTemplate,
            Spacer,
            Table,
            TableStyle,
        )

        buffer = BytesIO()

        doc = SimpleDocTemplate(
            buffer,
            pagesize=letter,
            rightMargin=MARGIN_INCH * inch,
            leftMargin=MARGIN_INCH * inch,
            topMargin=MARGIN_INCH * inch,
            bottomMargin=MARGIN_INCH * inch,
        )

        styles = getSampleStyleSheet()

        # Custom styles
        style_title = ParagraphStyle(
            "TitleMain",
            parent=styles["Heading1"],
            fontSize=22,
            spaceAfter=6,
            alignment=TA_CENTER,
            textColor=rl_colors.HexColor(DARK_GRAY),
        )
        style_subtitle = ParagraphStyle(
            "Subtitle",
            parent=styles["Normal"],
            fontSize=13,
            spaceAfter=16,
            alignment=TA_CENTER,
            textColor=rl_colors.HexColor(MEDIUM_GRAY),
        )
        style_section = ParagraphStyle(
            "SectionHeader",
            parent=styles["Heading2"],
            fontSize=15,
            spaceBefore=18,
            spaceAfter=10,
            textColor=rl_colors.HexColor(NVIDIA_GREEN),
            borderColor=rl_colors.HexColor(NVIDIA_GREEN),
            borderWidth=1,
            borderPadding=5,
        )
        style_subsection = ParagraphStyle(
            "SubSection",
            parent=styles["Heading3"],
            fontSize=12,
            spaceBefore=12,
            spaceAfter=6,
            textColor=rl_colors.HexColor(DARK_GRAY),
        )
        style_body = ParagraphStyle(
            "BodyJustified",
            parent=styles["Normal"],
            fontSize=10,
            spaceAfter=8,
            alignment=TA_JUSTIFY,
            leading=14,
        )
        style_small = ParagraphStyle(
            "SmallText",
            parent=styles["Normal"],
            fontSize=8,
            textColor=rl_colors.HexColor(MEDIUM_GRAY),
        )
        style_code = ParagraphStyle(
            "CodeText",
            parent=styles["Normal"],
            fontSize=7,
            fontName="Courier",
            backColor=rl_colors.HexColor(LIGHT_GRAY),
            borderColor=rl_colors.HexColor("#DDDDDD"),
            borderWidth=1,
            borderPadding=4,
        )

        # Color objects for table styling
        c_nvidia_green = rl_colors.HexColor(NVIDIA_GREEN)
        c_dark_gray = rl_colors.HexColor(DARK_GRAY)
        c_table_header = rl_colors.HexColor(TABLE_HEADER_BG)
        c_light_gray = rl_colors.HexColor(LIGHT_GRAY)
        c_alt_row = rl_colors.HexColor(TABLE_ALT_ROW)
        c_red = rl_colors.HexColor(PATHOGENIC_RED)
        c_blue = rl_colors.HexColor(STRUCTURE_BLUE)
        c_white = rl_colors.white
        c_border = rl_colors.HexColor("#DDDDDD")
        c_inner_border = rl_colors.HexColor("#EEEEEE")
        c_green_bg = rl_colors.HexColor("#F0F8E8")

        # Common table style base
        def _base_table_style(header_bg=c_table_header):
            return [
                ("BACKGROUND", (0, 0), (-1, 0), header_bg),
                ("TEXTCOLOR", (0, 0), (-1, 0), c_white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 9),
                ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [c_white, c_light_gray]),
                ("BOX", (0, 0), (-1, -1), 1, c_border),
                ("INNERGRID", (0, 0), (-1, -1), 0.5, c_inner_border),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
            ]

        elements: list[Any] = []

        # ---- Title Page ----
        elements.append(Spacer(1, 0.5 * inch))

        # Green header bar
        header_bar = Table(
            [[""]],
            colWidths=[6.5 * inch],
            rowHeights=[4],
        )
        header_bar.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), c_nvidia_green),
        ]))
        elements.append(header_bar)
        elements.append(Spacer(1, 16))

        elements.append(Paragraph(REPORT_TITLE.upper(), style_title))
        elements.append(Paragraph(REPORT_SUBTITLE, style_subtitle))
        elements.append(Spacer(1, 12))

        # Report info box
        info_data = [[
            Paragraph(
                f"<b>Patient ID:</b> {_esc(report.patient_id)}",
                style_body,
            ),
            Paragraph(
                f"<b>Target Gene:</b> {_esc(report.target_gene)}",
                style_body,
            ),
            Paragraph(
                f"<b>Date:</b> {datetime.now().strftime('%B %d, %Y')}",
                style_body,
            ),
        ]]
        info_table = Table(info_data, colWidths=[2.1 * inch, 2.1 * inch, 2.3 * inch])
        info_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), c_green_bg),
            ("BOX", (0, 0), (-1, -1), 1, c_nvidia_green),
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("TOPPADDING", (0, 0), (-1, -1), 8),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
        ]))
        elements.append(info_table)
        elements.append(Spacer(1, 12))

        # Report ID line
        elements.append(
            Paragraph(
                f"Report ID: {report.report_id} | Pipeline version: {report.pipeline_version}",
                style_small,
            )
        )
        elements.append(Spacer(1, 8))

        # Green divider
        elements.append(
            HRFlowable(width="100%", thickness=2, color=c_nvidia_green)
        )
        elements.append(Spacer(1, 8))

        # Table of Contents
        elements.append(Paragraph("Table of Contents", style_subsection))
        toc_items = [
            "1. Genomic Profile",
            "2. Target Analysis",
            "3. CAR-T Evaluation",
            "4. Small Molecule Candidates",
            "5. Treatment Recommendation",
        ]
        for item in toc_items:
            elements.append(Paragraph(item, style_body))
        elements.append(Spacer(1, 10))

        # ---- Section 1: Genomic Profile ----
        elements.append(
            Paragraph("1. GENOMIC PROFILE", style_section)
        )
        elements.append(Spacer(1, 6))

        if report.total_variants_analyzed > 0:
            elements.append(
                Paragraph(
                    f"<b>Total variants analyzed:</b> {report.total_variants_analyzed:,}",
                    style_body,
                )
            )
        elements.append(
            Paragraph(
                f"<b>Pathogenic variants:</b> {report.pathogenic_count} | "
                f"<b>High-impact variants:</b> {report.high_impact_count}",
                style_body,
            )
        )
        elements.append(Spacer(1, 8))

        if report.variants:
            elements.append(Paragraph("Variants in Target Gene", style_subsection))
            variant_header = [
                "Gene", "Position", "Ref/Alt", "rsID", "Consequence",
                "Impact", "ClinVar", "AM Score",
            ]
            variant_rows = [variant_header]
            for v in report.variants[:MAX_VARIANTS_IN_REPORT]:
                am_str = (
                    f"{v.am_pathogenicity:.2f}" if v.am_pathogenicity is not None else "N/A"
                )
                variant_rows.append([
                    _esc(v.gene),
                    f"{_esc(v.chrom)}:{v.pos}",
                    f"{_esc(v.ref)}/{_esc(v.alt)}",
                    _esc(v.rsid) or "N/A",
                    _esc(v.consequence)[:25],
                    _esc(v.impact),
                    _esc(v.clinical_significance)[:20],
                    am_str,
                ])

            col_widths = [0.6*inch, 1.0*inch, 0.6*inch, 0.8*inch, 1.0*inch, 0.6*inch, 0.9*inch, 0.6*inch]
            variant_table = Table(variant_rows, colWidths=col_widths)
            variant_table.setStyle(TableStyle(_base_table_style()))
            elements.append(variant_table)
            elements.append(Spacer(1, 12))

            # Pathogenicity highlights
            pathogenic = [v for v in report.variants if v.is_pathogenic]
            if pathogenic:
                elements.append(
                    Paragraph("Pathogenic Variant Highlights", style_subsection)
                )
                for v in pathogenic[:5]:
                    elements.append(
                        Paragraph(
                            f"<b>{_esc(v.rsid or v.variant_id)}</b> ({_esc(v.gene)}): "
                            f"{_esc(v.consequence)}, {_esc(v.clinical_significance)}"
                            + (f" -- {_esc(v.disease_associations)}" if v.disease_associations else ""),
                            style_body,
                        )
                    )
                elements.append(Spacer(1, 8))

        # ---- Section 2: Target Analysis ----
        elements.append(PageBreak())
        elements.append(
            Paragraph("2. TARGET ANALYSIS", style_section)
        )
        elements.append(Spacer(1, 6))

        if report.targets:
            for target in report.targets:
                elements.append(
                    Paragraph(
                        f"{_esc(target.gene)} ({_esc(target.protein)})",
                        style_subsection,
                    )
                )

                target_rows = [
                    ["Property", "Value"],
                    ["Gene", _esc(target.gene)],
                    ["Protein", _esc(target.protein)],
                ]
                if target.function:
                    target_rows.append(["Function", _esc(target.function)])
                if target.pathway:
                    target_rows.append(["Pathway", _esc(target.pathway)])
                target_rows.append(
                    ["Druggable", "Yes" if target.druggable else "No"]
                )
                target_rows.append(["Drug Status", _esc(target.drug_status)])
                if target.known_drugs:
                    target_rows.append(
                        ["Known Drugs", _esc(", ".join(target.known_drugs))]
                    )
                if target.diseases:
                    target_rows.append(
                        ["Diseases", _esc(", ".join(target.diseases[:4]))]
                    )

                target_table = Table(
                    target_rows,
                    colWidths=[2.0 * inch, 4.5 * inch],
                )
                target_table.setStyle(TableStyle(_base_table_style()))
                elements.append(target_table)
                elements.append(Spacer(1, 10))

                if target.rationale:
                    elements.append(
                        Paragraph(
                            f"<b>Rationale:</b> {_esc(target.rationale)}",
                            style_body,
                        )
                    )
                    elements.append(Spacer(1, 6))

        # Structures table
        if report.structures:
            elements.append(
                Paragraph("Protein Structures", style_subsection)
            )
            struct_header = ["PDB ID", "Method", "Resolution", "Description", "Primary"]
            struct_rows = [struct_header]
            for s in report.structures[:MAX_STRUCTURES_IN_REPORT]:
                res_str = f"{s.resolution:.1f} A" if s.resolution else "N/A"
                struct_rows.append([
                    s.pdb_id,
                    s.method or "N/A",
                    res_str,
                    _esc(s.description)[:35] or "",
                    "Yes" if s.is_primary else "",
                ])

            struct_table = Table(
                struct_rows,
                colWidths=[0.8*inch, 0.8*inch, 0.8*inch, 2.5*inch, 0.7*inch],
            )
            struct_table.setStyle(TableStyle(_base_table_style(c_blue)))
            elements.append(struct_table)
            elements.append(Spacer(1, 12))

        # ---- Section 3: CAR-T Evaluation ----
        elements.append(
            Paragraph("3. CAR-T EVALUATION", style_section)
        )
        elements.append(Spacer(1, 6))

        cart = report.cart_evaluation
        if not cart.target_antigen:
            elements.append(
                Paragraph(
                    "<i>CAR-T evaluation not performed.</i>",
                    style_body,
                )
            )
        elif not cart.applicable:
            elements.append(
                Paragraph(
                    f"<b>Status:</b> Not applicable for {_esc(cart.target_antigen)}",
                    style_body,
                )
            )
            elements.append(
                Paragraph(
                    f"<b>Rationale:</b> {_esc(cart.rationale)}",
                    style_body,
                )
            )
        else:
            elements.append(
                Paragraph(
                    f"<b>Target Antigen:</b> {_esc(cart.target_antigen)} | "
                    f"<b>Status:</b> CAR-T therapy applicable",
                    style_body,
                )
            )

            if cart.approved_products:
                elements.append(Spacer(1, 6))
                elements.append(
                    Paragraph("Approved Products", style_subsection)
                )
                for product in cart.approved_products:
                    elements.append(
                        Paragraph(f"  - {_esc(product)}", style_body)
                    )

            if cart.safety_profile:
                elements.append(Spacer(1, 6))
                elements.append(
                    Paragraph("Safety Profile", style_subsection)
                )
                safety_rows = [["Toxicity", "Details"]]
                for key, value in cart.safety_profile.items():
                    safety_rows.append([
                        key.upper().replace("_", " "),
                        _esc(str(value))[:80],
                    ])
                safety_table = Table(
                    safety_rows,
                    colWidths=[1.5 * inch, 5.0 * inch],
                )
                safety_table.setStyle(TableStyle(_base_table_style(c_red)))
                elements.append(safety_table)

        elements.append(Spacer(1, 12))

        # ---- Section 4: Drug Candidates ----
        elements.append(PageBreak())
        elements.append(
            Paragraph("4. SMALL MOLECULE CANDIDATES", style_section)
        )
        elements.append(Spacer(1, 6))

        if report.reference_compound:
            elements.append(
                Paragraph(
                    f"<b>Reference Compound:</b> {_esc(report.reference_compound)}",
                    style_body,
                )
            )
        if report.reference_smiles:
            elements.append(
                Paragraph(
                    f"<b>Reference SMILES:</b> <font face='Courier' size='7'>"
                    f"{_esc(report.reference_smiles)}</font>",
                    style_body,
                )
            )
        elements.append(
            Paragraph(
                f"<b>Generation:</b> {_esc(report.generation_method)} | "
                f"<b>Docking:</b> {_esc(report.docking_method)}",
                style_body,
            )
        )
        elements.append(Spacer(1, 8))

        if report.drug_candidates:
            # Ranking table
            elements.append(
                Paragraph("Top Ranked Candidates", style_subsection)
            )
            cand_header = [
                "Rank", "Docking\n(kcal/mol)", "QED", "MW (Da)",
                "LogP", "Lipinski\nViol.", "Composite",
            ]
            cand_rows = [cand_header]
            for c in report.drug_candidates[:MAX_DRUG_CANDIDATES_IN_REPORT]:
                qed_str = f"{c.qed_score:.3f}" if c.qed_score is not None else "N/A"
                cand_rows.append([
                    str(c.rank),
                    f"{c.docking_score:.2f}",
                    qed_str,
                    f"{c.molecular_weight:.1f}",
                    f"{c.logp:.2f}",
                    str(c.lipinski_violations),
                    f"{c.composite_score:.3f}",
                ])

            cand_col_widths = [0.5*inch, 0.9*inch, 0.6*inch, 0.7*inch, 0.6*inch, 0.7*inch, 0.8*inch]
            cand_table = Table(cand_rows, colWidths=cand_col_widths)
            cand_style = _base_table_style(c_nvidia_green)
            cand_style.append(("ALIGN", (0, 0), (-1, -1), "CENTER"))
            cand_table.setStyle(TableStyle(cand_style))
            elements.append(cand_table)
            elements.append(Spacer(1, 12))

            # Individual candidate details
            for c in report.drug_candidates[:10]:
                elements.append(
                    Paragraph(
                        f"<b>Candidate #{c.rank}</b>",
                        style_subsection,
                    )
                )
                # SMILES (monospace, wrapped)
                smiles_display = _esc(c.smiles)
                if len(smiles_display) > 90:
                    smiles_display = smiles_display[:90] + "..."
                elements.append(
                    Paragraph(
                        f"<font face='Courier' size='7'>{smiles_display}</font>",
                        style_small,
                    )
                )
                elements.append(Spacer(1, 4))

                # Properties row
                prop_data = [[
                    f"Docking: {c.docking_score:.2f}",
                    f"QED: {c.qed_score:.3f}" if c.qed_score is not None else "QED: N/A",
                    f"MW: {c.molecular_weight:.1f}",
                    f"LogP: {c.logp:.2f}",
                    f"Score: {c.composite_score:.3f}",
                ]]
                prop_table = Table(prop_data, colWidths=[1.1*inch] * 5)
                prop_table.setStyle(TableStyle([
                    ("FONTSIZE", (0, 0), (-1, -1), 8),
                    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                    ("BACKGROUND", (0, 0), (-1, -1), c_light_gray),
                    ("BOX", (0, 0), (-1, -1), 0.5, c_border),
                    ("TOPPADDING", (0, 0), (-1, -1), 4),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                ]))
                elements.append(prop_table)
                elements.append(Spacer(1, 6))
        else:
            elements.append(
                Paragraph(
                    "<i>No drug candidates available. Run the drug discovery pipeline.</i>",
                    style_body,
                )
            )

        # ---- Section 5: Treatment Recommendation ----
        elements.append(PageBreak())
        elements.append(
            Paragraph("5. TREATMENT RECOMMENDATION", style_section)
        )
        elements.append(Spacer(1, 6))

        if report.recommendation:
            # Convert markdown to paragraphs
            for line in report.recommendation.split("\n"):
                line = line.strip()
                if not line:
                    elements.append(Spacer(1, 4))
                elif line.startswith("###"):
                    elements.append(
                        Paragraph(_esc(line.lstrip("#").strip()), style_subsection)
                    )
                elif line.startswith("##"):
                    elements.append(
                        Paragraph(_esc(line.lstrip("#").strip()), style_subsection)
                    )
                elif line.startswith("#"):
                    elements.append(
                        Paragraph(_esc(line.lstrip("#").strip()), style_subsection)
                    )
                elif line.startswith("- ") or line.startswith("* "):
                    elements.append(
                        Paragraph(
                            f"  - {_esc(line[2:])}",
                            style_body,
                        )
                    )
                elif line.startswith(("1.", "2.", "3.", "4.", "5.", "6.", "7.", "8.", "9.")):
                    elements.append(Paragraph(_esc(line), style_body))
                else:
                    # Convert bold markdown
                    converted = re.sub(
                        r'\*\*(.+?)\*\*', r'<b>\1</b>', _esc(line)
                    )
                    elements.append(Paragraph(converted, style_body))
        else:
            elements.append(
                Paragraph(
                    "<i>No treatment recommendation generated.</i>",
                    style_body,
                )
            )

        if report.confidence:
            elements.append(Spacer(1, 8))
            elements.append(
                Paragraph(
                    f"<b>Confidence:</b> {report.confidence:.0%}",
                    style_body,
                )
            )

        if report.follow_up_questions:
            elements.append(Spacer(1, 10))
            elements.append(
                Paragraph("Suggested Follow-Up Questions", style_subsection)
            )
            for q in report.follow_up_questions:
                elements.append(
                    Paragraph(f"  - {_esc(q)}", style_body)
                )

        # ---- Footer ----
        elements.append(Spacer(1, 30))
        elements.append(
            HRFlowable(width="100%", thickness=2, color=c_nvidia_green)
        )
        elements.append(Spacer(1, 8))
        elements.append(
            Paragraph(
                f"<b>{REPORT_FOOTER}</b><br/>"
                f"Generated: {datetime.now().strftime('%B %d, %Y at %H:%M')}<br/>"
                f"Platform: NVIDIA DGX Spark | Parabricks 4.6 | BioNeMo NIM | Claude<br/>"
                f"Report ID: {report.report_id}<br/>"
                f"<br/>"
                f"<font color='{NVIDIA_GREEN}'>Powered by NVIDIA Accelerated Computing</font>",
                style_small,
            )
        )

        # Green footer bar
        footer_bar = Table(
            [[""]],
            colWidths=[6.5 * inch],
            rowHeights=[4],
        )
        footer_bar.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), c_nvidia_green),
        ]))
        elements.append(Spacer(1, 8))
        elements.append(footer_bar)

        # Build the PDF
        doc.build(elements)
        return buffer.getvalue()

    def save_pdf(self, report: ReportData, output_path: str | Path) -> Path:
        """
        Save report as a branded PDF file.

        Args:
            report: Populated ReportData.
            output_path: Output file path.

        Returns:
            Path to the generated PDF file.
        """
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)

        pdf_bytes = self.to_pdf(report)
        with open(path, "wb") as f:
            f.write(pdf_bytes)

        logger.info("PDF report saved: %s (%d bytes)", path, len(pdf_bytes))
        return path

    # ------------------------------------------------------------------
    # Convenience: generate + save all formats
    # ------------------------------------------------------------------

    def generate_and_save(
        self,
        patient_id: str,
        target_gene: str,
        output_dir: str | Path | None = None,
        formats: Sequence[str] = ("markdown", "json", "pdf"),
        **kwargs: Any,
    ) -> dict[str, Path]:
        """
        Generate report and save in requested formats.

        Args:
            patient_id: Patient identifier.
            target_gene: Primary target gene symbol.
            output_dir: Output directory (defaults to self.output_base_dir / reports).
            formats: Tuple of format names to generate ("markdown", "json", "pdf").
            **kwargs: Additional arguments passed to generate_report().

        Returns:
            Dict mapping format name to output file path.
        """
        output_dir = Path(output_dir) if output_dir else self.output_base_dir / "reports"
        output_dir.mkdir(parents=True, exist_ok=True)

        report = self.generate_report(
            patient_id=patient_id,
            target_gene=target_gene,
            **kwargs,
        )

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = f"{target_gene}_{patient_id}_{timestamp}"
        saved: dict[str, Path] = {}

        if "markdown" in formats:
            md_path = output_dir / f"{base_name}.md"
            with open(md_path, "w") as f:
                f.write(self.to_markdown(report))
            saved["markdown"] = md_path
            logger.info("Markdown report saved: %s", md_path)

        if "json" in formats:
            json_path = self.save_json(report, output_dir / f"{base_name}.json")
            saved["json"] = json_path

        if "pdf" in formats:
            try:
                pdf_path = self.save_pdf(report, output_dir / f"{base_name}.pdf")
                saved["pdf"] = pdf_path
            except ImportError:
                logger.warning(
                    "reportlab not installed -- skipping PDF generation. "
                    "Install with: pip install reportlab"
                )
            except Exception as exc:
                logger.error("PDF generation failed: %s", exc)

        return saved


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _esc(text: str) -> str:
    """Escape text for safe inclusion in reportlab paragraphs."""
    if not text:
        return ""
    text = str(text)
    text = text.replace("&", "&amp;")
    text = text.replace("<", "&lt;")
    text = text.replace(">", "&gt;")
    return text


# ---------------------------------------------------------------------------
# Factory function
# ---------------------------------------------------------------------------

def create_report_generator(
    milvus_host: str | None = None,
    milvus_port: int | None = None,
    output_base_dir: str | None = None,
    **kwargs: Any,
) -> ClinicalReportGenerator:
    """
    Factory function to create a ClinicalReportGenerator with settings from HCLSSettings.

    Args:
        milvus_host: Milvus host (falls back to HCLS_MILVUS_HOST).
        milvus_port: Milvus port (falls back to HCLS_MILVUS_PORT).
        output_base_dir: Drug discovery output directory.
        **kwargs: Additional arguments passed to ClinicalReportGenerator.

    Returns:
        Configured ClinicalReportGenerator instance.
    """
    try:
        from hcls_common.config import get_settings

        settings = get_settings()
        milvus_host = milvus_host or settings.milvus_host
        milvus_port = milvus_port or settings.milvus_port
        embedding_model = kwargs.pop("embedding_model", settings.embedding_model)
    except Exception:
        milvus_host = milvus_host or "localhost"
        milvus_port = milvus_port or 19530
        embedding_model = kwargs.pop("embedding_model", "BAAI/bge-small-en-v1.5")

    return ClinicalReportGenerator(
        milvus_host=milvus_host,
        milvus_port=milvus_port,
        embedding_model=embedding_model,
        output_base_dir=output_base_dir or "outputs",
        **kwargs,
    )
