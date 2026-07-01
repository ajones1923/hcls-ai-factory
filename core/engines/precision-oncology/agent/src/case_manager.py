"""
Oncology Intelligence Agent - Case Manager
VCF parsing, case lifecycle, and Molecular Tumor Board (MTB) packet generation.

This module is the core oncology-specific logic for managing patient cases,
parsing genomic variants from VCF files, and generating structured MTB packets
for clinical decision support.

Author: Adam Jones
Date: February 2026
License: Apache 2.0
"""

import re
import uuid
import logging
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from src.models import CaseSnapshot, MTBPacket
from src.knowledge import ACTIONABLE_TARGETS, get_target_context, classify_variant_actionability

logger = logging.getLogger(__name__)

# Regex for validating filter values (prevents Milvus injection)
_SAFE_FILTER_RE = re.compile(r"^[A-Za-z0-9 _.\-/]+$")


class OncologyCaseManager:
    """Manages oncology case lifecycle: creation, retrieval, and MTB packet generation.

    Responsibilities:
      - Parse VCF text or accept pre-parsed variant lists
      - Create and persist CaseSnapshot documents in the onco_cases collection
      - Generate structured MTB packets with actionable variants, evidence,
        therapy rankings, trial matches, and open questions
    """

    COLLECTION_NAME = "onco_cases"

    def __init__(self, collection_manager, embedder, knowledge, rag_engine):
        """Initialize the case manager.

        Args:
            collection_manager: Milvus collection manager for vector storage.
            embedder: Embedding model wrapper (e.g., BGE-small-en-v1.5).
            knowledge: Knowledge module providing ACTIONABLE_TARGETS and helpers.
            rag_engine: RAG engine for evidence retrieval and citation generation.
        """
        self.collection_manager = collection_manager
        self.embedder = embedder
        self.knowledge = knowledge
        self.rag_engine = rag_engine

    # ------------------------------------------------------------------
    # Case creation
    # ------------------------------------------------------------------

    def create_case(
        self,
        patient_id: str,
        cancer_type: str,
        stage: str,
        vcf_content_or_variants: Any,
        biomarkers: Optional[Dict[str, Any]] = None,
        prior_therapies: Optional[List[str]] = None,
    ) -> CaseSnapshot:
        """Create a new oncology case from patient data and genomic variants.

        Args:
            patient_id: Unique patient identifier.
            cancer_type: Cancer type string (e.g. "NSCLC", "CRC", "Melanoma").
            stage: Clinical stage (e.g. "IIIB", "IV").
            vcf_content_or_variants: Either raw VCF text (str) or a pre-parsed
                list of variant dicts with keys: gene, variant, chrom, pos, ref, alt.
            biomarkers: Optional dict of biomarker results, e.g.
                {"MSI": "MSI-H", "TMB": 14.2, "PD-L1_TPS": 80}.
            prior_therapies: Optional list of previously administered drug names.

        Returns:
            A fully populated CaseSnapshot.
        """
        biomarkers = biomarkers or {}
        prior_therapies = prior_therapies or []

        # Parse variants from VCF text or accept pre-parsed list
        if isinstance(vcf_content_or_variants, str):
            variants = self._parse_vcf_text(vcf_content_or_variants)
            logger.info("Parsed %d PASS variants from VCF text", len(variants))
        elif isinstance(vcf_content_or_variants, list):
            variants = vcf_content_or_variants
        else:
            raise ValueError(
                "vcf_content_or_variants must be a VCF string or list of variant dicts"
            )

        # Classify actionability for each variant
        for v in variants:
            gene = v.get("gene", "")
            variant_str = v.get("variant", v.get("hgvs", ""))
            v["actionability"] = self._classify_variant_actionability(gene, variant_str)

        case_id = str(uuid.uuid4())
        now = datetime.now(timezone.utc).isoformat()

        # Build a text summary for the case
        variant_genes = [v.get("gene", "unknown") for v in variants[:20]] if variants else []
        text_summary = (
            f"Patient {patient_id} with {cancer_type} stage {stage}. "
            f"Variants: {', '.join(variant_genes) if variant_genes else 'none'}. "
            f"Biomarkers: {biomarkers}. "
            f"Prior therapies: {', '.join(prior_therapies) if prior_therapies else 'none'}."
        )

        snapshot = CaseSnapshot(
            case_id=case_id,
            patient_id=patient_id,
            cancer_type=cancer_type,
            stage=stage,
            variants=variants,
            biomarkers=biomarkers,
            prior_therapies=prior_therapies,
            text_summary=text_summary,
            created_at=now,
            updated_at=now,
        )

        # Embed and store in onco_cases collection
        self._store_case(snapshot)
        logger.info("Created case %s for patient %s (%s %s)", case_id, patient_id, cancer_type, stage)
        return snapshot

    # ------------------------------------------------------------------
    # Case retrieval
    # ------------------------------------------------------------------

    def get_case(self, case_id: str) -> Optional[CaseSnapshot]:
        """Retrieve a case from the onco_cases collection by case_id.

        Args:
            case_id: The UUID of the case to retrieve.

        Returns:
            CaseSnapshot if found, otherwise None.
        """
        try:
            safe_case_id = case_id.strip()
            if not _SAFE_FILTER_RE.match(safe_case_id):
                logger.warning("Rejected unsafe case_id filter value: %r", case_id)
                return None
            results = self.collection_manager.query(
                collection_name=self.COLLECTION_NAME,
                filter_expr=f'case_id == "{safe_case_id}"',
                output_fields=["case_id", "patient_id", "cancer_type", "stage",
                               "variants", "biomarkers", "prior_therapies",
                               "created_at", "updated_at"],
                limit=1,
            )
            if not results:
                logger.warning("Case %s not found", case_id)
                return None

            rec = results[0]
            return CaseSnapshot(
                case_id=rec["case_id"],
                patient_id=rec["patient_id"],
                cancer_type=rec["cancer_type"],
                stage=rec["stage"],
                variants=rec.get("variants", []),
                biomarkers=rec.get("biomarkers", {}),
                prior_therapies=rec.get("prior_therapies", []),
                text_summary=rec.get("text_summary", ""),
                created_at=rec.get("created_at", ""),
                updated_at=rec.get("updated_at", ""),
            )
        except Exception as exc:
            logger.error("Failed to retrieve case %s: %s", case_id, exc)
            return None

    # ------------------------------------------------------------------
    # MTB packet generation
    # ------------------------------------------------------------------

    def generate_mtb_packet(self, case_id_or_snapshot) -> MTBPacket:
        """Generate a Molecular Tumor Board packet for clinical review.

        The MTB packet contains:
          - variant_table: all variants with actionability classification
          - evidence_table: RAG-retrieved evidence and citations per actionable variant
          - therapy_ranking: therapies ranked by evidence level with resistance flags
          - trial_matches: matching clinical trials from onco_trials
          - open_questions: VUS, missing biomarkers, uncertain evidence items

        Args:
            case_id_or_snapshot: Either a case_id string or a CaseSnapshot object.

        Returns:
            A fully populated MTBPacket.
        """
        if isinstance(case_id_or_snapshot, str):
            snapshot = self.get_case(case_id_or_snapshot)
            if snapshot is None:
                raise ValueError(f"Case {case_id_or_snapshot} not found")
        else:
            snapshot = case_id_or_snapshot

        # 1. Extract actionable variants
        actionable_variants = [
            v for v in snapshot.variants
            if v.get("actionability", "VUS") != "VUS"
        ]
        logger.info(
            "Case %s: %d total variants, %d actionable",
            snapshot.case_id, len(snapshot.variants), len(actionable_variants),
        )

        # 2. Build variant table
        variant_table = self._build_variant_table(snapshot.variants, snapshot.cancer_type)

        # 3. Build evidence table via RAG retrieval
        evidence_table = self._build_evidence_table(actionable_variants, snapshot.cancer_type)

        # 4. Build therapy ranking
        therapy_ranking = self._build_therapy_ranking(
            actionable_variants, snapshot.biomarkers,
            snapshot.cancer_type, snapshot.prior_therapies,
        )

        # 5. Build trial matches
        trial_matches = self._build_trial_matches(snapshot)

        # 6. Identify open questions
        open_questions = self._build_open_questions(snapshot)

        packet = MTBPacket(
            case_id=snapshot.case_id,
            patient_id=snapshot.patient_id,
            cancer_type=snapshot.cancer_type,
            stage=snapshot.stage,
            variant_table=variant_table,
            evidence_table=evidence_table,
            therapy_ranking=therapy_ranking,
            trial_matches=trial_matches,
            open_questions=open_questions,
            generated_at=datetime.now(timezone.utc).isoformat(),
        )
        logger.info("Generated MTB packet for case %s", snapshot.case_id)
        return packet

    # ------------------------------------------------------------------
    # VCF parsing
    # ------------------------------------------------------------------

    def _parse_vcf_text(self, vcf_text: str) -> List[Dict]:
        """Parse a VCF text string and extract PASS-only variants.

        Delegates to src.utils.vcf_parser for robust VCF parsing that
        supports SnpEff (ANN=), VEP (CSQ=), and GENE=/GENEINFO= formats.

        Args:
            vcf_text: Raw VCF file content as a string.

        Returns:
            List of variant dicts with keys: chrom, pos, ref, alt, gene,
            consequence, filter, variant (human-readable summary).
        """
        from src.utils.vcf_parser import (
            parse_vcf_text,
            filter_pass_variants,
            extract_gene_from_info,
            extract_consequence_from_info,
        )

        raw_variants = parse_vcf_text(vcf_text)
        pass_variants = filter_pass_variants(raw_variants)

        results = []
        for v in pass_variants:
            info = v.get("info", "")
            gene = extract_gene_from_info(info) if info else v.get("gene", "")
            consequence = extract_consequence_from_info(info) if info else v.get("consequence", "")
            chrom = v.get("chrom", "")
            pos = v.get("pos", "")
            ref_allele = v.get("ref", "")
            alt_allele = v.get("alt", "")

            variant_str = (
                f"{gene} {chrom}:{pos} {ref_allele}>{alt_allele}"
                if gene else f"{chrom}:{pos} {ref_allele}>{alt_allele}"
            )

            results.append({
                "chrom": chrom,
                "pos": int(pos) if isinstance(pos, str) and pos.isdigit() else pos,
                "ref": ref_allele,
                "alt": alt_allele,
                "gene": gene,
                "consequence": consequence,
                "filter": v.get("filter", ""),
                "variant": variant_str,
            })

        return results

    # ------------------------------------------------------------------
    # Actionability classification
    # ------------------------------------------------------------------

    def _classify_variant_actionability(self, gene: str, variant: str) -> str:
        """Classify a variant's actionability against known targets.

        Delegates to the shared classify_variant_actionability function
        in src.knowledge to keep classification logic in one place.

        Args:
            gene: Gene symbol (e.g. "EGFR", "BRAF").
            variant: Variant description (e.g. "L858R", "V600E").

        Returns:
            Evidence level string ("A", "B", "C", "D") or "VUS" if not recognized.
        """
        return classify_variant_actionability(gene, variant)

    # ------------------------------------------------------------------
    # Internal helpers for MTB packet building
    # ------------------------------------------------------------------

    def _store_case(self, snapshot: CaseSnapshot) -> None:
        """Embed and store a case snapshot in the onco_cases collection."""
        summary_text = (
            f"Patient {snapshot.patient_id} with {snapshot.cancer_type} stage {snapshot.stage}. "
            f"Variants: {', '.join(v.get('gene', 'unknown') for v in snapshot.variants[:20])}. "
            f"Biomarkers: {snapshot.biomarkers}. "
            f"Prior therapies: {', '.join(snapshot.prior_therapies) if snapshot.prior_therapies else 'none'}."
        )
        embedding = self.embedder.embed(summary_text)

        # Serialize complex types to strings for Milvus VARCHAR fields
        variants_str = ", ".join(
            f"{v.get('gene', '')} {v.get('variant', '')}" for v in snapshot.variants[:20]
        )
        biomarkers_str = ", ".join(
            f"{k}: {v}" for k, v in snapshot.biomarkers.items()
        ) if isinstance(snapshot.biomarkers, dict) else str(snapshot.biomarkers)
        therapies_str = ", ".join(snapshot.prior_therapies) if snapshot.prior_therapies else ""

        self.collection_manager.insert(
            collection_name=self.COLLECTION_NAME,
            data={
                "id": str(snapshot.case_id)[:100],
                "patient_id": str(snapshot.patient_id)[:100],
                "cancer_type": str(snapshot.cancer_type)[:50],
                "stage": str(snapshot.stage or "")[:20],
                "variants": variants_str[:1000],
                "biomarkers": biomarkers_str[:1000],
                "prior_therapies": therapies_str[:500],
                "embedding": embedding,
                "text_summary": summary_text[:3000],
            },
        )

    def _build_variant_table(self, variants: List[Dict], cancer_type: str) -> List[Dict]:
        """Build the variant table section of the MTB packet."""
        table = []
        for v in variants:
            gene = v.get("gene", "")
            variant_str = v.get("variant", "")
            actionability = v.get("actionability", "VUS")

            # Look up drugs if actionable
            drugs = []
            if actionability != "VUS" and gene.upper() in ACTIONABLE_TARGETS:
                target_info = ACTIONABLE_TARGETS[gene.upper()]
                drugs = target_info.get("drugs", [])

            table.append({
                "gene": gene,
                "variant": variant_str,
                "consequence": v.get("consequence", ""),
                "type": v.get("consequence", "unknown"),
                "evidence_level": actionability,
                "drugs": drugs,
                "chrom": v.get("chrom", ""),
                "pos": v.get("pos", ""),
            })
        return table

    def _build_evidence_table(self, actionable_variants: List[Dict], cancer_type: str) -> List[Dict]:
        """RAG-retrieve evidence for each actionable variant and collect citations."""
        evidence_table = []
        for v in actionable_variants:
            gene = v.get("gene", "")
            variant_str = v.get("variant", "")
            query_text = f"{gene} {variant_str} {cancer_type} targeted therapy clinical evidence"

            try:
                rag_results = self.rag_engine.retrieve(
                    query=query_text,
                    collection_names=["onco_literature", "onco_therapies"],
                    top_k=5,
                )
                citations = [
                    {
                        "source": r.get("source", ""),
                        "text": r.get("text", ""),
                        "score": r.get("score", 0.0),
                    }
                    for r in rag_results
                ]
            except Exception as exc:
                logger.warning("RAG retrieval failed for %s %s: %s", gene, variant_str, exc)
                citations = []

            evidence_table.append({
                "gene": gene,
                "variant": variant_str,
                "evidence_level": v.get("actionability", "VUS"),
                "citations": citations,
                "target_context": get_target_context(gene),
            })
        return evidence_table

    def _build_therapy_ranking(
        self,
        actionable_variants: List[Dict],
        biomarkers: Dict,
        cancer_type: str,
        prior_therapies: List[str],
    ) -> List[Dict]:
        """Rank therapies by evidence level, check resistance, flag contraindications."""
        therapies = []
        seen_drugs = set()

        for v in actionable_variants:
            gene = v.get("gene", "").upper()
            if gene in ACTIONABLE_TARGETS:
                for drug_info in ACTIONABLE_TARGETS[gene].get("drugs", []):
                    drug_name = drug_info if isinstance(drug_info, str) else drug_info.get("name", "")
                    if drug_name and drug_name not in seen_drugs:
                        seen_drugs.add(drug_name)
                        resistance_flag = drug_name.lower() in [t.lower() for t in prior_therapies]
                        therapies.append({
                            "drug_name": drug_name,
                            "evidence_level": v.get("actionability", "C"),
                            "target_gene": gene,
                            "resistance_flag": resistance_flag,
                            "contraindication_flag": False,
                        })

        # Sort by evidence level (A > B > C > D)
        level_order = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "VUS": 5}
        therapies.sort(key=lambda t: level_order.get(t.get("evidence_level", "E"), 4))

        # Assign ranks
        for i, t in enumerate(therapies):
            t["rank"] = i + 1

        return therapies

    def _build_trial_matches(self, snapshot: CaseSnapshot) -> List[Dict]:
        """Search onco_trials for matching cancer_type + biomarkers."""
        query_text = (
            f"{snapshot.cancer_type} stage {snapshot.stage} clinical trial "
            f"biomarkers: {snapshot.biomarkers}"
        )
        try:
            embedding = self.embedder.embed(query_text)
            results = self.collection_manager.search(
                collection_name="onco_trials",
                query_vector=embedding,
                top_k=10,
                output_fields=["trial_id", "title", "phase", "status", "criteria", "text"],
            )
            return [
                {
                    "trial_id": r.get("trial_id", ""),
                    "title": r.get("title", ""),
                    "phase": r.get("phase", ""),
                    "status": r.get("status", ""),
                    "match_score": r.get("score", 0.0),
                }
                for r in results
            ]
        except Exception as exc:
            logger.warning("Trial match search failed: %s", exc)
            return []

    def _build_open_questions(self, snapshot: CaseSnapshot) -> List[str]:
        """Identify VUS, missing biomarkers, and uncertain evidence items."""
        questions = []

        # Flag VUS variants
        vus_variants = [v for v in snapshot.variants if v.get("actionability") == "VUS"]
        if vus_variants:
            vus_genes = list({v.get("gene", "unknown") for v in vus_variants})
            questions.append(
                f"Variants of uncertain significance in: {', '.join(vus_genes[:10])}. "
                "Consider functional assays or updated databases."
            )

        # Check for missing key biomarkers
        expected_biomarkers = {"MSI", "TMB", "PD-L1_TPS"}
        present = set(snapshot.biomarkers.keys())
        missing = expected_biomarkers - present
        if missing:
            questions.append(f"Missing biomarker data: {', '.join(sorted(missing))}. Recommend testing.")

        # Flag if no actionable variants found
        actionable = [v for v in snapshot.variants if v.get("actionability", "VUS") != "VUS"]
        if not actionable:
            questions.append(
                "No actionable variants identified. Consider expanded panel testing, "
                "liquid biopsy, or whole-genome sequencing."
            )

        return questions
