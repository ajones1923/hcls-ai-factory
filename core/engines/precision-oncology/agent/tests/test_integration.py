"""
Integration tests for the Oncology Intelligence Agent.
====================================================
Exercises the full agent pipeline without external dependencies (no Milvus,
no LLM API). Uses realistic oncology patient profiles with actual gene names
and validates cross-module consistency.

These tests verify:
  - Full analysis pipeline with realistic patient data
  - Search planning with real gene/cancer type combinations
  - Evidence evaluation across multiple evidence profiles
  - Export functions (Markdown, JSON, FHIR R4)
  - Cross-module consistency (agent -> export round-trip)

Author: Adam Jones
Date: March 2026
"""

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

_AGENT_ROOT = Path(__file__).resolve().parents[1]
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))

from src.agent import OncoIntelligenceAgent
from src.export import export_fhir_r4, export_json, export_markdown
from src.models import (
    AgentResponse,
    CancerType,
    CrossCollectionResult,
    EvidenceLevel,
    MTBPacket,
    OncologyVariant,
    SearchHit,
    VariantType,
)


# ═══════════════════════════════════════════════════════════════════════════
# Realistic Patient Profiles
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def nsclc_egfr_patient():
    """Stage IV NSCLC patient with EGFR L858R and TP53 co-mutation."""
    return {
        "patient_id": "INTEG-NSCLC-001",
        "cancer_type": "nsclc",
        "stage": "IV",
        "sample_id": "S-NSCLC-001",
        "summary": (
            "65-year-old never-smoker female with stage IV lung adenocarcinoma. "
            "NGS panel reveals EGFR L858R (VAF 35%) and TP53 R175H (VAF 42%). "
            "PD-L1 TPS 80%. No prior systemic therapy."
        ),
        "variants": [
            {
                "gene": "EGFR",
                "variant_name": "L858R",
                "variant_type": "missense",
                "vaf": 0.35,
                "consequence": "missense_variant",
                "tier": "I",
            },
            {
                "gene": "TP53",
                "variant_name": "R175H",
                "variant_type": "missense",
                "vaf": 0.42,
                "consequence": "missense_variant",
                "tier": "III",
            },
        ],
        "biomarkers": {
            "tmb": 8.2,
            "msi": "MSS",
            "pdl1": 80,
        },
        "therapies": [
            {
                "name": "osimertinib",
                "targets": ["EGFR"],
                "evidence_level": "A",
                "line_of_therapy": "1L",
                "notes": "FLAURA trial: PFS 18.9 mo vs 10.2 mo for 1st-gen TKIs",
            },
            {
                "name": "osimertinib + chemotherapy",
                "targets": ["EGFR"],
                "evidence_level": "A",
                "line_of_therapy": "1L",
                "notes": "FLAURA2 trial: PFS 25.5 mo",
            },
        ],
        "clinical_trials": [
            {
                "nct_id": "NCT04487080",
                "title": "MARIPOSA: Amivantamab + Lazertinib vs Osimertinib",
                "phase": "Phase 3",
                "status": "Active",
                "match_rationale": "EGFR exon 21 L858R eligible",
            },
        ],
        "resistance_mechanisms": [
            {
                "mechanism": "C797S",
                "drug": "osimertinib",
                "description": "Cysteine-to-serine gatekeeper mutation in EGFR exon 20.",
            },
            {
                "mechanism": "MET amplification",
                "drug": "osimertinib",
                "description": "Bypass pathway activation via MET amplification.",
            },
        ],
        "open_questions": [
            "TP53 R175H co-mutation may confer worse prognosis on EGFR TKI therapy.",
            "Consider serial ctDNA monitoring for early resistance detection.",
        ],
    }


@pytest.fixture
def brca1_breast_patient():
    """Triple-negative breast cancer patient with BRCA1 germline mutation."""
    return {
        "patient_id": "INTEG-BREAST-001",
        "cancer_type": "breast",
        "stage": "IIIA",
        "sample_id": "S-BREAST-001",
        "summary": (
            "42-year-old female with triple-negative breast cancer (TNBC). "
            "Germline BRCA1 5382insC pathogenic variant detected. "
            "Somatic TP53 Y220C (VAF 55%). PD-L1 CPS 15."
        ),
        "variants": [
            {
                "gene": "BRCA1",
                "variant_name": "5382insC",
                "variant_type": "frameshift",
                "vaf": 0.50,
                "consequence": "frameshift_variant",
                "tier": "I",
            },
            {
                "gene": "TP53",
                "variant_name": "Y220C",
                "variant_type": "missense",
                "vaf": 0.55,
                "consequence": "missense_variant",
                "tier": "III",
            },
        ],
        "biomarkers": {
            "tmb": 6.1,
            "msi": "MSS",
            "pdl1": 15,
        },
        "therapies": [
            {
                "name": "olaparib",
                "targets": ["PARP"],
                "evidence_level": "A",
                "line_of_therapy": "2L+",
                "notes": "OlympiAD: BRCA1/2 germline mutation, HER2-negative",
            },
            {
                "name": "talazoparib",
                "targets": ["PARP"],
                "evidence_level": "A",
                "line_of_therapy": "2L+",
                "notes": "EMBRACA trial",
            },
        ],
        "clinical_trials": [
            {
                "nct_id": "NCT04191135",
                "title": "Olaparib + Durvalumab in BRCA-mutated TNBC",
                "phase": "Phase 2",
                "status": "Recruiting",
                "match_rationale": "BRCA1 pathogenic germline variant",
            },
        ],
        "open_questions": [
            "Consider genetic counseling referral for family members.",
            "Evaluate platinum sensitivity given BRCA1 status.",
        ],
    }


@pytest.fixture
def colorectal_kras_patient():
    """Colorectal cancer patient with KRAS G12D and MSI-H."""
    return {
        "patient_id": "INTEG-CRC-001",
        "cancer_type": "colorectal",
        "stage": "IV",
        "sample_id": "S-CRC-001",
        "summary": (
            "58-year-old male with metastatic colorectal adenocarcinoma. "
            "KRAS G12D (VAF 28%), MSI-High, TMB 42 mut/Mb. "
            "Prior FOLFOX + bevacizumab (PD after 6 months)."
        ),
        "variants": [
            {
                "gene": "KRAS",
                "variant_name": "G12D",
                "variant_type": "missense",
                "vaf": 0.28,
                "consequence": "missense_variant",
                "tier": "I",
            },
        ],
        "biomarkers": {
            "tmb": 42,
            "msi": "MSI-H",
            "pdl1": 60,
        },
        "therapies": [
            {
                "name": "pembrolizumab",
                "targets": ["PD-1"],
                "evidence_level": "A",
                "line_of_therapy": "2L",
                "notes": "KEYNOTE-177: MSI-H/dMMR CRC",
            },
            {
                "name": "nivolumab + ipilimumab",
                "targets": ["PD-1", "CTLA-4"],
                "evidence_level": "A",
                "line_of_therapy": "2L",
                "notes": "CheckMate-142: MSI-H CRC",
            },
        ],
        "open_questions": [
            "Lynch syndrome evaluation recommended given MSI-H and age < 60.",
            "KRAS G12D is not targetable with approved agents (sotorasib is G12C-specific).",
        ],
    }


@pytest.fixture
def melanoma_braf_patient():
    """BRAF V600E melanoma patient."""
    return {
        "patient_id": "INTEG-MEL-001",
        "cancer_type": "melanoma",
        "stage": "IIIC",
        "sample_id": "S-MEL-001",
        "summary": (
            "51-year-old male with BRAF V600E-mutant cutaneous melanoma, "
            "unresectable stage IIIC. LDH elevated. No brain metastases."
        ),
        "variants": [
            {
                "gene": "BRAF",
                "variant_name": "V600E",
                "variant_type": "missense",
                "vaf": 0.40,
                "consequence": "missense_variant",
                "tier": "I",
            },
        ],
        "biomarkers": {"tmb": 15, "pdl1": 30},
        "therapies": [
            {
                "name": "dabrafenib + trametinib",
                "targets": ["BRAF", "MEK"],
                "evidence_level": "A",
                "line_of_therapy": "1L",
                "notes": "COMBI-d/v: BRAF V600E/K melanoma",
            },
            {
                "name": "nivolumab + ipilimumab",
                "targets": ["PD-1", "CTLA-4"],
                "evidence_level": "A",
                "line_of_therapy": "1L",
                "notes": "CheckMate-067: 5-year OS 52%",
            },
        ],
    }


# ═══════════════════════════════════════════════════════════════════════════
# Mock RAG engine that returns realistic evidence
# ═══════════════════════════════════════════════════════════════════════════


def _make_mock_rag_engine(gene: str = "EGFR", cancer: str = "NSCLC"):
    """Create a mock RAG engine that returns realistic cross-collection evidence."""
    rag = MagicMock()

    evidence = CrossCollectionResult(
        query=f"{gene} {cancer} treatment",
        hits=[
            SearchHit(
                collection="onco_variants",
                id=f"CIViC:{gene}:001",
                score=0.91,
                text=f"{gene} mutation is a validated therapeutic target in {cancer}.",
                metadata={"gene": gene, "evidence_level": "A"},
            ),
            SearchHit(
                collection="onco_literature",
                id="PMID:33096080",
                score=0.88,
                text=f"Targeted therapy against {gene} in {cancer} improves PFS.",
                metadata={"year": 2024, "cancer_type": cancer.lower()},
            ),
            SearchHit(
                collection="onco_therapies",
                id=f"THERAPY:{gene}_inhibitor",
                score=0.85,
                text=f"FDA-approved {gene} inhibitor for {cancer}.",
                metadata={"category": "targeted", "evidence_level": "A"},
            ),
            SearchHit(
                collection="onco_guidelines",
                id=f"NCCN:{cancer}:2025",
                score=0.82,
                text=f"NCCN recommends {gene} testing for all {cancer} patients.",
                metadata={"org": "NCCN", "year": 2025},
            ),
            SearchHit(
                collection="onco_trials",
                id="NCT04000001",
                score=0.79,
                text=f"Phase 3 trial of next-gen {gene} inhibitor in {cancer}.",
                metadata={"phase": "Phase 3", "status": "Recruiting"},
            ),
        ],
        total_collections_searched=11,
        search_time_ms=85.3,
    )

    rag.cross_collection_search.return_value = [evidence]

    def _synthesize(question, evidence, plan):
        return AgentResponse(
            question=question,
            answer=f"Based on the evidence, {gene} is a key therapeutic target in {cancer}.",
            evidence=CrossCollectionResult(
                query=question,
                hits=evidence[0].hits if evidence else [],
                total_collections_searched=11,
                search_time_ms=85.3,
            ),
            knowledge_used=[f"Gene: {gene}", f"Cancer: {cancer}"],
        )

    rag.synthesize.side_effect = _synthesize
    return rag


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Full Agent Pipeline
# ═══════════════════════════════════════════════════════════════════════════


class TestAgentPipelineIntegration:
    """Test the full plan -> search -> evaluate -> synthesize pipeline."""

    def test_egfr_nsclc_full_pipeline(self):
        """Full pipeline for EGFR L858R in NSCLC."""
        rag = _make_mock_rag_engine("EGFR", "NSCLC")
        agent = OncoIntelligenceAgent(rag_engine=rag)

        response = agent.run("What is the best treatment for EGFR L858R in NSCLC?")

        assert isinstance(response, AgentResponse)
        assert "EGFR" in response.answer
        assert response.evidence.hit_count > 0
        assert response.report is not None
        assert "# Oncology Intelligence Report" in response.report

    def test_braf_melanoma_full_pipeline(self):
        """Full pipeline for BRAF V600E melanoma."""
        rag = _make_mock_rag_engine("BRAF", "MELANOMA")
        agent = OncoIntelligenceAgent(rag_engine=rag)

        response = agent.run("BRAF V600E targeted therapy in melanoma")

        assert isinstance(response, AgentResponse)
        assert response.evidence.hit_count >= 3
        assert response.report is not None

    def test_kras_colorectal_pipeline(self):
        """Full pipeline for KRAS in colorectal cancer."""
        rag = _make_mock_rag_engine("KRAS", "COLORECTAL")
        agent = OncoIntelligenceAgent(rag_engine=rag)

        response = agent.run("KRAS G12C treatment options in colorectal cancer")

        assert isinstance(response, AgentResponse)
        assert "KRAS" in response.answer

    def test_brca1_breast_pipeline(self):
        """Full pipeline for BRCA1 in breast cancer."""
        rag = _make_mock_rag_engine("BRCA1", "BREAST")
        agent = OncoIntelligenceAgent(rag_engine=rag)

        response = agent.run("BRCA1 mutation PARP inhibitor therapy in breast cancer")

        assert isinstance(response, AgentResponse)
        assert response.report is not None

    def test_comparative_query_pipeline(self):
        """Full pipeline for a comparative query."""
        rag = _make_mock_rag_engine("EGFR", "NSCLC")
        agent = OncoIntelligenceAgent(rag_engine=rag)

        response = agent.run("Compare osimertinib vs erlotinib for EGFR NSCLC")

        assert isinstance(response, AgentResponse)
        # Should have used comparative strategy
        assert response.plan is not None

    def test_multi_gene_query_pipeline(self):
        """Full pipeline with multiple genes in query."""
        rag = _make_mock_rag_engine("EGFR", "NSCLC")
        agent = OncoIntelligenceAgent(rag_engine=rag)

        response = agent.run("EGFR and ALK testing guidelines in NSCLC")

        assert isinstance(response, AgentResponse)
        assert response.report is not None


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Search Plan -> Evidence Evaluation Consistency
# ═══════════════════════════════════════════════════════════════════════════


class TestPlanEvidenceConsistency:
    """Verify that search plans produce consistent evidence evaluation."""

    @pytest.fixture
    def agent(self):
        return OncoIntelligenceAgent(rag_engine=MagicMock())

    @pytest.mark.parametrize("query,expected_genes,expected_cancers", [
        ("EGFR L858R treatment in NSCLC", ["EGFR"], ["NSCLC"]),
        ("BRAF V600E in melanoma", ["BRAF"], ["MELANOMA"]),
        ("KRAS mutations in colorectal cancer", ["KRAS"], ["COLORECTAL"]),
        ("BRCA1 breast cancer PARP inhibitor", ["BRCA1"], ["BREAST"]),
        ("ALK fusion in lung cancer", ["ALK"], ["NSCLC"]),
        ("TP53 mutations in pancreatic cancer", ["TP53"], ["PANCREATIC"]),
        ("PIK3CA breast cancer", ["PIK3CA"], ["BREAST"]),
        ("IDH1 mutation in glioblastoma", ["IDH1"], ["GLIOBLASTOMA"]),
    ])
    def test_plan_gene_cancer_extraction(self, agent, query, expected_genes, expected_cancers):
        """Search plan correctly identifies genes and cancer types."""
        plan = agent.search_plan(query)

        for gene in expected_genes:
            assert gene in plan.target_genes, (
                f"Expected gene {gene} in plan for query: {query}"
            )
        for cancer in expected_cancers:
            assert cancer in plan.relevant_cancer_types, (
                f"Expected cancer {cancer} in plan for query: {query}"
            )

    def test_targeted_plan_leads_to_targeted_strategy(self, agent):
        """Specific gene + cancer -> targeted strategy."""
        plan = agent.search_plan("EGFR L858R in NSCLC")
        assert plan.search_strategy == "targeted"

    def test_broad_plan_for_vague_query(self, agent):
        """Vague query -> broad strategy."""
        plan = agent.search_plan("What are the latest advances in oncology?")
        assert plan.search_strategy == "broad"

    def test_sufficient_evidence_evaluation(self, agent):
        """Multiple high-quality hits from multiple collections -> sufficient."""
        evidence_items = []
        for col in ["onco_variants", "onco_literature", "onco_therapies",
                     "onco_guidelines", "onco_trials"]:
            item = MagicMock()
            item.collection = col
            item.score = 0.85
            evidence_items.append(item)
        result = agent.evaluate_evidence(evidence_items)
        assert result == "sufficient"

    def test_insufficient_evidence_for_empty(self, agent):
        """No evidence -> insufficient."""
        assert agent.evaluate_evidence([]) == "insufficient"


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Export Functions (Markdown, JSON, FHIR)
# ═══════════════════════════════════════════════════════════════════════════


class TestExportMarkdownIntegration:
    """Test Markdown export with realistic patient data."""

    def test_nsclc_patient_markdown(self, nsclc_egfr_patient):
        """NSCLC patient export produces complete Markdown report."""
        md = export_markdown(nsclc_egfr_patient)

        assert isinstance(md, str)
        assert "EGFR" in md
        assert "L858R" in md
        assert "osimertinib" in md
        assert "INTEG-NSCLC-001" in md
        assert "Therapy Ranking" in md
        assert "Clinical Trial" in md
        assert "NCT04487080" in md
        assert "C797S" in md
        assert "research use only" in md.lower()

    def test_breast_patient_markdown(self, brca1_breast_patient):
        """BRCA1 breast cancer patient Markdown export."""
        md = export_markdown(brca1_breast_patient)

        assert "BRCA1" in md
        assert "olaparib" in md
        assert "5382insC" in md
        assert "PARP" in md

    def test_colorectal_patient_markdown(self, colorectal_kras_patient):
        """MSI-H colorectal patient Markdown export."""
        md = export_markdown(colorectal_kras_patient)

        assert "KRAS" in md
        assert "MSI-H" in md
        assert "pembrolizumab" in md
        assert "42" in md  # TMB value

    def test_melanoma_patient_markdown(self, melanoma_braf_patient):
        """BRAF melanoma patient Markdown export."""
        md = export_markdown(melanoma_braf_patient)

        assert "BRAF" in md
        assert "V600E" in md
        assert "dabrafenib" in md


class TestExportJSONIntegration:
    """Test JSON export with realistic patient data."""

    def test_nsclc_patient_json(self, nsclc_egfr_patient):
        """NSCLC patient JSON export has all expected fields."""
        result = export_json(nsclc_egfr_patient)

        assert isinstance(result, dict)
        assert result["patient_id"] == "INTEG-NSCLC-001"
        assert result["cancer_type"] == "nsclc"
        assert len(result["variants"]) == 2
        assert len(result["therapy_ranking"]) == 2
        assert result["biomarkers"]["tmb"] == 8.2

    def test_json_serializable(self, nsclc_egfr_patient):
        """JSON export produces serializable data."""
        result = export_json(nsclc_egfr_patient)
        serialized = json.dumps(result)
        reparsed = json.loads(serialized)
        assert reparsed["patient_id"] == "INTEG-NSCLC-001"

    def test_colorectal_patient_json(self, colorectal_kras_patient):
        """Colorectal patient JSON export."""
        result = export_json(colorectal_kras_patient)

        assert result["cancer_type"] == "colorectal"
        assert result["biomarkers"]["msi"] == "MSI-H"
        assert result["biomarkers"]["tmb"] == 42


class TestExportFHIRIntegration:
    """Test FHIR R4 export with realistic patient data."""

    def test_nsclc_fhir_bundle_structure(self, nsclc_egfr_patient):
        """NSCLC patient generates valid FHIR R4 Bundle."""
        bundle = export_fhir_r4(nsclc_egfr_patient, patient_id="INTEG-NSCLC-001")

        assert bundle["resourceType"] == "Bundle"
        assert bundle["type"] == "collection"
        assert "entry" in bundle
        assert len(bundle["entry"]) > 0

    def test_nsclc_fhir_has_patient(self, nsclc_egfr_patient):
        """FHIR bundle contains a Patient resource."""
        bundle = export_fhir_r4(nsclc_egfr_patient, patient_id="INTEG-NSCLC-001")
        patients = [
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "Patient"
        ]
        assert len(patients) == 1
        assert patients[0]["resource"]["identifier"][0]["value"] == "INTEG-NSCLC-001"

    def test_nsclc_fhir_has_diagnostic_report(self, nsclc_egfr_patient):
        """FHIR bundle contains a DiagnosticReport."""
        bundle = export_fhir_r4(nsclc_egfr_patient, patient_id="INTEG-NSCLC-001")
        reports = [
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        ]
        assert len(reports) == 1

    def test_nsclc_fhir_has_variant_observations(self, nsclc_egfr_patient):
        """FHIR bundle contains Observations for variants."""
        bundle = export_fhir_r4(nsclc_egfr_patient, patient_id="INTEG-NSCLC-001")
        observations = [
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "Observation"
        ]
        # At least 2 variant observations + TMB + MSI
        assert len(observations) >= 4

    def test_nsclc_fhir_snomed_cancer_code(self, nsclc_egfr_patient):
        """DiagnosticReport has SNOMED CT cancer type code."""
        bundle = export_fhir_r4(nsclc_egfr_patient, patient_id="INTEG-NSCLC-001")
        report = next(
            e["resource"] for e in bundle["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        )
        assert "conclusionCode" in report
        coding = report["conclusionCode"][0]["coding"][0]
        assert coding["system"] == "http://snomed.info/sct"

    def test_breast_fhir_bundle(self, brca1_breast_patient):
        """Breast cancer patient FHIR bundle."""
        bundle = export_fhir_r4(brca1_breast_patient, patient_id="INTEG-BREAST-001")
        assert bundle["resourceType"] == "Bundle"
        assert len(bundle["entry"]) >= 3  # Patient + report + observations

    def test_colorectal_fhir_bundle(self, colorectal_kras_patient):
        """Colorectal MSI-H patient FHIR bundle."""
        bundle = export_fhir_r4(colorectal_kras_patient, patient_id="INTEG-CRC-001")

        # Should have MSI observation
        msi_obs = [
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "Observation"
            and e["resource"]["code"]["coding"][0].get("code") == "81695-9"
        ]
        assert len(msi_obs) == 1

    def test_fhir_entry_fullurls(self, nsclc_egfr_patient):
        """All FHIR entries have urn:uuid: fullUrls."""
        bundle = export_fhir_r4(nsclc_egfr_patient, patient_id="INTEG-NSCLC-001")
        for entry in bundle["entry"]:
            assert entry["fullUrl"].startswith("urn:uuid:")


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Cross-Module Consistency
# ═══════════════════════════════════════════════════════════════════════════


class TestCrossModuleConsistency:
    """Verify consistency when data flows across agent modules."""

    def test_agent_report_contains_plan_info(self):
        """Agent report includes search plan details."""
        rag = _make_mock_rag_engine("EGFR", "NSCLC")
        agent = OncoIntelligenceAgent(rag_engine=rag)

        response = agent.run("EGFR resistance mechanisms in NSCLC")

        assert response.report is not None
        # Report should reference the strategy
        assert "Strategy" in response.report or "strategy" in response.report

    def test_json_export_matches_markdown_export(self, nsclc_egfr_patient):
        """JSON and Markdown exports produce consistent data for the same patient."""
        md = export_markdown(nsclc_egfr_patient)
        js = export_json(nsclc_egfr_patient)

        # Both should reference the patient ID
        assert "INTEG-NSCLC-001" in md
        assert js["patient_id"] == "INTEG-NSCLC-001"

        # Both should reference EGFR
        assert "EGFR" in md
        egfr_in_variants = any(
            v.get("gene") == "EGFR" for v in js.get("variants", [])
        )
        assert egfr_in_variants

    def test_fhir_export_matches_json_export(self, nsclc_egfr_patient):
        """FHIR bundle variant count matches JSON variant count."""
        js = export_json(nsclc_egfr_patient)
        bundle = export_fhir_r4(nsclc_egfr_patient, patient_id="INTEG-NSCLC-001")

        json_variant_count = len(js.get("variants", []))
        fhir_variant_obs = [
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "Observation"
            and e["resource"]["code"]["coding"][0].get("code") == "69548-6"
        ]
        assert len(fhir_variant_obs) == json_variant_count

    def test_all_export_formats_handle_same_patient(self, melanoma_braf_patient):
        """All three export formats work for the same patient without errors."""
        md = export_markdown(melanoma_braf_patient)
        js = export_json(melanoma_braf_patient)
        fhir = export_fhir_r4(melanoma_braf_patient, patient_id="INTEG-MEL-001")

        assert isinstance(md, str) and len(md) > 0
        assert isinstance(js, dict)
        assert isinstance(fhir, dict)
        assert fhir["resourceType"] == "Bundle"


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Model Construction
# ═══════════════════════════════════════════════════════════════════════════


class TestModelIntegration:
    """Test Pydantic model construction with realistic data."""

    def test_oncology_variant_creation(self):
        """Create OncologyVariant with real gene data."""
        variant = OncologyVariant(
            id="CIViC:EGFR:L858R",
            gene="EGFR",
            variant_name="L858R",
            variant_type=VariantType.SNV,
            cancer_type=CancerType.NSCLC,
            evidence_level=EvidenceLevel.LEVEL_A,
            drugs=["osimertinib", "erlotinib", "gefitinib", "afatinib"],
            text_summary="EGFR L858R is a sensitizing mutation in exon 21.",
            clinical_significance="Pathogenic",
            allele_frequency=0.35,
        )
        assert variant.gene == "EGFR"
        assert variant.evidence_level == EvidenceLevel.LEVEL_A
        assert len(variant.drugs) == 4
        embedding_text = variant.to_embedding_text()
        assert "EGFR" in embedding_text
        assert "L858R" in embedding_text

    def test_mtb_packet_creation(self):
        """Create MTBPacket with realistic tumor board data."""
        packet = MTBPacket(
            case_id="MTB-2026-0042",
            patient_summary="65F, Stage IV NSCLC, EGFR L858R, PD-L1 80%",
            patient_id="PT-NSCLC-042",
            cancer_type="nsclc",
            stage="IV",
            variant_table=[
                {"gene": "EGFR", "variant": "L858R", "tier": "I"},
                {"gene": "TP53", "variant": "R175H", "tier": "III"},
            ],
            therapy_ranking=[
                {"drug": "osimertinib", "evidence": "A", "line": "1L"},
            ],
            trial_matches=[
                {"nct_id": "NCT04487080", "title": "MARIPOSA"},
            ],
            open_questions=[
                "TP53 co-mutation prognostic impact?",
            ],
        )
        assert packet.case_id == "MTB-2026-0042"
        assert len(packet.variant_table) == 2
        assert len(packet.therapy_ranking) == 1

    def test_cross_collection_result_grouping(self):
        """CrossCollectionResult.hits_by_collection groups correctly."""
        result = CrossCollectionResult(
            query="EGFR NSCLC",
            hits=[
                SearchHit(collection="onco_variants", id="1", score=0.9, text="A"),
                SearchHit(collection="onco_variants", id="2", score=0.8, text="B"),
                SearchHit(collection="onco_literature", id="3", score=0.7, text="C"),
                SearchHit(collection="onco_trials", id="4", score=0.6, text="D"),
            ],
            total_collections_searched=3,
        )
        grouped = result.hits_by_collection()
        assert len(grouped["onco_variants"]) == 2
        assert len(grouped["onco_literature"]) == 1
        assert len(grouped["onco_trials"]) == 1
        assert result.hit_count == 4
