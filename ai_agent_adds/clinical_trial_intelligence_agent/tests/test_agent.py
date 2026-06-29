"""Tests for agent creation, workflow detection, and entity detection.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.models import (
    TrialWorkflowType,
    TrialQuery,
    TrialResponse,
    PatientProfile,
    SearchPlan,
)


# ===================================================================
# WORKFLOW DETECTION LOGIC
# ===================================================================

# Keyword-based workflow detection (simplified agent logic)
WORKFLOW_KEYWORDS = {
    TrialWorkflowType.PROTOCOL_DESIGN: [
        "protocol", "design", "study design", "arms", "randomization",
    ],
    TrialWorkflowType.PATIENT_MATCHING: [
        "match", "eligible", "patient matching", "find trial",
        "suitable trial", "enroll",
    ],
    TrialWorkflowType.SITE_SELECTION: [
        "site", "location", "center", "facility", "investigator site",
    ],
    TrialWorkflowType.ELIGIBILITY_OPTIMIZATION: [
        "eligibility", "inclusion", "exclusion", "criteria",
        "broaden", "restrict",
    ],
    TrialWorkflowType.ADAPTIVE_DESIGN: [
        "adaptive", "bayesian", "interim analysis", "dose finding",
        "sample size re-estimation",
    ],
    TrialWorkflowType.SAFETY_SIGNAL: [
        "safety", "adverse event", "toxicity", "side effect",
        "signal detection", "DILI", "QT",
    ],
    TrialWorkflowType.REGULATORY_DOCS: [
        "regulatory", "IND", "NDA", "BLA", "FDA", "EMA",
        "submission", "guidance",
    ],
    TrialWorkflowType.COMPETITIVE_INTEL: [
        "competitor", "competitive", "landscape", "rival",
        "market", "pipeline",
    ],
    TrialWorkflowType.DIVERSITY_ASSESSMENT: [
        "diversity", "representation", "equity", "inclusion",
        "underrepresented",
    ],
    TrialWorkflowType.DECENTRALIZED_PLANNING: [
        "decentralized", "DCT", "remote", "telemedicine",
        "wearable", "eConsent",
    ],
    TrialWorkflowType.ELIGIBILITY_ANALYSIS: [
        "eligibility analysis", "criteria analysis", "criterion",
    ],
    TrialWorkflowType.ENDPOINT_STRATEGY: [
        "endpoint strategy", "primary endpoint", "surrogate endpoint",
        "composite endpoint",
    ],
    TrialWorkflowType.SAFETY_MONITORING: [
        "safety monitoring", "DSMB", "pharmacovigilance",
        "risk management",
    ],
    TrialWorkflowType.REGULATORY_STRATEGY: [
        "regulatory strategy", "regulatory pathway", "breakthrough",
        "accelerated approval",
    ],
    TrialWorkflowType.COMPETITIVE_INTELLIGENCE: [
        "competitive intelligence", "threat assessment",
        "competing trials",
    ],
    TrialWorkflowType.BIOMARKER_STRATEGY: [
        "biomarker strategy", "companion diagnostic", "CDx",
        "predictive biomarker",
    ],
    TrialWorkflowType.RWE_ANALYSIS: [
        "real-world evidence", "RWE", "real-world data", "RWD",
        "retrospective",
    ],
    TrialWorkflowType.RECRUITMENT_OPTIMIZATION: [
        "recruitment optimization", "enrollment", "retention",
        "digital recruitment",
    ],
}


def detect_workflow(question: str) -> TrialWorkflowType:
    """Simple keyword-based workflow detection."""
    question_lower = question.lower()
    scores = {}
    for wtype, keywords in WORKFLOW_KEYWORDS.items():
        score = sum(1 for kw in keywords if kw.lower() in question_lower)
        if score > 0:
            scores[wtype] = score
    if scores:
        return max(scores, key=scores.get)
    return TrialWorkflowType.GENERAL


ENTITY_PATTERNS = {
    "therapeutic_area": {
        "oncology": ["cancer", "tumor", "carcinoma", "sarcoma", "lymphoma", "melanoma", "NSCLC"],
        "cardiology": ["heart", "cardiac", "cardiovascular", "HF", "coronary", "arrhythmia"],
        "neurology": ["alzheimer", "parkinson", "neurological", "brain", "dementia"],
        "immunology": ["autoimmune", "rheumatoid", "lupus", "immune"],
    },
    "drug": {
        "pembrolizumab": ["pembrolizumab", "keytruda"],
        "nivolumab": ["nivolumab", "opdivo"],
        "osimertinib": ["osimertinib", "tagrisso"],
        "dapagliflozin": ["dapagliflozin", "farxiga"],
    },
    "biomarker": {
        "PD-L1": ["PD-L1", "pd-l1", "PDL1"],
        "EGFR": ["EGFR", "egfr"],
        "HER2": ["HER2", "her2", "ERBB2"],
        "BRCA": ["BRCA", "brca1", "brca2"],
    },
}


def detect_entities(question: str) -> dict:
    """Simple keyword-based entity detection."""
    results = {"therapeutic_areas": [], "drugs": [], "biomarkers": []}

    for area, keywords in ENTITY_PATTERNS["therapeutic_area"].items():
        if any(kw.lower() in question.lower() for kw in keywords):
            results["therapeutic_areas"].append(area)

    for drug, keywords in ENTITY_PATTERNS["drug"].items():
        if any(kw.lower() in question.lower() for kw in keywords):
            results["drugs"].append(drug)

    for marker, keywords in ENTITY_PATTERNS["biomarker"].items():
        if any(kw.lower() in question.lower() for kw in keywords):
            results["biomarkers"].append(marker)

    return results


# ===================================================================
# TESTS
# ===================================================================


class TestWorkflowDetection:
    """Test workflow detection from query text."""

    def test_protocol_design(self):
        wtype = detect_workflow("Help me design a Phase III protocol")
        assert wtype == TrialWorkflowType.PROTOCOL_DESIGN

    def test_patient_matching(self):
        wtype = detect_workflow("Match this patient to suitable trials")
        assert wtype == TrialWorkflowType.PATIENT_MATCHING

    def test_site_selection(self):
        wtype = detect_workflow("Recommend sites for this trial")
        assert wtype == TrialWorkflowType.SITE_SELECTION

    def test_eligibility(self):
        wtype = detect_workflow("Optimize eligibility criteria")
        assert wtype == TrialWorkflowType.ELIGIBILITY_OPTIMIZATION

    def test_adaptive_design(self):
        wtype = detect_workflow("Consider Bayesian adaptive design for dose finding")
        assert wtype == TrialWorkflowType.ADAPTIVE_DESIGN

    def test_safety(self):
        wtype = detect_workflow("Assess adverse event safety signals")
        assert wtype == TrialWorkflowType.SAFETY_SIGNAL

    def test_regulatory(self):
        wtype = detect_workflow("Prepare FDA IND submission documents")
        assert wtype == TrialWorkflowType.REGULATORY_DOCS

    def test_competitive(self):
        wtype = detect_workflow("Analyze the competitive landscape")
        assert wtype == TrialWorkflowType.COMPETITIVE_INTEL

    def test_diversity(self):
        wtype = detect_workflow("Assess enrollment diversity and representation")
        assert wtype == TrialWorkflowType.DIVERSITY_ASSESSMENT

    def test_dct(self):
        wtype = detect_workflow("Plan a decentralized trial with telemedicine")
        assert wtype == TrialWorkflowType.DECENTRALIZED_PLANNING

    def test_general_fallback(self):
        wtype = detect_workflow("Tell me about clinical research")
        assert wtype == TrialWorkflowType.GENERAL

    def test_explicit_workflow_type(self):
        """If workflow_type is set explicitly, it should be used."""
        q = TrialQuery(
            question="General question",
            workflow_type=TrialWorkflowType.SAFETY_SIGNAL,
        )
        assert q.workflow_type == TrialWorkflowType.SAFETY_SIGNAL


class TestEntityDetection:
    """Test entity detection from query text."""

    def test_oncology_detection(self):
        entities = detect_entities("Find breast cancer clinical trials")
        assert "oncology" in entities["therapeutic_areas"]

    def test_cardiology_detection(self):
        entities = detect_entities("Heart failure trial outcomes")
        assert "cardiology" in entities["therapeutic_areas"]

    def test_drug_detection(self):
        entities = detect_entities("Trials using pembrolizumab for melanoma")
        assert "pembrolizumab" in entities["drugs"]
        assert "oncology" in entities["therapeutic_areas"]

    def test_biomarker_detection(self):
        entities = detect_entities("EGFR-mutant NSCLC with PD-L1 expression")
        assert "EGFR" in entities["biomarkers"]
        assert "PD-L1" in entities["biomarkers"]

    def test_no_entities(self):
        entities = detect_entities("General clinical research question")
        assert entities["therapeutic_areas"] == []
        assert entities["drugs"] == []
        assert entities["biomarkers"] == []

    def test_multiple_entities(self):
        entities = detect_entities(
            "Compare pembrolizumab and nivolumab in melanoma"
        )
        assert "pembrolizumab" in entities["drugs"]
        assert "nivolumab" in entities["drugs"]
        assert "oncology" in entities["therapeutic_areas"]


class TestAgentCreation:
    """Test agent-like object creation patterns."""

    def test_query_creation(self):
        q = TrialQuery(question="Test query")
        assert q.question == "Test query"

    def test_response_creation(self):
        r = TrialResponse(answer="Test answer", confidence=0.5)
        assert r.answer == "Test answer"

    def test_patient_profile_for_matching(self):
        p = PatientProfile(
            age=62,
            sex="male",
            diagnosis="Non-small cell lung cancer",
            biomarkers=["EGFR T790M", "PD-L1 80%"],
            medications=["osimertinib"],
        )
        assert p.age == 62
        assert len(p.biomarkers) == 2

    def test_search_plan_creation(self):
        plan = SearchPlan(
            question="Find NSCLC trials",
            therapeutic_areas=["oncology"],
            relevant_workflows=[
                TrialWorkflowType.PATIENT_MATCHING,
                TrialWorkflowType.COMPETITIVE_INTEL,
            ],
        )
        assert len(plan.relevant_workflows) == 2


class TestAllWorkflowKeywords:
    """Verify all workflow types have detection keywords."""

    @pytest.mark.parametrize("wtype", [
        wt for wt in TrialWorkflowType if wt != TrialWorkflowType.GENERAL
    ])
    def test_workflow_has_keywords(self, wtype):
        assert wtype in WORKFLOW_KEYWORDS, f"No keywords for {wtype}"
        assert len(WORKFLOW_KEYWORDS[wtype]) >= 2, (
            f"Insufficient keywords for {wtype}"
        )
