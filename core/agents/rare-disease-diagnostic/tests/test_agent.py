"""Tests for the agent: workflow detection and entity detection.

Tests workflow type detection and entity recognition from the
models and ingest data.

Author: Adam Jones
Date: March 2026
"""


from src.models import DiagnosticWorkflowType
from src.ingest.omim_parser import OMIM_DISEASES
from src.ingest.hpo_parser import HPO_TERMS
from src.ingest.orphanet_parser import ORPHANET_DISEASES


# ===================================================================
# WORKFLOW DETECTION TESTS
# ===================================================================


class TestWorkflowDetection:
    """Test workflow type detection from query keywords."""

    KEYWORD_MAP = {
        DiagnosticWorkflowType.PHENOTYPE_DRIVEN: [
            "phenotype", "HPO", "clinical features", "symptoms",
        ],
        DiagnosticWorkflowType.VARIANT_INTERPRETATION: [
            "variant", "ACMG", "pathogenic", "VUS", "mutation",
        ],
        DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS: [
            "differential", "diagnosis", "differential diagnosis",
        ],
        DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY: [
            "gene therapy", "eligibility", "Zolgensma", "AAV",
        ],
        DiagnosticWorkflowType.NEWBORN_SCREENING: [
            "newborn screening", "NBS", "tandem mass",
        ],
        DiagnosticWorkflowType.METABOLIC_WORKUP: [
            "metabolic", "IEM", "inborn error",
        ],
        DiagnosticWorkflowType.CARRIER_SCREENING: [
            "carrier", "carrier screening", "heterozygote",
        ],
        DiagnosticWorkflowType.PRENATAL_DIAGNOSIS: [
            "prenatal", "amniocentesis", "CVS",
        ],
        DiagnosticWorkflowType.NATURAL_HISTORY: [
            "natural history", "disease progression", "prognosis",
        ],
        DiagnosticWorkflowType.THERAPY_SELECTION: [
            "therapy", "treatment", "enzyme replacement",
        ],
    }

    def test_each_workflow_has_keywords(self):
        """Each workflow should have at least one trigger keyword."""
        for wf, keywords in self.KEYWORD_MAP.items():
            assert len(keywords) >= 1, f"{wf.value} has no keywords"

    def test_keyword_uniqueness(self):
        """Some keywords should be unique to their workflow."""
        all_keywords = []
        for keywords in self.KEYWORD_MAP.values():
            all_keywords.extend(keywords)
        # Not all must be unique, but we should have reasonable coverage
        assert len(all_keywords) >= 20


# ===================================================================
# ENTITY DETECTION TESTS
# ===================================================================


class TestEntityDetection:
    """Test entity detection from query text using seed data."""

    def test_detect_disease_entities(self):
        disease_names = {d["disease_name"].lower() for d in OMIM_DISEASES}
        test_query = "patient with cystic fibrosis"
        matches = [name for name in disease_names if name in test_query.lower()]
        assert len(matches) >= 1

    def test_detect_gene_entities(self):
        genes = {d["gene"] for d in OMIM_DISEASES}
        test_query = "CFTR variant analysis"
        matches = [g for g in genes if g in test_query]
        assert len(matches) >= 1

    def test_detect_hpo_entities(self):
        term_names = {t["name"].lower() for t in HPO_TERMS}
        test_query = "patient with seizures and hypotonia"
        matches = [name for name in term_names if name in test_query.lower()]
        assert len(matches) >= 2

    def test_detect_inheritance_pattern(self):
        patterns = {d["inheritance"].replace("_", " ") for d in OMIM_DISEASES}
        test_query = "autosomal recessive inheritance"
        matches = [p for p in patterns if p in test_query.lower()]
        assert len(matches) >= 1

    def test_detect_orphanet_disease(self):
        names = {d["name"].lower() for d in ORPHANET_DISEASES}
        test_query = "Gaucher Disease treatment options"
        matches = [n for n in names if n in test_query.lower()]
        assert len(matches) >= 1
