"""
Production readiness tests for the Precision Autoimmune Agent.

Tests document processor, knowledge base completeness, settings validation,
collection schemas, and model edge cases.

Author: Adam Jones
Date: March 2026
"""

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


class TestDocumentProcessor:
    """Tests for document_processor.py."""

    def setup_method(self):
        from src.document_processor import DocumentProcessor
        self.proc = DocumentProcessor()

    def test_chunk_text_basic(self):
        text = "First sentence. Second sentence. Third sentence."
        chunks = self.proc.chunk_text(text)
        assert len(chunks) >= 1
        assert "First sentence" in chunks[0]

    def test_chunk_text_respects_max_size(self):
        # Create text larger than max_chunk_size
        sentences = [f"Sentence number {i} with some additional text for padding." for i in range(200)]
        text = " ".join(sentences)
        chunks = self.proc.chunk_text(text)
        assert len(chunks) > 1
        for chunk in chunks:
            # Allow some slack for overlap
            assert len(chunk) < self.proc.max_chunk_size * 1.5

    def test_chunk_text_empty(self):
        chunks = self.proc.chunk_text("")
        assert len(chunks) <= 1

    def test_classify_document_type_lab(self):
        text = "Laboratory Report: CBC results showing WBC 8.5, reference range 4.5-11.0"
        doc_type = self.proc.classify_document_type(text)
        assert doc_type == "lab_report"

    def test_classify_document_type_progress_note(self):
        text = "Progress Note: Chief complaint of joint pain. History of present illness..."
        doc_type = self.proc.classify_document_type(text)
        assert doc_type == "progress_note"

    def test_classify_document_type_genetic(self):
        text = "Genetic Test Report: HLA typing results. Allele DRB1*04:01 detected."
        doc_type = self.proc.classify_document_type(text)
        assert doc_type == "genetic_report"

    def test_classify_document_type_unknown(self):
        text = "Random text with no medical markers xyz abc."
        doc_type = self.proc.classify_document_type(text)
        assert doc_type == "clinical_note"  # default fallback

    def test_detect_specialty_rheumatology(self):
        text = "Rheumatology consultation for rheumatoid arthritis management."
        specialty = self.proc.detect_specialty(text)
        assert specialty == "rheumatology"

    def test_detect_specialty_neurology(self):
        text = "Neurology referral for multiple sclerosis evaluation."
        specialty = self.proc.detect_specialty(text)
        assert specialty == "neurology"

    def test_detect_specialty_unknown(self):
        text = "General checkup without specific markers."
        specialty = self.proc.detect_specialty(text)
        assert specialty == "general"

    def test_extract_date_iso(self):
        text = "Visit date: 2026-03-10 follow-up appointment"
        date = self.proc.extract_date(text)
        assert date == "2026-03-10"

    def test_extract_date_us_format(self):
        text = "Date collected: 03/10/2026"
        date = self.proc.extract_date(text)
        assert "03" in date and "10" in date

    def test_extract_date_none(self):
        text = "No date information here."
        date = self.proc.extract_date(text)
        assert date == ""

    def test_extract_provider(self):
        text = "Seen by Dr. Sarah Mitchell, Rheumatology"
        provider = self.proc.extract_provider(text)
        assert "Sarah Mitchell" in provider

    def test_extract_autoantibodies_positive(self):
        text = "ANA positive 1:320 homogeneous pattern. anti-dsDNA elevated 120 IU/mL."
        found = self.proc.extract_autoantibodies(text)
        names = [f["antibody"] for f in found]
        assert "ANA" in names
        assert "anti-dsDNA" in names

    def test_extract_autoantibodies_with_titer(self):
        text = "ANA positive 1:640 speckled pattern"
        found = self.proc.extract_autoantibodies(text)
        ana = [f for f in found if f["antibody"] == "ANA"][0]
        assert ana["titer"] == "1:640"
        assert ana["positive"] is True

    def test_extract_lab_values(self):
        text = "CRP: 12.5 mg/L (elevated). ESR 45 mm/hr."
        found = self.proc.extract_lab_values(text)
        names = [f["test_name"] for f in found]
        assert "CRP" in names

    def test_extract_patient_id_from_path(self):
        path = Path("/data/demo_data/sarah_mitchell/lab_report.pdf")
        pid = self.proc.extract_patient_id_from_path(path)
        assert pid == "sarah_mitchell"

    def test_extract_patient_id_from_path_no_demo(self):
        path = Path("/data/patients/john_doe/report.pdf")
        pid = self.proc.extract_patient_id_from_path(path)
        assert pid == "john_doe"  # falls back to parent name


class TestKnowledgeBaseCompleteness:
    """Verify knowledge base is production-complete."""

    def test_all_biologics_have_monitoring(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        for therapy in BIOLOGIC_THERAPIES:
            assert "monitoring_requirements" in therapy, \
                f"{therapy['drug_name']} missing monitoring_requirements"
            assert len(therapy["monitoring_requirements"]) >= 3, \
                f"{therapy['drug_name']} has too few monitoring requirements"

    def test_all_biologics_have_contraindications(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        for therapy in BIOLOGIC_THERAPIES:
            assert len(therapy.get("contraindications", [])) >= 1, \
                f"{therapy['drug_name']} missing contraindications"

    def test_all_biologics_have_pgx(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        for therapy in BIOLOGIC_THERAPIES:
            assert len(therapy.get("pgx_considerations", [])) >= 1, \
                f"{therapy['drug_name']} missing pgx_considerations"

    def test_hla_associations_have_pmid(self):
        from src.knowledge import HLA_DISEASE_ASSOCIATIONS
        for allele, assocs in HLA_DISEASE_ASSOCIATIONS.items():
            for a in assocs:
                assert "pmid" in a, f"{allele} association missing PMID"
                assert "odds_ratio" in a, f"{allele} association missing odds_ratio"

    def test_autoantibody_map_has_sensitivity_specificity(self):
        from src.knowledge import AUTOANTIBODY_DISEASE_MAP
        for ab, assocs in AUTOANTIBODY_DISEASE_MAP.items():
            for a in assocs:
                assert "sensitivity" in a, f"{ab} missing sensitivity"
                assert "specificity" in a, f"{ab} missing specificity"
                assert 0 <= a["sensitivity"] <= 1, f"{ab} sensitivity out of range"
                assert 0 <= a["specificity"] <= 1, f"{ab} specificity out of range"

    def test_flare_patterns_have_required_fields(self):
        from src.knowledge import FLARE_BIOMARKER_PATTERNS
        for disease, pattern in FLARE_BIOMARKER_PATTERNS.items():
            assert "early_warning_biomarkers" in pattern, f"{disease} missing biomarkers"
            assert "thresholds" in pattern, f"{disease} missing thresholds"
            assert "protective_signals" in pattern, f"{disease} missing protective signals"
            assert len(pattern["early_warning_biomarkers"]) >= 3, \
                f"{disease} has too few biomarkers"

    def test_disease_activity_thresholds_complete(self):
        from src.knowledge import DISEASE_ACTIVITY_THRESHOLDS
        for score_name, info in DISEASE_ACTIVITY_THRESHOLDS.items():
            thresholds = info["thresholds"]
            assert "remission" in thresholds
            assert "low" in thresholds
            assert "moderate" in thresholds
            assert "high" in thresholds
            assert "components" in info
            assert "reference" in info

    def test_autoantibody_names_match_knowledge_base(self):
        """Verify document_processor AUTOANTIBODY_NAMES covers knowledge base."""
        from src.document_processor import AUTOANTIBODY_NAMES
        names_lower = {n.lower() for n in AUTOANTIBODY_NAMES}
        # Check key antibodies from the knowledge base are present
        key_antibodies = ["ANA", "anti-dsDNA", "anti-Smith", "RF", "anti-CCP"]
        for ab in key_antibodies:
            assert ab.lower() in names_lower, f"{ab} not in AUTOANTIBODY_NAMES"


class TestSettingsValidation:
    """Test settings configuration."""

    def test_settings_instantiates(self):
        from config.settings import AutoimmuneSettings
        s = AutoimmuneSettings()
        assert s.EMBEDDING_DIM == 384

    def test_weight_sum_approximately_one(self):
        from config.settings import settings
        total = (
            settings.WEIGHT_CLINICAL_DOCUMENTS + settings.WEIGHT_PATIENT_LABS +
            settings.WEIGHT_AUTOANTIBODY_PANELS + settings.WEIGHT_HLA_ASSOCIATIONS +
            settings.WEIGHT_DISEASE_CRITERIA + settings.WEIGHT_DISEASE_ACTIVITY +
            settings.WEIGHT_FLARE_PATTERNS + settings.WEIGHT_BIOLOGIC_THERAPIES +
            settings.WEIGHT_PGX_RULES + settings.WEIGHT_CLINICAL_TRIALS +
            settings.WEIGHT_LITERATURE + settings.WEIGHT_PATIENT_TIMELINES +
            settings.WEIGHT_CROSS_DISEASE + settings.WEIGHT_GENOMIC_EVIDENCE
        )
        assert abs(total - 1.0) < 0.05, f"Weights sum to {total}, expected ~1.0"

    def test_collection_config_has_all_collections(self):
        from config.settings import settings
        config = settings.collection_config
        assert len(config) == 14
        for name, info in config.items():
            assert "weight" in info
            assert "label" in info
            assert "name" in info

    def test_all_collection_names(self):
        from config.settings import settings
        names = settings.all_collection_names
        assert len(names) == 14
        assert "autoimmune_clinical_documents" in names

    def test_paths(self):
        from config.settings import settings
        assert settings.DATA_DIR.name == "data"
        assert settings.CACHE_DIR.name == "cache"
        assert settings.REFERENCE_DIR.name == "reference"

    def test_ports_distinct(self):
        from config.settings import settings
        assert settings.STREAMLIT_PORT != settings.API_PORT


class TestModelEdgeCases:
    """Test model validation edge cases."""

    def test_patient_profile_lowercase_sex(self):
        from src.models import AutoimmunePatientProfile
        p = AutoimmunePatientProfile(
            patient_id="test", age=30, sex="f",
            biomarkers={"CRP": 5.0},
        )
        assert p.sex == "f"

    def test_patient_profile_uppercase_sex(self):
        from src.models import AutoimmunePatientProfile
        p = AutoimmunePatientProfile(
            patient_id="test", age=30, sex="M",
            biomarkers={"CRP": 5.0},
        )
        assert p.sex == "M"

    def test_patient_profile_requires_data(self):
        from src.models import AutoimmunePatientProfile
        with pytest.raises(Exception):  # ValueError from model_validator
            AutoimmunePatientProfile(
                patient_id="test", age=30, sex="M",
            )

    def test_flare_risk_score_bounds(self):
        from src.models import AutoimmuneDisease, DiseaseActivityLevel, FlarePredictor, FlareRisk
        with pytest.raises(Exception):
            FlarePredictor(
                disease=AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
                current_activity=DiseaseActivityLevel.HIGH,
                predicted_risk=FlareRisk.HIGH,
                risk_score=1.5,  # > 1.0, should fail
            )

    def test_autoantibody_panel_positive_filtering(self):
        from src.models import AutoantibodyPanel, AutoantibodyResult
        panel = AutoantibodyPanel(
            patient_id="test",
            collection_date="2026-03-10",
            results=[
                AutoantibodyResult(antibody="ANA", value=320, positive=True),
                AutoantibodyResult(antibody="RF", value=5, positive=False),
                AutoantibodyResult(antibody="anti-CCP", value=200, positive=True),
            ],
        )
        assert panel.positive_count == 2
        assert "ANA" in panel.positive_antibodies
        assert "RF" not in panel.positive_antibodies

    def test_hla_profile_all_alleles(self):
        from src.models import HLAProfile
        profile = HLAProfile(
            hla_a=["A*02:01"],
            hla_b=["B*27:05"],
            hla_drb1=["DRB1*04:01"],
        )
        all_alleles = profile.all_alleles
        assert len(all_alleles) == 3
        assert "B*27:05" in all_alleles

    def test_disease_enum_values(self):
        from src.models import AutoimmuneDisease
        assert len(AutoimmuneDisease) == 13
        assert AutoimmuneDisease.RHEUMATOID_ARTHRITIS.value == "rheumatoid_arthritis"
        assert AutoimmuneDisease.CELIAC_DISEASE.value == "celiac_disease"


class TestExpandedKnowledgeBase:
    """Tests for the expanded knowledge base (v2.0)."""

    def test_knowledge_version_2(self):
        from src.knowledge import KNOWLEDGE_VERSION
        assert KNOWLEDGE_VERSION["version"] == "2.0.0"
        assert "stats" in KNOWLEDGE_VERSION

    def test_biologic_count_expanded(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        assert len(BIOLOGIC_THERAPIES) >= 20, f"Expected 20+ biologics, got {len(BIOLOGIC_THERAPIES)}"

    def test_biologic_drug_classes_diverse(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        classes = {t["drug_class"] for t in BIOLOGIC_THERAPIES}
        assert len(classes) >= 8, f"Expected 8+ drug classes, got {len(classes)}: {classes}"
        assert any("TNF" in c for c in classes)
        assert any("JAK" in c for c in classes)
        assert any("IL-17" in c for c in classes)

    def test_multiple_tnf_inhibitors(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        tnf_drugs = [t for t in BIOLOGIC_THERAPIES if "TNF" in t["drug_class"]]
        assert len(tnf_drugs) >= 3, "Should have multiple TNF inhibitors"

    def test_ms_biologics_present(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        ms_drugs = [t for t in BIOLOGIC_THERAPIES
                     if "multiple_sclerosis" in t.get("indicated_diseases", [])]
        assert len(ms_drugs) >= 2, "MS should have biologics (natalizumab, ocrelizumab)"

    def test_ibd_biologics_present(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        ibd_drugs = [t for t in BIOLOGIC_THERAPIES
                      if "inflammatory_bowel_disease" in t.get("indicated_diseases", [])]
        assert len(ibd_drugs) >= 4, f"IBD should have 4+ biologics, got {len(ibd_drugs)}"

    def test_disease_activity_scores_expanded(self):
        from src.knowledge import DISEASE_ACTIVITY_THRESHOLDS
        assert len(DISEASE_ACTIVITY_THRESHOLDS) >= 15, \
            f"Expected 15+ scoring systems, got {len(DISEASE_ACTIVITY_THRESHOLDS)}"

    def test_pasi_score_exists(self):
        from src.knowledge import DISEASE_ACTIVITY_THRESHOLDS
        assert "PASI" in DISEASE_ACTIVITY_THRESHOLDS
        assert DISEASE_ACTIVITY_THRESHOLDS["PASI"]["disease"] == "psoriasis"

    def test_edss_score_exists(self):
        from src.knowledge import DISEASE_ACTIVITY_THRESHOLDS
        assert "EDSS" in DISEASE_ACTIVITY_THRESHOLDS
        assert DISEASE_ACTIVITY_THRESHOLDS["EDSS"]["disease"] == "multiple_sclerosis"

    def test_essdai_score_exists(self):
        from src.knowledge import DISEASE_ACTIVITY_THRESHOLDS
        assert "ESSDAI" in DISEASE_ACTIVITY_THRESHOLDS
        assert DISEASE_ACTIVITY_THRESHOLDS["ESSDAI"]["disease"] == "sjogrens_syndrome"

    def test_mrss_score_exists(self):
        from src.knowledge import DISEASE_ACTIVITY_THRESHOLDS
        assert "mRSS" in DISEASE_ACTIVITY_THRESHOLDS
        assert DISEASE_ACTIVITY_THRESHOLDS["mRSS"]["disease"] == "systemic_sclerosis"

    def test_all_13_diseases_have_flare_patterns(self):
        from src.knowledge import FLARE_BIOMARKER_PATTERNS
        assert len(FLARE_BIOMARKER_PATTERNS) == 13
        expected = [
            "rheumatoid_arthritis", "systemic_lupus_erythematosus",
            "inflammatory_bowel_disease", "ankylosing_spondylitis",
            "psoriasis", "sjogrens_syndrome", "systemic_sclerosis",
            "multiple_sclerosis", "type_1_diabetes", "myasthenia_gravis",
            "celiac_disease", "graves_disease", "hashimoto_thyroiditis",
        ]
        for disease in expected:
            assert disease in FLARE_BIOMARKER_PATTERNS, f"{disease} missing flare pattern"

    def test_lab_test_patterns_expanded(self):
        from src.document_processor import LAB_TEST_PATTERNS
        assert len(LAB_TEST_PATTERNS) >= 40, \
            f"Expected 40+ lab patterns, got {len(LAB_TEST_PATTERNS)}"
        # Verify key new additions
        assert "ALT" in LAB_TEST_PATTERNS
        assert "ferritin" in LAB_TEST_PATTERNS
        assert "HbA1c" in LAB_TEST_PATTERNS
        assert "IgG" in LAB_TEST_PATTERNS
        assert "free_T4" in LAB_TEST_PATTERNS
        assert "CK" in LAB_TEST_PATTERNS
        assert "NT_proBNP" in LAB_TEST_PATTERNS

    def test_classification_criteria_expanded(self):
        from src.diagnostic_engine import CLASSIFICATION_CRITERIA
        assert len(CLASSIFICATION_CRITERIA) >= 8, \
            f"Expected 8+ classification criteria, got {len(CLASSIFICATION_CRITERIA)}"

    def test_overlap_syndromes_expanded(self):
        from src.diagnostic_engine import OVERLAP_SYNDROMES
        assert len(OVERLAP_SYNDROMES) >= 7, \
            f"Expected 7+ overlap syndromes, got {len(OVERLAP_SYNDROMES)}"
        assert "t1d_celiac_overlap" in OVERLAP_SYNDROMES
        assert "lupus_aps_overlap" in OVERLAP_SYNDROMES

    def test_every_biologic_has_all_fields(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        required = ["drug_name", "drug_class", "mechanism", "indicated_diseases",
                     "pgx_considerations", "contraindications", "monitoring_requirements"]
        for t in BIOLOGIC_THERAPIES:
            for field in required:
                assert field in t, f"{t['drug_name']} missing {field}"
            assert len(t["indicated_diseases"]) >= 1, f"{t['drug_name']} has no indications"

    def test_certolizumab_pregnancy_safe(self):
        """Verify PGx note about no Fc region = pregnancy-safe."""
        from src.knowledge import BIOLOGIC_THERAPIES
        cert = [t for t in BIOLOGIC_THERAPIES if t["drug_name"] == "Certolizumab pegol"]
        if cert:
            pgx_text = " ".join(cert[0]["pgx_considerations"])
            assert "pregnancy" in pgx_text.lower() or "placental" in pgx_text.lower()

    def test_il23_inhibitors_present(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        il23 = [t for t in BIOLOGIC_THERAPIES if "IL-23" in t["drug_class"]]
        assert len(il23) >= 1, "Should have IL-23 specific inhibitors"

    def test_tyk2_inhibitor_present(self):
        from src.knowledge import BIOLOGIC_THERAPIES
        tyk2 = [t for t in BIOLOGIC_THERAPIES if "TYK2" in t["drug_class"]]
        assert len(tyk2) >= 1, "Should have TYK2 inhibitor (deucravacitinib)"
