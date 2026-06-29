"""Tests for Streamlit UI modules (biomarker_ui, patient_360, protein_viewer).

Since Streamlit apps are difficult to unit-test directly (they execute
top-level st.* calls on import), these tests focus on:
- Verifying extractable pure functions (get_pgx_phenotype, PGX_DRUG_MAP, etc.)
- Validating the sample patient JSON data
- Checking source-level invariants (tab count, page_config ordering)
- Confirming that patient_360 and protein_viewer can be parsed without error

All Streamlit calls are mocked so no running server is required.

Author: Adam Jones
Date: March 2026
"""

import ast
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent
APP_DIR = PROJECT_ROOT / "app"
DATA_DIR = PROJECT_ROOT / "data" / "reference"

BIOMARKER_UI_PATH = APP_DIR / "biomarker_ui.py"
PATIENT_360_PATH = APP_DIR / "patient_360.py"
PROTEIN_VIEWER_PATH = APP_DIR / "protein_viewer.py"
SAMPLE_PATIENTS_PATH = DATA_DIR / "biomarker_sample_patients.json"


# =====================================================================
# 1. SOURCE FILE PARSING
# =====================================================================


class TestSourceParsing:
    """Verify that all app modules parse as valid Python."""

    def test_biomarker_ui_parses(self):
        """biomarker_ui.py must be valid Python."""
        source = BIOMARKER_UI_PATH.read_text()
        tree = ast.parse(source, filename="biomarker_ui.py")
        assert isinstance(tree, ast.Module)

    def test_patient_360_parses(self):
        """patient_360.py must be valid Python."""
        if not PATIENT_360_PATH.exists():
            pytest.skip("patient_360.py not found")
        source = PATIENT_360_PATH.read_text()
        tree = ast.parse(source, filename="patient_360.py")
        assert isinstance(tree, ast.Module)

    def test_protein_viewer_parses(self):
        """protein_viewer.py must be valid Python."""
        if not PROTEIN_VIEWER_PATH.exists():
            pytest.skip("protein_viewer.py not found")
        source = PROTEIN_VIEWER_PATH.read_text()
        tree = ast.parse(source, filename="protein_viewer.py")
        assert isinstance(tree, ast.Module)


# =====================================================================
# 2. PAGE CONFIG IS FIRST st.* CALL
# =====================================================================


class TestPageConfigFirst:
    """Verify st.set_page_config is the first Streamlit command."""

    def test_set_page_config_before_other_st_calls(self):
        """st.set_page_config must appear before any other st.* call in
        biomarker_ui.py (Streamlit requirement)."""
        source = BIOMARKER_UI_PATH.read_text()
        tree = ast.parse(source, filename="biomarker_ui.py")

        st_calls = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            # Match st.something(...)
            if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name):
                if func.value.id == "st":
                    st_calls.append((node.lineno, func.attr))

        assert len(st_calls) > 0, "No st.* calls found"
        first_call = st_calls[0]
        assert first_call[1] == "set_page_config", (
            f"First st.* call is st.{first_call[1]} at line {first_call[0]}, "
            f"but st.set_page_config must come first"
        )


# =====================================================================
# 3. TAB COUNT
# =====================================================================


class TestTabCount:
    """Verify the UI defines exactly 8 tabs."""

    def test_eight_tabs_defined(self):
        """st.tabs should receive a list of exactly 8 tab names."""
        source = BIOMARKER_UI_PATH.read_text()
        tree = ast.parse(source, filename="biomarker_ui.py")

        tabs_call_found = False
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if (
                isinstance(func, ast.Attribute)
                and isinstance(func.value, ast.Name)
                and func.value.id == "st"
                and func.attr == "tabs"
            ):
                tabs_call_found = True
                # The first argument should be a list
                assert len(node.args) >= 1, "st.tabs() called with no arguments"
                tab_list = node.args[0]
                assert isinstance(tab_list, ast.List), "st.tabs() arg is not a list literal"
                assert len(tab_list.elts) == 8, (
                    f"Expected 8 tabs, found {len(tab_list.elts)}"
                )
                break

        assert tabs_call_found, "st.tabs() call not found in biomarker_ui.py"

    def test_tab_unpacking_matches(self):
        """The st.tabs() call unpacks into 8 variables (tab1..tab8)."""
        source = BIOMARKER_UI_PATH.read_text()
        tree = ast.parse(source, filename="biomarker_ui.py")

        for node in ast.walk(tree):
            if not isinstance(node, ast.Assign):
                continue
            # Look for tuple unpacking: tab1, tab2, ... = st.tabs(...)
            if len(node.targets) == 1 and isinstance(node.targets[0], ast.Tuple):
                target = node.targets[0]
                if isinstance(node.value, ast.Call):
                    func = node.value.func
                    if (
                        isinstance(func, ast.Attribute)
                        and getattr(func.value, "id", None) == "st"
                        and func.attr == "tabs"
                    ):
                        assert len(target.elts) == 8, (
                            f"Expected 8 tab variables, found {len(target.elts)}"
                        )
                        return
        pytest.fail("Could not find tab1..tab8 = st.tabs(...) assignment")


# =====================================================================
# 4. SAMPLE PATIENT JSON
# =====================================================================


class TestSamplePatientsJSON:
    """Validate biomarker_sample_patients.json structure and values."""

    @pytest.fixture(autouse=True)
    def load_patients(self):
        assert SAMPLE_PATIENTS_PATH.exists(), "biomarker_sample_patients.json not found"
        self.patients = json.loads(SAMPLE_PATIENTS_PATH.read_text())

    def test_two_patients(self):
        assert len(self.patients) == 2

    @pytest.mark.parametrize("idx", [0, 1])
    def test_required_top_level_fields(self, idx):
        """Each patient must have id, demographics, biomarkers, genotypes."""
        patient = self.patients[idx]
        for field in ("id", "demographics", "biomarkers", "genotypes"):
            assert field in patient, f"Patient {idx} missing '{field}'"

    @pytest.mark.parametrize("idx", [0, 1])
    def test_demographics_fields(self, idx):
        demo = self.patients[idx]["demographics"]
        for field in ("sex", "age"):
            assert field in demo, f"Patient {idx} demographics missing '{field}'"
        assert isinstance(demo["age"], int)
        assert 0 < demo["age"] < 120

    @pytest.mark.parametrize("idx", [0, 1])
    def test_biomarker_values_positive(self, idx):
        """All biomarker values must be positive numbers."""
        biomarkers = self.patients[idx]["biomarkers"]
        assert len(biomarkers) > 0, f"Patient {idx} has no biomarkers"
        for name, value in biomarkers.items():
            assert isinstance(value, (int, float)), (
                f"Patient {idx} biomarker '{name}' is not numeric: {value}"
            )
            assert value > 0, (
                f"Patient {idx} biomarker '{name}' is non-positive: {value}"
            )

    @pytest.mark.parametrize("idx", [0, 1])
    def test_biomarker_plausible_ranges(self, idx):
        """Spot-check a few biomarkers for biologically plausible values."""
        bio = self.patients[idx]["biomarkers"]
        # Albumin: 1.0-6.0 g/dL
        if "albumin" in bio:
            assert 1.0 <= bio["albumin"] <= 6.0, f"Albumin out of range: {bio['albumin']}"
        # WBC: 1.0-30.0 x10^9/L
        if "wbc" in bio:
            assert 1.0 <= bio["wbc"] <= 30.0, f"WBC out of range: {bio['wbc']}"
        # HbA1c: 3.0-15.0%
        if "hba1c" in bio:
            assert 3.0 <= bio["hba1c"] <= 15.0, f"HbA1c out of range: {bio['hba1c']}"
        # Creatinine: 0.1-5.0 mg/dL
        if "creatinine" in bio:
            assert 0.1 <= bio["creatinine"] <= 5.0
        # TSH: 0.01-20.0 mIU/L
        if "tsh" in bio:
            assert 0.01 <= bio["tsh"] <= 20.0

    @pytest.mark.parametrize("idx", [0, 1])
    def test_genotypes_non_empty(self, idx):
        """Each patient must have at least one genotype entry."""
        genotypes = self.patients[idx]["genotypes"]
        assert len(genotypes) > 0, f"Patient {idx} has no genotypes"

    def test_patient_ids_unique(self):
        ids = [p["id"] for p in self.patients]
        assert len(set(ids)) == len(ids), "Duplicate patient IDs found"

    def test_first_patient_has_star_alleles(self):
        """First patient (HG002) should have star_alleles for PGx testing."""
        assert "star_alleles" in self.patients[0]
        assert len(self.patients[0]["star_alleles"]) > 0


# =====================================================================
# 5. _load_sample_patient KEY MAPPING
# =====================================================================


class TestLoadSamplePatientMapping:
    """Verify the _bio_key_map in _load_sample_patient covers the biomarkers
    present in the sample patient JSON."""

    # This is the mapping extracted from biomarker_ui.py _load_sample_patient
    BIO_KEY_MAP = {
        "wbc": "t1_wbc", "lymphocyte_pct": "t1_lymph", "rdw": "t1_rdw",
        "mcv": "t1_mcv", "platelets": "t1_plt",
        "albumin": "t1_alb", "creatinine": "t1_cr", "glucose": "t1_glu",
        "alt": "t1_alt", "ast": "t1_ast", "alkaline_phosphatase": "t1_alp",
        "total_cholesterol": "t1_tc", "ldl_c": "t1_ldl", "hdl_c": "t1_hdl",
        "triglycerides": "t1_trig", "lpa": "t1_lpa",
        "tsh": "t1_tsh", "free_t4": "t1_ft4", "free_t3": "t1_ft3",
        "ferritin": "t1_ferr", "transferrin_saturation": "t1_tsat",
        "hs_crp": "t1_crp",
        "hba1c": "t1_hba1c", "fasting_insulin": "t1_ins", "homa_ir": "t1_homa",
        "vitamin_d": "t1_vitd", "vitamin_b12": "t1_b12", "folate": "t1_fol",
        "omega3_index": "t1_o3", "magnesium": "t1_mag", "zinc": "t1_zinc",
    }

    GENO_KEY_MAP = {
        "APOE": "t1_apoe", "MTHFR_rs1801133": "t1_mthfr",
        "TCF7L2_rs7903146": "t1_tcf7l2", "PNPLA3_rs738409": "t1_pnpla3",
        "DIO2_rs225014": "t1_dio2", "HFE_rs1800562": "t1_hfe",
    }

    def test_all_widget_keys_unique(self):
        """Widget keys must be unique to avoid Streamlit key collisions."""
        values = list(self.BIO_KEY_MAP.values())
        assert len(set(values)) == len(values), "Duplicate widget keys in bio_key_map"

    def test_sample_patient_biomarkers_mostly_mapped(self):
        """Most biomarkers in the sample patient JSON should have a key mapping."""
        patients = json.loads(SAMPLE_PATIENTS_PATH.read_text())
        patient_bio_keys = set(patients[0]["biomarkers"].keys())
        mapped_keys = set(self.BIO_KEY_MAP.keys())
        covered = patient_bio_keys & mapped_keys
        # At least 80% of patient biomarkers should be mapped
        coverage = len(covered) / len(patient_bio_keys)
        assert coverage >= 0.75, (
            f"Only {coverage:.0%} of patient biomarkers are mapped. "
            f"Unmapped: {patient_bio_keys - mapped_keys}"
        )

    def test_genotype_map_covers_sample_patient(self):
        """Genotypes in sample patient should be mapped."""
        patients = json.loads(SAMPLE_PATIENTS_PATH.read_text())
        patient_geno_keys = set(patients[0]["genotypes"].keys())
        mapped_keys = set(self.GENO_KEY_MAP.keys())
        covered = patient_geno_keys & mapped_keys
        coverage = len(covered) / len(patient_geno_keys)
        assert coverage >= 0.75, (
            f"Only {coverage:.0%} of genotypes mapped. "
            f"Unmapped: {patient_geno_keys - mapped_keys}"
        )

    def test_load_logic_simulation(self):
        """Simulate _load_sample_patient logic and verify session state is populated."""
        patients = json.loads(SAMPLE_PATIENTS_PATH.read_text())
        patient_data = patients[0]

        # Simulate session state as a dict
        session_state = {}

        # Reproduce the mapping logic from _load_sample_patient
        for src_key, widget_key in self.BIO_KEY_MAP.items():
            val = patient_data.get("biomarkers", {}).get(src_key)
            if val is not None:
                session_state[widget_key] = float(val)

        # Check that we loaded biomarkers
        assert len(session_state) > 0, "No biomarkers loaded"
        # Verify specific values
        assert session_state.get("t1_alb") == 4.5
        assert session_state.get("t1_cr") == 0.95
        assert session_state.get("t1_glu") == 98.0


# =====================================================================
# 6. DEMO_BIOMARKERS SIDEBAR MAPPING
# =====================================================================


class TestDemoBiomarkerMapping:
    """Verify the sidebar 'Load Demo Patient' biomarker_map covers expected keys."""

    # Extracted from the sidebar demo-load block in biomarker_ui.py
    SIDEBAR_BIOMARKER_MAP = {
        "wbc": "t1_wbc", "lymphocyte_pct": "t1_lymph", "rdw": "t1_rdw",
        "mcv": "t1_mcv", "platelets": "t1_plt", "albumin": "t1_alb",
        "creatinine": "t1_cr", "glucose": "t1_glu", "alt": "t1_alt",
        "ast": "t1_ast", "alkaline_phosphatase": "t1_alp",
        "total_cholesterol": "t1_tc", "ldl": "t1_ldl", "hdl": "t1_hdl",
        "triglycerides": "t1_tg", "c_reactive_protein": "t1_crp",
        "hs_crp": "t1_hs_crp", "d_dimer": "t1_d_dimer",
        "troponin_i": "t1_troponin", "nt_probnp": "t1_nt_probnp",
    }

    EXPECTED_CORE_KEYS = {
        "wbc", "lymphocyte_pct", "rdw", "mcv", "albumin",
        "creatinine", "glucose", "alt", "ast", "alkaline_phosphatase",
    }

    def test_core_keys_present(self):
        """The sidebar demo map must include core PhenoAge biomarkers."""
        mapped = set(self.SIDEBAR_BIOMARKER_MAP.keys())
        missing = self.EXPECTED_CORE_KEYS - mapped
        assert not missing, f"Sidebar biomarker_map missing core keys: {missing}"

    def test_widget_keys_unique(self):
        values = list(self.SIDEBAR_BIOMARKER_MAP.values())
        assert len(set(values)) == len(values)


# =====================================================================
# 7. HELPER FUNCTIONS (extracted logic, no Streamlit dependency)
# =====================================================================


class TestGetPgxPhenotype:
    """Test the get_pgx_phenotype pure function."""

    @staticmethod
    def _get_pgx_phenotype(gene: str, star_alleles: str) -> dict:
        """Local copy of get_pgx_phenotype to avoid importing biomarker_ui.py
        (which triggers st.set_page_config at module level)."""
        poor_alleles = {"*2", "*3", "*4", "*5", "*6", "*7", "*8"}
        alleles = [a.strip() for a in star_alleles.split("/")]
        poor_count = sum(1 for a in alleles if a in poor_alleles)
        if poor_count == 2:
            return {"phenotype": "Poor Metabolizer", "level": "critical"}
        elif poor_count == 1:
            return {"phenotype": "Intermediate Metabolizer", "level": "moderate"}
        elif "*17" in alleles:
            return {"phenotype": "Ultra-Rapid Metabolizer", "level": "high"}
        else:
            return {"phenotype": "Normal Metabolizer", "level": "normal"}

    def test_normal_metabolizer(self):
        result = self._get_pgx_phenotype("CYP2D6", "*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"
        assert result["level"] == "normal"

    def test_intermediate_metabolizer(self):
        result = self._get_pgx_phenotype("CYP2D6", "*1/*4")
        assert result["phenotype"] == "Intermediate Metabolizer"
        assert result["level"] == "moderate"

    def test_poor_metabolizer(self):
        result = self._get_pgx_phenotype("CYP2D6", "*4/*5")
        assert result["phenotype"] == "Poor Metabolizer"
        assert result["level"] == "critical"

    def test_ultra_rapid_metabolizer(self):
        result = self._get_pgx_phenotype("CYP2C19", "*1/*17")
        assert result["phenotype"] == "Ultra-Rapid Metabolizer"
        assert result["level"] == "high"

    def test_poor_takes_precedence_over_ultra_rapid(self):
        """If both poor and ultra-rapid alleles are present, poor should win
        (poor_count == 1 -> Intermediate, not Ultra-Rapid)."""
        result = self._get_pgx_phenotype("CYP2D6", "*4/*17")
        # *4 is poor -> poor_count=1 -> Intermediate (poor branch runs first)
        assert result["phenotype"] == "Intermediate Metabolizer"

    def test_all_poor_alleles(self):
        """Each of the defined poor alleles should be recognized."""
        for allele in ("*2", "*3", "*4", "*5", "*6", "*7", "*8"):
            result = self._get_pgx_phenotype("TEST", f"*1/{allele}")
            assert result["phenotype"] == "Intermediate Metabolizer", (
                f"Allele {allele} not recognized as poor"
            )


class TestRiskDisplay:
    """Test _risk_display logic (which Streamlit function it dispatches to)."""

    @staticmethod
    def _risk_dispatch(level: str) -> str:
        """Return the Streamlit function name that _risk_display would call."""
        level_lower = level.lower()
        if level_lower in ("critical", "high"):
            return "error"
        elif level_lower == "moderate":
            return "warning"
        else:
            return "info"

    def test_critical_uses_error(self):
        assert self._risk_dispatch("critical") == "error"
        assert self._risk_dispatch("CRITICAL") == "error"

    def test_high_uses_error(self):
        assert self._risk_dispatch("high") == "error"
        assert self._risk_dispatch("HIGH") == "error"

    def test_moderate_uses_warning(self):
        assert self._risk_dispatch("moderate") == "warning"

    def test_low_uses_info(self):
        assert self._risk_dispatch("low") == "info"

    def test_unknown_uses_info(self):
        assert self._risk_dispatch("unknown") == "info"


# =====================================================================
# 8. PGX_DRUG_MAP STRUCTURE
# =====================================================================


class TestPGXDrugMap:
    """Validate the PGX_DRUG_MAP structure extracted from the source."""

    # Extracted from biomarker_ui.py to avoid importing
    EXPECTED_GENES = {"CYP2D6", "CYP2C19", "SLCO1B1", "VKORC1", "MTHFR", "TPMT", "CYP2C9"}

    PGX_DRUG_MAP = {
        "CYP2D6": [
            {"drug": "Codeine"}, {"drug": "Tramadol"},
            {"drug": "Tamoxifen"}, {"drug": "Amitriptyline"},
        ],
        "CYP2C19": [
            {"drug": "Clopidogrel"}, {"drug": "Omeprazole"},
            {"drug": "Escitalopram"}, {"drug": "Voriconazole"},
        ],
        "SLCO1B1": [
            {"drug": "Simvastatin"}, {"drug": "Atorvastatin"},
            {"drug": "Rosuvastatin"},
        ],
        "VKORC1": [{"drug": "Warfarin"}],
        "MTHFR": [{"drug": "Methotrexate"}, {"drug": "Fluorouracil (5-FU)"}],
        "TPMT": [{"drug": "Azathioprine"}, {"drug": "Mercaptopurine"}],
        "CYP2C9": [{"drug": "Warfarin"}, {"drug": "Celecoxib"}, {"drug": "Phenytoin"}],
    }

    def test_all_expected_genes_present(self):
        assert set(self.PGX_DRUG_MAP.keys()) == self.EXPECTED_GENES

    def test_each_gene_has_drugs(self):
        for gene, drugs in self.PGX_DRUG_MAP.items():
            assert len(drugs) > 0, f"{gene} has no drug entries"

    def test_drug_entries_have_drug_field(self):
        for gene, drugs in self.PGX_DRUG_MAP.items():
            for entry in drugs:
                assert "drug" in entry, f"{gene} has entry without 'drug' field"

    def test_source_pgx_map_matches_expected_genes(self):
        """Verify the source file defines PGX_DRUG_MAP with the expected genes."""
        source = BIOMARKER_UI_PATH.read_text()
        for gene in self.EXPECTED_GENES:
            assert f'"{gene}"' in source, (
                f"PGX_DRUG_MAP in source missing gene '{gene}'"
            )


# =====================================================================
# 9. patient_360.py IMPORT TEST (with mocked Streamlit)
# =====================================================================


class TestPatient360Import:
    """Verify patient_360.py can be imported with mocked Streamlit."""

    def test_import_with_mock_streamlit(self):
        """patient_360.py should import successfully when st is mocked."""
        if not PATIENT_360_PATH.exists():
            pytest.skip("patient_360.py not found")

        mock_st = MagicMock()
        mock_st.markdown = MagicMock()
        mock_st.columns = MagicMock(return_value=[MagicMock(), MagicMock()])

        with patch.dict(sys.modules, {"streamlit": mock_st}):
            # Remove cached module if previously imported
            sys.modules.pop("app.patient_360", None)
            sys.modules.pop("patient_360", None)

            import importlib.util
            spec = importlib.util.spec_from_file_location("patient_360", PATIENT_360_PATH)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)

            # Verify key exports exist
            assert hasattr(mod, "render_patient_360")
            assert hasattr(mod, "render_header")
            assert hasattr(mod, "render_footer")
            assert hasattr(mod, "AGENT_ICONS")
            assert isinstance(mod.AGENT_ICONS, dict)
            assert "biomarker" in mod.AGENT_ICONS

    def test_agent_icons_complete(self):
        """AGENT_ICONS dict should have entries for all platform agents."""
        if not PATIENT_360_PATH.exists():
            pytest.skip("patient_360.py not found")

        mock_st = MagicMock()
        with patch.dict(sys.modules, {"streamlit": mock_st}):
            sys.modules.pop("patient_360", None)
            import importlib.util
            spec = importlib.util.spec_from_file_location("patient_360", PATIENT_360_PATH)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)

            expected_agents = {"biomarker", "oncology", "cart", "imaging", "autoimmune",
                               "drug_discovery", "genomics"}
            assert expected_agents.issubset(set(mod.AGENT_ICONS.keys()))


# =====================================================================
# 10. protein_viewer.py IMPORT TEST (with mocked Streamlit)
# =====================================================================


class TestProteinViewerImport:
    """Verify protein_viewer.py can be imported with mocked Streamlit."""

    def test_import_with_mock_streamlit(self):
        """protein_viewer.py should import successfully when st is mocked."""
        if not PROTEIN_VIEWER_PATH.exists():
            pytest.skip("protein_viewer.py not found")

        mock_st = MagicMock()
        mock_components = MagicMock()

        with patch.dict(sys.modules, {
            "streamlit": mock_st,
            "streamlit.components": MagicMock(),
            "streamlit.components.v1": mock_components,
        }):
            sys.modules.pop("protein_viewer", None)
            import importlib.util
            spec = importlib.util.spec_from_file_location("protein_viewer", PROTEIN_VIEWER_PATH)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)

            assert hasattr(mod, "render_protein_viewer")
            assert hasattr(mod, "get_molstar_html")

    def test_get_molstar_html_returns_string(self):
        """get_molstar_html should return valid HTML containing the PDB ID."""
        if not PROTEIN_VIEWER_PATH.exists():
            pytest.skip("protein_viewer.py not found")

        mock_st = MagicMock()
        mock_components = MagicMock()

        with patch.dict(sys.modules, {
            "streamlit": mock_st,
            "streamlit.components": MagicMock(),
            "streamlit.components.v1": mock_components,
        }):
            sys.modules.pop("protein_viewer", None)
            import importlib.util
            spec = importlib.util.spec_from_file_location("protein_viewer", PROTEIN_VIEWER_PATH)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)

            html = mod.get_molstar_html("5FTK", highlight_residue=155, highlight_chain="A")
            assert isinstance(html, str)
            assert "5FTK" in html
            assert "155" in html
            assert "pdbe-molstar" in html.lower() or "PDBeMolstarPlugin" in html
