"""Biomarker Intelligence Agent -- Streamlit UI v1.0.

Full-featured UI with 8 tabs:
- Biomarker Analysis: full patient analysis pipeline (with sample patient quick-load)
- Biological Age: PhenoAge calculator
- Disease Risk: focused disease trajectory analysis
- PGx Profile: pharmacogenomic drug interaction mapping
- Evidence Explorer: RAG Q&A with collection filtering
- Reports: PDF and FHIR R4 export
- Patient 360: unified cross-agent intelligence dashboard
- Longitudinal: biomarker trend tracking across multiple visits

Port: 8528 (assigned to Biomarker Intelligence Agent)

Usage:
    streamlit run app/biomarker_ui.py --server.port 8528

Author: Adam Jones
Date: March 2026
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import streamlit as st

# Page config MUST be the first Streamlit command
st.set_page_config(
    page_title="Biomarker Intelligence Agent -- HCLS AI Factory",
    page_icon="\U0001fa78",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.warning(
    "**Clinical Decision Support Tool** — This system provides evidence-based guidance "
    "for research and clinical decision support only. All recommendations must be verified "
    "by a qualified healthcare professional. Not FDA-cleared. Not a substitute for professional "
    "clinical judgment."
)

# Add project root to path (must happen before src imports)
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Load API key: prefer env vars, then settings, then rag-chat-pipeline .env file
_api_key = (
    os.environ.get("ANTHROPIC_API_KEY")
    or os.environ.get("BIOMARKER_ANTHROPIC_API_KEY")
)
if not _api_key:
    from config.settings import settings

    if settings.ANTHROPIC_API_KEY:
        _api_key = settings.ANTHROPIC_API_KEY
    else:
        env_path = settings.RAG_PIPELINE_ROOT / ".env"
        if env_path.exists():
            for line in env_path.read_text().splitlines():
                if line.startswith("ANTHROPIC_API_KEY="):
                    _api_key = line.split("=", 1)[1].strip().strip('"')
                    break
if _api_key:
    os.environ["ANTHROPIC_API_KEY"] = _api_key


# =====================================================================
# ENGINE INITIALIZATION
# =====================================================================


@st.cache_resource(ttl=300)
def init_engine():
    """Initialize the Biomarker Intelligence analysis engine (cached across reruns)."""
    try:
        from src.biological_age import BiologicalAgeCalculator
        from src.collections import BiomarkerCollectionManager
        from src.critical_values import CriticalValueEngine
        from src.discordance_detector import DiscordanceDetector
        from src.disease_trajectory import DiseaseTrajectoryAnalyzer
        from src.lab_range_interpreter import LabRangeInterpreter
        from src.pharmacogenomics import PharmacogenomicMapper

        manager = BiomarkerCollectionManager()
        manager.connect()

        bio_age_calc = BiologicalAgeCalculator()
        disease_analyzer = DiseaseTrajectoryAnalyzer()
        critical_engine = CriticalValueEngine()
        discordance_engine = DiscordanceDetector()
        lab_range_engine = LabRangeInterpreter()
        pgx_mapper = PharmacogenomicMapper()

        try:
            from sentence_transformers import SentenceTransformer

            class SimpleEmbedder:
                def __init__(self):
                    self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")

                def embed_text(self, text):
                    return self.model.encode(text).tolist()

                def encode(self, texts):
                    return self.model.encode(texts).tolist()

            embedder = SimpleEmbedder()
        except ImportError:
            embedder = None

        try:
            import anthropic

            from config.settings import settings as cfg

            class SimpleLLMClient:
                def __init__(self):
                    self.client = anthropic.Anthropic()

                def generate(
                    self,
                    prompt,
                    system_prompt="",
                    max_tokens=2048,
                    temperature=0.7,
                ):
                    msg = self.client.messages.create(
                        model=cfg.LLM_MODEL,
                        max_tokens=max_tokens,
                        temperature=temperature,
                        system=system_prompt,
                        messages=[{"role": "user", "content": prompt}],
                    )
                    return msg.content[0].text

                def generate_stream(
                    self,
                    prompt,
                    system_prompt="",
                    max_tokens=2048,
                    temperature=0.7,
                ):
                    with self.client.messages.stream(
                        model=cfg.LLM_MODEL,
                        max_tokens=max_tokens,
                        temperature=temperature,
                        system=system_prompt,
                        messages=[{"role": "user", "content": prompt}],
                    ) as stream:
                        for text in stream.text_stream:
                            yield text

            llm_client = SimpleLLMClient()
        except (ImportError, Exception):
            llm_client = None

        return {
            "manager": manager,
            "bio_age_calc": bio_age_calc,
            "disease_analyzer": disease_analyzer,
            "critical_engine": critical_engine,
            "discordance_engine": discordance_engine,
            "lab_range_engine": lab_range_engine,
            "pgx_mapper": pgx_mapper,
            "embedder": embedder,
            "llm_client": llm_client,
        }
    except Exception as e:
        st.error(f"Failed to initialize engine: {e}")
        return None


engine = init_engine()


# =====================================================================
# PAGE HEADER
# =====================================================================

st.caption("⚕️ Research Use Only — Not for clinical decision-making without healthcare provider review.")

# =====================================================================
# CUSTOM CSS -- NVIDIA Black + Green theme
# =====================================================================

st.markdown(
    """
<style>
    .stApp { background-color: #0a0a0f; }
    .stApp, .stApp p, .stApp span, .stApp li, .stApp td, .stApp th,
    .stApp label, .stApp .stMarkdown {
        color: #ffffff;
    }
    .stApp h1, .stApp h2, .stApp h3, .stApp h4, .stApp h5, .stApp h6 {
        color: #ffffff;
    }
    section[data-testid="stSidebar"] {
        background-color: #12121a;
    }
    section[data-testid="stSidebar"] p,
    section[data-testid="stSidebar"] span,
    section[data-testid="stSidebar"] label {
        color: #e0e0e8;
    }
    section[data-testid="stSidebar"] .stMarkdown h3,
    section[data-testid="stSidebar"] .stMarkdown h2 {
        color: #76B900;
    }
    .stTextInput input, .stTextArea textarea {
        background-color: #1a1a24 !important;
        color: #ffffff !important;
        border: 1px solid #333 !important;
    }
    .stTextInput input:focus, .stTextArea textarea:focus {
        border-color: #76B900 !important;
        box-shadow: 0 0 0 1px #76B900 !important;
    }
    .stSelectbox > div > div {
        background-color: #1a1a24 !important;
        color: #ffffff !important;
    }
    .stNumberInput input {
        background-color: #1a1a24 !important;
        color: #ffffff !important;
    }
    .stButton > button {
        background-color: #76B900 !important;
        color: #000000 !important;
        font-weight: 600;
        border: none;
    }
    .stButton > button:hover {
        background-color: #5a9100 !important;
    }
    .risk-critical { color: #ff4444; font-weight: bold; }
    .risk-high { color: #ff8800; font-weight: bold; }
    .risk-moderate { color: #ffcc00; font-weight: bold; }
    .risk-low { color: #76B900; font-weight: bold; }
    .risk-normal { color: #76B900; font-weight: bold; }
</style>
""",
    unsafe_allow_html=True,
)


# =====================================================================
# SIDEBAR
# =====================================================================

with st.sidebar:
    st.markdown("## Biomarker Intelligence Agent")
    st.caption("HCLS AI Factory | v1.0")

    st.markdown("---")
    st.markdown("### Collection Stats")

    if engine and engine.get("manager"):
        try:
            stats = engine["manager"].get_collection_stats()
            total = sum(stats.values())
            for name, count in stats.items():
                short = name.replace("biomarker_", "").replace("_", " ").title()
                st.text(f"{short}: {count:,}")
            st.markdown(f"**Total: {total:,} vectors**")
        except Exception:
            st.warning("Milvus not connected")
    else:
        st.warning("Engine not initialized")

    st.markdown("---")
    st.markdown("### Service Status")

    status_items = {
        "Milvus": engine is not None and engine.get("manager") is not None,
        "Embedder": engine is not None and engine.get("embedder") is not None,
        "LLM": engine is not None and engine.get("llm_client") is not None,
    }
    for svc, ok in status_items.items():
        icon = "+" if ok else "x"
        st.text(f"[{icon}] {svc}")

    if st.sidebar.button("🔄 Reconnect", key="reconnect_btn"):
        st.cache_resource.clear()
        st.rerun()

    st.markdown("---")
    st.markdown("### About")
    st.caption(
        "Genotype-aware biomarker interpretation, biological age estimation, "
        "disease trajectory detection, and pharmacogenomic profiling. "
        "Part of the HCLS AI Factory precision medicine platform."
    )

    st.sidebar.markdown("---")
    st.sidebar.markdown("### \U0001f3af Demo Mode")
    if st.sidebar.button("Load Demo Patient", key="load_demo"):
        import sys as _sys
        _lib_path = os.environ.get("HCLS_LIB_PATH", "/app/lib")
        _sys.path.insert(0, _lib_path)
        from hcls_common.demo_data import (
            DEMO_BIOMARKERS,
            DEMO_PATIENT_AGE,
            DEMO_PATIENT_ID,
            DEMO_PATIENT_SEX,
        )
        st.session_state["t1_patient_id"] = DEMO_PATIENT_ID
        st.session_state["t1_age"] = DEMO_PATIENT_AGE
        st.session_state["t1_sex"] = DEMO_PATIENT_SEX
        # Map biomarkers to session state keys
        biomarker_map = {
            "wbc": "t1_wbc", "lymphocyte_pct": "t1_lymph", "rdw": "t1_rdw",
            "mcv": "t1_mcv", "platelets": "t1_plt", "albumin": "t1_alb",
            "creatinine": "t1_cr", "glucose": "t1_glu", "alt": "t1_alt",
            "ast": "t1_ast", "alkaline_phosphatase": "t1_alp",
            "total_cholesterol": "t1_tc", "ldl": "t1_ldl", "hdl": "t1_hdl",
            "triglycerides": "t1_tg", "c_reactive_protein": "t1_crp",
            "hs_crp": "t1_hs_crp", "d_dimer": "t1_d_dimer",
            "troponin_i": "t1_troponin", "nt_probnp": "t1_nt_probnp",
        }
        for src, dst in biomarker_map.items():
            if src in DEMO_BIOMARKERS:
                st.session_state[dst] = DEMO_BIOMARKERS[src]
        st.toast("\u2705 Demo patient loaded! Switch to any tab to analyze.", icon="\U0001f3af")
        st.session_state["auto_run_demo"] = True
        st.rerun()


# =====================================================================
# HELPER FUNCTIONS
# =====================================================================


def risk_badge(level: str) -> str:
    """Return an HTML span with risk-level coloring."""
    css_class = f"risk-{level.lower()}"
    return f'<span class="{css_class}">{level.upper()}</span>'


def _risk_display(level: str, prefix: str = "", suffix: str = "") -> None:
    """Render a risk level using safe Streamlit components instead of raw HTML.

    Maps risk levels to appropriate Streamlit callout functions to avoid
    XSS via unsafe_allow_html with user-controlled data.
    """
    label = level.upper()
    text = f"{prefix}: {label} {suffix}".strip() if prefix else f"{label} {suffix}".strip()
    level_lower = level.lower()
    if level_lower in ("critical", "high"):
        st.error(text)
    elif level_lower == "moderate":
        st.warning(text)
    else:
        st.info(text)


def get_pgx_phenotype(gene: str, star_alleles: str) -> dict:
    """Determine metabolizer phenotype from star alleles (simplified)."""
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


# PGx drug mapping (gene -> list of drugs with recommendations)
PGX_DRUG_MAP = {
    "CYP2D6": [
        {"drug": "Codeine", "poor": "AVOID - no analgesic effect", "ultra_rapid": "AVOID - risk of toxicity", "normal": "Standard dosing"},
        {"drug": "Tramadol", "poor": "AVOID - reduced efficacy", "ultra_rapid": "AVOID - risk of respiratory depression", "normal": "Standard dosing"},
        {"drug": "Tamoxifen", "poor": "Consider aromatase inhibitor", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
        {"drug": "Amitriptyline", "poor": "50% dose reduction", "ultra_rapid": "Avoid or increase monitoring", "normal": "Standard dosing"},
    ],
    "CYP2C19": [
        {"drug": "Clopidogrel", "poor": "AVOID - use prasugrel/ticagrelor", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
        {"drug": "Omeprazole", "poor": "50% dose reduction", "ultra_rapid": "Increase dose 2-3x", "normal": "Standard dosing"},
        {"drug": "Escitalopram", "poor": "50% dose reduction", "ultra_rapid": "Consider dose increase", "normal": "Standard dosing"},
        {"drug": "Voriconazole", "poor": "Use alternative antifungal", "ultra_rapid": "Increase dose", "normal": "Standard dosing"},
    ],
    "SLCO1B1": [
        {"drug": "Simvastatin", "poor": "AVOID high doses - myopathy risk", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
        {"drug": "Atorvastatin", "poor": "Use lower dose", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
        {"drug": "Rosuvastatin", "poor": "Use lower starting dose", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
    ],
    "VKORC1": [
        {"drug": "Warfarin", "poor": "Reduce dose 25-50%", "ultra_rapid": "Increase dose", "normal": "Standard dosing"},
    ],
    "MTHFR": [
        {"drug": "Methotrexate", "poor": "Dose reduction; supplement with leucovorin", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
        {"drug": "Fluorouracil (5-FU)", "poor": "Monitor for toxicity", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
    ],
    "TPMT": [
        {"drug": "Azathioprine", "poor": "AVOID or 90% dose reduction", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
        {"drug": "Mercaptopurine", "poor": "AVOID or 90% dose reduction", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
    ],
    "CYP2C9": [
        {"drug": "Warfarin", "poor": "Reduce dose 25-50%", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
        {"drug": "Celecoxib", "poor": "50% dose reduction", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
        {"drug": "Phenytoin", "poor": "Reduce dose; monitor levels", "ultra_rapid": "Standard dosing", "normal": "Standard dosing"},
    ],
}


# =====================================================================
# TABS
# =====================================================================

st.markdown("# Biomarker Intelligence Agent")
st.caption("Genotype-aware biomarker interpretation | HCLS AI Factory")

tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs(
    [
        "\U0001f52c Biomarker Analysis",
        "\U0001f9ec Biological Age",
        "\U0001f4ca Disease Risk",
        "\U0001f48a PGx Profile",
        "\U0001f50d Evidence Explorer",
        "\U0001f4cb Reports",
        "\U0001f310 Patient 360",
        "\U0001f4c8 Longitudinal",
    ]
)


# =====================================================================
# TAB 1: BIOMARKER ANALYSIS
# =====================================================================

with tab1:
    # Auto-run demo analysis if triggered (flag consumed at line ~705 via pop)
    if st.session_state.get("auto_run_demo"):
        st.toast("🚀 Auto-running demo analysis...", icon="⚡")

    st.markdown("## Full Patient Biomarker Analysis")
    st.caption(
        "Enter patient information, biomarker values, and genotype data for "
        "a comprehensive precision medicine analysis."
    )

    # Sample Patient Quick Load
    st.markdown("#### Quick Start: Load Sample Patient")
    sample_col1, sample_col2 = st.columns(2)

    _sample_patients_path = PROJECT_ROOT / "data" / "reference" / "biomarker_sample_patients.json"
    if _sample_patients_path.exists():
        _sample_data = json.loads(_sample_patients_path.read_text())
        def _load_sample_patient(patient_data):
            """Load sample patient biomarkers, genotypes, and star alleles into session state."""
            st.session_state["demo_patient"] = patient_data
            # Map biomarker values to widget session state keys
            _bio_key_map = {
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
            for src_key, widget_key in _bio_key_map.items():
                val = patient_data.get("biomarkers", {}).get(src_key)
                if val is not None:
                    st.session_state[widget_key] = float(val)
            # Load genotypes
            for gt_key, gt_val in patient_data.get("genotypes", {}).items():
                widget_key = f"t1_{gt_key.lower().replace('_rs', '_rs').split('_')[0]}"
                # Map to known widget keys
                _geno_key_map = {
                    "APOE": "t1_apoe", "MTHFR_rs1801133": "t1_mthfr",
                    "TCF7L2_rs7903146": "t1_tcf7l2", "PNPLA3_rs738409": "t1_pnpla3",
                    "DIO2_rs225014": "t1_dio2", "HFE_rs1800562": "t1_hfe",
                }
                if gt_key in _geno_key_map:
                    st.session_state[_geno_key_map[gt_key]] = gt_val

        with sample_col1:
            if st.button("Load Male Patient (HG002)", key="load_male"):
                _load_sample_patient(_sample_data[0])
                st.rerun()
        with sample_col2:
            if st.button("Load Female Patient", key="load_female"):
                _load_sample_patient(_sample_data[1])
                st.rerun()

    st.markdown("---")

    # Pre-fill from demo patient if loaded
    _demo = st.session_state.get("demo_patient")
    _default_pid = _demo["id"] if _demo else "PATIENT-001"
    _default_age = _demo["demographics"]["age"] if _demo else 45
    _default_sex_val = _demo["demographics"]["sex"] if _demo else "M"
    # Map full sex name to M/F abbreviation used by selectbox
    _sex_options = ["M", "F"]
    if _default_sex_val in ("Male", "male"):
        _default_sex_idx = 0
    elif _default_sex_val in ("Female", "female"):
        _default_sex_idx = 1
    else:
        _default_sex_idx = _sex_options.index(_default_sex_val) if _default_sex_val in _sex_options else 0

    patient_id = st.text_input("Patient ID", value=_default_pid, key="t1_patient_id")

    col_info1, col_info2 = st.columns(2)
    with col_info1:
        patient_age = st.number_input("Age", min_value=0, max_value=150, value=_default_age, key="t1_age")
    with col_info2:
        patient_sex = st.selectbox("Sex", _sex_options, index=_default_sex_idx, key="t1_sex")

    # -- Biomarker value inputs --
    biomarkers = {}

    with st.expander("CBC (Complete Blood Count)", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.number_input("WBC (10^3/uL)", min_value=0.0, value=0.0, step=0.1, key="t1_wbc")
            if v > 0:
                biomarkers["wbc"] = v
        with c2:
            v = st.number_input("Lymphocyte %", min_value=0.0, max_value=100.0, value=0.0, step=0.1, key="t1_lymph")
            if v > 0:
                biomarkers["lymphocyte_pct"] = v
        with c3:
            v = st.number_input("RDW %", min_value=0.0, value=0.0, step=0.1, key="t1_rdw")
            if v > 0:
                biomarkers["rdw"] = v
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.number_input("MCV (fL)", min_value=0.0, value=0.0, step=0.1, key="t1_mcv")
            if v > 0:
                biomarkers["mcv"] = v
        with c2:
            v = st.number_input("Platelets (10^3/uL)", min_value=0.0, value=0.0, step=1.0, key="t1_plt")
            if v > 0:
                biomarkers["platelets"] = v

    with st.expander("CMP (Comprehensive Metabolic Panel)", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.number_input("Albumin (g/dL)", min_value=0.0, value=0.0, step=0.1, key="t1_alb")
            if v > 0:
                biomarkers["albumin"] = v
        with c2:
            v = st.number_input("Creatinine (mg/dL)", min_value=0.0, value=0.0, step=0.01, key="t1_cr")
            if v > 0:
                biomarkers["creatinine"] = v
        with c3:
            v = st.number_input("Glucose (mg/dL)", min_value=0.0, value=0.0, step=1.0, key="t1_glu")
            if v > 0:
                biomarkers["glucose"] = v
                biomarkers["fasting_glucose"] = v
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.number_input("ALT (U/L)", min_value=0.0, value=0.0, step=1.0, key="t1_alt")
            if v > 0:
                biomarkers["alt"] = v
        with c2:
            v = st.number_input("AST (U/L)", min_value=0.0, value=0.0, step=1.0, key="t1_ast")
            if v > 0:
                biomarkers["ast"] = v
        with c3:
            v = st.number_input("Alk Phos (U/L)", min_value=0.0, value=0.0, step=1.0, key="t1_alp")
            if v > 0:
                biomarkers["alkaline_phosphatase"] = v

    with st.expander("Lipids", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.number_input("Total Cholesterol (mg/dL)", min_value=0.0, value=0.0, step=1.0, key="t1_tc")
            if v > 0:
                biomarkers["total_cholesterol"] = v
        with c2:
            v = st.number_input("LDL-C (mg/dL)", min_value=0.0, value=0.0, step=1.0, key="t1_ldl")
            if v > 0:
                biomarkers["ldl_c"] = v
        with c3:
            v = st.number_input("HDL-C (mg/dL)", min_value=0.0, value=0.0, step=1.0, key="t1_hdl")
            if v > 0:
                biomarkers["hdl_c"] = v
        c1, c2, _ = st.columns(3)
        with c1:
            v = st.number_input("Triglycerides (mg/dL)", min_value=0.0, value=0.0, step=1.0, key="t1_trig")
            if v > 0:
                biomarkers["triglycerides"] = v
        with c2:
            v = st.number_input("Lp(a) (nmol/L)", min_value=0.0, value=0.0, step=1.0, key="t1_lpa")
            if v > 0:
                biomarkers["lpa"] = v

    with st.expander("Thyroid", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.number_input("TSH (mIU/L)", min_value=0.0, value=0.0, step=0.01, key="t1_tsh")
            if v > 0:
                biomarkers["tsh"] = v
        with c2:
            v = st.number_input("Free T4 (ng/dL)", min_value=0.0, value=0.0, step=0.01, key="t1_ft4")
            if v > 0:
                biomarkers["free_t4"] = v
        with c3:
            v = st.number_input("Free T3 (pg/mL)", min_value=0.0, value=0.0, step=0.01, key="t1_ft3")
            if v > 0:
                biomarkers["free_t3"] = v

    with st.expander("Iron Studies", expanded=False):
        c1, c2 = st.columns(2)
        with c1:
            v = st.number_input("Ferritin (ng/mL)", min_value=0.0, value=0.0, step=1.0, key="t1_ferr")
            if v > 0:
                biomarkers["ferritin"] = v
        with c2:
            v = st.number_input("Transferrin Saturation (%)", min_value=0.0, max_value=100.0, value=0.0, step=0.1, key="t1_tsat")
            if v > 0:
                biomarkers["transferrin_saturation"] = v

    with st.expander("Inflammation", expanded=False):
        v = st.number_input("hs-CRP (mg/L)", min_value=0.0, value=0.0, step=0.01, key="t1_crp")
        if v > 0:
            biomarkers["hs_crp"] = v

    with st.expander("Diabetes", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.number_input("HbA1c (%)", min_value=0.0, value=0.0, step=0.1, key="t1_hba1c")
            if v > 0:
                biomarkers["hba1c"] = v
        with c2:
            v = st.number_input("Fasting Insulin (uIU/mL)", min_value=0.0, value=0.0, step=0.1, key="t1_ins")
            if v > 0:
                biomarkers["fasting_insulin"] = v
        with c3:
            v = st.number_input("HOMA-IR", min_value=0.0, value=0.0, step=0.1, key="t1_homa")
            if v > 0:
                biomarkers["homa_ir"] = v

    with st.expander("Vitamins & Minerals", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.number_input("Vitamin D (ng/mL)", min_value=0.0, value=0.0, step=0.1, key="t1_vitd")
            if v > 0:
                biomarkers["vitamin_d"] = v
        with c2:
            v = st.number_input("Vitamin B12 (pg/mL)", min_value=0.0, value=0.0, step=1.0, key="t1_b12")
            if v > 0:
                biomarkers["vitamin_b12"] = v
        with c3:
            v = st.number_input("Folate (ng/mL)", min_value=0.0, value=0.0, step=0.1, key="t1_fol")
            if v > 0:
                biomarkers["folate"] = v
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.number_input("Omega-3 Index (%)", min_value=0.0, max_value=100.0, value=0.0, step=0.1, key="t1_o3")
            if v > 0:
                biomarkers["omega3_index"] = v
        with c2:
            v = st.number_input("Magnesium (mg/dL)", min_value=0.0, value=0.0, step=0.1, key="t1_mag")
            if v > 0:
                biomarkers["magnesium"] = v
        with c3:
            v = st.number_input("Zinc (mcg/dL)", min_value=0.0, value=0.0, step=1.0, key="t1_zinc")
            if v > 0:
                biomarkers["zinc"] = v

    # -- Genotype inputs --
    genotypes = {}
    with st.expander("Genotype Inputs (rsIDs)", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.text_input("TCF7L2 rs7903146", placeholder="e.g., CT", key="t1_tcf7l2")
            if v:
                genotypes["TCF7L2_rs7903146"] = v
        with c2:
            v = st.text_input("PNPLA3 rs738409", placeholder="e.g., CG", key="t1_pnpla3")
            if v:
                genotypes["PNPLA3_rs738409"] = v
        with c3:
            v = st.text_input("DIO2 rs225014", placeholder="e.g., GA", key="t1_dio2")
            if v:
                genotypes["DIO2_rs225014"] = v
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.text_input("APOE", placeholder="e.g., E3/E4", key="t1_apoe")
            if v:
                genotypes["APOE"] = v
        with c2:
            v = st.text_input("HFE rs1800562", placeholder="e.g., GA", key="t1_hfe")
            if v:
                genotypes["HFE_rs1800562"] = v
        with c3:
            v = st.text_input("MTHFR rs1801133", placeholder="e.g., CT", key="t1_mthfr")
            if v:
                genotypes["MTHFR_rs1801133"] = v
        c1, c2, _ = st.columns(3)
        with c1:
            v = st.text_input("FADS1 rs174546", placeholder="e.g., CT", key="t1_fads1")
            if v:
                genotypes["FADS1_rs174546"] = v
        with c2:
            v = st.text_input("HSD17B13 rs72613567", placeholder="e.g., T/TA", key="t1_hsd")
            if v:
                genotypes["HSD17B13_rs72613567"] = v

    # -- Star allele inputs --
    star_alleles = {}
    with st.expander("Star Allele Inputs (PGx)", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.text_input("CYP2D6", placeholder="e.g., *1/*2", key="t1_cyp2d6")
            if v:
                star_alleles["CYP2D6"] = v
        with c2:
            v = st.text_input("CYP2C19", placeholder="e.g., *1/*1", key="t1_cyp2c19")
            if v:
                star_alleles["CYP2C19"] = v
        with c3:
            v = st.text_input("SLCO1B1", placeholder="e.g., *1/*1", key="t1_slco1b1")
            if v:
                star_alleles["SLCO1B1"] = v
        c1, c2, c3 = st.columns(3)
        with c1:
            v = st.text_input("VKORC1", placeholder="e.g., A/G", key="t1_vkorc1")
            if v:
                star_alleles["VKORC1"] = v
        with c2:
            v = st.text_input("MTHFR (star)", placeholder="e.g., CT", key="t1_mthfr_star")
            if v:
                star_alleles["MTHFR"] = v
        with c3:
            v = st.text_input("TPMT", placeholder="e.g., *1/*1", key="t1_tpmt")
            if v:
                star_alleles["TPMT"] = v

    # -- Run Analysis --
    auto_run = st.session_state.pop("auto_run_demo", False)
    if st.button("Run Full Analysis", type="primary", key="t1_run") or auto_run:
        # Clear stale results from previous runs
        for key in ["analysis_results", "patient_info", "report_markdown"]:
            st.session_state.pop(key, None)

        if not biomarkers:
            st.warning("Please enter at least one biomarker value.")
        elif engine is None:
            st.error("Engine not initialized. Check Milvus connection.")
        else:
            with st.spinner("Running comprehensive analysis..."):
                results = {}

                # Biological age
                try:
                    bio_age = engine["bio_age_calc"].calculate(
                        chronological_age=float(patient_age),
                        biomarkers=biomarkers,
                    )
                    results["biological_age"] = bio_age
                except Exception as e:
                    st.warning(f"Biological age calculation error: {e}")

                # Disease trajectories
                try:
                    trajectories = engine["disease_analyzer"].analyze_all(
                        biomarkers=biomarkers,
                        genotypes=genotypes,
                        age=float(patient_age),
                        sex="male" if patient_sex == "M" else "female",
                    )
                    results["disease_trajectories"] = trajectories
                except Exception as e:
                    st.warning(f"Disease trajectory error: {e}")

                # Critical value checking
                try:
                    critical_alerts = engine["critical_engine"].check(biomarkers)
                    results["critical_alerts"] = critical_alerts
                except Exception as e:
                    st.warning(f"Critical value check error: {e}")

                # Discordance detection
                try:
                    discordance_findings = engine["discordance_engine"].check(biomarkers)
                    results["discordance_findings"] = discordance_findings
                except Exception as e:
                    st.warning(f"Discordance detection error: {e}")

                # Lab range interpretation
                try:
                    lab_sex = "Male" if patient_sex == "M" else "Female"
                    lab_interpretations = engine["lab_range_engine"].interpret(biomarkers, sex=lab_sex)
                    lab_discrepancies = engine["lab_range_engine"].get_discrepancies(biomarkers, sex=lab_sex)
                    results["lab_interpretations"] = lab_interpretations
                    results["lab_discrepancies"] = lab_discrepancies
                except Exception as e:
                    st.warning(f"Lab range interpretation error: {e}")

                # Store results in session state
                st.session_state["analysis_results"] = results
                st.session_state["patient_info"] = {
                    "patient_id": patient_id,
                    "age": patient_age,
                    "sex": patient_sex,
                    "biomarkers": biomarkers,
                    "genotypes": genotypes,
                    "star_alleles": star_alleles,
                }

            # Display results
            if "biological_age" in results:
                ba = results["biological_age"]
                pheno = ba.get("phenoage", ba)
                st.markdown("### Biological Age")
                c1, c2, c3 = st.columns(3)
                with c1:
                    st.metric("Chronological Age", f"{patient_age} years")
                with c2:
                    bio_val = pheno.get("biological_age", ba.get("biological_age", "N/A"))
                    st.metric("Biological Age", f"{bio_val} years")
                with c3:
                    accel = pheno.get("age_acceleration", ba.get("age_acceleration", 0))
                    risk = pheno.get("mortality_risk", ba.get("mortality_risk", "NORMAL"))
                    st.metric("Acceleration", f"{accel:+.1f} years")
                    _risk_display(risk, prefix="Risk")

            if "disease_trajectories" in results:
                st.markdown("### Disease Trajectories")
                for traj in results["disease_trajectories"]:
                    level = traj.get("risk_level", "LOW")
                    name = traj.get("display_name", traj.get("disease", ""))
                    stage = traj.get("stage", "")
                    _risk_display(level, prefix=f"**{name}**", suffix=f"| Stage: {stage}")
                    if traj.get("findings"):
                        for f in traj["findings"]:
                            st.markdown(f"- {f}")
                    if traj.get("recommendations"):
                        with st.expander(f"Recommendations for {name}"):
                            for r in traj["recommendations"]:
                                st.markdown(f"- {r}")

            if results.get("critical_alerts"):
                st.markdown("### Critical Value Alerts")
                for alert in results["critical_alerts"]:
                    if alert.severity.lower() == "critical":
                        st.error(alert.to_alert_string())
                    else:
                        st.warning(alert.to_alert_string())

            if results.get("discordance_findings"):
                st.markdown("### Cross-Biomarker Discordance")
                for finding in results["discordance_findings"]:
                    st.warning(finding.to_alert_string())

            if results.get("lab_discrepancies"):
                st.markdown("### Lab Range Optimization Opportunities")
                for comp in results["lab_discrepancies"]:
                    st.info(comp.to_interpretation())


# =====================================================================
# TAB 2: BIOLOGICAL AGE
# =====================================================================

with tab2:
    st.markdown("## Biological Age Calculator")
    st.caption(
        "Quick PhenoAge estimation using 9 routine blood biomarkers "
        "(Levine et al., Aging 2018)."
    )

    t2_age = st.number_input("Chronological Age", min_value=1, max_value=150, value=45, key="t2_age")

    st.markdown("#### PhenoAge Biomarkers")
    c1, c2, c3 = st.columns(3)
    with c1:
        t2_alb = st.number_input("Albumin (g/dL)", min_value=0.01, value=4.2, step=0.1, key="t2_alb")
        t2_cr = st.number_input("Creatinine (mg/dL)", min_value=0.01, value=0.9, step=0.01, key="t2_cr")
        t2_glu = st.number_input("Glucose (mg/dL)", min_value=0.01, value=95.0, step=1.0, key="t2_glu")
    with c2:
        t2_crp = st.number_input("hs-CRP (mg/L)", min_value=0.01, value=1.5, step=0.1, key="t2_crp")
        t2_lymph = st.number_input("Lymphocyte %", min_value=0.01, max_value=100.0, value=30.0, step=0.1, key="t2_lymph")
        t2_mcv = st.number_input("MCV (fL)", min_value=0.01, value=89.0, step=0.1, key="t2_mcv")
    with c3:
        t2_rdw = st.number_input("RDW %", min_value=0.01, value=13.0, step=0.1, key="t2_rdw")
        t2_alp = st.number_input("Alk Phos (U/L)", min_value=0.01, value=65.0, step=1.0, key="t2_alp")
        t2_wbc = st.number_input("WBC (10^3/uL)", min_value=0.01, value=6.0, step=0.1, key="t2_wbc")

    if st.button("Calculate Biological Age", type="primary", key="t2_run"):
        if engine is None:
            st.error("Engine not initialized.")
        else:
            pheno_biomarkers = {
                "albumin": t2_alb,
                "creatinine": t2_cr,
                "glucose": t2_glu,
                "hs_crp": t2_crp,
                "lymphocyte_pct": t2_lymph,
                "mcv": t2_mcv,
                "rdw": t2_rdw,
                "alkaline_phosphatase": t2_alp,
                "wbc": t2_wbc,
            }

            result = engine["bio_age_calc"].calculate(
                chronological_age=float(t2_age),
                biomarkers=pheno_biomarkers,
            )

            c1, c2, c3 = st.columns(3)
            with c1:
                st.metric("Biological Age", f"{result['biological_age']} years")
            with c2:
                st.metric("Age Acceleration", f"{result['age_acceleration']:+.1f} years")
            with c3:
                _risk_display(result['mortality_risk'], prefix="Risk Level")

            st.markdown("#### Top Aging Drivers")
            driver_data = []
            for d in result.get("phenoage", {}).get("top_aging_drivers", result.get("top_aging_drivers", [])):
                driver_data.append({
                    "Biomarker": d["biomarker"],
                    "Value": d["value"],
                    "Contribution": d["contribution"],
                    "Direction": d["direction"],
                })
            if driver_data:
                st.table(driver_data)


# =====================================================================
# TAB 3: DISEASE RISK
# =====================================================================

with tab3:
    st.markdown("## Disease Risk Analysis")
    st.caption("Focused disease trajectory analysis with genotype-stratified thresholds.")

    t3_age = st.number_input("Patient Age", min_value=1, max_value=120, value=45, key="t3_age_input")

    disease_options = {
        "Type 2 Diabetes": "type2_diabetes",
        "Cardiovascular Disease": "cardiovascular",
        "Liver Disease (NAFLD/Fibrosis)": "liver",
        "Thyroid Dysfunction": "thyroid",
        "Iron Metabolism Disorder": "iron",
        "Nutritional Deficiency": "nutritional",
    }

    selected_disease = st.selectbox(
        "Disease Category",
        list(disease_options.keys()),
        key="t3_disease",
    )
    disease_key = disease_options[selected_disease]

    st.markdown(f"#### Relevant Biomarkers for {selected_disease}")
    t3_biomarkers = {}
    t3_genotypes = {}

    if disease_key == "type2_diabetes":
        c1, c2 = st.columns(2)
        with c1:
            v = st.number_input("HbA1c (%)", min_value=0.0, value=0.0, step=0.1, key="t3_hba1c")
            if v > 0:
                t3_biomarkers["hba1c"] = v
            v = st.number_input("Fasting Glucose (mg/dL)", min_value=0.0, value=0.0, step=1.0, key="t3_glu")
            if v > 0:
                t3_biomarkers["fasting_glucose"] = v
        with c2:
            v = st.number_input("Fasting Insulin (uIU/mL)", min_value=0.0, value=0.0, step=0.1, key="t3_ins")
            if v > 0:
                t3_biomarkers["fasting_insulin"] = v
            gt = st.text_input("TCF7L2 rs7903146", placeholder="e.g., CT", key="t3_tcf7l2")
            if gt:
                t3_genotypes["TCF7L2_rs7903146"] = gt

    elif disease_key == "cardiovascular":
        c1, c2 = st.columns(2)
        with c1:
            v = st.number_input("LDL-C (mg/dL)", min_value=0.0, value=0.0, step=1.0, key="t3_ldl")
            if v > 0:
                t3_biomarkers["ldl_c"] = v
            v = st.number_input("Lp(a) (nmol/L)", min_value=0.0, value=0.0, step=1.0, key="t3_lpa")
            if v > 0:
                t3_biomarkers["lpa"] = v
            v = st.number_input("hs-CRP (mg/L)", min_value=0.0, value=0.0, step=0.1, key="t3_crp")
            if v > 0:
                t3_biomarkers["hs_crp"] = v
        with c2:
            v = st.number_input("ApoB (mg/dL)", min_value=0.0, value=0.0, step=1.0, key="t3_apob")
            if v > 0:
                t3_biomarkers["apob"] = v
            gt = st.text_input("APOE", placeholder="e.g., E3/E4", key="t3_apoe")
            if gt:
                t3_genotypes["APOE"] = gt

    elif disease_key == "liver":
        c1, c2 = st.columns(2)
        with c1:
            v = st.number_input("ALT (U/L)", min_value=0.0, value=0.0, step=1.0, key="t3_alt")
            if v > 0:
                t3_biomarkers["alt"] = v
            v = st.number_input("AST (U/L)", min_value=0.0, value=0.0, step=1.0, key="t3_ast")
            if v > 0:
                t3_biomarkers["ast"] = v
            v = st.number_input("Platelets (10^3/uL)", min_value=0.0, value=0.0, step=1.0, key="t3_plt")
            if v > 0:
                t3_biomarkers["platelets"] = v
        with c2:
            v = st.number_input("Ferritin (ng/mL)", min_value=0.0, value=0.0, step=1.0, key="t3_ferr")
            if v > 0:
                t3_biomarkers["ferritin"] = v
            gt = st.text_input("PNPLA3 rs738409", placeholder="e.g., CG", key="t3_pnpla3")
            if gt:
                t3_genotypes["PNPLA3_rs738409"] = gt

    elif disease_key == "thyroid":
        c1, c2 = st.columns(2)
        with c1:
            v = st.number_input("TSH (mIU/L)", min_value=0.0, value=0.0, step=0.01, key="t3_tsh")
            if v > 0:
                t3_biomarkers["tsh"] = v
            v = st.number_input("Free T4 (ng/dL)", min_value=0.0, value=0.0, step=0.01, key="t3_ft4")
            if v > 0:
                t3_biomarkers["free_t4"] = v
        with c2:
            v = st.number_input("Free T3 (pg/mL)", min_value=0.0, value=0.0, step=0.01, key="t3_ft3")
            if v > 0:
                t3_biomarkers["free_t3"] = v
            gt = st.text_input("DIO2 rs225014", placeholder="e.g., GA", key="t3_dio2")
            if gt:
                t3_genotypes["DIO2_rs225014"] = gt

    elif disease_key == "iron":
        c1, c2 = st.columns(2)
        with c1:
            v = st.number_input("Ferritin (ng/mL)", min_value=0.0, value=0.0, step=1.0, key="t3_ferr_iron")
            if v > 0:
                t3_biomarkers["ferritin"] = v
            v = st.number_input("Transferrin Sat (%)", min_value=0.0, value=0.0, step=0.1, key="t3_tsat")
            if v > 0:
                t3_biomarkers["transferrin_saturation"] = v
        with c2:
            gt = st.text_input("HFE rs1800562 (C282Y)", placeholder="e.g., GA", key="t3_hfe")
            if gt:
                t3_genotypes["HFE_rs1800562"] = gt

    elif disease_key == "nutritional":
        c1, c2 = st.columns(2)
        with c1:
            v = st.number_input("Vitamin D (ng/mL)", min_value=0.0, value=0.0, step=0.1, key="t3_vitd")
            if v > 0:
                t3_biomarkers["vitamin_d"] = v
            v = st.number_input("Vitamin B12 (pg/mL)", min_value=0.0, value=0.0, step=1.0, key="t3_b12")
            if v > 0:
                t3_biomarkers["vitamin_b12"] = v
            v = st.number_input("Folate (ng/mL)", min_value=0.0, value=0.0, step=0.1, key="t3_fol")
            if v > 0:
                t3_biomarkers["folate"] = v
        with c2:
            v = st.number_input("Omega-3 Index (%)", min_value=0.0, value=0.0, step=0.1, key="t3_o3")
            if v > 0:
                t3_biomarkers["omega3_index"] = v
            gt = st.text_input("MTHFR rs1801133", placeholder="e.g., CT", key="t3_mthfr")
            if gt:
                t3_genotypes["MTHFR_rs1801133"] = gt
            gt2 = st.text_input("FADS1 rs174546", placeholder="e.g., CT", key="t3_fads1")
            if gt2:
                t3_genotypes["FADS1_rs174546"] = gt2

    if st.button("Analyze Risk", type="primary", key="t3_run"):
        if not t3_biomarkers:
            st.warning("Please enter at least one biomarker value.")
        elif engine is None:
            st.error("Engine not initialized.")
        else:
            analyzer = engine["disease_analyzer"]
            method_map = {
                "type2_diabetes": analyzer.analyze_type2_diabetes,
                "cardiovascular": analyzer.analyze_cardiovascular,
                "liver": lambda b, g: analyzer.analyze_liver(b, g, age=float(t3_age)),
                "thyroid": analyzer.analyze_thyroid,
                "iron": lambda b, g: analyzer.analyze_iron(b, g, sex="male"),
                "nutritional": analyzer.analyze_nutritional,
            }

            result = method_map[disease_key](t3_biomarkers, t3_genotypes)

            st.markdown(f"### {result['display_name']}")
            c1, c2 = st.columns(2)
            with c1:
                _risk_display(result['risk_level'], prefix="Risk Level")
            with c2:
                st.markdown(f"Stage: **{result['stage']}**")

            if result.get("findings"):
                st.markdown("#### Findings")
                for f in result["findings"]:
                    st.markdown(f"- {f}")

            if result.get("recommendations"):
                st.markdown("#### Recommendations")
                for r in result["recommendations"]:
                    st.markdown(f"- {r}")

            # Store disease results for pathway map
            st.session_state["disease_results"] = [result]
            st.session_state["bm_results"] = t3_biomarkers

    # Interactive Pathway Network Visualization
    if st.session_state.get("disease_results"):
        with st.expander("🗺️ Interactive Pathway Map", expanded=False):
            try:
                import plotly.graph_objects as go

                # Build network from disease trajectory results
                diseases = st.session_state.get("disease_results", [])
                biomarkers_used = st.session_state.get("bm_results", {})

                # Define pathway network nodes
                nodes = []
                node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
                edges_x, edges_y = [], []

                # Center: Patient
                nodes.append("Patient")
                node_x.append(0.5); node_y.append(0.5)  # noqa: E702
                node_text.append(f"Patient {st.session_state.get('t1_patient_id', '')}")
                node_color.append("#76B900"); node_size.append(30)  # noqa: E702

                # Disease nodes in a circle around center
                import math
                disease_names = []
                if isinstance(diseases, list):
                    for d in diseases:
                        if isinstance(d, dict):
                            disease_names.append(d.get("disease", d.get("name", d.get("display_name", "Unknown"))))
                        elif hasattr(d, "disease"):
                            disease_names.append(d.disease)

                risk_colors = {"LOW": "#22c55e", "MODERATE": "#f59e0b", "HIGH": "#ef4444", "CRITICAL": "#dc2626", "NORMAL": "#22c55e"}

                n_diseases = max(len(disease_names), 1)
                for i, disease in enumerate(disease_names):
                    angle = 2 * math.pi * i / n_diseases
                    x = 0.5 + 0.35 * math.cos(angle)
                    y = 0.5 + 0.35 * math.sin(angle)

                    nodes.append(disease)
                    node_x.append(x); node_y.append(y)  # noqa: E702

                    # Get risk level
                    risk = "UNKNOWN"
                    if isinstance(diseases, list) and i < len(diseases):
                        d = diseases[i]
                        if isinstance(d, dict):
                            risk = d.get("risk_level", d.get("stage", "UNKNOWN"))
                        elif hasattr(d, "risk_level"):
                            risk = d.risk_level

                    node_text.append(f"{disease}<br>Risk: {risk}")
                    node_color.append(risk_colors.get(str(risk).upper(), "#6b7280"))
                    node_size.append(22)

                    # Edge from patient to disease
                    edges_x.extend([0.5, x, None])
                    edges_y.extend([0.5, y, None])

                # Biomarker nodes in outer ring
                key_biomarkers = {
                    "albumin": "Albumin", "creatinine": "Creatinine", "glucose": "Glucose",
                    "c_reactive_protein": "CRP", "hba1c": "HbA1c", "alt": "ALT",
                    "tsh": "TSH", "ferritin": "Ferritin", "homocysteine": "Homocysteine",
                }

                bm_vals = {}
                if isinstance(biomarkers_used, dict):
                    bm_vals = biomarkers_used

                bm_i = 0
                for key, label in key_biomarkers.items():
                    val = bm_vals.get(key) or st.session_state.get(f"t1_{key}")
                    if val:
                        angle = 2 * math.pi * bm_i / len(key_biomarkers) + 0.3
                        x = 0.5 + 0.18 * math.cos(angle)
                        y = 0.5 + 0.18 * math.sin(angle)

                        nodes.append(label)
                        node_x.append(x); node_y.append(y)  # noqa: E702
                        node_text.append(f"{label}: {val}")
                        node_color.append("#3b82f6")
                        node_size.append(14)

                        # Edge from biomarker to patient
                        edges_x.extend([x, 0.5, None])
                        edges_y.extend([y, 0.5, None])
                        bm_i += 1

                fig = go.Figure()

                # Add edges
                fig.add_trace(go.Scatter(
                    x=edges_x, y=edges_y, mode='lines',
                    line=dict(width=1, color='#4b5563'),
                    hoverinfo='none',
                ))

                # Add nodes
                fig.add_trace(go.Scatter(
                    x=node_x, y=node_y, mode='markers+text',
                    marker=dict(size=node_size, color=node_color, line=dict(width=1, color='#1f2937')),
                    text=nodes, textposition="top center", textfont=dict(size=10, color='white'),
                    hovertext=node_text, hoverinfo='text',
                ))

                fig.update_layout(
                    showlegend=False,
                    plot_bgcolor='#0a0a0f',
                    paper_bgcolor='#0a0a0f',
                    font=dict(color='white'),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    margin=dict(l=20, r=20, t=30, b=20),
                    height=500,
                    title=dict(text="Biomarker → Disease Pathway Network", font=dict(size=14, color="#76B900")),
                )

                st.plotly_chart(fig, use_container_width=True)
                st.caption("🟢 Green = Low risk | 🟡 Yellow = Moderate | 🔴 Red = High/Critical | 🔵 Blue = Biomarkers")
            except ImportError:
                st.info("Install plotly for interactive pathway visualization: `pip install plotly`")


# =====================================================================
# TAB 4: PGx PROFILE
# =====================================================================

with tab4:
    st.markdown("## Pharmacogenomic Profile")
    st.caption("Star allele interpretation for 7 pharmacogenes + HLA-B*57:01 status.")

    t4_star_alleles = {}
    c1, c2, c3 = st.columns(3)
    with c1:
        v = st.text_input("CYP2D6", placeholder="e.g., *1/*4", key="t4_cyp2d6")
        if v:
            t4_star_alleles["CYP2D6"] = v
        v = st.text_input("CYP2C19", placeholder="e.g., *2/*2", key="t4_cyp2c19")
        if v:
            t4_star_alleles["CYP2C19"] = v
        v = st.text_input("CYP2C9", placeholder="e.g., *1/*1", key="t4_cyp2c9")
        if v:
            t4_star_alleles["CYP2C9"] = v
    t4_genotypes = {}
    with c2:
        v = st.text_input("SLCO1B1 rs4149056", placeholder="e.g., TC", key="t4_slco1b1")
        if v:
            t4_genotypes["SLCO1B1_rs4149056"] = v
        v = st.text_input("VKORC1 rs9923231", placeholder="e.g., AG", key="t4_vkorc1")
        if v:
            t4_genotypes["VKORC1_rs9923231"] = v
        v = st.text_input("MTHFR rs1801133", placeholder="e.g., CT", key="t4_mthfr")
        if v:
            t4_genotypes["MTHFR_rs1801133"] = v
    with c3:
        v = st.text_input("TPMT", placeholder="e.g., *3A/*3A", key="t4_tpmt")
        if v:
            t4_star_alleles["TPMT"] = v
        t4_hla = st.checkbox("HLA-B*57:01 Positive", key="t4_hla")

    if st.button("Map Drug Interactions", type="primary", key="t4_run"):
        if not t4_star_alleles and not t4_genotypes and not t4_hla:
            st.warning("Please enter at least one gene or check HLA-B*57:01.")
        else:
            st.markdown("### Drug Interaction Map")

            table_data = []
            critical_alerts = []
            used_mapper = False

            # Try the proper PharmacogenomicMapper first
            mapper = engine.get("pgx_mapper") if engine else None
            if mapper is not None:
                try:
                    genotypes_for_mapper = dict(t4_genotypes)
                    if t4_hla:
                        genotypes_for_mapper["HLA_B5701"] = "positive"
                    mapper_result = mapper.map_all(
                        star_alleles=t4_star_alleles,
                        genotypes=genotypes_for_mapper,
                    )
                    used_mapper = True

                    for gene_res in mapper_result.get("gene_results", []):
                        gene = gene_res.get("gene", "")
                        phenotype = gene_res.get("phenotype") or "Unknown"
                        star_al = gene_res.get("star_alleles") or gene_res.get("genotype") or ""
                        for drug_rec in gene_res.get("affected_drugs", []):
                            table_data.append({
                                "Gene": gene,
                                "Star Alleles": star_al,
                                "Phenotype": phenotype,
                                "Drug": drug_rec["drug"],
                                "Action": drug_rec.get("action", ""),
                                "Recommendation": drug_rec["recommendation"],
                            })

                    for alert in mapper_result.get("critical_alerts", []):
                        if isinstance(alert, dict):
                            msg = alert.get("message", str(alert))
                            gene = alert.get("gene", "")
                            drug = alert.get("drug", "")
                            critical_alerts.append(f"{gene} + {drug}: {msg}" if gene else msg)
                        else:
                            critical_alerts.append(str(alert))

                except Exception as e:
                    st.warning(f"PharmacogenomicMapper error, falling back to simplified logic: {e}")
                    used_mapper = False
                    table_data = []
                    critical_alerts = []

            # Fallback: simplified local logic
            if not used_mapper:
                for gene, alleles in t4_star_alleles.items():
                    pheno_info = get_pgx_phenotype(gene, alleles)
                    drugs = PGX_DRUG_MAP.get(gene, [])

                    for drug_info in drugs:
                        phenotype = pheno_info["phenotype"]
                        if "Poor" in phenotype:
                            rec = drug_info.get("poor", "Standard dosing")
                        elif "Ultra" in phenotype:
                            rec = drug_info.get("ultra_rapid", "Standard dosing")
                        else:
                            rec = drug_info.get("normal", "Standard dosing")

                        if "AVOID" in rec:
                            critical_alerts.append(f"{gene} {alleles} + {drug_info['drug']}: {rec}")

                        table_data.append({
                            "Gene": gene,
                            "Star Alleles": alleles,
                            "Phenotype": phenotype,
                            "Drug": drug_info["drug"],
                            "Recommendation": rec,
                        })

                if t4_hla:
                    critical_alerts.append("HLA-B*57:01 Positive + Abacavir: CONTRAINDICATED")
                    table_data.append({
                        "Gene": "HLA-B",
                        "Star Alleles": "*57:01",
                        "Phenotype": "Positive",
                        "Drug": "Abacavir",
                        "Recommendation": "CONTRAINDICATED - risk of hypersensitivity reaction",
                    })

            if critical_alerts:
                st.markdown("#### Critical Alerts")
                for alert in critical_alerts:
                    st.error(f"WARNING: {alert}")

            if table_data:
                st.table(table_data)


# =====================================================================
# TAB 5: EVIDENCE EXPLORER
# =====================================================================

with tab5:
    st.markdown("## Evidence Explorer")
    st.caption("RAG-powered Q&A across all biomarker knowledge collections.")

    # Sidebar collection filters (using columns inside the tab)
    collection_names = [
        "biomarker_reference",
        "biomarker_genetic_variants",
        "biomarker_pgx_rules",
        "biomarker_disease_trajectories",
        "biomarker_clinical_evidence",
        "biomarker_nutrition",
        "biomarker_drug_interactions",
        "biomarker_aging_markers",
        "biomarker_genotype_adjustments",
        "biomarker_monitoring",
        "biomarker_critical_values",
        "biomarker_discordance_rules",
        "biomarker_aj_carrier_screening",
    ]

    with st.expander("Collection Filters", expanded=False):
        selected_collections = []
        c1, c2 = st.columns(2)
        for i, coll in enumerate(collection_names):
            short = coll.replace("biomarker_", "").replace("_", " ").title()
            col = c1 if i % 2 == 0 else c2
            with col:
                if st.checkbox(short, value=True, key=f"t5_coll_{coll}"):
                    selected_collections.append(coll)

    # Initialize conversation history
    if "evidence_history" not in st.session_state:
        st.session_state["evidence_history"] = []

    question = st.text_input(
        "Ask a question about biomarkers, genetics, or pharmacogenomics:",
        placeholder="e.g., How does MTHFR C677T affect folate metabolism?",
        key="t5_question",
    )

    if st.button("Search Evidence", type="primary", key="t5_run"):
        if not question:
            st.warning("Please enter a question.")
        elif engine is None or engine.get("embedder") is None:
            st.error("Engine or embedder not initialized.")
        else:
            with st.spinner("Searching evidence..."):
                # Embed the question
                prefix = "Represent this sentence for searching relevant passages: "
                query_vec = engine["embedder"].embed_text(prefix + question)

                # Search selected collections
                all_hits = []
                for coll in selected_collections:
                    try:
                        hits = engine["manager"].search(
                            collection_name=coll,
                            query_embedding=query_vec,
                            top_k=5,
                            score_threshold=0.4,
                        )
                        all_hits.extend(hits)
                    except Exception:
                        pass

                # Sort by score
                all_hits.sort(key=lambda x: x.get("score", 0), reverse=True)

                # Build evidence context
                evidence_text = ""
                for i, hit in enumerate(all_hits[:15]):
                    coll_name = hit.get("collection", "unknown")
                    score = hit.get("score", 0)
                    text = hit.get("text_chunk", hit.get("text_summary", hit.get("finding", "")))
                    evidence_text += f"\n[{i+1}] ({coll_name}, score={score:.3f}): {text}\n"

                # Generate LLM response if available
                if engine.get("llm_client") and evidence_text:
                    system_prompt = (
                        "You are a precision medicine biomarker expert. Answer the question "
                        "using the provided evidence. Cite sources by their number [1], [2], etc. "
                        "Be concise and clinically accurate."
                    )
                    prompt = f"Question: {question}\n\nEvidence:\n{evidence_text}\n\nAnswer:"

                    response_container = st.empty()
                    full_response = ""
                    for chunk in engine["llm_client"].generate_stream(
                        prompt, system_prompt=system_prompt
                    ):
                        full_response += chunk
                        response_container.markdown(full_response)

                    # Store in history
                    st.session_state["evidence_history"].append({
                        "question": question,
                        "answer": full_response,
                        "hits": len(all_hits),
                        "timestamp": datetime.now().strftime("%H:%M:%S"),
                    })
                elif evidence_text:
                    st.markdown("### Evidence Found")
                    st.markdown(evidence_text)
                else:
                    st.info("No matching evidence found in selected collections.")

    # Show conversation history
    if st.session_state.get("evidence_history"):
        with st.expander("Conversation History", expanded=False):
            for entry in reversed(st.session_state["evidence_history"]):
                st.markdown(f"**Q** ({entry['timestamp']}): {entry['question']}")
                st.markdown(f"**A**: {entry['answer'][:200]}...")
                st.markdown("---")


# =====================================================================
# TAB 6: REPORTS
# =====================================================================

with tab6:
    st.markdown("## Report Generation")
    st.caption("Generate and export comprehensive biomarker analysis reports.")

    results = st.session_state.get("analysis_results")
    patient = st.session_state.get("patient_info")

    if not results or not patient:
        st.info(
            "No analysis results available. Run a Full Analysis in the "
            "Biomarker Analysis tab first."
        )
    else:
        if st.button("Generate Full Report", type="primary", key="t6_gen"):
            with st.spinner("Generating 12-section report..."):
                from src.export import (
                    export_csv,
                    export_fhir_diagnostic_report,
                    export_pdf,
                )
                from src.models import (
                    AnalysisResult,
                    BiologicalAgeResult,
                    DiseaseCategory,
                    DiseaseTrajectoryResult,
                    GenotypeAdjustmentResult,
                    MetabolizerPhenotype,
                    PatientProfile,
                    PGxResult,
                    RiskLevel,
                )
                from src.report_generator import ReportGenerator

                # Build PatientProfile model from session state
                profile = PatientProfile(
                    patient_id=patient.get("patient_id", "PATIENT"),
                    age=int(patient["age"]),
                    sex=patient["sex"],
                    biomarkers=patient.get("biomarkers", {}),
                    genotypes=patient.get("genotypes", {}),
                    star_alleles=patient.get("star_alleles", {}),
                )

                # Build BiologicalAgeResult from raw dict
                ba_raw = results.get("biological_age", {})
                pheno = ba_raw.get("phenoage", ba_raw)
                bio_age_result = BiologicalAgeResult(
                    chronological_age=int(patient["age"]),
                    biological_age=float(pheno.get("biological_age", patient["age"])),
                    age_acceleration=float(pheno.get("age_acceleration", 0)),
                    phenoage_score=float(pheno.get("phenoage_score",
                                                    pheno.get("mortality_score", 0))),
                    grimage_score=pheno.get("grimage_score"),
                    mortality_risk=float(pheno.get("mortality_risk", 0)),
                    aging_drivers=pheno.get("aging_drivers",
                                            ba_raw.get("aging_drivers", [])),
                )

                # Build DiseaseTrajectoryResult list
                disease_trajectories = []
                disease_map = {
                    "diabetes": DiseaseCategory.DIABETES,
                    "cardiovascular": DiseaseCategory.CARDIOVASCULAR,
                    "liver": DiseaseCategory.LIVER,
                    "thyroid": DiseaseCategory.THYROID,
                    "iron": DiseaseCategory.IRON,
                    "nutritional": DiseaseCategory.NUTRITIONAL,
                    "kidney": DiseaseCategory.KIDNEY,
                    "bone_health": DiseaseCategory.BONE_HEALTH,
                    "cognitive": DiseaseCategory.COGNITIVE,
                }
                for traj_raw in results.get("disease_trajectories", []):
                    disease_key = traj_raw.get("disease", "").lower()
                    disease_cat = disease_map.get(disease_key)
                    if not disease_cat:
                        continue
                    risk_val = traj_raw.get("risk_level", "normal").lower()
                    try:
                        risk_level = RiskLevel(risk_val)
                    except ValueError:
                        risk_level = RiskLevel.NORMAL
                    disease_trajectories.append(DiseaseTrajectoryResult(
                        disease=disease_cat,
                        risk_level=risk_level,
                        current_markers=traj_raw.get("current_markers", {}),
                        genetic_risk_factors=traj_raw.get(
                            "genetic_risk_factors",
                            traj_raw.get("genetic_factors", []),
                        ),
                        years_to_onset_estimate=traj_raw.get(
                            "years_to_onset_estimate"),
                        intervention_recommendations=traj_raw.get(
                            "intervention_recommendations",
                            traj_raw.get("recommendations", []),
                        ),
                    ))

                # Build PGx results from session (may come from PGx tab)
                pgx_results = []
                for pgx_raw in results.get("pgx_results", []):
                    phenotype_val = pgx_raw.get("phenotype", "normal").lower()
                    try:
                        phenotype = MetabolizerPhenotype(phenotype_val)
                    except ValueError:
                        phenotype = MetabolizerPhenotype.NORMAL
                    pgx_results.append(PGxResult(
                        gene=pgx_raw.get("gene", ""),
                        star_alleles=pgx_raw.get("star_alleles", ""),
                        phenotype=phenotype,
                        drugs_affected=pgx_raw.get("drugs_affected", []),
                    ))

                # Build genotype adjustments
                genotype_adjustments = []
                for adj_raw in results.get("genotype_adjustments", []):
                    genotype_adjustments.append(GenotypeAdjustmentResult(
                        biomarker=adj_raw.get("biomarker", ""),
                        standard_range=adj_raw.get("standard_range", ""),
                        adjusted_range=adj_raw.get("adjusted_range", ""),
                        genotype=adj_raw.get("genotype", ""),
                        gene=adj_raw.get("gene", ""),
                        rationale=adj_raw.get("rationale", ""),
                    ))

                # Assemble AnalysisResult
                analysis = AnalysisResult(
                    patient_profile=profile,
                    biological_age=bio_age_result,
                    disease_trajectories=disease_trajectories,
                    pgx_results=pgx_results,
                    genotype_adjustments=genotype_adjustments,
                    critical_alerts=[a.to_alert_string() if hasattr(a, 'to_alert_string') else str(a) for a in results.get("critical_alerts", [])],
                )

                # Generate the full 12-section report
                generator = ReportGenerator()
                report_md = generator.generate(analysis, profile)

                # Store for display and export
                st.session_state["full_report"] = report_md
                st.session_state["report_analysis"] = analysis
                st.session_state["report_profile"] = profile

        # Display report and export buttons
        if "full_report" in st.session_state:
            st.markdown(st.session_state["full_report"])

            st.markdown("---")
            st.markdown("### Export Options")

            analysis = st.session_state.get("report_analysis")
            profile = st.session_state.get("report_profile")

            c1, c2, c3, c4 = st.columns(4)

            with c1:
                st.download_button(
                    "Markdown",
                    data=st.session_state["full_report"].encode("utf-8"),
                    file_name=f"biomarker_report_{datetime.now().strftime('%Y%m%d_%H%M')}.md",
                    mime="text/markdown",
                    key="t6_md",
                )

            with c2:
                try:
                    import reportlab  # noqa: F401
                except ImportError:
                    st.warning("reportlab not installed — PDF will be plain text")
                from src.export import export_pdf
                pdf_bytes = export_pdf(st.session_state["full_report"])
                st.download_button(
                    "PDF Report",
                    data=pdf_bytes,
                    file_name=f"biomarker_report_{datetime.now().strftime('%Y%m%d_%H%M')}.pdf",
                    mime="application/pdf",
                    key="t6_pdf",
                )

            with c3:
                if analysis and profile:
                    from src.export import export_fhir_diagnostic_report
                    fhir_json = export_fhir_diagnostic_report(analysis, profile)
                    st.download_button(
                        "FHIR R4",
                        data=fhir_json,
                        file_name=f"biomarker_fhir_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                        mime="application/json",
                        key="t6_fhir",
                    )

            with c4:
                if analysis:
                    from src.export import export_csv
                    csv_bytes = export_csv(analysis)
                    st.download_button(
                        "CSV",
                        data=csv_bytes,
                        file_name=f"biomarker_report_{datetime.now().strftime('%Y%m%d_%H%M')}.csv",
                        mime="text/csv",
                        key="t6_csv",
                    )


# =====================================================================
# TAB 7: PATIENT 360
# =====================================================================

with tab7:
    st.markdown("### \U0001f310 Unified Patient 360 Dashboard")
    st.caption("Cross-agent intelligence view combining genomics, biomarkers, drug candidates, and clinical evidence")

    # Check if we have analysis results to display
    _analysis = st.session_state.get("analysis_results")
    _patient = st.session_state.get("patient_info")

    if _analysis and _patient:
        from patient_360 import (
            render_biomarker_panel as render_360_biomarker_panel,
        )
        from patient_360 import (
            render_concordance_matrix,
            render_drug_candidates,
            render_evidence_chain,
            render_header,
            render_pipeline_status,
        )

        _pid = _patient.get("patient_id", "UNKNOWN")
        render_header(_pid, f"Patient {_pid}")

        # Pipeline status -- determine completed stages from available results
        _completed_stages = ["Biomarker"]
        _ba = _analysis.get("biological_age")
        if _ba:
            _completed_stages.append("Biomarker")  # already present, but confirms
        if _analysis.get("disease_trajectories"):
            _completed_stages.append("Biomarker")
        if _analysis.get("pgx_results"):
            _completed_stages.append("Drug Discovery")
        render_pipeline_status(list(set(_completed_stages)))

        # Biomarker panel -- unpack from analysis_results
        if _ba:
            _pheno = _ba.get("phenoage", _ba)
            render_360_biomarker_panel(
                biological_age=float(_pheno.get("biological_age", _patient.get("age", 0))),
                chronological_age=int(_patient.get("age", 0)),
                age_acceleration=float(_pheno.get("age_acceleration", 0)),
                mortality_risk=str(_pheno.get("mortality_risk", "UNKNOWN")).upper(),
                disease_trajectories=_analysis.get("disease_trajectories", []),
                top_drivers=_pheno.get("aging_drivers", _ba.get("aging_drivers", [])),
            )

        # Drug candidates (if any exist in session)
        render_drug_candidates(
            candidates=_analysis.get("drug_candidates", []),
            pgx_filtered=_analysis.get("pgx_filtered_candidates", []),
        )

        # Evidence chain
        render_evidence_chain(st.session_state.get("evidence_history", []))

        # Concordance
        render_concordance_matrix(_analysis.get("concordance_scores", {}))
    else:
        st.info("Run an analysis in the other tabs first, then return here for the unified Patient 360 view.")
        st.markdown("---")
        st.markdown("**Or explore with demo data:**")
        if st.button("\U0001f3af Load Demo Patient 360", key="demo_360"):
            from patient_360 import render_demo_patient_360
            render_demo_patient_360()

# =====================================================================
# TAB 8: LONGITUDINAL TRACKING
# =====================================================================

with tab8:
    st.header("Longitudinal Biomarker Tracking")
    st.markdown("Track biomarker trends across multiple visits to identify improving, stable, and crisis patterns.")

    _long_path = PROJECT_ROOT / "data" / "reference" / "biomarker_longitudinal_tracking.json"
    if _long_path.exists():
        _long_data = json.loads(_long_path.read_text())

        _long_patient_id = _long_data.get("patient_id", "Unknown")
        visits = _long_data.get("visits", [])
        tracked = _long_data.get("tracked_biomarkers", [])

        st.subheader(f"Patient: {_long_patient_id}")

        # Visit timeline
        visit_cols = st.columns(len(visits))
        for i, visit in enumerate(visits):
            with visit_cols[i]:
                st.metric(visit.get("label", f"Visit {i+1}"), visit.get("date", ""))

        st.markdown("---")

        # Trend filter
        trend_filter = st.selectbox(
            "Filter by trend",
            ["All", "Improving", "Crisis", "Stable"],
            key="trend_filter",
        )

        # Filter biomarkers by trend
        if trend_filter != "All":
            filtered = [b for b in tracked if b.get("trend", "").lower() == trend_filter.lower()]
        else:
            filtered = tracked

        # Summary metrics
        improving = len([b for b in tracked if b.get("trend") == "improving"])
        crisis = len([b for b in tracked if b.get("trend") == "crisis"])
        stable = len([b for b in tracked if b.get("trend") == "stable"])

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total Tracked", len(tracked))
        m2.metric("Improving", improving)
        m3.metric("Crisis", crisis)
        m4.metric("Stable", stable)

        st.markdown("---")

        # Display biomarkers as expandable cards
        for bio in filtered:
            name = bio.get("biomarker", "Unknown")
            trend = bio.get("trend", "unknown")
            values = bio.get("values", [])
            unit = bio.get("unit", "")
            pct_change = bio.get("pct_change_total", 0)
            note = bio.get("clinical_note", "")

            # Color-code by trend
            if trend == "improving":
                icon = "arrow_down_small" if pct_change < 0 else "arrow_up_small"
                color = "green"
            elif trend == "crisis":
                icon = "warning"
                color = "red"
            else:
                icon = "heavy_minus_sign"
                color = "gray"

            _expander_label = (
                f"{'🟢' if trend == 'improving' else '🔴' if trend == 'crisis' else '⚪'} "
                f"{name} — {trend.upper()} ({pct_change:+.1f}%)"
                if pct_change
                else f"{'🟢' if trend == 'improving' else '🔴' if trend == 'crisis' else '⚪'} "
                     f"{name} — {trend.upper()}"
            )

            with st.expander(_expander_label):
                # Sparkline using st.line_chart
                if len(values) >= 2:
                    import pandas as pd
                    visit_labels = [v.get("label", f"V{i+1}") for i, v in enumerate(visits)]
                    # Handle None values
                    chart_values = [v if v is not None else 0 for v in values[:len(visit_labels)]]
                    df = pd.DataFrame({"Visit": visit_labels[:len(chart_values)], name: chart_values})
                    df = df.set_index("Visit")
                    st.line_chart(df)

                col1, col2, col3 = st.columns(3)
                for i, visit in enumerate(visits):
                    if i < len(values):
                        val = values[i]
                        label = visit.get("label", f"Visit {i+1}")
                        [col1, col2, col3][i].metric(label, f"{val} {unit}" if val is not None else "N/A")

                if note:
                    st.info(note)
    else:
        st.warning("Longitudinal tracking data not found. Run the data generation scripts first.")
