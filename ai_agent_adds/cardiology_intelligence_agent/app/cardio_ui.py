"""Cardiology Intelligence Agent -- 10-Tab Streamlit UI.

NVIDIA dark-themed cardiovascular clinical decision support interface with
RAG-powered queries, six risk calculators, GDMT optimization, eight
clinical workflows, multi-modal imaging integration, cardio-oncology
surveillance, evidence exploration, and multi-format report generation.

Usage:
    streamlit run app/cardio_ui.py --server.port 8536

Author: Adam Jones
Date: March 2026
"""

import json
import os
from typing import Optional

import requests
import streamlit as st

# =====================================================================
# Configuration
# =====================================================================

API_BASE = os.environ.get("CARDIO_API_BASE", "http://localhost:8126")

NVIDIA_THEME = {
    "bg_primary": "#1a1a2e",
    "bg_secondary": "#16213e",
    "bg_card": "#0f3460",
    "text_primary": "#e0e0e0",
    "text_secondary": "#a0a0b0",
    "accent": "#76b900",
    "accent_hover": "#8ed100",
    "danger": "#e74c3c",
    "warning": "#f39c12",
    "info": "#3498db",
    "success": "#76b900",
}


# =====================================================================
# Page Config & Custom CSS
# =====================================================================

st.set_page_config(
    page_title="Cardiology Intelligence Agent",
    page_icon="🫀",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(f"""
<style>
    /* Main background */
    .stApp {{
        background-color: {NVIDIA_THEME['bg_primary']};
        color: {NVIDIA_THEME['text_primary']};
    }}

    /* Sidebar */
    section[data-testid="stSidebar"] {{
        background-color: {NVIDIA_THEME['bg_secondary']};
    }}
    section[data-testid="stSidebar"] .stMarkdown {{
        color: {NVIDIA_THEME['text_primary']};
    }}

    /* Cards */
    div[data-testid="stMetric"] {{
        background-color: {NVIDIA_THEME['bg_card']};
        border: 1px solid {NVIDIA_THEME['accent']};
        border-radius: 8px;
        padding: 12px;
    }}
    div[data-testid="stMetric"] label {{
        color: {NVIDIA_THEME['text_secondary']};
    }}
    div[data-testid="stMetric"] div[data-testid="stMetricValue"] {{
        color: {NVIDIA_THEME['accent']};
    }}

    /* Headers */
    h1, h2, h3 {{
        color: {NVIDIA_THEME['accent']} !important;
    }}

    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {{
        gap: 4px;
        background-color: {NVIDIA_THEME['bg_secondary']};
        border-radius: 8px;
        padding: 4px;
    }}
    .stTabs [data-baseweb="tab"] {{
        color: {NVIDIA_THEME['text_secondary']};
        border-radius: 6px;
    }}
    .stTabs [aria-selected="true"] {{
        background-color: {NVIDIA_THEME['bg_card']};
        color: {NVIDIA_THEME['accent']} !important;
    }}

    /* Buttons */
    .stButton > button {{
        background-color: {NVIDIA_THEME['accent']};
        color: #000000;
        font-weight: 600;
        border: none;
        border-radius: 6px;
    }}
    .stButton > button:hover {{
        background-color: {NVIDIA_THEME['accent_hover']};
        color: #000000;
    }}

    /* Expanders */
    .streamlit-expanderHeader {{
        background-color: {NVIDIA_THEME['bg_card']};
        color: {NVIDIA_THEME['text_primary']};
        border-radius: 6px;
    }}

    /* Dataframes */
    .stDataFrame {{
        border: 1px solid {NVIDIA_THEME['accent']};
        border-radius: 6px;
    }}

    /* Text inputs */
    .stTextInput > div > div > input,
    .stTextArea > div > div > textarea,
    .stSelectbox > div > div {{
        background-color: {NVIDIA_THEME['bg_secondary']};
        color: {NVIDIA_THEME['text_primary']};
        border: 1px solid {NVIDIA_THEME['bg_card']};
    }}

    /* Success/warning/error boxes */
    .success-box {{
        background-color: rgba(118, 185, 0, 0.15);
        border-left: 4px solid {NVIDIA_THEME['accent']};
        padding: 12px;
        border-radius: 4px;
        margin: 8px 0;
    }}
    .warning-box {{
        background-color: rgba(243, 156, 18, 0.15);
        border-left: 4px solid {NVIDIA_THEME['warning']};
        padding: 12px;
        border-radius: 4px;
        margin: 8px 0;
    }}
    .info-box {{
        background-color: rgba(52, 152, 219, 0.15);
        border-left: 4px solid {NVIDIA_THEME['info']};
        padding: 12px;
        border-radius: 4px;
        margin: 8px 0;
    }}
</style>
""", unsafe_allow_html=True)

st.warning(
    "**Clinical Decision Support Tool** — This system provides evidence-based guidance "
    "for research and clinical decision support only. All recommendations must be verified "
    "by a qualified healthcare professional. Not FDA-cleared. Not a substitute for professional "
    "clinical judgment."
)


# =====================================================================
# API Helpers
# =====================================================================

def api_get(endpoint: str, params: dict = None) -> Optional[dict]:
    """GET request to cardiology API."""
    try:
        resp = requests.get(f"{API_BASE}{endpoint}", params=params, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        st.error(f"API error: {e}")
        return None


def api_post(endpoint: str, data: dict) -> Optional[dict]:
    """POST request to cardiology API."""
    try:
        resp = requests.post(f"{API_BASE}{endpoint}", json=data, timeout=60)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        st.error(f"API error: {e}")
        return None


def risk_badge(category: str) -> str:
    """Return colored HTML badge for risk category."""
    colors = {
        "low": NVIDIA_THEME["success"],
        "borderline": NVIDIA_THEME["info"],
        "low-moderate": NVIDIA_THEME["info"],
        "moderate": NVIDIA_THEME["warning"],
        "moderate-high": NVIDIA_THEME["warning"],
        "intermediate": NVIDIA_THEME["warning"],
        "high": NVIDIA_THEME["danger"],
        "very_high": NVIDIA_THEME["danger"],
    }
    color = colors.get(category, NVIDIA_THEME["text_secondary"])
    return f'<span style="background-color:{color};color:#000;padding:4px 12px;border-radius:12px;font-weight:600;">{category.upper()}</span>'


# =====================================================================
# Sidebar
# =====================================================================

with st.sidebar:
    st.image("https://developer.nvidia.com/favicon.ico", width=40)
    st.markdown("### Cardiology Intelligence Agent")
    st.markdown(f'<span style="color:{NVIDIA_THEME["text_secondary"]};">RAG-powered cardiovascular CDS</span>', unsafe_allow_html=True)
    st.divider()

    # Workflow selector
    workflow_options = [
        "None", "CAD Assessment", "Heart Failure", "Valvular Disease",
        "Arrhythmia & EP", "Cardiac MRI", "Stress Testing",
        "Prevention", "Cardio-Oncology",
    ]
    selected_workflow = st.selectbox("Active Workflow", workflow_options, index=0)
    workflow_map = {
        "CAD Assessment": "cad", "Heart Failure": "heart_failure",
        "Valvular Disease": "valvular", "Arrhythmia & EP": "arrhythmia",
        "Cardiac MRI": "cardiac_mri", "Stress Testing": "stress_test",
        "Prevention": "prevention", "Cardio-Oncology": "cardio_oncology",
    }
    active_workflow = workflow_map.get(selected_workflow)

    # Collection filter
    st.markdown("#### Collection Filter")
    all_collections = [
        "cardio_guidelines", "cardio_trials", "cardio_pharmacology",
        "cardio_imaging", "cardio_electrophysiology", "cardio_genetics",
        "cardio_biomarkers", "cardio_interventional", "cardio_prevention",
        "cardio_heart_failure", "cardio_valvular", "cardio_oncology",
    ]
    selected_collections = st.multiselect(
        "Search collections",
        all_collections,
        default=all_collections,
        label_visibility="collapsed",
    )

    # Patient context toggle
    st.divider()
    use_patient_ctx = st.toggle("Enable Patient Context", value=False)
    patient_context = {}
    if use_patient_ctx:
        patient_context["age"] = st.number_input("Age", 18, 120, 65)
        patient_context["sex"] = st.selectbox("Sex", ["male", "female"])
        patient_context["chief_complaint"] = st.text_input("Chief Complaint", "")

    # About
    st.divider()
    with st.expander("About"):
        st.markdown("""
        **Cardiology Intelligence Agent** v1.0.0

        Part of the HCLS AI Factory precision medicine platform.
        Provides RAG-powered clinical decision support across 12
        cardiovascular knowledge collections, 6 validated risk
        calculators, GDMT optimization, and 8 clinical workflows.

        Built on NVIDIA DGX Spark with Claude AI.
        """)
        st.markdown(f'<span style="color:{NVIDIA_THEME["text_secondary"]};">Author: Adam Jones | March 2026</span>', unsafe_allow_html=True)


# =====================================================================
# Main Content -- 10 Tabs
# =====================================================================

tab_dashboard, tab_query, tab_risk, tab_hf, tab_cad, tab_arrhythmia, \
    tab_imaging, tab_onc, tab_evidence, tab_reports = st.tabs([
        "Dashboard",
        "Clinical Query",
        "Risk Calculator",
        "Heart Failure",
        "CAD Assessment",
        "Arrhythmia",
        "Imaging",
        "Cardio-Oncology",
        "Evidence Explorer",
        "Report Generator",
    ])


# ── Tab 1: Dashboard ──

with tab_dashboard:
    st.header("Cardiovascular Decision Support Dashboard")

    # System health metrics
    health = api_get("/health")
    if health:
        col1, col2, col3, col4, col5, col6 = st.columns(6)
        col1.metric("Status", health.get("status", "unknown").upper())
        col2.metric("Collections", health.get("collections", 0))
        col3.metric("Vectors", f"{health.get('total_vectors', 0):,}")
        col4.metric("Workflows", health.get("workflows", 0))
        col5.metric("Risk Calcs", health.get("risk_calculators", 0))
        col6.metric("Engine", "Ready" if health.get("components", {}).get("rag_engine") == "ready" else "Offline")
    else:
        st.warning("API unavailable -- showing static dashboard")
        col1, col2, col3, col4, col5, col6 = st.columns(6)
        col1.metric("Status", "OFFLINE")
        col2.metric("Collections", 12)
        col3.metric("Vectors", "18,965")
        col4.metric("Workflows", 8)
        col5.metric("Risk Calcs", 6)
        col6.metric("Engine", "Offline")

    st.divider()

    # Collections overview
    col_left, col_right = st.columns(2)

    with col_left:
        st.subheader("Knowledge Collections")
        collections = api_get("/collections")
        if collections:
            for coll in collections.get("collections", []):
                name = coll.get("name", "Unknown")
                count = coll.get("count", 0)
                source = coll.get("source", "")
                st.markdown(f"**{name}** -- {count:,} records  \n_{source}_")
        else:
            st.info("Connect to API to view collection details.")

    with col_right:
        st.subheader("Available Workflows")
        workflows = api_get("/workflows")
        if workflows:
            for wf in workflows.get("workflows", []):
                with st.expander(f"{wf['name']} ({wf['id']})"):
                    st.write(wf.get("description", ""))
                    if wf.get("risk_calculators"):
                        st.write(f"**Risk calculators:** {', '.join(wf['risk_calculators'])}")
        else:
            for name in workflow_options[1:]:
                st.markdown(f"- {name}")

    st.divider()

    # Quick stats
    st.subheader("Agent Capabilities")
    cap1, cap2, cap3, cap4 = st.columns(4)
    cap1.metric("Conditions Covered", 45)
    cap2.metric("Drug Classes", 32)
    cap3.metric("Cardio Genes", 56)
    cap4.metric("Guidelines", 63)


# ── Tab 2: Clinical Query ──

with tab_query:
    st.header("Clinical Query")
    st.markdown("Ask any cardiovascular clinical question. The RAG engine searches across all collections and cites relevant guidelines.")

    query_text = st.text_area(
        "Enter your clinical question",
        placeholder="e.g., What is the recommended initial therapy for a 58-year-old male with newly diagnosed HFrEF (LVEF 30%) and NYHA Class III symptoms?",
        height=120,
    )

    col_q1, col_q2, col_q3 = st.columns(3)
    with col_q1:
        query_workflow = st.selectbox(
            "Workflow context",
            ["Auto-detect"] + list(workflow_map.keys()),
            key="query_workflow",
        )
    with col_q2:
        query_top_k = st.slider("Evidence passages", 1, 20, 5, key="query_top_k")
    with col_q3:
        include_guidelines = st.checkbox("Include guideline citations", value=True, key="query_guidelines")

    if st.button("Submit Query", key="btn_query", type="primary"):
        if query_text.strip():
            with st.spinner("Searching cardiovascular knowledge base..."):
                wf_type = workflow_map.get(query_workflow) if query_workflow != "Auto-detect" else None
                payload = {
                    "question": query_text,
                    "workflow_type": wf_type or active_workflow,
                    "patient_context": patient_context if use_patient_ctx else None,
                    "top_k": query_top_k,
                    "include_guidelines": include_guidelines,
                }
                result = api_post("/v1/cardio/query", payload)

            if result:
                st.markdown("### Answer")
                st.markdown(result.get("answer", "No answer returned."))

                if result.get("confidence"):
                    st.progress(result["confidence"], text=f"Confidence: {result['confidence']:.0%}")

                if result.get("guidelines_cited"):
                    st.markdown("### Guidelines Cited")
                    for g in result["guidelines_cited"]:
                        st.markdown(f"- {g}")

                if result.get("evidence"):
                    with st.expander(f"Evidence ({len(result['evidence'])} passages)"):
                        for i, ev in enumerate(result["evidence"], 1):
                            st.markdown(f"**[{i}]** _{ev.get('collection', '')}_ (score: {ev.get('score', 'N/A')})")
                            st.markdown(ev.get("text", ""))
                            st.divider()
        else:
            st.warning("Please enter a question.")


# ── Tab 3: Risk Calculator ──

with tab_risk:
    st.header("Cardiovascular Risk Calculators")

    calc_choice = st.selectbox("Select Calculator", [
        "ASCVD (10-year risk)",
        "HEART Score (ACS triage)",
        "CHA2DS2-VASc (AF stroke risk)",
        "HAS-BLED (Bleeding risk)",
        "MAGGIC (HF mortality)",
        "EuroSCORE II (Surgical risk)",
    ])

    st.divider()

    if calc_choice == "ASCVD (10-year risk)":
        st.subheader("Pooled Cohort Equations -- ASCVD 10-Year Risk")
        col_a1, col_a2 = st.columns(2)
        with col_a1:
            ascvd_age = st.number_input("Age", 40, 79, 55, key="ascvd_age")
            ascvd_sex = st.selectbox("Sex", ["male", "female"], key="ascvd_sex")
            ascvd_race = st.selectbox("Race", ["white", "african_american", "other"], key="ascvd_race")
            ascvd_tc = st.number_input("Total Cholesterol (mg/dL)", 100.0, 400.0, 220.0, key="ascvd_tc")
            ascvd_hdl = st.number_input("HDL Cholesterol (mg/dL)", 20.0, 150.0, 50.0, key="ascvd_hdl")
        with col_a2:
            ascvd_sbp = st.number_input("Systolic BP (mmHg)", 80.0, 250.0, 130.0, key="ascvd_sbp")
            ascvd_bp_tx = st.checkbox("On BP treatment", key="ascvd_bp_tx")
            ascvd_dm = st.checkbox("Diabetes", key="ascvd_dm")
            ascvd_smoke = st.checkbox("Current smoker", key="ascvd_smoke")

        if st.button("Calculate ASCVD Risk", key="btn_ascvd", type="primary"):
            with st.spinner("Calculating..."):
                result = api_post("/v1/cardio/risk/ascvd", {
                    "age": ascvd_age, "sex": ascvd_sex, "race": ascvd_race,
                    "total_cholesterol": ascvd_tc, "hdl_cholesterol": ascvd_hdl,
                    "systolic_bp": ascvd_sbp, "bp_treatment": ascvd_bp_tx,
                    "diabetes": ascvd_dm, "smoker": ascvd_smoke,
                })
            if result:
                st.markdown(risk_badge(result.get("risk_category", "")), unsafe_allow_html=True)
                r1, r2 = st.columns(2)
                r1.metric("10-Year ASCVD Risk", f"{result.get('score', 0)}%")
                r2.metric("Category", result.get("risk_category", "").upper())
                st.markdown(f"**Interpretation:** {result.get('interpretation', '')}")
                st.markdown("**Recommendations:**")
                for rec in result.get("recommendations", []):
                    st.markdown(f"- {rec}")

    elif calc_choice == "HEART Score (ACS triage)":
        st.subheader("HEART Score -- Acute Chest Pain Evaluation")
        col_h1, col_h2 = st.columns(2)
        with col_h1:
            heart_history = st.selectbox("History", [0, 1, 2], format_func=lambda x: ["Slightly suspicious", "Moderately suspicious", "Highly suspicious"][x], key="heart_hist")
            heart_ecg = st.selectbox("ECG", [0, 1, 2], format_func=lambda x: ["Normal", "Non-specific changes", "Significant deviation"][x], key="heart_ecg")
            heart_age = st.number_input("Age", 18, 120, 55, key="heart_age")
        with col_h2:
            heart_rf = st.number_input("Risk factors (count)", 0, 5, 2, key="heart_rf")
            heart_trop = st.selectbox("Troponin", [0, 1, 2], format_func=lambda x: ["Normal", "1-3x ULN", ">3x ULN"][x], key="heart_trop")

        if st.button("Calculate HEART Score", key="btn_heart", type="primary"):
            with st.spinner("Calculating..."):
                result = api_post("/v1/cardio/risk/heart-score", {
                    "history_score": heart_history, "ecg_score": heart_ecg,
                    "age": heart_age, "risk_factors": heart_rf, "troponin_score": heart_trop,
                })
            if result:
                st.markdown(risk_badge(result.get("risk_category", "")), unsafe_allow_html=True)
                r1, r2 = st.columns(2)
                r1.metric("HEART Score", int(result.get("score", 0)))
                r2.metric("Category", result.get("risk_category", "").upper())
                st.markdown(f"**Interpretation:** {result.get('interpretation', '')}")
                st.markdown("**Recommendations:**")
                for rec in result.get("recommendations", []):
                    st.markdown(f"- {rec}")

    elif calc_choice == "CHA2DS2-VASc (AF stroke risk)":
        st.subheader("CHA2DS2-VASc -- Stroke Risk in Atrial Fibrillation")
        col_c1, col_c2 = st.columns(2)
        with col_c1:
            cha_chf = st.checkbox("CHF / LV dysfunction", key="cha_chf")
            cha_htn = st.checkbox("Hypertension", key="cha_htn")
            cha_age = st.number_input("Age", 18, 120, 70, key="cha_age")
            cha_dm = st.checkbox("Diabetes", key="cha_dm")
        with col_c2:
            cha_stroke = st.checkbox("Prior Stroke / TIA", key="cha_stroke")
            cha_vasc = st.checkbox("Vascular disease", key="cha_vasc")
            cha_sex = st.selectbox("Sex", ["male", "female"], key="cha_sex")

        if st.button("Calculate CHA2DS2-VASc", key="btn_cha", type="primary"):
            with st.spinner("Calculating..."):
                result = api_post("/v1/cardio/risk/cha2ds2-vasc", {
                    "chf": cha_chf, "hypertension": cha_htn, "age": cha_age,
                    "diabetes": cha_dm, "stroke_tia": cha_stroke,
                    "vascular_disease": cha_vasc, "sex": cha_sex,
                })
            if result:
                st.markdown(risk_badge(result.get("risk_category", "")), unsafe_allow_html=True)
                r1, r2 = st.columns(2)
                r1.metric("CHA2DS2-VASc", int(result.get("score", 0)))
                r2.metric("Annual Stroke Rate", f"{result.get('details', {}).get('annual_stroke_rate_pct', 'N/A')}%")
                st.markdown(f"**Interpretation:** {result.get('interpretation', '')}")
                st.markdown("**Recommendations:**")
                for rec in result.get("recommendations", []):
                    st.markdown(f"- {rec}")

    elif calc_choice == "HAS-BLED (Bleeding risk)":
        st.subheader("HAS-BLED -- Bleeding Risk Assessment")
        col_b1, col_b2 = st.columns(2)
        with col_b1:
            hb_htn = st.checkbox("Uncontrolled hypertension (>160 mmHg)", key="hb_htn")
            hb_renal = st.checkbox("Renal disease (dialysis, transplant, Cr >2.3)", key="hb_renal")
            hb_liver = st.checkbox("Liver disease (cirrhosis, bilirubin >2x, AST/ALT >3x)", key="hb_liver")
            hb_stroke = st.checkbox("Stroke history", key="hb_stroke")
        with col_b2:
            hb_bleed = st.checkbox("Prior major bleeding", key="hb_bleed")
            hb_inr = st.checkbox("Labile INR (TTR <60%)", key="hb_inr")
            hb_age = st.checkbox("Age > 65", key="hb_age")
            hb_drugs = st.selectbox("Drugs/alcohol", [0, 1, 2], format_func=lambda x: ["None", "One (antiplatelet/NSAID or alcohol)", "Both"][x], key="hb_drugs")

        if st.button("Calculate HAS-BLED", key="btn_hb", type="primary"):
            with st.spinner("Calculating..."):
                result = api_post("/v1/cardio/risk/has-bled", {
                    "hypertension_uncontrolled": hb_htn, "renal_disease": hb_renal,
                    "liver_disease": hb_liver, "stroke_history": hb_stroke,
                    "bleeding_history": hb_bleed, "labile_inr": hb_inr,
                    "age_over_65": hb_age, "drugs_alcohol": hb_drugs,
                })
            if result:
                st.markdown(risk_badge(result.get("risk_category", "")), unsafe_allow_html=True)
                st.metric("HAS-BLED Score", int(result.get("score", 0)))
                st.markdown(f"**Interpretation:** {result.get('interpretation', '')}")
                st.info(result.get("details", {}).get("note", ""))
                st.markdown("**Recommendations:**")
                for rec in result.get("recommendations", []):
                    st.markdown(f"- {rec}")

    elif calc_choice == "MAGGIC (HF mortality)":
        st.subheader("MAGGIC -- Heart Failure Mortality Risk")
        col_m1, col_m2, col_m3 = st.columns(3)
        with col_m1:
            mag_age = st.number_input("Age", 18, 100, 65, key="mag_age")
            mag_sex = st.selectbox("Sex", ["male", "female"], key="mag_sex")
            mag_lvef = st.number_input("LVEF (%)", 5.0, 80.0, 30.0, key="mag_lvef")
            mag_nyha = st.selectbox("NYHA Class", [1, 2, 3, 4], index=2, key="mag_nyha")
            mag_sbp = st.number_input("Systolic BP (mmHg)", 60.0, 250.0, 120.0, key="mag_sbp")
        with col_m2:
            mag_bmi = st.number_input("BMI (kg/m2)", 10.0, 60.0, 27.0, key="mag_bmi")
            mag_cr = st.number_input("Creatinine (mg/dL)", 0.3, 15.0, 1.2, key="mag_cr")
            mag_dm = st.checkbox("Diabetes", key="mag_dm")
            mag_copd = st.checkbox("COPD", key="mag_copd")
            mag_hf18 = st.checkbox("HF duration > 18 months", key="mag_hf18")
        with col_m3:
            mag_smoke = st.checkbox("Current smoker", key="mag_smoke")
            mag_bb = st.checkbox("On beta-blocker", key="mag_bb")
            mag_acei = st.checkbox("On ACEi/ARB", key="mag_acei")

        if st.button("Calculate MAGGIC", key="btn_mag", type="primary"):
            with st.spinner("Calculating..."):
                result = api_post("/v1/cardio/risk/maggic", {
                    "age": mag_age, "sex": mag_sex, "lvef": mag_lvef,
                    "nyha_class": mag_nyha, "systolic_bp": mag_sbp,
                    "bmi": mag_bmi, "creatinine": mag_cr, "diabetes": mag_dm,
                    "copd": mag_copd, "hf_duration_18m": mag_hf18,
                    "smoker": mag_smoke, "beta_blocker": mag_bb, "acei_arb": mag_acei,
                })
            if result:
                st.markdown(risk_badge(result.get("risk_category", "")), unsafe_allow_html=True)
                r1, r2, r3 = st.columns(3)
                r1.metric("MAGGIC Score", int(result.get("score", 0)))
                r2.metric("1-Year Mortality", f"{result.get('details', {}).get('one_year_mortality_pct', 'N/A')}%")
                r3.metric("3-Year Mortality", f"{result.get('details', {}).get('three_year_mortality_pct', 'N/A')}%")
                st.markdown("**Recommendations:**")
                for rec in result.get("recommendations", []):
                    st.markdown(f"- {rec}")

    elif calc_choice == "EuroSCORE II (Surgical risk)":
        st.subheader("EuroSCORE II -- Cardiac Surgical Mortality Risk")
        st.markdown("_Enter patient and procedural factors for operative risk estimation._")
        col_e1, col_e2, col_e3 = st.columns(3)
        with col_e1:
            es_age = st.number_input("Age", 18, 100, 70, key="es_age")
            es_sex = st.selectbox("Sex", ["male", "female"], key="es_sex")
            es_crcl = st.number_input("Creatinine clearance (mL/min)", 0.0, 200.0, 65.0, key="es_crcl")
            es_extra = st.checkbox("Extracardiac arteriopathy", key="es_extra")
            es_mob = st.checkbox("Poor mobility", key="es_mob")
            es_prev = st.checkbox("Previous cardiac surgery", key="es_prev")
        with col_e2:
            es_lung = st.checkbox("Chronic lung disease", key="es_lung")
            es_endo = st.checkbox("Active endocarditis", key="es_endo")
            es_crit = st.checkbox("Critical preop state", key="es_crit")
            es_dm = st.checkbox("Diabetes on insulin", key="es_dm")
            es_nyha = st.selectbox("NYHA Class", [1, 2, 3, 4], index=1, key="es_nyha")
            es_ccs4 = st.checkbox("CCS Class 4 angina", key="es_ccs4")
        with col_e3:
            es_lvef = st.number_input("LVEF (%)", 5.0, 80.0, 55.0, key="es_lvef")
            es_mi = st.checkbox("Recent MI (<90 days)", key="es_mi")
            es_ph = st.selectbox("Pulmonary hypertension", ["none", "moderate", "severe"], key="es_ph")
            es_urg = st.selectbox("Urgency", ["elective", "urgent", "emergency", "salvage"], key="es_urg")
            es_aorta = st.checkbox("Surgery on thoracic aorta", key="es_aorta")

        if st.button("Calculate EuroSCORE II", key="btn_es", type="primary"):
            with st.spinner("Calculating..."):
                result = api_post("/v1/cardio/risk/euroscore", {
                    "age": es_age, "sex": es_sex, "creatinine_clearance": es_crcl,
                    "extracardiac_arteriopathy": es_extra, "poor_mobility": es_mob,
                    "previous_cardiac_surgery": es_prev, "chronic_lung_disease": es_lung,
                    "active_endocarditis": es_endo, "critical_preop_state": es_crit,
                    "diabetes_on_insulin": es_dm, "nyha_class": es_nyha,
                    "ccs_class_4_angina": es_ccs4, "lvef": es_lvef,
                    "recent_mi": es_mi, "pulmonary_hypertension": es_ph,
                    "urgency": es_urg, "thoracic_aorta": es_aorta,
                })
            if result:
                st.markdown(risk_badge(result.get("risk_category", "")), unsafe_allow_html=True)
                r1, r2 = st.columns(2)
                r1.metric("Predicted Mortality", f"{result.get('details', {}).get('predicted_mortality_pct', 'N/A')}%")
                r2.metric("Category", result.get("risk_category", "").upper())
                st.markdown("**Recommendations:**")
                for rec in result.get("recommendations", []):
                    st.markdown(f"- {rec}")


# ── Tab 4: Heart Failure ──

with tab_hf:
    st.header("Heart Failure -- GDMT Optimizer")
    st.markdown("Optimize Guideline-Directed Medical Therapy based on HF phenotype, LVEF, NYHA class, and current medications.")

    col_hf1, col_hf2 = st.columns(2)

    with col_hf1:
        st.subheader("Patient Parameters")
        hf_lvef = st.number_input("LVEF (%)", 5.0, 80.0, 28.0, key="hf_lvef")
        hf_nyha = st.selectbox("NYHA Class", ["I", "II", "III", "IV"], index=2, key="hf_nyha")

        st.subheader("Current Medications")
        hf_meds = []
        n_meds = st.number_input("Number of current medications", 0, 10, 2, key="hf_n_meds")
        for i in range(int(n_meds)):
            col_med_name, col_med_dose = st.columns(2)
            with col_med_name:
                med_name = st.text_input(f"Med {i+1} name", key=f"hf_med_name_{i}",
                                          placeholder="e.g., carvedilol")
            with col_med_dose:
                med_dose = st.text_input(f"Med {i+1} dose", key=f"hf_med_dose_{i}",
                                          placeholder="e.g., 12.5 mg BID")
            if med_name:
                hf_meds.append({"name": med_name, "dose": med_dose})

    with col_hf2:
        st.subheader("Additional Data")
        hf_bp = st.number_input("Systolic BP (mmHg)", 60.0, 250.0, 110.0, key="hf_bp")
        hf_hr = st.number_input("Heart rate (bpm)", 30, 200, 78, key="hf_hr")
        hf_k = st.number_input("Potassium (mEq/L)", 2.5, 7.0, 4.2, key="hf_k")
        hf_cr = st.number_input("Creatinine (mg/dL)", 0.3, 15.0, 1.3, key="hf_cr")
        hf_egfr = st.number_input("eGFR (mL/min/1.73m2)", 5, 150, 55, key="hf_egfr")

    if st.button("Optimize GDMT", key="btn_gdmt", type="primary"):
        with st.spinner("Analyzing GDMT status..."):
            payload = {
                "lvef": hf_lvef,
                "nyha_class": hf_nyha,
                "current_medications": hf_meds,
                "patient_data": {
                    "systolic_bp": hf_bp, "heart_rate": hf_hr,
                    "potassium": hf_k, "creatinine": hf_cr, "egfr": hf_egfr,
                },
            }
            result = api_post("/v1/cardio/gdmt/optimize", payload)

        if result:
            st.markdown(f"### {result.get('hf_phenotype', 'Unknown')} -- LVEF {result.get('lvef')}%, NYHA {result.get('nyha_class')}")

            # Four pillars status
            st.subheader("Four Pillars of GDMT")
            pillar_cols = st.columns(4)
            pillars = result.get("four_pillars_status", {})
            for i, (pillar, status) in enumerate(pillars.items()):
                with pillar_cols[i % 4]:
                    color = NVIDIA_THEME["success"] if "TARGET" in status or "ON" in status else NVIDIA_THEME["danger"]
                    st.markdown(f'<div style="text-align:center;padding:10px;background:{NVIDIA_THEME["bg_card"]};border-radius:8px;border:2px solid {color};">'
                                f'<strong style="color:{NVIDIA_THEME["text_primary"]};">{pillar}</strong><br>'
                                f'<span style="color:{color};font-weight:bold;">{status}</span></div>',
                                unsafe_allow_html=True)

            # Recommendations
            st.subheader("Optimization Recommendations")
            for rec in result.get("recommendations", []):
                if isinstance(rec, dict):
                    with st.expander(f"{rec.get('recommendation', '')} -- {rec.get('medication_class', '')}"):
                        st.write(f"**Current status:** {rec.get('current_status', '')}")
                        if rec.get("target_dose"):
                            st.write(f"**Target dose:** {rec['target_dose']}")
                        if rec.get("evidence_level"):
                            st.write(f"**Evidence:** {rec['evidence_level']}")
                        if rec.get("guideline_ref"):
                            st.write(f"**Guideline:** {rec['guideline_ref']}")
                        if rec.get("caution"):
                            st.warning(rec["caution"])

            st.markdown(f"**Summary:** {result.get('summary', '')}")


# ── Tab 5: CAD Assessment ──

with tab_cad:
    st.header("Coronary Artery Disease Assessment")
    st.markdown("Calcium score interpretation, CAD-RADS classification, and plaque characterization.")

    col_cad1, col_cad2 = st.columns(2)

    with col_cad1:
        st.subheader("Coronary Calcium Score")
        cad_cac = st.number_input("Agatston Calcium Score", 0, 10000, 150, key="cad_cac")
        cad_age = st.number_input("Age", 30, 100, 58, key="cad_age")
        cad_sex = st.selectbox("Sex", ["male", "female"], key="cad_sex")

        # Calcium score interpretation
        if cad_cac == 0:
            cac_interp = "No identifiable plaque -- very low ASCVD risk"
            cac_color = NVIDIA_THEME["success"]
        elif cad_cac <= 100:
            cac_interp = "Mild plaque burden -- low-to-moderate risk"
            cac_color = NVIDIA_THEME["info"]
        elif cad_cac <= 400:
            cac_interp = "Moderate plaque burden -- consider statin therapy"
            cac_color = NVIDIA_THEME["warning"]
        else:
            cac_interp = "Extensive plaque burden -- high risk, statin + aggressive risk factor management"
            cac_color = NVIDIA_THEME["danger"]

        st.markdown(f'<div class="info-box" style="border-left-color:{cac_color};">{cac_interp}</div>', unsafe_allow_html=True)

    with col_cad2:
        st.subheader("CAD-RADS Classification")
        cad_rads = st.selectbox("CAD-RADS", [
            "0 - No stenosis (0%)",
            "1 - Minimal (1-24%)",
            "2 - Mild (25-49%)",
            "3 - Moderate (50-69%)",
            "4A - Severe (70-99%)",
            "4B - Left main >50% or 3-vessel >70%",
            "5 - Total occlusion",
        ], key="cad_rads")

        st.subheader("Plaque Features")
        cad_napkin = st.checkbox("Napkin-ring sign (vulnerable plaque)", key="cad_napkin")
        cad_remodel = st.checkbox("Positive remodeling", key="cad_remodel")
        cad_low_att = st.checkbox("Low-attenuation plaque", key="cad_low_att")
        cad_spotty = st.checkbox("Spotty calcification", key="cad_spotty")

    high_risk_features = sum([cad_napkin, cad_remodel, cad_low_att, cad_spotty])

    if st.button("Run CAD Assessment", key="btn_cad", type="primary"):
        with st.spinner("Assessing..."):
            result = api_post("/v1/cardio/workflow/cad", {
                "calcium_score": cad_cac, "age": cad_age, "sex": cad_sex,
                "cad_rads": cad_rads.split(" - ")[0],
                "high_risk_plaque_features": high_risk_features,
                "napkin_ring": cad_napkin, "positive_remodeling": cad_remodel,
                "low_attenuation": cad_low_att, "spotty_calcification": cad_spotty,
            })
        if result:
            st.json(result)
        else:
            st.info("CAD workflow requires API connection for full analysis.")

    if high_risk_features >= 2:
        st.warning(f"{high_risk_features} high-risk plaque features detected -- consider further evaluation regardless of stenosis severity.")


# ── Tab 6: Arrhythmia ──

with tab_arrhythmia:
    st.header("Arrhythmia & Electrophysiology")
    st.markdown("AF management, stroke/bleeding risk stratification, and anticoagulation decisions.")

    col_arr1, col_arr2 = st.columns(2)

    with col_arr1:
        st.subheader("Rhythm Assessment")
        arr_type = st.selectbox("Arrhythmia Type", [
            "Atrial Fibrillation (paroxysmal)",
            "Atrial Fibrillation (persistent)",
            "Atrial Fibrillation (permanent)",
            "Atrial Flutter",
            "SVT",
            "Ventricular Tachycardia",
            "Other",
        ], key="arr_type")

        st.subheader("ECG Parameters")
        arr_rate = st.number_input("Ventricular rate (bpm)", 30, 300, 88, key="arr_rate")
        arr_qtc = st.number_input("QTc (ms)", 300, 700, 440, key="arr_qtc")
        arr_pr = st.number_input("PR interval (ms)", 0, 400, 160, key="arr_pr")
        arr_qrs = st.number_input("QRS duration (ms)", 60, 250, 96, key="arr_qrs")

    with col_arr2:
        st.subheader("Anticoagulation Decision")
        st.markdown("Use CHA2DS2-VASc and HAS-BLED calculators in the **Risk Calculator** tab for detailed scoring.")

        arr_cha = st.number_input("CHA2DS2-VASc score (pre-calculated)", 0, 9, 3, key="arr_cha")
        arr_hb = st.number_input("HAS-BLED score (pre-calculated)", 0, 9, 2, key="arr_hb")

        if "fibrillation" in arr_type.lower() or "flutter" in arr_type.lower():
            if arr_cha == 0:
                st.markdown('<div class="success-box">No anticoagulation needed (CHA2DS2-VASc = 0)</div>', unsafe_allow_html=True)
            elif arr_cha == 1:
                st.markdown('<div class="info-box">Consider anticoagulation -- shared decision-making</div>', unsafe_allow_html=True)
            else:
                st.markdown('<div class="warning-box">Anticoagulation recommended -- prefer DOAC over warfarin</div>', unsafe_allow_html=True)

            if arr_hb >= 3:
                st.markdown('<div class="warning-box">Elevated bleeding risk (HAS-BLED >= 3) -- address modifiable factors but do NOT withhold anticoagulation</div>', unsafe_allow_html=True)

        st.subheader("Rate vs. Rhythm Control")
        arr_strategy = st.radio("Management Strategy", ["Rate control", "Rhythm control"], key="arr_strategy")

    if st.button("Run Arrhythmia Workflow", key="btn_arr", type="primary"):
        with st.spinner("Analyzing..."):
            result = api_post("/v1/cardio/workflow/arrhythmia", {
                "arrhythmia_type": arr_type,
                "ventricular_rate": arr_rate, "qtc": arr_qtc,
                "pr_interval": arr_pr, "qrs_duration": arr_qrs,
                "cha2ds2_vasc": arr_cha, "has_bled": arr_hb,
                "management_strategy": arr_strategy,
            })
        if result:
            st.json(result)


# ── Tab 7: Imaging ──

with tab_imaging:
    st.header("Multi-Modal Cardiac Imaging")
    st.markdown("Integrate findings across echocardiography, CT, MRI, and nuclear imaging.")

    img_modality = st.selectbox("Imaging Modality", [
        "Echocardiography", "Cardiac CT", "Cardiac MRI", "Nuclear (SPECT/PET)",
    ], key="img_modality")

    if img_modality == "Echocardiography":
        col_e1, col_e2, col_e3 = st.columns(3)
        with col_e1:
            echo_lvef = st.number_input("LVEF (%)", 5.0, 80.0, 55.0, key="echo_lvef")
            echo_lvidd = st.number_input("LVIDd (cm)", 2.0, 10.0, 4.8, key="echo_lvidd")
            echo_la = st.number_input("LA volume index (mL/m2)", 10.0, 100.0, 28.0, key="echo_la")
        with col_e2:
            echo_gls = st.number_input("GLS (%)", -30.0, 0.0, -19.0, key="echo_gls")
            echo_e_prime = st.number_input("e' septal (cm/s)", 2.0, 20.0, 8.0, key="echo_eprime")
            echo_ee = st.number_input("E/e' ratio", 2.0, 30.0, 10.0, key="echo_ee")
        with col_e3:
            echo_trv = st.number_input("TR velocity (m/s)", 1.0, 5.0, 2.4, key="echo_trv")
            echo_tapse = st.number_input("TAPSE (mm)", 5.0, 35.0, 22.0, key="echo_tapse")
            echo_rwma = st.checkbox("Regional wall motion abnormality", key="echo_rwma")

        if st.button("Analyze Echo", key="btn_echo", type="primary"):
            # Diastolic function grading
            if echo_ee > 14 and echo_la > 34 and echo_trv > 2.8:
                dd_grade = "Grade II-III (elevated filling pressures)"
            elif echo_ee > 14 or echo_la > 34:
                dd_grade = "Grade I-II (indeterminate -- correlate clinically)"
            else:
                dd_grade = "Normal diastolic function"

            st.subheader("Echo Interpretation")
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("LVEF", f"{echo_lvef}%")
            m2.metric("GLS", f"{echo_gls}%")
            m3.metric("E/e'", f"{echo_ee}")
            m4.metric("TAPSE", f"{echo_tapse} mm")
            st.markdown(f"**Diastolic function:** {dd_grade}")
            if echo_rwma:
                st.warning("Regional wall motion abnormality detected -- consider ischemic evaluation.")
            if echo_gls > -16:
                st.warning(f"Reduced GLS ({echo_gls}%) -- subclinical LV dysfunction. Consider cardiotoxicity or early cardiomyopathy.")

    elif img_modality == "Cardiac CT":
        st.info("Enter findings in the CAD Assessment tab for calcium score and CAD-RADS analysis.")

    elif img_modality == "Cardiac MRI":
        col_mri1, col_mri2 = st.columns(2)
        with col_mri1:
            mri_lvef = st.number_input("LVEF (%)", 5.0, 80.0, 50.0, key="mri_lvef")
            mri_rvef = st.number_input("RVEF (%)", 5.0, 80.0, 52.0, key="mri_rvef")
            mri_lge = st.selectbox("LGE Pattern", ["None", "Subendocardial (ischemic)", "Mid-wall (non-ischemic)", "Epicardial", "Diffuse", "Patchy"], key="mri_lge")
        with col_mri2:
            mri_t1 = st.number_input("Native T1 (ms)", 800, 1500, 1050, key="mri_t1")
            mri_t2 = st.number_input("T2 (ms)", 30, 80, 48, key="mri_t2")
            mri_ecv = st.number_input("ECV (%)", 15.0, 60.0, 28.0, key="mri_ecv")

        if st.button("Analyze Cardiac MRI", key="btn_mri", type="primary"):
            with st.spinner("Analyzing..."):
                result = api_post("/v1/cardio/workflow/cardiac-mri", {
                    "lvef": mri_lvef, "rvef": mri_rvef,
                    "lge_pattern": mri_lge, "native_t1": mri_t1,
                    "t2": mri_t2, "ecv": mri_ecv,
                })
            if result:
                st.json(result)
            # Local interpretation
            st.subheader("MRI Findings Summary")
            m1, m2, m3 = st.columns(3)
            m1.metric("LVEF", f"{mri_lvef}%")
            m2.metric("Native T1", f"{mri_t1} ms")
            m3.metric("ECV", f"{mri_ecv}%")
            if mri_lge != "None":
                st.warning(f"LGE pattern: {mri_lge}")
            if mri_t2 > 55:
                st.warning("Elevated T2 -- suggests active inflammation or edema.")
            if mri_ecv > 32:
                st.warning("Elevated ECV -- suggests diffuse fibrosis or infiltrative disease.")

    elif img_modality == "Nuclear (SPECT/PET)":
        col_nuc1, col_nuc2 = st.columns(2)
        with col_nuc1:
            nuc_type = st.selectbox("Study Type", ["SPECT MPI", "PET MPI", "PET Viability", "Pyrophosphate scan"], key="nuc_type")
            nuc_stress = st.selectbox("Stress Protocol", ["Exercise", "Regadenoson", "Adenosine", "Dobutamine"], key="nuc_stress")
            nuc_lvef = st.number_input("Post-stress LVEF (%)", 5.0, 80.0, 55.0, key="nuc_lvef")
        with col_nuc2:
            nuc_sss = st.number_input("Summed Stress Score (SSS)", 0, 56, 4, key="nuc_sss")
            nuc_srs = st.number_input("Summed Rest Score (SRS)", 0, 56, 2, key="nuc_srs")
            nuc_sds = nuc_sss - nuc_srs
            st.metric("Summed Difference Score (SDS)", nuc_sds)

        if st.button("Analyze Nuclear Study", key="btn_nuc", type="primary"):
            if nuc_sds < 2:
                st.markdown('<div class="success-box">Normal perfusion -- low risk for ischemia</div>', unsafe_allow_html=True)
            elif nuc_sds <= 6:
                st.markdown(f'<div class="warning-box">Mild-moderate ischemia (SDS={nuc_sds}) -- consider optimizing medical therapy</div>', unsafe_allow_html=True)
            else:
                st.markdown(f'<div class="warning-box" style="border-left-color:{NVIDIA_THEME["danger"]};">Significant ischemia (SDS={nuc_sds}) -- consider invasive angiography</div>', unsafe_allow_html=True)


# ── Tab 8: Cardio-Oncology ──

with tab_onc:
    st.header("Cardio-Oncology Surveillance")
    st.markdown("Monitor cardiac function during cardiotoxic chemotherapy with GLS tracking and biomarker surveillance.")

    col_onc1, col_onc2 = st.columns(2)

    with col_onc1:
        st.subheader("Chemotherapy Agent")
        onc_agent = st.selectbox("Agent Class", [
            "Anthracycline (doxorubicin, epirubicin)",
            "HER2-targeted (trastuzumab, pertuzumab)",
            "Checkpoint inhibitor (nivolumab, pembrolizumab)",
            "TKI (sunitinib, sorafenib)",
            "Proteasome inhibitor (carfilzomib)",
            "VEGF inhibitor (bevacizumab)",
            "Radiation (mediastinal)",
        ], key="onc_agent")

        onc_cum_dose = st.number_input("Cumulative anthracycline dose (mg/m2)", 0, 1000, 240, key="onc_cum")
        onc_cycles = st.number_input("Treatment cycles completed", 0, 30, 4, key="onc_cycles")

    with col_onc2:
        st.subheader("Cardiac Monitoring")
        onc_lvef_base = st.number_input("Baseline LVEF (%)", 30.0, 80.0, 62.0, key="onc_lvef_base")
        onc_lvef_curr = st.number_input("Current LVEF (%)", 20.0, 80.0, 55.0, key="onc_lvef_curr")
        onc_gls_base = st.number_input("Baseline GLS (%)", -30.0, 0.0, -20.0, key="onc_gls_base")
        onc_gls_curr = st.number_input("Current GLS (%)", -30.0, 0.0, -17.0, key="onc_gls_curr")
        onc_trop = st.number_input("hs-Troponin (ng/L)", 0.0, 1000.0, 8.0, key="onc_trop")
        onc_bnp = st.number_input("NT-proBNP (pg/mL)", 0.0, 10000.0, 125.0, key="onc_bnp")

    if st.button("Assess Cardiotoxicity", key="btn_onc", type="primary"):
        lvef_drop = onc_lvef_base - onc_lvef_curr
        gls_delta = abs(onc_gls_curr) - abs(onc_gls_base)  # Negative = worsening
        gls_rel_change = (gls_delta / abs(onc_gls_base)) * 100 if onc_gls_base != 0 else 0

        st.subheader("Cardiotoxicity Assessment")

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("LVEF Change", f"{-lvef_drop:+.1f}%", delta=f"{-lvef_drop:.1f}%")
        m2.metric("GLS Relative Change", f"{gls_rel_change:+.1f}%", delta=f"{gls_rel_change:.1f}%")
        m3.metric("hs-Troponin", f"{onc_trop} ng/L")
        m4.metric("NT-proBNP", f"{onc_bnp} pg/mL")

        # CTRCD definitions (2022 ESC Cardio-Oncology guidelines)
        if onc_lvef_curr < 40:
            st.markdown(f'<div class="warning-box" style="border-left-color:{NVIDIA_THEME["danger"]};">'
                        f'<strong>Severe CTRCD:</strong> LVEF < 40% -- Hold cardiotoxic therapy. Initiate HF treatment (GDMT). Urgent cardiology referral.</div>',
                        unsafe_allow_html=True)
        elif lvef_drop >= 10 and onc_lvef_curr < 50:
            st.markdown('<div class="warning-box"><strong>Moderate CTRCD:</strong> LVEF drop >= 10% to below 50% -- '
                        'Consider holding therapy. Start cardioprotective agents. Close monitoring.</div>',
                        unsafe_allow_html=True)
        elif gls_rel_change < -15:
            st.markdown('<div class="warning-box"><strong>Subclinical CTRCD:</strong> GLS relative decrease > 15% -- '
                        'Consider cardioprotective agents (ACEi + beta-blocker). Repeat in 2-4 weeks.</div>',
                        unsafe_allow_html=True)
        else:
            st.markdown('<div class="success-box">No significant cardiotoxicity detected. Continue surveillance per protocol.</div>',
                        unsafe_allow_html=True)

        if onc_trop > 14:
            st.warning(f"Elevated hs-Troponin ({onc_trop} ng/L) -- subclinical myocardial injury. Consider cardioprotection.")
        if onc_bnp > 300:
            st.warning(f"Elevated NT-proBNP ({onc_bnp} pg/mL) -- volume overload or myocardial stress.")

        # Cumulative anthracycline risk
        if "anthracycline" in onc_agent.lower():
            if onc_cum_dose > 550:
                st.error(f"Cumulative anthracycline dose ({onc_cum_dose} mg/m2) exceeds safety threshold (550 mg/m2). High risk of irreversible cardiotoxicity.")
            elif onc_cum_dose > 400:
                st.warning(f"Cumulative anthracycline dose ({onc_cum_dose} mg/m2) approaching limit. Consider dexrazoxane or liposomal formulation.")


# ── Tab 9: Evidence Explorer ──

with tab_evidence:
    st.header("Evidence Explorer")
    st.markdown("Search and browse individual cardiology knowledge collections.")

    search_text = st.text_input("Search query", placeholder="e.g., SGLT2 inhibitor heart failure mortality benefit", key="ev_search")

    col_ev1, col_ev2 = st.columns([2, 1])
    with col_ev1:
        ev_collections = st.multiselect(
            "Collections to search",
            all_collections,
            default=selected_collections[:3] if selected_collections else all_collections[:3],
            key="ev_collections",
        )
    with col_ev2:
        ev_top_k = st.slider("Results per collection", 1, 20, 5, key="ev_top_k")

    if st.button("Search Evidence", key="btn_ev_search", type="primary"):
        if search_text.strip():
            with st.spinner("Searching collections..."):
                result = api_post("/v1/cardio/search", {
                    "question": search_text,
                    "collections": ev_collections,
                    "top_k": ev_top_k,
                })
            if result:
                st.metric("Total Results", result.get("total", 0))
                st.markdown(f"**Collections searched:** {', '.join(result.get('collections_searched', []))}")

                for i, res in enumerate(result.get("results", []), 1):
                    with st.expander(f"[{i}] {res.get('collection', 'Unknown')} (score: {res.get('score', 'N/A'):.3f})"):
                        st.markdown(res.get("text", ""))
                        if res.get("metadata"):
                            st.json(res["metadata"])
            else:
                st.info("Search requires API connection.")
        else:
            st.warning("Enter a search query.")

    # Browse reference data
    st.divider()
    st.subheader("Reference Catalogues")
    ref_choice = st.selectbox("Browse", ["Conditions", "Biomarkers", "Drug Classes", "Genes"], key="ev_ref")

    if st.button("Load Reference", key="btn_ev_ref"):
        endpoint_map = {
            "Conditions": "/v1/cardio/conditions",
            "Biomarkers": "/v1/cardio/biomarkers",
            "Drug Classes": "/v1/cardio/drugs",
            "Genes": "/v1/cardio/genes",
        }
        result = api_get(endpoint_map.get(ref_choice, "/v1/cardio/conditions"))
        if result:
            # Extract the main list from the result
            for key in ("conditions", "biomarkers", "drug_classes", "genes"):
                if key in result:
                    st.dataframe(result[key], use_container_width=True)
                    break


# ── Tab 10: Report Generator ──

with tab_reports:
    st.header("Report Generator")
    st.markdown("Generate and export clinical reports in multiple formats.")

    col_rpt1, col_rpt2 = st.columns(2)

    with col_rpt1:
        rpt_type = st.selectbox("Report Type", [
            "clinical_summary",
            "risk_assessment",
            "gdmt_report",
            "workflow_report",
            "imaging_report",
        ], key="rpt_type")

        rpt_format = st.selectbox("Export Format", [
            "markdown", "json", "pdf", "fhir",
        ], key="rpt_format")

        rpt_title = st.text_input("Report Title (optional)", key="rpt_title")
        rpt_patient = st.text_input("Patient ID (optional)", key="rpt_patient")

    with col_rpt2:
        st.subheader("Report Data")
        st.markdown("Paste JSON data from a risk calculation, GDMT optimization, or workflow result.")
        rpt_data_str = st.text_area(
            "Report data (JSON)",
            placeholder='{"calculator": "ASCVD", "score": 12.5, "risk_category": "intermediate", ...}',
            height=200,
            key="rpt_data",
        )

    include_evidence = st.checkbox("Include evidence citations", value=True, key="rpt_evidence")
    include_recs = st.checkbox("Include recommendations", value=True, key="rpt_recs")

    if st.button("Generate Report", key="btn_rpt", type="primary"):
        try:
            rpt_data = json.loads(rpt_data_str) if rpt_data_str.strip() else {}
        except json.JSONDecodeError:
            st.error("Invalid JSON in report data field.")
            rpt_data = None

        if rpt_data is not None:
            with st.spinner("Generating report..."):
                payload = {
                    "report_type": rpt_type,
                    "format": rpt_format,
                    "patient_id": rpt_patient or None,
                    "title": rpt_title or None,
                    "data": rpt_data,
                    "include_evidence": include_evidence,
                    "include_recommendations": include_recs,
                }
                result = api_post("/v1/reports/generate", payload)

            if result:
                st.success(f"Report generated: {result.get('report_id', 'N/A')}")

                col_meta1, col_meta2, col_meta3 = st.columns(3)
                col_meta1.metric("Report ID", result.get("report_id", ""))
                col_meta2.metric("Format", result.get("format", "").upper())
                col_meta3.metric("Generated", result.get("generated_at", "")[:10])

                content = result.get("content", "")
                if rpt_format == "markdown":
                    st.markdown(content)
                elif rpt_format in ("json", "fhir"):
                    try:
                        st.json(json.loads(content))
                    except (json.JSONDecodeError, TypeError):
                        st.code(content, language="json")
                else:
                    st.code(content)

                # Download button
                ext_map = {"markdown": ".md", "json": ".json", "pdf": ".pdf", "fhir": ".json"}
                mime_map = {"markdown": "text/markdown", "json": "application/json", "pdf": "application/pdf", "fhir": "application/fhir+json"}
                st.download_button(
                    "Download Report",
                    data=content,
                    file_name=f"cardio_report_{result.get('report_id', 'unknown')}{ext_map.get(rpt_format, '.txt')}",
                    mime=mime_map.get(rpt_format, "text/plain"),
                )

    # Supported formats reference
    st.divider()
    formats = api_get("/v1/reports/formats")
    if formats:
        st.subheader("Supported Formats")
        for fmt in formats.get("formats", []):
            st.markdown(f"- **{fmt['name']}** ({fmt['extension']}) -- {fmt['description']}")
