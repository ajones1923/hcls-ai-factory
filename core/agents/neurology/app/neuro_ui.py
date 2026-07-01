"""Neurology Intelligence Agent -- 5-Tab Streamlit UI.

NVIDIA dark-themed neurology clinical decision support interface with
RAG-powered queries, clinical scale calculators, workflow runners,
and real-time dashboard monitoring.

Usage:
    streamlit run app/neuro_ui.py --server.port 8529

Author: Adam Jones
Date: March 2026
"""

import json
import os
from datetime import datetime
from typing import Optional

import requests
import streamlit as st

# =====================================================================
# Configuration
# =====================================================================

API_BASE = os.environ.get("NEURO_API_BASE", "http://localhost:8528")

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
    page_title="Neurology Intelligence Agent",
    page_icon="\U0001F9E0",
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
    div[data-testid="stMetric"] [data-testid="stMetricValue"] {{
        color: {NVIDIA_THEME['accent']};
    }}

    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {{
        background-color: {NVIDIA_THEME['bg_secondary']};
        border-radius: 8px;
        padding: 4px;
    }}
    .stTabs [data-baseweb="tab"] {{
        color: {NVIDIA_THEME['text_secondary']};
    }}
    .stTabs [aria-selected="true"] {{
        color: {NVIDIA_THEME['accent']};
        border-bottom-color: {NVIDIA_THEME['accent']};
    }}

    /* Buttons */
    .stButton > button {{
        background-color: {NVIDIA_THEME['accent']};
        color: #000000;
        border: none;
        border-radius: 6px;
        font-weight: 600;
    }}
    .stButton > button:hover {{
        background-color: {NVIDIA_THEME['accent_hover']};
        color: #000000;
    }}

    /* Expanders */
    details {{
        background-color: {NVIDIA_THEME['bg_card']};
        border: 1px solid {NVIDIA_THEME['accent']}40;
        border-radius: 6px;
    }}

    /* Text inputs */
    .stTextInput > div > div > input,
    .stTextArea > div > div > textarea {{
        background-color: {NVIDIA_THEME['bg_secondary']};
        color: {NVIDIA_THEME['text_primary']};
        border: 1px solid {NVIDIA_THEME['accent']}60;
    }}

    /* Select boxes */
    .stSelectbox > div > div {{
        background-color: {NVIDIA_THEME['bg_secondary']};
        color: {NVIDIA_THEME['text_primary']};
    }}

    /* Status indicators */
    .status-healthy {{ color: {NVIDIA_THEME['success']}; font-weight: bold; }}
    .status-degraded {{ color: {NVIDIA_THEME['warning']}; font-weight: bold; }}
    .status-error {{ color: {NVIDIA_THEME['danger']}; font-weight: bold; }}

    /* Agent header */
    .agent-header {{
        background: linear-gradient(135deg, {NVIDIA_THEME['bg_card']}, {NVIDIA_THEME['bg_secondary']});
        border-left: 4px solid {NVIDIA_THEME['accent']};
        padding: 16px 20px;
        border-radius: 0 8px 8px 0;
        margin-bottom: 20px;
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

def api_get(path: str, timeout: int = 15) -> Optional[dict]:
    """GET request to neuro API with error handling."""
    try:
        resp = requests.get(f"{API_BASE}{path}", timeout=timeout)
        resp.raise_for_status()
        return resp.json()
    except requests.exceptions.ConnectionError:
        st.error(f"Cannot connect to API at {API_BASE}. Is the server running?")
        return None
    except requests.exceptions.Timeout:
        st.error(f"API request timed out: {path}")
        return None
    except Exception as exc:
        st.error(f"API error: {exc}")
        return None


def api_post(path: str, data: dict, timeout: int = 60) -> Optional[dict]:
    """POST request to neuro API with error handling."""
    try:
        resp = requests.post(
            f"{API_BASE}{path}",
            json=data,
            timeout=timeout,
            headers={"Content-Type": "application/json"},
        )
        resp.raise_for_status()
        return resp.json()
    except requests.exceptions.ConnectionError:
        st.error(f"Cannot connect to API at {API_BASE}. Is the server running?")
        return None
    except requests.exceptions.Timeout:
        st.error(f"API request timed out: {path}")
        return None
    except requests.exceptions.HTTPError as exc:
        try:
            detail = exc.response.json().get("detail", str(exc))
        except Exception:
            detail = str(exc)
        st.error(f"API error ({exc.response.status_code}): {detail}")
        return None
    except Exception as exc:
        st.error(f"API error: {exc}")
        return None


# =====================================================================
# Sidebar
# =====================================================================

with st.sidebar:
    st.markdown(f"""
    <div class="agent-header">
        <h2 style="color: {NVIDIA_THEME['accent']}; margin: 0;">Neurology Intelligence</h2>
        <p style="color: {NVIDIA_THEME['text_secondary']}; margin: 4px 0 0 0; font-size: 0.85em;">
            HCLS AI Factory Agent
        </p>
    </div>
    """, unsafe_allow_html=True)

    # Health status
    health = api_get("/health")
    if health:
        status = health.get("status", "unknown")
        status_class = "status-healthy" if status == "healthy" else "status-degraded"
        st.markdown(f'<p class="{status_class}">Status: {status.upper()}</p>', unsafe_allow_html=True)

        components = health.get("components", {})
        for comp, state in components.items():
            icon = "+" if state in ("connected", "ready") else "-"
            st.text(f"  {icon} {comp}: {state}")

        st.markdown("---")
        st.metric("Collections", health.get("collections", 0))
        st.metric("Vectors", f"{health.get('total_vectors', 0):,}")
        st.metric("Workflows", health.get("workflows", 0))
        st.metric("Scales", health.get("scales", 0))
    else:
        st.warning("API unavailable")

    st.markdown("---")
    st.caption(f"API: {API_BASE}")
    st.caption(f"v1.0.0 | {datetime.now().strftime('%Y-%m-%d')}")


# =====================================================================
# Main Content - Tabs
# =====================================================================

tab_dashboard, tab_explorer, tab_scales, tab_workflows, tab_reports = st.tabs([
    "Dashboard",
    "Evidence Explorer",
    "Clinical Scales",
    "Workflow Runner",
    "Reports & Export",
])


# =====================================================================
# Tab 1: Dashboard
# =====================================================================

with tab_dashboard:
    st.header("Neurology Intelligence Dashboard")

    # Health overview
    col1, col2, col3, col4 = st.columns(4)

    if health:
        with col1:
            st.metric("Service Status", health.get("status", "unknown").upper())
        with col2:
            st.metric("Collections", health.get("collections", 0))
        with col3:
            st.metric("Total Vectors", f"{health.get('total_vectors', 0):,}")
        with col4:
            st.metric("Clinical Scales", health.get("scales", 10))
    else:
        st.info("Connect to the API to view dashboard metrics.")

    st.markdown("---")

    # Domains overview
    st.subheader("Neurology Domains")
    domains = api_get("/v1/neuro/domains")
    if domains:
        domain_list = domains.get("domains", [])
        cols = st.columns(2)
        for i, domain in enumerate(domain_list):
            with cols[i % 2]:
                with st.expander(f"{domain.get('name', 'Unknown')}", expanded=False):
                    st.write(domain.get("description", "No description available."))
                    st.caption(f"ID: {domain.get('id', 'N/A')}")

    # Scales overview
    st.subheader("Available Clinical Scales")
    scales = api_get("/v1/neuro/scales")
    if scales:
        scale_list = scales.get("scales", [])
        cols = st.columns(3)
        for i, scale in enumerate(scale_list):
            with cols[i % 3]:
                st.text(f"  {scale.get('name', 'N/A')} (0-{scale.get('max_score', '?')})")

    # Metrics
    st.subheader("Service Metrics")
    try:
        resp = requests.get(f"{API_BASE}/metrics", timeout=10)
        if resp.status_code == 200:
            st.code(resp.text, language="text")
    except Exception:
        st.info("Metrics unavailable.")


# =====================================================================
# Tab 2: Evidence Explorer (RAG Q&A)
# =====================================================================

with tab_explorer:
    st.header("Evidence Explorer")
    st.write("RAG-powered neurology Q&A across all knowledge collections.")

    # Domain selector
    domain_options = [
        "auto", "stroke", "dementia", "epilepsy", "tumors", "ms",
        "parkinsons", "headache", "neuromuscular", "general",
    ]
    selected_domain = st.selectbox(
        "Domain Focus",
        domain_options,
        index=0,
        help="Select a neurology domain to guide the query, or leave as 'auto'.",
    )

    # Query input
    question = st.text_area(
        "Neurology Question",
        placeholder="e.g., What are the 2024 updates to mechanical thrombectomy eligibility criteria for large vessel occlusion stroke?",
        height=100,
    )

    col_topk, col_guidelines = st.columns(2)
    with col_topk:
        top_k = st.slider("Evidence passages (top_k)", 1, 20, 5)
    with col_guidelines:
        include_guidelines = st.checkbox("Include guideline citations", value=True)

    if st.button("Search", key="explorer_search"):
        if question.strip():
            with st.spinner("Searching neurology knowledge base..."):
                payload = {
                    "question": question.strip(),
                    "top_k": top_k,
                    "include_guidelines": include_guidelines,
                }
                if selected_domain != "auto":
                    payload["domain"] = selected_domain

                result = api_post("/v1/neuro/query", payload)

            if result:
                st.subheader("Answer")
                st.markdown(result.get("answer", "No answer generated."))

                if result.get("guidelines_cited"):
                    st.subheader("Guidelines Cited")
                    for g in result["guidelines_cited"]:
                        st.write(f"- {g}")

                confidence = result.get("confidence", 0)
                st.progress(confidence, text=f"Confidence: {confidence:.0%}")

                evidence = result.get("evidence", [])
                if evidence:
                    st.subheader(f"Evidence ({len(evidence)} passages)")
                    for i, ev in enumerate(evidence):
                        with st.expander(f"[{ev.get('collection', 'unknown')}] Score: {ev.get('score', 0):.3f}"):
                            st.write(ev.get("text", ""))
                            if ev.get("metadata"):
                                st.json(ev["metadata"])
        else:
            st.warning("Please enter a question.")


# =====================================================================
# Tab 3: Clinical Scales
# =====================================================================

with tab_scales:
    st.header("Clinical Scale Calculators")
    st.write("Validated neurological assessment scales with automated scoring and interpretation.")

    scale_choice = st.selectbox(
        "Select Scale",
        ["NIHSS (Stroke)", "GCS (Consciousness)", "MoCA (Cognition)",
         "UPDRS-III (Parkinson's Motor)", "EDSS (MS Disability)",
         "mRS (Rankin Outcome)", "HIT-6 (Headache Impact)",
         "ALSFRS-R (ALS Function)", "ASPECTS (Stroke CT)",
         "Hoehn-Yahr (PD Staging)"],
        key="scale_selector",
    )

    scale_map = {
        "NIHSS (Stroke)": "nihss",
        "GCS (Consciousness)": "gcs",
        "MoCA (Cognition)": "moca",
        "UPDRS-III (Parkinson's Motor)": "updrs",
        "EDSS (MS Disability)": "edss",
        "mRS (Rankin Outcome)": "mrs",
        "HIT-6 (Headache Impact)": "hit6",
        "ALSFRS-R (ALS Function)": "alsfrs",
        "ASPECTS (Stroke CT)": "aspects",
        "Hoehn-Yahr (PD Staging)": "hoehn_yahr",
    }
    scale_id = scale_map.get(scale_choice, "nihss")

    items = {}

    if scale_id == "nihss":
        st.subheader("NIH Stroke Scale (0-42)")
        c1, c2, c3 = st.columns(3)
        with c1:
            items["1a_loc"] = st.selectbox("1a. LOC", [0, 1, 2, 3], key="nihss_1a")
            items["1b_loc_questions"] = st.selectbox("1b. LOC Questions", [0, 1, 2], key="nihss_1b")
            items["1c_loc_commands"] = st.selectbox("1c. LOC Commands", [0, 1, 2], key="nihss_1c")
            items["2_gaze"] = st.selectbox("2. Best Gaze", [0, 1, 2], key="nihss_2")
            items["3_visual"] = st.selectbox("3. Visual Fields", [0, 1, 2, 3], key="nihss_3")
        with c2:
            items["4_facial_palsy"] = st.selectbox("4. Facial Palsy", [0, 1, 2, 3], key="nihss_4")
            items["5a_motor_left_arm"] = st.selectbox("5a. L Arm Motor", [0, 1, 2, 3, 4], key="nihss_5a")
            items["5b_motor_right_arm"] = st.selectbox("5b. R Arm Motor", [0, 1, 2, 3, 4], key="nihss_5b")
            items["6a_motor_left_leg"] = st.selectbox("6a. L Leg Motor", [0, 1, 2, 3, 4], key="nihss_6a")
            items["6b_motor_right_leg"] = st.selectbox("6b. R Leg Motor", [0, 1, 2, 3, 4], key="nihss_6b")
        with c3:
            items["7_limb_ataxia"] = st.selectbox("7. Limb Ataxia", [0, 1, 2], key="nihss_7")
            items["8_sensory"] = st.selectbox("8. Sensory", [0, 1, 2], key="nihss_8")
            items["9_language"] = st.selectbox("9. Best Language", [0, 1, 2, 3], key="nihss_9")
            items["10_dysarthria"] = st.selectbox("10. Dysarthria", [0, 1, 2], key="nihss_10")
            items["11_extinction"] = st.selectbox("11. Extinction/Inattention", [0, 1, 2], key="nihss_11")

    elif scale_id == "gcs":
        st.subheader("Glasgow Coma Scale (3-15)")
        c1, c2, c3 = st.columns(3)
        with c1:
            items["eye"] = st.selectbox("Eye Opening", [1, 2, 3, 4], format_func=lambda x: {1: "1 - None", 2: "2 - To pain", 3: "3 - To voice", 4: "4 - Spontaneous"}[x], key="gcs_e")
        with c2:
            items["verbal"] = st.selectbox("Verbal Response", [1, 2, 3, 4, 5], format_func=lambda x: {1: "1 - None", 2: "2 - Incomprehensible", 3: "3 - Inappropriate", 4: "4 - Confused", 5: "5 - Oriented"}[x], key="gcs_v")
        with c3:
            items["motor"] = st.selectbox("Motor Response", [1, 2, 3, 4, 5, 6], format_func=lambda x: {1: "1 - None", 2: "2 - Extension", 3: "3 - Flexion", 4: "4 - Withdrawal", 5: "5 - Localizing", 6: "6 - Obeys"}[x], key="gcs_m")

    elif scale_id == "moca":
        st.subheader("Montreal Cognitive Assessment (0-30)")
        c1, c2 = st.columns(2)
        with c1:
            items["visuospatial"] = st.number_input("Visuospatial/Executive (0-5)", 0, 5, 5, key="moca_vs")
            items["naming"] = st.number_input("Naming (0-3)", 0, 3, 3, key="moca_nm")
            items["attention"] = st.number_input("Attention (0-6)", 0, 6, 6, key="moca_at")
            items["language"] = st.number_input("Language (0-3)", 0, 3, 3, key="moca_lg")
        with c2:
            items["abstraction"] = st.number_input("Abstraction (0-2)", 0, 2, 2, key="moca_ab")
            items["delayed_recall"] = st.number_input("Delayed Recall (0-5)", 0, 5, 5, key="moca_dr")
            items["orientation"] = st.number_input("Orientation (0-6)", 0, 6, 6, key="moca_or")
            items["education_years"] = st.number_input("Education (years)", 0, 30, 16, key="moca_ed")

    elif scale_id == "updrs":
        st.subheader("MDS-UPDRS Part III Motor (0-132)")
        st.write("Enter item scores (each 0-4). Simplified to key domains:")
        c1, c2 = st.columns(2)
        with c1:
            items["speech"] = st.selectbox("3.1 Speech", [0, 1, 2, 3, 4], key="updrs_speech")
            items["facial_expression"] = st.selectbox("3.2 Facial Expression", [0, 1, 2, 3, 4], key="updrs_face")
            items["rigidity_neck"] = st.selectbox("3.3a Rigidity - Neck", [0, 1, 2, 3, 4], key="updrs_rig_n")
            items["rigidity_rue"] = st.selectbox("3.3b Rigidity - RUE", [0, 1, 2, 3, 4], key="updrs_rig_rue")
            items["rigidity_lue"] = st.selectbox("3.3c Rigidity - LUE", [0, 1, 2, 3, 4], key="updrs_rig_lue")
            items["finger_tapping_r"] = st.selectbox("3.4a Finger Tap - R", [0, 1, 2, 3, 4], key="updrs_ft_r")
            items["finger_tapping_l"] = st.selectbox("3.4b Finger Tap - L", [0, 1, 2, 3, 4], key="updrs_ft_l")
        with c2:
            items["hand_movements_r"] = st.selectbox("3.5a Hand Mvt - R", [0, 1, 2, 3, 4], key="updrs_hm_r")
            items["hand_movements_l"] = st.selectbox("3.5b Hand Mvt - L", [0, 1, 2, 3, 4], key="updrs_hm_l")
            items["3.10"] = st.selectbox("3.10 Gait", [0, 1, 2, 3, 4], key="updrs_gait")
            items["3.12"] = st.selectbox("3.12 Postural Stability", [0, 1, 2, 3, 4], key="updrs_post")
            items["3.15a"] = st.selectbox("3.15a Postural Tremor - R", [0, 1, 2, 3, 4], key="updrs_pt_r")
            items["3.15b"] = st.selectbox("3.15b Postural Tremor - L", [0, 1, 2, 3, 4], key="updrs_pt_l")
            items["3.17a"] = st.selectbox("3.17a Rest Tremor - RUE", [0, 1, 2, 3, 4], key="updrs_rt_r")

    elif scale_id == "edss":
        st.subheader("Expanded Disability Status Scale (0-10)")
        c1, c2 = st.columns(2)
        with c1:
            items["pyramidal"] = st.selectbox("Pyramidal FS", [0, 1, 2, 3, 4, 5, 6], key="edss_pyr")
            items["cerebellar"] = st.selectbox("Cerebellar FS", [0, 1, 2, 3, 4, 5], key="edss_cer")
            items["brainstem"] = st.selectbox("Brainstem FS", [0, 1, 2, 3, 4, 5], key="edss_bs")
            items["sensory"] = st.selectbox("Sensory FS", [0, 1, 2, 3, 4, 5, 6], key="edss_sen")
        with c2:
            items["bowel_bladder"] = st.selectbox("Bowel/Bladder FS", [0, 1, 2, 3, 4, 5, 6], key="edss_bb")
            items["visual"] = st.selectbox("Visual FS", [0, 1, 2, 3, 4, 5, 6], key="edss_vis")
            items["cerebral"] = st.selectbox("Cerebral FS", [0, 1, 2, 3, 4, 5], key="edss_cer2")
            items["ambulation"] = st.number_input("Ambulation (EDSS 0-10, 0=use FS only)", 0.0, 10.0, 0.0, 0.5, key="edss_amb")

    elif scale_id == "mrs":
        st.subheader("Modified Rankin Scale (0-6)")
        items["score"] = st.selectbox(
            "mRS Score", [0, 1, 2, 3, 4, 5, 6],
            format_func=lambda x: {
                0: "0 - No symptoms",
                1: "1 - No significant disability",
                2: "2 - Slight disability",
                3: "3 - Moderate disability",
                4: "4 - Moderately severe disability",
                5: "5 - Severe disability",
                6: "6 - Dead",
            }[x],
            key="mrs_score",
        )

    elif scale_id == "hit6":
        st.subheader("Headache Impact Test-6 (36-78)")
        st.write("Score each item: 6=Never, 8=Rarely, 10=Sometimes, 11=Very Often, 13=Always")
        options = [6, 8, 10, 11, 13]
        labels = {6: "Never", 8: "Rarely", 10: "Sometimes", 11: "Very Often", 13: "Always"}
        c1, c2 = st.columns(2)
        with c1:
            items["pain_severity"] = st.selectbox("Pain severity limits activities", options, format_func=lambda x: f"{x} - {labels[x]}", key="hit6_1")
            items["daily_activities"] = st.selectbox("Limits daily activities", options, format_func=lambda x: f"{x} - {labels[x]}", key="hit6_2")
            items["lie_down"] = st.selectbox("Wish to lie down", options, format_func=lambda x: f"{x} - {labels[x]}", key="hit6_3")
        with c2:
            items["too_tired"] = st.selectbox("Too tired for work/daily activities", options, format_func=lambda x: f"{x} - {labels[x]}", key="hit6_4")
            items["fed_up"] = st.selectbox("Fed up or irritated", options, format_func=lambda x: f"{x} - {labels[x]}", key="hit6_5")
            items["concentration"] = st.selectbox("Difficulty concentrating", options, format_func=lambda x: f"{x} - {labels[x]}", key="hit6_6")

    elif scale_id == "alsfrs":
        st.subheader("ALS Functional Rating Scale - Revised (0-48)")
        st.write("Score each item 0-4 (4 = normal function)")
        c1, c2, c3 = st.columns(3)
        with c1:
            items["speech"] = st.selectbox("Speech", [4, 3, 2, 1, 0], key="als_speech")
            items["salivation"] = st.selectbox("Salivation", [4, 3, 2, 1, 0], key="als_saliva")
            items["swallowing"] = st.selectbox("Swallowing", [4, 3, 2, 1, 0], key="als_swallow")
            items["handwriting"] = st.selectbox("Handwriting", [4, 3, 2, 1, 0], key="als_write")
        with c2:
            items["cutting_food"] = st.selectbox("Cutting Food", [4, 3, 2, 1, 0], key="als_cut")
            items["dressing"] = st.selectbox("Dressing/Hygiene", [4, 3, 2, 1, 0], key="als_dress")
            items["turning_in_bed"] = st.selectbox("Turning in Bed", [4, 3, 2, 1, 0], key="als_turn")
            items["walking"] = st.selectbox("Walking", [4, 3, 2, 1, 0], key="als_walk")
        with c3:
            items["climbing_stairs"] = st.selectbox("Climbing Stairs", [4, 3, 2, 1, 0], key="als_stairs")
            items["dyspnea"] = st.selectbox("Dyspnea", [4, 3, 2, 1, 0], key="als_dysp")
            items["orthopnea"] = st.selectbox("Orthopnea", [4, 3, 2, 1, 0], key="als_ortho")
            items["respiratory_insufficiency"] = st.selectbox("Respiratory Insufficiency", [4, 3, 2, 1, 0], key="als_resp")

    elif scale_id == "aspects":
        st.subheader("ASPECTS (0-10)")
        st.write("Mark each region: 1 = Normal, 0 = Early ischemic change")
        c1, c2 = st.columns(2)
        regions = ["C (caudate)", "L (lentiform)", "IC (internal capsule)", "I (insular ribbon)", "M1", "M2", "M3", "M4", "M5", "M6"]
        region_keys = ["C", "L", "IC", "I", "M1", "M2", "M3", "M4", "M5", "M6"]
        for i, (name, key) in enumerate(zip(regions, region_keys)):
            with c1 if i < 5 else c2:
                items[key] = st.selectbox(name, [1, 0], format_func=lambda x: "Normal" if x == 1 else "Ischemic change", key=f"aspects_{key}")

    elif scale_id == "hoehn_yahr":
        st.subheader("Hoehn and Yahr Scale (0-5)")
        items["stage"] = st.selectbox(
            "Stage", [0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0],
            format_func=lambda x: {
                0.0: "0 - No signs",
                1.0: "1 - Unilateral only",
                1.5: "1.5 - Unilateral + axial",
                2.0: "2 - Bilateral, no balance impairment",
                2.5: "2.5 - Mild bilateral, recovery on pull",
                3.0: "3 - Mild-moderate bilateral, postural instability",
                4.0: "4 - Severe, still ambulatory",
                5.0: "5 - Wheelchair/bedridden",
            }[x],
            key="hy_stage",
        )

    if st.button("Calculate Score", key="scale_calculate"):
        with st.spinner("Calculating..."):
            payload = {
                "scale_name": scale_id,
                "items": {k: v for k, v in items.items()},
            }
            result = api_post("/v1/neuro/scale/calculate", payload)

        if result:
            st.subheader("Results")
            rc1, rc2, rc3 = st.columns(3)
            with rc1:
                st.metric("Score", f"{result.get('total_score', 0)} / {result.get('max_score', 0)}")
            with rc2:
                st.metric("Severity", result.get("severity_category", "N/A").replace("_", " ").title())
            with rc3:
                st.metric("Items Scored", result.get("items_received", 0))

            st.write(f"**Interpretation:** {result.get('interpretation', 'N/A')}")

            recs = result.get("recommendations", [])
            if recs:
                st.subheader("Recommendations")
                for rec in recs:
                    st.write(f"- {rec}")

            # Export button
            if st.button("Export Scale Report", key="scale_export"):
                report_data = {
                    "report_type": "scale_summary",
                    "format": "markdown",
                    "title": f"{scale_choice} Assessment",
                    "data": result,
                }
                report = api_post("/v1/reports/generate", report_data)
                if report:
                    st.download_button(
                        "Download Report",
                        data=report.get("content", ""),
                        file_name=f"neuro_scale_{scale_id}_{report.get('report_id', 'report')}.md",
                        mime="text/markdown",
                    )


# =====================================================================
# Tab 4: Workflow Runner
# =====================================================================

with tab_workflows:
    st.header("Neurology Workflow Runner")
    st.write("Execute specialized neurology clinical workflows.")

    workflow_choice = st.selectbox(
        "Select Workflow",
        [
            "Acute Stroke Triage",
            "Dementia Evaluation",
            "Epilepsy Classification",
            "Brain Tumor Grading",
            "MS Assessment",
            "Parkinson's Assessment",
            "Headache Classification",
            "Neuromuscular Evaluation",
        ],
        key="workflow_selector",
    )

    workflow_endpoints = {
        "Acute Stroke Triage": "/v1/neuro/stroke/triage",
        "Dementia Evaluation": "/v1/neuro/dementia/evaluate",
        "Epilepsy Classification": "/v1/neuro/epilepsy/classify",
        "Brain Tumor Grading": "/v1/neuro/tumor/grade",
        "MS Assessment": "/v1/neuro/ms/assess",
        "Parkinson's Assessment": "/v1/neuro/parkinsons/assess",
        "Headache Classification": "/v1/neuro/headache/classify",
        "Neuromuscular Evaluation": "/v1/neuro/neuromuscular/evaluate",
    }

    wf_payload = {}

    if workflow_choice == "Acute Stroke Triage":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["nihss_total"] = st.number_input("NIHSS Total (0-42)", 0, 42, 0, key="wf_nihss")
            wf_payload["onset_time_hours"] = st.number_input("Hours since onset", 0.0, 48.0, 2.0, 0.5, key="wf_onset")
            wf_payload["ct_aspects_score"] = st.number_input("ASPECTS Score (0-10)", 0, 10, 10, key="wf_aspects")
        with c2:
            wf_payload["lvo_suspected"] = st.checkbox("LVO Suspected", key="wf_lvo")
            wf_payload["age"] = st.number_input("Age", 0, 120, 65, key="wf_age")
            wf_payload["anticoagulant_use"] = st.checkbox("Anticoagulant Use", key="wf_ac")
        wf_payload["clinical_notes"] = st.text_area("Clinical Notes (optional)", key="wf_stroke_notes")

    elif workflow_choice == "Dementia Evaluation":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["moca_total"] = st.number_input("MoCA Total (0-30)", 0, 30, 26, key="wf_moca")
            wf_payload["age"] = st.number_input("Age", 0, 120, 72, key="wf_dem_age")
            wf_payload["education_years"] = st.number_input("Education (years)", 0, 30, 12, key="wf_edu")
        with c2:
            symptoms_text = st.text_area("Dominant Symptoms (one per line)", placeholder="memory_loss\nbehavioral_change", key="wf_dem_sym")
            wf_payload["dominant_symptoms"] = [s.strip() for s in symptoms_text.split("\n") if s.strip()]
            wf_payload["motor_features"] = st.text_input("Motor Features", key="wf_dem_motor")
        wf_payload["clinical_notes"] = st.text_area("Clinical Notes (optional)", key="wf_dem_notes")

    elif workflow_choice == "Epilepsy Classification":
        wf_payload["seizure_description"] = st.text_area("Seizure Description", placeholder="Describe the seizure semiology...", key="wf_epi_desc")
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["age_at_onset"] = st.number_input("Age at Onset", 0, 120, 20, key="wf_epi_age")
            wf_payload["eeg_findings"] = st.text_input("EEG Findings", key="wf_epi_eeg")
        with c2:
            wf_payload["mri_findings"] = st.text_input("MRI Findings", key="wf_epi_mri")
            aeds = st.text_input("Current AEDs (comma-separated)", key="wf_epi_aeds")
            wf_payload["current_aeds"] = [a.strip() for a in aeds.split(",") if a.strip()]

    elif workflow_choice == "Brain Tumor Grading":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["histology"] = st.text_input("Histological Diagnosis", placeholder="e.g., Diffuse astrocytoma", key="wf_tumor_hist")
            markers = st.text_input("Molecular Markers (comma-separated)", placeholder="idh_mutant, mgmt_methylated", key="wf_tumor_mol")
            wf_payload["molecular_markers"] = [m.strip() for m in markers.split(",") if m.strip()]
        with c2:
            wf_payload["location"] = st.text_input("Tumor Location", key="wf_tumor_loc")
            wf_payload["kps"] = st.number_input("KPS (0-100)", 0, 100, 80, key="wf_tumor_kps")
            wf_payload["extent_of_resection"] = st.selectbox("Extent of Resection", ["GTR", "STR", "biopsy"], key="wf_tumor_eor")

    elif workflow_choice == "MS Assessment":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["ms_course"] = st.selectbox("MS Course", ["rrms", "spms", "ppms", "cis", "ris"], key="wf_ms_course")
            wf_payload["edss_score"] = st.number_input("EDSS Score (0-10)", 0.0, 10.0, 2.0, 0.5, key="wf_ms_edss")
            wf_payload["relapse_count_2yr"] = st.number_input("Relapses (past 2yr)", 0, 20, 0, key="wf_ms_relapse")
        with c2:
            wf_payload["new_t2_lesions"] = st.number_input("New T2 Lesions", 0, 50, 0, key="wf_ms_t2")
            wf_payload["gad_enhancing_lesions"] = st.number_input("Gd+ Lesions", 0, 20, 0, key="wf_ms_gad")
            wf_payload["current_dmt"] = st.text_input("Current DMT", key="wf_ms_dmt")
            wf_payload["jcv_status"] = st.selectbox("JCV Status", ["unknown", "positive", "negative"], key="wf_ms_jcv")

    elif workflow_choice == "Parkinson's Assessment":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["updrs_total"] = st.number_input("UPDRS-III Total (0-132)", 0, 132, 20, key="wf_pd_updrs")
            wf_payload["hoehn_yahr"] = st.selectbox("Hoehn-Yahr Stage", [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0], key="wf_pd_hy")
            wf_payload["disease_duration_years"] = st.number_input("Disease Duration (years)", 0.0, 40.0, 5.0, key="wf_pd_dur")
        with c2:
            wf_payload["motor_fluctuations"] = st.checkbox("Motor Fluctuations", key="wf_pd_fluct")
            wf_payload["dyskinesia"] = st.checkbox("Dyskinesia", key="wf_pd_dysk")
            meds = st.text_input("Current Medications (comma-separated)", key="wf_pd_meds")
            wf_payload["current_medications"] = [m.strip() for m in meds.split(",") if m.strip()]
            nms = st.text_input("Non-motor Symptoms (comma-separated)", key="wf_pd_nms")
            wf_payload["non_motor_symptoms"] = [s.strip() for s in nms.split(",") if s.strip()]

    elif workflow_choice == "Headache Classification":
        wf_payload["headache_description"] = st.text_area("Headache Description", placeholder="Describe headache characteristics...", key="wf_ha_desc")
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["duration_hours"] = st.number_input("Duration (hours)", 0.0, 168.0, 8.0, key="wf_ha_dur")
            wf_payload["frequency_per_month"] = st.number_input("Frequency/month", 0, 30, 4, key="wf_ha_freq")
            wf_payload["location"] = st.selectbox("Location", ["unilateral", "bilateral", "occipital", "diffuse"], key="wf_ha_loc")
        with c2:
            wf_payload["quality"] = st.selectbox("Quality", ["throbbing", "pressure", "stabbing", "dull"], key="wf_ha_qual")
            wf_payload["aura"] = st.checkbox("Aura present", key="wf_ha_aura")
            symptoms = st.text_input("Associated symptoms (comma-separated)", placeholder="nausea, photophobia", key="wf_ha_sym")
            wf_payload["associated_symptoms"] = [s.strip() for s in symptoms.split(",") if s.strip()]

    elif workflow_choice == "Neuromuscular Evaluation":
        wf_payload["presentation"] = st.text_area("Clinical Presentation", placeholder="Describe weakness pattern and symptoms...", key="wf_nm_pres")
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["weakness_pattern"] = st.selectbox("Weakness Pattern", ["proximal", "distal", "bulbar", "respiratory", "diffuse"], key="wf_nm_weak")
            wf_payload["sensory_involvement"] = st.checkbox("Sensory Involvement", key="wf_nm_sens")
            wf_payload["emg_findings"] = st.text_input("EMG/NCS Findings", key="wf_nm_emg")
        with c2:
            wf_payload["ck_level"] = st.number_input("CK Level (IU/L)", 0.0, 50000.0, 100.0, key="wf_nm_ck")
            abs_text = st.text_input("Antibodies (comma-separated)", key="wf_nm_abs")
            wf_payload["antibodies"] = [a.strip() for a in abs_text.split(",") if a.strip()]
            wf_payload["progression_rate"] = st.selectbox("Progression Rate", ["slow", "moderate", "rapid"], key="wf_nm_prog")

    if st.button("Run Workflow", key="wf_run"):
        endpoint = workflow_endpoints.get(workflow_choice, "")
        if endpoint:
            # Clean payload: remove empty strings and None values for optional fields
            clean_payload = {k: v for k, v in wf_payload.items() if v is not None and v != "" and v != []}
            with st.spinner(f"Running {workflow_choice}..."):
                result = api_post(endpoint, clean_payload)

            if result:
                st.subheader("Results")
                st.json(result)

                # Recommendations
                recs = result.get("recommendations", [])
                if recs:
                    st.subheader("Recommendations")
                    for rec in recs:
                        st.write(f"- {rec}")

                # Guidelines
                guidelines = result.get("guidelines_cited", [])
                if guidelines:
                    st.subheader("Guidelines Cited")
                    for gl in guidelines:
                        st.write(f"- {gl}")


# =====================================================================
# Tab 5: Reports & Export
# =====================================================================

with tab_reports:
    st.header("Reports & Export")
    st.write("Generate and export structured neurology reports.")

    report_type = st.selectbox(
        "Report Type",
        [
            "stroke_triage", "cognitive_assessment", "epilepsy_classification",
            "tumor_grading", "ms_assessment", "parkinsons_assessment",
            "headache_report", "neuromuscular_evaluation", "scale_summary",
            "general",
        ],
        key="report_type",
    )

    export_format = st.selectbox(
        "Export Format",
        ["markdown", "json", "fhir", "pdf"],
        key="report_format",
    )

    report_title = st.text_input("Report Title (optional)", key="report_title")
    patient_id = st.text_input("Patient ID (optional)", key="report_patient")
    encounter_id = st.text_input("Encounter ID (optional)", key="report_encounter")

    report_data_raw = st.text_area(
        "Report Data (JSON)",
        value='{\n  "summary": "Clinical assessment results"\n}',
        height=150,
        key="report_data",
    )

    if st.button("Generate Report", key="report_generate"):
        try:
            data_dict = json.loads(report_data_raw)
        except json.JSONDecodeError:
            st.error("Invalid JSON in Report Data field.")
            data_dict = None

        if data_dict is not None:
            with st.spinner("Generating report..."):
                payload = {
                    "report_type": report_type,
                    "format": export_format,
                    "data": data_dict,
                    "include_evidence": True,
                    "include_recommendations": True,
                }
                if report_title:
                    payload["title"] = report_title
                if patient_id:
                    payload["patient_id"] = patient_id
                if encounter_id:
                    payload["encounter_id"] = encounter_id

                result = api_post("/v1/reports/generate", payload)

            if result:
                st.subheader(result.get("title", "Report"))
                st.caption(f"Report ID: {result.get('report_id', 'N/A')} | Generated: {result.get('generated_at', 'N/A')}")

                content = result.get("content", "")
                if export_format == "markdown":
                    st.markdown(content)
                elif export_format in ("json", "fhir"):
                    st.json(json.loads(content) if isinstance(content, str) else content)
                else:
                    st.code(content)

                # Download button
                ext = {"markdown": ".md", "json": ".json", "fhir": ".json", "pdf": ".pdf"}.get(export_format, ".txt")
                mime = {"markdown": "text/markdown", "json": "application/json", "fhir": "application/fhir+json", "pdf": "application/pdf"}.get(export_format, "text/plain")
                st.download_button(
                    "Download Report",
                    data=content,
                    file_name=f"neuro_{report_type}_{result.get('report_id', 'report')}{ext}",
                    mime=mime,
                )

    # Available formats reference
    st.markdown("---")
    st.subheader("Supported Formats")
    formats = api_get("/v1/reports/formats")
    if formats:
        for fmt in formats.get("formats", []):
            st.write(f"- **{fmt.get('name', 'N/A')}** ({fmt.get('extension', '')}) -- {fmt.get('description', '')}")
