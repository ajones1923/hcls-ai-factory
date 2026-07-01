"""Clinical Trial Intelligence Agent -- 5-Tab Streamlit UI.

NVIDIA dark-themed clinical trial decision support interface with
RAG-powered queries, patient-trial matching, protocol optimization,
competitive landscape analysis, and real-time dashboard monitoring.

Usage:
    streamlit run app/trial_ui.py --server.port 8128

Author: Adam Jones
Date: March 2026
"""

import os
from datetime import datetime
from typing import Optional

import requests
import streamlit as st

# =====================================================================
# Configuration
# =====================================================================

API_BASE = os.environ.get("TRIAL_API_BASE", "http://localhost:8538")

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
    page_title="Clinical Trial Intelligence Agent",
    page_icon="\\U0001F9EA",
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
    """GET request to trial API with error handling."""
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
    """POST request to trial API with error handling."""
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
        <h2 style="color: {NVIDIA_THEME['accent']}; margin: 0;">Clinical Trial Intelligence</h2>
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
    else:
        st.warning("API unavailable")

    st.markdown("---")
    st.caption(f"API: {API_BASE}")
    st.caption(f"v1.0.0 | {datetime.now().strftime('%Y-%m-%d')}")


# =====================================================================
# Main Content - Tabs
# =====================================================================

tab_dashboard, tab_explorer, tab_matcher, tab_protocol, tab_competitive = st.tabs([
    "Dashboard",
    "Trial Explorer",
    "Patient Matcher",
    "Protocol Analyzer",
    "Competitive Dashboard",
])


# =====================================================================
# Tab 1: Dashboard
# =====================================================================

with tab_dashboard:
    st.header("Clinical Trial Intelligence Dashboard")

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
            st.metric("Workflows", health.get("workflows", 11))
    else:
        st.info("Connect to the API to view dashboard metrics.")

    st.markdown("---")

    # Workflows overview
    st.subheader("Available Workflows")
    workflows = api_get("/workflows")
    if workflows:
        wf_list = workflows.get("workflows", [])
        cols = st.columns(2)
        for i, wf in enumerate(wf_list):
            with cols[i % 2]:
                with st.expander(f"{wf.get('name', wf.get('id', 'Unknown'))}", expanded=False):
                    st.write(wf.get("description", "No description available."))
                    st.caption(f"ID: {wf.get('id', 'N/A')}")

    # Collections overview
    st.subheader("Knowledge Collections")
    collections = api_get("/collections")
    if collections:
        coll_list = collections.get("collections", [])
        if coll_list:
            cols = st.columns(3)
            for i, cname in enumerate(coll_list):
                with cols[i % 3]:
                    st.text(f"  {cname}")
        else:
            st.info("No collections loaded.")

    # Metrics
    st.subheader("Service Metrics")
    try:
        resp = requests.get(f"{API_BASE}/metrics", timeout=10)
        if resp.status_code == 200:
            st.code(resp.text, language="text")
    except Exception:
        st.info("Metrics unavailable.")


# =====================================================================
# Tab 2: Trial Explorer (RAG Q&A)
# =====================================================================

with tab_explorer:
    st.header("Trial Explorer")
    st.write("RAG-powered clinical trial Q&A across all knowledge collections.")

    # Workflow selector
    workflow_options = [
        "auto", "protocol_design", "patient_matching", "site_selection",
        "eligibility_optimization", "adaptive_design", "safety_signal",
        "regulatory_docs", "competitive_intel", "diversity_assessment",
        "decentralized_planning", "general",
    ]
    selected_workflow = st.selectbox(
        "Workflow Focus",
        workflow_options,
        index=0,
        help="Select a workflow to guide the query, or leave as 'auto' for automatic routing.",
    )

    # Query input
    question = st.text_area(
        "Clinical Trial Question",
        placeholder="e.g., What are the latest adaptive design approaches for Phase II oncology trials?",
        height=100,
    )

    col_topk, col_guidelines = st.columns(2)
    with col_topk:
        top_k = st.slider("Evidence passages (top_k)", 1, 20, 5)
    with col_guidelines:
        include_guidelines = st.checkbox("Include guideline citations", value=True)

    if st.button("Search", key="explorer_search"):
        if question.strip():
            with st.spinner("Searching trial knowledge base..."):
                payload = {
                    "question": question.strip(),
                    "top_k": top_k,
                    "include_guidelines": include_guidelines,
                }
                if selected_workflow != "auto":
                    payload["workflow_type"] = selected_workflow

                result = api_post("/v1/trial/query", payload)

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
# Tab 3: Patient Matcher
# =====================================================================

with tab_matcher:
    st.header("Patient-Trial Matcher")
    st.write("Input a patient profile to find matching clinical trials.")

    col_left, col_right = st.columns(2)

    with col_left:
        st.subheader("Patient Profile")
        pt_age = st.number_input("Age", min_value=0, max_value=120, value=55)
        pt_sex = st.selectbox("Sex", ["male", "female"])
        pt_diagnosis = st.text_input(
            "Primary Diagnosis",
            placeholder="e.g., Non-small cell lung cancer, EGFR-mutated",
        )
        pt_biomarkers = st.text_area(
            "Biomarkers (one per line)",
            placeholder="HER2+\nEGFR T790M\nPD-L1 > 50%",
            height=80,
        )
        pt_medications = st.text_area(
            "Current Medications (one per line)",
            placeholder="Osimertinib 80mg daily\nPembrolizumab 200mg q3w",
            height=80,
        )

    with col_right:
        st.subheader("Additional Context")
        pt_variants = st.text_area(
            "Genomic Variants (one per line)",
            placeholder="EGFR L858R\nTP53 R248W",
            height=80,
        )
        pt_comorbidities = st.text_area(
            "Comorbidities (one per line)",
            placeholder="Type 2 diabetes\nHypertension",
            height=80,
        )
        pt_location = st.text_input(
            "Geographic Location",
            placeholder="e.g., Boston, MA, USA",
        )
        pt_therapeutic_area = st.selectbox(
            "Therapeutic Area",
            ["", "oncology", "cardiology", "neurology", "immunology",
             "rare_disease", "infectious_disease", "metabolic", "respiratory",
             "dermatology", "ophthalmology", "hematology", "gastroenterology", "other"],
        )
        pt_max_results = st.slider("Maximum Results", 1, 50, 10)

    if st.button("Find Matching Trials", key="match_search"):
        if pt_diagnosis.strip():
            with st.spinner("Matching patient to trials..."):
                payload = {
                    "age": pt_age,
                    "sex": pt_sex,
                    "diagnosis": pt_diagnosis.strip(),
                    "biomarkers": [b.strip() for b in pt_biomarkers.strip().split("\n") if b.strip()],
                    "medications": [m.strip() for m in pt_medications.strip().split("\n") if m.strip()],
                    "genomic_variants": [v.strip() for v in pt_variants.strip().split("\n") if v.strip()],
                    "comorbidities": [c.strip() for c in pt_comorbidities.strip().split("\n") if c.strip()],
                    "geographic_location": pt_location.strip() or None,
                    "therapeutic_area": pt_therapeutic_area or None,
                    "max_results": pt_max_results,
                }

                result = api_post("/v1/trial/match", payload)

            if result:
                matches = result.get("matches", [])
                st.subheader(f"Matching Trials ({len(matches)} found)")
                st.caption(f"Screened {result.get('total_screened', 0)} trials")

                for match in matches:
                    score = match.get("overall_score", 0)
                    score_color = NVIDIA_THEME["success"] if score >= 0.7 else (
                        NVIDIA_THEME["warning"] if score >= 0.4 else NVIDIA_THEME["danger"]
                    )
                    with st.expander(
                        f"{match.get('trial_id', 'N/A')} | "
                        f"{match.get('trial_title', 'Unknown')[:60]} | "
                        f"Score: {score:.2f}"
                    ):
                        mcol1, mcol2, mcol3 = st.columns(3)
                        with mcol1:
                            st.metric("Match Score", f"{score:.2f}")
                        with mcol2:
                            st.metric("Phase", match.get("phase", "N/A"))
                        with mcol3:
                            st.metric("Status", match.get("status", "N/A"))

                        inc_met = match.get("inclusion_met", 0)
                        inc_total = match.get("inclusion_total", 0)
                        exc_clear = match.get("exclusion_clear", 0)
                        exc_total = match.get("exclusion_total", 0)
                        if inc_total > 0 or exc_total > 0:
                            st.write(f"Inclusion: {inc_met}/{inc_total} met | Exclusion: {exc_clear}/{exc_total} clear")

                        st.write(f"Confidence: {match.get('confidence', 0):.2f}")
        else:
            st.warning("Please enter a primary diagnosis.")


# =====================================================================
# Tab 4: Protocol Analyzer
# =====================================================================

with tab_protocol:
    st.header("Protocol Analyzer")
    st.write("Optimize clinical trial protocols with complexity scoring and AI-driven recommendations.")

    protocol_summary = st.text_area(
        "Protocol Summary",
        placeholder="Paste protocol synopsis or summary here...",
        height=200,
    )

    pcol1, pcol2, pcol3 = st.columns(3)
    with pcol1:
        proto_therapeutic = st.selectbox(
            "Therapeutic Area",
            ["", "oncology", "cardiology", "neurology", "immunology",
             "rare_disease", "infectious_disease", "metabolic", "respiratory",
             "dermatology", "ophthalmology", "hematology", "gastroenterology"],
            key="proto_ta",
        )
        proto_phase = st.selectbox(
            "Phase",
            ["", "phase_i", "phase_i_ii", "phase_ii", "phase_ii_iii", "phase_iii", "phase_iv"],
            key="proto_phase",
        )
    with pcol2:
        proto_indication = st.text_input("Indication", placeholder="e.g., Advanced melanoma")
        proto_visits = st.number_input("Number of Visits", min_value=0, max_value=200, value=0)
    with pcol3:
        proto_procedures = st.number_input("Number of Procedures", min_value=0, max_value=200, value=0)

    proto_endpoints = st.text_area(
        "Endpoints (one per line)",
        placeholder="Overall survival\nProgression-free survival\nObjective response rate",
        height=80,
    )
    proto_criteria = st.text_area(
        "Key Eligibility Criteria (one per line)",
        placeholder="Age >= 18\nECOG 0-1\nAdequate organ function",
        height=80,
    )

    if st.button("Analyze Protocol", key="protocol_analyze"):
        if protocol_summary.strip():
            with st.spinner("Analyzing protocol complexity..."):
                payload = {
                    "protocol_summary": protocol_summary.strip(),
                    "therapeutic_area": proto_therapeutic or None,
                    "phase": proto_phase or None,
                    "indication": proto_indication.strip() or None,
                    "endpoints": [e.strip() for e in proto_endpoints.strip().split("\n") if e.strip()],
                    "eligibility_criteria": [c.strip() for c in proto_criteria.strip().split("\n") if c.strip()],
                    "visit_count": proto_visits if proto_visits > 0 else None,
                    "procedure_count": proto_procedures if proto_procedures > 0 else None,
                }

                result = api_post("/v1/trial/protocol/optimize", payload)

            if result:
                st.subheader("Complexity Assessment")

                rcol1, rcol2, rcol3 = st.columns(3)
                with rcol1:
                    complexity = result.get("complexity_score", 0)
                    st.metric("Complexity Score", f"{complexity:.3f}")
                with rcol2:
                    percentile = result.get("percentile_rank", 0)
                    st.metric("Percentile Rank", f"{percentile:.1f}%")
                with rcol3:
                    impact = result.get("estimated_enrollment_impact", 0)
                    st.metric("Enrollment Impact", f"{impact:.2f}")

                st.progress(complexity, text=f"Complexity: {complexity:.1%}")

                recs = result.get("optimization_recommendations", [])
                if recs:
                    st.subheader("Optimization Recommendations")
                    for rec in recs:
                        st.write(f"- {rec}")

                risks = result.get("risk_factors", [])
                if risks:
                    st.subheader("Risk Factors")
                    for risk in risks:
                        if isinstance(risk, dict):
                            st.write(f"- **{risk.get('factor', '')}**: {risk.get('description', '')}")
                        else:
                            st.write(f"- {risk}")

                # Generate report
                if st.button("Generate Report", key="proto_report"):
                    report_data = {
                        "report_type": "protocol_summary",
                        "format": "markdown",
                        "title": f"Protocol Analysis: {proto_indication or 'Clinical Trial'}",
                        "data": result,
                    }
                    report = api_post("/v1/reports/generate", report_data)
                    if report:
                        st.download_button(
                            "Download Report",
                            data=report.get("content", ""),
                            file_name=f"protocol_analysis_{report.get('report_id', 'report')}.md",
                            mime="text/markdown",
                        )
        else:
            st.warning("Please enter a protocol summary.")


# =====================================================================
# Tab 5: Competitive Dashboard
# =====================================================================

with tab_competitive:
    st.header("Competitive Landscape Dashboard")
    st.write("Analyze the competitive landscape by therapeutic area and indication.")

    ccol1, ccol2 = st.columns(2)
    with ccol1:
        comp_therapeutic = st.selectbox(
            "Therapeutic Area",
            ["oncology", "cardiology", "neurology", "immunology",
             "rare_disease", "infectious_disease", "metabolic", "respiratory",
             "dermatology", "ophthalmology", "hematology", "gastroenterology"],
            key="comp_ta",
        )
        comp_indication = st.text_input(
            "Indication (optional)",
            placeholder="e.g., Non-small cell lung cancer",
            key="comp_indication",
        )
    with ccol2:
        comp_mechanism = st.text_input(
            "Mechanism of Action (optional)",
            placeholder="e.g., PD-1 inhibitor",
            key="comp_mechanism",
        )
        comp_include_completed = st.checkbox("Include Completed Trials", value=False)
        comp_max = st.slider("Maximum Competitors", 5, 50, 20, key="comp_max")

    if st.button("Analyze Landscape", key="comp_analyze"):
        with st.spinner("Analyzing competitive landscape..."):
            payload = {
                "therapeutic_area": comp_therapeutic,
                "indication": comp_indication.strip() or None,
                "mechanism": comp_mechanism.strip() or None,
                "include_completed": comp_include_completed,
                "max_competitors": comp_max,
            }

            result = api_post("/v1/trial/competitive/landscape", payload)

        if result:
            st.subheader("Landscape Summary")
            st.markdown(result.get("landscape_summary", "No summary available."))

            competitors = result.get("competitors", [])
            if competitors:
                st.subheader(f"Competitors ({len(competitors)})")
                for comp in competitors:
                    if isinstance(comp, dict):
                        with st.expander(f"{comp.get('trial_id', 'N/A')} - {comp.get('sponsor', 'Unknown')}"):
                            st.write(f"**Phase:** {comp.get('phase', 'N/A')}")
                            st.write(f"**Indication:** {comp.get('indication', 'N/A')}")
                            st.write(f"**Mechanism:** {comp.get('mechanism', 'N/A')}")
                            progress = comp.get("enrollment_progress", 0)
                            st.progress(progress, text=f"Enrollment: {progress:.0%}")

            threat = result.get("threat_assessment", {})
            if threat:
                st.subheader("Threat Assessment")
                st.json(threat)

            enrollment_race = result.get("enrollment_race", [])
            if enrollment_race:
                st.subheader("Enrollment Race")
                for er in enrollment_race:
                    if isinstance(er, dict):
                        st.write(f"- {er.get('trial_id', 'N/A')}: {er.get('progress', 'N/A')}")

    # Reference data
    st.markdown("---")
    st.subheader("Reference Data")

    ref_col1, ref_col2, ref_col3 = st.columns(3)
    with ref_col1:
        if st.button("Therapeutic Areas", key="ref_ta"):
            data = api_get("/v1/trial/therapeutic-areas")
            if data:
                for ta in data.get("therapeutic_areas", []):
                    st.write(f"- {ta.get('name', ta.get('id', 'N/A'))}")

    with ref_col2:
        if st.button("Trial Phases", key="ref_phases"):
            data = api_get("/v1/trial/phases")
            if data:
                for phase in data.get("phases", []):
                    st.write(f"- **{phase.get('name', 'N/A')}:** {phase.get('description', '')}")

    with ref_col3:
        if st.button("Guidelines", key="ref_guidelines"):
            data = api_get("/v1/trial/guidelines")
            if data:
                for gl in data.get("guidelines", []):
                    st.write(f"- **{gl.get('name', 'N/A')}:** {gl.get('description', '')}")
