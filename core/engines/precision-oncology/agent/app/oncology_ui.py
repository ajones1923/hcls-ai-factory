"""
Oncology Intelligence Agent - MTB Workbench UI
=================================================
Streamlit application providing a Molecular Tumor Board workbench with
five tabs: Case Workbench, Evidence Explorer, Trial Finder, Therapy Ranker,
and Outcomes Dashboard.

Author: Adam Jones
Date: February 2026
"""

import os
from typing import Optional

import requests
import streamlit as st

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
API_BASE = os.environ.get("ONCO_API_BASE_URL", "http://localhost:8527")
PAGE_TITLE = "Oncology Intelligence MTB Workbench"
PAGE_ICON = "🧬"

CANCER_TYPES = [
    "B-Cell Acute Lymphoblastic Leukemia (B-ALL)",
    "T-Cell Acute Lymphoblastic Leukemia (T-ALL)",
    "Acute Myeloid Leukemia (AML)",
    "Chronic Myeloid Leukemia (CML)",
    "Non-Small Cell Lung Cancer (NSCLC)",
    "Small Cell Lung Cancer (SCLC)",
    "Breast Cancer",
    "Colorectal Cancer",
    "Pancreatic Cancer",
    "Melanoma",
    "Glioblastoma",
    "Medulloblastoma",
    "Neuroblastoma",
    "Rhabdomyosarcoma",
    "Wilms Tumor (Nephroblastoma)",
    "Ewing Sarcoma",
    "Retinoblastoma",
    "Ovarian Cancer",
    "Prostate Cancer",
    "Bladder Cancer",
    "Hepatocellular Carcinoma",
    "Renal Cell Carcinoma",
    "Head and Neck Squamous Cell Carcinoma",
    "Gastric Cancer",
    "Esophageal Cancer",
    "Cholangiocarcinoma",
    "Thyroid Cancer",
    "Other",
]

STAGES = ["I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IIIC", "IV", "IVA", "IVB"]

BIOMARKER_OPTIONS = [
    "EGFR+", "ALK+", "ROS1+", "BRAF V600E", "KRAS G12C",
    "MSI-H", "TMB-H", "PD-L1>=50%", "HER2+", "BRCA+",
    "NTRK fusion", "RET fusion", "MET amplification",
    "PIK3CA mutation", "FGFR alteration",
]

THERAPY_OPTIONS = [
    "Platinum-based chemotherapy",
    "Pembrolizumab",
    "Nivolumab",
    "Atezolizumab",
    "Osimertinib",
    "Erlotinib",
    "Gefitinib",
    "Crizotinib",
    "Alectinib",
    "Vemurafenib",
    "Dabrafenib + Trametinib",
    "Sotorasib",
    "Adagrasib",
    "Trastuzumab",
    "Bevacizumab",
    "Olaparib",
    "Larotrectinib",
    "Entrectinib",
    "Capmatinib",
    "Radiation therapy",
    "Surgery",
]


# ---------------------------------------------------------------------------
# API Client
# ---------------------------------------------------------------------------
@st.cache_resource
def get_api_session():
    """Create a reusable requests session for the oncology API."""
    session = requests.Session()
    session.headers.update({"Content-Type": "application/json"})
    return session


def api_get(endpoint: str, params: Optional[dict] = None) -> Optional[dict]:
    """GET request to the oncology API."""
    session = get_api_session()
    try:
        resp = session.get(f"{API_BASE}{endpoint}", params=params, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except requests.exceptions.ConnectionError:
        st.error(f"Cannot connect to API at {API_BASE}. Is the service running?")
        return None
    except Exception as exc:
        st.error(f"API error: {exc}")
        return None


def api_post(endpoint: str, payload: dict) -> Optional[dict]:
    """POST request to the oncology API."""
    session = get_api_session()
    try:
        resp = session.post(f"{API_BASE}{endpoint}", json=payload, timeout=120)
        resp.raise_for_status()
        return resp.json()
    except requests.exceptions.ConnectionError:
        st.error(f"Cannot connect to API at {API_BASE}. Is the service running?")
        return None
    except Exception as exc:
        st.error(f"API error: {exc}")
        return None


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
def render_sidebar():
    """Render sidebar with branding, service status, and navigation."""
    with st.sidebar:
        st.markdown(
            '<div style="background:#76b900;color:white;padding:12px 16px;'
            'border-radius:6px;text-align:center;font-size:18px;font-weight:bold;">'
            'HCLS AI Factory</div>',
            unsafe_allow_html=True,
        )
        st.markdown("## Oncology Intelligence Agent")
        st.markdown("*Molecular Tumor Board Workbench*")
        st.divider()

        # Service status
        st.markdown("### Service Status")
        health = api_get("/health")
        if health:
            status = health.get("status", "unknown")
            status_icon = "green" if status == "healthy" else "orange"
            st.markdown(
                f"**API:** :{status_icon}[{status}]"
            )
            services = health.get("services", {})
            for svc, ok in services.items():
                icon = "green" if ok else "red"
                label = svc.replace("_", " ").title()
                st.markdown(f"- :{icon}[{label}]")

            total_v = health.get("total_vectors", 0)
            st.metric("Total Vectors", f"{total_v:,}")
        else:
            st.markdown("**API:** :red[Offline]")

        st.divider()

        # Navigation links
        st.markdown("### Links")
        st.markdown("- [HCLS AI Factory](https://github.com/aj-geddes/hcls-ai-factory)")
        st.markdown("- [Milvus Attu](http://localhost:8000)")
        st.markdown("- [Grafana](http://localhost:3000)")
        st.markdown("- [Portal](http://localhost:8510)")

        st.divider()

        st.markdown("### \U0001f3af Demo Mode")
        if st.button("Load Demo Patient", key="load_demo"):
            import sys as _sys
            _lib_path = os.environ.get("HCLS_LIB_PATH", "/app/lib")
            _sys.path.insert(0, _lib_path)
            from hcls_common.demo_data import DEMO_ONCOLOGY
            st.session_state["demo_loaded"] = True
            st.session_state["demo_cancer_type"] = DEMO_ONCOLOGY["cancer_type"]
            st.session_state["demo_stage"] = DEMO_ONCOLOGY["stage"]
            st.session_state["demo_variants"] = DEMO_ONCOLOGY["variants"]
            st.toast("\u2705 Demo patient loaded! Use Case Workbench tab.", icon="\U0001f3af")
            st.rerun()

        st.divider()
        st.caption("Author: Adam Jones | February 2026")


# ---------------------------------------------------------------------------
# Tab 1: Case Workbench
# ---------------------------------------------------------------------------
def render_case_workbench():
    """Case creation and MTB packet generation."""
    st.header("Case Workbench")
    st.markdown("Create a patient case and generate a Molecular Tumor Board packet.")

    col1, col2 = st.columns(2)
    with col1:
        patient_id = st.text_input("Patient ID", value="", placeholder="PT-001")
        cancer_type = st.selectbox("Cancer Type", CANCER_TYPES)
        stage = st.selectbox("Stage", [""] + STAGES)
    with col2:
        st.markdown("**Biomarkers**")
        tmb = st.number_input("TMB (mut/Mb)", min_value=0.0, max_value=500.0, value=0.0, step=0.1)
        msi_status = st.selectbox("MSI Status", ["MSS", "MSI-L", "MSI-H"])
        pdl1_tps = st.number_input("PD-L1 TPS (%)", min_value=0, max_value=100, value=0)
        hrd_score = st.number_input("HRD Score", min_value=0.0, max_value=100.0, value=0.0, step=0.1)

    # Prior therapies
    prior_therapies = st.multiselect("Prior Therapies", THERAPY_OPTIONS)

    # Variant input
    st.subheader("Variants")
    input_mode = st.radio("Input Mode", ["Manual Entry", "Paste VCF"], horizontal=True)

    vcf_text = None
    variants_list = []

    if input_mode == "Paste VCF":
        vcf_text = st.text_area(
            "VCF Content",
            height=200,
            placeholder="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO ...",
        )
    else:
        st.markdown("Add variants below:")
        if "variant_rows" not in st.session_state:
            st.session_state.variant_rows = [{"gene": "", "variant": "", "type": "SNV"}]

        for i, vrow in enumerate(st.session_state.variant_rows):
            c1, c2, c3, c4 = st.columns([2, 3, 2, 1])
            with c1:
                vrow["gene"] = st.text_input(f"Gene##{i}", value=vrow.get("gene", ""), key=f"vg_{i}")
            with c2:
                vrow["variant"] = st.text_input(f"Variant##{i}", value=vrow.get("variant", ""), key=f"vv_{i}")
            with c3:
                vrow["type"] = st.selectbox(
                    f"Type##{i}",
                    ["SNV", "Insertion", "Deletion", "CNV", "Fusion", "Rearrangement"],
                    key=f"vt_{i}",
                )
            with c4:
                if st.button("X", key=f"vdel_{i}"):
                    st.session_state.variant_rows.pop(i)
                    st.rerun()

        if st.button("+ Add Variant"):
            st.session_state.variant_rows.append({"gene": "", "variant": "", "type": "SNV"})
            st.rerun()

        variants_list = [
            {"gene": v["gene"], "variant": v["variant"], "variant_type": v["type"]}
            for v in st.session_state.variant_rows
            if v.get("gene") and v.get("variant")
        ]

    # Create case
    st.divider()
    col_a, col_b = st.columns(2)

    with col_a:
        if st.button("Create Case", type="primary", use_container_width=True):
            if not patient_id:
                st.warning("Please enter a Patient ID.")
                return

            biomarkers = {
                "tmb": tmb,
                "msi_status": msi_status,
                "pdl1_tps": pdl1_tps,
                "hrd_score": hrd_score,
            }

            payload = {
                "patient_id": patient_id,
                "cancer_type": cancer_type,
                "stage": stage if stage else None,
                "biomarkers": biomarkers,
                "prior_therapies": prior_therapies,
            }
            if vcf_text and vcf_text.strip():
                payload["vcf_text"] = vcf_text
            elif variants_list:
                payload["variants"] = variants_list

            with st.spinner("Creating case ..."):
                result = api_post("/api/cases", payload)

            if result:
                st.success(f"Case created: **{result.get('case_id', 'N/A')}**")
                st.session_state["current_case_id"] = result.get("case_id")
                st.json(result)

                # Publish cross-agent event for case creation
                try:
                    import sys
                    _lib_path = os.environ.get("HCLS_LIB_PATH", "/app/lib")
                    sys.path.insert(0, _lib_path)
                    from hcls_common.event_bus import publish_event, EventType, PipelineStage
                    publish_event(
                        EventType.ONCOLOGY_CASE_CREATED,
                        source_stage=PipelineStage.ANNOTATION,
                        payload={
                            "case_id": result.get("case_id"),
                            "cancer_type": cancer_type,
                            "variant_count": len(variants_list),
                        },
                        patient_id=patient_id,
                    )
                except Exception:
                    pass

    with col_b:
        case_id = st.session_state.get("current_case_id", "")
        if st.button("Generate MTB Packet", use_container_width=True, disabled=not case_id):
            with st.spinner("Generating MTB packet ..."):
                mtb = api_post(f"/api/cases/{case_id}/mtb", {})

            if mtb:
                st.success("MTB Packet generated successfully.")

                # Publish cross-agent event for MTB generation
                try:
                    import sys
                    _lib_path = os.environ.get("HCLS_LIB_PATH", "/app/lib")
                    sys.path.insert(0, _lib_path)
                    from hcls_common.event_bus import publish_event, EventType, PipelineStage
                    publish_event(
                        EventType.THERAPY_RANKED,
                        source_stage=PipelineStage.ANNOTATION,
                        payload={
                            "case_id": case_id,
                            "therapy_count": len(mtb.get("therapy_ranking", [])),
                            "trial_matches": len(mtb.get("trial_matches", [])),
                            "variant_count": len(mtb.get("variants", [])),
                        },
                        patient_id=st.session_state.get("current_case_id", ""),
                    )
                except Exception:
                    pass

                _render_mtb_packet(mtb)


def _render_mtb_packet(mtb: dict):
    """Display an MTB packet in structured sections."""
    st.subheader("MTB Packet Summary")

    # Variant table
    variants = mtb.get("variants", [])
    if variants:
        st.markdown("**Actionable Variants**")
        st.dataframe(variants, use_container_width=True)

    # Evidence table
    evidence = mtb.get("evidence", [])
    if evidence:
        st.markdown("**Evidence Summary**")
        st.dataframe(evidence, use_container_width=True)

    # Therapy ranking
    therapies = mtb.get("therapy_ranking", [])
    if therapies:
        st.markdown("**Therapy Ranking**")
        for i, tx in enumerate(therapies, 1):
            name = tx.get("name", "Unknown")
            level = tx.get("evidence_level", "N/A")
            score = tx.get("score", 0)
            st.markdown(f"{i}. **{name}** - Evidence Level: `{level}` | Score: {score:.2f}")

    # Trial matches
    trial_matches = mtb.get("trial_matches", [])
    if trial_matches:
        st.markdown("**Matched Clinical Trials**")
        for t in trial_matches:
            nct = t.get("nct_id", "N/A")
            title = t.get("title", "Untitled")
            match_score = t.get("match_score", 0)
            st.markdown(f"- **{nct}**: {title} (match: {match_score:.2f})")

    # Open questions
    questions = mtb.get("open_questions", [])
    if questions:
        st.markdown("**Open Questions for Discussion**")
        for q in questions:
            st.markdown(f"- {q}")


# ---------------------------------------------------------------------------
# Tab 2: Evidence Explorer
# ---------------------------------------------------------------------------
def render_evidence_explorer():
    """RAG-powered evidence Q&A interface."""
    st.header("Evidence Explorer")
    st.markdown("Ask clinical questions powered by oncology knowledge collections.")

    # Filters in sidebar area
    with st.expander("Filters", expanded=False):
        filter_cancer = st.selectbox("Cancer Type Filter", [""] + CANCER_TYPES, key="ev_cancer")
        filter_gene = st.text_input("Gene Filter", key="ev_gene", placeholder="e.g. EGFR, BRAF")

    question = st.text_input(
        "Ask a question",
        placeholder="What are the approved therapies for EGFR-mutant NSCLC?",
        key="ev_question",
    )

    if st.button("Ask", type="primary", key="ev_ask"):
        if not question:
            st.warning("Please enter a question.")
            return

        payload = {
            "question": question,
            "cancer_type": filter_cancer if filter_cancer else None,
            "gene": filter_gene if filter_gene else None,
            "top_k": 10,
            "include_follow_ups": True,
        }

        with st.spinner("Searching evidence ..."):
            result = api_post("/api/ask", payload)

        if result:
            st.session_state["ev_result"] = result

    # Display result
    result = st.session_state.get("ev_result")
    if result:
        st.markdown("### Answer")
        st.markdown(result.get("answer", "No answer."))

        confidence = result.get("confidence", 0)
        st.progress(min(confidence, 1.0), text=f"Confidence: {confidence:.1%}")

        # Sources
        sources = result.get("sources", [])
        if sources:
            st.markdown(f"### Evidence Sources ({len(sources)})")
            for src in sources:
                coll = src.get("collection", "unknown")
                score = src.get("score", 0)
                text = src.get("text", "")
                badge_color = {
                    "onco_variants": "blue",
                    "onco_therapies": "green",
                    "onco_resistance": "red",
                    "onco_pathways": "violet",
                    "onco_biomarkers": "orange",
                    "onco_trials": "gray",
                }.get(coll, "gray")

                st.markdown(f":{badge_color}[{coll}] Score: {score:.3f}")
                st.markdown(f"> {text[:300]}{'...' if len(text) > 300 else ''}")
                st.divider()

        timing = result.get("processing_time_ms", 0)
        st.caption(f"Processing time: {timing:.0f} ms")

        # Follow-up questions
        follow_ups = result.get("follow_up_questions", [])
        if follow_ups:
            st.markdown("### Follow-up Questions")
            for fq in follow_ups:
                if st.button(fq, key=f"fq_{hash(fq)}"):
                    st.session_state["ev_question"] = fq
                    st.rerun()


# ---------------------------------------------------------------------------
# Tab 3: Trial Finder
# ---------------------------------------------------------------------------
def render_trial_finder():
    """Clinical trial matching interface."""
    st.header("Trial Finder")
    st.markdown("Find matching clinical trials based on cancer type and biomarker profile.")

    col1, col2 = st.columns(2)
    with col1:
        cancer_type = st.selectbox("Cancer Type", CANCER_TYPES, key="tf_cancer")
        stage = st.selectbox("Stage", [""] + STAGES, key="tf_stage")
        age = st.number_input("Patient Age", min_value=0, max_value=120, value=60, key="tf_age")
    with col2:
        st.markdown("**Biomarkers**")
        selected_biomarkers = []
        for bm in BIOMARKER_OPTIONS:
            if st.checkbox(bm, key=f"tf_bm_{bm}"):
                selected_biomarkers.append(bm)

    biomarker_dict = {bm: True for bm in selected_biomarkers}

    if st.button("Find Trials", type="primary", key="tf_search"):
        payload = {
            "cancer_type": cancer_type,
            "biomarkers": biomarker_dict,
            "stage": stage if stage else None,
            "age": age,
            "top_k": 15,
        }

        with st.spinner("Matching trials ..."):
            result = api_post("/api/trials/match", payload)

        if result:
            st.session_state["tf_result"] = result

    # Display results
    result = st.session_state.get("tf_result")
    if result:
        matches = result.get("matches", [])
        st.markdown(f"### Matched Trials ({len(matches)})")

        if not matches:
            st.info("No matching trials found for the given criteria.")
            return

        for i, trial in enumerate(matches, 1):
            nct_id = trial.get("nct_id", "N/A")
            title = trial.get("title", "Untitled")
            match_score = trial.get("match_score", 0)
            phase = trial.get("phase", "N/A")
            status = trial.get("status", "Unknown")
            explanation = trial.get("explanation", "")
            eligibility = trial.get("eligibility_status", "Unknown")

            elig_icon = "green" if eligibility == "Eligible" else (
                "orange" if eligibility == "Potentially Eligible" else "red"
            )

            with st.container(border=True):
                st.markdown(f"**{i}. {nct_id}** | Phase: {phase} | Status: {status}")
                st.markdown(f"**{title}**")
                st.markdown(f"Match Score: **{match_score:.2f}** | Eligibility: :{elig_icon}[{eligibility}]")
                if explanation:
                    st.markdown(f"*{explanation}*")

        timing = result.get("processing_time_ms", 0)
        st.caption(f"Processing time: {timing:.0f} ms")


# ---------------------------------------------------------------------------
# Tab 4: Therapy Ranker
# ---------------------------------------------------------------------------
def render_therapy_ranker():
    """Therapy ranking interface based on molecular profile."""
    st.header("Therapy Ranker")
    st.markdown("Rank therapy options based on molecular variants, biomarkers, and treatment history.")

    cancer_type = st.selectbox("Cancer Type", CANCER_TYPES, key="tr_cancer")

    # Variant entry
    st.subheader("Variants")
    if "tr_variant_rows" not in st.session_state:
        st.session_state.tr_variant_rows = [{"gene": "", "variant": "", "type": "SNV"}]

    for i, vrow in enumerate(st.session_state.tr_variant_rows):
        c1, c2, c3, c4 = st.columns([2, 3, 2, 1])
        with c1:
            vrow["gene"] = st.text_input(f"Gene##{i}", value=vrow.get("gene", ""), key=f"trg_{i}")
        with c2:
            vrow["variant"] = st.text_input(f"Variant##{i}", value=vrow.get("variant", ""), key=f"trv_{i}")
        with c3:
            vrow["type"] = st.selectbox(
                f"Type##{i}",
                ["SNV", "Insertion", "Deletion", "CNV", "Fusion", "Rearrangement"],
                key=f"trt_{i}",
            )
        with c4:
            if st.button("X", key=f"trdel_{i}"):
                st.session_state.tr_variant_rows.pop(i)
                st.rerun()

    if st.button("+ Add Variant", key="tr_add_var"):
        st.session_state.tr_variant_rows.append({"gene": "", "variant": "", "type": "SNV"})
        st.rerun()

    # Biomarkers
    col1, col2 = st.columns(2)
    with col1:
        tr_tmb = st.number_input("TMB (mut/Mb)", min_value=0.0, value=0.0, step=0.1, key="tr_tmb")
        tr_msi = st.selectbox("MSI Status", ["MSS", "MSI-L", "MSI-H"], key="tr_msi")
    with col2:
        tr_pdl1 = st.number_input("PD-L1 TPS (%)", min_value=0, max_value=100, value=0, key="tr_pdl1")
        tr_hrd = st.number_input("HRD Score", min_value=0.0, value=0.0, step=0.1, key="tr_hrd")

    # Prior therapies
    tr_prior = st.multiselect("Prior Therapies", THERAPY_OPTIONS, key="tr_prior")

    if st.button("Rank Therapies", type="primary", key="tr_rank"):
        variants = [
            {"gene": v["gene"], "variant": v["variant"], "variant_type": v["type"]}
            for v in st.session_state.tr_variant_rows
            if v.get("gene") and v.get("variant")
        ]

        biomarkers = {
            "tmb": tr_tmb,
            "msi_status": tr_msi,
            "pdl1_tps": tr_pdl1,
            "hrd_score": tr_hrd,
        }

        payload = {
            "cancer_type": cancer_type,
            "variants": variants,
            "biomarkers": biomarkers,
            "prior_therapies": tr_prior,
        }

        with st.spinner("Ranking therapies ..."):
            result = api_post("/api/therapies/rank", payload)

        if result:
            st.session_state["tr_result"] = result

    # Display results
    result = st.session_state.get("tr_result")
    if result:
        therapies = result.get("therapies", [])
        st.markdown(f"### Ranked Therapies ({len(therapies)})")

        if not therapies:
            st.info("No therapy recommendations available for the given profile.")
            return

        for i, tx in enumerate(therapies, 1):
            name = tx.get("name", "Unknown")
            score = tx.get("score", 0)
            evidence_level = tx.get("evidence_level", "N/A")
            mechanism = tx.get("mechanism", "")
            guideline_ref = tx.get("guideline_reference", "")
            resistance_warnings = tx.get("resistance_warnings", [])

            # Evidence level badge color
            ev_color = {
                "1": "green", "1A": "green", "1B": "green",
                "2": "blue", "2A": "blue", "2B": "blue",
                "3": "orange", "3A": "orange", "3B": "orange",
                "4": "red", "R1": "red", "R2": "red",
            }.get(str(evidence_level).upper(), "gray")

            with st.container(border=True):
                st.markdown(
                    f"**{i}. {name}** | :{ev_color}[Level {evidence_level}] | "
                    f"Score: {score:.2f}"
                )
                if mechanism:
                    st.markdown(f"Mechanism: {mechanism}")
                if guideline_ref:
                    st.caption(f"Guideline: {guideline_ref}")
                if resistance_warnings:
                    for rw in resistance_warnings:
                        st.warning(f"Resistance: {rw}", icon="!")

        timing = result.get("processing_time_ms", 0)
        st.caption(f"Processing time: {timing:.0f} ms")


# ---------------------------------------------------------------------------
# Tab 5: Outcomes Dashboard
# ---------------------------------------------------------------------------
def render_outcomes_dashboard():
    """Knowledge stats and event log overview."""
    st.header("Outcomes Dashboard")

    # Knowledge stats
    st.subheader("Knowledge Base Statistics")
    stats = api_get("/knowledge/stats")
    if stats:
        col1, col2, col3, col4, col5 = st.columns(5)
        col1.metric("Targets", f"{stats.get('target_count', 0):,}")
        col2.metric("Therapies", f"{stats.get('therapy_count', 0):,}")
        col3.metric("Resistance", f"{stats.get('resistance_count', 0):,}")
        col4.metric("Pathways", f"{stats.get('pathway_count', 0):,}")
        col5.metric("Biomarkers", f"{stats.get('biomarker_count', 0):,}")

        # Collection sizes chart
        st.subheader("Collection Sizes")
        coll_counts = stats.get("collection_counts", {})
        if coll_counts:
            import pandas as pd

            df = pd.DataFrame(
                [{"Collection": k, "Vectors": v} for k, v in coll_counts.items()]
            )
            df = df.sort_values("Vectors", ascending=True)
            st.bar_chart(df.set_index("Collection"), horizontal=True)

    # Recent events
    st.subheader("Recent Events")
    events = api_get("/api/events", params={"limit": 20})
    if events and events.get("events"):
        for evt in events["events"]:
            ts = evt.get("timestamp", "")
            etype = evt.get("event_type", "unknown")
            details = evt.get("details", {})
            user = evt.get("user", "system")


            st.markdown(
                f"**{etype}** | {ts[:19]} | User: {user}"
            )
            if details:
                st.json(details)
            st.divider()
    else:
        st.info("No events recorded yet.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    st.set_page_config(
        page_title=PAGE_TITLE,
        page_icon=PAGE_ICON,
        layout="wide",
        initial_sidebar_state="expanded",
    )

    render_sidebar()

    st.title(PAGE_TITLE)
    st.markdown(
        "Clinical decision support for molecular tumor boards, powered by "
        "RAG-driven oncology intelligence."
    )

    st.warning(
        "**Clinical Decision Support Tool** — This system provides evidence-based guidance "
        "for research and clinical decision support only. All recommendations must be verified "
        "by a qualified healthcare professional. Not FDA-cleared. Not a substitute for professional "
        "clinical judgment."
    )

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Case Workbench",
        "Evidence Explorer",
        "Trial Finder",
        "Therapy Ranker",
        "Outcomes Dashboard",
    ])

    with tab1:
        render_case_workbench()

    with tab2:
        render_evidence_explorer()

    with tab3:
        render_trial_finder()

    with tab4:
        render_therapy_ranker()

    with tab5:
        render_outcomes_dashboard()


if __name__ == "__main__":
    main()
