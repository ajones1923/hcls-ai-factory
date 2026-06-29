"""Streamlit UI for Imaging Intelligence Agent.

Medical imaging RAG chat interface with NIM integration,
comparative analysis, workflow demos, and report export.
Runs on port 8525.

Author: Adam Jones
Date: February 2026
"""

import json
import time
from datetime import datetime
from typing import Dict, List, Optional

import streamlit as st
from loguru import logger

# =====================================================================
# Page Configuration
# =====================================================================

st.set_page_config(
    page_title="Imaging Intelligence Agent",
    page_icon="🩻",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "About": (
            "Imaging Intelligence Agent -- HCLS AI Factory\n\n"
            "Multi-collection RAG engine for medical imaging intelligence.\n"
            "Apache 2.0 Licensed."
        ),
    },
)


# =====================================================================
# Engine Initialization (cached)
# =====================================================================


@st.cache_resource
def init_engine():
    """Initialize all agent components once and cache across sessions."""
    from config.settings import settings
    from sentence_transformers import SentenceTransformer
    from src.collections import ImagingCollectionManager
    from src.nim.service_manager import NIMServiceManager
    from src.rag_engine import ImagingRAGEngine

    # 1. Collection manager
    manager = ImagingCollectionManager(
        host=settings.MILVUS_HOST,
        port=settings.MILVUS_PORT,
    )
    try:
        manager.connect()
        manager.ensure_collections()
        logger.info("Milvus connected and collections ensured")
    except Exception as e:
        logger.warning(f"Milvus connection deferred: {e}")

    # 2. Embedding model
    embedder = SentenceTransformer(settings.EMBEDDING_MODEL)
    logger.info(f"Loaded embedding model: {settings.EMBEDDING_MODEL}")

    # 3. NIM service manager (includes LLM with Anthropic fallback)
    nim_manager = NIMServiceManager(settings)

    # 4. RAG engine
    engine = ImagingRAGEngine(
        collection_manager=manager,
        embedder=embedder,
        llm_client=nim_manager.llm,
        nim_service_manager=nim_manager,
    )

    return {
        "manager": manager,
        "embedder": embedder,
        "nim_manager": nim_manager,
        "engine": engine,
        "settings": settings,
    }


def safe_init():
    """Attempt initialization with error handling for display."""
    try:
        return init_engine()
    except Exception as e:
        st.error(f"Initialization error: {e}")
        logger.error(f"Engine init failed: {e}")
        return None


# =====================================================================
# Session State Defaults
# =====================================================================


def init_session_state():
    """Initialize Streamlit session state variables."""
    defaults = {
        "conversation_history": [],
        "last_evidence": None,
        "active_collections": [
            "imaging_literature", "imaging_trials", "imaging_findings",
            "imaging_protocols", "imaging_devices", "imaging_anatomy",
            "imaging_benchmarks", "imaging_guidelines",
            "imaging_report_templates", "imaging_datasets",
        ],
        "modality_filter": "All",
        "body_region_filter": "All",
        "year_range": (2015, 2026),
        "search_top_k": 5,
        "nim_mode": "auto",
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


init_session_state()


# =====================================================================
# Helpers
# =====================================================================

MODALITIES = ["All", "ct", "mri", "xray", "cxr", "ultrasound", "pet", "pet_ct", "mammography"]
BODY_REGIONS = [
    "All", "head", "neck", "chest", "abdomen", "pelvis",
    "spine", "extremity", "brain", "cardiac", "breast",
]
COLLECTION_LABELS = {
    "imaging_literature": "Literature",
    "imaging_trials": "Clinical Trials",
    "imaging_findings": "Findings",
    "imaging_protocols": "Protocols",
    "imaging_devices": "AI Devices",
    "imaging_anatomy": "Anatomy",
    "imaging_benchmarks": "Benchmarks",
    "imaging_guidelines": "Guidelines",
    "imaging_report_templates": "Report Templates",
    "imaging_datasets": "Datasets",
}

NIM_STATUS_COLORS = {
    "available": "🟢",
    "mock": "🟡",
    "unavailable": "🔴",
}

WORKFLOW_OPTIONS = {
    "ct_head_hemorrhage": "CT Head -- Hemorrhage Detection",
    "ct_chest_lung_nodule": "CT Chest -- Lung Nodule Analysis",
    "cxr_rapid_findings": "CXR -- Rapid Findings Triage",
    "mri_brain_ms_lesion": "MRI Brain -- MS Lesion Quantification",
}


def format_score_color(score: float, settings) -> str:
    """Return a colored indicator based on citation relevance score."""
    if score >= settings.CITATION_HIGH_THRESHOLD:
        return "🟢"
    elif score >= settings.CITATION_MEDIUM_THRESHOLD:
        return "🟡"
    return "⚪"


def build_conversation_context(history: list, max_turns: int = 3) -> str:
    """Build conversation context string from recent history."""
    if not history:
        return ""
    recent = history[-max_turns * 2:]
    parts = []
    for msg in recent:
        role = msg.get("role", "user")
        content = msg.get("content", "")
        parts.append(f"{role.upper()}: {content[:500]}")
    return "\n".join(parts)


# =====================================================================
# Sidebar
# =====================================================================


def render_sidebar(ctx: Optional[Dict]):
    """Render the sidebar with NIM status, collection stats, and filters."""
    with st.sidebar:
        st.markdown("## Imaging Intelligence Agent")
        st.caption("Multi-collection RAG for medical imaging")

        if ctx is None:
            st.warning("Engine not initialized")
            return

        ctx["settings"]
        nim_manager = ctx["nim_manager"]
        manager = ctx["manager"]

        # -- NIM Status Panel --
        st.markdown("### NIM Services")
        try:
            nim_status = nim_manager.check_all_services()
        except Exception:
            nim_status = {
                "vista3d": "unavailable",
                "maisi": "unavailable",
                "vila_m3": "unavailable",
                "llm": "unavailable",
            }

        nim_labels = {"vista3d": "VISTA-3D", "maisi": "MAISI", "vila_m3": "VILA-M3", "llm": "Llama-3 / Claude"}
        cols = st.columns(2)
        for i, (key, label) in enumerate(nim_labels.items()):
            status = nim_status.get(key, "unavailable")
            icon = NIM_STATUS_COLORS.get(status, "🔴")
            cols[i % 2].markdown(f"{icon} **{label}**")

        st.divider()

        # -- Collection Stats --
        st.markdown("### Collection Stats")
        try:
            stats = manager.get_collection_stats()
        except Exception:
            stats = {}

        if stats:
            for coll, label in COLLECTION_LABELS.items():
                count = stats.get(coll, 0)
                st.metric(label=label, value=f"{count:,}")
        else:
            st.info("Milvus not connected -- stats unavailable")

        st.divider()

        # -- Filters --
        st.markdown("### Filters")
        st.session_state.modality_filter = st.selectbox(
            "Modality", MODALITIES, index=0,
        )
        st.session_state.body_region_filter = st.selectbox(
            "Body Region", BODY_REGIONS, index=0,
        )
        st.session_state.year_range = st.slider(
            "Year Range", 2000, 2026, st.session_state.year_range,
        )

        st.divider()

        # -- Collection Toggles --
        st.markdown("### Collections to Search")
        active = []
        for coll, label in COLLECTION_LABELS.items():
            checked = st.checkbox(label, value=(coll in st.session_state.active_collections), key=f"coll_{coll}")
            if checked:
                active.append(coll)
        st.session_state.active_collections = active


# =====================================================================
# Tab: Ask (Chat Interface)
# =====================================================================


def render_ask_tab(ctx: Optional[Dict]):
    """Render the chat-based RAG question-answer interface."""
    if ctx is None:
        st.warning("Engine not initialized. Check the sidebar for connection status.")
        return

    engine = ctx["engine"]
    settings = ctx["settings"]

    # Display conversation history
    for msg in st.session_state.conversation_history:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])

    # Chat input
    user_question = st.chat_input("Ask about medical imaging AI...")
    if not user_question:
        return

    # Display user message
    st.session_state.conversation_history.append({"role": "user", "content": user_question})
    with st.chat_message("user"):
        st.markdown(user_question)

    # Build query kwargs
    kwargs = {
        "top_k_per_collection": st.session_state.search_top_k,
        "collections_filter": st.session_state.active_collections or None,
    }
    if st.session_state.modality_filter != "All":
        kwargs["modality_filter"] = st.session_state.modality_filter
    if st.session_state.year_range[0] > 2000:
        kwargs["year_min"] = st.session_state.year_range[0]
    if st.session_state.year_range[1] < 2026:
        kwargs["year_max"] = st.session_state.year_range[1]

    conversation_ctx = build_conversation_context(st.session_state.conversation_history)

    # Retrieve evidence first for the sidebar
    with st.spinner("Searching collections..."):
        evidence = engine.retrieve(user_question, **kwargs)
        st.session_state.last_evidence = evidence

    # Stream the LLM response
    with st.chat_message("assistant"):
        message_placeholder = st.empty()
        full_response = ""
        try:
            for chunk in engine.query_stream(user_question, conversation_context=conversation_ctx, **kwargs):
                full_response += chunk
                message_placeholder.markdown(full_response + "▌")
            message_placeholder.markdown(full_response)
        except Exception as e:
            full_response = f"Error generating response: {e}"
            message_placeholder.error(full_response)

    st.session_state.conversation_history.append({"role": "assistant", "content": full_response})

    # Evidence sidebar
    if evidence and evidence.hits:
        with st.expander(f"Evidence ({evidence.hit_count} hits from {evidence.total_collections_searched} collections, {evidence.search_time_ms:.0f}ms)", expanded=False):
            grouped = evidence.hits_by_collection()
            for coll_name, hits in grouped.items():
                label = COLLECTION_LABELS.get(coll_name, coll_name)
                st.markdown(f"**{label}** ({len(hits)} hits)")
                for hit in hits:
                    color = format_score_color(hit.score, settings)
                    st.markdown(f"  {color} `{hit.id}` (score: {hit.score:.3f})")
                    st.caption(hit.text[:300])
                st.markdown("---")


# =====================================================================
# Tab: Comparative Analysis
# =====================================================================


def render_comparative_tab(ctx: Optional[Dict]):
    """Render the two-entity comparative analysis interface."""
    if ctx is None:
        st.warning("Engine not initialized.")
        return

    engine = ctx["engine"]
    settings = ctx["settings"]

    st.markdown("### Comparative Analysis")
    st.caption("Compare two imaging modalities, techniques, devices, or pathologies side by side.")

    col1, col2 = st.columns(2)
    with col1:
        entity_a = st.text_input("Entity A", placeholder="e.g., CT head", key="comp_a")
    with col2:
        entity_b = st.text_input("Entity B", placeholder="e.g., MRI brain", key="comp_b")

    comparison_question = st.text_input(
        "Comparison question (optional)",
        placeholder="e.g., Which is better for detecting acute hemorrhage?",
        key="comp_question",
    )

    if st.button("Compare", type="primary", key="comp_run"):
        if not entity_a or not entity_b:
            st.warning("Please enter both entities.")
            return

        query = comparison_question or f"{entity_a} vs {entity_b}"
        with st.spinner("Running comparative retrieval..."):
            try:
                result = engine.retrieve_comparative(query)
            except Exception as e:
                st.error(f"Comparison failed: {e}")
                return

        col_a, col_b = st.columns(2)
        with col_a:
            st.markdown(f"#### {result.entity_a}")
            st.metric("Evidence hits", result.evidence_a.hit_count)
            for hit in result.evidence_a.hits[:10]:
                color = format_score_color(hit.score, settings)
                st.markdown(f"{color} **{COLLECTION_LABELS.get(hit.collection, hit.collection)}** ({hit.score:.3f})")
                st.caption(hit.text[:250])

        with col_b:
            st.markdown(f"#### {result.entity_b}")
            st.metric("Evidence hits", result.evidence_b.hit_count)
            for hit in result.evidence_b.hits[:10]:
                color = format_score_color(hit.score, settings)
                st.markdown(f"{color} **{COLLECTION_LABELS.get(hit.collection, hit.collection)}** ({hit.score:.3f})")
                st.caption(hit.text[:250])

        if result.comparison_context:
            with st.expander("Domain Knowledge Context", expanded=True):
                st.markdown(result.comparison_context)

        # Synthesize comparison via LLM
        with st.spinner("Synthesizing comparison..."):
            try:
                answer = engine.llm_client.generate([
                    {"role": "system", "content": "You are a medical imaging comparison expert. Provide a structured comparison based on the evidence."},
                    {"role": "user", "content": f"Compare {result.entity_a} vs {result.entity_b}. Question: {query}\n\nContext:\n{result.comparison_context}"},
                ])
                st.markdown("### Synthesized Comparison")
                st.markdown(answer)
            except Exception as e:
                st.warning(f"LLM synthesis unavailable: {e}")


# =====================================================================
# Tab: Workflow Demo
# =====================================================================


def render_workflow_tab(ctx: Optional[Dict]):
    """Render the imaging workflow demo interface."""
    if ctx is None:
        st.warning("Engine not initialized.")
        return

    from src.workflows import WORKFLOW_REGISTRY

    st.markdown("### Imaging Workflow Demo")
    st.caption("Run reference imaging analysis workflows in mock mode to demonstrate clinical pipeline capabilities.")

    selected_workflow = st.selectbox(
        "Select workflow",
        options=list(WORKFLOW_OPTIONS.keys()),
        format_func=lambda x: WORKFLOW_OPTIONS[x],
    )

    # Show workflow info
    workflow_class = WORKFLOW_REGISTRY.get(selected_workflow)
    if workflow_class:
        info_instance = workflow_class(mock_mode=True)
        info = info_instance.get_workflow_info()
        col1, col2, col3 = st.columns(3)
        col1.metric("Modality", info.get("modality", "N/A").upper())
        col2.metric("Body Region", info.get("body_region", "N/A").title())
        col3.metric("Target Latency", f"{info.get('target_latency_sec', 0):.0f}s")

        st.markdown(f"**Models used:** {', '.join(info.get('models_used', []))}")

    if st.button("Run Demo", type="primary", key="workflow_run"):
        with st.spinner(f"Running {WORKFLOW_OPTIONS.get(selected_workflow, selected_workflow)}..."):
            try:
                workflow = workflow_class(mock_mode=True)
                result = workflow.run()
            except Exception as e:
                st.error(f"Workflow failed: {e}")
                return

        # Status header
        status_color = "🟢" if result.status.value == "completed" else "🔴"
        st.markdown(f"### {status_color} Workflow Result: {result.workflow_name}")
        st.caption(f"Completed in {result.inference_time_ms:.1f}ms {'(mock)' if result.is_mock else '(live)'}")

        # Severity badge
        severity_colors = {
            "critical": "🔴", "urgent": "🟠", "significant": "🟡",
            "routine": "🟢", "normal": "⚪",
        }
        sev = result.severity.value
        st.markdown(f"**Severity:** {severity_colors.get(sev, '⚪')} {sev.upper()}")
        st.markdown(f"**Classification:** `{result.classification}`")

        # Findings
        st.markdown("#### Findings")
        for i, finding in enumerate(result.findings):
            sev_icon = severity_colors.get(finding.get("severity", "routine"), "⚪")
            st.markdown(f"{sev_icon} **Finding {i + 1}:** {finding.get('description', 'N/A')}")
            if finding.get("recommendation"):
                st.info(f"Recommendation: {finding['recommendation']}")

        # Measurements
        if result.measurements:
            st.markdown("#### Measurements")
            meas_cols = st.columns(min(len(result.measurements), 4))
            for i, (key, val) in enumerate(result.measurements.items()):
                meas_cols[i % len(meas_cols)].metric(
                    label=key.replace("_", " ").title(),
                    value=f"{val:.2f}" if isinstance(val, float) else str(val),
                )

        # Raw JSON
        with st.expander("Raw Result (JSON)"):
            st.json(result.model_dump())


# =====================================================================
# Tab: Reports
# =====================================================================


def render_reports_tab(ctx: Optional[Dict]):
    """Render the report export interface."""
    st.markdown("### Report Export")
    st.caption("Generate downloadable reports from your conversation and evidence.")

    history = st.session_state.conversation_history
    evidence = st.session_state.last_evidence

    if not history:
        st.info("No conversation history to export. Ask a question in the Ask tab first.")
        return

    col1, col2, col3 = st.columns(3)

    # -- Markdown Export --
    with col1:
        if st.button("Export Markdown", key="export_md"):
            md = _generate_markdown_report(history, evidence)
            st.download_button(
                label="Download .md",
                data=md,
                file_name=f"imaging_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md",
                mime="text/markdown",
                key="dl_md",
            )

    # -- JSON Export --
    with col2:
        if st.button("Export JSON", key="export_json"):
            report_data = {
                "generated_at": datetime.now().isoformat(),
                "conversation": history,
                "evidence_count": evidence.hit_count if evidence else 0,
                "collections_searched": evidence.total_collections_searched if evidence else 0,
            }
            json_str = json.dumps(report_data, indent=2, default=str)
            st.download_button(
                label="Download .json",
                data=json_str,
                file_name=f"imaging_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                mime="application/json",
                key="dl_json",
            )

    # -- PDF Export --
    with col3:
        if st.button("Export PDF", key="export_pdf"):
            try:
                pdf_bytes = _generate_pdf_report(history, evidence)
                st.download_button(
                    label="Download .pdf",
                    data=pdf_bytes,
                    file_name=f"imaging_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
                    mime="application/pdf",
                    key="dl_pdf",
                )
            except ImportError:
                st.warning("PDF export requires the `reportlab` package. Install with: pip install reportlab")
            except Exception as e:
                st.error(f"PDF generation failed: {e}")


def _generate_markdown_report(history: list, evidence) -> str:
    """Generate a Markdown report from conversation and evidence."""
    lines = [
        "# Imaging Intelligence Agent Report",
        f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "---",
        "",
        "## Conversation",
        "",
    ]
    for msg in history:
        role = msg.get("role", "user").title()
        content = msg.get("content", "")
        lines.append(f"### {role}")
        lines.append(content)
        lines.append("")

    if evidence and evidence.hits:
        lines.append("---")
        lines.append("")
        lines.append("## Evidence")
        lines.append(f"Total hits: {evidence.hit_count} across {evidence.total_collections_searched} collections")
        lines.append("")
        for hit in evidence.hits[:20]:
            lines.append(f"- **[{hit.collection}]** `{hit.id}` (score: {hit.score:.3f}): {hit.text[:200]}")
        lines.append("")

    lines.append("---")
    lines.append("*Generated by HCLS AI Factory -- Imaging Intelligence Agent*")
    return "\n".join(lines)


def _generate_pdf_report(history: list, evidence) -> bytes:
    """Generate a PDF report using reportlab."""
    from io import BytesIO

    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer

    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    styles = getSampleStyleSheet()
    story = []

    title_style = ParagraphStyle("Title", parent=styles["Title"], fontSize=18)
    story.append(Paragraph("Imaging Intelligence Agent Report", title_style))
    story.append(Spacer(1, 12))
    story.append(Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", styles["Normal"]))
    story.append(Spacer(1, 24))

    heading_style = ParagraphStyle("Heading", parent=styles["Heading2"], fontSize=14)
    for msg in history:
        role = msg.get("role", "user").title()
        content = msg.get("content", "")[:2000]
        story.append(Paragraph(role, heading_style))
        story.append(Paragraph(content.replace("\n", "<br/>"), styles["Normal"]))
        story.append(Spacer(1, 12))

    if evidence and evidence.hits:
        story.append(Paragraph("Evidence", heading_style))
        for hit in evidence.hits[:15]:
            text = f"[{hit.collection}] {hit.id} (score: {hit.score:.3f}): {hit.text[:200]}"
            story.append(Paragraph(text, styles["Normal"]))
            story.append(Spacer(1, 4))

    doc.build(story)
    return buffer.getvalue()


# =====================================================================
# Tab: Settings
# =====================================================================


def render_settings_tab(ctx: Optional[Dict]):
    """Render the settings / configuration interface."""
    st.markdown("### Settings")

    if ctx is None:
        st.warning("Engine not initialized.")
        return

    settings = ctx["settings"]

    # Search configuration
    st.markdown("#### Search Configuration")
    st.session_state.search_top_k = st.slider(
        "Results per collection (top_k)",
        min_value=1, max_value=20, value=st.session_state.search_top_k,
    )

    st.markdown("#### NIM Mode")
    st.session_state.nim_mode = st.radio(
        "NIM service mode",
        options=["auto", "local", "mock"],
        index=0,
        help=(
            "auto: Use local NIM if available, fallback to mock. "
            "local: Require live NIM endpoints. "
            "mock: Always use synthetic responses."
        ),
    )

    st.markdown("#### Citation Thresholds")
    col1, col2 = st.columns(2)
    with col1:
        st.number_input(
            "High relevance threshold",
            min_value=0.0, max_value=1.0,
            value=settings.CITATION_HIGH_THRESHOLD,
            step=0.05, key="cite_high",
        )
    with col2:
        st.number_input(
            "Medium relevance threshold",
            min_value=0.0, max_value=1.0,
            value=settings.CITATION_MEDIUM_THRESHOLD,
            step=0.05, key="cite_med",
        )

    st.markdown("#### Collection Search Weights")
    st.caption("Weights control how results from each collection influence final ranking.")
    from src.rag_engine import COLLECTION_CONFIG

    weight_cols = st.columns(3)
    for i, (coll, config) in enumerate(COLLECTION_CONFIG.items()):
        label = COLLECTION_LABELS.get(coll, coll)
        weight_cols[i % 3].slider(
            label, min_value=0.0, max_value=0.5,
            value=config["weight"], step=0.01,
            key=f"weight_{coll}",
        )

    # Clear conversation
    st.markdown("---")
    if st.button("Clear Conversation History", key="clear_history"):
        st.session_state.conversation_history = []
        st.session_state.last_evidence = None
        st.rerun()


# =====================================================================
# Main
# =====================================================================


def main():
    """Main application entry point."""
    ctx = safe_init()
    render_sidebar(ctx)

    st.markdown("# Imaging Intelligence Agent")
    st.caption("Multi-collection RAG engine for CT, MRI, X-ray, and medical imaging AI")

    tab_ask, tab_comp, tab_workflow, tab_reports, tab_settings = st.tabs([
        "Ask", "Comparative", "Workflow Demo", "Reports", "Settings",
    ])

    with tab_ask:
        render_ask_tab(ctx)

    with tab_comp:
        render_comparative_tab(ctx)

    with tab_workflow:
        render_workflow_tab(ctx)

    with tab_reports:
        render_reports_tab(ctx)

    with tab_settings:
        render_settings_tab(ctx)


if __name__ == "__main__":
    main()
