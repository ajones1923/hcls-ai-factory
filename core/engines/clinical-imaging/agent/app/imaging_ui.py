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
# Tab Module Imports (graceful fallbacks)
# =====================================================================

try:
    from app.tabs import (
        evidence_explorer,
        workflow_runner,
        protocol_advisor,
        device_ecosystem,
        dose_intelligence,
        reports_export,
        patient_360,
        benchmarks,
        image_gallery,
        analytics,
    )
except ImportError:
    from tabs import (
        evidence_explorer,
        workflow_runner,
        protocol_advisor,
        device_ecosystem,
        dose_intelligence,
        reports_export,
        patient_360,
        benchmarks,
        image_gallery,
        analytics,
    )

# =====================================================================
# Page Configuration
# =====================================================================

st.set_page_config(
    page_title="Clinical Imaging Engine",
    page_icon="🩻",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "About": (
            "Clinical Imaging Engine -- Engine 4 of the HCLS AI Factory\n\n"
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
    "imaging_radiomics": "Radiomics Features",
    "imaging_reports": "Radiology Reports",
}

NIM_STATUS_COLORS = {
    "available": "🟢",
    "mock": "🟡",
    "unavailable": "🔴",
}


def init_session_state():
    """Initialize Streamlit session state variables."""
    defaults = {
        "conversation_history": [],
        "last_evidence": None,
        "last_workflow_result": None,
        "last_query_result": None,
        "active_collections": list(COLLECTION_LABELS.keys()),
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
# Sidebar
# =====================================================================


def render_sidebar(ctx: Optional[Dict]):
    """Render the sidebar with NIM status, collection stats, and filters."""
    with st.sidebar:
        st.markdown("## Clinical Imaging Engine")
        st.caption("Multi-collection RAG for medical imaging")

        # -- Guided Tour: Start Here --
        with st.expander("🎯 **Start Here — Demo Guide**", expanded=not st.session_state.get("_tour_dismissed", False)):
            st.markdown("""
**Recommended Demo Flow (28 min):**

**1️⃣ Workflow Runner** *(5 min)*
> Run 3 pre-loaded clinical cases:
> - 🧠 Emergency stroke (ICH detection)
> - 🫁 Lung cancer screening (nodule tracking)
> - ❤️ Cardiac workup (coronary stenosis)

**2️⃣ Image Gallery** *(3 min)*
> AI-annotated medical imaging showcase:
> - Toggle between raw ↔ AI-analyzed views
> - Scroll through 3D volume slices
> - See detection overlays with measurements

**3️⃣ Evidence Explorer** *(3 min)*
> Multi-collection RAG search:
> - Try: *"What AI models detect hemorrhage?"*
> - Try: *"Compare CT vs MRI for stroke"*

**4️⃣ Protocol Advisor** *(2 min)*
> Patient-specific imaging protocols:
> - Enter: *"chest pain, rule out PE"*
> - See AI-optimized dose reduction

**5️⃣ Dose Intelligence** *(2 min)*
> 36% average dose reduction with DLIR

**6️⃣ Device Ecosystem** *(1 min)*
> 50 FDA-cleared AI devices catalog

**7️⃣ Reports & Export** *(2 min)*
> Generate PDF with NVIDIA branding
> Export FHIR R4 DiagnosticReport

**8️⃣ Patient 360** *(3 min)*
> Cross-modal imaging ↔ genomics integration
> Interactive gene-finding network graph

**9️⃣ Analytics** *(3 min)*
> Population-level imaging analytics
> GPU-accelerated with NVIDIA RAPIDS

**🔟 Benchmarks** *(1 min)*
> AI model performance validation
            """)

            if st.button("✅ Got it, dismiss tour", use_container_width=True):
                st.session_state["_tour_dismissed"] = True
                st.rerun()

        if ctx is None:
            st.warning("Engine not initialized")
            return

        settings = ctx["settings"]
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
        modalities = ["All", "ct", "mri", "xray", "cxr", "ultrasound", "pet", "pet_ct", "mammography"]
        body_regions = [
            "All", "head", "neck", "chest", "abdomen", "pelvis",
            "spine", "extremity", "brain", "cardiac", "breast",
        ]
        st.session_state.modality_filter = st.selectbox("Modality", modalities, index=0)
        st.session_state.body_region_filter = st.selectbox("Body Region", body_regions, index=0)
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

        st.divider()
        st.markdown("### \U0001f3af Demo Mode")
        if st.button("Load Demo Patient", key="load_demo"):
            import sys as _sys
            import os as _os
            _lib_path = _os.environ.get("HCLS_LIB_PATH", "/app/lib")
            if _lib_path not in _sys.path:
                _sys.path.insert(0, _lib_path)
            from hcls_common.demo_data import DEMO_IMAGING
            st.session_state["modality_filter"] = DEMO_IMAGING["modality"]
            st.session_state["body_region_filter"] = DEMO_IMAGING["body_region"]
            st.session_state["demo_question"] = DEMO_IMAGING["question"]
            st.toast("\u2705 Demo patient loaded! Ask the demo question in Evidence Explorer tab.", icon="\U0001f3af")
            st.rerun()

        st.divider()
        st.markdown("### DICOM Viewer")
        ohif_url = getattr(settings, "OHIF_URL", "http://localhost:8526")
        st.link_button("Open OHIF Viewer", ohif_url)
        try:
            import requests
            resp = requests.get(
                f"{settings.ORTHANC_URL}/studies",
                auth=(settings.ORTHANC_USERNAME, settings.ORTHANC_PASSWORD) if settings.ORTHANC_PASSWORD else None,
                timeout=2,
            )
            if resp.status_code == 200:
                st.caption(f"Orthanc: {len(resp.json())} studies")
        except Exception:
            st.caption("Orthanc: offline")


# =====================================================================
# Main
# =====================================================================


def main():
    """Main application entry point."""
    # Inject custom theme
    try:
        from app.custom_theme import inject_theme, render_hero_header, render_clinical_disclaimer, render_footer
    except ImportError:
        from custom_theme import inject_theme, render_hero_header, render_clinical_disclaimer, render_footer
    inject_theme()

    ctx = safe_init()

    # Store engine in session state for tab modules to access
    if ctx:
        st.session_state["rag_engine"] = ctx["engine"]
        st.session_state["nim_manager"] = ctx["nim_manager"]
        st.session_state["settings"] = ctx["settings"]

    render_sidebar(ctx)

    # Branded hero header
    render_hero_header()

    # Styled clinical disclaimer
    render_clinical_disclaimer()

    # Main tabs
    tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10 = st.tabs([
        "Evidence Explorer",
        "Workflow Runner",
        "Image Gallery",
        "Protocol Advisor",
        "Device & AI Ecosystem",
        "Dose Intelligence",
        "Reports & Export",
        "Patient 360",
        "Analytics",
        "Benchmarks & Validation",
    ])

    with tab1:
        evidence_explorer.render()

    with tab2:
        workflow_runner.render()

    with tab3:
        image_gallery.render()

    with tab4:
        protocol_advisor.render()

    with tab5:
        device_ecosystem.render()

    with tab6:
        dose_intelligence.render()

    with tab7:
        reports_export.render()

    with tab8:
        patient_360.render()

    with tab9:
        analytics.render()

    with tab10:
        benchmarks.render()

    # Footer
    render_footer()


if __name__ == "__main__":
    main()
