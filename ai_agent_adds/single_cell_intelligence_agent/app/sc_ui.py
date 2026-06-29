"""Single-Cell Intelligence Agent -- 5-Tab Streamlit UI.

NVIDIA dark-themed single-cell genomics clinical decision support
interface with RAG-powered queries, TME profiling, workflow runners,
and real-time dashboard monitoring.

Usage:
    streamlit run app/sc_ui.py --server.port 8130

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

API_BASE = os.environ.get("SC_API_BASE", "http://localhost:8540")

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
    page_title="Single-Cell Intelligence Agent",
    page_icon="\U0001F9EC",
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
    """GET request to single-cell API with error handling."""
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
    """POST request to single-cell API with error handling."""
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
        <h2 style="color: {NVIDIA_THEME['accent']}; margin: 0;">Single-Cell Intelligence</h2>
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

tab_dashboard, tab_explorer, tab_tme, tab_workflows, tab_reports = st.tabs([
    "Dashboard",
    "Evidence Explorer",
    "TME Profiler",
    "Workflow Runner",
    "Reports & Export",
])


# =====================================================================
# Tab 1: Dashboard
# =====================================================================

with tab_dashboard:
    st.header("Single-Cell Intelligence Dashboard")

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
            st.metric("Workflows", health.get("workflows", 10))
    else:
        st.info("Connect to the API to view dashboard metrics.")

    st.markdown("---")

    # Cell types overview
    st.subheader("Cell Type Catalogue")
    cell_types = api_get("/v1/sc/cell-types")
    if cell_types:
        ct_list = cell_types.get("cell_types", [])
        cols = st.columns(2)
        compartments = {}
        for ct in ct_list:
            comp = ct.get("compartment", "other")
            compartments.setdefault(comp, []).append(ct)

        for i, (comp, cts) in enumerate(compartments.items()):
            with cols[i % 2]:
                with st.expander(f"{comp.title()} ({len(cts)} types)", expanded=False):
                    for ct in cts:
                        markers = ", ".join(ct.get("canonical_markers", [])[:5])
                        st.text(f"  {ct.get('name', 'N/A')}: {markers}")

    # GPU status placeholder
    st.subheader("GPU Status")
    st.info("GPU monitoring available when NVIDIA drivers are detected.")

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
    st.write("RAG-powered single-cell biology Q&A across all knowledge collections.")

    # Domain selector
    domain_options = [
        "auto", "transcriptomics", "immunology", "oncology",
        "spatial", "multiomics", "pharmacology", "cell_therapy", "general",
    ]
    selected_domain = st.selectbox(
        "Domain Focus",
        domain_options,
        index=0,
        help="Select a domain to guide the query, or leave as 'auto'.",
    )

    # Query input
    question = st.text_area(
        "Single-Cell Biology Question",
        placeholder="e.g., What marker genes distinguish M1 vs M2 macrophages in the tumor microenvironment?",
        height=100,
    )

    col_topk, col_guidelines = st.columns(2)
    with col_topk:
        top_k = st.slider("Evidence passages (top_k)", 1, 20, 5)
    with col_guidelines:
        include_guidelines = st.checkbox("Include guideline citations", value=True)

    if st.button("Search", key="explorer_search"):
        if question.strip():
            with st.spinner("Searching single-cell knowledge base..."):
                payload = {
                    "question": question.strip(),
                    "top_k": top_k,
                    "include_guidelines": include_guidelines,
                }
                if selected_domain != "auto":
                    payload["domain"] = selected_domain

                result = api_post("/v1/sc/query", payload)

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
# Tab 3: TME Profiler
# =====================================================================

with tab_tme:
    st.header("Tumor Microenvironment Profiler")
    st.write("Upload tumor cell composition data to classify the microenvironment and predict therapy response.")

    tumor_type = st.selectbox(
        "Tumor Type",
        ["NSCLC", "Melanoma", "Breast", "Colorectal", "Ovarian",
         "Pancreatic", "Glioblastoma", "RCC", "HNSCC", "Other"],
        key="tme_tumor_type",
    )

    st.subheader("Cell Type Proportions")
    st.write("Enter estimated proportions (0.0 - 1.0) for each cell type:")

    c1, c2, c3 = st.columns(3)
    proportions = {}
    with c1:
        proportions["cd8_t_cells"] = st.number_input("CD8+ T Cells", 0.0, 1.0, 0.10, 0.01, key="tme_cd8")
        proportions["cd4_t_cells"] = st.number_input("CD4+ T Cells", 0.0, 1.0, 0.08, 0.01, key="tme_cd4")
        proportions["tregs"] = st.number_input("Tregs", 0.0, 1.0, 0.03, 0.01, key="tme_treg")
        proportions["nk_cells"] = st.number_input("NK Cells", 0.0, 1.0, 0.02, 0.01, key="tme_nk")
    with c2:
        proportions["macrophages_m1"] = st.number_input("Macrophages (M1)", 0.0, 1.0, 0.05, 0.01, key="tme_m1")
        proportions["macrophages_m2"] = st.number_input("Macrophages (M2)", 0.0, 1.0, 0.08, 0.01, key="tme_m2")
        proportions["dendritic_cells"] = st.number_input("Dendritic Cells", 0.0, 1.0, 0.02, 0.01, key="tme_dc")
        proportions["b_cells"] = st.number_input("B Cells", 0.0, 1.0, 0.03, 0.01, key="tme_b")
    with c3:
        proportions["fibroblasts"] = st.number_input("Fibroblasts/CAFs", 0.0, 1.0, 0.15, 0.01, key="tme_fib")
        proportions["endothelial"] = st.number_input("Endothelial", 0.0, 1.0, 0.05, 0.01, key="tme_endo")
        proportions["malignant"] = st.number_input("Malignant", 0.0, 1.0, 0.35, 0.01, key="tme_mal")
        proportions["other"] = st.number_input("Other", 0.0, 1.0, 0.04, 0.01, key="tme_other")

    immune_markers_text = st.text_area(
        "Immune Marker Expression (JSON, optional)",
        placeholder='{"PD_L1": 0.45, "CTLA4": 0.12, "LAG3": 0.08, "TIM3": 0.15}',
        height=80,
        key="tme_markers",
    )

    clinical_notes = st.text_area("Clinical Notes (optional)", key="tme_notes", height=80)

    if st.button("Profile TME", key="tme_run"):
        immune_markers = None
        if immune_markers_text.strip():
            try:
                immune_markers = json.loads(immune_markers_text)
            except json.JSONDecodeError:
                st.error("Invalid JSON in Immune Marker Expression field.")
                immune_markers = None

        with st.spinner("Profiling tumor microenvironment..."):
            payload = {
                "cell_type_proportions": {k: v for k, v in proportions.items() if v > 0},
                "tumor_type": tumor_type.lower(),
            }
            if immune_markers:
                payload["immune_markers"] = immune_markers
            if clinical_notes:
                payload["clinical_notes"] = clinical_notes

            result = api_post("/v1/sc/tme-profile", payload)

        if result:
            st.subheader("TME Classification")
            rc1, rc2, rc3 = st.columns(3)
            with rc1:
                st.metric("Classification", result.get("tme_classification", "N/A").replace("_", " ").title())
            with rc2:
                st.metric("Immune Score", f"{result.get('immune_score', 0):.2f}")
            with rc3:
                st.metric("Stromal Score", f"{result.get('stromal_score', 0):.2f}")

            st.write(f"**Immune Phenotype:** {result.get('immune_phenotype', 'N/A')}")

            checkpoint = result.get("checkpoint_expression", {})
            if checkpoint:
                st.subheader("Checkpoint Expression")
                for k, v in checkpoint.items():
                    st.write(f"- **{k}:** {v}")

            therapy = result.get("therapy_prediction", {})
            if therapy:
                st.subheader("Therapy Predictions")
                for k, v in therapy.items():
                    st.write(f"- **{k}:** {v}")

            recs = result.get("recommendations", [])
            if recs:
                st.subheader("Recommendations")
                for rec in recs:
                    st.write(f"- {rec}")

            st.subheader("Full Results")
            st.json(result)


# =====================================================================
# Tab 4: Workflow Runner
# =====================================================================

with tab_workflows:
    st.header("Single-Cell Workflow Runner")
    st.write("Execute specialized single-cell analysis workflows.")

    workflow_choice = st.selectbox(
        "Select Workflow",
        [
            "Cell Type Annotation",
            "TME Profiling",
            "Drug Response Prediction",
            "Subclonal Analysis",
            "Spatial Niche Mapping",
            "Trajectory Inference",
            "Ligand-Receptor Interaction",
            "Biomarker Discovery",
            "CAR-T Target Validation",
            "Treatment Monitoring",
        ],
        key="workflow_selector",
    )

    workflow_endpoints = {
        "Cell Type Annotation": "/v1/sc/annotate",
        "TME Profiling": "/v1/sc/tme-profile",
        "Drug Response Prediction": "/v1/sc/drug-response",
        "Subclonal Analysis": "/v1/sc/subclonal",
        "Spatial Niche Mapping": "/v1/sc/spatial-niche",
        "Trajectory Inference": "/v1/sc/trajectory",
        "Ligand-Receptor Interaction": "/v1/sc/ligand-receptor",
        "Biomarker Discovery": "/v1/sc/biomarker",
        "CAR-T Target Validation": "/v1/sc/cart-validate",
        "Treatment Monitoring": "/v1/sc/treatment-monitor",
    }

    wf_payload = {}

    if workflow_choice == "Cell Type Annotation":
        c1, c2 = st.columns(2)
        with c1:
            markers_text = st.text_area(
                "Marker Genes (one per line)",
                placeholder="CD3D\nCD8A\nGZMB\nPRF1",
                key="wf_markers",
            )
            wf_payload["marker_genes"] = [m.strip() for m in markers_text.split("\n") if m.strip()]
            wf_payload["tissue_type"] = st.text_input("Tissue Type", placeholder="e.g., lung, PBMC, breast", key="wf_tissue")
        with c2:
            wf_payload["strategy"] = st.selectbox("Strategy", ["marker_based", "reference_based", "llm_based"], key="wf_strategy")
            wf_payload["num_cells"] = st.number_input("Number of Cells", 1, 1000000, 500, key="wf_ncells")
            wf_payload["cluster_id"] = st.text_input("Cluster ID (optional)", key="wf_cluster")
            wf_payload["top_n"] = st.slider("Max Cell Types", 1, 20, 5, key="wf_topn")

    elif workflow_choice == "Drug Response Prediction":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["drug_name"] = st.text_input("Drug Name", placeholder="e.g., Pembrolizumab", key="wf_drug")
            wf_payload["drug_class"] = st.selectbox("Drug Class", ["checkpoint_inhibitor", "tki", "chemotherapy", "antibody_drug_conjugate", "car_t", "other"], key="wf_drug_class")
            wf_payload["cell_type"] = st.text_input("Target Cell Type", key="wf_drug_ct")
        with c2:
            wf_payload["tumor_type"] = st.text_input("Tumor Type", key="wf_drug_tumor")
            alt_text = st.text_area("Genomic Alterations (one per line)", key="wf_drug_alt")
            wf_payload["genomic_alterations"] = [a.strip() for a in alt_text.split("\n") if a.strip()]
        wf_payload["clinical_notes"] = st.text_area("Clinical Notes (optional)", key="wf_drug_notes")

    elif workflow_choice == "Subclonal Analysis":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["cell_count"] = st.number_input("Total Cells", 1, 1000000, 5000, key="wf_sub_cells")
            wf_payload["tumor_type"] = st.text_input("Tumor Type", key="wf_sub_tumor")
        with c2:
            wf_payload["target_antigen"] = st.text_input("Target Antigen (for escape analysis)", placeholder="e.g., CD19, HER2", key="wf_sub_target")
        wf_payload["clinical_notes"] = st.text_area("Clinical Notes (optional)", key="wf_sub_notes")

    elif workflow_choice == "Spatial Niche Mapping":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["platform"] = st.selectbox("Platform", ["visium", "merfish", "cosmx", "xenium", "slide_seq", "stereo_seq", "seqfish"], key="wf_sp_plat")
            wf_payload["tissue_type"] = st.text_input("Tissue Type", key="wf_sp_tissue")
        with c2:
            ct_text = st.text_area("Cell Types (JSON: spot_id -> cell_type)", placeholder='{"spot_1": "T_cell", "spot_2": "Macrophage"}', key="wf_sp_ct")
            if ct_text.strip():
                try:
                    wf_payload["cell_types"] = json.loads(ct_text)
                except json.JSONDecodeError:
                    st.error("Invalid JSON for cell types.")
        wf_payload["clinical_notes"] = st.text_area("Clinical Notes (optional)", key="wf_sp_notes")

    elif workflow_choice == "Trajectory Inference":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["method"] = st.selectbox("Method", ["monocle3", "paga", "rna_velocity", "scvelo", "cellrank", "palantir"], key="wf_traj_method")
            wf_payload["root_cell_type"] = st.text_input("Root Cell Type", placeholder="e.g., HSC, progenitor", key="wf_traj_root")
        with c2:
            ct_text = st.text_area("Cell Types (one per line)", placeholder="HSC\nProgenitor\nMature_T_cell", key="wf_traj_ct")
            wf_payload["cell_types"] = [c.strip() for c in ct_text.split("\n") if c.strip()]
            wf_payload["tissue_type"] = st.text_input("Tissue Type", key="wf_traj_tissue")

    elif workflow_choice == "Ligand-Receptor Interaction":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["source_cell_type"] = st.text_input("Source Cell Type (ligand)", key="wf_lr_src")
            wf_payload["target_cell_type"] = st.text_input("Target Cell Type (receptor)", key="wf_lr_tgt")
        with c2:
            wf_payload["database"] = st.selectbox("Database", ["cellphonedb", "nichenet", "cellchat", "celltalkdb"], key="wf_lr_db")
            wf_payload["tissue_type"] = st.text_input("Tissue Type", key="wf_lr_tissue")
            ct_text = st.text_area("All Cell Types (one per line)", key="wf_lr_all_ct")
            wf_payload["cell_types"] = [c.strip() for c in ct_text.split("\n") if c.strip()]

    elif workflow_choice == "Biomarker Discovery":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["cell_type"] = st.text_input("Cell Type of Interest", key="wf_bio_ct")
            wf_payload["condition_a"] = st.text_input("Condition A (e.g., responder)", key="wf_bio_a")
            wf_payload["condition_b"] = st.text_input("Condition B (e.g., non-responder)", key="wf_bio_b")
        with c2:
            wf_payload["tumor_type"] = st.text_input("Tumor Type", key="wf_bio_tumor")
            wf_payload["min_log2fc"] = st.number_input("Min Log2FC", 0.0, 10.0, 1.0, 0.1, key="wf_bio_fc")
            wf_payload["max_pval"] = st.number_input("Max Adj P-value", 0.001, 1.0, 0.05, 0.005, key="wf_bio_pval")

    elif workflow_choice == "CAR-T Target Validation":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["target_gene"] = st.text_input("Target Gene", placeholder="e.g., CD19, BCMA, HER2", key="wf_cart_gene")
            wf_payload["tumor_type"] = st.text_input("Tumor Type", key="wf_cart_tumor")
        with c2:
            normal_text = st.text_area(
                "Normal Tissue Expression (JSON)",
                placeholder='{"heart": 0.01, "lung": 0.05, "liver": 0.02}',
                key="wf_cart_normal",
            )
            if normal_text.strip():
                try:
                    wf_payload["normal_tissue_expression"] = json.loads(normal_text)
                except json.JSONDecodeError:
                    st.error("Invalid JSON for normal tissue expression.")
        wf_payload["clinical_notes"] = st.text_area("Clinical Notes (optional)", key="wf_cart_notes")

    elif workflow_choice == "Treatment Monitoring":
        c1, c2 = st.columns(2)
        with c1:
            wf_payload["treatment"] = st.text_input("Treatment Regimen", key="wf_tm_treat")
            wf_payload["tumor_type"] = st.text_input("Tumor Type", key="wf_tm_tumor")
            wf_payload["target_antigen"] = st.text_input("Target Antigen", key="wf_tm_antigen")
        with c2:
            baseline_text = st.text_area(
                "Baseline Composition (JSON)",
                placeholder='{"cd8_t": 0.10, "malignant": 0.60, "treg": 0.05}',
                key="wf_tm_base",
            )
            if baseline_text.strip():
                try:
                    wf_payload["baseline_composition"] = json.loads(baseline_text)
                except json.JSONDecodeError:
                    st.error("Invalid JSON for baseline composition.")
            current_text = st.text_area(
                "Current Composition (JSON)",
                placeholder='{"cd8_t": 0.25, "malignant": 0.30, "treg": 0.08}',
                key="wf_tm_curr",
            )
            if current_text.strip():
                try:
                    wf_payload["current_composition"] = json.loads(current_text)
                except json.JSONDecodeError:
                    st.error("Invalid JSON for current composition.")

    if st.button("Run Workflow", key="wf_run"):
        endpoint = workflow_endpoints.get(workflow_choice, "")
        if endpoint:
            # Clean payload: remove empty strings and None values
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
    st.write("Generate and export structured single-cell analysis reports.")

    report_type = st.selectbox(
        "Report Type",
        [
            "cell_type_annotation", "tme_profile", "drug_response",
            "subclonal_analysis", "spatial_niche", "trajectory",
            "ligand_receptor", "biomarker", "cart_validation",
            "treatment_monitoring", "general",
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
        value='{\n  "summary": "Single-cell analysis results"\n}',
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
                    file_name=f"sc_{report_type}_{result.get('report_id', 'report')}{ext}",
                    mime=mime,
                )

    # Available formats reference
    st.markdown("---")
    st.subheader("Supported Formats")
    formats = api_get("/v1/reports/formats")
    if formats:
        for fmt in formats.get("formats", []):
            st.write(f"- **{fmt.get('name', 'N/A')}** ({fmt.get('extension', '')}) -- {fmt.get('description', '')}")
