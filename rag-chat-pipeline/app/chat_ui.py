"""
Streamlit Chat Interface for Genomic Evidence RAG.
Includes VCF preview, LLM performance metrics, and Target Hypothesis management.
"""
import os
import streamlit as st
import time
import json
from typing import Optional, List, Dict, Any
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.rag_engine import create_rag_engine, RAGEngine
from src.target_hypothesis import (
    TargetHypothesis,
    TargetHypothesisManager,
    create_hypothesis_from_chat
)
from src.knowledge import KNOWLEDGE_CONNECTIONS, get_gene_reference_data
from config.settings import settings


# ═══════════════════════════════════════════════════════════════════════════════
# KNOWLEDGE CONNECTION DATABASE now imported from src/knowledge.py
# Maps: Gene → Protein → Pathway → Disease → Drug
# ═══════════════════════════════════════════════════════════════════════════════


@st.cache_data(ttl=300)  # Cache for 5 minutes
def get_variant_stats() -> Dict[str, int]:
    """Get variant counts from Milvus database."""
    try:
        from pymilvus import connections, Collection
        connections.connect("default", host=os.environ.get("MILVUS_HOST", "localhost"), port=os.environ.get("MILVUS_PORT", "19530"))
        collection = Collection("genomic_evidence")
        collection.load()

        # Get VCP variant count
        vcp_results = collection.query(
            expr='gene == "VCP"',
            output_fields=["id"],
            limit=1000
        )
        vcp_count = len(vcp_results)

        # Get total neurodegeneration gene count
        neuro_genes = ["VCP", "C9orf72", "GRN", "TBK1", "FUS", "TARDBP", "MAPT", "SOD1"]
        total_neuro = 0
        for gene in neuro_genes:
            results = collection.query(
                expr=f'gene == "{gene}"',
                output_fields=["id"],
                limit=1000
            )
            total_neuro += len(results)

        connections.disconnect("default")
        return {"vcp": vcp_count, "neuro_total": total_neuro}
    except Exception:
        # Fallback to reasonable estimates if Milvus unavailable
        return {"vcp": 13, "neuro_total": 35}


# Page configuration
st.set_page_config(
    page_title="Genomic RAG | NVIDIA DGX Spark Demo",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS for polished 10/10 UI - Light Theme
st.markdown("""
<style>
    /* Global improvements - Light Theme */
    .stApp {
        background: linear-gradient(180deg, #ffffff 0%, #f8f9fa 100%);
    }

    /* Sidebar styling - Light */
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #f8f9fa 0%, #ffffff 100%);
        border-right: 1px solid #e1e4e8;
    }
    [data-testid="stSidebar"] .stMarkdown h3 {
        color: #2ea44f !important;
        font-size: 1.1rem;
        letter-spacing: 0.5px;
    }

    /* Header styling - Light with NVIDIA green accent */
    .main-header {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        padding: 1.5rem 2rem;
        border-radius: 12px;
        margin-bottom: 1.5rem;
        border: 1px solid #30363d;
        border-left: 4px solid #76b900;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.15);
    }
    .main-header h1 {
        color: #ffffff;
        margin: 0;
        font-size: 1.8rem;
        font-weight: 600;
        letter-spacing: -0.5px;
    }
    .main-header p {
        color: #a0aec0;
        margin: 0.5rem 0 0 0;
        font-size: 0.95rem;
    }

    /* Demo badge with glow */
    .demo-badge {
        background: linear-gradient(135deg, #76b900 0%, #5a8c00 100%);
        color: #000;
        padding: 0.3rem 0.9rem;
        border-radius: 20px;
        font-size: 0.7rem;
        font-weight: 700;
        display: inline-block;
        margin-left: 1rem;
        text-transform: uppercase;
        letter-spacing: 1px;
        box-shadow: 0 0 15px rgba(118, 185, 0, 0.4);
        animation: pulse-glow 2s ease-in-out infinite;
    }
    @keyframes pulse-glow {
        0%, 100% { box-shadow: 0 0 15px rgba(118, 185, 0, 0.4); }
        50% { box-shadow: 0 0 25px rgba(118, 185, 0, 0.6); }
    }

    /* Stats metrics styling - Light */
    [data-testid="stMetricValue"] {
        font-size: 1.6rem !important;
        font-weight: 700 !important;
        color: #1f2937 !important;
    }
    [data-testid="stMetricLabel"] {
        color: #6b7280 !important;
        font-size: 0.85rem !important;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }
    [data-testid="stMetricDelta"] {
        font-size: 0.75rem !important;
    }

    /* Demo question buttons */
    .demo-btn {
        background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);
        border: 1px solid #e1e4e8;
        border-radius: 8px;
        padding: 0.6rem 1rem;
        margin: 0.3rem 0;
        cursor: pointer;
        transition: all 0.2s ease;
        width: 100%;
        text-align: left;
    }
    .demo-btn:hover {
        border-color: #76b900;
        background: linear-gradient(135deg, #f0fff4 0%, #ffffff 100%);
        transform: translateX(4px);
        box-shadow: 0 2px 8px rgba(118, 185, 0, 0.2);
    }

    /* Sidebar buttons - Light */
    [data-testid="stSidebar"] .stButton > button {
        background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%) !important;
        border: 1px solid #e1e4e8 !important;
        border-radius: 8px !important;
        color: #1f2937 !important;
        font-size: 0.85rem !important;
        padding: 0.6rem 1rem !important;
        transition: all 0.2s ease !important;
        text-align: left !important;
    }
    [data-testid="stSidebar"] .stButton > button:hover {
        border-color: #76b900 !important;
        background: linear-gradient(135deg, #f0fff4 0%, #ffffff 100%) !important;
        transform: translateX(4px) !important;
        box-shadow: 0 2px 8px rgba(118, 185, 0, 0.2) !important;
    }

    /* Chat container - Light */
    .stChatMessage {
        border-radius: 12px;
        background: #ffffff;
        border: 1px solid #e1e4e8;
    }

    /* Chat input - Light */
    .stChatInput > div {
        border-radius: 12px !important;
        border: 1px solid #e1e4e8 !important;
        background: #ffffff !important;
    }
    .stChatInput > div:focus-within {
        border-color: #76b900 !important;
        box-shadow: 0 0 0 2px rgba(118, 185, 0, 0.2) !important;
    }

    /* Tabs styling - Light */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
        background: transparent;
    }
    .stTabs [data-baseweb="tab"] {
        background: #f8f9fa;
        border-radius: 8px 8px 0 0;
        border: 1px solid #e1e4e8;
        border-bottom: none;
        color: #6b7280;
        padding: 0.5rem 1.5rem;
    }
    .stTabs [aria-selected="true"] {
        background: #ffffff !important;
        border-color: #76b900 !important;
        color: #1f2937 !important;
    }

    /* Expander styling - Light */
    .streamlit-expanderHeader {
        background: #f8f9fa !important;
        border-radius: 8px !important;
        border: 1px solid #e1e4e8 !important;
    }
    .streamlit-expanderHeader:hover {
        border-color: #76b900 !important;
    }

    /* Section header - Light */
    .section-header {
        color: #2ea44f;
        font-size: 0.75rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 1.5px;
        margin: 1rem 0 0.5rem 0;
        padding-bottom: 0.3rem;
        border-bottom: 1px solid #e1e4e8;
    }

    /* Model selector badge */
    .model-badge {
        display: inline-block;
        padding: 0.15rem 0.5rem;
        border-radius: 4px;
        font-size: 0.65rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }
    .model-badge.cloud {
        background: #1f6feb;
        color: #ffffff;
    }
    .model-badge.local {
        background: #2ea44f;
        color: #ffffff;
    }

    /* Welcome card - Light */
    .welcome-card {
        background: linear-gradient(135deg, #f8f9fa 0%, #ffffff 100%);
        border: 1px solid #e1e4e8;
        border-radius: 12px;
        padding: 1.5rem;
        margin: 1rem 0;
    }
    .welcome-card h3 {
        color: #2ea44f;
        margin-bottom: 0.5rem;
    }

    /* Divider - Light */
    .sidebar-divider {
        border: none;
        border-top: 1px solid #e1e4e8;
        margin: 1rem 0;
    }

    /* Hide default Streamlit elements */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}

    /* Smooth scrolling */
    html {
        scroll-behavior: smooth;
    }

    /* ═══════════════════════════════════════════════════════════
       KNOWLEDGE CONNECTION VISUALIZATION (Clinker-style)
       ═══════════════════════════════════════════════════════════ */
    .knowledge-container {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        border-radius: 12px;
        padding: 1.5rem;
        margin: 1rem 0;
        border: 1px solid #30363d;
        border-left: 4px solid #76b900;
    }
    .knowledge-header {
        color: #76b900;
        font-size: 1rem;
        font-weight: 600;
        margin-bottom: 1rem;
        display: flex;
        align-items: center;
        gap: 0.5rem;
    }
    .knowledge-path {
        display: flex;
        align-items: center;
        justify-content: space-between;
        flex-wrap: wrap;
        gap: 0.5rem;
        padding: 1rem 0;
    }
    .knowledge-node {
        background: linear-gradient(135deg, #ffffff 0%, #f0f0f0 100%);
        border-radius: 8px;
        padding: 0.75rem 1rem;
        text-align: center;
        min-width: 120px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.3);
        transition: all 0.2s ease;
    }
    .knowledge-node:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(118, 185, 0, 0.4);
    }
    .knowledge-node.variant {
        background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%);
        border: 2px solid #f59e0b;
    }
    .knowledge-node.gene {
        background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%);
        border: 2px solid #3b82f6;
    }
    .knowledge-node.protein {
        background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%);
        border: 2px solid #10b981;
    }
    .knowledge-node.pathway {
        background: linear-gradient(135deg, #e0e7ff 0%, #c7d2fe 100%);
        border: 2px solid #6366f1;
    }
    .knowledge-node.disease {
        background: linear-gradient(135deg, #fee2e2 0%, #fecaca 100%);
        border: 2px solid #ef4444;
    }
    .knowledge-node.drug {
        background: linear-gradient(135deg, #76b900 0%, #5a8c00 100%);
        border: 2px solid #4a7c00;
        color: white;
    }
    .node-label {
        font-size: 0.65rem;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        color: #6b7280;
        margin-bottom: 0.25rem;
    }
    .knowledge-node.drug .node-label {
        color: rgba(255,255,255,0.8);
    }
    .node-value {
        font-size: 0.85rem;
        font-weight: 600;
        color: #1f2937;
    }
    .knowledge-node.drug .node-value {
        color: white;
    }
    .knowledge-arrow {
        color: #76b900;
        font-size: 1.5rem;
        font-weight: bold;
    }
    .knowledge-details {
        margin-top: 1rem;
        padding-top: 1rem;
        border-top: 1px solid #30363d;
    }
    .knowledge-detail-row {
        display: flex;
        justify-content: space-between;
        padding: 0.5rem 0;
        border-bottom: 1px solid rgba(255,255,255,0.1);
    }
    .knowledge-detail-label {
        color: #a0aec0;
        font-size: 0.8rem;
    }
    .knowledge-detail-value {
        color: #ffffff;
        font-size: 0.8rem;
        font-weight: 500;
    }
    .druggable-badge {
        background: #76b900;
        color: black;
        padding: 0.2rem 0.5rem;
        border-radius: 4px;
        font-size: 0.7rem;
        font-weight: 700;
    }
    .not-druggable-badge {
        background: #6b7280;
        color: white;
        padding: 0.2rem 0.5rem;
        border-radius: 4px;
        font-size: 0.7rem;
    }
</style>
""", unsafe_allow_html=True)


# Available models for selection
AVAILABLE_MODELS = {
    # Local models (Ollama) - Optimized options first
    "llama3.1:8b-fast-ctx2k": {"name": "Llama 8B Fast", "description": "Optimized Q4_0, 2K context", "size": "4.3GB", "provider": "ollama"},
    "llama3.1:8b-instruct-q4_0": {"name": "Llama 8B Q4_0", "description": "Smaller quantization", "size": "4.3GB", "provider": "ollama"},
    "llama3.1:8b": {"name": "Llama 3.1 8B", "description": "Standard Q4_K_M", "size": "4.9GB", "provider": "ollama"},
    "llama3.1:70b": {"name": "Llama 3.1 70B", "description": "Best local quality", "size": "40GB", "provider": "ollama"},
    # Cloud models (Anthropic)
    "claude-sonnet-4-20250514": {"name": "Claude Sonnet 4", "description": "Fast cloud (2-5s)", "size": "Cloud", "provider": "anthropic"},
    "claude-opus-4-20250514": {"name": "Claude Opus 4", "description": "Best quality (5-10s)", "size": "Cloud", "provider": "anthropic"},
}


def load_shared_model() -> tuple:
    """Load model and provider from shared file (set by portal)."""
    model_file = settings.DATA_DIR / 'current_model.json'
    try:
        if model_file.exists():
            with open(model_file, 'r') as f:
                data = json.load(f)
                model = data.get('model', settings.LLM_MODEL)
                provider = data.get('provider', settings.LLM_PROVIDER)
                return model, provider
    except (json.JSONDecodeError, OSError, IOError):
        pass
    return settings.LLM_MODEL, settings.LLM_PROVIDER


def get_provider_for_model(model: str) -> str:
    """Get the provider for a given model."""
    if model in AVAILABLE_MODELS:
        return AVAILABLE_MODELS[model].get('provider', 'ollama')
    elif model.startswith('claude-'):
        return 'anthropic'
    elif model.startswith('gpt-'):
        return 'openai'
    return 'ollama'


def get_rag_engine(model_name: str = None) -> RAGEngine:
    """Initialize RAG engine with specified model."""
    # Priority: explicit model_name > session state > shared file > settings
    if model_name:
        model = model_name
        provider = get_provider_for_model(model)
    elif st.session_state.get('selected_model'):
        model = st.session_state.get('selected_model')
        provider = get_provider_for_model(model)
    else:
        model, provider = load_shared_model()

    # Check if we need to reinitialize
    cache_key = f"rag_engine_{model}"
    if cache_key not in st.session_state:
        st.session_state[cache_key] = create_rag_engine(
            milvus_host=settings.MILVUS_HOST,
            milvus_port=settings.MILVUS_PORT,
            collection_name=settings.MILVUS_COLLECTION,
            embedding_model=settings.EMBEDDING_MODEL,
            llm_provider=provider,
            llm_model=model,
            top_k=settings.RAG_TOP_K,
            score_threshold=settings.RAG_SCORE_THRESHOLD,
        )
    return st.session_state[cache_key]


@st.cache_resource
def get_target_manager() -> TargetHypothesisManager:
    """Initialize target hypothesis manager (cached)."""
    return TargetHypothesisManager(storage_dir=settings.DATA_DIR / "targets")


def get_vcf_preview(vcf_path: str, limit: int = 50):
    """Get VCF preview data."""
    variants = []
    try:
        try:
            from cyvcf2 import VCF
            vcf = VCF(vcf_path)
            for i, record in enumerate(vcf):
                if i >= limit:
                    break
                variants.append({
                    'CHROM': record.CHROM,
                    'POS': record.POS,
                    'REF': record.REF[:15] + ('...' if len(record.REF) > 15 else ''),
                    'ALT': ','.join([str(a)[:15] + ('...' if len(str(a)) > 15 else '') for a in record.ALT]),
                    'QUAL': f"{record.QUAL:.1f}" if record.QUAL else 'N/A',
                    'FILTER': record.FILTER or 'PASS',
                })
            vcf.close()
        except ImportError:
            import gzip
            opener = gzip.open if vcf_path.endswith('.gz') else open
            with opener(vcf_path, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    if len(variants) >= limit:
                        break
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        variants.append({
                            'CHROM': parts[0],
                            'POS': int(parts[1]),
                            'REF': parts[3][:15] + ('...' if len(parts[3]) > 15 else ''),
                            'ALT': parts[4][:15] + ('...' if len(parts[4]) > 15 else ''),
                            'QUAL': parts[5] if parts[5] != '.' else 'N/A',
                            'FILTER': parts[6] if parts[6] != '.' else 'PASS',
                        })
    except Exception as e:
        st.error(f"Error reading VCF: {e}")
        return []

    return variants


def render_sidebar():
    """Render sidebar with settings and info."""
    with st.sidebar:
        # Header with branding
        st.markdown("""
        <div style="text-align: center; padding: 0.5rem 0 1rem 0;">
            <span style="font-size: 1.5rem;">🧬</span>
            <span style="font-size: 1.1rem; font-weight: 600; color: #76b900; margin-left: 0.3rem;">RAG Config</span>
        </div>
        """, unsafe_allow_html=True)

        # ═══════════════════════════════════════════════════════════
        # SECTION 1: LLM Model Selection
        # ═══════════════════════════════════════════════════════════
        st.markdown('<div class="section-header">LLM Model</div>', unsafe_allow_html=True)
        current_model = st.session_state.get('selected_model', settings.LLM_MODEL)

        # Group models by provider
        cloud_models = {k: v for k, v in AVAILABLE_MODELS.items() if v['provider'] == 'anthropic'}
        local_models = {k: v for k, v in AVAILABLE_MODELS.items() if v['provider'] == 'ollama'}
        model_options = list(cloud_models.keys()) + list(local_models.keys())

        selected_model = st.selectbox(
            "Select Model",
            options=model_options,
            index=model_options.index(current_model) if current_model in model_options else 0,
            format_func=lambda x: f"{'☁️' if AVAILABLE_MODELS[x]['provider'] == 'anthropic' else '🖥️'} {AVAILABLE_MODELS[x]['name']}",
            help="☁️ Cloud (Claude) = Best quality | 🖥️ Local (Ollama) = On-device",
            label_visibility="collapsed"
        )

        # Show model info with styled badge
        model_info = AVAILABLE_MODELS[selected_model]
        badge_class = "cloud" if model_info['provider'] == 'anthropic' else "local"
        badge_text = "CLOUD" if model_info['provider'] == 'anthropic' else "LOCAL"
        st.markdown(f"""
        <div style="display: flex; align-items: center; gap: 0.5rem; margin-top: -0.5rem;">
            <span class="model-badge {badge_class}">{badge_text}</span>
            <span style="color: #8b949e; font-size: 0.8rem;">{model_info['description']} · {model_info['size']}</span>
        </div>
        """, unsafe_allow_html=True)

        # Update session state if model changed
        if selected_model != current_model:
            st.session_state.selected_model = selected_model
            st.toast(f"Switched to {model_info['name']}", icon="🔄")

        # ═══════════════════════════════════════════════════════════
        # SECTION 2: Demo Questions (Moved here - right after LLM)
        # ═══════════════════════════════════════════════════════════
        st.markdown('<div class="section-header">Demo Questions</div>', unsafe_allow_html=True)
        st.caption("Click to ask showcase questions")

        # FTD/VCP Demo Questions (Primary narrative)
        st.markdown('<p style="color: #2ea44f; font-size: 0.7rem; margin: 0.5rem 0 0.25rem 0;">FTD/VCP DEMO</p>', unsafe_allow_html=True)

        if st.button("🧠 FTD-Associated Variants", use_container_width=True, key="demo_ftd_variants"):
            st.session_state.quick_query = "What variants are associated with frontotemporal dementia?"

        if st.button("🔬 VCP Gene Analysis", use_container_width=True, key="demo_vcp"):
            st.session_state.quick_query = "Tell me about VCP variants and their disease associations. Is VCP a potential drug target?"

        if st.button("🧬 rs188935092 Details", use_container_width=True, key="demo_rsid"):
            st.session_state.quick_query = "What is the clinical significance of rs188935092? What gene is it in and what conditions is it associated with?"

        if st.button("🧪 Neurodegeneration Genes", use_container_width=True, key="demo_neuro"):
            st.session_state.quick_query = "What genes in this genome are linked to neurodegeneration? Include VCP, C9orf72, GRN, and any others found."

        # Other Demo Questions
        st.markdown('<p style="color: #6b7280; font-size: 0.7rem; margin: 0.5rem 0 0.25rem 0;">OTHER QUERIES</p>', unsafe_allow_html=True)

        if st.button("🔴 Pathogenic Variants", use_container_width=True, key="demo_pathogenic"):
            st.session_state.quick_query = "Are there any pathogenic variants? What diseases are they associated with?"

        if st.button("💊 Drug Targets", use_container_width=True, key="demo_drug"):
            st.session_state.quick_query = "Based on the variants found, what genes could be potential drug targets?"

        if st.button("💊 PPI/Acid Reflux", use_container_width=True, key="demo_ppi"):
            st.session_state.quick_query = "Would a proton pump inhibitor work well for this patient? Check CYP2C19 variants."

        # ═══════════════════════════════════════════════════════════
        # SECTION 3: Search Filters (Now after Demo Questions)
        # ═══════════════════════════════════════════════════════════
        st.markdown('<div class="section-header">Search Filters</div>', unsafe_allow_html=True)

        # Tabs for different sidebar sections
        sidebar_tab = st.radio(
            "Navigation",
            ["Filters", "Targets", "Files", "VCF Preview", "Metrics"],
            horizontal=True,
            label_visibility="collapsed"
        )

        if sidebar_tab == "Filters":
            render_search_filters()
        elif sidebar_tab == "Targets":
            render_targets_sidebar()
        elif sidebar_tab == "Files":
            render_file_manager()
        elif sidebar_tab == "VCF Preview":
            render_vcf_preview_sidebar()
        elif sidebar_tab == "Metrics":
            render_metrics_sidebar()

        # ═══════════════════════════════════════════════════════════
        # Footer with status
        # ═══════════════════════════════════════════════════════════
        st.markdown("---")
        current = st.session_state.get('selected_model', settings.LLM_MODEL)
        model_name = AVAILABLE_MODELS.get(current, {}).get('name', current)
        st.markdown(f"""
        <div style="text-align: center; padding: 0.5rem; background: #ffffff; border-radius: 8px; border: 1px solid #e1e4e8;">
            <span style="color: #2ea44f; font-size: 0.75rem;">●</span>
            <span style="color: #6b7280; font-size: 0.75rem;"> {model_name} Active</span>
        </div>
        """, unsafe_allow_html=True)

        return get_filter_expression()


def render_search_filters():
    """Render search filter controls."""
    st.subheader("Search Filters")

    st.session_state.filter_gene = st.text_input(
        "Filter by Gene",
        value=st.session_state.get('filter_gene', ''),
        placeholder="e.g., EGFR, BRCA1",
        help="Only search variants in this gene"
    )

    st.session_state.filter_chrom = st.selectbox(
        "Filter by Chromosome",
        ["All"] + settings.INCLUDE_CHROMOSOMES,
        help="Only search variants on this chromosome"
    )

    st.session_state.filter_impact = st.multiselect(
        "Filter by Impact",
        ["HIGH", "MODERATE", "LOW", "MODIFIER"],
        help="Only include variants with these impact levels"
    )

    # Database stats
    st.markdown("---")
    st.subheader("Database Info")
    try:
        engine = get_rag_engine()
        stats = engine.milvus.get_stats()
        st.metric("Total Variants", f"{stats['num_entities']:,}")
    except Exception as e:
        st.error(f"Database not connected: {e}")


def render_targets_sidebar():
    """Render target hypotheses management in sidebar."""
    st.subheader("🎯 Target Hypotheses")

    manager = get_target_manager()
    summary = manager.get_summary()

    # Summary metrics
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Total Targets", summary['total'])
    with col2:
        st.metric("High Priority", summary['high_priority'])

    st.markdown("---")

    # List existing targets
    targets = manager.list_all()
    if targets:
        st.caption("Saved Targets:")
        for target in targets[:5]:  # Show first 5
            priority_icon = "🔴" if target.priority >= 4 else "🟡" if target.priority >= 2 else "⚪"
            with st.expander(f"{priority_icon} {target.gene}", expanded=False):
                st.write(f"**Confidence:** {target.confidence}")
                st.write(f"**Priority:** {target.priority}/5")
                st.write(f"**Variants:** {target.variant_count}")
                st.write(f"**Status:** {target.status}")
                if target.rationale:
                    st.write(f"**Rationale:** {target.rationale[:100]}...")

                if st.button(f"Delete", key=f"del_{target.id}"):
                    manager.delete(target.id)
                    st.rerun()

        if len(targets) > 5:
            st.caption(f"... and {len(targets) - 5} more")
    else:
        st.info("No targets saved yet. Use the chat to identify potential drug targets.")

    st.markdown("---")

    # Export button
    if targets:
        if st.button("📤 Export for Phase 5"):
            output_file = manager.export_for_phase5()
            st.success(f"Exported to {output_file}")

            # Show download link
            with open(output_file, 'r') as f:
                st.download_button(
                    "Download JSON",
                    f.read(),
                    file_name="targets_for_phase5.json",
                    mime="application/json"
                )


def render_vcf_preview_sidebar():
    """Render VCF preview in sidebar."""
    st.subheader("VCF Preview")

    vcf_path = str(settings.VCF_INPUT_PATH)
    if Path(vcf_path).exists():
        st.success(f"VCF: {Path(vcf_path).name}")

        # Show variant count option
        limit = st.slider("Preview rows", 10, 100, 25)

        if st.button("Load Preview"):
            with st.spinner("Loading VCF..."):
                variants = get_vcf_preview(vcf_path, limit)
                st.session_state.vcf_preview = variants

        if 'vcf_preview' in st.session_state and st.session_state.vcf_preview:
            st.dataframe(
                st.session_state.vcf_preview,
                use_container_width=True,
                height=300
            )
    else:
        st.warning(f"VCF not found: {vcf_path}")


def render_file_manager():
    """Render file manager for browsing and uploading VCF files."""
    import os
    import datetime

    st.subheader("File Manager")

    # Define directories
    input_dir = settings.DATA_DIR / "input"
    output_dir = settings.DATA_DIR / "output"

    # Ensure directories exist
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Directory selector
    selected_dir = st.selectbox(
        "Browse Directory",
        ["input", "output"],
        format_func=lambda x: f"📁 {x.upper()}"
    )

    current_dir = input_dir if selected_dir == "input" else output_dir

    # ═══════════════════════════════════════════════════════════
    # FILE UPLOAD SECTION (only for input directory)
    # ═══════════════════════════════════════════════════════════
    if selected_dir == "input":
        st.markdown("---")
        st.markdown('<p style="color: #2ea44f; font-size: 0.75rem; font-weight: 600; text-transform: uppercase; letter-spacing: 1px;">Upload VCF Files</p>', unsafe_allow_html=True)

        uploaded_file = st.file_uploader(
            "Choose VCF file",
            type=["vcf", "vcf.gz", "gz"],
            help="Upload VCF or compressed VCF.gz files",
            label_visibility="collapsed"
        )

        if uploaded_file is not None:
            # Validate file extension
            valid_extensions = ['.vcf', '.vcf.gz', '.gz']
            file_ext = ''.join(Path(uploaded_file.name).suffixes).lower()

            if any(uploaded_file.name.lower().endswith(ext) for ext in valid_extensions):
                # Show file info
                file_size = len(uploaded_file.getvalue())
                st.info(f"📄 **{uploaded_file.name}** ({format_file_size(file_size)})")

                if st.button("⬆️ Upload File", use_container_width=True):
                    try:
                        # Save file
                        dest_path = input_dir / uploaded_file.name
                        with open(dest_path, 'wb') as f:
                            f.write(uploaded_file.getvalue())
                        st.success(f"Uploaded: {uploaded_file.name}")
                        st.rerun()
                    except Exception as e:
                        st.error(f"Upload failed: {e}")
            else:
                st.warning("Please upload a .vcf or .vcf.gz file")

    # ═══════════════════════════════════════════════════════════
    # FILE LISTING SECTION
    # ═══════════════════════════════════════════════════════════
    st.markdown("---")
    st.markdown(f'<p style="color: #2ea44f; font-size: 0.75rem; font-weight: 600; text-transform: uppercase; letter-spacing: 1px;">Files in {selected_dir.upper()}</p>', unsafe_allow_html=True)

    # Get files in directory
    try:
        files = []
        for f in current_dir.iterdir():
            if f.is_file():
                stat = f.stat()
                files.append({
                    'name': f.name,
                    'path': f,
                    'size': stat.st_size,
                    'modified': datetime.datetime.fromtimestamp(stat.st_mtime),
                    'type': get_file_type(f.name)
                })

        # Sort by modification time (newest first)
        files.sort(key=lambda x: x['modified'], reverse=True)

        if files:
            for file_info in files:
                # File type icon
                icon = "🧬" if file_info['type'] == 'vcf' else "📄"

                # Create expandable file entry
                with st.expander(f"{icon} {file_info['name']}", expanded=False):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.caption(f"**Size:** {format_file_size(file_info['size'])}")
                    with col2:
                        st.caption(f"**Modified:** {file_info['modified'].strftime('%Y-%m-%d %H:%M')}")

                    # Action buttons
                    btn_col1, btn_col2 = st.columns(2)

                    with btn_col1:
                        # Download button
                        with open(file_info['path'], 'rb') as f:
                            st.download_button(
                                "⬇️ Download",
                                f.read(),
                                file_name=file_info['name'],
                                mime="application/octet-stream",
                                use_container_width=True,
                                key=f"dl_{file_info['name']}"
                            )

                    with btn_col2:
                        # Delete button
                        if st.button("🗑️ Delete", use_container_width=True, key=f"del_{file_info['name']}"):
                            try:
                                file_info['path'].unlink()
                                st.success(f"Deleted: {file_info['name']}")
                                st.rerun()
                            except Exception as e:
                                st.error(f"Delete failed: {e}")

            # Summary stats
            st.markdown("---")
            total_size = sum(f['size'] for f in files)
            vcf_count = sum(1 for f in files if f['type'] == 'vcf')
            st.caption(f"**{len(files)} files** ({format_file_size(total_size)}) | **{vcf_count} VCF files**")
        else:
            st.info(f"No files in {selected_dir}/ directory")

    except Exception as e:
        st.error(f"Error reading directory: {e}")


def format_file_size(size_bytes: int) -> str:
    """Format file size in human-readable format."""
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024 * 1024:
        return f"{size_bytes / 1024:.1f} KB"
    elif size_bytes < 1024 * 1024 * 1024:
        return f"{size_bytes / (1024 * 1024):.1f} MB"
    else:
        return f"{size_bytes / (1024 * 1024 * 1024):.2f} GB"


def get_file_type(filename: str) -> str:
    """Determine file type from filename."""
    lower_name = filename.lower()
    if lower_name.endswith('.vcf') or lower_name.endswith('.vcf.gz'):
        return 'vcf'
    elif lower_name.endswith('.json'):
        return 'json'
    elif lower_name.endswith('.tsv') or lower_name.endswith('.csv'):
        return 'tabular'
    else:
        return 'other'


def render_metrics_sidebar():
    """Render LLM performance metrics."""
    st.subheader("LLM Performance")

    # Get metrics from session state
    metrics = st.session_state.get('llm_metrics', {})

    col1, col2 = st.columns(2)
    with col1:
        ttft = metrics.get('ttft', '--')
        st.metric(
            "TTFT",
            f"{ttft} ms" if isinstance(ttft, (int, float)) else ttft,
            help="Time to First Token"
        )
    with col2:
        tps = metrics.get('tokens_per_sec', '--')
        st.metric(
            "Tokens/sec",
            f"{tps:.1f}" if isinstance(tps, (int, float)) else tps,
            help="Token generation speed"
        )

    st.markdown("---")

    # Cache stats
    st.caption("Cache Statistics")
    cache_hits = metrics.get('cache_hits', 0)
    cache_misses = metrics.get('cache_misses', 0)
    total = cache_hits + cache_misses
    hit_rate = (cache_hits / total * 100) if total > 0 else 0

    st.progress(hit_rate / 100)
    st.caption(f"Hit Rate: {hit_rate:.1f}% ({cache_hits}/{total})")

    # Query stats
    st.markdown("---")
    st.caption("Session Stats")
    st.text(f"Queries: {metrics.get('total_queries', 0)}")
    st.text(f"Avg Latency: {metrics.get('avg_latency', '--')} ms")


def get_filter_expression():
    """Build filter expression from session state."""
    filter_parts = []

    filter_gene = st.session_state.get('filter_gene', '')
    filter_chrom = st.session_state.get('filter_chrom', 'All')
    filter_impact = st.session_state.get('filter_impact', [])

    if filter_gene:
        filter_parts.append(f'gene == "{filter_gene}"')
    if filter_chrom != "All":
        filter_parts.append(f'chrom == "{filter_chrom}"')
    if filter_impact:
        impact_expr = " or ".join([f'impact == "{i}"' for i in filter_impact])
        filter_parts.append(f"({impact_expr})")

    return " and ".join(filter_parts) if filter_parts else None


def render_evidence(evidence_list: list):
    """Render retrieved evidence in expandable sections."""
    if not evidence_list:
        st.info("No evidence retrieved for this query.")
        return

    with st.expander(f"📚 Retrieved Evidence ({len(evidence_list)} items)", expanded=False):
        for i, ev in enumerate(evidence_list, 1):
            score_color = "green" if ev['score'] > 0.7 else "orange" if ev['score'] > 0.5 else "red"

            col1, col2, col3 = st.columns([3, 1, 1])
            with col1:
                st.markdown(f"**{i}. {ev.get('gene', 'Unknown gene')}** - {ev['chrom']}:{ev['pos']}")
            with col2:
                st.markdown(f":{score_color}[Score: {ev['score']:.2f}]")
            with col3:
                # AlphaMissense badge
                am_score = ev.get('am_pathogenicity')
                if am_score is not None:
                    am_class = ev.get('am_class', '')
                    if am_score >= 0.564:
                        st.markdown(f":red[AM: {am_score:.2f}]")
                    elif am_score <= 0.34:
                        st.markdown(f":green[AM: {am_score:.2f}]")
                    else:
                        st.markdown(f":orange[AM: {am_score:.2f}]")

            st.markdown(f"_{ev['text_summary']}_")

            # Details in columns - now with 5 columns for AlphaMissense
            detail_cols = st.columns(5)
            with detail_cols[0]:
                st.caption(f"Ref: {ev['ref']}")
            with detail_cols[1]:
                st.caption(f"Alt: {ev['alt']}")
            with detail_cols[2]:
                st.caption(f"Consequence: {ev.get('consequence', 'N/A')}")
            with detail_cols[3]:
                st.caption(f"Impact: {ev.get('impact', 'N/A')}")
            with detail_cols[4]:
                am_score = ev.get('am_pathogenicity')
                if am_score is not None:
                    st.caption(f"AlphaMissense: {ev.get('am_class', 'N/A')}")
                else:
                    st.caption("AlphaMissense: N/A")

            st.markdown("---")


def render_knowledge_connections(evidence_list: list):
    """
    Render knowledge connection visualization (Clinker-style semantic layer).
    Shows: Variant → Gene → Protein → Pathway → Disease → Drug
    """
    if not evidence_list:
        return

    # Extract unique genes from evidence that have knowledge connections
    genes_with_connections = []
    seen_genes = set()

    for ev in evidence_list:
        gene = ev.get('gene', '').upper()
        if gene and gene in KNOWLEDGE_CONNECTIONS and gene not in seen_genes:
            genes_with_connections.append({
                'gene': gene,
                'variant': ev.get('rsid') or f"{ev.get('chrom')}:{ev.get('pos')}",
                'consequence': ev.get('consequence', 'variant'),
                'am_pathogenicity': ev.get('am_pathogenicity'),
                'am_class': ev.get('am_class'),
                'clinical_significance': ev.get('clinical_significance'),
                **KNOWLEDGE_CONNECTIONS[gene]
            })
            seen_genes.add(gene)

    if not genes_with_connections:
        return

    # Render the knowledge connection visualization
    with st.expander("🔗 Knowledge Connections (Clinker)", expanded=True):
        st.caption("Semantic connections from genomic variants to therapeutic opportunities")

        for conn in genes_with_connections[:3]:  # Limit to top 3 for clarity
            # Get the primary disease and drug
            primary_disease = conn['diseases'][0] if conn['diseases'] else 'Unknown'
            primary_drug = conn['drugs'][0] if conn['drugs'] else 'None'

            # Druggable badge
            if conn['druggable']:
                st.success(f"🧬 **{conn['gene']} Knowledge Path** — DRUGGABLE TARGET")
            else:
                st.info(f"🧬 **{conn['gene']} Knowledge Path** — Research Target")

            # Create the path visualization using columns
            cols = st.columns(6)

            with cols[0]:
                # Build variant box with AlphaMissense score if available
                am_score = conn.get('am_pathogenicity')
                am_line = ""
                if am_score is not None:
                    am_color = "#dc2626" if am_score >= 0.564 else "#16a34a" if am_score <= 0.34 else "#d97706"
                    am_line = f'<div style="font-size: 0.6rem; color: {am_color}; font-weight: 600;">AM: {am_score:.2f}</div>'

                st.markdown(f"""
                <div style="background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%); border: 2px solid #f59e0b; border-radius: 8px; padding: 8px; text-align: center;">
                    <div style="font-size: 0.65rem; color: #92400e; text-transform: uppercase;">Variant</div>
                    <div style="font-size: 0.8rem; font-weight: 600; color: #78350f;">{conn['variant'][:12]}</div>
                    {am_line}
                </div>
                """, unsafe_allow_html=True)

            with cols[1]:
                st.markdown(f"""
                <div style="background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%); border: 2px solid #3b82f6; border-radius: 8px; padding: 8px; text-align: center;">
                    <div style="font-size: 0.65rem; color: #1e40af; text-transform: uppercase;">Gene</div>
                    <div style="font-size: 0.8rem; font-weight: 600; color: #1e3a8a;">{conn['gene']}</div>
                </div>
                """, unsafe_allow_html=True)

            with cols[2]:
                st.markdown(f"""
                <div style="background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%); border: 2px solid #10b981; border-radius: 8px; padding: 8px; text-align: center;">
                    <div style="font-size: 0.65rem; color: #065f46; text-transform: uppercase;">Protein</div>
                    <div style="font-size: 0.8rem; font-weight: 600; color: #064e3b;">{conn['protein'][:12]}</div>
                </div>
                """, unsafe_allow_html=True)

            with cols[3]:
                st.markdown(f"""
                <div style="background: linear-gradient(135deg, #e0e7ff 0%, #c7d2fe 100%); border: 2px solid #6366f1; border-radius: 8px; padding: 8px; text-align: center;">
                    <div style="font-size: 0.65rem; color: #3730a3; text-transform: uppercase;">Pathway</div>
                    <div style="font-size: 0.8rem; font-weight: 600; color: #312e81;">{conn['pathway'][:12]}</div>
                </div>
                """, unsafe_allow_html=True)

            with cols[4]:
                st.markdown(f"""
                <div style="background: linear-gradient(135deg, #fee2e2 0%, #fecaca 100%); border: 2px solid #ef4444; border-radius: 8px; padding: 8px; text-align: center;">
                    <div style="font-size: 0.65rem; color: #991b1b; text-transform: uppercase;">Disease</div>
                    <div style="font-size: 0.8rem; font-weight: 600; color: #7f1d1d;">{primary_disease.split('(')[0].strip()[:12]}</div>
                </div>
                """, unsafe_allow_html=True)

            with cols[5]:
                st.markdown(f"""
                <div style="background: linear-gradient(135deg, #76b900 0%, #5a8c00 100%); border: 2px solid #4a7c00; border-radius: 8px; padding: 8px; text-align: center;">
                    <div style="font-size: 0.65rem; color: rgba(255,255,255,0.9); text-transform: uppercase;">Drug</div>
                    <div style="font-size: 0.8rem; font-weight: 600; color: white;">{primary_drug.split('(')[0].strip()[:12]}</div>
                </div>
                """, unsafe_allow_html=True)

            # Details section
            st.markdown("---")
            detail_cols = st.columns(2)
            with detail_cols[0]:
                st.markdown(f"**Function:** {conn['function']}")
                st.markdown(f"**Diseases:** {', '.join(conn['diseases'][:2])}")
            with detail_cols[1]:
                st.markdown(f"**Drugs:** {', '.join(conn['drugs'][:2])}")
                st.markdown(f"**Status:** {conn['drug_status']}")
                if conn['pdb_ids']:
                    st.markdown(f"**PDB:** {', '.join(conn['pdb_ids'])}")

            st.markdown("---")


def render_performance_panel():
    """Render real-time performance metrics panel."""
    metrics = st.session_state.get('llm_metrics', {})

    # Create metrics row
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        ttft = metrics.get('ttft', '--')
        st.metric(
            "Time to First Token",
            f"{ttft} ms" if isinstance(ttft, (int, float)) else ttft,
            delta=None
        )

    with col2:
        tps = metrics.get('tokens_per_sec', '--')
        st.metric(
            "Tokens/sec",
            f"{tps:.1f}" if isinstance(tps, (int, float)) else tps,
            delta=None
        )

    with col3:
        cache_hits = metrics.get('cache_hits', 0)
        cache_misses = metrics.get('cache_misses', 0)
        total = cache_hits + cache_misses
        hit_rate = (cache_hits / total * 100) if total > 0 else 0
        st.metric(
            "Cache Hit Rate",
            f"{hit_rate:.1f}%",
            delta=None
        )

    with col4:
        st.metric(
            "Total Queries",
            metrics.get('total_queries', 0),
            delta=None
        )


def render_target_panel():
    """Render target hypothesis management panel."""
    manager = get_target_manager()

    col1, col2 = st.columns([2, 1])

    with col1:
        st.subheader("🎯 Target Hypotheses")

    with col2:
        if st.button("➕ Add Target"):
            st.session_state.show_add_target = True

    # Add target form
    if st.session_state.get('show_add_target', False):
        with st.form("add_target_form"):
            st.markdown("### Add New Target Hypothesis")

            gene = st.text_input("Gene Symbol*", placeholder="e.g., EGFR, BRCA1")
            protein = st.text_input("Protein Name", placeholder="e.g., Epidermal growth factor receptor")
            rationale = st.text_area("Rationale*", placeholder="Why is this a potential drug target?")

            col1, col2 = st.columns(2)
            with col1:
                confidence = st.selectbox("Confidence", ["low", "medium", "high"])
            with col2:
                priority = st.slider("Priority", 1, 5, 3)

            therapeutic_area = st.text_input("Therapeutic Area", placeholder="e.g., Oncology, Cardiology")
            mechanism = st.text_input("Mechanism", placeholder="e.g., Kinase inhibition")
            pdb_ids = st.text_input("PDB IDs (comma-separated)", placeholder="e.g., 7SYE, 1M17")

            col1, col2 = st.columns(2)
            with col1:
                submitted = st.form_submit_button("Save Target")
            with col2:
                if st.form_submit_button("Cancel"):
                    st.session_state.show_add_target = False
                    st.rerun()

            if submitted and gene and rationale:
                # Get recent evidence from session state
                recent_evidence = []
                if st.session_state.messages:
                    for msg in st.session_state.messages[-3:]:
                        if 'evidence' in msg:
                            for ev in msg['evidence']:
                                if ev.get('gene', '').upper() == gene.upper():
                                    recent_evidence.append(ev)

                # Auto-populate from knowledge base if gene is known
                ref_data = get_gene_reference_data(gene)
                parsed_pdb_ids = [p.strip() for p in pdb_ids.split(',') if p.strip()] if pdb_ids else []

                target = TargetHypothesis(
                    gene=gene.upper(),
                    protein=protein or ref_data.get('protein') or None,
                    rationale=rationale,
                    confidence=confidence,
                    priority=priority,
                    therapeutic_area=therapeutic_area or None,
                    mechanism=mechanism or ref_data.get('mechanism') or None,
                    pdb_ids=parsed_pdb_ids or ref_data.get('pdb_ids', []),
                    uniprot_id=ref_data.get('uniprot_id'),
                    reference_smiles=ref_data.get('reference_smiles'),
                    reference_drug=ref_data.get('reference_drug'),
                    druggability='high' if ref_data.get('druggable') else 'low',
                    variants=recent_evidence,
                    source_query=st.session_state.messages[-1]['content'] if st.session_state.messages else None
                )

                manager.add(target)
                st.success(f"Target {gene} saved!")
                st.session_state.show_add_target = False
                st.rerun()

    # Display existing targets
    st.markdown("---")
    targets = manager.list_all()

    if targets:
        # Sort by priority
        targets_sorted = sorted(targets, key=lambda x: x.priority, reverse=True)

        for target in targets_sorted:
            priority_color = "red" if target.priority >= 4 else "orange" if target.priority >= 2 else "gray"

            with st.container():
                col1, col2, col3, col4 = st.columns([3, 1, 1, 1])

                with col1:
                    st.markdown(f"**{target.gene}** - {target.protein or 'Unknown protein'}")
                with col2:
                    st.markdown(f":{priority_color}[Priority: {target.priority}/5]")
                with col3:
                    st.markdown(f"Confidence: {target.confidence}")
                with col4:
                    st.markdown(f"Status: {target.status}")

                with st.expander("Details", expanded=False):
                    st.write(f"**Rationale:** {target.rationale}")
                    if target.therapeutic_area:
                        st.write(f"**Therapeutic Area:** {target.therapeutic_area}")
                    if target.mechanism:
                        st.write(f"**Mechanism:** {target.mechanism}")
                    if target.pdb_ids:
                        st.write(f"**PDB IDs:** {', '.join(target.pdb_ids)}")
                    st.write(f"**Variant Evidence:** {target.variant_count} variants")
                    st.write(f"**Created:** {target.created_at}")

                    # Status update
                    new_status = st.selectbox(
                        "Update Status",
                        ["hypothesis", "validated", "selected", "rejected"],
                        index=["hypothesis", "validated", "selected", "rejected"].index(target.status),
                        key=f"status_{target.id}"
                    )
                    if new_status != target.status:
                        manager.update(target.id, {'status': new_status})
                        st.rerun()

                    if st.button("🗑️ Delete", key=f"delete_{target.id}"):
                        manager.delete(target.id)
                        st.rerun()

                st.markdown("---")

        # Export section
        st.markdown("### Export Targets")
        col1, col2 = st.columns(2)
        with col1:
            if st.button("📤 Export for Phase 5 (Cryo-EM)"):
                output_file = manager.export_for_phase5()
                st.success(f"Exported {len([t for t in targets if t.status in ('validated', 'selected') or t.priority >= 4])} targets")

                with open(output_file, 'r') as f:
                    export_data = f.read()

                st.download_button(
                    "⬇️ Download JSON",
                    export_data,
                    file_name="targets_for_phase5.json",
                    mime="application/json"
                )

        with col2:
            st.caption("Exports validated/selected targets or those with priority >= 4")

    else:
        st.info("No target hypotheses saved yet. Use the chat to identify potential drug targets, then add them here.")


def update_metrics(ttft: float, total_tokens: int, duration: float):
    """Update LLM performance metrics and save to shared file."""
    if 'llm_metrics' not in st.session_state:
        st.session_state.llm_metrics = {
            'ttft': None,
            'tokens_per_sec': None,
            'cache_hits': 0,
            'cache_misses': 0,
            'total_queries': 0,
            'avg_latency': None,
            'latencies': [],
        }

    metrics = st.session_state.llm_metrics
    metrics['ttft'] = int(ttft * 1000)  # Convert to ms
    metrics['tokens_per_sec'] = total_tokens / duration if duration > 0 else 0
    metrics['total_queries'] += 1
    metrics['latencies'].append(duration * 1000)

    # Calculate average latency
    latencies = metrics['latencies'][-100:]  # Keep last 100
    metrics['avg_latency'] = int(sum(latencies) / len(latencies)) if latencies else 0

    # Save to shared metrics file for portal access
    save_shared_metrics(metrics)


def _extract_knowledge(evidence: list) -> list:
    """Extract unique knowledge connections from evidence genes."""
    seen = set()
    connections = []
    for ev in (evidence or []):
        gene = ev.get('gene', '').upper()
        if gene and gene not in seen and gene in KNOWLEDGE_CONNECTIONS:
            seen.add(gene)
            kc = KNOWLEDGE_CONNECTIONS[gene]
            connections.append({
                'gene': gene,
                'protein': kc.get('protein', ''),
                'pathway': kc.get('pathway', ''),
                'diseases': kc.get('diseases', []),
                'drugs': kc.get('drugs', []),
                'druggable': kc.get('druggable', False),
                'pdb_ids': kc.get('pdb_ids', []),
            })
    return connections


def trigger_report_refresh(query: str, evidence: list, response: str):
    """Trigger background PDF report regeneration with latest query context."""
    import subprocess
    import threading

    def _regenerate():
        try:
            # Write context for the report generator
            context = {
                'query': query,
                'response': response,
                'evidence': [
                    {
                        'gene': ev.get('gene', ''),
                        'rsid': ev.get('rsid', ''),
                        'chrom': ev.get('chrom', ''),
                        'pos': ev.get('pos', ''),
                        'consequence': ev.get('consequence', ''),
                        'impact': ev.get('impact', ''),
                        'clinical_significance': ev.get('clinical_significance', ''),
                        'score': ev.get('score', 0),
                        'am_pathogenicity': ev.get('am_pathogenicity', ''),
                        'am_class': ev.get('am_class', ''),
                        'text_summary': ev.get('text_summary', ''),
                    }
                    for ev in (evidence or [])
                ],
                'knowledge_connections': _extract_knowledge(evidence),
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                'total_queries': st.session_state.get('llm_metrics', {}).get('total_queries', 0),
                'model': st.session_state.get('selected_model', 'unknown'),
            }

            context_file = Path(__file__).parent.parent / 'data' / 'report_context.json'
            context_file.parent.mkdir(parents=True, exist_ok=True)
            with open(context_file, 'w') as f:
                json.dump(context, f, indent=2, default=str)

            # Call the report generator via subprocess (uses drug-discovery venv with reportlab)
            base_dir = Path(__file__).parent.parent.parent
            generator_script = base_dir / 'drug-discovery-pipeline' / 'generate_vcp_report_enhanced.py'
            generator_python = base_dir / 'drug-discovery-pipeline' / 'venv' / 'bin' / 'python'

            if generator_script.exists() and generator_python.exists():
                subprocess.run(
                    [str(generator_python), str(generator_script), '--context', str(context_file)],
                    cwd=str(generator_script.parent),
                    timeout=60,
                    capture_output=True,
                )
        except Exception:
            pass  # Report refresh is best-effort, never block the chat

    # Run in background thread so the chat UI is not blocked
    thread = threading.Thread(target=_regenerate, daemon=True)
    thread.start()


def save_shared_metrics(metrics: dict):
    """Save metrics to shared file for portal access."""
    metrics_file = settings.DATA_DIR / 'llm_metrics.json'
    try:
        # Create a serializable copy (exclude latencies list for file)
        shared_data = {
            'ttft': metrics.get('ttft'),
            'tokens_per_sec': metrics.get('tokens_per_sec'),
            'cache_hits': metrics.get('cache_hits', 0),
            'cache_misses': metrics.get('cache_misses', 0),
            'total_queries': metrics.get('total_queries', 0),
            'avg_latency': metrics.get('avg_latency'),
            'last_updated': time.strftime('%Y-%m-%d %H:%M:%S'),
        }
        with open(metrics_file, 'w') as f:
            json.dump(shared_data, f)
    except (OSError, IOError) as e:
        # Silently fail - metrics storage is optional
        pass


def main():
    """Main chat interface."""
    # Initialize session state
    if "messages" not in st.session_state:
        st.session_state.messages = []
    if "quick_query" not in st.session_state:
        st.session_state.quick_query = None
    if "selected_model" not in st.session_state:
        # Load from shared file (portal setting) or use default
        model, provider = load_shared_model()
        st.session_state.selected_model = model
    if "llm_metrics" not in st.session_state:
        st.session_state.llm_metrics = {
            'ttft': None,
            'tokens_per_sec': None,
            'cache_hits': 0,
            'cache_misses': 0,
            'total_queries': 0,
            'avg_latency': None,
            'latencies': [],
        }

    # Render sidebar and get filter
    filter_expr = render_sidebar()

    # Main content - Hero header
    st.markdown("""
    <div class="main-header">
        <h1>🧬 Genomic Evidence RAG Chat <span class="demo-badge">NVIDIA DGX Spark</span></h1>
        <p>AI-powered genomic variant analysis with ClinVar clinical annotations | HG002 (GIAB) WGS Sample</p>
    </div>
    """, unsafe_allow_html=True)

    # Stats row with better styling
    try:
        engine = get_rag_engine()
        stats = engine.milvus.get_stats()
        variant_count = stats.get('num_entities', 0)
    except Exception:
        variant_count = 35616  # Fallback to known count

    # Stats in styled columns - Light Theme
    st.markdown("""
    <style>
        .stats-container {
            display: flex;
            gap: 1rem;
            margin-bottom: 1.5rem;
        }
        .stat-box {
            flex: 1;
            background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);
            border: 1px solid #e1e4e8;
            border-radius: 10px;
            padding: 1rem;
            text-align: center;
            transition: all 0.2s ease;
        }
        .stat-box:hover {
            border-color: #76b900;
            box-shadow: 0 4px 12px rgba(118, 185, 0, 0.15);
            transform: translateY(-2px);
        }
        .stat-value {
            font-size: 1.5rem;
            font-weight: 700;
            color: #1f2937;
        }
        .stat-value.highlight {
            color: #2ea44f;
        }
        .stat-value.warning {
            color: #cf222e;
        }
        .stat-value.info {
            color: #0969da;
        }
        .stat-label {
            font-size: 0.75rem;
            color: #6b7280;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin-top: 0.25rem;
        }
        .stat-delta {
            font-size: 0.65rem;
            color: #cf222e;
            margin-top: 0.2rem;
        }
    </style>
    """, unsafe_allow_html=True)

    current_model = st.session_state.get('selected_model', settings.LLM_MODEL)
    model_name = AVAILABLE_MODELS.get(current_model, {}).get('name', current_model)

    st.markdown(f"""
    <div class="stats-container">
        <div class="stat-box">
            <div class="stat-value highlight">{variant_count:,}</div>
            <div class="stat-label">Annotated Variants</div>
        </div>
        <div class="stat-box">
            <div class="stat-value" style="color: #e65100;">13</div>
            <div class="stat-label">VCP Variants</div>
            <div class="stat-delta" style="color: #e65100;">FTD/ALS target</div>
        </div>
        <div class="stat-box">
            <div class="stat-value info">35</div>
            <div class="stat-label">FTD/ALS Genes</div>
            <div class="stat-delta">neurodegeneration</div>
        </div>
        <div class="stat-box">
            <div class="stat-value warning">3</div>
            <div class="stat-label">Pathogenic</div>
            <div class="stat-delta">disease-linked</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">55</div>
            <div class="stat-label">Drug Response</div>
        </div>
        <div class="stat-box">
            <div class="stat-value" style="font-size: 1rem;">{model_name}</div>
            <div class="stat-label">Active LLM</div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # Tabs for main content
    tab1, tab2, tab3 = st.tabs(["💬 Chat", "🎯 Target Hypotheses", "📊 Metrics"])

    with tab1:
        # Welcome message if no chat history
        if not st.session_state.messages:
            st.markdown("""
            <div class="welcome-card">
                <h3>Welcome to the Genomic Evidence RAG Chat</h3>
                <p style="color: #8b949e; margin-bottom: 1rem;">
                    AI-powered genomic variant analysis using <strong>Retrieval-Augmented Generation (RAG)</strong>
                    with ClinVar clinical annotations on HG002 (GIAB) whole genome sequencing data.
                </p>
            </div>
            """, unsafe_allow_html=True)

            # FTD Demo Highlight Box
            st.markdown("""
            <div style="background: linear-gradient(135deg, #fff3e0 0%, #ffffff 100%); border: 2px solid #ff9800; border-radius: 10px; padding: 1rem; margin-bottom: 1rem;">
                <h4 style="color: #e65100; margin: 0 0 0.5rem 0;">🧠 Featured Demo: Frontotemporal Dementia → Drug Discovery</h4>
                <p style="color: #4b5563; font-size: 0.85rem; margin: 0;">
                    Explore <strong>VCP gene variants</strong> linked to FTD/ALS, view Cryo-EM structures, and generate drug candidates.
                    Start with the <strong>FTD/VCP Demo</strong> buttons in the sidebar.
                </p>
            </div>
            """, unsafe_allow_html=True)

            col1, col2 = st.columns(2)

            with col1:
                st.markdown("""
                <div style="background: #ffffff; border: 1px solid #e1e4e8; border-radius: 10px; padding: 1rem; height: 100%;">
                    <h4 style="color: #2ea44f; margin-bottom: 0.75rem;">🧠 FTD Demo Questions</h4>
                    <p style="color: #1f2937; font-size: 0.9rem; margin-bottom: 0.5rem;">Click <strong>FTD/VCP Demo</strong> buttons or ask:</p>
                    <ul style="color: #4b5563; font-size: 0.85rem; margin: 0; padding-left: 1.2rem;">
                        <li>"What variants are associated with frontotemporal dementia?"</li>
                        <li>"Tell me about VCP variants and disease associations"</li>
                        <li>"What is the clinical significance of rs188935092?"</li>
                        <li>"What genes are linked to neurodegeneration?"</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)

            with col2:
                # Get dynamic variant stats from database
                variant_stats = get_variant_stats()
                vcp_count = variant_stats.get("vcp", 13)
                neuro_total = variant_stats.get("neuro_total", 35)

                st.markdown(f"""
                <div style="background: #ffffff; border: 1px solid #e1e4e8; border-radius: 10px; padding: 1rem; height: 100%;">
                    <h4 style="color: #2ea44f; margin-bottom: 0.75rem;">📊 FTD/Neurodegeneration Data</h4>
                    <ul style="color: #4b5563; font-size: 0.85rem; margin: 0; padding-left: 1.2rem;">
                        <li><span style="color: #e65100; font-weight: bold;">VCP:</span> {vcp_count} variants (rs188935092 = HIGH impact)</li>
                        <li><span style="color: #0969da;">FTD/ALS genes:</span> C9orf72, GRN, TBK1, FUS, TARDBP</li>
                        <li><span style="color: #2ea44f;">Total:</span> {neuro_total} neurodegeneration gene variants</li>
                        <li><span style="color: #6b7280;">Drug Target:</span> VCP inhibitors in clinical development</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)

            st.markdown("<br>", unsafe_allow_html=True)

        # Display chat history
        for message in st.session_state.messages:
            with st.chat_message(message["role"]):
                st.markdown(message["content"])
                if message["role"] == "assistant" and "evidence" in message:
                    render_knowledge_connections(message["evidence"])
                    render_evidence(message["evidence"])

        # Handle quick query from sidebar
        if st.session_state.quick_query:
            prompt = st.session_state.quick_query
            st.session_state.quick_query = None
        else:
            prompt = st.chat_input("Ask about genetic variants or drug targets...")

        if prompt:
            # Add user message
            st.session_state.messages.append({"role": "user", "content": prompt})
            with st.chat_message("user"):
                st.markdown(prompt)

            # Generate response
            with st.chat_message("assistant"):
                try:
                    engine = get_rag_engine()

                    with st.spinner("Searching evidence..."):
                        # Track timing
                        start_time = time.time()
                        first_token_time = None
                        token_count = 0

                        # Stream response
                        response_placeholder = st.empty()

                        full_response = ""
                        evidence = []

                        for chunk in engine.query_stream(prompt, filter_expr=filter_expr):
                            if chunk["type"] == "evidence":
                                evidence = chunk["content"]
                            elif chunk["type"] == "token":
                                if first_token_time is None:
                                    first_token_time = time.time()
                                token_count += 1
                                full_response += chunk["content"]
                                response_placeholder.markdown(full_response + "▌")
                            elif chunk["type"] == "done":
                                response_placeholder.markdown(full_response)

                        end_time = time.time()

                        # Update metrics
                        if first_token_time:
                            ttft = first_token_time - start_time
                            duration = end_time - start_time
                            update_metrics(ttft, token_count, duration)

                        # Show knowledge connections and evidence
                        render_knowledge_connections(evidence)
                        render_evidence(evidence)

                        # Save to history
                        st.session_state.messages.append({
                            "role": "assistant",
                            "content": full_response,
                            "evidence": evidence,
                        })

                        # Trigger pipeline report refresh (background, non-blocking)
                        trigger_report_refresh(prompt, evidence, full_response)

                except Exception as e:
                    st.error(f"Error: {e}")
                    st.info("Make sure Milvus is running and the database is populated.")

    with tab2:
        render_target_panel()

    with tab3:
        st.subheader("📊 LLM Performance Metrics")
        render_performance_panel()

        st.markdown("---")
        st.subheader("Session Information")
        col1, col2 = st.columns(2)
        current_model = st.session_state.get('selected_model', settings.LLM_MODEL)
        model_info = AVAILABLE_MODELS.get(current_model, {"name": current_model, "size": "?"})
        with col1:
            st.write(f"**Messages:** {len(st.session_state.messages)}")
            st.write(f"**Model:** {model_info['name']} ({model_info['size']})")
        with col2:
            st.write(f"**Provider:** {settings.LLM_PROVIDER}")
            st.write(f"**Embedding:** {settings.EMBEDDING_MODEL}")


if __name__ == "__main__":
    main()
