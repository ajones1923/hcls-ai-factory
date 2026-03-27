"""Pharmacogenomics Intelligence Agent — Streamlit UI v1.0.

Full-featured UI with:
- All 15 collection stats in sidebar
- 10 specialized PGx tabs (Dashboard, Drug Check, Medication Review,
  Warfarin Dosing, Chemo Safety, HLA Screening, Report Generator,
  Evidence Explorer, Phenoconversion Modeler, Population Analytics)
- Deep Research mode (autonomous agent pipeline)
- Conversation memory for follow-up queries
- Collection-specific filtering
- Citation relevance scoring (high/medium/low)
- Knowledge graph stats
- NVIDIA dark theme (green #76B900 on dark background)

Port: 8507 (assigned to Pharmacogenomics Intelligence Agent)

Usage:
    streamlit run app/pgx_ui.py --server.port 8507

Author: Adam Jones
Date: March 2026
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import streamlit as st

# Add project root to path (must happen before src imports)
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Load API key from rag-chat-pipeline .env if not already set
if not os.environ.get("ANTHROPIC_API_KEY"):
    from config.settings import settings
    if settings.ANTHROPIC_API_KEY:
        os.environ["ANTHROPIC_API_KEY"] = settings.ANTHROPIC_API_KEY
    else:
        env_path = settings.RAG_PIPELINE_ROOT / ".env"
        if env_path.exists():
            for line in env_path.read_text().splitlines():
                if line.startswith("ANTHROPIC_API_KEY="):
                    os.environ["ANTHROPIC_API_KEY"] = line.split("=", 1)[1].strip().strip('"')
                    break

from src.export import export_markdown, export_json, export_pdf, generate_filename


# ═══════════════════════════════════════════════════════════════════════
# PAGE CONFIG (must be the first Streamlit command)
# ═══════════════════════════════════════════════════════════════════════

st.set_page_config(
    page_title="PGx Intelligence | HCLS AI Factory",
    page_icon="\U0001f48a",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.warning(
    "**Clinical Decision Support Tool** — This system provides evidence-based guidance "
    "for research and clinical decision support only. All recommendations must be verified "
    "by a qualified healthcare professional. Not FDA-cleared. Not a substitute for professional "
    "clinical judgment."
)


# ═══════════════════════════════════════════════════════════════════════
# ENGINE INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════


@st.cache_resource
def init_engine():
    """Initialize the PGx RAG engine (cached across reruns)."""
    try:
        from src.collections import PGxCollectionManager
        from src.rag_engine import PGxRAGEngine
        from src import knowledge as kg
        from src import query_expansion as qe

        manager = PGxCollectionManager()
        manager.connect()

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
            from config.settings import settings as _settings

            class SimpleLLMClient:
                def __init__(self):
                    self.client = anthropic.Anthropic()

                def generate(self, prompt, system_prompt="", max_tokens=2048, temperature=0.7):
                    msg = self.client.messages.create(
                        model=_settings.LLM_MODEL,
                        max_tokens=max_tokens,
                        temperature=temperature,
                        system=system_prompt,
                        messages=[{"role": "user", "content": prompt}],
                    )
                    return msg.content[0].text

                def generate_stream(self, prompt, system_prompt="", max_tokens=2048, temperature=0.7):
                    with self.client.messages.stream(
                        model=_settings.LLM_MODEL,
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

        engine = PGxRAGEngine(
            collection_manager=manager,
            embedder=embedder,
            llm_client=llm_client,
            knowledge=kg,
            query_expander=qe,
        )
        return engine, manager
    except Exception as e:
        st.error(f"Failed to initialize: {e}")
        return None, None


@st.cache_resource
def init_agent(_engine):
    """Initialize the autonomous PGx Intelligence Agent."""
    if not _engine:
        return None
    try:
        from src.agent import PGxIntelligenceAgent
        return PGxIntelligenceAgent(_engine)
    except Exception:
        return None


@st.cache_resource
def init_pipeline_components():
    """Initialize PGx pipeline components (cached)."""
    components = {}
    try:
        from src.pgx_pipeline import StarAlleleCaller, PhenotypeTranslator, DrugGeneMatcher
        components["star_allele_caller"] = StarAlleleCaller()
        components["phenotype_translator"] = PhenotypeTranslator()
        components["drug_gene_matcher"] = DrugGeneMatcher()
    except Exception:
        pass
    try:
        from src.phenoconversion import PhenoconversionDetector
        components["phenoconversion_detector"] = PhenoconversionDetector()
    except Exception:
        pass
    try:
        from src.hla_screener import HLAScreener
        components["hla_screener"] = HLAScreener()
    except Exception:
        pass
    try:
        from src.dosing import DosingCalculator
        components["dosing_calculator"] = DosingCalculator()
    except Exception:
        pass
    return components


engine, manager = init_engine()
agent = init_agent(engine)
pipeline = init_pipeline_components()


# ═══════════════════════════════════════════════════════════════════════
# KNOWLEDGE DATA (lazy import for sidebar stats)
# ═══════════════════════════════════════════════════════════════════════


@st.cache_data
def load_knowledge():
    """Load PGx knowledge graph data for UI display."""
    try:
        from src.knowledge import (
            get_knowledge_stats, PHARMACOGENES, DRUG_CATEGORIES,
            METABOLIZER_PHENOTYPES, HLA_DRUG_ASSOCIATIONS,
            POPULATION_FREQUENCIES, DRUG_ALTERNATIVES,
            ACTIVITY_SCORE_TABLES, CLINICAL_DECISION_NOTES,
        )
        stats = get_knowledge_stats()
        all_drugs = []
        for cat_data in DRUG_CATEGORIES.values():
            all_drugs.extend(cat_data.get("drugs", []))
        all_drugs = sorted(set(all_drugs))
        return {
            "stats": stats,
            "pharmacogenes": PHARMACOGENES,
            "drug_categories": DRUG_CATEGORIES,
            "phenotypes": METABOLIZER_PHENOTYPES,
            "hla_associations": HLA_DRUG_ASSOCIATIONS,
            "population_frequencies": POPULATION_FREQUENCIES,
            "drug_alternatives": DRUG_ALTERNATIVES,
            "activity_scores": ACTIVITY_SCORE_TABLES,
            "clinical_notes": CLINICAL_DECISION_NOTES,
            "all_drugs": all_drugs,
        }
    except Exception:
        return None


knowledge_data = load_knowledge()


# ═══════════════════════════════════════════════════════════════════════
# CUSTOM CSS — NVIDIA Black + Green theme
# ═══════════════════════════════════════════════════════════════════════

st.markdown("""
<style>
    /* -- Base App -------------------------------------------------- */
    .stApp { background-color: #0a0a0f; }

    /* -- Streamlit Native Element Overrides ------------------------ */
    .stApp, .stApp p, .stApp span, .stApp li, .stApp td, .stApp th,
    .stApp label, .stApp .stMarkdown {
        color: #ffffff;
    }
    .stApp h1, .stApp h2, .stApp h3, .stApp h4, .stApp h5, .stApp h6 {
        color: #ffffff;
    }

    /* Sidebar */
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

    /* Text input and text area */
    .stTextInput input, .stTextArea textarea {
        background-color: #1a1a24 !important;
        color: #ffffff !important;
        border: 1px solid #333 !important;
    }
    .stTextInput input:focus, .stTextArea textarea:focus {
        border-color: #76B900 !important;
        box-shadow: 0 0 0 1px #76B900 !important;
    }

    /* Select boxes */
    .stSelectbox > div > div,
    .stMultiSelect > div > div {
        background-color: #1a1a24 !important;
        color: #ffffff !important;
        border-color: #333 !important;
    }

    /* Buttons */
    .stButton > button {
        background-color: #1a1a24;
        color: #e0e0e8;
        border: 1px solid #333;
        transition: all 0.2s ease;
    }
    .stButton > button:hover {
        background-color: #222230;
        border-color: #76B900;
        color: #ffffff;
    }

    /* Primary button */
    .stButton > button[kind="primary"] {
        background-color: #76B900;
        color: #0a0a0f;
        border: none;
        font-weight: 600;
    }
    .stButton > button[kind="primary"]:hover {
        background-color: #8DD100;
    }

    /* Expander */
    .streamlit-expanderHeader {
        background-color: #1a1a24 !important;
        color: #ffffff !important;
        border-radius: 8px;
    }
    .streamlit-expanderContent {
        background-color: #12121a !important;
        border: 1px solid #222230 !important;
    }

    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        background-color: transparent;
    }
    .stTabs [data-baseweb="tab"] {
        color: #a0a0b0;
    }
    .stTabs [aria-selected="true"] {
        color: #76B900 !important;
    }

    /* Metrics */
    [data-testid="stMetricValue"] {
        color: #76B900 !important;
    }
    [data-testid="stMetricLabel"] {
        color: #a0a0b0 !important;
    }

    /* Dividers */
    hr {
        border-color: #222230 !important;
    }

    /* Scrollbar */
    ::-webkit-scrollbar {
        width: 8px;
        height: 8px;
    }
    ::-webkit-scrollbar-track {
        background: #0a0a0f;
    }
    ::-webkit-scrollbar-thumb {
        background: #333;
        border-radius: 4px;
    }
    ::-webkit-scrollbar-thumb:hover {
        background: #555;
    }

    /* -- Custom Components ---------------------------------------- */
    .main-title {
        font-size: 2rem;
        font-weight: 700;
        color: #76B900;
        margin-bottom: 0.2rem;
    }
    .sub-title {
        font-size: 1rem;
        color: #a0a0b0;
        margin-bottom: 1.5rem;
    }

    .collection-badge {
        display: inline-block;
        padding: 2px 8px;
        border-radius: 4px;
        font-size: 0.75rem;
        font-weight: 600;
        margin-right: 4px;
    }
    .badge-genereference { background: #1D6FA4; color: white; }
    .badge-drugguideline { background: #76B900; color: white; }
    .badge-druginteraction { background: #DF6500; color: white; }
    .badge-hla { background: #E53935; color: white; }
    .badge-phenoconversion { background: #952FC6; color: white; }
    .badge-dosingalgorithm { background: #F9C500; color: #0a0a0f; }
    .badge-evidence { background: #1DBFA4; color: white; }
    .badge-populationdata { background: #00897B; color: white; }
    .badge-trial { background: #5C6BC0; color: white; }
    .badge-fdalabel { background: #F06292; color: white; }
    .badge-drugalternative { background: #8D6E63; color: white; }
    .badge-patientprofile { background: #26A69A; color: white; }
    .badge-implementation { background: #7E57C2; color: white; }
    .badge-education { background: #42A5F5; color: white; }
    .badge-genomic { background: #78909C; color: white; }

    .relevance-high { color: #76B900; font-weight: 700; }
    .relevance-medium { color: #F9C500; font-weight: 600; }
    .relevance-low { color: #6a6a7a; font-weight: 400; }

    .evidence-card {
        background: #1a1a24;
        border: 1px solid #222230;
        border-radius: 8px;
        padding: 12px;
        margin: 6px 0;
    }
    .evidence-card .score {
        color: #76B900;
        font-weight: 600;
        font-size: 0.85rem;
    }
    .evidence-card .snippet {
        color: #e0e0e8;
        font-size: 0.85rem;
        margin-top: 6px;
        line-height: 1.5;
    }
    .evidence-card a {
        color: #76B900;
        text-decoration: none;
    }
    .evidence-card a:hover {
        text-decoration: underline;
    }

    .alert-critical {
        background: #3a1010;
        border: 1px solid #E53935;
        border-radius: 8px;
        padding: 12px;
        margin: 6px 0;
        color: #ff8a80;
    }
    .alert-warning {
        background: #3a2a10;
        border: 1px solid #F9C500;
        border-radius: 8px;
        padding: 12px;
        margin: 6px 0;
        color: #ffe082;
    }
    .alert-info {
        background: #102a3a;
        border: 1px solid #42A5F5;
        border-radius: 8px;
        padding: 12px;
        margin: 6px 0;
        color: #90caf9;
    }

    .mode-badge {
        display: inline-block;
        padding: 3px 10px;
        border-radius: 12px;
        font-size: 0.7rem;
        font-weight: 700;
        margin-left: 8px;
    }
    .mode-deep { background: #952FC6; color: white; }
    .mode-quick { background: #1D6FA4; color: white; }

    .dose-result {
        background: #1a2a1a;
        border: 1px solid #76B900;
        border-radius: 8px;
        padding: 16px;
        margin: 8px 0;
    }
    .dose-value {
        font-size: 2rem;
        font-weight: 700;
        color: #76B900;
    }
    .dose-label {
        font-size: 0.9rem;
        color: #a0a0b0;
    }

    .phenoconv-arrow {
        font-size: 1.5rem;
        color: #F9C500;
        font-weight: 700;
    }
</style>
""", unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════════
# HEADER
# ═══════════════════════════════════════════════════════════════════════

st.markdown(
    '<div class="main-title">Pharmacogenomics Intelligence Agent</div>',
    unsafe_allow_html=True,
)
st.markdown(
    '<div class="sub-title">'
    'Precision Prescribing Intelligence Across 15 PGx Collections '
    '| HCLS AI Factory'
    '</div>',
    unsafe_allow_html=True,
)


# ═══════════════════════════════════════════════════════════════════════
# SESSION STATE INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════

if "pgx_messages" not in st.session_state:
    st.session_state["pgx_messages"] = []
if "pgx_conversation_context" not in st.session_state:
    st.session_state["pgx_conversation_context"] = ""
if "pgx_last_evidence" not in st.session_state:
    st.session_state["pgx_last_evidence"] = None
if "pgx_last_response" not in st.session_state:
    st.session_state["pgx_last_response"] = ""


# ═══════════════════════════════════════════════════════════════════════
# TABS
# ═══════════════════════════════════════════════════════════════════════

(tab_dashboard, tab_drug_check, tab_med_review, tab_warfarin,
 tab_chemo, tab_hla, tab_report, tab_evidence, tab_phenoconv,
 tab_population) = st.tabs([
    "PGx Dashboard",
    "Drug Check",
    "Medication Review",
    "Warfarin Dosing",
    "Chemotherapy Safety",
    "HLA Screening",
    "PGx Report Generator",
    "Evidence Explorer",
    "Phenoconversion Modeler",
    "Population Analytics",
])


# ═══════════════════════════════════════════════════════════════════════
# SIDEBAR
# ═══════════════════════════════════════════════════════════════════════

with st.sidebar:
    st.markdown("### Configuration")

    # Deep Research toggle
    deep_research = st.toggle(
        "Deep Research Mode", value=False,
        help="Uses autonomous agent with sub-question decomposition and evidence evaluation",
    )
    if deep_research:
        st.markdown('<span class="mode-badge mode-deep">DEEP RESEARCH</span>', unsafe_allow_html=True)
    else:
        st.markdown('<span class="mode-badge mode-quick">QUICK RAG</span>', unsafe_allow_html=True)

    # Collection selection
    st.markdown("---")
    st.markdown("### Collections")

    COLLECTION_LABELS = {
        "pgx_gene_reference": "Gene Reference",
        "pgx_drug_guidelines": "Drug Guidelines",
        "pgx_drug_interactions": "Drug Interactions",
        "pgx_hla_hypersensitivity": "HLA Hypersensitivity",
        "pgx_phenoconversion": "Phenoconversion",
        "pgx_dosing_algorithms": "Dosing Algorithms",
        "pgx_clinical_evidence": "Clinical Evidence",
        "pgx_population_data": "Population Data",
        "pgx_clinical_trials": "Clinical Trials",
        "pgx_fda_labels": "FDA Labels",
        "pgx_drug_alternatives": "Drug Alternatives",
        "pgx_patient_profiles": "Patient Profiles",
        "pgx_implementation": "Implementation",
        "pgx_education": "Education",
        "genomic_evidence": "Genomic Evidence",
    }

    # Get live stats
    collection_stats = {}
    if manager:
        try:
            collection_stats = manager.get_collection_stats()
        except Exception:
            pass

    # Collection toggles with record counts
    selected_collections = []
    for coll_id, label in COLLECTION_LABELS.items():
        count = collection_stats.get(coll_id, 0)
        checked = st.checkbox(f"{label} ({count:,})", value=True, key=f"coll_{coll_id}")
        if checked:
            selected_collections.append(coll_id)

    total_vectors = sum(collection_stats.values())
    st.markdown(f"**Total: {total_vectors:,} vectors across {len(selected_collections)} collections**")

    # Knowledge graph stats
    st.markdown("---")
    st.markdown("### Knowledge Graph")
    if knowledge_data and knowledge_data.get("stats"):
        ks = knowledge_data["stats"]
        st.markdown(f"- **Pharmacogenes:** {ks.get('pharmacogenes', 25)}")
        st.markdown(f"- **Phenotypes:** {ks.get('metabolizer_phenotypes', 5)}")
        st.markdown(f"- **Drug Categories:** {ks.get('drug_categories', 12)}")
        st.markdown(f"- **Drugs Tracked:** {ks.get('drugs_tracked', 100)}")
        st.markdown(f"- **HLA Associations:** {ks.get('hla_drug_associations', 12)}")
        st.markdown(f"- **Activity Score Tables:** {ks.get('activity_score_tables', 0)}")
        st.markdown(f"- **Entity Aliases:** {ks.get('entity_aliases', 0)}")
    else:
        st.markdown("- **Pharmacogenes:** 25")
        st.markdown("- **Phenotypes:** 5")
        st.markdown("- **Drug Categories:** 12")
        st.markdown("- **Drugs Tracked:** 100+")
        st.markdown("- **HLA Associations:** 12")

    # System Health
    st.markdown("---")
    st.markdown("### System Health")
    milvus_ok = manager is not None
    embedder_ok = engine is not None and engine.embedder is not None if engine else False
    llm_ok = engine is not None and engine.llm is not None if engine else False
    st.markdown(
        f"- Milvus: {'Connected' if milvus_ok else 'Disconnected'}\n"
        f"- Embedder: {'Ready' if embedder_ok else 'Unavailable'}\n"
        f"- LLM: {'Ready' if llm_ok else 'Unavailable'}"
    )

    # Demo Queries
    st.markdown("---")
    st.markdown("### Demo Queries")
    demo_queries = [
        "What is the CYP2D6 genotype impact on codeine efficacy?",
        "Warfarin dosing for CYP2C9 *1/*3, VKORC1 -1639 A/A",
        "DPYD screening before fluorouracil chemotherapy",
        "HLA-B*57:01 testing before abacavir prescription",
        "Compare CYP2D6 vs CYP2C19 clinical impact",
        "Phenoconversion risk with fluoxetine + codeine",
        "CYP3A5 expresser tacrolimus dosing in transplant",
        "TPMT and NUDT15 testing for thiopurine therapy",
        "Population frequency of CYP2D6*4 across ethnicities",
        "Clopidogrel alternatives for CYP2C19 poor metabolizer",
    ]
    for q in demo_queries:
        if st.button(q, key=f"demo_{q[:25]}"):
            st.session_state["demo_query"] = q

    # Clear conversation
    st.markdown("---")
    if st.button("Clear Conversation", key="clear_conv"):
        st.session_state["pgx_messages"] = []
        st.session_state["pgx_conversation_context"] = ""
        st.session_state["pgx_last_evidence"] = None
        st.session_state["pgx_last_response"] = ""
        st.rerun()


# ═══════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════


def render_evidence_cards(evidence):
    """Render evidence cards with relevance indicators."""
    if not evidence or not hasattr(evidence, "hits_by_collection"):
        return
    by_coll = evidence.hits_by_collection()
    for coll_name, hits in by_coll.items():
        badge_class = f"badge-{coll_name.lower()}"
        for hit in hits[:5]:
            source_link = ""
            if hit.collection == "Evidence" and hit.id.isdigit():
                source_link = (
                    f' <a href="https://pubmed.ncbi.nlm.nih.gov/{hit.id}/"'
                    f' target="_blank">PubMed</a>'
                )
            elif hit.collection == "Trial" and hit.id.upper().startswith("NCT"):
                source_link = (
                    f' <a href="https://clinicaltrials.gov/study/{hit.id}"'
                    f' target="_blank">ClinicalTrials.gov</a>'
                )

            relevance = hit.metadata.get("relevance", "")
            relevance_html = ""
            if relevance:
                relevance_html = f' <span class="relevance-{relevance}">[{relevance}]</span>'

            snippet = hit.text[:200].replace("<", "&lt;").replace(">", "&gt;")
            st.markdown(
                f'<div class="evidence-card">'
                f'<span class="collection-badge {badge_class}">{hit.collection}</span>'
                f' <strong>{hit.id}</strong>'
                f' <span class="score">{hit.score:.3f}</span>'
                f'{relevance_html}'
                f'{source_link}'
                f'<div class="snippet">{snippet}...</div>'
                f'</div>',
                unsafe_allow_html=True,
            )


def run_query(question, collections_filter=None):
    """Run a RAG query and return (response_text, evidence)."""
    if not engine:
        return "Engine not initialized. Check Milvus and embedder connections.", None

    conv_ctx = st.session_state.get("pgx_conversation_context", "")

    try:
        from src.models import AgentQuery
        agent_query = AgentQuery(question=question)
        evidence = engine.retrieve(
            agent_query,
            collections_filter=collections_filter or selected_collections,
            conversation_context=conv_ctx if conv_ctx else None,
        )

        if engine.llm:
            prompt = engine._build_prompt(question, evidence)
            from src.rag_engine import PGX_SYSTEM_PROMPT
            response = engine.llm.generate(
                prompt=prompt,
                system_prompt=PGX_SYSTEM_PROMPT,
                max_tokens=2048,
                temperature=0.7,
            )
        else:
            response = _format_evidence_only(evidence)

        # Update conversation memory
        st.session_state["pgx_conversation_context"] = (
            f"Previous Q: {question}\nPrevious A: {response[:500]}"
        )
        st.session_state["pgx_last_evidence"] = evidence
        st.session_state["pgx_last_response"] = response

        return response, evidence
    except Exception as e:
        return f"Query failed: {e}", None


def _format_evidence_only(evidence):
    """Format evidence results when LLM is unavailable."""
    if not evidence or not evidence.hits:
        return "No evidence found for this query."
    lines = [f"Found {len(evidence.hits)} results across {evidence.total_collections_searched} collections:\n"]
    for hit in evidence.hits[:10]:
        lines.append(f"- **[{hit.collection}:{hit.id}]** (score: {hit.score:.3f})")
        lines.append(f"  {hit.text[:200]}...")
    return "\n".join(lines)


def render_export_buttons(query, response, evidence=None, key_prefix="export"):
    """Render Markdown / JSON / PDF export buttons."""
    col_md, col_json, col_pdf = st.columns(3)
    with col_md:
        try:
            md_content = export_markdown(query, response, evidence)
            st.download_button(
                "Export Markdown",
                data=md_content,
                file_name=generate_filename("md"),
                mime="text/markdown",
                key=f"{key_prefix}_md",
            )
        except Exception:
            st.button("Export Markdown (unavailable)", disabled=True, key=f"{key_prefix}_md_dis")
    with col_json:
        try:
            json_content = export_json(query, response, evidence)
            st.download_button(
                "Export JSON",
                data=json_content,
                file_name=generate_filename("json"),
                mime="application/json",
                key=f"{key_prefix}_json",
            )
        except Exception:
            st.button("Export JSON (unavailable)", disabled=True, key=f"{key_prefix}_json_dis")
    with col_pdf:
        try:
            pdf_content = export_pdf(query, response, evidence)
            st.download_button(
                "Export PDF",
                data=pdf_content,
                file_name=generate_filename("pdf"),
                mime="application/pdf",
                key=f"{key_prefix}_pdf",
            )
        except Exception:
            st.button("Export PDF (unavailable)", disabled=True, key=f"{key_prefix}_pdf_dis")


def get_all_drug_names():
    """Return a sorted list of all drug names from knowledge."""
    if knowledge_data and knowledge_data.get("all_drugs"):
        return knowledge_data["all_drugs"]
    return [
        "abacavir", "allopurinol", "amitriptyline", "aripiprazole", "atomoxetine",
        "atorvastatin", "azathioprine", "capecitabine", "carbamazepine", "carvedilol",
        "celecoxib", "citalopram", "clopidogrel", "clozapine", "codeine",
        "cyclosporine", "dapsone", "desipramine", "duloxetine", "efavirenz",
        "escitalopram", "fentanyl", "flecainide", "fluconazole", "fluorouracil",
        "fluoxetine", "fluvastatin", "fluvoxamine", "haloperidol", "hydrocodone",
        "ibuprofen", "imipramine", "irinotecan", "isoniazid", "ivacaftor",
        "lamotrigine", "lansoprazole", "levetiracetam", "lovastatin", "mercaptopurine",
        "methadone", "metoprolol", "morphine", "nortriptyline", "olanzapine",
        "omeprazole", "ondansetron", "oxcarbazepine", "oxycodone", "pantoprazole",
        "paroxetine", "phenytoin", "pitavastatin", "pravastatin", "propafenone",
        "propranolol", "quetiapine", "rasburicase", "risperidone", "rosuvastatin",
        "sertraline", "simvastatin", "sirolimus", "tacrolimus", "tamoxifen",
        "tegafur", "thioguanine", "tramadol", "trimipramine", "valproic acid",
        "venlafaxine", "voriconazole", "warfarin",
    ]


def get_gene_list():
    """Return a list of pharmacogene names."""
    if knowledge_data and knowledge_data.get("pharmacogenes"):
        return sorted(knowledge_data["pharmacogenes"].keys())
    return [
        "ABCB1", "CACNA1S", "CFTR", "CYP1A2", "CYP2B6", "CYP2C19", "CYP2C8",
        "CYP2C9", "CYP2D6", "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "F5",
        "G6PD", "HLA", "IFNL3", "MTHFR", "NAT2", "NUDT15", "RYR1",
        "SLCO1B1", "TPMT", "UGT1A1", "VKORC1",
    ]


# ═══════════════════════════════════════════════════════════════════════
# TAB 1: PGx Dashboard
# ═══════════════════════════════════════════════════════════════════════

with tab_dashboard:
    st.markdown("## PGx Intelligence Dashboard")

    # Quick metrics row
    col1, col2, col3, col4, col5 = st.columns(5)
    with col1:
        st.metric("Collections", len(collection_stats) if collection_stats else 15)
    with col2:
        st.metric("Total Vectors", f"{total_vectors:,}" if total_vectors else "N/A")
    with col3:
        st.metric("Pharmacogenes", 25)
    with col4:
        drug_count = len(get_all_drug_names())
        st.metric("Drugs Tracked", drug_count)
    with col5:
        st.metric("HLA Pairs", 12)

    st.markdown("---")

    # Collection overview table
    st.markdown("### Collection Overview")
    if collection_stats:
        coll_data = []
        for coll_id, label in COLLECTION_LABELS.items():
            count = collection_stats.get(coll_id, 0)
            coll_data.append({"Collection": label, "Vectors": count})
        st.dataframe(coll_data, use_container_width=True, hide_index=True)
    else:
        st.info("Connect to Milvus to see live collection statistics.")

    st.markdown("---")

    # Quick-start query box
    st.markdown("### Quick Query")
    dash_query = st.text_input(
        "Ask a pharmacogenomics question",
        placeholder="e.g., What genes affect warfarin dosing?",
        key="dashboard_query",
    )
    # Check for demo query
    if st.session_state.get("demo_query"):
        dash_query = st.session_state.pop("demo_query")

    if dash_query:
        with st.spinner("Searching across PGx collections..."):
            response, evidence = run_query(dash_query)
        st.markdown("### Response")
        st.markdown(response)
        if evidence:
            with st.expander("Evidence Details", expanded=False):
                render_evidence_cards(evidence)
            render_export_buttons(dash_query, response, evidence, key_prefix="dash")

    # Pharmacogene overview
    st.markdown("---")
    st.markdown("### Key Pharmacogenes")
    gene_list = get_gene_list()
    gene_cols = st.columns(5)
    for i, gene in enumerate(gene_list):
        with gene_cols[i % 5]:
            if knowledge_data and knowledge_data.get("pharmacogenes", {}).get(gene):
                gdata = knowledge_data["pharmacogenes"][gene]
                cpic_count = len(gdata.get("cpic_guidelines", []))
                st.markdown(
                    f"**{gene}**  \n"
                    f"<span style='color:#a0a0b0;font-size:0.8rem;'>"
                    f"{gdata.get('full_name', '')[:30]}"
                    f" | {cpic_count} CPIC guideline{'s' if cpic_count != 1 else ''}"
                    f"</span>",
                    unsafe_allow_html=True,
                )
            else:
                st.markdown(f"**{gene}**")


# ═══════════════════════════════════════════════════════════════════════
# TAB 2: Drug Check
# ═══════════════════════════════════════════════════════════════════════

with tab_drug_check:
    st.markdown("## Single Drug PGx Lookup")
    st.markdown(
        "Select a drug and enter patient genotype/phenotype information to "
        "receive pharmacogenomic prescribing recommendations."
    )

    dc_col1, dc_col2 = st.columns(2)

    with dc_col1:
        all_drugs = get_all_drug_names()
        selected_drug = st.selectbox(
            "Select Drug",
            options=[""] + all_drugs,
            key="drug_check_drug",
            help="Choose from 100+ drugs with PGx annotations",
        )

    with dc_col2:
        gene_input = st.text_input(
            "Gene / Phenotype",
            placeholder="e.g., CYP2D6 poor metabolizer or CYP2C19 *1/*2",
            key="drug_check_gene",
        )

    if selected_drug and gene_input:
        query_text = (
            f"What are the pharmacogenomic prescribing recommendations for "
            f"{selected_drug} in a patient who is {gene_input}? "
            f"Include CPIC guidelines, dose adjustments, and alternative drugs."
        )

        if st.button("Get PGx Recommendation", key="drug_check_go", type="primary"):
            with st.spinner(f"Searching PGx evidence for {selected_drug}..."):
                response, evidence = run_query(query_text)

            st.markdown("### Recommendation")
            st.markdown(response)

            if evidence:
                with st.expander("Supporting Evidence", expanded=False):
                    render_evidence_cards(evidence)
                render_export_buttons(query_text, response, evidence, key_prefix="drugcheck")

    elif selected_drug:
        # Show drug category info
        if knowledge_data and knowledge_data.get("drug_categories"):
            for cat_name, cat_data in knowledge_data["drug_categories"].items():
                if selected_drug in cat_data.get("drugs", []):
                    st.info(
                        f"**Category:** {cat_name.replace('_', ' ').title()}  \n"
                        f"**Description:** {cat_data.get('description', '')}  \n"
                        f"**Primary Genes:** {', '.join(cat_data.get('primary_genes', []))}"
                    )
                    break

    # Drug categories overview
    st.markdown("---")
    st.markdown("### Drug Categories with PGx Relevance")
    if knowledge_data and knowledge_data.get("drug_categories"):
        for cat_name, cat_data in knowledge_data["drug_categories"].items():
            with st.expander(f"{cat_name.replace('_', ' ').title()} ({len(cat_data.get('drugs', []))} drugs)"):
                st.markdown(f"**{cat_data.get('description', '')}**")
                st.markdown(f"Primary genes: **{', '.join(cat_data.get('primary_genes', []))}**")
                st.markdown(f"Drugs: {', '.join(cat_data.get('drugs', []))}")


# ═══════════════════════════════════════════════════════════════════════
# TAB 3: Medication Review
# ═══════════════════════════════════════════════════════════════════════

with tab_med_review:
    st.markdown("## Polypharmacy Medication Review")
    st.markdown(
        "Enter a patient's medication list to check for drug-drug-gene "
        "interactions and phenoconversion risks. The phenoconversion "
        "detector identifies when concomitant medications alter effective "
        "metabolizer phenotypes."
    )

    med_list_input = st.text_area(
        "Medication List (one per line or comma-separated)",
        placeholder=(
            "fluoxetine\n"
            "codeine\n"
            "omeprazole\n"
            "metoprolol\n"
            "warfarin"
        ),
        height=150,
        key="med_review_list",
    )

    mr_col1, mr_col2 = st.columns(2)
    with mr_col1:
        mr_gene = st.selectbox(
            "Patient CYP2D6 Phenotype (genetic)",
            ["Normal Metabolizer", "Intermediate Metabolizer",
             "Poor Metabolizer", "Ultrarapid Metabolizer"],
            key="mr_cyp2d6_pheno",
        )
    with mr_col2:
        mr_gene2 = st.selectbox(
            "Patient CYP2C19 Phenotype (genetic)",
            ["Normal Metabolizer", "Intermediate Metabolizer",
             "Poor Metabolizer", "Ultrarapid Metabolizer", "Rapid Metabolizer"],
            key="mr_cyp2c19_pheno",
        )

    if st.button("Run Medication Review", key="med_review_go", type="primary"):
        if med_list_input.strip():
            # Parse medication list
            meds = [
                m.strip().lower()
                for m in med_list_input.replace(",", "\n").split("\n")
                if m.strip()
            ]

            st.markdown(f"### Reviewing {len(meds)} medications")
            st.markdown(f"Medications: **{', '.join(meds)}**")

            # Phenoconversion detection
            detector = pipeline.get("phenoconversion_detector")
            if detector:
                genetic_profile = {
                    "CYP2D6": {"phenotype": mr_gene, "diplotype": "unknown"},
                    "CYP2C19": {"phenotype": mr_gene2, "diplotype": "unknown"},
                    "CYP2C9": {"phenotype": "Normal Metabolizer", "diplotype": "unknown"},
                    "CYP3A4": {"phenotype": "Normal Metabolizer", "diplotype": "unknown"},
                    "CYP3A5": {"phenotype": "CYP3A5 Non-expresser", "diplotype": "unknown"},
                }

                alerts = detector.detect(genetic_profile, meds)

                if alerts:
                    st.markdown("### Phenoconversion Alerts")
                    for alert in alerts:
                        severity_class = "alert-critical" if "CLINICAL ACTION" in alert.clinical_significance else "alert-warning"
                        st.markdown(
                            f'<div class="{severity_class}">'
                            f'<strong>{alert.gene}</strong>: '
                            f'{alert.genetic_phenotype} &rarr; '
                            f'{alert.effective_phenotype}<br>'
                            f'<em>Precipitant:</em> {alert.precipitant_drug} '
                            f'({alert.inhibitor_strength} effect)<br>'
                            f'<em>Affected substrates:</em> '
                            f'{", ".join(alert.affected_substrates) if alert.affected_substrates else "None co-prescribed"}<br>'
                            f'<em>{alert.clinical_significance}</em>'
                            f'</div>',
                            unsafe_allow_html=True,
                        )
                else:
                    st.success("No phenoconversion events detected with the current medication list.")

                # Show effective phenotypes
                effective = detector.get_all_effective_phenotypes(genetic_profile, meds)
                st.markdown("### Effective Phenotypes After Phenoconversion")
                eff_data = []
                for gene, data in effective.items():
                    eff_data.append({
                        "Gene": gene,
                        "Genetic Phenotype": data.get("phenotype", ""),
                        "Effective Phenotype": data.get("effective_phenotype", ""),
                        "Phenoconverted": "Yes" if data.get("phenoconverted") else "No",
                    })
                st.dataframe(eff_data, use_container_width=True, hide_index=True)
            else:
                st.warning("Phenoconversion detector not available. Running RAG query instead.")

            # RAG query for comprehensive review
            st.markdown("---")
            st.markdown("### Comprehensive Interaction Analysis")
            review_query = (
                f"Review the following medication list for pharmacogenomic interactions, "
                f"drug-drug-gene interactions, and phenoconversion risks: {', '.join(meds)}. "
                f"Patient is {mr_gene} for CYP2D6 and {mr_gene2} for CYP2C19."
            )
            with st.spinner("Querying PGx evidence..."):
                response, evidence = run_query(review_query)
            st.markdown(response)
            if evidence:
                with st.expander("Evidence Details", expanded=False):
                    render_evidence_cards(evidence)
                render_export_buttons(review_query, response, evidence, key_prefix="medreview")
        else:
            st.warning("Please enter at least one medication.")


# ═══════════════════════════════════════════════════════════════════════
# TAB 4: Warfarin Dosing
# ═══════════════════════════════════════════════════════════════════════

with tab_warfarin:
    st.markdown("## IWPC Warfarin Dose Calculator")
    st.markdown(
        "Pharmacogenomic warfarin dosing using the International Warfarin "
        "Pharmacogenetics Consortium (IWPC) algorithm. Incorporates CYP2C9, "
        "VKORC1, and CYP4F2 genotypes with clinical variables."
    )

    wd_col1, wd_col2, wd_col3 = st.columns(3)

    with wd_col1:
        st.markdown("#### Genotype Inputs")
        cyp2c9_dip = st.selectbox(
            "CYP2C9 Diplotype",
            ["*1/*1", "*1/*2", "*1/*3", "*2/*2", "*2/*3", "*3/*3"],
            key="wd_cyp2c9",
        )
        vkorc1_geno = st.selectbox(
            "VKORC1 -1639 Genotype",
            ["G/G (wild-type)", "G/A (heterozygous)", "A/A (homozygous variant)"],
            key="wd_vkorc1",
        )
        cyp4f2_geno = st.selectbox(
            "CYP4F2 Genotype",
            ["*1/*1", "*1/*3", "*3/*3"],
            key="wd_cyp4f2",
        )

    with wd_col2:
        st.markdown("#### Clinical Variables")
        wd_age = st.number_input("Age (years)", min_value=18, max_value=100, value=60, key="wd_age")
        wd_weight = st.number_input("Weight (kg)", min_value=30.0, max_value=200.0, value=75.0, key="wd_weight")
        wd_height = st.number_input("Height (cm)", min_value=120.0, max_value=220.0, value=170.0, key="wd_height")

    with wd_col3:
        st.markdown("#### Additional Factors")
        wd_race = st.selectbox(
            "Race/Ethnicity",
            ["European", "African American", "Asian", "Other"],
            key="wd_race",
        )
        wd_amiodarone = st.checkbox("Amiodarone use", key="wd_amio")
        wd_smoker = st.checkbox("Current smoker", key="wd_smoker")
        wd_enzyme_inducer = st.checkbox(
            "Enzyme inducer use (rifampin, phenytoin, carbamazepine)",
            key="wd_inducer",
        )

    if st.button("Calculate Warfarin Dose", key="wd_calc", type="primary"):
        calculator = pipeline.get("dosing_calculator")
        if calculator:
            # Map VKORC1 selection to algorithm input
            vkorc1_map = {
                "G/G (wild-type)": "*1/*1",
                "G/A (heterozygous)": "*1/-1639G>A",
                "A/A (homozygous variant)": "-1639G>A/-1639G>A",
            }
            vkorc1_input = vkorc1_map.get(vkorc1_geno, "*1/*1")

            # Map race
            race_map = {
                "European": "european",
                "African American": "african_american",
                "Asian": "asian",
                "Other": "other",
            }
            race_input = race_map.get(wd_race, "other")

            result = calculator.warfarin_dose(
                cyp2c9_diplotype=cyp2c9_dip,
                vkorc1_genotype=vkorc1_input,
                cyp4f2_genotype=cyp4f2_geno,
                age=wd_age,
                weight=wd_weight,
                height=wd_height,
                race=race_input,
                amiodarone=wd_amiodarone,
                smoker=wd_smoker,
                enzyme_inducer=wd_enzyme_inducer,
            )

            # Display results
            st.markdown("### Predicted Warfarin Dose")
            res_c1, res_c2, res_c3 = st.columns(3)
            with res_c1:
                st.markdown(
                    f'<div class="dose-result">'
                    f'<div class="dose-value">{result["predicted_weekly_dose"]} mg</div>'
                    f'<div class="dose-label">Weekly Dose</div>'
                    f'</div>',
                    unsafe_allow_html=True,
                )
            with res_c2:
                st.markdown(
                    f'<div class="dose-result">'
                    f'<div class="dose-value">{result["predicted_daily_dose"]} mg</div>'
                    f'<div class="dose-label">Daily Dose</div>'
                    f'</div>',
                    unsafe_allow_html=True,
                )
            with res_c3:
                ci = result.get("confidence_interval", (0, 0))
                st.markdown(
                    f'<div class="dose-result">'
                    f'<div class="dose-value">{result.get("dose_category", "N/A")}</div>'
                    f'<div class="dose-label">95% CI: {ci[0]}-{ci[1]} mg/week</div>'
                    f'</div>',
                    unsafe_allow_html=True,
                )

            st.markdown(f"**Algorithm:** {result.get('algorithm_name', 'IWPC')}")

            # Clinical notes
            notes = result.get("clinical_notes", [])
            if notes:
                st.markdown("### Clinical Notes")
                for note in notes:
                    st.markdown(f"- {note}")

            # Variables used
            with st.expander("Variables Used in Calculation"):
                st.json(result.get("variables_used", {}))

            # Export
            dose_report = (
                f"## Warfarin Dose Prediction\n\n"
                f"**Weekly Dose:** {result['predicted_weekly_dose']} mg\n\n"
                f"**Daily Dose:** {result['predicted_daily_dose']} mg\n\n"
                f"**Category:** {result.get('dose_category', 'N/A')}\n\n"
                f"**Algorithm:** {result.get('algorithm_name', 'IWPC')}\n\n"
                f"### Clinical Notes\n\n"
                + "\n".join(f"- {n}" for n in notes)
            )
            render_export_buttons(
                f"Warfarin dosing: CYP2C9 {cyp2c9_dip}, VKORC1 {vkorc1_geno}",
                dose_report, key_prefix="warfarin",
            )
        else:
            st.warning("Dosing calculator not available. Running RAG query instead.")
            query_text = (
                f"Calculate warfarin dose for patient with CYP2C9 {cyp2c9_dip}, "
                f"VKORC1 {vkorc1_geno}, CYP4F2 {cyp4f2_geno}, age {wd_age}, "
                f"weight {wd_weight}kg, height {wd_height}cm, race {wd_race}, "
                f"amiodarone={'yes' if wd_amiodarone else 'no'}"
            )
            with st.spinner("Querying warfarin dosing evidence..."):
                response, evidence = run_query(query_text)
            st.markdown(response)


# ═══════════════════════════════════════════════════════════════════════
# TAB 5: Chemotherapy Safety
# ═══════════════════════════════════════════════════════════════════════

with tab_chemo:
    st.markdown("## Chemotherapy PGx Safety Screening")
    st.markdown(
        "Pre-treatment screening for DPYD (fluoropyrimidines), TPMT, and "
        "NUDT15 (thiopurines). Deficiency in these enzymes can cause "
        "life-threatening toxicity."
    )

    chemo_type = st.radio(
        "Select Chemotherapy Type",
        ["Fluoropyrimidines (5-FU, capecitabine, tegafur)",
         "Thiopurines (azathioprine, mercaptopurine, thioguanine)"],
        key="chemo_type",
    )

    calculator = pipeline.get("dosing_calculator")

    if "Fluoropyrimidine" in chemo_type:
        st.markdown("### DPYD Genotype-Guided Dosing")
        st.markdown(
            "Dihydropyrimidine dehydrogenase (DPD) is the rate-limiting enzyme "
            "for fluoropyrimidine catabolism. DPD deficiency causes severe, "
            "potentially fatal toxicity."
        )

        dpyd_dip = st.selectbox(
            "DPYD Diplotype",
            ["*1/*1 (Normal)", "*1/*2A (One no-function allele)",
             "*1/c.2846A>T (One decreased-function allele)",
             "*1/*13 (One no-function allele)",
             "*2A/*2A (Two no-function alleles)",
             "*1/HapB3 (One decreased-function allele)"],
            key="chemo_dpyd",
        )

        if st.button("Screen DPYD", key="chemo_dpyd_go", type="primary"):
            # Extract diplotype string
            dpyd_clean = dpyd_dip.split(" (")[0]

            if calculator:
                result = calculator.fluoropyrimidine_dose(dpyd_clean)

                # Color-code by severity
                if result["dose_reduction_percent"] >= 100:
                    alert_class = "alert-critical"
                elif result["dose_reduction_percent"] >= 50:
                    alert_class = "alert-warning"
                else:
                    alert_class = "alert-info"

                st.markdown(
                    f'<div class="{alert_class}">'
                    f'<strong>{result["phenotype"]}</strong><br>'
                    f'Activity Score: {result["activity_score"]}<br>'
                    f'Dose Reduction: {result["dose_reduction_percent"]}%<br><br>'
                    f'{result["recommendation"]}'
                    f'</div>',
                    unsafe_allow_html=True,
                )

                st.markdown("#### Affected Drugs")
                for drug in result.get("affected_drugs", []):
                    st.markdown(f"- {drug}")

                st.markdown("#### Monitoring")
                st.markdown(result.get("monitoring", ""))

                with st.expander("Allele Activity Scores"):
                    st.json(result.get("allele_scores", {}))
            else:
                with st.spinner("Running RAG query..."):
                    response, evidence = run_query(
                        f"DPYD genotype {dpyd_clean} fluoropyrimidine dose recommendation"
                    )
                st.markdown(response)

    else:
        st.markdown("### TPMT + NUDT15 Genotype-Guided Dosing")
        st.markdown(
            "Both TPMT and NUDT15 must be considered for thiopurine dosing. "
            "The more restrictive (lower activity) score determines the dose."
        )

        tc1, tc2 = st.columns(2)
        with tc1:
            tpmt_dip = st.selectbox(
                "TPMT Diplotype",
                ["*1/*1 (Normal)", "*1/*3A (One no-function)",
                 "*1/*3C (One no-function)", "*1/*2 (One no-function)",
                 "*3A/*3A (Two no-function)", "*3A/*3C (Two no-function)"],
                key="chemo_tpmt",
            )
        with tc2:
            nudt15_dip = st.selectbox(
                "NUDT15 Diplotype",
                ["*1/*1 (Normal)", "*1/*3 (One no-function)",
                 "*1/*2 (One no-function)", "*3/*3 (Two no-function)"],
                key="chemo_nudt15",
            )

        if st.button("Screen TPMT/NUDT15", key="chemo_thio_go", type="primary"):
            tpmt_clean = tpmt_dip.split(" (")[0]
            nudt15_clean = nudt15_dip.split(" (")[0]

            if calculator:
                result = calculator.thiopurine_dose(tpmt_clean, nudt15_clean)

                if result["dose_percent_of_standard"] <= 10:
                    alert_class = "alert-critical"
                elif result["dose_percent_of_standard"] <= 50:
                    alert_class = "alert-warning"
                else:
                    alert_class = "alert-info"

                st.markdown(
                    f'<div class="{alert_class}">'
                    f'<strong>TPMT:</strong> {result["tpmt_phenotype"]} '
                    f'(AS {result["tpmt_activity_score"]})<br>'
                    f'<strong>NUDT15:</strong> {result["nudt15_phenotype"]} '
                    f'(AS {result["nudt15_activity_score"]})<br>'
                    f'<strong>Dose:</strong> {result["dose_percent_of_standard"]}% of standard<br><br>'
                    f'{result["dose_recommendation"]}'
                    f'</div>',
                    unsafe_allow_html=True,
                )

                st.markdown("#### Affected Drugs")
                for drug in result.get("affected_drugs", []):
                    st.markdown(f"- {drug}")

                st.markdown("#### Monitoring")
                st.markdown(result.get("monitoring", ""))
            else:
                with st.spinner("Running RAG query..."):
                    response, evidence = run_query(
                        f"TPMT {tpmt_clean} NUDT15 {nudt15_clean} thiopurine dose recommendation"
                    )
                st.markdown(response)


# ═══════════════════════════════════════════════════════════════════════
# TAB 6: HLA Screening
# ═══════════════════════════════════════════════════════════════════════

with tab_hla:
    st.markdown("## HLA-Drug Hypersensitivity Screening")
    st.markdown(
        "Screen for HLA alleles associated with severe adverse drug reactions "
        "(SJS/TEN, DRESS, DILI, hypersensitivity syndrome). Pre-emptive "
        "testing prevents life-threatening reactions."
    )

    hla_mode = st.radio(
        "Screening Mode",
        ["Pre-Prescription (check one drug)", "Panel Screen (check all drugs)"],
        key="hla_mode",
    )

    screener = pipeline.get("hla_screener")

    if "Pre-Prescription" in hla_mode:
        hla_c1, hla_c2 = st.columns(2)
        with hla_c1:
            hla_drug = st.selectbox(
                "Drug to Screen",
                ["", "abacavir", "allopurinol", "carbamazepine", "clozapine",
                 "dapsone", "flucloxacillin", "lamotrigine", "methazolamide",
                 "minocycline", "nevirapine", "oxcarbazepine", "phenytoin",
                 "sulfasalazine", "ticlopidine", "trimethoprim-sulfamethoxazole"],
                key="hla_drug_select",
            )
        with hla_c2:
            st.markdown("#### Patient HLA Typing")
            hla_b_alleles = st.text_input(
                "HLA-B alleles (comma-separated)",
                placeholder="e.g., *57:01, *44:02",
                key="hla_b_input",
            )
            hla_a_alleles = st.text_input(
                "HLA-A alleles (comma-separated)",
                placeholder="e.g., *31:01, *02:01",
                key="hla_a_input",
            )
            hla_dqb1_alleles = st.text_input(
                "HLA-DQB1 alleles (comma-separated, for clozapine)",
                placeholder="e.g., *05:02, *03:01",
                key="hla_dqb1_input",
            )

        if hla_drug and st.button("Screen HLA", key="hla_screen_go", type="primary"):
            # Parse HLA typing
            hla_typing = {}
            if hla_b_alleles.strip():
                hla_typing["HLA-B"] = [a.strip() for a in hla_b_alleles.split(",") if a.strip()]
            if hla_a_alleles.strip():
                hla_typing["HLA-A"] = [a.strip() for a in hla_a_alleles.split(",") if a.strip()]
            if hla_dqb1_alleles.strip():
                hla_typing["HLA-DQB1"] = [a.strip() for a in hla_dqb1_alleles.split(",") if a.strip()]

            if not hla_typing:
                st.warning("Please enter at least one HLA allele.")
            elif screener:
                results = screener.screen_drug(hla_drug, hla_typing)

                st.markdown(f"### Screening Results for {hla_drug.title()}")
                for r in results:
                    if r.status.value == "CONTRAINDICATED":
                        alert_class = "alert-critical"
                    elif r.status.value == "HIGH_RISK":
                        alert_class = "alert-warning"
                    elif r.status.value == "SAFE":
                        alert_class = "alert-info"
                    else:
                        alert_class = "alert-info"

                    st.markdown(
                        f'<div class="{alert_class}">'
                        f'<strong>Status: {r.status.value}</strong><br>'
                        f'<strong>HLA Allele:</strong> {r.hla_allele}<br>'
                        f'<strong>Reaction:</strong> {r.reaction_type}<br>'
                        f'<strong>Severity:</strong> {r.severity.value}<br><br>'
                        f'{r.recommendation}<br><br>'
                        f'<em>Evidence:</em> {r.evidence_level}<br>'
                        f'<em>Population Risk:</em> {r.population_risk}'
                        f'</div>',
                        unsafe_allow_html=True,
                    )

                    if r.alternatives:
                        st.markdown("**Alternative Drugs:**")
                        for alt in r.alternatives:
                            st.markdown(f"- {alt}")
            else:
                with st.spinner("Running RAG query..."):
                    response, evidence = run_query(
                        f"HLA screening for {hla_drug} with patient HLA typing {hla_typing}"
                    )
                st.markdown(response)

    else:
        # Panel screen mode
        st.markdown("#### Enter Complete HLA Typing")
        panel_hla_b = st.text_input(
            "HLA-B alleles", placeholder="*57:01, *44:02", key="panel_hla_b",
        )
        panel_hla_a = st.text_input(
            "HLA-A alleles", placeholder="*31:01, *02:01", key="panel_hla_a",
        )
        panel_hla_drb1 = st.text_input(
            "HLA-DRB1 alleles (optional)", placeholder="*01:01, *07:01", key="panel_hla_drb1",
        )
        panel_hla_dqb1 = st.text_input(
            "HLA-DQB1 alleles (optional, for clozapine)", placeholder="*05:02, *03:01", key="panel_hla_dqb1",
        )

        if st.button("Run Panel Screen", key="hla_panel_go", type="primary"):
            hla_typing = {}
            if panel_hla_b.strip():
                hla_typing["HLA-B"] = [a.strip() for a in panel_hla_b.split(",") if a.strip()]
            if panel_hla_a.strip():
                hla_typing["HLA-A"] = [a.strip() for a in panel_hla_a.split(",") if a.strip()]
            if panel_hla_drb1.strip():
                hla_typing["HLA-DRB1"] = [a.strip() for a in panel_hla_drb1.split(",") if a.strip()]
            if panel_hla_dqb1.strip():
                hla_typing["HLA-DQB1"] = [a.strip() for a in panel_hla_dqb1.split(",") if a.strip()]

            if not hla_typing:
                st.warning("Please enter at least one HLA allele.")
            elif screener:
                results = screener.screen_all_drugs(hla_typing)

                if results:
                    st.markdown(f"### Panel Results: {len(results)} risk(s) identified")
                    for r in results:
                        alert_class = "alert-critical" if r.status.value == "CONTRAINDICATED" else "alert-warning"
                        st.markdown(
                            f'<div class="{alert_class}">'
                            f'<strong>{r.drug.title()}</strong> — {r.status.value}<br>'
                            f'HLA: {r.hla_allele} | Reaction: {r.reaction_type}<br>'
                            f'{r.recommendation}'
                            f'</div>',
                            unsafe_allow_html=True,
                        )
                else:
                    st.success(
                        "No HLA-drug contraindications or high-risk associations "
                        "identified for the provided HLA typing."
                    )
            else:
                st.warning("HLA screener not available.")

    # HLA association reference table
    st.markdown("---")
    st.markdown("### HLA-Drug Association Reference")
    if knowledge_data and knowledge_data.get("hla_associations"):
        hla_ref = []
        for key, data in knowledge_data["hla_associations"].items():
            hla_ref.append({
                "HLA Allele": data.get("hla_allele", ""),
                "Drug": data.get("drug", key.split("_")[-1] if "_" in key else key),
                "Reaction": data.get("reaction_type", data.get("reaction", "")),
                "Severity": data.get("severity", ""),
                "Evidence": data.get("evidence_level", ""),
            })
        if hla_ref:
            st.dataframe(hla_ref, use_container_width=True, hide_index=True)


# ═══════════════════════════════════════════════════════════════════════
# TAB 7: PGx Report Generator
# ═══════════════════════════════════════════════════════════════════════

with tab_report:
    st.markdown("## Comprehensive PGx Report Generator")
    st.markdown(
        "Generate a complete pharmacogenomic report for a patient. Enter "
        "genotyping results and medication list to produce an actionable "
        "clinical report with dosing recommendations, interaction alerts, "
        "and HLA screening."
    )

    rpt_col1, rpt_col2 = st.columns(2)

    with rpt_col1:
        st.markdown("#### Patient Genotype Panel")
        rpt_cyp2d6 = st.selectbox(
            "CYP2D6 Diplotype", ["*1/*1", "*1/*4", "*4/*4", "*1/*2", "*1/*41",
                                  "*2/*2", "*41/*41", "*1xN/*1xN"],
            key="rpt_cyp2d6",
        )
        rpt_cyp2c19 = st.selectbox(
            "CYP2C19 Diplotype", ["*1/*1", "*1/*2", "*2/*2", "*1/*17",
                                   "*17/*17", "*2/*17", "*1/*3"],
            key="rpt_cyp2c19",
        )
        rpt_cyp2c9 = st.selectbox(
            "CYP2C9 Diplotype", ["*1/*1", "*1/*2", "*1/*3", "*2/*2",
                                  "*2/*3", "*3/*3"],
            key="rpt_cyp2c9",
        )
        rpt_cyp3a5 = st.selectbox(
            "CYP3A5 Diplotype", ["*1/*1", "*1/*3", "*3/*3"],
            key="rpt_cyp3a5",
        )
        rpt_dpyd = st.selectbox(
            "DPYD Diplotype", ["*1/*1", "*1/*2A", "*1/c.2846A>T", "*2A/*2A"],
            key="rpt_dpyd",
        )
        rpt_tpmt = st.selectbox(
            "TPMT Diplotype", ["*1/*1", "*1/*3A", "*3A/*3A"],
            key="rpt_tpmt",
        )
        rpt_slco1b1 = st.selectbox(
            "SLCO1B1 Genotype", ["*1/*1 (521TT)", "*1/*5 (521TC)", "*5/*5 (521CC)"],
            key="rpt_slco1b1",
        )
        rpt_vkorc1 = st.selectbox(
            "VKORC1 -1639 Genotype", ["G/G", "G/A", "A/A"],
            key="rpt_vkorc1",
        )

    with rpt_col2:
        st.markdown("#### Patient Information")
        rpt_age = st.number_input("Age", min_value=1, max_value=120, value=55, key="rpt_age")
        rpt_sex = st.selectbox("Sex", ["Male", "Female"], key="rpt_sex")
        rpt_ethnicity = st.selectbox(
            "Ethnicity",
            ["European", "African American", "East Asian", "South Asian",
             "Latino/Hispanic", "Other"],
            key="rpt_ethnicity",
        )

        st.markdown("#### Current Medications")
        rpt_meds = st.text_area(
            "Medications (one per line)",
            placeholder="metoprolol\nomeprazole\nclopidogrel\nsimvastatin",
            height=150,
            key="rpt_meds",
        )

    if st.button("Generate PGx Report", key="rpt_generate", type="primary"):
        meds = [m.strip() for m in rpt_meds.split("\n") if m.strip()] if rpt_meds else []

        report_query = (
            f"Generate a comprehensive pharmacogenomic report for a {rpt_age}-year-old "
            f"{rpt_sex.lower()} ({rpt_ethnicity}) patient with the following genotypes: "
            f"CYP2D6 {rpt_cyp2d6}, CYP2C19 {rpt_cyp2c19}, CYP2C9 {rpt_cyp2c9}, "
            f"CYP3A5 {rpt_cyp3a5}, DPYD {rpt_dpyd}, TPMT {rpt_tpmt}, "
            f"SLCO1B1 {rpt_slco1b1.split(' ')[0]}, VKORC1 {rpt_vkorc1}. "
            f"Current medications: {', '.join(meds) if meds else 'none listed'}. "
            f"Include: phenotype assignments, drug recommendations, interaction alerts, "
            f"HLA screening needs, and clinical decision support."
        )

        with st.spinner("Generating comprehensive PGx report..."):
            response, evidence = run_query(report_query)

        st.markdown("### Patient PGx Report")
        st.markdown(response)

        if evidence:
            with st.expander("Supporting Evidence", expanded=False):
                render_evidence_cards(evidence)

        # Phenoconversion check if medications provided
        if meds and pipeline.get("phenoconversion_detector"):
            detector = pipeline["phenoconversion_detector"]
            translator = pipeline.get("phenotype_translator")

            genetic_profile = {
                "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": rpt_cyp2d6},
                "CYP2C19": {"phenotype": "Normal Metabolizer", "diplotype": rpt_cyp2c19},
                "CYP2C9": {"phenotype": "Normal Metabolizer", "diplotype": rpt_cyp2c9},
                "CYP3A5": {"phenotype": "CYP3A5 Non-expresser", "diplotype": rpt_cyp3a5},
            }

            # Translate diplotypes if translator available
            if translator:
                try:
                    for gene in ["CYP2D6", "CYP2C19", "CYP2C9"]:
                        dip = genetic_profile[gene]["diplotype"]
                        result = translator.translate(gene, dip)
                        if result:
                            genetic_profile[gene]["phenotype"] = result.get("phenotype", genetic_profile[gene]["phenotype"])
                except Exception:
                    pass

            alerts = detector.detect(genetic_profile, meds)
            if alerts:
                st.markdown("### Phenoconversion Alerts")
                for alert in alerts:
                    st.markdown(
                        f'<div class="alert-warning">'
                        f'<strong>{alert.gene}:</strong> '
                        f'{alert.genetic_phenotype} &rarr; {alert.effective_phenotype} '
                        f'(precipitant: {alert.precipitant_drug})'
                        f'</div>',
                        unsafe_allow_html=True,
                    )

        # Export buttons
        st.markdown("---")
        render_export_buttons(report_query, response, evidence, key_prefix="report")


# ═══════════════════════════════════════════════════════════════════════
# TAB 8: Evidence Explorer
# ═══════════════════════════════════════════════════════════════════════

with tab_evidence:
    st.markdown("## Evidence Explorer")
    st.markdown(
        "Search across all 15 PGx collections with optional Deep Research "
        "mode for autonomous multi-step evidence synthesis."
    )

    # Conversation history
    for msg in st.session_state["pgx_messages"]:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])

    # Chat input
    ev_query = st.chat_input("Ask a pharmacogenomics question...", key="evidence_chat")

    if ev_query:
        # Add user message
        st.session_state["pgx_messages"].append({"role": "user", "content": ev_query})
        with st.chat_message("user"):
            st.markdown(ev_query)

        # Run query
        with st.chat_message("assistant"):
            if deep_research and agent:
                with st.spinner("Deep Research in progress (multi-step evidence synthesis)..."):
                    try:
                        result = agent.run(ev_query)
                        response = result.get("answer", str(result))
                        evidence = result.get("evidence", None)
                    except Exception as e:
                        response = f"Deep Research failed: {e}. Falling back to standard RAG."
                        evidence = None
                        with st.spinner("Running standard RAG..."):
                            response, evidence = run_query(ev_query)
            else:
                with st.spinner("Searching across PGx collections..."):
                    response, evidence = run_query(ev_query)

            st.markdown(response)

            if evidence:
                with st.expander("Evidence Details", expanded=False):
                    render_evidence_cards(evidence)
                    st.markdown(
                        f"*{evidence.hit_count} results from "
                        f"{evidence.total_collections_searched} collections "
                        f"in {evidence.search_time_ms:.0f}ms*"
                    )

            # Export
            render_export_buttons(ev_query, response, evidence, key_prefix=f"ev_{len(st.session_state['pgx_messages'])}")

        # Add assistant message
        st.session_state["pgx_messages"].append({"role": "assistant", "content": response})

    # Entity explorer
    st.markdown("---")
    st.markdown("### Entity Explorer")
    st.markdown("Find all evidence related to a specific gene, drug, or concept across all 15 collections.")

    entity_input = st.text_input(
        "Enter entity (gene, drug, allele, etc.)",
        placeholder="e.g., CYP2D6, codeine, HLA-B*57:01",
        key="entity_explorer",
    )

    if entity_input and st.button("Find Related Evidence", key="entity_go"):
        if engine:
            with st.spinner(f"Searching all collections for '{entity_input}'..."):
                try:
                    related = engine.find_related(entity_input, top_k=5)
                    if related:
                        for coll_name, hits in related.items():
                            label = COLLECTION_LABELS.get(coll_name, coll_name)
                            with st.expander(f"{label} ({len(hits)} results)"):
                                for hit in hits:
                                    st.markdown(
                                        f"**{hit.id}** (score: {hit.score:.3f})  \n"
                                        f"{hit.text[:300]}..."
                                    )
                    else:
                        st.info(f"No evidence found for '{entity_input}'.")
                except Exception as e:
                    st.error(f"Entity search failed: {e}")
        else:
            st.warning("Engine not initialized.")


# ═══════════════════════════════════════════════════════════════════════
# TAB 9: Phenoconversion Modeler
# ═══════════════════════════════════════════════════════════════════════

with tab_phenoconv:
    st.markdown("## Interactive Phenoconversion Modeler")
    st.markdown(
        "Visualize how concomitant medications alter the effective "
        "metabolizer phenotype through enzyme inhibition or induction. "
        "Enter a patient's genetic phenotype and concurrent medications "
        "to see the effective phenotype."
    )

    pc_col1, pc_col2 = st.columns(2)

    with pc_col1:
        st.markdown("#### Genetic Phenotype")
        pc_gene = st.selectbox(
            "Select Gene",
            ["CYP2D6", "CYP2C19", "CYP2C9", "CYP3A4", "CYP3A5"],
            key="pc_gene",
        )
        pc_phenotype = st.selectbox(
            "Genetic Phenotype",
            ["Ultrarapid Metabolizer", "Normal Metabolizer",
             "Intermediate Metabolizer", "Poor Metabolizer"],
            key="pc_phenotype",
        )

    with pc_col2:
        st.markdown("#### Concurrent Medications")
        pc_meds = st.text_area(
            "Medications (one per line)",
            placeholder="fluoxetine\nbupropion\nparoxetine",
            height=120,
            key="pc_meds",
        )

    if st.button("Model Phenoconversion", key="pc_model", type="primary"):
        meds = [m.strip().lower() for m in pc_meds.split("\n") if m.strip()] if pc_meds else []

        detector = pipeline.get("phenoconversion_detector")
        if detector and meds:
            effective = detector.get_effective_phenotype(pc_gene, pc_phenotype, meds)

            # Visualization
            is_converted = effective != pc_phenotype
            st.markdown("### Phenoconversion Result")

            pc_r1, pc_r2, pc_r3 = st.columns([2, 1, 2])
            with pc_r1:
                st.markdown(
                    f'<div class="dose-result" style="text-align:center;">'
                    f'<div class="dose-label">Genetic Phenotype</div>'
                    f'<div class="dose-value" style="font-size:1.2rem;">{pc_phenotype}</div>'
                    f'</div>',
                    unsafe_allow_html=True,
                )
            with pc_r2:
                if is_converted:
                    st.markdown(
                        '<div style="text-align:center; margin-top:20px;">'
                        '<span class="phenoconv-arrow">&rarr;</span>'
                        '</div>',
                        unsafe_allow_html=True,
                    )
                else:
                    st.markdown(
                        '<div style="text-align:center; margin-top:20px;">'
                        '<span style="color:#76B900; font-size:1.5rem;">=</span>'
                        '</div>',
                        unsafe_allow_html=True,
                    )
            with pc_r3:
                border_color = "#F9C500" if is_converted else "#76B900"
                bg_color = "#3a2a1a" if is_converted else "#1a2a1a"
                st.markdown(
                    f'<div style="background:{bg_color}; border:1px solid {border_color}; '
                    f'border-radius:8px; padding:16px; text-align:center;">'
                    f'<div class="dose-label">Effective Phenotype</div>'
                    f'<div class="dose-value" style="font-size:1.2rem; color:{border_color};">'
                    f'{effective}</div>'
                    f'</div>',
                    unsafe_allow_html=True,
                )

            if is_converted:
                st.warning(
                    f"Phenoconversion detected: {pc_gene} {pc_phenotype} is "
                    f"effectively converted to {effective} due to concomitant "
                    f"medication(s). Adjust dosing of {pc_gene} substrates "
                    f"accordingly."
                )
            else:
                st.success(
                    f"No phenoconversion detected. The effective phenotype for "
                    f"{pc_gene} remains {pc_phenotype}."
                )

            # Show which medications are causing the effect
            st.markdown("#### Medication Analysis")
            from src.phenoconversion import CYP_INHIBITORS, CYP_INDUCERS
            for med in meds:
                inhibitor_info = CYP_INHIBITORS.get(med, {}).get(pc_gene)
                inducer_info = CYP_INDUCERS.get(med, {}).get(pc_gene)
                if inhibitor_info:
                    st.markdown(
                        f"- **{med}**: {inhibitor_info.value} **inhibitor** of {pc_gene}"
                    )
                elif inducer_info:
                    st.markdown(
                        f"- **{med}**: {inducer_info.value} **inducer** of {pc_gene}"
                    )
                else:
                    st.markdown(
                        f"- **{med}**: No known effect on {pc_gene}"
                    )

            # Show affected substrates
            from src.phenoconversion import CYP_SUBSTRATES
            substrates = CYP_SUBSTRATES.get(pc_gene, [])
            affected_subs = [m for m in meds if m in substrates]
            if affected_subs:
                st.markdown("#### Co-prescribed Substrates Affected")
                for sub in affected_subs:
                    st.markdown(
                        f'<div class="alert-warning">'
                        f'<strong>{sub}</strong> is a {pc_gene} substrate. '
                        f'Consider dosing based on effective phenotype ({effective}) '
                        f'rather than genetic phenotype ({pc_phenotype}).'
                        f'</div>',
                        unsafe_allow_html=True,
                    )
        elif not meds:
            st.warning("Please enter at least one medication.")
        else:
            st.warning("Phenoconversion detector not available.")

    # Reference: CYP Inhibitors and Inducers
    st.markdown("---")
    st.markdown("### Reference: Known CYP Inhibitors and Inducers")
    if knowledge_data and knowledge_data.get("pharmacogenes"):
        try:
            from src.phenoconversion import CYP_INHIBITORS, CYP_INDUCERS

            ref_tab1, ref_tab2 = st.tabs(["Inhibitors", "Inducers"])
            with ref_tab1:
                for gene in ["CYP2D6", "CYP2C19", "CYP2C9", "CYP3A4", "CYP3A5"]:
                    inhibitors = {
                        drug: info[gene].value
                        for drug, info in CYP_INHIBITORS.items()
                        if gene in info
                    }
                    if inhibitors:
                        with st.expander(f"{gene} Inhibitors ({len(inhibitors)})"):
                            for drug, strength in sorted(inhibitors.items()):
                                color = {"strong": "#E53935", "moderate": "#F9C500", "weak": "#a0a0b0"}.get(strength, "#a0a0b0")
                                st.markdown(
                                    f"- <span style='color:{color};font-weight:600;'>[{strength}]</span> {drug}",
                                    unsafe_allow_html=True,
                                )

            with ref_tab2:
                for gene in ["CYP2C19", "CYP3A4", "CYP3A5"]:
                    inducers = {
                        drug: info[gene].value
                        for drug, info in CYP_INDUCERS.items()
                        if gene in info
                    }
                    if inducers:
                        with st.expander(f"{gene} Inducers ({len(inducers)})"):
                            for drug, strength in sorted(inducers.items()):
                                color = {"strong": "#E53935", "moderate": "#F9C500", "weak": "#a0a0b0"}.get(strength, "#a0a0b0")
                                st.markdown(
                                    f"- <span style='color:{color};font-weight:600;'>[{strength}]</span> {drug}",
                                    unsafe_allow_html=True,
                                )
        except ImportError:
            st.info("Phenoconversion module not available for reference display.")


# ═══════════════════════════════════════════════════════════════════════
# TAB 10: Population Analytics
# ═══════════════════════════════════════════════════════════════════════

with tab_population:
    st.markdown("## Population Pharmacogenetics Analytics")
    st.markdown(
        "Explore population-specific allele frequencies for key pharmacogenes. "
        "Understanding ancestry-based allele distribution is critical for "
        "equitable precision prescribing."
    )

    pop_gene = st.selectbox(
        "Select Pharmacogene",
        get_gene_list(),
        key="pop_gene",
    )

    if pop_gene and knowledge_data:
        gene_data = knowledge_data.get("pharmacogenes", {}).get(pop_gene)
        if gene_data:
            st.markdown(f"### {pop_gene} — {gene_data.get('full_name', '')}")
            st.markdown(f"**Chromosome:** {gene_data.get('chromosome', 'N/A')}")
            st.markdown(f"**Function:** {gene_data.get('function', 'N/A')}")
            st.markdown(f"**Complexity:** {gene_data.get('complexity_level', 'N/A')}")

            if gene_data.get("star_alleles_defined"):
                st.markdown(f"**Star alleles defined:** {gene_data['star_alleles_defined']}")
            if gene_data.get("cpic_guidelines"):
                st.markdown(f"**CPIC guidelines:** {', '.join(gene_data['cpic_guidelines'])}")

            # Key variants
            variants = gene_data.get("key_variants", [])
            if variants:
                st.markdown("#### Key Variants")
                for v in variants:
                    st.markdown(f"- {v}")

            # Population frequency data
            pop_freq = knowledge_data.get("population_frequencies", {}).get(pop_gene, {})
            if pop_freq:
                st.markdown("#### Allele Frequencies by Population")

                # Build frequency table
                freq_data = []
                populations = set()
                alleles = set()
                for allele, pop_data in pop_freq.items():
                    alleles.add(allele)
                    for pop, freq in pop_data.items():
                        populations.add(pop)

                for allele in sorted(alleles):
                    row = {"Allele": allele}
                    for pop in sorted(populations):
                        row[pop] = pop_freq.get(allele, {}).get(pop, "N/A")
                    freq_data.append(row)

                if freq_data:
                    st.dataframe(freq_data, use_container_width=True, hide_index=True)
                else:
                    st.info(f"No detailed population frequency data available for {pop_gene}.")

                # Visualization with bar chart
                try:
                    import pandas as pd
                    if freq_data:
                        df = pd.DataFrame(freq_data)
                        df = df.set_index("Allele")
                        # Convert string percentages to float if needed
                        for col in df.columns:
                            try:
                                df[col] = df[col].astype(str).str.rstrip("%").astype(float) / 100
                            except (ValueError, TypeError):
                                pass
                        numeric_cols = df.select_dtypes(include=["float64", "int64"])
                        if not numeric_cols.empty:
                            st.markdown("#### Frequency Visualization")
                            st.bar_chart(numeric_cols)
                except Exception:
                    pass
            else:
                # Show inline allele frequencies from gene data if available
                inline_freqs = {}
                for key in ["allele_frequency_global", "allele_frequency_european",
                             "allele_frequency_african", "allele_frequency_east_asian",
                             "allele_frequency_south_asian", "allele_frequency_latino"]:
                    val = gene_data.get(key)
                    if val is not None:
                        label = key.replace("allele_frequency_", "").replace("_", " ").title()
                        inline_freqs[label] = val

                if inline_freqs:
                    st.markdown("#### Reference Allele Frequencies")
                    freq_items = [{"Population": k, "Frequency": v} for k, v in inline_freqs.items()]
                    st.dataframe(freq_items, use_container_width=True, hide_index=True)
                else:
                    st.info(f"No population frequency data available for {pop_gene} in knowledge graph.")

        # RAG-based population query
        st.markdown("---")
        pop_query = st.text_input(
            "Ask about population-specific PGx data",
            placeholder=f"e.g., What is the frequency of {pop_gene} poor metabolizers in East Asians?",
            key="pop_query",
        )
        if pop_query:
            with st.spinner("Searching population data..."):
                response, evidence = run_query(
                    pop_query,
                    collections_filter=["pgx_population_data", "pgx_gene_reference", "pgx_clinical_evidence"],
                )
            st.markdown(response)
            if evidence:
                with st.expander("Evidence", expanded=False):
                    render_evidence_cards(evidence)

    # Population health equity section
    st.markdown("---")
    st.markdown("### Health Equity in Pharmacogenomics")
    st.markdown(
        "Pharmacogenomic allele frequencies vary significantly across "
        "populations. Clinical implementation must account for these "
        "differences to ensure equitable access to precision prescribing."
    )

    equity_data = [
        {
            "Gene": "CYP2D6",
            "Clinical Impact": "Codeine, tramadol, tamoxifen efficacy",
            "Key Population Consideration": "*10 is very common in East Asians (~40%); *17 is common in Africans (~20%)",
        },
        {
            "Gene": "CYP2C19",
            "Clinical Impact": "Clopidogrel efficacy",
            "Key Population Consideration": "*2 is very common in East Asians (~30% vs ~15% European); *17 is more common in Europeans",
        },
        {
            "Gene": "DPYD",
            "Clinical Impact": "Fluoropyrimidine toxicity",
            "Key Population Consideration": "*2A is predominantly European (~1%); HapB3 also primarily European",
        },
        {
            "Gene": "NUDT15",
            "Clinical Impact": "Thiopurine toxicity",
            "Key Population Consideration": "*3 is common in East Asians (~10%) but rare in Europeans (<1%)",
        },
        {
            "Gene": "HLA-B*15:02",
            "Clinical Impact": "Carbamazepine SJS/TEN",
            "Key Population Consideration": "~8% in Southeast Asians vs <1% in Europeans — FDA mandates testing in at-risk populations",
        },
        {
            "Gene": "HLA-B*58:01",
            "Clinical Impact": "Allopurinol SJS/TEN",
            "Key Population Consideration": "~6-8% in Southeast Asians and Koreans; ~3-4% in African Americans; ~1-2% in Europeans",
        },
        {
            "Gene": "CYP3A5",
            "Clinical Impact": "Tacrolimus dosing",
            "Key Population Consideration": "*1 (expresser) is common in Africans (~70%) vs rare in Europeans (~10%)",
        },
    ]
    st.dataframe(equity_data, use_container_width=True, hide_index=True)


# ═══════════════════════════════════════════════════════════════════════
# FOOTER
# ═══════════════════════════════════════════════════════════════════════

st.markdown("---")
st.markdown(
    '<div style="text-align:center; color:#6a6a7a; font-size:0.8rem;">'
    'Pharmacogenomics Intelligence Agent v1.0 | HCLS AI Factory | '
    'Port 8507 | Powered by NVIDIA DGX Spark + Claude AI<br>'
    '15 Collections | 25 Pharmacogenes | 12 Drug Categories | '
    '12 HLA Associations | BGE-small-en-v1.5 Embeddings<br>'
    '<em>For research and clinical decision support only. '
    'Always verify recommendations with current CPIC/DPWG guidelines.</em>'
    '</div>',
    unsafe_allow_html=True,
)
