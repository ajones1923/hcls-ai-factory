"""
Precision Autoimmune Intelligence Agent — Streamlit UI

10-tab interface:
  1. Clinical Query        — RAG-powered Q&A with evidence citations
  2. Patient Analysis      — Full autoimmune analysis pipeline
  3. Document Ingest       — Upload clinical PDFs for patients
  4. Diagnostic Odyssey    — Timeline visualization and odyssey analysis
  5. Autoantibody Panel    — Interactive antibody interpretation
  6. HLA Analysis          — HLA-disease association lookup
  7. Disease Activity      — Activity scoring dashboards
  8. Flare Prediction      — Biomarker-based flare risk
  9. Therapy Advisor       — Biologic therapy recommendations with PGx
 10. Knowledge Base        — Collection stats, evidence explorer

Port: 8531
Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import os
import sys
import tempfile
import time
from pathlib import Path

import streamlit as st

# ── Page config (must be first Streamlit call) ────────────────────────────
st.set_page_config(
    page_title="Precision Autoimmune Intelligence",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── NVIDIA theme CSS ──────────────────────────────────────────────────────
st.markdown("""
<style>
    .stApp { background-color: #1a1a2e; }
    .stTabs [data-baseweb="tab-list"] { gap: 2px; }
    .stTabs [data-baseweb="tab"] {
        background-color: #16213e; border-radius: 4px 4px 0 0;
        padding: 8px 16px; color: #76B900;
    }
    .stTabs [aria-selected="true"] {
        background-color: #0f3460; border-bottom: 2px solid #76B900;
    }
    .stMetric { background-color: #16213e; border-radius: 8px; padding: 12px; }
    div[data-testid="stMetricValue"] { color: #76B900; }
    .stAlert { border-radius: 8px; }
    .stExpander { background-color: #16213e; border-radius: 8px; }
</style>
""", unsafe_allow_html=True)

st.warning(
    "**Clinical Decision Support Tool** — This system provides evidence-based guidance "
    "for research and clinical decision support only. All recommendations must be verified "
    "by a qualified healthcare professional. Not FDA-cleared. Not a substitute for professional "
    "clinical judgment."
)

# ── Path setup ────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.logging import configure_logging
from config.settings import settings
from src.models import AutoimmuneDisease, AutoimmunePatientProfile

# ── Logging setup ────────────────────────────────────────────────────────
configure_logging(
    log_level=settings.LOG_LEVEL,
    log_dir=settings.LOG_DIR or None,
    service_name="autoimmune-ui",
)

# ── Cached initialization ─────────────────────────────────────────────────

@st.cache_resource(ttl=600)
def init_services():
    """Initialize all backend services (cached)."""
    from src.agent import AutoimmuneAgent
    from src.collections import AutoimmuneCollectionManager
    from src.diagnostic_engine import DiagnosticEngine
    from src.document_processor import DocumentProcessor
    from src.knowledge import KNOWLEDGE_VERSION
    from src.rag_engine import AutoimmuneRAGEngine
    from src.timeline_builder import TimelineBuilder

    # Collection manager
    cm = AutoimmuneCollectionManager(
        host=settings.MILVUS_HOST,
        port=settings.MILVUS_PORT,
    )
    try:
        cm.connect()
        milvus_ok = True
    except Exception:
        milvus_ok = False

    # Embedder
    embedder = None
    try:
        from sentence_transformers import SentenceTransformer
        embedder = SentenceTransformer(settings.EMBEDDING_MODEL)
    except Exception:
        pass

    # LLM
    llm = None
    api_key = settings.ANTHROPIC_API_KEY or os.environ.get("ANTHROPIC_API_KEY", "")
    if api_key:
        try:
            import anthropic
            llm = anthropic.Anthropic(api_key=api_key)
        except Exception:
            pass

    # RAG engine
    rag = AutoimmuneRAGEngine(
        collection_manager=cm,
        embedder=embedder,
        llm_client=llm,
        settings=settings,
        knowledge=KNOWLEDGE_VERSION,
    )

    # Other engines
    agent = AutoimmuneAgent()
    doc_proc = DocumentProcessor(collection_manager=cm, embedder=embedder)
    diag = DiagnosticEngine(agent=agent, rag_engine=rag)
    timeline = TimelineBuilder(rag_engine=rag)

    return {
        "cm": cm,
        "embedder": embedder,
        "llm": llm,
        "rag": rag,
        "agent": agent,
        "doc_processor": doc_proc,
        "diagnostic": diag,
        "timeline": timeline,
        "milvus_ok": milvus_ok,
    }


services = init_services()


# ── Sidebar ───────────────────────────────────────────────────────────────

def render_sidebar():
    with st.sidebar:
        st.image(
            "https://upload.wikimedia.org/wikipedia/sco/thumb/2/21/Nvidia_logo.svg/351px-Nvidia_logo.svg.png",
            width=120,
        )
        st.markdown("### Precision Autoimmune Intelligence")
        st.markdown("---")

        # Service status
        st.markdown("**Service Status**")
        milvus_ok = services.get("milvus_ok", False)
        embedder_ok = services.get("embedder") is not None
        llm_ok = services.get("llm") is not None

        st.markdown(f"{'🟢' if milvus_ok else '🔴'} Milvus Vector DB")
        st.markdown(f"{'🟢' if embedder_ok else '🔴'} BGE Embedder")
        st.markdown(f"{'🟢' if llm_ok else '🔴'} Claude LLM")

        if milvus_ok:
            try:
                stats = services["cm"].get_collection_stats()
                total = sum(v for v in stats.values() if v > 0)
                st.metric("Collections", len(stats))
                st.metric("Total Vectors", f"{total:,}")
            except Exception:
                pass

        st.markdown("---")
        st.markdown("**Demo Patients**")
        patient = st.selectbox("Load patient", [
            "", "sarah_mitchell", "maya_rodriguez",
            "linda_chen", "david_park", "rachel_thompson",
            "emma_williams", "james_cooper", "karen_foster",
            "michael_torres",
        ])
        if patient:
            st.session_state["active_patient"] = patient
            st.toast(f"Loaded patient: {patient.replace('_', ' ').title()}", icon="🔬")

        # Case summary card
        active = st.session_state.get("active_patient", "")
        if active:
            case_summaries = {
                "sarah_mitchell": {
                    "age": 34, "sex": "F",
                    "dx": "Systemic Lupus (SLE)",
                    "hla": "DRB1*03:01",
                    "key_abs": "ANA+ 1:640, anti-dsDNA+, anti-Smith+",
                    "duration": "3.5 years",
                    "current_tx": "Belimumab + HCQ",
                },
                "maya_rodriguez": {
                    "age": 28, "sex": "F",
                    "dx": "POTS / hEDS / MCAS",
                    "hla": "—",
                    "key_abs": "Elevated tryptase, joint hypermobility",
                    "duration": "4 years",
                    "current_tx": "Fludrocortisone + Cromolyn",
                },
                "linda_chen": {
                    "age": 45, "sex": "F",
                    "dx": "Sjogren's Syndrome",
                    "hla": "DRB1*03:01",
                    "key_abs": "ANA+, anti-SSA/Ro+, RF+",
                    "duration": "7 years",
                    "current_tx": "Pilocarpine + HCQ",
                },
                "david_park": {
                    "age": 45, "sex": "M",
                    "dx": "Ankylosing Spondylitis",
                    "hla": "B*27:05",
                    "key_abs": "HLA-B27+",
                    "duration": "6 years",
                    "current_tx": "Secukinumab",
                },
                "rachel_thompson": {
                    "age": 38, "sex": "F",
                    "dx": "Mixed Connective Tissue Disease (MCTD)",
                    "hla": "DRB1*04:01",
                    "key_abs": "ANA+ speckled, anti-U1 RNP+",
                    "duration": "5 years",
                    "current_tx": "HCQ + low-dose prednisone",
                },
                "emma_williams": {
                    "age": 34, "sex": "F",
                    "dx": "Multiple Sclerosis (RRMS)",
                    "hla": "DRB1*15:01",
                    "key_abs": "Oligoclonal bands+, elevated NfL",
                    "duration": "7 months",
                    "current_tx": "Ocrelizumab",
                },
                "james_cooper": {
                    "age": 19, "sex": "M",
                    "dx": "Type 1 Diabetes + Celiac Disease",
                    "hla": "DQ2, DQ8",
                    "key_abs": "anti-GAD65+, anti-IA-2+, tTG-IgA+",
                    "duration": "6 months",
                    "current_tx": "Insulin pump + GF diet",
                },
                "karen_foster": {
                    "age": 48, "sex": "F",
                    "dx": "Systemic Sclerosis (dcSSc)",
                    "hla": "DPB1*13:01",
                    "key_abs": "ANA+, anti-Scl-70+",
                    "duration": "1 year",
                    "current_tx": "Tocilizumab + MMF",
                },
                "michael_torres": {
                    "age": 41, "sex": "M",
                    "dx": "Graves' Disease",
                    "hla": "—",
                    "key_abs": "TRAb+, TSI+, anti-TPO+",
                    "duration": "1 year",
                    "current_tx": "Levothyroxine (post-RAI)",
                },
            }
            summary = case_summaries.get(active)
            if summary:
                st.markdown("---")
                st.markdown(f"**{active.replace('_', ' ').title()}**")
                st.markdown(f"**{summary['age']}{summary['sex']}** | {summary['dx']}")
                st.caption(f"HLA: {summary['hla']}")
                st.caption(f"Abs: {summary['key_abs']}")
                st.caption(f"Duration: {summary['duration']} | Tx: {summary['current_tx']}")

        # Export buttons
        if active:
            st.markdown("---")
            st.markdown("**Export**")
            export_col1, export_col2 = st.columns(2)
            with export_col1:
                if st.button("📄 Markdown", key="export_md_btn"):
                    st.session_state["export_format"] = "markdown"
            with export_col2:
                if st.button("📋 FHIR R4", key="export_fhir_btn"):
                    st.session_state["export_format"] = "fhir"

        st.markdown("---")
        st.caption("HCLS AI Factory | DGX Spark")
        st.caption(f"API: localhost:{settings.API_PORT}")


render_sidebar()


# ── Tab definitions ───────────────────────────────────────────────────────

tabs = st.tabs([
    "Clinical Query",
    "Patient Analysis",
    "Document Ingest",
    "Diagnostic Odyssey",
    "Autoantibody Panel",
    "HLA Analysis",
    "Disease Activity",
    "Flare Prediction",
    "Therapy Advisor",
    "Knowledge Base",
])


# ════════════════════════════════════════════════════════════════════════════
# TAB 1: Clinical Query
# ════════════════════════════════════════════════════════════════════════════
with tabs[0]:
    st.header("Clinical Query")
    st.markdown("Ask questions about autoimmune diseases, treatment options, "
                "diagnostic criteria, and patient-specific data.")

    col1, col2 = st.columns([3, 1])
    with col2:
        patient_filter = st.text_input("Patient ID filter", value=st.session_state.get("active_patient", ""))
        collection_filter = st.multiselect(
            "Filter collections",
            options=settings.all_collection_names,
            default=[],
        )

    with col1:
        question = st.text_area("Enter your clinical question:", height=100,
                                placeholder="e.g., What is the significance of anti-CCP positivity "
                                "with HLA-DRB1*04:01 in a patient with symmetric polyarthritis?")

        if "query_history" not in st.session_state:
            st.session_state["query_history"] = []

        if st.button("Search & Analyze", type="primary", key="query_btn"):
            if not question:
                st.warning("Please enter a question.")
            elif not services.get("rag"):
                st.error("RAG engine not available.")
            else:
                with st.spinner("Searching 14 collections..."):
                    # Evidence retrieval
                    result = services["rag"].retrieve(
                        question,
                        patient_id=patient_filter or None,
                        collections_filter=collection_filter or None,
                    )

                    st.markdown(f"*Found {len(result.hits)} evidence items from "
                                f"{result.total_collections_searched} collections "
                                f"in {result.search_time_ms:.0f}ms*")

                with st.spinner("Synthesizing response with Claude..."):
                    st.markdown("---")
                    try:
                        response_container = st.empty()
                        full_response = []
                        for chunk in services["rag"].query_stream(
                            question,
                            patient_id=patient_filter or None,
                            collections_filter=collection_filter or None,
                        ):
                            full_response.append(chunk)
                            response_container.markdown("".join(full_response))
                        answer = "".join(full_response)
                    except Exception:
                        # Fallback to non-streaming
                        answer = services["rag"].query(
                            question,
                            patient_id=patient_filter or None,
                            collections_filter=collection_filter or None,
                        )
                        st.markdown(answer)

                    st.session_state["last_answer"] = answer
                    st.session_state["query_history"].append({
                        "question": question,
                        "answer": answer[:500],
                        "evidence_count": len(result.hits),
                        "timestamp": time.strftime("%H:%M:%S"),
                    })

                # Show evidence sources
                with st.expander("Evidence Sources", expanded=False):
                    for i, hit in enumerate(result.hits[:15], 1):
                        badge = {"high": "🟢", "medium": "🟡", "low": "🔴"}.get(hit.relevance, "⚪")
                        st.markdown(
                            f"**{i}. {badge} [{hit.collection}]** (score: {hit.score:.3f})\n\n"
                            f"{hit.text[:500]}"
                        )
                        st.markdown("---")

    # Conversation history
    if st.session_state.get("query_history"):
        with st.expander("Query History", expanded=False):
            for i, entry in enumerate(reversed(st.session_state["query_history"][-10:]), 1):
                st.markdown(f"**{entry['timestamp']}** — {entry['question'][:100]}")
                st.caption(f"{entry['evidence_count']} evidence items | {entry['answer'][:200]}...")
                st.markdown("---")


# ════════════════════════════════════════════════════════════════════════════
# TAB 2: Patient Analysis
# ════════════════════════════════════════════════════════════════════════════
with tabs[1]:
    st.header("Patient Analysis Pipeline")
    st.markdown("Run the full autoimmune analysis: autoantibody interpretation, "
                "HLA associations, disease activity, flare risk, therapy recommendations.")

    with st.expander("Sample Patient: RA (52F)", expanded=False):
        if st.button("Load RA Sample Patient"):
            st.session_state["sample_profile"] = {
                "patient_id": "demo_ra_52f",
                "age": 52, "sex": "F",
                "antibodies": [
                    {"antibody": "anti-CCP", "value": 340, "positive": True},
                    {"antibody": "RF", "value": 89, "positive": True},
                    {"antibody": "ANA", "value": 1.2, "positive": False},
                ],
                "hla": {"drb1": ["DRB1*04:01", "DRB1*07:01"]},
                "biomarkers": {"CRP": 12.5, "ESR": 38},
                "conditions": ["rheumatoid_arthritis"],
            }

    with st.expander("Sample Patient: SLE (28F)", expanded=False):
        if st.button("Load SLE Sample Patient"):
            st.session_state["sample_profile"] = {
                "patient_id": "demo_sle_28f",
                "age": 28, "sex": "F",
                "antibodies": [
                    {"antibody": "ANA", "value": 640, "positive": True, "titer": "1:640", "pattern": "homogeneous"},
                    {"antibody": "anti-dsDNA", "value": 120, "positive": True},
                    {"antibody": "anti-Smith", "value": 45, "positive": True},
                ],
                "hla": {"drb1": ["DRB1*03:01", "DRB1*15:01"]},
                "biomarkers": {"CRP": 8.2, "ESR": 55, "complement_C3": 65, "complement_C4": 8},
                "conditions": ["systemic_lupus_erythematosus"],
            }

    # Manual entry
    st.subheader("Patient Profile")
    col1, col2, col3 = st.columns(3)
    with col1:
        pid = st.text_input("Patient ID", value=st.session_state.get("sample_profile", {}).get("patient_id", ""))
    with col2:
        age = st.number_input("Age", 0, 150, value=st.session_state.get("sample_profile", {}).get("age", 45))
    with col3:
        sex = st.selectbox("Sex", ["F", "M"])

    st.subheader("Biomarkers")
    bc1, bc2, bc3, bc4 = st.columns(4)
    sample_bio = st.session_state.get("sample_profile", {}).get("biomarkers", {})
    with bc1:
        crp = st.number_input("CRP (mg/L)", 0.0, 500.0, value=float(sample_bio.get("CRP", 0.0)))
    with bc2:
        esr = st.number_input("ESR (mm/hr)", 0.0, 200.0, value=float(sample_bio.get("ESR", 0.0)))
    with bc3:
        c3 = st.number_input("Complement C3", 0.0, 300.0, value=float(sample_bio.get("complement_C3", 0.0)))
    with bc4:
        c4 = st.number_input("Complement C4", 0.0, 100.0, value=float(sample_bio.get("complement_C4", 0.0)))

    st.subheader("Conditions")
    conditions = st.multiselect(
        "Diagnosed conditions",
        [d.value for d in AutoimmuneDisease],
        default=st.session_state.get("sample_profile", {}).get("conditions", []),
        key="analysis_conditions",
    )

    if st.button("Run Analysis", type="primary", key="analyze_btn"):
        if not pid:
            st.warning("Please enter a patient ID.")
        else:
            biomarkers = {}
            if crp > 0: biomarkers["CRP"] = crp  # noqa: E701
            if esr > 0: biomarkers["ESR"] = esr  # noqa: E701
            if c3 > 0: biomarkers["complement_C3"] = c3  # noqa: E701
            if c4 > 0: biomarkers["complement_C4"] = c4  # noqa: E701

            try:
                profile = AutoimmunePatientProfile(
                    patient_id=pid,
                    age=age,
                    sex=sex,
                    biomarkers=biomarkers,
                    diagnosed_conditions=[AutoimmuneDisease(c) for c in conditions] if conditions else [],
                )

                with st.spinner("Running analysis pipeline..."):
                    result = services["agent"].analyze_patient(profile)

                # Display results
                if result.critical_alerts:
                    for alert in result.critical_alerts:
                        st.error(f"⚠️ {alert}")

                col1, col2 = st.columns(2)
                with col1:
                    st.subheader("Disease Activity Scores")
                    for score in result.disease_activity_scores:
                        color = {"remission": "green", "low": "blue", "moderate": "orange",
                                 "high": "red", "very_high": "red"}.get(score.level.value, "gray")
                        st.markdown(
                            f"**{score.score_name}**: {score.score_value} "
                            f"(:{color}[{score.level.value.upper()}])"
                        )

                with col2:
                    st.subheader("Flare Predictions")
                    for pred in result.flare_predictions:
                        color = {"low": "green", "moderate": "orange", "high": "red",
                                 "imminent": "red"}.get(pred.predicted_risk.value, "gray")
                        st.markdown(
                            f"**{pred.disease.value}**: :{color}[{pred.predicted_risk.value.upper()}] "
                            f"(score: {pred.risk_score:.2f})"
                        )
                        if pred.contributing_factors:
                            st.caption("Contributing: " + ", ".join(pred.contributing_factors))

                if result.biologic_recommendations:
                    st.subheader("Biologic Therapy Recommendations")
                    for therapy in result.biologic_recommendations:
                        with st.expander(f"💊 {therapy.drug_name} ({therapy.drug_class})"):
                            st.markdown(f"**Mechanism:** {therapy.mechanism}")
                            if therapy.pgx_considerations:
                                st.markdown("**PGx Considerations:**")
                                for pgx in therapy.pgx_considerations:
                                    st.markdown(f"- {pgx}")
                            if therapy.contraindications:
                                st.markdown("**Contraindications:**")
                                for ci in therapy.contraindications:
                                    st.markdown(f"- ⚠️ {ci}")
                            if therapy.monitoring_requirements:
                                st.markdown("**Monitoring:**")
                                for m in therapy.monitoring_requirements:
                                    st.markdown(f"- 📋 {m}")

            except Exception as exc:
                st.error(f"Analysis failed: {exc}")


# ════════════════════════════════════════════════════════════════════════════
# TAB 3: Document Ingest
# ════════════════════════════════════════════════════════════════════════════
with tabs[2]:
    st.header("Document Ingestion")
    st.markdown("Upload clinical PDFs to build patient knowledge bases. "
                "Documents are automatically classified, chunked, and embedded.")

    col1, col2 = st.columns([2, 1])

    with col1:
        patient_id = st.text_input(
            "Patient ID",
            value=st.session_state.get("active_patient", ""),
            key="ingest_patient_id",
        )
        uploaded_files = st.file_uploader(
            "Upload clinical PDFs",
            type=["pdf"],
            accept_multiple_files=True,
        )

        if st.button("Ingest Uploaded Files", type="primary", key="ingest_btn"):
            if not uploaded_files:
                st.warning("Please upload at least one PDF.")
            elif not patient_id:
                st.warning("Please enter a patient ID.")
            else:
                doc_proc = services.get("doc_processor")
                if not doc_proc:
                    st.error("Document processor not available.")
                else:
                    progress = st.progress(0)
                    total_chunks = 0
                    for i, f in enumerate(uploaded_files):
                        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
                            tmp.write(f.read())
                            tmp_path = Path(tmp.name)

                        try:
                            records = doc_proc.process_pdf(tmp_path, patient_id)
                            if records:
                                records = doc_proc.embed_records(records)
                                count = doc_proc.collection_manager.insert_batch(
                                    settings.COLL_CLINICAL_DOCUMENTS, records
                                )
                                total_chunks += count
                                st.success(f"✅ {f.name}: {count} chunks ({records[0].get('doc_type', '?')})")
                        except Exception as exc:
                            st.error(f"❌ {f.name}: {exc}")
                        finally:
                            tmp_path.unlink(missing_ok=True)

                        progress.progress((i + 1) / len(uploaded_files))

                    st.info(f"Total: {total_chunks} chunks ingested for patient {patient_id}")

    with col2:
        st.subheader("Demo Data")
        st.markdown("Pre-built patient cases with 130 clinical PDFs across 5 patients.")

        demo_patients = ["sarah_mitchell", "maya_rodriguez", "linda_chen", "david_park", "rachel_thompson"]
        for p in demo_patients:
            demo_dir = settings.DEMO_DATA_DIR / p
            if demo_dir.exists():
                pdf_count = len(list(demo_dir.glob("*.pdf")))
                st.markdown(f"📁 **{p.replace('_', ' ').title()}** ({pdf_count} PDFs)")

        if st.button("Ingest All Demo Data", type="secondary", key="demo_ingest_btn"):
            doc_proc = services.get("doc_processor")
            if doc_proc:
                with st.spinner("Ingesting 130 demo PDFs... this may take a few minutes."):
                    results = doc_proc.ingest_demo_data(settings.DEMO_DATA_DIR)
                    for patient, count in results.items():
                        st.success(f"✅ {patient}: {count} chunks")
                    st.info(f"Total: {sum(results.values())} chunks ingested")
            else:
                st.error("Document processor not available.")


# ════════════════════════════════════════════════════════════════════════════
# TAB 4: Diagnostic Odyssey
# ════════════════════════════════════════════════════════════════════════════
with tabs[3]:
    st.header("Diagnostic Odyssey Timeline")
    st.markdown("Visualize a patient's journey from first symptoms through diagnosis. "
                "Identify diagnostic delays, specialist visits, and missed opportunities.")

    patient_id = st.text_input(
        "Patient ID",
        value=st.session_state.get("active_patient", ""),
        key="odyssey_patient",
    )

    if st.button("Build Timeline", type="primary", key="timeline_btn"):
        if not patient_id:
            st.warning("Please enter a patient ID.")
        else:
            rag = services.get("rag")
            timeline_builder = services.get("timeline")
            diag = services.get("diagnostic")

            with st.spinner("Retrieving patient records..."):
                # Search clinical documents for this patient
                hits = rag.search(
                    f"patient history timeline events for {patient_id}",
                    patient_id=patient_id,
                    collections_filter=[settings.COLL_CLINICAL_DOCUMENTS],
                    top_k_per_collection=50,
                )

                chunks = [{"text_chunk": h.text, **h.metadata} for h in hits]

                # Build timeline
                timeline_data = timeline_builder.build_timeline(patient_id, chunks)

                # Analyze odyssey
                odyssey = diag.analyze_diagnostic_odyssey(timeline_data["events"])

            # Display
            col1, col2, col3 = st.columns(3)
            with col1:
                delay = odyssey.get("diagnostic_delay")
                if delay:
                    st.metric("Diagnostic Delay", f"{delay['years']:.1f} years")
                else:
                    st.metric("Diagnostic Delay", "N/A")
            with col2:
                st.metric("Specialists Seen", odyssey.get("num_specialists_seen", 0))
            with col3:
                st.metric("Misdiagnoses", odyssey.get("num_misdiagnoses", 0))

            st.subheader("Timeline Events")
            for event in timeline_data["events"]:
                icon = {
                    "symptom_onset": "🟡", "diagnosis": "✅", "misdiagnosis": "❌",
                    "lab_result": "🧪", "imaging": "📷", "biopsy": "🔬",
                    "genetic_test": "🧬", "treatment_start": "💊",
                    "treatment_change": "🔄", "flare": "🔥",
                    "referral": "📋", "er_visit": "🚑",
                }.get(event["event_type"], "📝")

                date = event.get("event_date", "Unknown")
                st.markdown(
                    f"{icon} **{date}** — {event['event_type'].replace('_', ' ').title()}\n\n"
                    f"{event.get('description', '')[:200]}"
                )
                if event.get("specialty"):
                    st.caption(f"Specialty: {event['specialty']} | Provider: {event.get('provider', 'Unknown')}")
                st.markdown("---")


# ════════════════════════════════════════════════════════════════════════════
# TAB 5: Autoantibody Panel
# ════════════════════════════════════════════════════════════════════════════
with tabs[4]:
    st.header("Autoantibody Panel Interpretation")
    st.markdown("Enter autoantibody results for disease association analysis.")

    from src.knowledge import AUTOANTIBODY_DISEASE_MAP
    from src.models import AutoantibodyPanel, AutoantibodyResult

    ab_names = list(AUTOANTIBODY_DISEASE_MAP.keys())

    # Input table
    st.subheader("Enter Results")
    ab_results = []
    cols = st.columns(4)
    with cols[0]:
        st.markdown("**Antibody**")
    with cols[1]:
        st.markdown("**Value**")
    with cols[2]:
        st.markdown("**Positive?**")
    with cols[3]:
        st.markdown("**Titer**")

    for i, ab_name in enumerate(ab_names):
        cols = st.columns(4)
        with cols[0]:
            st.text(ab_name)
        with cols[1]:
            val = st.number_input(f"val_{ab_name}", 0.0, 10000.0, 0.0, key=f"ab_val_{i}", label_visibility="collapsed")
        with cols[2]:
            pos = st.checkbox("Pos", key=f"ab_pos_{i}", label_visibility="collapsed")
        with cols[3]:
            titer = st.text_input(f"titer_{ab_name}", "", key=f"ab_titer_{i}", label_visibility="collapsed")

        if pos or val > 0:
            ab_results.append(AutoantibodyResult(
                antibody=ab_name, value=val, positive=pos,
                titer=titer if titer else None,
            ))

    if st.button("Interpret Panel", type="primary", key="ab_interpret_btn"):
        if not ab_results:
            st.warning("Please mark at least one antibody as positive.")
        else:
            panel = AutoantibodyPanel(
                patient_id=st.session_state.get("active_patient", "manual"),
                collection_date="2026-03-10",
                results=ab_results,
            )

            agent = services.get("agent")
            findings = agent.interpret_autoantibodies(panel)

            st.subheader(f"Findings ({len(findings)} disease associations)")
            for f in findings:
                spec = f.get("specificity", 0)
                badge = "🟢 High" if spec >= 0.95 else "🟡 Moderate" if spec >= 0.85 else "🔴 Low"
                st.markdown(
                    f"**{f['antibody']}** → {f['disease'].replace('_', ' ').title()} "
                    f"(sens={f.get('sensitivity', '?')}, spec={spec}) — Specificity: {badge}"
                )
                if f.get("note"):
                    st.caption(f.get("note"))

            # Differential diagnosis
            diag = services.get("diagnostic")
            positive = [r.antibody for r in ab_results if r.positive]
            differential = diag.generate_differential(positive)

            st.subheader("Differential Diagnosis (ranked)")
            for d in differential[:8]:
                st.markdown(f"**{d['rank']}. {d['disease'].replace('_', ' ').title()}** (score: {d['score']})")
                for ev in d["evidence"]:
                    st.caption(f"  • {ev}")


# ════════════════════════════════════════════════════════════════════════════
# TAB 6: HLA Analysis
# ════════════════════════════════════════════════════════════════════════════
with tabs[5]:
    st.header("HLA-Disease Association Analysis")
    st.markdown("Enter HLA alleles to identify autoimmune disease susceptibility.")

    from src.knowledge import HLA_DISEASE_ASSOCIATIONS
    from src.models import HLAProfile

    col1, col2 = st.columns(2)
    with col1:
        hla_b = st.text_input("HLA-B alleles (comma-separated)", placeholder="B*27:05")
        hla_drb1 = st.text_input("HLA-DRB1 alleles (comma-separated)", placeholder="DRB1*04:01, DRB1*03:01")
    with col2:
        hla_c = st.text_input("HLA-C alleles (comma-separated)", placeholder="C*06:02")
        hla_dqb1 = st.text_input("HLA-DQB1 alleles (comma-separated)", placeholder="DQB1*02:01")

    def parse_alleles(text: str) -> list:
        return [a.strip() for a in text.split(",") if a.strip()]

    if st.button("Analyze HLA Associations", type="primary", key="hla_btn"):
        profile = HLAProfile(
            hla_b=parse_alleles(hla_b),
            hla_c=parse_alleles(hla_c),
            hla_drb1=parse_alleles(hla_drb1),
            hla_dqb1=parse_alleles(hla_dqb1),
        )

        if not profile.all_alleles:
            st.warning("Please enter at least one HLA allele.")
        else:
            agent = services.get("agent")
            assocs = agent.analyze_hla_associations(profile)

            st.subheader(f"Found {len(assocs)} Disease Associations")
            for a in assocs:
                or_val = a["odds_ratio"]
                color = "red" if or_val > 10 else "orange" if or_val > 3 else "blue"
                st.markdown(
                    f"**HLA-{a['allele']}** → {a['disease'].replace('_', ' ').title()} — "
                    f":{color}[OR = {or_val}]"
                )
                if a.get("pmid"):
                    st.caption(f"PMID: {a['pmid']}")
                if a.get("note"):
                    st.caption(a["note"])

            # Reference table
            with st.expander("Full HLA Reference Database"):
                import pandas as pd
                rows = []
                for allele, assocs_list in HLA_DISEASE_ASSOCIATIONS.items():
                    for a in assocs_list:
                        rows.append({
                            "Allele": allele,
                            "Disease": a["disease"].replace("_", " ").title(),
                            "OR": a["odds_ratio"],
                            "PMID": a.get("pmid", ""),
                        })
                df = pd.DataFrame(rows).sort_values("OR", ascending=False)
                st.dataframe(df, use_container_width=True)


# ════════════════════════════════════════════════════════════════════════════
# TAB 7: Disease Activity
# ════════════════════════════════════════════════════════════════════════════
with tabs[6]:
    st.header("Disease Activity Scoring")
    st.markdown("Calculate standardized disease activity scores from biomarkers.")

    from src.knowledge import DISEASE_ACTIVITY_THRESHOLDS

    score_system = st.selectbox(
        "Scoring System",
        list(DISEASE_ACTIVITY_THRESHOLDS.keys()),
    )

    info = DISEASE_ACTIVITY_THRESHOLDS.get(score_system, {})
    st.caption(f"Disease: {info.get('disease', '?').replace('_', ' ').title()} | "
               f"Reference: {info.get('reference', 'N/A')}")

    col1, col2 = st.columns(2)
    with col1:
        crp_val = st.number_input("CRP (mg/L)", 0.0, 500.0, 5.0, key="act_crp")
        esr_val = st.number_input("ESR (mm/hr)", 0.0, 200.0, 20.0, key="act_esr")
    with col2:
        st.markdown("**Thresholds:**")
        thresholds = info.get("thresholds", {})
        for level, val in thresholds.items():
            st.markdown(f"- {level.title()}: ≤ {val}")

    if st.button("Calculate Score", type="primary", key="activity_btn"):
        agent = services.get("agent")
        disease_enum = None
        for d in AutoimmuneDisease:
            if d.value == info.get("disease"):
                disease_enum = d
                break

        if disease_enum:
            biomarkers = {"CRP": crp_val, "ESR": esr_val}
            scores = agent.calculate_disease_activity(biomarkers, [disease_enum])

            for score in scores:
                color = {"remission": "green", "low": "blue", "moderate": "orange",
                         "high": "red"}.get(score.level.value, "gray")
                st.markdown(f"### {score.score_name}: {score.score_value}")
                st.markdown(f"Activity Level: :{color}[**{score.level.value.upper()}**]")

                # Visual gauge
                import plotly.graph_objects as go
                fig = go.Figure(go.Indicator(
                    mode="gauge+number",
                    value=score.score_value,
                    gauge={
                        "axis": {"range": info.get("range", [0, 10])},
                        "bar": {"color": {"remission": "#00ff00", "low": "#76B900",
                                          "moderate": "#ff8800", "high": "#ff0000"}.get(score.level.value, "#888")},
                        "steps": [
                            {"range": [0, thresholds.get("low", 3)], "color": "#1a3a1a"},
                            {"range": [thresholds.get("low", 3), thresholds.get("moderate", 5)], "color": "#3a3a1a"},
                            {"range": [thresholds.get("moderate", 5), info.get("range", [0, 10])[1]], "color": "#3a1a1a"},
                        ],
                    },
                    title={"text": score.score_name},
                ))
                fig.update_layout(height=300, paper_bgcolor="rgba(0,0,0,0)", font_color="white")
                st.plotly_chart(fig, use_container_width=True)


# ════════════════════════════════════════════════════════════════════════════
# TAB 8: Flare Prediction
# ════════════════════════════════════════════════════════════════════════════
with tabs[7]:
    st.header("Flare Risk Prediction")
    st.markdown("Predict flare risk from biomarker patterns and disease trajectory.")

    from src.knowledge import FLARE_BIOMARKER_PATTERNS

    disease = st.selectbox(
        "Disease",
        list(FLARE_BIOMARKER_PATTERNS.keys()),
        format_func=lambda x: x.replace("_", " ").title(),
    )

    pattern = FLARE_BIOMARKER_PATTERNS.get(disease, {})
    markers = pattern.get("early_warning_biomarkers", [])

    st.subheader("Enter Current Biomarker Values")
    biomarker_values = {}
    cols = st.columns(min(len(markers), 5) or 1)
    for i, marker in enumerate(markers):
        with cols[i % len(cols)]:
            val = st.number_input(marker, 0.0, 10000.0, 0.0, key=f"flare_{marker}")
            if val > 0:
                biomarker_values[marker] = val

    if st.button("Predict Flare Risk", type="primary", key="flare_btn"):
        if not biomarker_values:
            st.warning("Please enter at least one biomarker value.")
        else:
            agent = services.get("agent")
            disease_enum = None
            for d in AutoimmuneDisease:
                if d.value == disease:
                    disease_enum = d
                    break

            if disease_enum:
                predictions = agent.predict_flares(biomarker_values, [disease_enum])

                for pred in predictions:
                    color = {"low": "green", "moderate": "orange", "high": "red",
                             "imminent": "red"}.get(pred.predicted_risk.value, "gray")

                    st.markdown(f"### Flare Risk: :{color}[{pred.predicted_risk.value.upper()}]")
                    st.markdown(f"Risk Score: **{pred.risk_score:.2f}** / 1.00")

                    # Visual bar
                    import plotly.graph_objects as go
                    fig = go.Figure(go.Indicator(
                        mode="gauge+number",
                        value=pred.risk_score,
                        gauge={
                            "axis": {"range": [0, 1]},
                            "bar": {"color": {"low": "#00ff00", "moderate": "#ff8800",
                                              "high": "#ff4400", "imminent": "#ff0000"}.get(pred.predicted_risk.value)},
                            "steps": [
                                {"range": [0, 0.4], "color": "#1a3a1a"},
                                {"range": [0.4, 0.6], "color": "#3a3a1a"},
                                {"range": [0.6, 0.8], "color": "#3a2a1a"},
                                {"range": [0.8, 1.0], "color": "#3a1a1a"},
                            ],
                        },
                        title={"text": "Flare Risk Score"},
                    ))
                    fig.update_layout(height=250, paper_bgcolor="rgba(0,0,0,0)", font_color="white")
                    st.plotly_chart(fig, use_container_width=True)

                    col1, col2 = st.columns(2)
                    with col1:
                        if pred.contributing_factors:
                            st.markdown("**Contributing Factors:**")
                            for f in pred.contributing_factors:
                                st.markdown(f"- ⚠️ {f}")
                    with col2:
                        if pred.protective_factors:
                            st.markdown("**Protective Factors:**")
                            for f in pred.protective_factors:
                                st.markdown(f"- ✅ {f}")

                    if pred.recommended_monitoring:
                        st.markdown("**Recommended Monitoring:**")
                        for m in pred.recommended_monitoring:
                            st.markdown(f"- 📋 {m}")


# ════════════════════════════════════════════════════════════════════════════
# TAB 9: Therapy Advisor
# ════════════════════════════════════════════════════════════════════════════
with tabs[8]:
    st.header("Biologic Therapy Advisor")
    st.markdown("Evidence-based biologic therapy recommendations with pharmacogenomic filtering.")

    from src.knowledge import BIOLOGIC_THERAPIES

    conditions = st.multiselect(
        "Select conditions",
        [d.value for d in AutoimmuneDisease],
        format_func=lambda x: x.replace("_", " ").title(),
        key="therapy_conditions",
    )

    if st.button("Get Recommendations", type="primary", key="therapy_btn"):
        if not conditions:
            st.warning("Please select at least one condition.")
        else:
            agent = services.get("agent")
            condition_enums = [AutoimmuneDisease(c) for c in conditions]
            therapies = agent.recommend_biologics(condition_enums, {})

            if not therapies:
                st.info("No biologic therapies found for selected conditions.")
            else:
                st.subheader(f"{len(therapies)} Recommended Therapies")
                for t in therapies:
                    with st.expander(f"💊 {t.drug_name} — {t.drug_class}", expanded=True):
                        col1, col2 = st.columns(2)
                        with col1:
                            st.markdown(f"**Mechanism:** {t.mechanism}")
                            if t.pgx_considerations:
                                st.markdown("**PGx Considerations:**")
                                for pgx in t.pgx_considerations:
                                    st.markdown(f"- 🧬 {pgx}")
                        with col2:
                            if t.contraindications:
                                st.markdown("**Contraindications:**")
                                for ci in t.contraindications:
                                    st.markdown(f"- ⚠️ {ci}")
                            if t.monitoring_requirements:
                                st.markdown("**Monitoring:**")
                                for m in t.monitoring_requirements:
                                    st.markdown(f"- 📋 {m}")

    # Full drug reference
    with st.expander("Full Biologic Therapy Reference"):
        import pandas as pd
        rows = []
        for t in BIOLOGIC_THERAPIES:
            rows.append({
                "Drug": t["drug_name"],
                "Class": t["drug_class"],
                "Diseases": ", ".join(d.replace("_", " ").title() for d in t["indicated_diseases"]),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True)


# ════════════════════════════════════════════════════════════════════════════
# TAB 10: Knowledge Base
# ════════════════════════════════════════════════════════════════════════════
with tabs[9]:
    st.header("Knowledge Base & Collections")

    # Collection stats
    if services.get("milvus_ok"):
        try:
            stats = services["cm"].get_collection_stats()

            st.subheader("Collection Statistics")
            import pandas as pd
            rows = []
            for name, count in sorted(stats.items()):
                label = settings.collection_config.get(name, {}).get("label", name)
                weight = settings.collection_config.get(name, {}).get("weight", 0)
                rows.append({"Collection": name, "Label": label, "Vectors": count, "Weight": weight})

            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True)

            total = sum(v for v in stats.values() if v > 0)
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Collections", len(stats))
            with col2:
                st.metric("Total Vectors", f"{total:,}")
            with col3:
                st.metric("Embedding Model", settings.EMBEDDING_MODEL.split("/")[-1])

        except Exception as exc:
            st.error(f"Failed to load collection stats: {exc}")
    else:
        st.warning("Milvus not connected. Collections unavailable.")

    # Knowledge base version
    st.subheader("Knowledge Base Version")
    from src.knowledge import KNOWLEDGE_VERSION
    st.json(KNOWLEDGE_VERSION)

    # Quick stats
    st.subheader("Reference Data")
    from src.knowledge import (
        AUTOANTIBODY_DISEASE_MAP,
        BIOLOGIC_THERAPIES,
        DISEASE_ACTIVITY_THRESHOLDS,
        FLARE_BIOMARKER_PATTERNS,
        HLA_DISEASE_ASSOCIATIONS,
    )

    col1, col2, col3, col4, col5 = st.columns(5)
    with col1:
        total_hla = sum(len(v) for v in HLA_DISEASE_ASSOCIATIONS.values())
        st.metric("HLA Associations", total_hla)
    with col2:
        st.metric("Autoantibody Types", len(AUTOANTIBODY_DISEASE_MAP))
    with col3:
        st.metric("Biologic Drugs", len(BIOLOGIC_THERAPIES))
    with col4:
        st.metric("Activity Scores", len(DISEASE_ACTIVITY_THRESHOLDS))
    with col5:
        st.metric("Flare Patterns", len(FLARE_BIOMARKER_PATTERNS))

    # Collection management
    st.subheader("Collection Management")
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Create All Collections", key="create_colls_btn"):
            with st.spinner("Creating collections..."):
                try:
                    colls = services["cm"].create_all_collections()
                    st.success(f"Created {len(colls)} collections")
                except Exception as exc:
                    st.error(f"Failed: {exc}")
    with col2:
        if st.button("Recreate Collections (DROP existing)", type="secondary", key="recreate_colls_btn"):
            with st.spinner("Dropping and recreating..."):
                try:
                    colls = services["cm"].create_all_collections(drop_existing=True)
                    st.success(f"Recreated {len(colls)} collections")
                except Exception as exc:
                    st.error(f"Failed: {exc}")

# ── Export handler ─────────────────────────────────────────────────────────
export_format = st.session_state.pop("export_format", None)
if export_format:
    try:
        from src.export import AutoimmuneExporter
        from src.knowledge import KNOWLEDGE_VERSION
        exporter = AutoimmuneExporter(knowledge_version=KNOWLEDGE_VERSION)
        active_patient = st.session_state.get("active_patient", "unknown")

        if export_format == "markdown":
            md = exporter.to_markdown(active_patient, query_answer=st.session_state.get("last_answer"))
            st.download_button(
                "Download Markdown Report",
                data=md,
                file_name=f"autoimmune_report_{active_patient}.md",
                mime="text/markdown",
            )
        elif export_format == "fhir":
            fhir_json = exporter.to_fhir_json(active_patient)
            st.download_button(
                "Download FHIR R4 Bundle",
                data=fhir_json,
                file_name=f"autoimmune_fhir_{active_patient}.json",
                mime="application/json",
            )
    except Exception as exc:
        st.error(f"Export failed: {exc}")
