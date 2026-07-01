"""Evidence Explorer tab — RAG search across imaging knowledge base."""
import streamlit as st
from typing import Optional


def _render_comparative_results(engine, result, query: str):
    """Render comparative analysis results with side-by-side layout."""
    entity_a_label = result.entity_a_resolved.get("full_name", result.entity_a)
    entity_b_label = result.entity_b_resolved.get("full_name", result.entity_b)
    entity_a_type = result.entity_a_resolved.get("type", "unknown")
    entity_b_type = result.entity_b_resolved.get("type", "unknown")

    st.subheader(f"Comparative Analysis: {entity_a_label} vs {entity_b_label}")

    # Entity resolution summary
    col_meta_a, col_meta_b = st.columns(2)
    with col_meta_a:
        st.metric(label=entity_a_label, value=entity_a_type.title())
        st.caption(f"Resolved key: `{result.entity_a_resolved.get('key', result.entity_a)}`")
    with col_meta_b:
        st.metric(label=entity_b_label, value=entity_b_type.title())
        st.caption(f"Resolved key: `{result.entity_b_resolved.get('key', result.entity_b)}`")

    # Search stats
    stats_cols = st.columns(4)
    with stats_cols[0]:
        st.metric("Evidence A", result.evidence_a.hit_count)
    with stats_cols[1]:
        st.metric("Evidence B", result.evidence_b.hit_count)
    with stats_cols[2]:
        st.metric("Shared Evidence", result.shared_hit_count)
    with stats_cols[3]:
        st.metric("Search Time", f"{result.total_search_time_ms:.0f} ms")

    st.divider()

    # AI-synthesized comparative analysis
    with st.spinner("Generating comparative analysis..."):
        try:
            synthesis = engine.query(query)
            st.subheader("AI-Synthesized Comparison")
            st.markdown(synthesis)
        except Exception as e:
            st.error(f"LLM synthesis error: {e}")

    st.divider()

    # Side-by-side evidence panels
    st.subheader("Retrieved Evidence (Side-by-Side)")
    col_a, col_b = st.columns(2)

    with col_a:
        st.markdown(f"#### {entity_a_label}")
        if result.evidence_a.hits:
            for i, hit in enumerate(result.evidence_a.hits[:10]):
                with st.expander(
                    f"[{hit.collection}] {hit.score:.3f} — {hit.text[:60]}..."
                ):
                    st.markdown(hit.text)
                    if hit.metadata:
                        st.json(hit.metadata)
        else:
            st.info("No evidence retrieved.")

    with col_b:
        st.markdown(f"#### {entity_b_label}")
        if result.evidence_b.hits:
            for i, hit in enumerate(result.evidence_b.hits[:10]):
                with st.expander(
                    f"[{hit.collection}] {hit.score:.3f} — {hit.text[:60]}..."
                ):
                    st.markdown(hit.text)
                    if hit.metadata:
                        st.json(hit.metadata)
        else:
            st.info("No evidence retrieved.")

    # Shared evidence section
    if result.shared_evidence:
        st.divider()
        st.subheader(f"Shared Evidence ({len(result.shared_evidence)} sources)")
        st.caption("These sources appeared in both entity searches — high relevance for comparison.")
        for hit in result.shared_evidence:
            with st.expander(
                f"[{hit.collection}] {hit.score:.3f} — {hit.text[:80]}..."
            ):
                st.markdown(hit.text)
                if hit.metadata:
                    st.json(hit.metadata)

    # Domain knowledge context
    if result.comparison_context:
        with st.expander("Domain Knowledge Context", expanded=False):
            st.markdown(result.comparison_context)


def render():
    """Render the Evidence Explorer tab."""
    st.header("Evidence Explorer")
    st.caption("Search across 10 imaging knowledge collections with AI-powered synthesis")

    # Example queries section
    st.markdown("**Try an example query:**")
    example_cols = st.columns(4)
    examples = [
        "What AI models detect intracranial hemorrhage?",
        "Lung-RADS guidelines for nodule management",
        "Compare CT vs MRI for stroke detection",
        "Deep learning for cardiac calcium scoring",
    ]
    for i, (col, example) in enumerate(zip(example_cols, examples)):
        with col:
            if st.button(example, key=f"example_{i}", use_container_width=True):
                st.session_state["evidence_query"] = example

    # Query input
    query = st.text_area("Ask a question about medical imaging",
                         value=st.session_state.get("evidence_query", ""),
                         placeholder="e.g., What are the current guidelines for lung nodule management?",
                         height=80)

    # Filters in expander
    with st.expander("Search Filters", expanded=False):
        col1, col2, col3 = st.columns(3)
        with col1:
            modality = st.selectbox("Modality", ["All", "CT", "MRI", "CXR", "Mammography", "Ultrasound", "PET/CT"])
        with col2:
            body_region = st.selectbox("Body Region", ["All", "Head", "Chest", "Abdomen", "Brain", "Cardiac", "Breast", "Spine", "Pelvis", "Extremity", "Neck", "Whole Body"])
        with col3:
            top_k = st.slider("Results per collection", 1, 10, 5)

    if st.button("Search", type="primary", use_container_width=True):
        if not query:
            st.warning("Please enter a question.")
            return

        with st.spinner("Searching across collections..."):
            try:
                engine = st.session_state.get("rag_engine")
                if engine:
                    mod_filter = None if modality == "All" else modality.lower()
                    region_filter = None if body_region == "All" else body_region.lower()

                    # Detect comparative query using the engine's regex-based detection
                    is_comparative = engine._is_comparative(query)

                    if is_comparative:
                        # Run comparative retrieval pipeline
                        comp_result = engine.retrieve_comparative(
                            query,
                            modality_filter=mod_filter,
                            body_region_filter=region_filter,
                            top_k_per_collection=top_k,
                        )
                        _render_comparative_results(engine, comp_result, query)
                    else:
                        # Standard retrieval + synthesis
                        result = engine.retrieve(
                            query=query,
                            top_k_per_collection=top_k,
                            modality_filter=mod_filter,
                            body_region_filter=region_filter,
                        )

                        # Synthesize with LLM
                        answer = engine.query(query)
                        st.subheader("AI-Synthesized Answer")
                        st.markdown(answer)

                        # Evidence panel
                        if hasattr(result, 'hits') and result.hits:
                            st.subheader(f"Evidence ({len(result.hits)} sources)")
                            for i, hit in enumerate(result.hits):
                                with st.expander(f"[{hit.collection}] Score: {hit.score:.3f} — {hit.text[:80]}..."):
                                    st.markdown(hit.text)
                                    if hit.metadata:
                                        st.json(hit.metadata)

                            # Collection contribution chart
                            collection_counts = {}
                            for hit in result.hits:
                                col_name = hit.collection
                                collection_counts[col_name] = collection_counts.get(col_name, 0) + 1

                            if collection_counts:
                                st.divider()
                                try:
                                    import plotly.express as px
                                    fig = px.pie(
                                        values=list(collection_counts.values()),
                                        names=list(collection_counts.keys()),
                                        title="Evidence Sources by Collection",
                                        color_discrete_sequence=[
                                            "#76B900", "#4CAF50", "#2196F3", "#FF9800",
                                            "#9C27B0", "#F44336", "#00BCD4", "#FFEB3B",
                                            "#795548", "#607D8B",
                                        ],
                                        hole=0.4,
                                    )
                                    fig.update_layout(
                                        plot_bgcolor="rgba(0,0,0,0)",
                                        paper_bgcolor="rgba(0,0,0,0)",
                                        font_color="white",
                                        showlegend=True,
                                        height=350,
                                    )
                                    st.plotly_chart(fig, use_container_width=True)
                                except ImportError:
                                    import pandas as pd
                                    df = pd.DataFrame({
                                        "Collection": list(collection_counts.keys()),
                                        "Count": list(collection_counts.values()),
                                    }).set_index("Collection")
                                    st.subheader("Evidence Sources by Collection")
                                    st.bar_chart(df)
                else:
                    st.info("RAG engine not initialized. Start the API server first.")
            except Exception as e:
                st.error(f"Search error: {e}")
