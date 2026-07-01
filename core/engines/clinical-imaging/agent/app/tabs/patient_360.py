"""Patient 360 / Cross-Modal tab — imaging + genomics integration."""
import streamlit as st
import json
import math
from pathlib import Path

DEMO_CASES_PATH = Path(__file__).parent.parent.parent / "data" / "reference" / "demo_cases.json"


def _try_live_genomic_query(genes: list, cross_modal_queries: list) -> dict:
    """Attempt live cross-modal query against genomic_evidence collection.
    Falls back to None if Milvus/collection unavailable."""
    try:
        from pymilvus import connections, Collection

        # Check if genomic_evidence collection exists
        connections.connect(alias="default", host="localhost", port="19530")

        if not Collection("genomic_evidence").num_entities:
            return None

        # Use the RAG engine if available
        engine = st.session_state.get("rag_engine")
        if engine is None:
            return None

        # Run cross-modal queries and collect results
        live_results = []
        for query in cross_modal_queries[:3]:  # Limit to 3 queries
            try:
                result = engine.retrieve(query, top_k_per_collection=3, collections=["genomic_evidence"])
                if hasattr(result, 'hits') and result.hits:
                    for hit in result.hits[:3]:
                        live_results.append({
                            "query": query,
                            "text": hit.text[:300] if hit.text else "",
                            "score": hit.score,
                            "id": hit.id,
                        })
            except Exception:
                continue

        return {"live_results": live_results, "query_count": len(cross_modal_queries)} if live_results else None
    except Exception:
        return None


def _render_cross_modal_graph(case_data: dict):
    """Render an interactive network graph of imaging-genomic relationships."""
    genomic = case_data.get("genomic_context", {})
    genes = genomic.get("genes", [])
    overrides = case_data.get("workflow_overrides", {})

    if not genes:
        return

    st.subheader("Cross-Modal Relationship Graph")

    try:
        import plotly.graph_objects as go

        # Build nodes and edges
        nodes = []
        edges = []
        node_colors = []
        node_sizes = []
        node_labels = []

        # Center node: the clinical case
        case_title_short = case_data.get("title", "Case").split(":")[0].strip()
        nodes.append({"id": "case", "label": case_title_short, "x": 0, "y": 0})
        node_colors.append("#76B900")  # NVIDIA green
        node_sizes.append(40)
        node_labels.append(case_title_short)

        # Imaging findings ring
        findings = []
        if "hemorrhage_type" in overrides:
            findings.extend([
                "Hemorrhage",
                f"Volume: {overrides.get('volume_ml', '?')} mL",
                f"Midline Shift: {overrides.get('midline_shift_mm', '?')} mm",
            ])
        elif "nodules" in overrides:
            findings.extend([
                "Lung Nodule",
                f"Lung-RADS: {overrides.get('nodules', [{}])[0].get('lung_rads', '?')}",
                f"Size: {overrides.get('nodules', [{}])[0].get('long_axis_mm', '?')} mm",
            ])
        elif "calcium_score" in overrides:
            findings.extend([
                "Coronary Stenosis",
                f"Ca Score: {overrides.get('calcium_score', '?')}",
                f"CAD-RADS: {case_data.get('expected_classification', '?')}",
            ])
        elif "findings" in overrides and isinstance(overrides["findings"], list):
            for f in overrides["findings"][:3]:
                if isinstance(f, dict):
                    findings.append(f.get("finding", "Finding").replace("_", " ").title())

        if not findings:
            findings = ["Finding 1", "Finding 2"]

        n_findings = len(findings)
        for i, finding in enumerate(findings):
            angle = (2 * math.pi * i / n_findings) - math.pi / 2
            x = 2.0 * math.cos(angle)
            y = 2.0 * math.sin(angle) - 1.0
            nodes.append({"id": f"f_{i}", "label": finding, "x": x, "y": y})
            node_colors.append("#FF8C00")  # Orange for findings
            node_sizes.append(25)
            node_labels.append(finding)
            edges.append(("case", f"f_{i}"))

        # Gene ring (outer)
        n_genes = len(genes)
        for i, gene in enumerate(genes):
            angle = (2 * math.pi * i / n_genes) + math.pi / 4
            x = 4.0 * math.cos(angle)
            y = 4.0 * math.sin(angle)
            nodes.append({"id": f"g_{i}", "label": gene, "x": x, "y": y})
            node_colors.append("#00D4FF")  # Cyan for genes
            node_sizes.append(20)
            node_labels.append(gene)
            # Connect genes to relevant findings (connect to nearest finding)
            nearest_finding = f"f_{i % n_findings}"
            edges.append((nearest_finding, f"g_{i}"))

        # Build plotly figure
        edge_x, edge_y = [], []
        node_map = {n["id"]: (n["x"], n["y"]) for n in nodes}
        for src, dst in edges:
            x0, y0 = node_map[src]
            x1, y1 = node_map[dst]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y, mode='lines',
            line=dict(width=1.5, color='rgba(150,150,150,0.5)'),
            hoverinfo='none'
        )

        node_x = [n["x"] for n in nodes]
        node_y = [n["y"] for n in nodes]

        node_trace = go.Scatter(
            x=node_x, y=node_y, mode='markers+text',
            marker=dict(size=node_sizes, color=node_colors, line=dict(width=2, color='white')),
            text=node_labels,
            textposition="top center",
            textfont=dict(size=11, color="white"),
            hoverinfo='text',
        )

        fig = go.Figure(data=[edge_trace, node_trace])
        fig.update_layout(
            showlegend=False,
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            font_color="white",
            height=450,
            margin=dict(l=20, r=20, t=20, b=20),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        )

        # Add legend annotations
        fig.add_annotation(
            x=0.02, y=0.98, xref="paper", yref="paper",
            text="<b>Case</b>  <b style='color:#FF8C00'>Findings</b>  <b style='color:#00D4FF'>Genes</b>",
            showarrow=False, font=dict(size=12, color="white"), align="left",
        )

        st.plotly_chart(fig, use_container_width=True)

    except ImportError:
        # Fallback: text-based representation
        st.markdown("**Cross-Modal Links:**")
        for gene in genes:
            st.markdown(f"- **{gene}** -- {case_data.get('title', 'Case').split(':')[0]}")


def _load_demo_cases():
    try:
        with open(DEMO_CASES_PATH) as f:
            return json.load(f)
    except Exception:
        return []

def render():
    st.header("Patient 360 — Cross-Modal Integration")
    st.caption("See how imaging findings connect to genomic intelligence across the HCLS AI Factory pipeline")

    cases = _load_demo_cases()
    if not cases:
        st.warning("No demo cases available.")
        return

    # Case selector
    case_options = {c["case_id"]: f"{c['case_id']}: {c['title']}" for c in cases}
    selected_id = st.selectbox("Select Patient Case", list(case_options.keys()), format_func=lambda x: case_options[x])
    case = next(c for c in cases if c["case_id"] == selected_id)

    # Patient header
    demo = case["patient_demographics"]
    st.markdown(f"### {case['title']}")

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Age / Sex", f"{demo.get('age')}y {demo.get('sex', 'N/A')}")
    with col2:
        st.metric("Weight", f"{demo.get('weight_kg', 'N/A')} kg")
    with col3:
        st.metric("Workflow", case["workflow_name"])
    with col4:
        sev = case.get("expected_severity", "N/A")
        colors = {"critical": "🔴", "urgent": "🟠", "significant": "🟡", "routine": "🟢"}
        st.metric("Severity", f"{colors.get(sev, '⚪')} {sev.upper()}")

    st.markdown(f"**Clinical Scenario:** {case['clinical_scenario']}")

    st.divider()

    # Two-column layout: Imaging | Genomics
    col_img, col_gen = st.columns(2)

    with col_img:
        st.subheader("🔬 Imaging Analysis")
        st.markdown(f"**Workflow:** `{case['workflow_name']}`")

        overrides = case.get("workflow_overrides", {})

        # Render findings based on workflow type
        if "hemorrhage" in case["workflow_name"]:
            st.markdown("**Hemorrhage Detection**")
            st.markdown(f"- Type: {overrides.get('hemorrhage_type', 'N/A')}")
            st.markdown(f"- Location: {overrides.get('location', 'N/A')}")
            st.markdown(f"- Volume: {overrides.get('volume_ml', 'N/A')} mL")
            st.markdown(f"- Midline Shift: {overrides.get('midline_shift_mm', 'N/A')} mm")
            st.markdown(f"- IVE: {'Yes' if overrides.get('intraventricular_extension') else 'No'}")
        elif "lung_nodule" in case["workflow_name"]:
            st.markdown("**Nodule Detection**")
            nodules = overrides.get("nodules", [])
            for i, n in enumerate(nodules):
                st.markdown(f"**Nodule {i+1}:** {n.get('location', 'N/A')}")
                st.markdown(f"- Size: {n.get('long_axis_mm', 'N/A')} × {n.get('short_axis_mm', 'N/A')} mm")
                st.markdown(f"- Type: {n.get('type', 'N/A')} | Margin: {n.get('margin', 'N/A')}")
                st.markdown(f"- Lung-RADS: **{n.get('lung_rads', 'N/A')}**")
                if n.get("doubling_time_days"):
                    st.markdown(f"- Volume Doubling: {n['doubling_time_days']} days")
        elif "coronary" in case["workflow_name"]:
            st.markdown("**Coronary CTA Analysis**")
            st.markdown(f"- Calcium Score (Agatston): **{overrides.get('calcium_score_agatston', 'N/A')}**")
            st.markdown(f"- Calcium Percentile: {overrides.get('calcium_percentile', 'N/A')}th")
            st.markdown(f"- EF: {overrides.get('ejection_fraction_pct', 'N/A')}%")
            vessels = overrides.get("vessels", {})
            for vessel, info in vessels.items():
                if isinstance(info, dict) and info.get("stenosis_pct", 0) > 0:
                    st.markdown(f"- **{vessel.replace('_', ' ').title()}:** {info['stenosis_pct']}% stenosis ({info.get('plaque_type', 'N/A')} plaque)")
            st.markdown(f"- CAD-RADS: **{overrides.get('cadrads', 'N/A')}**")
        else:
            for k, v in overrides.items():
                if isinstance(v, (str, int, float, bool)):
                    st.markdown(f"- **{k.replace('_', ' ').title()}:** {v}")

        st.markdown(f"**Classification:** {case.get('expected_classification', 'N/A')}")

    with col_gen:
        st.subheader("🧬 Genomic Enrichment")
        genomic = case.get("genomic_context")
        if genomic:
            genes = genomic.get("genes", [])
            cross_modal_queries = genomic.get("cross_modal_queries", [])
            st.markdown(f"**Linked Genes:** {', '.join(f'`{g}`' for g in genes)}")
            st.markdown(f"**Clinical Relevance:**")
            st.markdown(genomic.get("relevance", ""))

            st.markdown("**Cross-Modal Queries:**")
            for q in cross_modal_queries:
                st.markdown(f"- _{q}_")

            # Live cross-modal genomic query with fallback
            st.markdown("---")
            live_data = _try_live_genomic_query(genes, cross_modal_queries)
            if live_data and live_data["live_results"]:
                st.markdown(
                    '<span style="background:#76B900;color:white;padding:4px 12px;'
                    'border-radius:4px;font-weight:bold;font-size:12px;">'
                    'LIVE CROSS-MODAL</span>',
                    unsafe_allow_html=True,
                )
                st.caption(
                    f"Queried genomic_evidence collection "
                    f"({len(live_data['live_results'])} results)"
                )
                for r in live_data["live_results"]:
                    with st.container(border=True):
                        c1, c2 = st.columns([3, 1])
                        with c1:
                            st.markdown(f"**Query:** {r['query']}")
                            display_text = (
                                r['text'][:200] + "..."
                                if len(r['text']) > 200
                                else r['text']
                            )
                            st.markdown(display_text)
                        with c2:
                            st.metric("Similarity", f"{r['score']:.3f}")
            else:
                st.markdown(
                    '<span style="background:#555;color:white;padding:4px 12px;'
                    'border-radius:4px;font-weight:bold;font-size:12px;">'
                    'DEMO MODE</span>',
                    unsafe_allow_html=True,
                )
                st.caption(
                    "Connect genomic_evidence collection for live cross-modal enrichment"
                )
        else:
            st.info("No genomic context for this case.")

    st.divider()

    # Cross-modal relationship network graph
    _render_cross_modal_graph(case)

    # Pipeline visualization
    st.subheader("HCLS AI Factory Pipeline")

    if "hemorrhage" in case["workflow_name"]:
        pipeline = "DICOM CT Head → AI Hemorrhage Detection → Severity Grading → Genomic Risk Factors → Structured ICH Report"
    elif "lung_nodule" in case["workflow_name"]:
        pipeline = "Low-Dose CT → AI Nodule Detection → Lung-RADS Classification → Cancer Genomics (EGFR/ALK/ROS1) → Targeted Therapy Planning"
    elif "coronary" in case["workflow_name"]:
        pipeline = "Coronary CTA → Calcium Scoring → AI Stenosis Grading → CAD-RADS → Cardiovascular Genomics (LDLR/PCSK9) → Risk Management"
    else:
        pipeline = "Imaging Study → AI Analysis → Cross-Modal Genomics → Clinical Report"

    steps = pipeline.split(" → ")
    cols = st.columns(len(steps))
    for i, (col, step) in enumerate(zip(cols, steps)):
        with col:
            stage_colors = ["🟦", "🟩", "🟨", "🟧", "🟥"]
            color = stage_colors[i % len(stage_colors)]
            st.markdown(f"**{color} Step {i+1}**")
            st.markdown(f"_{step}_")
            if i < len(steps) - 1:
                st.markdown("→")

    # Talking points
    talking_points = case.get("talking_points", [])
    if talking_points:
        st.divider()
        st.subheader("💬 Demo Talking Points")
        for tp in talking_points:
            st.success(tp)
