"""Benchmarks & Validation tab — AI model performance metrics."""
import streamlit as st
import json
from pathlib import Path
from collections import Counter

BENCHMARK_DATA_PATH = Path(__file__).parent.parent.parent / "data" / "reference" / "benchmark_seed_data.json"

def _load_benchmarks():
    with open(BENCHMARK_DATA_PATH) as f:
        return json.load(f)

def render():
    st.header("Benchmarks & Validation")
    st.caption("AI model performance metrics across architectures, modalities, and deployment targets")

    benchmarks = _load_benchmarks()

    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Benchmarks", len(benchmarks))
    with col2:
        architectures = set(b.get("model_architecture", "") for b in benchmarks)
        st.metric("Architectures", len(architectures))
    with col3:
        dice_scores = [b["metric_value"] for b in benchmarks if b.get("metric_name") == "Dice"]
        if dice_scores:
            st.metric("Best Dice Score", f"{max(dice_scores):.3f}")
    with col4:
        times = [float(b.get("inference_time_ms", "0").replace("ms", "").strip()) for b in benchmarks if b.get("inference_time_ms")]
        if times:
            st.metric("Fastest Inference", f"{min(t for t in times if t > 0):.0f} ms")

    st.divider()

    # Filters
    col1, col2, col3 = st.columns(3)
    with col1:
        all_arch = sorted(set(b.get("model_architecture", "") for b in benchmarks))
        sel_arch = st.multiselect("Architecture", all_arch, default=[])
    with col2:
        all_tasks = sorted(set(b.get("ai_task", "") for b in benchmarks))
        sel_tasks = st.multiselect("AI Task", all_tasks, default=[])
    with col3:
        all_hw = sorted(set(b.get("hardware", "") for b in benchmarks))
        sel_hw = st.multiselect("Hardware", all_hw, default=[])

    filtered = benchmarks
    if sel_arch:
        filtered = [b for b in filtered if b.get("model_architecture") in sel_arch]
    if sel_tasks:
        filtered = [b for b in filtered if b.get("ai_task") in sel_tasks]
    if sel_hw:
        filtered = [b for b in filtered if b.get("hardware") in sel_hw]

    st.markdown(f"**Showing {len(filtered)} of {len(benchmarks)} benchmarks**")

    # Performance chart
    st.subheader("Performance by Model")
    import pandas as pd
    if filtered:
        df_data = []
        for b in filtered:
            df_data.append({
                "Model": b.get("model_name", "")[:25],
                "Metric": b.get("metric_name", ""),
                "Value": b.get("metric_value", 0),
                "Task": b.get("ai_task", ""),
                "Architecture": b.get("model_architecture", ""),
                "Hardware": b.get("hardware", ""),
            })
        df = pd.DataFrame(df_data)
        st.dataframe(df, use_container_width=True, hide_index=True)

    # Distribution charts
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**By Architecture**")
        arch_counts = Counter(b.get("model_architecture", "unknown") for b in filtered)
        st.bar_chart(dict(sorted(arch_counts.items(), key=lambda x: -x[1])))
    with col2:
        st.markdown("**By Hardware**")
        hw_counts = Counter(b.get("hardware", "unknown") for b in filtered)
        st.bar_chart(dict(sorted(hw_counts.items(), key=lambda x: -x[1])))

    # Detail expanders
    st.subheader("Benchmark Details")
    for b in filtered:
        val = b.get("metric_value", 0)
        quality = "🟢" if val >= 0.9 else "🟡" if val >= 0.8 else "🔵"
        with st.expander(f"{quality} **{b.get('model_name', 'Unknown')}** — {b.get('metric_name', '')}: {val:.3f}"):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"**Architecture:** {b.get('model_architecture', 'N/A')}")
                st.markdown(f"**AI Task:** {b.get('ai_task', 'N/A')}")
                st.markdown(f"**Modality:** {b.get('modality', 'N/A').upper()}")
                st.markdown(f"**Body Region:** {b.get('body_region', 'N/A')}")
            with col2:
                st.markdown(f"**Dataset:** {b.get('dataset_name', 'N/A')}")
                st.markdown(f"**Training Data:** {b.get('training_data_size', 'N/A')}")
                st.markdown(f"**Inference Time:** {b.get('inference_time_ms', 'N/A')}")
                st.markdown(f"**Hardware:** {b.get('hardware', 'N/A')}")
            st.markdown(f"_{b.get('text_summary', '')}_")
