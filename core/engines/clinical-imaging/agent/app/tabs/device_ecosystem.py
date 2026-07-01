"""Device & AI Ecosystem tab — FDA-cleared imaging AI device dashboard."""
import streamlit as st
import json
from pathlib import Path
from collections import Counter

DEVICE_DATA_PATH = Path(__file__).parent.parent.parent / "data" / "reference" / "device_seed_data.json"

def _load_devices():
    with open(DEVICE_DATA_PATH) as f:
        return json.load(f)

def render():
    st.header("Device & AI Ecosystem")
    st.caption("Explore 50 FDA-cleared medical imaging AI devices across modalities, body regions, and AI tasks")

    devices = _load_devices()

    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Devices", len(devices))
    with col2:
        modalities = set(d.get("modality", "") for d in devices)
        st.metric("Modalities", len(modalities))
    with col3:
        tasks = set(d.get("ai_task", "") for d in devices)
        st.metric("AI Task Types", len(tasks))
    with col4:
        cleared = sum(1 for d in devices if d.get("regulatory_status") == "510k_cleared")
        st.metric("510(k) Cleared", cleared)

    st.divider()

    # Filters
    col1, col2, col3 = st.columns(3)
    with col1:
        all_modalities = sorted(set(d.get("modality", "") for d in devices))
        sel_modalities = st.multiselect("Filter by Modality", all_modalities, default=[])
    with col2:
        all_tasks = sorted(set(d.get("ai_task", "") for d in devices))
        sel_tasks = st.multiselect("Filter by AI Task", all_tasks, default=[])
    with col3:
        all_regions = sorted(set(d.get("body_region", "") for d in devices))
        sel_regions = st.multiselect("Filter by Body Region", all_regions, default=[])

    # Apply filters
    filtered = devices
    if sel_modalities:
        filtered = [d for d in filtered if d.get("modality") in sel_modalities]
    if sel_tasks:
        filtered = [d for d in filtered if d.get("ai_task") in sel_tasks]
    if sel_regions:
        filtered = [d for d in filtered if d.get("body_region") in sel_regions]

    st.markdown(f"**Showing {len(filtered)} of {len(devices)} devices**")

    # Distribution charts
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**By Modality**")
        mod_counts = Counter(d.get("modality", "unknown") for d in filtered)
        st.bar_chart(dict(sorted(mod_counts.items(), key=lambda x: -x[1])))
    with col2:
        st.markdown("**By AI Task**")
        task_counts = Counter(d.get("ai_task", "unknown") for d in filtered)
        st.bar_chart(dict(sorted(task_counts.items(), key=lambda x: -x[1])))

    # Device table
    st.subheader("Device Details")
    for d in filtered:
        with st.expander(f"**{d.get('device_name', 'Unknown')}** — {d.get('manufacturer', '')} | {d.get('modality', '').upper()} | {d.get('ai_task', '')}"):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"**Regulatory:** {d.get('regulatory_status', 'N/A')}")
                st.markdown(f"**Clearance Date:** {d.get('clearance_date', 'N/A')}")
                st.markdown(f"**Body Region:** {d.get('body_region', 'N/A')}")
                st.markdown(f"**Architecture:** {d.get('model_architecture', 'N/A')}")
            with col2:
                st.markdown(f"**Intended Use:** {d.get('intended_use', 'N/A')}")
                st.markdown(f"**Performance:** {d.get('performance_summary', 'N/A')}")
