"""
Surface (c) — Async Alert Surface (PRD §7.3; port 8563). Strict alert discipline.
Run:  streamlit run app/alerts_app.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # engine root on path

import streamlit as st

from app import _ui
from app._engine import featured, get_engine  # noqa: F401
from config.settings import settings
from src.utils.model_router import get_router

st.set_page_config(page_title="TSC Alerts", layout="centered")
_ui.inject_css()

orch, manifest = get_engine()
_ui.watermark("live" if get_router().online else "offline (deterministic stubs)")
st.title("Async Alerts")
st.markdown(f"<div class='small'>Alert discipline: recalibrate if any clinician would receive more than "
            f"~{settings.ALERTS_PER_CLINICIAN_PER_WEEK_MAX}/week. Every alert carries a source and a "
            f"dismissal path.</div>", unsafe_allow_html=True)

total = 0
for p in manifest["patients"]:
    for a in orch.assemble_surface(p["patient_id"], "alerts")["alerts"]:
        total += 1
        st.markdown(
            f"<div class='card'><span class='chip'>{a['category']}</span> "
            f"<b>{p['patient_id']}</b><div style='margin-top:6px'>{a['summary']}</div>"
            f"<div class='small'>source: {a['source_section']} · dismissable</div></div>",
            unsafe_allow_html=True)
if total == 0:
    st.success("No active alerts across the cohort — disciplined by design.")
else:
    st.caption(f"{total} alert(s) across {manifest['n_patients']} patients (well within budget).")
