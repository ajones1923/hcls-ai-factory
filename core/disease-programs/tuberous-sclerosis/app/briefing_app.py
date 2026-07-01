"""
Surface (a) — Pre-Visit Briefing (PRD §7.1; port 8561). One-screen, mobile-readable.
Run:  streamlit run app/briefing_app.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # engine root on path

import streamlit as st

from app import _ui
from app._engine import featured, get_engine
from src.utils.model_router import get_router

st.set_page_config(page_title="TSC Pre-Visit Briefing", layout="centered")
_ui.inject_css()

orch, manifest = get_engine()
fmap = featured()
_ui.watermark("live" if get_router().online else "offline (deterministic stubs)")

labels = {pid: f"{tag} · {pid}" for tag, pid in fmap.items()}
labels.update({p["patient_id"]: p["patient_id"] for p in manifest["patients"] if p["patient_id"] not in labels})
pid = st.selectbox("Patient", options=list(labels), format_func=lambda p: labels[p],
                   index=list(labels).index(fmap.get("A", list(labels)[0])))

s = orch.assemble_surface(pid, "briefing")
h = s["header"]
st.markdown(
    f"<div class='card'><h4>Pre-visit briefing · {pid}</h4>"
    f"<div class='mut'>{h.get('genotype') or ''} {h.get('variant') or ''}</div>"
    f"<div style='margin-top:8px'>{_ui.class_badge(h.get('classification'))}</div></div>",
    unsafe_allow_html=True)

st.markdown("<div class='card'><h4>Action items</h4></div>", unsafe_allow_html=True)
if s["action_items"]:
    for a in s["action_items"]:
        st.markdown(f"<div class='hl'>{a}</div>", unsafe_allow_html=True)
else:
    st.caption("No action items this visit.")

ok = sum(1 for v in s["staleness"].values() if v["status"] == "ok")
st.markdown(f"<div class='small' style='margin-top:10px'>{ok}/{len(s['staleness'])} sections current · "
            f"open the dashboard for detail and the full audit trail.</div>", unsafe_allow_html=True)
st.caption("SYNTHETIC demonstration data · decision support · clinician review required")
