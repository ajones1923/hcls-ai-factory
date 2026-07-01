"""
Surface (b) — In-Visit Dashboard (PRD §7.2; port 8562). Polished four-quadrant clinical UI.
Run:  streamlit run app/dashboard_app.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # engine root on path

import streamlit as st

from app import _ui
from app._engine import featured, get_engine
from src.utils.model_router import get_router

st.set_page_config(page_title="TSC In-Visit Dashboard", layout="wide")
_ui.inject_css()

orch, manifest = get_engine()
fmap = featured()
_ui.watermark("live" if get_router().online else "offline (deterministic stubs)")

labels = {pid: f"{tag} · {pid}" for tag, pid in fmap.items()}
labels.update({p["patient_id"]: p["patient_id"] for p in manifest["patients"] if p["patient_id"] not in labels})
pid = st.selectbox("Patient", options=list(labels), format_func=lambda p: labels[p],
                   index=list(labels).index(fmap.get("B", list(labels)[0])))

proj = orch.store.projection(pid)
st.title(f"In-Visit Dashboard — {pid}")

c1, c2 = st.columns(2)
with c1:
    _ui.variant_card((proj.get("variant_interp") or {}).get("primary"))
    st.markdown("<div class='card'><h4>Trajectory forecast</h4></div>", unsafe_allow_html=True)
    traj = proj.get("trajectory") or {}
    quants = traj.get("quantities") or traj.get("lesions", [])
    if quants:
        for tab, q in zip(st.tabs([f"{l['lesion']}" for l in quants]), quants):
            with tab:
                st.altair_chart(_ui.trajectory_chart(q), width="stretch")
                if q.get("crosses_in_12_18mo_window"):
                    note = f"forecast crosses threshold in ~{q['months_to_threshold']} months"
                elif q.get("at_or_above_threshold"):
                    note = "at or beyond the intervention threshold now"
                else:
                    note = "stable relative to threshold"
                grade = q.get("crossing_grade", "")
                p18 = (q.get("crossing_probability") or {}).get("m18")
                sr = q.get("surveillance_recommendation") or {}
                unit = q.get("unit", "cm")
                st.caption(f"{q['lesion']} @ {q['location']} · threshold {q['threshold_cm']} {unit} · {note}")
                if grade:
                    st.markdown("**Crossing risk:** " + _ui.chips([
                        f"{grade} (P₁₈={p18})",
                        f"surveillance {sr.get('itsc_floor_months')}→{sr.get('recommended_interval_months')}mo "
                        f"({sr.get('modality','')})",
                    ]), unsafe_allow_html=True)
    else:
        st.caption("No longitudinal series for this patient.")

with c2:
    st.markdown("<div class='card'><h4>Longitudinal phenotype (HPO)</h4></div>", unsafe_allow_html=True)
    _ui.hpo_timeline(proj.get("hpo_profile"))
    st.markdown("<div class='card'><h4>TAND surveillance · briefing material, never an alert</h4></div>", unsafe_allow_html=True)
    _ui.tand_block(proj.get("tand_briefing"))

th = proj.get("therapeutics")
if th:
    st.markdown("<div class='card'><h4>Therapeutic options brief · decision support, not a recommendation</h4></div>", unsafe_allow_html=True)
    cols = st.columns(2)
    for i, (k, v) in enumerate(th.get("sections", {}).items()):
        with cols[i % 2]:
            with st.expander(k.replace("_", " ").title(), expanded=(k == "trial_matching")):
                st.write(v)
    tms = th.get("trial_matches", [])
    if tms:
        st.markdown("**Trial matches:** " + _ui.chips(f"{m['nct']}: {m['match']}" for m in tms), unsafe_allow_html=True)
    if th.get("sources"):
        st.caption("Cited sources: " + " · ".join(s["source_uri"] for s in th["sources"]))

with st.expander("🔍  Audit trail / provenance — every output is traceable to model, tier, latency"):
    recs = [{"event": e.event_type.value, "records": (e.provenance or {}).get("records", [])}
            for e in orch.store.events_for(pid) if (e.provenance or {}).get("records")]
    _ui.provenance_table(recs)

st.caption(f"SYNTHETIC · {manifest['n_patients']} patients · staleness: "
           f"{ {k: v['status'] for k, v in proj['staleness'].items()} }")
