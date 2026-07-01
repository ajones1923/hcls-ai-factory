"""Shared UI helpers for the clinician surfaces — clinical-grade styling, badges,
the trajectory chart (altair), HPO timeline, TAND highlights, and the audit trail."""
from __future__ import annotations

import altair as alt
import pandas as pd
import streamlit as st

NAVY = "#0f1830"
ACCENT = "#7fb0ff"
_CLASS_COLOR = {
    "Pathogenic": "#e5484d", "Likely Pathogenic": "#e8833a",
    "Variant of Uncertain Significance": "#8b8f9a",
    "Likely Benign": "#3aa675", "Benign": "#2f9e64",
}


def inject_css() -> None:
    st.markdown("""<style>
      .block-container{padding-top:1.4rem;max-width:1200px}
      .card{background:#16213f;border:1px solid #26345c;border-radius:12px;padding:16px 18px;margin-bottom:10px}
      .card h4{margin:0 0 8px;color:#cfe0ff;font-size:14px;text-transform:uppercase;letter-spacing:.04em}
      .badge{display:inline-block;padding:3px 11px;border-radius:20px;color:#fff;font-weight:600;font-size:13px}
      .chip{display:inline-block;background:#22325c;color:#bcd2ff;border:1px solid #34467a;
            border-radius:6px;padding:2px 9px;margin:2px 4px 2px 0;font-size:12px}
      .mut{font-family:'JetBrains Mono',monospace;color:#9fc1ff}
      .hl{background:#3a2d12;border-left:3px solid #e8b54a;padding:8px 12px;border-radius:6px;margin:6px 0;color:#f0e6cf}
      .mk{display:inline-block;background:#5a4413;color:#ffd98a;border-radius:4px;padding:1px 7px;margin-right:5px;font-size:11px}
      .small{color:#9db4dc;font-size:12px}
    </style>""", unsafe_allow_html=True)


def watermark(reasoning: str = "offline") -> None:
    st.markdown(
        f"<div style='background:#3b2330;color:#ffb3c1;border:1px solid #6b3a4d;padding:6px 12px;"
        f"border-radius:8px;font-size:12px;margin-bottom:10px'>SYNTHETIC demonstration data · "
        f"reasoning: {reasoning} · decision support, clinician review required · not FDA-cleared</div>",
        unsafe_allow_html=True)


def class_badge(cls: str | None) -> str:
    c = _CLASS_COLOR.get(cls or "", "#8b8f9a")
    return f"<span class='badge' style='background:{c}'>{cls or '—'}</span>"


def chips(items) -> str:
    return "".join(f"<span class='chip'>{i}</span>" for i in items)


def variant_card(vi: dict) -> None:
    if not vi:
        st.info("No reportable variant on the available sample.")
        return
    crit = [f"{c['code']} ({c['bucket']})" for c in vi.get("acmg_criteria", [])]
    mosaic = " <span class='chip'>MOSAIC recovered</span>" if vi.get("recovered") else ""
    st.markdown(
        f"<div class='card'><h4>Variant interpretation</h4>"
        f"<div class='mut'>{vi.get('gene')} {vi.get('hgvsc') or ''} {vi.get('hgvsp') or ''}</div>"
        f"<div style='margin:8px 0'>{class_badge(vi.get('acmg_classification'))}"
        f" <span class='small'>VAF {vi.get('vaf')} · {vi.get('acmg_rule','')}</span>{mosaic}</div>"
        f"<div>{chips(crit)}</div>"
        f"<div class='small' style='margin-top:8px'>ddPCR: {vi.get('ddpcr_recommended')} · "
        f"review: {vi.get('review_status')}</div></div>", unsafe_allow_html=True)
    nar = vi.get("synthesis_narrative")
    if nar and nar.get("summary"):
        with st.expander("ACMG per-criterion reasoning (live model narrative)"):
            for pc in nar.get("per_criterion", []):
                st.markdown(f"**{pc.get('code')}** — {pc.get('reasoning')}")
            st.markdown(f"_{nar['summary']}_")


def trajectory_chart(lesion: dict):
    obs = pd.DataFrame({"month": lesion.get("observed_months", []),
                        "cm": lesion.get("observed_cm", [])})
    rows = []
    for h in (6, 12, 18):
        d = (lesion.get("forecast") or {}).get(f"m{h}")
        if d:
            rows.append({"month": h, "mean": d["mean_cm"], "lo50": d["pi50"][0], "hi50": d["pi50"][1],
                         "lo90": d["pi90"][0], "hi90": d["pi90"][1]})
    fc = pd.DataFrame(rows)
    thr = lesion.get("threshold_cm")
    layers = []
    if not fc.empty:
        layers += [
            alt.Chart(fc).mark_area(opacity=0.13, color=ACCENT).encode(x=alt.X("month:Q", title="months (0 = now)"), y=alt.Y("lo90:Q", title=f"{lesion['lesion']} size (cm)"), y2="hi90:Q"),
            alt.Chart(fc).mark_area(opacity=0.28, color=ACCENT).encode(x="month:Q", y="lo50:Q", y2="hi50:Q"),
            alt.Chart(fc).mark_line(point=True, color=ACCENT).encode(x="month:Q", y="mean:Q"),
        ]
    if not obs.empty:
        layers += [
            alt.Chart(obs).mark_line(strokeDash=[4, 2], color="#ffd166").encode(x="month:Q", y="cm:Q"),
            alt.Chart(obs).mark_point(filled=True, size=95, color="#ffd166").encode(x="month:Q", y="cm:Q"),
        ]
    if thr is not None:
        layers.append(alt.Chart(pd.DataFrame({"y": [thr]})).mark_rule(color="#ff6b6b", strokeDash=[6, 3]).encode(y="y:Q"))
    return alt.layer(*layers).properties(height=240).configure_view(strokeOpacity=0)


_VAL_MARK = {"verified": "✓ verified", "relabeled": "✓ relabeled", "remapped": "✓ remapped",
             "recovered": "✓ recovered", "unknown": "⚠ unknown", "unverified": "– unverified"}


def hpo_timeline(profile: dict) -> None:
    profile = profile or {}
    terms = profile.get("hpo_terms", [])
    if terms:
        df = pd.DataFrame([{"HPO": t["hpo_id"], "Phenotype": t["label"],
                            "Onset": t.get("onset") or "—", "Source": t.get("source", ""),
                            "Ontology": _VAL_MARK.get(t.get("validation", "unverified"), t.get("validation"))}
                           for t in terms])
        st.dataframe(df, hide_index=True, width="stretch")
    rel = profile.get("ontology_release")
    if rel and rel != "unavailable":
        dropped = profile.get("n_dropped_unverified", 0)
        extra = f" · {dropped} model-emitted code(s) rejected as non-ontology" if dropped else ""
        st.caption(f"🔗 Every term grounded against the {rel}{extra}")
    gaps = profile.get("surveillance_gaps", [])
    overdue = [g for g in gaps if g.get("status") == "overdue"]
    if overdue:
        st.markdown("**Surveillance gaps:** " + chips(f"{g['type']} (overdue)" for g in overdue), unsafe_allow_html=True)
    for d in profile.get("discordances", []):
        st.markdown(
            f"<div class='hl'><span class='mk'>discordance</span> {d.get('label')} — present in "
            f"{', '.join(d.get('present_sources', []))}; negated in {', '.join(d.get('absent_sources', []))}"
            f"<div class='small'>{d.get('note','')}</div></div>", unsafe_allow_html=True)


def tand_block(tand: dict) -> None:
    if not tand:
        return
    flagged = tand.get("flagged_clusters", [])
    st.markdown(f"**TAND clusters flagged:** {chips(flagged) if flagged else '<span class=small>none above threshold</span>'}",
                unsafe_allow_html=True)
    for h in tand.get("discourse_highlights", []):
        marks = "".join(f"<span class='mk'>{m}</span>" for m in h.get("markers", []))
        st.markdown(f"<div class='hl'>{marks}<br>{h.get('passage')}<div class='small'>— {h.get('source')} · cluster: {h.get('cluster')}</div></div>",
                    unsafe_allow_html=True)
    st.caption(tand.get("recommendation", ""))


def provenance_table(records: list[dict]) -> None:
    rows = []
    for ev in records:
        for r in ev.get("records", []):
            rows.append({"event": ev["event"], "step": r.get("step"), "tier": r.get("tier"),
                         "model": r.get("model_id"), "latency_ms": r.get("latency_ms"),
                         "prompt_v": r.get("prompt_template_version")})
    if rows:
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    else:
        st.caption("No model-backed provenance records (deterministic-only patient).")
