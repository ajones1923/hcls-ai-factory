"""
Patient 360 — Unified Patient Intelligence Dashboard.

A single Streamlit page showing a patient's complete precision medicine profile
across all HCLS AI Factory agents. This is the visualization of the closed-loop
architecture: from DNA to drug candidates, personalized to this patient.

Part of the HCLS AI Factory: Patient DNA → Drug Candidates in <5 hours on DGX Spark.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

from typing import Any

import streamlit as st

# ══════════════════════════════════════════════════════════════════════
# NVIDIA-branded styling (consistent across all agents)
# ══════════════════════════════════════════════════════════════════════

NVIDIA_GREEN = "#76B900"
NVIDIA_DARK = "#1B1B2F"
NVIDIA_LIGHT_BG = "#f8f9fa"
NVIDIA_ACCENT = "#333333"

CUSTOM_CSS = f"""
<style>
    /* Global NVIDIA branding */
    .stApp {{
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
    }}

    /* Header with green accent bar */
    .nvidia-header {{
        background: linear-gradient(135deg, {NVIDIA_DARK} 0%, #2d2d4e 100%);
        border-left: 4px solid {NVIDIA_GREEN};
        padding: 1.5rem 2rem;
        border-radius: 8px;
        margin-bottom: 1.5rem;
        color: white;
    }}
    .nvidia-header h1 {{
        color: white;
        margin: 0;
        font-size: 1.8rem;
    }}
    .nvidia-header p {{
        color: #b0b0b0;
        margin: 0.5rem 0 0 0;
    }}

    /* Agent cards */
    .agent-card {{
        background: white;
        border: 1px solid #e0e0e0;
        border-left: 4px solid {NVIDIA_GREEN};
        border-radius: 8px;
        padding: 1.2rem;
        margin-bottom: 1rem;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }}
    .agent-card h3 {{
        color: {NVIDIA_DARK};
        margin: 0 0 0.5rem 0;
        font-size: 1.1rem;
    }}

    /* Metric cards */
    .metric-card {{
        background: {NVIDIA_LIGHT_BG};
        border-radius: 8px;
        padding: 1rem;
        text-align: center;
        border: 1px solid #e0e0e0;
    }}
    .metric-value {{
        font-size: 2rem;
        font-weight: 700;
        color: {NVIDIA_DARK};
    }}
    .metric-label {{
        font-size: 0.85rem;
        color: #666;
        margin-top: 0.25rem;
    }}

    /* Risk indicators */
    .risk-low {{ color: {NVIDIA_GREEN}; font-weight: 600; }}
    .risk-moderate {{ color: #f0a500; font-weight: 600; }}
    .risk-high {{ color: #e74c3c; font-weight: 600; }}
    .risk-critical {{ color: #c0392b; font-weight: 700; background: #fce4e4; padding: 2px 8px; border-radius: 4px; }}

    /* Evidence chain */
    .evidence-item {{
        background: #f8f9fa;
        border-left: 3px solid {NVIDIA_GREEN};
        padding: 0.75rem 1rem;
        margin: 0.5rem 0;
        border-radius: 0 4px 4px 0;
        font-size: 0.9rem;
    }}
    .evidence-source {{
        color: #888;
        font-size: 0.8rem;
    }}

    /* Timeline */
    .timeline-step {{
        display: flex;
        align-items: flex-start;
        margin-bottom: 1rem;
    }}
    .timeline-dot {{
        width: 12px;
        height: 12px;
        border-radius: 50%;
        background: {NVIDIA_GREEN};
        margin-top: 4px;
        margin-right: 12px;
        flex-shrink: 0;
    }}
    .timeline-dot.inactive {{
        background: #ccc;
    }}

    /* Pipeline flow */
    .pipeline-flow {{
        display: flex;
        justify-content: space-between;
        align-items: center;
        padding: 1rem;
        background: {NVIDIA_DARK};
        border-radius: 8px;
        margin: 1rem 0;
    }}
    .pipeline-stage {{
        text-align: center;
        color: white;
        padding: 0.5rem;
    }}
    .pipeline-stage.active {{
        color: {NVIDIA_GREEN};
    }}
    .pipeline-arrow {{
        color: {NVIDIA_GREEN};
        font-size: 1.5rem;
    }}

    /* Footer */
    .nvidia-footer {{
        text-align: center;
        color: #888;
        font-size: 0.8rem;
        margin-top: 2rem;
        padding: 1rem;
        border-top: 1px solid #e0e0e0;
    }}
</style>
"""

# Agent icons (Unicode for broad compatibility)
AGENT_ICONS = {
    "biomarker": "\U0001F9EC",      # DNA helix
    "oncology": "\U0001F3AF",       # Target
    "cart": "\U0001F9EA",           # Test tube (T-cell)
    "imaging": "\U0001F9E0",        # Brain
    "autoimmune": "\U0001F6E1\uFE0F",  # Shield
    "drug_discovery": "\U0001F48A",  # Pill
    "genomics": "\U0001F52C",       # Microscope
}


def render_header(patient_id: str, patient_name: str = ""):
    """Render the NVIDIA-branded page header."""
    display = patient_name or patient_id
    st.markdown(CUSTOM_CSS, unsafe_allow_html=True)
    st.markdown(f"""
    <div class="nvidia-header">
        <h1>{AGENT_ICONS['biomarker']} Patient 360 — {display}</h1>
        <p>Unified Precision Medicine Intelligence | HCLS AI Factory on NVIDIA DGX Spark</p>
    </div>
    """, unsafe_allow_html=True)


def render_pipeline_status(stages_complete: list[str]):
    """Render the pipeline flow visualization showing which stages are complete."""
    all_stages = [
        ("Genomics", AGENT_ICONS["genomics"]),
        ("RAG/Chat", "\U0001F4AC"),
        ("Biomarker", AGENT_ICONS["biomarker"]),
        ("Oncology", AGENT_ICONS["oncology"]),
        ("Imaging", AGENT_ICONS["imaging"]),
        ("Drug Discovery", AGENT_ICONS["drug_discovery"]),
    ]

    cols = st.columns(len(all_stages) * 2 - 1)
    for i, (stage_name, icon) in enumerate(all_stages):
        col_idx = i * 2
        is_active = stage_name.lower().replace("/", "_") in [s.lower() for s in stages_complete]
        color = NVIDIA_GREEN if is_active else "#666"
        cols[col_idx].markdown(
            f"<div style='text-align:center;'>"
            f"<span style='font-size:1.5rem;'>{icon}</span><br>"
            f"<span style='color:{color};font-size:0.75rem;font-weight:{'700' if is_active else '400'};'>"
            f"{stage_name}</span>"
            f"{'<br><span style=\"color:' + NVIDIA_GREEN + ';font-size:0.7rem;\">✓</span>' if is_active else ''}"
            f"</div>",
            unsafe_allow_html=True,
        )
        if i < len(all_stages) - 1:
            cols[col_idx + 1].markdown(
                f"<div style='text-align:center;padding-top:0.5rem;color:{NVIDIA_GREEN};'>→</div>",
                unsafe_allow_html=True,
            )


def render_genomic_summary(
    variants: list[dict],
    key_genes: list[str],
    ancestry: str = "",
    pgx_profile: dict = None,
):
    """Left column: Patient genomic summary."""
    st.markdown(f"### {AGENT_ICONS['genomics']} Genomic Profile")

    if ancestry:
        st.markdown(f"**Ancestry:** {ancestry}")

    if key_genes:
        st.markdown("**Key Genes:**")
        for gene in key_genes[:10]:
            st.markdown(f"- `{gene}`")

    if variants:
        st.markdown(f"**Pathogenic Variants:** {len(variants)}")
        for v in variants[:5]:
            gene = v.get("gene", "Unknown")
            change = v.get("protein_change", v.get("hgvs", ""))
            sig = v.get("significance", "")
            st.markdown(f"- **{gene}** {change} — {sig}")

    if pgx_profile and pgx_profile.get("gene_results"):
        st.markdown("**Pharmacogenomic Profile:**")
        for gr in pgx_profile["gene_results"]:
            gene = gr.get("gene", "")
            phenotype = gr.get("phenotype", "")
            alleles = gr.get("star_alleles", "")
            risk_class = "risk-high" if "poor" in phenotype.lower() or "ultra" in phenotype.lower() else "risk-low"
            st.markdown(
                f"- **{gene}** ({alleles}): "
                f"<span class='{risk_class}'>{phenotype}</span>",
                unsafe_allow_html=True,
            )


def render_biomarker_panel(
    biological_age: float,
    chronological_age: int,
    age_acceleration: float,
    mortality_risk: str,
    disease_trajectories: list[dict],
    top_drivers: list[dict] = None,
):
    """Center: Biomarker radar chart + disease trajectory."""
    st.markdown(f"### {AGENT_ICONS['biomarker']} Biomarker Intelligence")

    # Key metrics
    col1, col2, col3 = st.columns(3)
    with col1:
        delta = round(age_acceleration, 1)
        st.metric("Biological Age", f"{biological_age:.1f}y", delta=f"{delta:+.1f}y",
                  delta_color="inverse")
    with col2:
        st.metric("Chronological Age", f"{chronological_age}y")
    with col3:
        risk_color = "🟢" if mortality_risk == "LOW" else "🟡" if mortality_risk == "MODERATE" else "🔴"
        st.metric("Mortality Risk", f"{risk_color} {mortality_risk}")

    # Top aging drivers
    if top_drivers:
        st.markdown("**Top Aging Drivers:**")
        for driver in top_drivers[:5]:
            marker = driver.get("biomarker", driver.get("marker", ""))
            contrib = driver.get("contribution_pct", driver.get("contribution", 0))
            direction = driver.get("direction", "")
            bar_width = min(abs(contrib), 100)
            color = "#e74c3c" if direction == "accelerating" else NVIDIA_GREEN
            st.markdown(
                f"<div style='margin:2px 0;'>"
                f"<span style='display:inline-block;width:100px;'>{marker}</span>"
                f"<span style='display:inline-block;width:{bar_width}%;max-width:200px;height:12px;"
                f"background:{color};border-radius:6px;'></span>"
                f" <span style='font-size:0.8rem;color:#666;'>{contrib:.0f}%</span>"
                f"</div>",
                unsafe_allow_html=True,
            )

    # Disease trajectories
    if disease_trajectories:
        st.markdown("**Disease Trajectories:**")
        for traj in disease_trajectories:
            disease = traj.get("disease", traj.get("category", ""))
            risk = traj.get("risk_level", traj.get("risk", ""))
            stage = traj.get("stage", "")
            risk_class = "risk-critical" if risk.lower() in ("high", "critical") else \
                         "risk-moderate" if risk.lower() == "moderate" else "risk-low"
            st.markdown(
                f"<div class='evidence-item'>"
                f"<strong>{disease}</strong>: <span class='{risk_class}'>{risk}</span>"
                f"{' — ' + stage if stage else ''}"
                f"</div>",
                unsafe_allow_html=True,
            )


def render_drug_candidates(
    candidates: list[dict],
    pgx_filtered: list[dict] = None,
):
    """Right column: Drug candidates ranked by PGx compatibility."""
    st.markdown(f"### {AGENT_ICONS['drug_discovery']} Drug Candidates")

    display_candidates = pgx_filtered if pgx_filtered else candidates

    if not display_candidates:
        st.info("No drug candidates available. Run Drug Discovery Pipeline to generate candidates.")
        return

    if pgx_filtered:
        st.success(f"PGx-filtered: {len(pgx_filtered)} candidates ranked by patient compatibility")

    for i, cand in enumerate(display_candidates[:10]):
        name = cand.get("name", cand.get("smiles", f"Candidate {i+1}"))
        score = cand.get("binding_affinity", cand.get("score", 0))
        risk = cand.get("metabolism_risk", "unknown")

        risk_icon = "🟢" if risk == "safe" else "🟡" if risk == "caution" else "🔴" if risk == "contraindicated" else "⚪"

        st.markdown(
            f"<div class='agent-card'>"
            f"<h3>{risk_icon} #{i+1}: {name}</h3>"
            f"<div>Binding: <strong>{score}</strong> kcal/mol</div>"
            f"<div>PGx Status: <strong>{risk}</strong></div>"
            f"{'<div style=\"color:#e74c3c;\">' + cand.get('recommendation', '') + '</div>' if risk == 'contraindicated' else ''}"
            f"</div>",
            unsafe_allow_html=True,
        )


def render_evidence_chain(evidence: list[dict]):
    """Bottom: Evidence provenance chain."""
    st.markdown("### 📋 Evidence Provenance Chain")
    st.caption("Every recommendation traced to its published source")

    if not evidence:
        st.info("Evidence chain will populate as agents process data.")
        return

    for item in evidence:
        agent = item.get("agent", "")
        finding = item.get("finding", "")
        source = item.get("source", "")
        pmid = item.get("pmid", "")
        confidence = item.get("confidence", "")

        source_link = f"[PMID:{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid})" if pmid else source

        st.markdown(
            f"<div class='evidence-item'>"
            f"<strong>{agent}</strong>: {finding}"
            f"<div class='evidence-source'>Source: {source_link}"
            f"{f' | Confidence: {confidence}' if confidence else ''}</div>"
            f"</div>",
            unsafe_allow_html=True,
        )


def render_concordance_matrix(concordance: dict[str, float]):
    """Show concordance scores between agent findings."""
    if not concordance:
        return

    st.markdown("### 🔗 Cross-Agent Concordance")
    st.caption("When multiple agents agree on a finding, confidence increases")

    for pair, score in concordance.items():
        color = NVIDIA_GREEN if score > 0.8 else "#f0a500" if score > 0.5 else "#e74c3c"
        bar_width = int(score * 100)
        st.markdown(
            f"<div style='margin:4px 0;'>"
            f"<span style='display:inline-block;width:200px;font-size:0.85rem;'>{pair}</span>"
            f"<span style='display:inline-block;width:{bar_width}%;max-width:300px;height:16px;"
            f"background:{color};border-radius:8px;'></span>"
            f" <span style='font-size:0.85rem;font-weight:600;'>{score:.0%}</span>"
            f"</div>",
            unsafe_allow_html=True,
        )


def render_footer():
    """Render the NVIDIA-branded footer."""
    st.markdown(
        "<div class='nvidia-footer'>"
        "HCLS AI Factory | Patient DNA → Drug Candidates in &lt;5 hours | "
        "NVIDIA DGX Spark ($4,699) | Apache 2.0"
        "</div>",
        unsafe_allow_html=True,
    )


# ══════════════════════════════════════════════════════════════════════
# Main page (can be used as a standalone Streamlit page or tab)
# ══════════════════════════════════════════════════════════════════════

def render_patient_360(patient_data: dict[str, Any] | None = None):
    """Main entry point for the Patient 360 dashboard.

    Parameters
    ----------
    patient_data : dict
        Complete patient data from PatientCase or session state.
        If None, uses st.session_state.
    """
    data = patient_data or {}

    # Try loading from session state if available
    if not data and hasattr(st, "session_state"):
        data = {
            "patient_id": getattr(st.session_state, "patient_id", "DEMO-001"),
            "biological_age": getattr(st.session_state, "biological_age", None),
            "age_acceleration": getattr(st.session_state, "age_acceleration", None),
            "disease_trajectories": getattr(st.session_state, "disease_trajectories", []),
            "pgx_results": getattr(st.session_state, "pgx_results", {}),
            "drug_candidates": getattr(st.session_state, "drug_candidates", []),
        }

    patient_id = data.get("patient_id", "DEMO-001")

    # Header
    render_header(patient_id)

    # Pipeline status
    stages = data.get("stages_complete", ["Genomics", "RAG/Chat", "Biomarker"])
    render_pipeline_status(stages)

    st.divider()

    # Three-column layout
    left, center, right = st.columns([1, 1.5, 1])

    with left:
        render_genomic_summary(
            variants=data.get("variants", []),
            key_genes=data.get("key_genes", []),
            ancestry=data.get("ancestry", ""),
            pgx_profile=data.get("pgx_results", {}),
        )

    with center:
        render_biomarker_panel(
            biological_age=data.get("biological_age", data.get("age", 0)),
            chronological_age=data.get("age", 0),
            age_acceleration=data.get("age_acceleration", 0),
            mortality_risk=data.get("mortality_risk", "UNKNOWN"),
            disease_trajectories=data.get("disease_trajectories", []),
            top_drivers=data.get("top_aging_drivers", []),
        )

    with right:
        render_drug_candidates(
            candidates=data.get("drug_candidates", []),
            pgx_filtered=data.get("pgx_filtered_candidates", []),
        )

    st.divider()

    # Evidence chain + concordance
    ev_col, conc_col = st.columns([2, 1])

    with ev_col:
        render_evidence_chain(data.get("evidence_chain", []))

    with conc_col:
        render_concordance_matrix(data.get("concordance_scores", {}))

    # Footer
    render_footer()


def render_demo_patient_360():
    """Render with demo data for showcasing the platform."""
    demo_data = {
        "patient_id": "DEMO-VCP-001",
        "age": 45,
        "sex": "M",
        "ancestry": "European",
        "stages_complete": ["Genomics", "RAG/Chat", "Biomarker", "Oncology", "Drug Discovery"],

        # Genomic profile
        "variants": [
            {"gene": "VCP", "protein_change": "p.R155H", "significance": "Pathogenic",
             "hgvs": "NM_007126.5:c.464G>A"},
            {"gene": "MTHFR", "protein_change": "p.C677T", "significance": "Risk Factor"},
        ],
        "key_genes": ["VCP", "MTHFR", "CYP2D6", "CYP2C19"],

        # PGx
        "pgx_results": {
            "gene_results": [
                {"gene": "CYP2D6", "star_alleles": "*1/*2", "phenotype": "Normal Metabolizer"},
                {"gene": "CYP2C19", "star_alleles": "*1/*1", "phenotype": "Normal Metabolizer"},
                {"gene": "SLCO1B1", "star_alleles": "*1/*5", "phenotype": "Intermediate Function"},
            ],
        },

        # Biomarker
        "biological_age": 48.2,
        "age_acceleration": 3.2,
        "mortality_risk": "MODERATE",
        "top_aging_drivers": [
            {"biomarker": "CRP", "contribution_pct": 32, "direction": "accelerating"},
            {"biomarker": "Glucose", "contribution_pct": 18, "direction": "accelerating"},
            {"biomarker": "Albumin", "contribution_pct": 15, "direction": "protective"},
            {"biomarker": "Creatinine", "contribution_pct": 12, "direction": "accelerating"},
        ],

        # Disease trajectories
        "disease_trajectories": [
            {"disease": "Type 2 Diabetes", "risk_level": "Moderate", "stage": "Pre-diabetic"},
            {"disease": "Cardiovascular", "risk_level": "Low", "stage": "Normal"},
            {"disease": "Neurodegeneration (VCP)", "risk_level": "High", "stage": "Pre-symptomatic"},
        ],

        # Drug candidates (PGx-filtered)
        "drug_candidates": [
            {"name": "CB-5083-A7", "binding_affinity": -9.2, "metabolism_risk": "safe", "recommendation": "PGx-compatible"},
            {"name": "CB-5083-B3", "binding_affinity": -8.8, "metabolism_risk": "safe", "recommendation": "PGx-compatible"},
            {"name": "CB-5083-C1", "binding_affinity": -8.5, "metabolism_risk": "caution", "recommendation": "SLCO1B1 intermediate: monitor levels"},
            {"name": "VCP-INH-12", "binding_affinity": -7.9, "metabolism_risk": "safe", "recommendation": "PGx-compatible"},
        ],
        "pgx_filtered_candidates": [
            {"name": "CB-5083-A7", "binding_affinity": -9.2, "metabolism_risk": "safe", "recommendation": "PGx-compatible: no metabolic concerns"},
            {"name": "CB-5083-B3", "binding_affinity": -8.8, "metabolism_risk": "safe", "recommendation": "PGx-compatible: no metabolic concerns"},
            {"name": "VCP-INH-12", "binding_affinity": -7.9, "metabolism_risk": "safe", "recommendation": "PGx-compatible: no metabolic concerns"},
            {"name": "CB-5083-C1", "binding_affinity": -8.5, "metabolism_risk": "caution", "recommendation": "SLCO1B1 intermediate function: dose adjustment may be needed"},
        ],

        # Evidence chain
        "evidence_chain": [
            {"agent": "Genomics Pipeline", "finding": "VCP p.R155H pathogenic variant identified",
             "source": "ClinVar", "pmid": "28622507", "confidence": "Pathogenic"},
            {"agent": "Biomarker Agent", "finding": "Biological age accelerated by 3.2 years, driven by CRP elevation",
             "source": "PhenoAge (Levine et al. 2018)", "pmid": "29676998", "confidence": "95% CI: [44.7, 51.7]"},
            {"agent": "Biomarker Agent", "finding": "CYP2D6 *1/*2: Normal Metabolizer — no drug restrictions",
             "source": "CPIC Guideline", "pmid": "31562822"},
            {"agent": "Drug Discovery", "finding": "CB-5083-A7 binds VCP D2 domain at -9.2 kcal/mol",
             "source": "DiffDock (PDB:5FTK)", "confidence": "Top 1/50 candidates"},
            {"agent": "Biomarker Agent", "finding": "Pre-diabetic trajectory detected: HbA1c 5.9%, fasting glucose 108",
             "source": "ADA Criteria 2024", "confidence": "Moderate risk"},
        ],

        # Concordance
        "concordance_scores": {
            "Biomarker ↔ Imaging": 0.85,
            "Biomarker ↔ Oncology": 0.72,
            "Genomics ↔ Biomarker": 0.94,
            "PGx ↔ Drug Discovery": 0.91,
        },
    }

    render_patient_360(demo_data)

    # ── Cohort Comparison ──────────────────────────────────────
    st.markdown("---")
    st.markdown("### 📊 Population Cohort Comparison")
    st.caption("How this patient compares to reference populations and VCP mutation carriers")

    col_a, col_b = st.columns(2)

    with col_a:
        st.markdown("**Biological Age vs. Population**")
        try:
            import plotly.graph_objects as go

            fig = go.Figure()
            # Population distribution (simulated normal curve)
            import numpy as np
            x_pop = np.linspace(40, 70, 100)
            y_pop = np.exp(-0.5 * ((x_pop - 52) / 5) ** 2)

            fig.add_trace(go.Scatter(
                x=x_pop, y=y_pop, fill='tozeroy', name='Population (age 52)',
                fillcolor='rgba(118,185,0,0.2)', line=dict(color='#76B900'),
            ))
            fig.add_vline(x=58.2, line_dash="dash", line_color="#ef4444",
                          annotation_text="This Patient: 58.2y", annotation_font_color="#ef4444")
            fig.add_vline(x=52, line_dash="dot", line_color="#76B900",
                          annotation_text="Pop. Average: 52y", annotation_font_color="#76B900")

            fig.update_layout(
                height=300, margin=dict(l=20, r=20, t=40, b=20),
                plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
                font=dict(color='#333'), showlegend=False,
                title=dict(text="Biological Age Distribution", font=dict(size=13)),
                xaxis_title="Biological Age (years)",
            )
            st.plotly_chart(fig, use_container_width=True)
        except ImportError:
            st.info("Install plotly and numpy for cohort visualization")

    with col_b:
        st.markdown("**Disease Risk Radar**")
        try:
            import plotly.graph_objects as go

            categories = ['Diabetes', 'Cardiovascular', 'Liver', 'Thyroid', 'Iron', 'Nutritional']
            patient_risk = [65, 45, 20, 15, 10, 35]  # demo patient percentiles
            pop_avg = [25, 30, 15, 10, 8, 20]  # population average

            fig = go.Figure()
            fig.add_trace(go.Scatterpolar(
                r=patient_risk, theta=categories, fill='toself', name='This Patient',
                fillcolor='rgba(239,68,68,0.2)', line=dict(color='#ef4444'),
            ))
            fig.add_trace(go.Scatterpolar(
                r=pop_avg, theta=categories, fill='toself', name='Population Avg',
                fillcolor='rgba(118,185,0,0.1)', line=dict(color='#76B900'),
            ))

            fig.update_layout(
                polar=dict(
                    radialaxis=dict(visible=True, range=[0, 100]),
                    bgcolor='rgba(0,0,0,0)',
                ),
                height=300, margin=dict(l=40, r=40, t=40, b=20),
                plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
                font=dict(color='#333'), showlegend=True,
                title=dict(text="Disease Risk Profile", font=dict(size=13)),
            )
            st.plotly_chart(fig, use_container_width=True)
        except ImportError:
            st.info("Install plotly for radar chart")

    # VCP Cohort stats
    st.markdown("**VCP p.R155H Carrier Cohort (N=47)**")
    cohort_cols = st.columns(4)
    with cohort_cols[0]:
        st.metric("Median Onset Age", "56y", delta="-4y vs. this patient", delta_color="inverse")
    with cohort_cols[1]:
        st.metric("5-Year FTD Risk", "34%", delta="+12% vs. population")
    with cohort_cols[2]:
        st.metric("IBM Prevalence", "62%", help="Inclusion body myopathy in VCP carriers")
    with cohort_cols[3]:
        st.metric("Paget Disease", "43%", help="Paget disease of bone in VCP carriers")


# Run as standalone page
if __name__ == "__main__":
    st.set_page_config(
        page_title="Patient 360 — HCLS AI Factory",
        page_icon=AGENT_ICONS["biomarker"],
        layout="wide",
    )
    render_demo_patient_360()
