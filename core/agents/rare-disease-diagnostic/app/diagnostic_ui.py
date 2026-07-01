"""Rare Disease Diagnostic Agent -- 5-Tab Streamlit UI.

NVIDIA dark-themed rare disease diagnostic support interface with
patient intake, differential diagnosis, variant review, therapeutic
options, and multi-format report generation.

Usage:
    streamlit run app/diagnostic_ui.py --server.port 8544

Author: Adam Jones
Date: March 2026
"""

import json
import os
from datetime import datetime
from typing import Optional

import requests
import streamlit as st

# =====================================================================
# Configuration
# =====================================================================

API_BASE = os.environ.get("RD_API_BASE", "http://localhost:8134")

NVIDIA_THEME = {
    "bg_primary": "#1a1a2e",
    "bg_secondary": "#16213e",
    "bg_card": "#0f3460",
    "text_primary": "#e0e0e0",
    "text_secondary": "#a0a0b0",
    "accent": "#76b900",
    "accent_hover": "#8ed100",
    "danger": "#e74c3c",
    "warning": "#f39c12",
    "info": "#3498db",
    "success": "#76b900",
}


# =====================================================================
# Page Config & Custom CSS
# =====================================================================

st.set_page_config(
    page_title="Rare Disease Diagnostic Agent",
    page_icon="\U0001F9EC",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(f"""
<style>
    /* Main background */
    .stApp {{
        background-color: {NVIDIA_THEME['bg_primary']};
        color: {NVIDIA_THEME['text_primary']};
    }}

    /* Sidebar */
    section[data-testid="stSidebar"] {{
        background-color: {NVIDIA_THEME['bg_secondary']};
    }}
    section[data-testid="stSidebar"] .stMarkdown {{
        color: {NVIDIA_THEME['text_primary']};
    }}

    /* Cards */
    div[data-testid="stMetric"] {{
        background-color: {NVIDIA_THEME['bg_card']};
        border: 1px solid {NVIDIA_THEME['accent']};
        border-radius: 8px;
        padding: 12px;
    }}
    div[data-testid="stMetric"] label {{
        color: {NVIDIA_THEME['text_secondary']};
    }}
    div[data-testid="stMetric"] [data-testid="stMetricValue"] {{
        color: {NVIDIA_THEME['accent']};
    }}

    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {{
        background-color: {NVIDIA_THEME['bg_secondary']};
        border-radius: 8px;
        padding: 4px;
    }}
    .stTabs [data-baseweb="tab"] {{
        color: {NVIDIA_THEME['text_secondary']};
    }}
    .stTabs [aria-selected="true"] {{
        color: {NVIDIA_THEME['accent']};
        border-bottom-color: {NVIDIA_THEME['accent']};
    }}

    /* Buttons */
    .stButton > button {{
        background-color: {NVIDIA_THEME['accent']};
        color: #000000;
        border: none;
        border-radius: 6px;
        font-weight: 600;
    }}
    .stButton > button:hover {{
        background-color: {NVIDIA_THEME['accent_hover']};
        color: #000000;
    }}

    /* Expanders */
    details {{
        background-color: {NVIDIA_THEME['bg_card']};
        border: 1px solid {NVIDIA_THEME['accent']}40;
        border-radius: 6px;
    }}

    /* Text inputs */
    .stTextInput > div > div > input,
    .stTextArea > div > div > textarea {{
        background-color: {NVIDIA_THEME['bg_secondary']};
        color: {NVIDIA_THEME['text_primary']};
        border: 1px solid {NVIDIA_THEME['accent']}60;
    }}

    /* Select boxes */
    .stSelectbox > div > div {{
        background-color: {NVIDIA_THEME['bg_secondary']};
        color: {NVIDIA_THEME['text_primary']};
    }}

    /* Status indicators */
    .status-healthy {{ color: {NVIDIA_THEME['success']}; font-weight: bold; }}
    .status-degraded {{ color: {NVIDIA_THEME['warning']}; font-weight: bold; }}
    .status-error {{ color: {NVIDIA_THEME['danger']}; font-weight: bold; }}

    /* Agent header */
    .agent-header {{
        background: linear-gradient(135deg, {NVIDIA_THEME['bg_card']}, {NVIDIA_THEME['bg_secondary']});
        border-left: 4px solid {NVIDIA_THEME['accent']};
        padding: 16px 20px;
        border-radius: 0 8px 8px 0;
        margin-bottom: 20px;
    }}
</style>
""", unsafe_allow_html=True)

st.warning(
    "**Clinical Decision Support Tool** — This system provides evidence-based guidance "
    "for research and clinical decision support only. All recommendations must be verified "
    "by a qualified healthcare professional. Not FDA-cleared. Not a substitute for professional "
    "clinical judgment."
)


# =====================================================================
# API Helpers
# =====================================================================

def api_get(path: str, timeout: int = 15) -> Optional[dict]:
    """GET request to rare disease API with error handling."""
    try:
        resp = requests.get(f"{API_BASE}{path}", timeout=timeout)
        resp.raise_for_status()
        return resp.json()
    except requests.exceptions.ConnectionError:
        st.error(f"Cannot connect to API at {API_BASE}. Is the server running?")
        return None
    except requests.exceptions.Timeout:
        st.error(f"API request timed out: {path}")
        return None
    except Exception as exc:
        st.error(f"API error: {exc}")
        return None


def api_post(path: str, data: dict, timeout: int = 60) -> Optional[dict]:
    """POST request to rare disease API with error handling."""
    try:
        resp = requests.post(
            f"{API_BASE}{path}",
            json=data,
            timeout=timeout,
            headers={"Content-Type": "application/json"},
        )
        resp.raise_for_status()
        return resp.json()
    except requests.exceptions.ConnectionError:
        st.error(f"Cannot connect to API at {API_BASE}. Is the server running?")
        return None
    except requests.exceptions.Timeout:
        st.error(f"API request timed out: {path}")
        return None
    except requests.exceptions.HTTPError as exc:
        try:
            detail = exc.response.json().get("detail", str(exc))
        except Exception:
            detail = str(exc)
        st.error(f"API error ({exc.response.status_code}): {detail}")
        return None
    except Exception as exc:
        st.error(f"API error: {exc}")
        return None


# =====================================================================
# Sidebar
# =====================================================================

with st.sidebar:
    st.markdown(f"""
    <div class="agent-header">
        <h2 style="color: {NVIDIA_THEME['accent']}; margin: 0;">Rare Disease Diagnostics</h2>
        <p style="color: {NVIDIA_THEME['text_secondary']}; margin: 4px 0 0 0; font-size: 0.85em;">
            HCLS AI Factory Agent
        </p>
    </div>
    """, unsafe_allow_html=True)

    # Health status
    health = api_get("/health")
    if health:
        status = health.get("status", "unknown")
        status_class = "status-healthy" if status == "healthy" else "status-degraded"
        st.markdown(f'<p class="{status_class}">Status: {status.upper()}</p>', unsafe_allow_html=True)

        components = health.get("components", {})
        for comp, state in components.items():
            icon = "+" if state in ("connected", "ready") else "-"
            st.text(f"  {icon} {comp}: {state}")

        st.markdown("---")
        st.metric("Collections", health.get("collections", 0))
        st.metric("Vectors", f"{health.get('total_vectors', 0):,}")
        st.metric("Workflows", health.get("workflows", 0))
    else:
        st.warning("API unavailable")

    st.markdown("---")
    st.caption(f"API: {API_BASE}")
    st.caption(f"v1.0.0 | {datetime.now().strftime('%Y-%m-%d')}")


# =====================================================================
# Main Content - Tabs
# =====================================================================

tab_intake, tab_diagnosis, tab_variants, tab_therapy, tab_reports = st.tabs([
    "Patient Intake",
    "Diagnostic Dashboard",
    "Variant Review",
    "Therapeutic Options",
    "Report Generator",
])


# =====================================================================
# Tab 1: Patient Intake
# =====================================================================

with tab_intake:
    st.header("Patient Intake")
    st.write("Enter patient phenotype, genotype, and clinical information for rare disease analysis.")

    col_left, col_right = st.columns(2)

    with col_left:
        st.subheader("Phenotype (HPO Terms)")
        hpo_input = st.text_area(
            "HPO Terms (one per line: ID | Label | Onset | Severity)",
            placeholder="HP:0001250 | Seizures | infantile | severe\nHP:0001263 | Global developmental delay | infantile | moderate\nHP:0000252 | Microcephaly | congenital | moderate",
            height=150,
            key="intake_hpo",
        )

        clinical_notes = st.text_area(
            "Clinical Notes",
            placeholder="Free-text clinical observations, history, prior workup results...",
            height=100,
            key="intake_notes",
        )

        family_history = st.text_area(
            "Family History",
            placeholder="Consanguinity, affected relatives, known genetic conditions...",
            height=80,
            key="intake_family",
        )

    with col_right:
        st.subheader("Genotype (Variants)")
        variant_input = st.text_area(
            "Genomic Variants (one per line: Gene | Variant | Zygosity)",
            placeholder="MECP2 | c.473C>T (p.Thr158Met) | heterozygous\nSCN1A | c.2836C>T (p.Arg946Cys) | heterozygous",
            height=120,
            key="intake_variants",
        )

        vcf_file = st.file_uploader(
            "Or upload VCF file",
            type=["vcf", "vcf.gz"],
            key="intake_vcf",
        )

        icol1, icol2 = st.columns(2)
        with icol1:
            pt_age = st.number_input("Age (years)", min_value=0.0, max_value=120.0, value=5.0, step=0.5, key="intake_age")
            pt_sex = st.selectbox("Sex", ["unknown", "male", "female"], key="intake_sex")
        with icol2:
            max_dx = st.slider("Max Diagnoses", 5, 50, 10, key="intake_max_dx")

    if st.button("Submit for Diagnosis", key="intake_submit"):
        # Parse HPO terms
        phenotypes = []
        for line in hpo_input.strip().split("\n"):
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split("|")]
            hpo = {"id": parts[0]}
            if len(parts) > 1:
                hpo["label"] = parts[1]
            if len(parts) > 2:
                hpo["onset"] = parts[2]
            if len(parts) > 3:
                hpo["severity"] = parts[3]
            phenotypes.append(hpo)

        if not phenotypes:
            st.warning("Please enter at least one HPO term.")
        else:
            # Parse variants
            variants = []
            for line in variant_input.strip().split("\n"):
                if not line.strip():
                    continue
                parts = [p.strip() for p in line.split("|")]
                if len(parts) >= 2:
                    v = {"gene": parts[0], "variant": parts[1]}
                    if len(parts) > 2:
                        v["zygosity"] = parts[2]
                    variants.append(v)

            with st.spinner("Analyzing patient data for rare disease differential diagnosis..."):
                payload = {
                    "phenotypes": phenotypes,
                    "variants": variants if variants else None,
                    "age_years": pt_age,
                    "sex": pt_sex,
                    "family_history": family_history.strip() or None,
                    "clinical_notes": clinical_notes.strip() or None,
                    "max_diagnoses": max_dx,
                }
                result = api_post("/v1/diagnostic/diagnose", payload)

            if result:
                st.success("Diagnostic analysis complete. See Diagnostic Dashboard tab for results.")
                st.session_state["diagnosis_result"] = result
                st.session_state["patient_phenotypes"] = phenotypes
                st.session_state["patient_variants"] = variants


# =====================================================================
# Tab 2: Diagnostic Dashboard
# =====================================================================

with tab_diagnosis:
    st.header("Diagnostic Dashboard")

    result = st.session_state.get("diagnosis_result")
    if not result:
        st.info("Submit patient data in the Patient Intake tab to see diagnostic results here.")

        # Also allow free-text queries
        st.subheader("Quick Query")
        question = st.text_area(
            "Rare Disease Question",
            placeholder="e.g., What are the genetic causes of infantile-onset epileptic encephalopathy?",
            height=80,
            key="dx_question",
        )
        top_k = st.slider("Evidence passages", 1, 20, 5, key="dx_topk")

        if st.button("Search", key="dx_search"):
            if question.strip():
                with st.spinner("Searching rare disease knowledge base..."):
                    payload = {
                        "question": question.strip(),
                        "top_k": top_k,
                    }
                    qresult = api_post("/v1/diagnostic/query", payload)

                if qresult:
                    st.subheader("Answer")
                    st.markdown(qresult.get("answer", "No answer generated."))

                    confidence = qresult.get("confidence", 0)
                    st.progress(confidence, text=f"Confidence: {confidence:.0%}")

                    evidence = qresult.get("evidence", [])
                    if evidence:
                        st.subheader(f"Evidence ({len(evidence)} passages)")
                        for i, ev in enumerate(evidence):
                            with st.expander(f"[{ev.get('collection', 'unknown')}] Score: {ev.get('score', 0):.3f}"):
                                st.write(ev.get("text", ""))
                                if ev.get("metadata"):
                                    st.json(ev["metadata"])
            else:
                st.warning("Please enter a question.")
    else:
        # Display differential diagnosis results
        st.subheader("Phenotype Summary")
        st.markdown(result.get("phenotype_summary", "No summary available."))

        odyssey_risk = result.get("diagnostic_odyssey_risk", "unknown")
        risk_color = NVIDIA_THEME["success"] if odyssey_risk == "low" else (
            NVIDIA_THEME["warning"] if odyssey_risk == "moderate" else NVIDIA_THEME["danger"]
        )
        st.markdown(f"**Diagnostic Odyssey Risk:** <span style='color:{risk_color};font-weight:bold'>{odyssey_risk.upper()}</span>", unsafe_allow_html=True)

        differential = result.get("differential", [])
        if differential:
            st.subheader(f"Differential Diagnosis ({len(differential)} candidates)")
            for i, dx in enumerate(differential):
                confidence = dx.get("confidence", 0)
                conf_color = NVIDIA_THEME["success"] if confidence >= 0.7 else (
                    NVIDIA_THEME["warning"] if confidence >= 0.4 else NVIDIA_THEME["danger"]
                )
                with st.expander(f"#{i+1} {dx.get('disease_name', 'Unknown')} | {dx.get('disease_id', '')} | Confidence: {confidence:.1%}"):
                    dcol1, dcol2, dcol3 = st.columns(3)
                    with dcol1:
                        st.metric("Confidence", f"{confidence:.1%}")
                    with dcol2:
                        st.metric("Phenotype Overlap", f"{dx.get('phenotype_overlap', 0):.1%}")
                    with dcol3:
                        st.metric("Gene Match", "Yes" if dx.get("gene_match") else "No")

                    st.write(f"**Inheritance:** {dx.get('inheritance_pattern', 'N/A')}")
                    st.write(f"**Source:** {dx.get('source', 'N/A')}")

                    evidence = dx.get("evidence", [])
                    if evidence:
                        st.write("**Evidence:**")
                        for ev in evidence:
                            st.write(f"- {ev}")

        recommendations = result.get("recommendations", [])
        if recommendations:
            st.subheader("Recommendations")
            for rec in recommendations:
                st.write(f"- {rec}")

        if st.button("Clear Results", key="dx_clear"):
            st.session_state.pop("diagnosis_result", None)
            st.rerun()


# =====================================================================
# Tab 3: Variant Review
# =====================================================================

with tab_variants:
    st.header("Variant Review")
    st.write("ACMG/AMP variant classification with ClinVar and in-silico predictor integration.")

    vcol1, vcol2 = st.columns(2)

    with vcol1:
        var_gene = st.text_input("Gene Symbol", placeholder="e.g., CFTR", key="var_gene")
        var_variant = st.text_input("Variant (HGVS)", placeholder="e.g., c.1521_1523delCTT", key="var_variant")
        var_transcript = st.text_input("Transcript (optional)", placeholder="e.g., NM_000492.4", key="var_transcript")

    with vcol2:
        var_zygosity = st.selectbox("Zygosity", ["", "heterozygous", "homozygous", "hemizygous"], key="var_zyg")
        var_inheritance = st.selectbox("Inheritance Pattern", ["", "AD", "AR", "XL", "XR", "MT"], key="var_inh")
        var_freq = st.text_input("gnomAD Allele Frequency (optional)", placeholder="e.g., 0.00001", key="var_freq")

    var_phenotypes = st.text_input(
        "Patient HPO IDs (comma-separated, optional)",
        placeholder="HP:0001250, HP:0001263",
        key="var_pheno",
    )

    if st.button("Classify Variant", key="var_classify"):
        if var_gene.strip() and var_variant.strip():
            with st.spinner("Classifying variant per ACMG/AMP guidelines..."):
                payload = {
                    "gene": var_gene.strip(),
                    "variant": var_variant.strip(),
                    "transcript": var_transcript.strip() or None,
                    "zygosity": var_zygosity or None,
                    "inheritance": var_inheritance or None,
                    "phenotypes": [p.strip() for p in var_phenotypes.split(",") if p.strip()] or None,
                    "population_frequency": float(var_freq) if var_freq.strip() else None,
                }
                result = api_post("/v1/diagnostic/variants/interpret", payload)

            if result:
                st.subheader("Classification Result")

                cls = result.get("classification", "VUS")
                cls_colors = {
                    "Pathogenic": NVIDIA_THEME["danger"],
                    "Likely Pathogenic": NVIDIA_THEME["warning"],
                    "VUS": NVIDIA_THEME["info"],
                    "Likely Benign": NVIDIA_THEME["success"],
                    "Benign": NVIDIA_THEME["success"],
                }
                cls_color = cls_colors.get(cls, NVIDIA_THEME["text_primary"])

                rcol1, rcol2, rcol3 = st.columns(3)
                with rcol1:
                    st.markdown(f"**Classification:** <span style='color:{cls_color};font-size:1.3em;font-weight:bold'>{cls}</span>", unsafe_allow_html=True)
                with rcol2:
                    st.metric("Pathogenicity Score", f"{result.get('pathogenicity_score', 0):.3f}")
                with rcol3:
                    clinvar = result.get("clinvar_significance", "Not found")
                    st.metric("ClinVar", clinvar or "Not found")

                criteria = result.get("acmg_criteria", [])
                if criteria:
                    st.subheader("ACMG Criteria Applied")
                    for c in criteria:
                        met_icon = "+" if c.get("met") else "-"
                        st.write(f"  {met_icon} **{c.get('code', '')}** ({c.get('strength', '')}): {c.get('evidence', '')}")

                recommendations = result.get("recommendations", [])
                if recommendations:
                    st.subheader("Recommendations")
                    for rec in recommendations:
                        st.write(f"- {rec}")
        else:
            st.warning("Please enter both gene symbol and variant notation.")

    # Saved variants from intake
    saved_variants = st.session_state.get("patient_variants", [])
    if saved_variants:
        st.markdown("---")
        st.subheader("Variants from Patient Intake")
        for v in saved_variants:
            st.write(f"- **{v.get('gene', '')}** {v.get('variant', '')} ({v.get('zygosity', 'N/A')})")

    # ACMG reference
    st.markdown("---")
    if st.button("Show ACMG Criteria Reference", key="var_acmg_ref"):
        data = api_get("/v1/diagnostic/acmg-criteria")
        if data:
            st.subheader("Pathogenic Criteria")
            for c in data.get("pathogenic_criteria", []):
                st.write(f"- **{c['code']}** ({c['strength']}): {c['description']}")
            st.subheader("Benign Criteria")
            for c in data.get("benign_criteria", []):
                st.write(f"- **{c['code']}** ({c['strength']}): {c['description']}")


# =====================================================================
# Tab 4: Therapeutic Options
# =====================================================================

with tab_therapy:
    st.header("Therapeutic Options")
    st.write("Search approved therapies, gene therapies, clinical trials, and investigational compounds.")

    tcol1, tcol2 = st.columns(2)

    with tcol1:
        th_disease_name = st.text_input(
            "Disease Name",
            placeholder="e.g., Spinal Muscular Atrophy",
            key="th_disease",
        )
        th_disease_id = st.text_input(
            "Disease ID (OMIM/Orphanet, optional)",
            placeholder="e.g., OMIM:253300 or ORPHA:70",
            key="th_id",
        )
        th_gene = st.text_input(
            "Gene (optional)",
            placeholder="e.g., SMN1",
            key="th_gene",
        )

    with tcol2:
        th_approved = st.checkbox("Include approved therapies", value=True, key="th_appr")
        th_investigational = st.checkbox("Include investigational", value=True, key="th_inv")
        th_gene_therapy = st.checkbox("Include gene therapies", value=True, key="th_gt")
        th_trials = st.checkbox("Include clinical trials", value=True, key="th_trials")
        th_max = st.slider("Maximum results", 5, 50, 20, key="th_max")

    if st.button("Search Therapies", key="th_search"):
        if th_disease_name.strip() or th_disease_id.strip() or th_gene.strip():
            with st.spinner("Searching therapeutic options..."):
                payload = {
                    "disease_name": th_disease_name.strip() or None,
                    "disease_id": th_disease_id.strip() or None,
                    "gene": th_gene.strip() or None,
                    "include_approved": th_approved,
                    "include_investigational": th_investigational,
                    "include_gene_therapy": th_gene_therapy,
                    "include_trials": th_trials,
                    "max_results": th_max,
                }
                result = api_post("/v1/diagnostic/therapy/search", payload)

            if result:
                tcol1, tcol2 = st.columns(2)
                with tcol1:
                    st.metric("Gene Therapy Eligible", "Yes" if result.get("gene_therapy_eligible") else "No")
                with tcol2:
                    st.metric("ERT Available", "Yes" if result.get("ert_available") else "No")

                therapies = result.get("therapies", [])
                if therapies:
                    st.subheader(f"Therapies Found ({len(therapies)})")
                    for th in therapies:
                        with st.expander(f"{th.get('therapy_name', 'Unknown')} | {th.get('therapy_type', '')} | {th.get('status', '')}"):
                            st.write(f"**Mechanism:** {th.get('mechanism', 'N/A')}")
                            st.write(f"**Target:** {th.get('target', 'N/A')}")
                            st.write(f"**Indication:** {th.get('indication', 'N/A')}")
                            st.write(f"**Evidence Level:** {th.get('evidence_level', 'N/A')}")

                trials = result.get("clinical_trials", [])
                if trials:
                    st.subheader(f"Clinical Trials ({len(trials)})")
                    for t in trials:
                        if isinstance(t, dict):
                            st.write(f"- **{t.get('trial_id', '')}**: {t.get('title', 'N/A')}")

                recommendations = result.get("recommendations", [])
                if recommendations:
                    st.subheader("Recommendations")
                    for rec in recommendations:
                        st.write(f"- {rec}")
        else:
            st.warning("Please enter a disease name, ID, or gene.")

    # Trial matching
    st.markdown("---")
    st.subheader("Clinical Trial Eligibility Matching")

    trial_location = st.text_input("Geographic Location (optional)", placeholder="e.g., Boston, MA, USA", key="th_location")

    if st.button("Find Matching Trials", key="th_trial_match"):
        if th_disease_name.strip() or th_gene.strip():
            with st.spinner("Matching to clinical trials..."):
                payload = {
                    "disease_name": th_disease_name.strip() or None,
                    "disease_id": th_disease_id.strip() or None,
                    "gene": th_gene.strip() or None,
                    "geographic_location": trial_location.strip() or None,
                    "max_results": 10,
                }
                result = api_post("/v1/diagnostic/trial/match", payload)

            if result:
                st.write(result.get("patient_summary", ""))
                matches = result.get("matches", [])
                if matches:
                    for m in matches:
                        with st.expander(f"{m.get('trial_id', 'N/A')} | {m.get('title', '')[:60]} | Score: {m.get('match_score', 0):.2f}"):
                            st.write(f"**Phase:** {m.get('phase', 'N/A')}")
                            st.write(f"**Status:** {m.get('status', 'N/A')}")
                            st.write(f"**Sponsor:** {m.get('sponsor', 'N/A')}")
                            st.write(f"**Eligibility:** {m.get('eligibility_summary', 'N/A')}")
        else:
            st.warning("Please enter a disease name or gene above.")

    # Gene therapy reference
    st.markdown("---")
    if st.button("Show Approved Gene Therapies", key="th_gt_ref"):
        data = api_get("/v1/diagnostic/gene-therapies")
        if data:
            st.subheader(f"FDA/EMA Approved Gene Therapies ({data.get('total_approved', 0)})")
            for gt in data.get("approved_therapies", []):
                st.write(f"- **{gt['name']}** -- {gt['disease']} ({gt['gene']}, {gt['year']})")
            st.caption(f"Pipeline estimate: ~{data.get('pipeline_estimate', 0)} gene therapies in development")


# =====================================================================
# Tab 5: Report Generator
# =====================================================================

with tab_reports:
    st.header("Report Generator")
    st.write("Generate structured rare disease diagnostic reports in multiple formats.")

    rcol1, rcol2 = st.columns(2)

    with rcol1:
        report_type = st.selectbox(
            "Report Type",
            [
                "differential_diagnosis",
                "variant_interpretation",
                "phenotype_analysis",
                "therapeutic_summary",
                "trial_eligibility",
                "case_summary",
                "gene_panel",
                "natural_history",
            ],
            key="rpt_type",
        )
        report_format = st.selectbox(
            "Export Format",
            ["markdown", "json", "pdf", "fhir", "phenopacket"],
            key="rpt_format",
        )

    with rcol2:
        report_title = st.text_input(
            "Report Title (optional)",
            placeholder="e.g., Diagnostic Evaluation - Patient 12345",
            key="rpt_title",
        )
        report_patient_id = st.text_input("Patient ID (optional)", key="rpt_patient")
        report_disease_id = st.text_input("Disease ID (optional)", key="rpt_disease")

    report_include_evidence = st.checkbox("Include evidence citations", value=True, key="rpt_evidence")
    report_include_recs = st.checkbox("Include recommendations", value=True, key="rpt_recs")

    # Use diagnosis result if available
    report_data = {}
    diagnosis_result = st.session_state.get("diagnosis_result")
    if diagnosis_result:
        st.info("Diagnosis results from Patient Intake are available and will be included.")
        report_data = diagnosis_result

    if st.button("Generate Report", key="rpt_generate"):
        with st.spinner("Generating report..."):
            payload = {
                "report_type": report_type,
                "format": report_format,
                "patient_id": report_patient_id.strip() or None,
                "disease_id": report_disease_id.strip() or None,
                "title": report_title.strip() or None,
                "data": report_data,
                "include_evidence": report_include_evidence,
                "include_recommendations": report_include_recs,
            }
            result = api_post("/v1/reports/generate", payload)

        if result:
            st.subheader(f"Report: {result.get('title', 'Untitled')}")
            st.caption(f"Format: {result.get('format', '')} | ID: {result.get('report_id', '')} | Generated: {result.get('generated_at', '')}")

            content = result.get("content", "")

            if report_format == "markdown":
                st.markdown(content)
            elif report_format in ("json", "fhir", "phenopacket"):
                try:
                    st.json(json.loads(content))
                except Exception:
                    st.code(content, language="json")
            else:
                st.code(content)

            # Download button
            ext_map = {"markdown": ".md", "json": ".json", "pdf": ".pdf", "fhir": ".json", "phenopacket": ".json"}
            mime_map = {"markdown": "text/markdown", "json": "application/json", "pdf": "application/pdf", "fhir": "application/fhir+json", "phenopacket": "application/json"}
            st.download_button(
                "Download Report",
                data=content,
                file_name=f"rd_report_{result.get('report_id', 'report')}{ext_map.get(report_format, '.txt')}",
                mime=mime_map.get(report_format, "text/plain"),
            )

    # Available formats reference
    st.markdown("---")
    if st.button("Show Available Formats", key="rpt_formats"):
        data = api_get("/v1/reports/formats")
        if data:
            st.subheader("Supported Formats")
            for f in data.get("formats", []):
                st.write(f"- **{f['name']}** ({f['extension']}): {f['description']}")
            st.subheader("Report Types")
            for rt in data.get("report_types", []):
                st.write(f"- {rt.replace('_', ' ').title()}")
