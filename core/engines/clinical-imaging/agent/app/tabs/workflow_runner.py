"""Workflow Runner tab — execute imaging AI workflows and demo cases."""
import streamlit as st
import json
import time
from pathlib import Path

DEMO_CASES_PATH = Path(__file__).parent.parent.parent / "data" / "reference" / "demo_cases.json"
SAMPLE_DIR = Path(__file__).parent.parent.parent / "data" / "sample_images"

# ---------------------------------------------------------------------------
# Image mappings — demo cases and standalone workflows
# ---------------------------------------------------------------------------
WORKFLOW_IMAGE_MAP = {
    "DEMO-001": {
        "primary": "annotated/ct_head_hemorrhage_annotated.png",
        "context": ["fullres/fullres_cxr_synth_000.png"],
        "modality_icon": "\U0001f9e0",
        "modality_label": "Non-Contrast CT Head",
    },
    "DEMO-002": {
        "primary": "annotated/ct_chest_nodule_annotated.png",
        "context": ["fullres/fullres_cxr_synth_001.png", "fullres/fullres_cxr_synth_002.png"],
        "modality_icon": "\U0001fac1",
        "modality_label": "Low-Dose CT Chest",
    },
    "DEMO-003": {
        "primary": "annotated/ct_coronary_annotated.png",
        "context": ["organ_ct_000.png", "organ_ct_003.png"],
        "modality_icon": "\u2764\ufe0f",
        "modality_label": "Coronary CT Angiography",
    },
    "DEMO-004": {
        "primary": "annotated/cxr_pneumonia_bilateral_annotated.png",
        "context": ["fullres/fullres_cxr_synth_002.png"],
        "modality_icon": "\U0001fac1",
        "modality_label": "Portable AP Chest X-Ray",
    },
}

WORKFLOW_IMAGES = {
    "ct_head_hemorrhage": "annotated/ct_head_hemorrhage_annotated.png",
    "ct_chest_lung_nodule": "annotated/ct_chest_nodule_annotated.png",
    "ct_coronary_angiography": "annotated/ct_coronary_annotated.png",
    "cxr_rapid_findings": "annotated/fullres_cxr_synth_001_annotated.png",
    "mri_brain_ms_lesion": "annotated/ct_head_hemorrhage_annotated.png",  # reuse CT head
    "mri_prostate_pirads": "organ_ct_001.png",  # reuse organ CT
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _load_demo_cases():
    try:
        with open(DEMO_CASES_PATH) as f:
            return json.load(f)
    except Exception:
        return []


def _severity_color(severity: str) -> str:
    colors = {"critical": "\U0001f534", "urgent": "\U0001f7e0", "significant": "\U0001f7e1", "routine": "\U0001f7e2", "normal": "\u26aa"}
    return colors.get(severity.lower(), "\u26aa")


def _severity_bg(severity: str) -> str:
    """Return a CSS background-color string for inline severity badges."""
    bg = {
        "critical": "#d32f2f",
        "urgent": "#e65100",
        "significant": "#f9a825",
        "routine": "#2e7d32",
        "normal": "#757575",
    }
    return bg.get(severity.lower(), "#757575")


def _safe_image(relative_path: str) -> Path | None:
    """Return the resolved path only if the image exists on disk."""
    p = SAMPLE_DIR / relative_path
    return p if p.exists() else None


# ---------------------------------------------------------------------------
# Main render
# ---------------------------------------------------------------------------
def render():
    st.header("Workflow Runner")

    # ==================================================================
    # Demo Cases section
    # ==================================================================
    st.subheader("\U0001f4cb Pre-loaded Demo Cases")
    cases = _load_demo_cases()

    if cases:
        case_options = {c["case_id"]: f"{c['case_id']}: {c['title']}" for c in cases}
        selected_case_id = st.selectbox("Select Demo Case", list(case_options.keys()), format_func=lambda x: case_options[x])
        selected_case = next(c for c in cases if c["case_id"] == selected_case_id)

        sev = selected_case.get("expected_severity", "N/A")
        img_info = WORKFLOW_IMAGE_MAP.get(selected_case_id, {})

        # --- Case details card ---
        with st.container(border=True):
            st.markdown(f"**{selected_case['title']}**")
            st.markdown(f"_{selected_case['clinical_scenario']}_")

            col1, col2, col3 = st.columns(3)
            with col1:
                demo = selected_case["patient_demographics"]
                st.metric("Patient", f"{demo.get('age', 'N/A')}y {demo.get('sex', 'N/A')}")
            with col2:
                st.metric("Workflow", selected_case["workflow_name"])
            with col3:
                st.metric("Expected Severity", f"{_severity_color(sev)} {sev.upper()}")

        # --- Visual preview (before running) ---
        if img_info:
            primary_path = _safe_image(img_info.get("primary", ""))
            if primary_path:
                modality_label = img_info.get("modality_label", "Medical Imaging")
                modality_icon = img_info.get("modality_icon", "")
                st.image(
                    str(primary_path),
                    caption=f"{modality_icon} {modality_label} — Representative imaging for this case",
                    use_container_width=True,
                )

        # --- Run button ---
        if st.button("\u25b6 Run Demo Case", type="primary", use_container_width=True):
            # Animated pipeline progress
            progress_bar = st.progress(0)
            status_text = st.empty()

            pipeline_steps = [
                ("\U0001f4e5 Loading DICOM data...", 0.15, 0.4),
                ("\U0001f527 Preprocessing & normalization...", 0.30, 0.6),
                ("\U0001f9e0 Running AI inference...", 0.55, 1.2),
                ("\U0001f50d Detecting findings...", 0.75, 0.5),
                ("\U0001f4ca Generating classifications...", 0.90, 0.3),
                ("\u2705 Analysis complete!", 1.0, 0.2),
            ]

            for step_text, progress, delay in pipeline_steps:
                status_text.markdown(f"**{step_text}**")
                progress_bar.progress(progress)
                time.sleep(delay)

            progress_bar.empty()
            status_text.empty()

            overrides = selected_case.get("workflow_overrides", {})
            workflow_name = selected_case["workflow_name"]

            st.subheader("Workflow Results")

            # ==================================================
            # Top section: Image + Key Findings side-by-side
            # ==================================================
            primary_path = _safe_image(img_info.get("primary", "")) if img_info else None
            has_image = primary_path is not None

            if has_image:
                img_col, findings_col = st.columns([0.55, 0.45])
            else:
                # Fall back to full width findings when no image
                img_col = None
                findings_col = st.container()

            # -- Left: annotated image --
            if has_image and img_col is not None:
                with img_col:
                    st.image(
                        str(primary_path),
                        caption=f"AI Analysis Overlay — {workflow_name}",
                        use_container_width=True,
                    )
                    # Show context images in a row if available
                    ctx_paths = [_safe_image(p) for p in img_info.get("context", [])]
                    ctx_paths = [p for p in ctx_paths if p is not None]
                    if ctx_paths:
                        ctx_cols = st.columns(len(ctx_paths))
                        for idx, cp in enumerate(ctx_paths):
                            with ctx_cols[idx]:
                                st.image(str(cp), use_container_width=True)

            # -- Right: key metrics & classification --
            with findings_col:
                # Extract a handful of important metrics to feature
                metric_keys = []
                for key, val in overrides.items():
                    if isinstance(val, (int, float)):
                        metric_keys.append((key, val))
                    elif isinstance(val, str) and any(ch.isdigit() for ch in val):
                        metric_keys.append((key, val))
                # Show up to 4 metrics
                shown_metrics = metric_keys[:4]
                if shown_metrics:
                    m_cols = st.columns(min(len(shown_metrics), 2))
                    for idx, (key, val) in enumerate(shown_metrics):
                        with m_cols[idx % len(m_cols)]:
                            st.metric(key.replace("_", " ").title(), val)

                # Severity badge
                sev_bg = _severity_bg(sev)
                st.markdown(
                    f'<span style="background:{sev_bg};color:#fff;padding:4px 12px;'
                    f'border-radius:4px;font-weight:600;">{_severity_color(sev)} {sev.upper()}</span>',
                    unsafe_allow_html=True,
                )

                # Classification
                classification = selected_case.get("expected_classification", "N/A")
                st.markdown(f"**Classification:** {classification}")

            # ==================================================
            # Middle section: Full findings in expander
            # ==================================================
            with st.expander("Full AI Findings", expanded=False):
                for key, val in overrides.items():
                    if isinstance(val, (str, int, float, bool)):
                        st.markdown(f"- **{key.replace('_', ' ').title()}:** {val}")

            st.divider()

            # ==================================================
            # Export Results
            # ==================================================
            genomic = selected_case.get("genomic_context")

            st.subheader("\U0001f4c4 Export Results")
            export_cols = st.columns(4)

            # Generate report content from the workflow results
            report_title = f"Imaging Intelligence Report: {selected_case['title']}"
            report_question = selected_case['clinical_scenario']

            with export_cols[0]:
                # Markdown download
                md_content = f"# {report_title}\n\n## Clinical Scenario\n{report_question}\n\n## Findings\n"
                for key, val in overrides.items():
                    if isinstance(val, (str, int, float, bool)):
                        md_content += f"- **{key.replace('_', ' ').title()}:** {val}\n"
                md_content += f"\n## Classification\n- **Classification:** {selected_case.get('expected_classification', 'N/A')}\n"
                md_content += f"- **Severity:** {sev.upper()}\n"
                if genomic:
                    md_content += f"\n## Genomic Context\n- **Genes:** {', '.join(genomic.get('genes', []))}\n"
                    md_content += f"- **Relevance:** {genomic.get('relevance', '')}\n"
                md_content += f"\n---\n*Generated by HCLS AI Factory \u2014 Imaging Intelligence Agent*\n"
                st.download_button("\U0001f4dd Markdown", md_content, f"{selected_case['case_id']}_report.md", "text/markdown", use_container_width=True)

            with export_cols[1]:
                # JSON download
                import json as _json
                json_report = {
                    "case_id": selected_case["case_id"],
                    "title": selected_case["title"],
                    "clinical_scenario": selected_case["clinical_scenario"],
                    "patient_demographics": selected_case["patient_demographics"],
                    "workflow_name": selected_case["workflow_name"],
                    "findings": overrides,
                    "classification": selected_case.get("expected_classification"),
                    "severity": sev,
                    "genomic_context": genomic,
                }
                st.download_button("\U0001f4ca JSON", _json.dumps(json_report, indent=2), f"{selected_case['case_id']}_report.json", "application/json", use_container_width=True)

            with export_cols[2]:
                # PDF download via API
                st.download_button(
                    "\U0001f4d5 PDF Report",
                    md_content.encode(),  # Fallback - ideally call the API
                    f"{selected_case['case_id']}_report.pdf",
                    "application/pdf",
                    use_container_width=True,
                    disabled=True,
                    help="Generate via Reports & Export tab for NVIDIA-branded PDF"
                )

            with export_cols[3]:
                # FHIR R4 placeholder
                fhir_bundle = {
                    "resourceType": "Bundle",
                    "type": "transaction",
                    "entry": [
                        {
                            "resource": {
                                "resourceType": "DiagnosticReport",
                                "status": "final",
                                "code": {"coding": [{"system": "http://loinc.org", "code": "18748-4", "display": "Diagnostic imaging study"}]},
                                "conclusion": f"{selected_case['title']} - {sev.upper()}",
                                "category": [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "RAD"}]}],
                            }
                        }
                    ]
                }
                st.download_button("\U0001f3e5 FHIR R4", _json.dumps(fhir_bundle, indent=2), f"{selected_case['case_id']}_fhir.json", "application/json", use_container_width=True)

            st.divider()

            # ==================================================
            # Genomic context
            # ==================================================
            if genomic:
                st.subheader("\U0001f9ec Cross-Modal Genomic Enrichment")
                st.markdown(f"**Genes:** {', '.join(genomic.get('genes', []))}")
                st.markdown(f"**Relevance:** {genomic.get('relevance', '')}")

            # ==================================================
            # Talking points
            # ==================================================
            talking_points = selected_case.get("talking_points", [])
            if talking_points:
                st.subheader("\U0001f4ac Demo Talking Points")
                for tp in talking_points:
                    st.markdown(f"- {tp}")

    st.divider()

    # ==================================================================
    # Manual Workflow Runner
    # ==================================================================
    st.subheader("\U0001f527 Manual Workflow Runner")
    workflows = ["ct_head_hemorrhage", "ct_chest_lung_nodule", "cxr_rapid_findings",
                  "mri_brain_ms_lesion", "ct_coronary_angiography", "mri_prostate_pirads"]

    selected_wf = st.selectbox("Select Workflow", workflows)

    # Show a thumbnail preview for the selected workflow
    wf_image_rel = WORKFLOW_IMAGES.get(selected_wf)
    wf_image_path = _safe_image(wf_image_rel) if wf_image_rel else None
    if wf_image_path:
        st.image(str(wf_image_path), caption=f"Reference image — {selected_wf}", width=300)

    mock_mode = st.checkbox("Mock Mode (demo data)", value=True)

    if st.button("\u25b6 Run Workflow"):
        with st.spinner(f"Running {selected_wf}..."):
            try:
                from src.workflows import WORKFLOW_REGISTRY
                wf_cls = WORKFLOW_REGISTRY.get(selected_wf)
                if wf_cls:
                    wf = wf_cls()
                    result = wf.run(mock_mode=mock_mode)
                    st.success("Workflow complete")

                    # Show image alongside JSON result
                    if wf_image_path:
                        res_img_col, res_json_col = st.columns([0.45, 0.55])
                        with res_img_col:
                            st.image(
                                str(wf_image_path),
                                caption=f"AI Analysis Overlay — {selected_wf}",
                                use_container_width=True,
                            )
                        with res_json_col:
                            if hasattr(result, "dict"):
                                st.json(result.dict())
                            elif hasattr(result, "model_dump"):
                                st.json(result.model_dump())
                            else:
                                st.json(result)
                    else:
                        if hasattr(result, "dict"):
                            st.json(result.dict())
                        elif hasattr(result, "model_dump"):
                            st.json(result.model_dump())
                        else:
                            st.json(result)
                else:
                    st.error(f"Workflow {selected_wf} not found in registry")
            except Exception as e:
                st.error(f"Workflow error: {e}")
