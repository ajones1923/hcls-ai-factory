"""Reports & Export tab -- generate clinical imaging reports in multiple formats.

Calls the FastAPI backend at /reports/generate for text-based reports
and provides DICOM SR export via the DICOMSRExporter.
"""

import json
import streamlit as st
import requests

API_BASE = "http://localhost:8524"

QUICK_QUERIES = [
    "What are the ACR guidelines for acute chest pain imaging?",
    "Describe standard CT head hemorrhage detection protocols",
    "Lung nodule management per Lung-RADS criteria",
    "Breast screening guidelines for high-risk patients",
    "Radiation dose optimization strategies for pediatric CT",
]


def render():
    st.header("Reports & Export")
    st.caption("Generate structured radiology reports in multiple formats")

    # ── Format selector ──────────────────────────────────────────────
    format_options = {
        "markdown": "Markdown (Human-readable)",
        "json": "JSON (Machine-readable)",
        "pdf": "PDF (NVIDIA-themed clinical report)",
        "fhir": "FHIR R4 DiagnosticReport Bundle",
        "dicom_sr": "DICOM Structured Report",
    }

    selected_format = st.selectbox(
        "Export Format",
        list(format_options.keys()),
        format_func=lambda x: format_options[x],
    )

    st.divider()

    # ── Workflow result export ────────────────────────────────────────
    last_result = st.session_state.get("last_workflow_result")
    if last_result:
        st.success("A workflow result is available for export.")
        with st.expander("Workflow Result Preview"):
            if isinstance(last_result, dict):
                st.json(last_result)
            else:
                st.write(last_result)

    # ── Quick report section ─────────────────────────────────────────
    st.subheader("Quick Reports")
    st.markdown("Select a pre-built query to generate a report instantly:")

    quick_cols = st.columns(2)
    for i, query in enumerate(QUICK_QUERIES):
        col = quick_cols[i % 2]
        with col:
            if st.button(query, key=f"quick_{i}", use_container_width=True):
                st.session_state["report_query"] = query

    st.divider()

    # ── Custom query input ───────────────────────────────────────────
    st.subheader("Generate Report")

    query = st.text_area(
        "Clinical Question",
        value=st.session_state.get("report_query", ""),
        height=100,
        placeholder="Enter a clinical imaging question for the report...",
    )

    col_a, col_b = st.columns(2)
    with col_a:
        report_title = st.text_input(
            "Report Title (optional)",
            placeholder="Custom title for the report",
        )
    with col_b:
        top_k = st.slider("Evidence items per collection", 1, 20, 5)

    # ── Generate button ──────────────────────────────────────────────
    generate_clicked = st.button(
        f"Generate {format_options[selected_format]} Report",
        type="primary",
        use_container_width=True,
    )

    if generate_clicked:
        if not query:
            st.warning("Please enter a clinical question.")
            return

        # ── DICOM SR export (handled separately) ─────────────────────
        if selected_format == "dicom_sr":
            _handle_dicom_sr_export(query)
            return

        # ── FHIR export (generated client-side from last workflow) ───
        if selected_format == "fhir":
            _handle_fhir_export(query, last_result)
            return

        # ── Standard formats (markdown, json, pdf) via API ───────────
        with st.spinner("Generating report..."):
            api_format = selected_format  # "markdown", "json", or "pdf"
            try:
                payload = {
                    "question": query,
                    "format": api_format,
                    "top_k": top_k,
                    "include_evidence": True,
                }
                if report_title:
                    payload["report_title"] = report_title

                resp = requests.post(
                    f"{API_BASE}/reports/generate",
                    json=payload,
                    timeout=60,
                )

                if resp.status_code != 200:
                    st.error(f"API error ({resp.status_code}): {resp.text[:500]}")
                    return

            except requests.ConnectionError:
                st.error(
                    f"Cannot connect to the Clinical Imaging Engine API at {API_BASE}. "
                    "Ensure the backend is running."
                )
                return
            except Exception as e:
                st.error(f"Report generation error: {e}")
                return

        # ── Display results ──────────────────────────────────────────
        if api_format == "pdf":
            # PDF comes back as binary
            st.success("PDF report generated successfully.")
            st.download_button(
                "Download PDF Report",
                data=resp.content,
                file_name="imaging_intelligence_report.pdf",
                mime="application/pdf",
                use_container_width=True,
            )
        elif api_format == "json":
            data = resp.json()
            st.success(
                f"JSON report generated -- {data.get('evidence_count', 0)} evidence items, "
                f"{data.get('generation_time_ms', 0):.0f}ms"
            )
            st.json(data)
            st.download_button(
                "Download JSON",
                data=json.dumps(data, indent=2),
                file_name="imaging_intelligence_report.json",
                mime="application/json",
                use_container_width=True,
            )
        else:
            # Markdown
            data = resp.json()
            content = data.get("content", data.get("answer", ""))
            st.success(
                f"Report generated -- {data.get('evidence_count', 0)} evidence items, "
                f"{data.get('generation_time_ms', 0):.0f}ms"
            )
            st.markdown(content)
            st.download_button(
                "Download Markdown",
                data=content,
                file_name="imaging_intelligence_report.md",
                mime="text/markdown",
                use_container_width=True,
            )

    # ── Report template browser ──────────────────────────────────────
    st.divider()
    st.subheader("Report Templates")
    try:
        from pathlib import Path

        templates_path = (
            Path(__file__).parent.parent.parent
            / "data"
            / "reference"
            / "report_template_seed_data.json"
        )
        with open(templates_path) as f:
            templates = json.load(f)

        for t in templates:
            with st.expander(
                f"**{t.get('template_name', 'Unknown')}** -- "
                f"{t.get('modality', '').upper()} | {t.get('body_region', '')}"
            ):
                st.markdown(f"**Finding Type:** {t.get('finding_type', 'N/A')}")
                st.markdown(f"**Coding System:** {t.get('coding_system', 'N/A')}")
                st.markdown(f"**Fields:** {t.get('structured_fields', 'N/A')}")
                st.text(t.get("example_report", "No example available"))
    except Exception as e:
        st.caption(f"Report templates not loaded: {e}")


def _handle_fhir_export(query: str, last_result):
    """Generate a FHIR R4 DiagnosticReport Bundle."""
    with st.spinner("Generating FHIR R4 bundle..."):
        # First get a markdown report for content
        try:
            resp = requests.post(
                f"{API_BASE}/reports/generate",
                json={"question": query, "format": "json", "top_k": 5},
                timeout=60,
            )
            if resp.status_code != 200:
                st.error(f"API error ({resp.status_code}): {resp.text[:500]}")
                return
            report_data = resp.json()
        except requests.ConnectionError:
            st.error(f"Cannot connect to the API at {API_BASE}.")
            return
        except Exception as e:
            st.error(f"Error: {e}")
            return

    # Build a minimal FHIR bundle client-side
    import uuid
    from datetime import datetime, timezone

    now_iso = datetime.now(timezone.utc).isoformat()
    patient_url = f"urn:uuid:{uuid.uuid4()}"

    bundle = {
        "resourceType": "Bundle",
        "id": str(uuid.uuid4()),
        "type": "collection",
        "timestamp": now_iso,
        "entry": [
            {
                "fullUrl": patient_url,
                "resource": {
                    "resourceType": "Patient",
                    "id": "anonymous",
                },
            },
            {
                "fullUrl": f"urn:uuid:{uuid.uuid4()}",
                "resource": {
                    "resourceType": "DiagnosticReport",
                    "id": str(uuid.uuid4()),
                    "status": "final",
                    "category": [
                        {
                            "coding": [
                                {
                                    "system": "http://loinc.org",
                                    "code": "LP29684-5",
                                    "display": "Radiology",
                                }
                            ]
                        }
                    ],
                    "code": {
                        "coding": [
                            {
                                "system": "http://loinc.org",
                                "code": "18748-4",
                                "display": "Diagnostic imaging study",
                            }
                        ],
                        "text": "Imaging Intelligence Agent Report",
                    },
                    "subject": {"reference": patient_url},
                    "effectiveDateTime": now_iso,
                    "issued": now_iso,
                    "conclusion": report_data.get("answer", ""),
                },
            },
        ],
    }

    fhir_json = json.dumps(bundle, indent=2)
    st.success("FHIR R4 DiagnosticReport Bundle generated.")
    st.json(bundle)
    st.download_button(
        "Download FHIR Bundle (JSON)",
        data=fhir_json,
        file_name="fhir_diagnostic_report.json",
        mime="application/fhir+json",
        use_container_width=True,
    )


def _handle_dicom_sr_export(query: str):
    """Generate a DICOM Structured Report via the DICOMSRExporter."""
    with st.spinner("Generating DICOM Structured Report..."):
        try:
            from src.export_dicom_sr import DICOMSRExporter

            exporter = DICOMSRExporter()
            if not exporter.available:
                st.warning(
                    "DICOM SR export requires `highdicom` and `pydicom`. "
                    "Install with: `pip install highdicom pydicom`"
                )
                return

            # Generate a mock SR based on the query
            workflow_name = "generic"
            query_lower = query.lower()
            if "head" in query_lower or "hemorrhage" in query_lower:
                workflow_name = "ct_head_hemorrhage"
            elif "chest" in query_lower or "cxr" in query_lower or "lung" in query_lower:
                if "nodule" in query_lower:
                    workflow_name = "ct_chest_lung_nodule"
                else:
                    workflow_name = "cxr_rapid_findings"
            elif "breast" in query_lower or "mammo" in query_lower:
                workflow_name = "breast_birads"

            import tempfile
            import os

            with tempfile.NamedTemporaryFile(suffix=".dcm", delete=False) as tmp:
                tmp_path = tmp.name

            result = exporter.export_mock(
                workflow_name=workflow_name,
                output_path=tmp_path,
            )

            if result.success:
                with open(tmp_path, "rb") as f:
                    dcm_bytes = f.read()
                os.unlink(tmp_path)

                st.success(
                    f"DICOM SR generated successfully. "
                    f"SOP Instance UID: `{result.sop_instance_uid}` | "
                    f"Content items: {result.content_items_count}"
                )
                st.download_button(
                    "Download DICOM SR (.dcm)",
                    data=dcm_bytes,
                    file_name=f"imaging_sr_{workflow_name}.dcm",
                    mime="application/dicom",
                    use_container_width=True,
                )
            else:
                st.error(f"DICOM SR generation failed: {result.error}")
                if tmp_path and os.path.exists(tmp_path):
                    os.unlink(tmp_path)

        except ImportError as e:
            st.warning(
                f"DICOM SR export dependencies not available: {e}. "
                "Install with: `pip install highdicom pydicom`"
            )
        except Exception as e:
            st.error(f"DICOM SR export error: {e}")
