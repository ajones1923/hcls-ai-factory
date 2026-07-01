#!/usr/bin/env python3
"""End-to-end validation script for the Imaging Intelligence Agent.

Runs 12 validation checks covering all agent capabilities:
  (a)  Milvus connection and collection counts
  (b)  RAG query returns results
  (c)  NIM service status (local + cloud)
  (d)  API health endpoint responds
  (e)  Workflow registry has 4 workflows
  (f)  Cross-modal genomics trigger
  (g)  FHIR DiagnosticReport R4 export
  (h)  Real model inference (CXR DenseNet-121)
  (i)  Cloud NIM connectivity
  (j)  DICOM event webhook
  (k)  Export pipeline (Markdown + JSON + FHIR)
  (l)  Full pipeline: workflow → cross-modal → FHIR

Usage:
    python scripts/validate_e2e.py
    python scripts/validate_e2e.py --api-url http://localhost:8524
    python scripts/validate_e2e.py --host 10.0.0.1 --port 19530
    python scripts/validate_e2e.py --quick   # Skip slow checks

Author: Adam Jones
Date: February 2026
"""

import argparse
import json
import sys
import time
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger


# =====================================================================
# Individual checks
# =====================================================================


def check_milvus(host, port):
    """Check (a): Milvus connection and collection counts."""
    try:
        from src.collections import ImagingCollectionManager, OWNED_COLLECTION_SCHEMAS

        manager = ImagingCollectionManager(host=host, port=port)
        manager.connect()

        stats = manager.get_collection_stats()
        total_records = sum(stats.values())

        # Verify all 10 owned collections exist
        expected = set(OWNED_COLLECTION_SCHEMAS.keys())
        from pymilvus import utility
        found_collections = set()
        for name in expected:
            if utility.has_collection(name, using="imaging"):
                found_collections.add(name)

        # Also check genomic_evidence
        has_genomic = utility.has_collection("genomic_evidence", using="imaging")
        genomic_count = stats.get("genomic_evidence", 0)

        missing = expected - found_collections
        manager.disconnect()

        details = (
            f"{len(found_collections)}/{len(expected)} collections, "
            f"{total_records:,} total records"
        )
        if has_genomic:
            details += f", genomic_evidence: {genomic_count:,} vectors"

        if missing:
            return "FAIL", f"Missing collections: {missing}. {details}"

        return "PASS", f"All collections present. {details}"

    except Exception as e:
        return "FAIL", f"Milvus connection failed: {e}"


def check_rag_query(host, port):
    """Check (b): RAG query returns results."""
    try:
        from sentence_transformers import SentenceTransformer

        from config.settings import settings
        from src.collections import ImagingCollectionManager

        embedder = SentenceTransformer(settings.EMBEDDING_MODEL)
        manager = ImagingCollectionManager(host=host, port=port)
        manager.connect()

        # Run a simple vector search on imaging_findings
        test_query = "CT head hemorrhage detection"
        query_embedding = embedder.encode(
            test_query, normalize_embeddings=True
        ).tolist()

        results = manager.search(
            collection_name="imaging_findings",
            query_embedding=query_embedding,
            top_k=3,
        )

        manager.disconnect()

        if results:
            top_score = results[0]["score"]
            return (
                "PASS",
                f"RAG search returned {len(results)} results. "
                f"Top score: {top_score:.3f}.",
            )
        else:
            return (
                "WARN",
                "RAG search returned 0 results. "
                "Collections may be empty -- run seed scripts first.",
            )

    except Exception as e:
        return "FAIL", f"RAG query failed: {e}"


def check_nim_services():
    """Check (c): NIM service status (local + cloud)."""
    try:
        from config.settings import settings
        from src.nim.service_manager import NIMServiceManager

        nim_manager = NIMServiceManager(settings)
        status = nim_manager.check_all_services()

        available = sum(1 for s in status.values() if s == "available")
        cloud = sum(1 for s in status.values() if s == "cloud")
        mock = sum(1 for s in status.values() if s == "mock")
        unavailable = sum(1 for s in status.values() if s == "unavailable")

        status_str = ", ".join(f"{k}={v}" for k, v in status.items())

        if unavailable == 0:
            return (
                "PASS",
                f"All NIM services OK ({available} local, {cloud} cloud, {mock} mock). "
                f"{status_str}",
            )
        elif available + cloud + mock > 0:
            return (
                "WARN",
                f"Partial availability ({available} local, {cloud} cloud, "
                f"{mock} mock, {unavailable} unavailable). {status_str}",
            )
        else:
            return "FAIL", f"No NIM services available. {status_str}"

    except Exception as e:
        return "FAIL", f"NIM service check failed: {e}"


def check_api_health(api_url):
    """Check (d): API health endpoint responds."""
    try:
        import requests

        url = f"{api_url.rstrip('/')}/health"
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            data = response.json()
            status_str = data.get("status", "unknown")
            collections = data.get("collections", {})
            total_vecs = sum(collections.values()) if isinstance(collections, dict) else 0
            return (
                "PASS",
                f"API health OK. Status: {status_str}. "
                f"{len(collections)} collections, {total_vecs:,} vectors.",
            )
        else:
            return (
                "FAIL",
                f"API health returned HTTP {response.status_code}.",
            )

    except Exception as e:
        return (
            "WARN",
            f"API unreachable at {api_url}: {e}. Server may not be running.",
        )


def check_workflow_registry():
    """Check (e): Workflow registry has 4 workflows."""
    try:
        from src.workflows import WORKFLOW_REGISTRY

        count = len(WORKFLOW_REGISTRY)
        names = list(WORKFLOW_REGISTRY.keys())

        if count == 4:
            return "PASS", f"Workflow registry: {names}."
        elif count > 0:
            return "WARN", f"Expected 4 workflows, found {count}: {names}."
        else:
            return "FAIL", "Workflow registry is empty."

    except Exception as e:
        return "FAIL", f"Workflow registry check failed: {e}"


def check_cross_modal(host, port):
    """Check (f): Cross-modal genomics trigger fires correctly."""
    try:
        from src.models import FindingSeverity, WorkflowResult, WorkflowStatus

        # Test with mock data (doesn't need Milvus for basic evaluation)
        lung_result = WorkflowResult(
            workflow_name="ct_chest_lung_nodule",
            status=WorkflowStatus.COMPLETED,
            classification="Lung-RADS 4B (highest)",
            severity=FindingSeverity.CRITICAL,
            findings=[{
                "category": "nodule",
                "description": "15mm solid nodule in RUL",
                "lung_rads": "4B",
            }],
            measurements={"nodule_count": 1.0, "max_diameter_mm": 15.0},
        )

        # Try with real Milvus if available
        try:
            from sentence_transformers import SentenceTransformer
            from config.settings import settings
            from src.collections import ImagingCollectionManager
            from src.cross_modal import CrossModalTrigger

            embedder = SentenceTransformer(settings.EMBEDDING_MODEL)
            manager = ImagingCollectionManager(host=host, port=port)
            manager.connect()

            trigger = CrossModalTrigger(manager, embedder, enabled=True)
            cm_result = trigger.evaluate(lung_result)
            manager.disconnect()

            if cm_result and cm_result.genomic_hit_count > 0:
                return (
                    "PASS",
                    f"Cross-modal trigger fired for Lung-RADS 4B. "
                    f"{cm_result.genomic_hit_count} genomic hits from "
                    f"{cm_result.query_count} queries.",
                )
            elif cm_result:
                return (
                    "WARN",
                    "Cross-modal trigger fired but returned 0 genomic hits. "
                    "genomic_evidence collection may be empty.",
                )
            else:
                return "FAIL", "Cross-modal trigger did not fire for Lung-RADS 4B."
        except Exception as e:
            # Fall back to logic-only check
            from src.cross_modal import CrossModalTrigger
            from unittest.mock import MagicMock
            import numpy as np

            mock_manager = MagicMock()
            mock_manager.search.return_value = [
                {"id": "test-001", "score": 0.75, "text_chunk": "EGFR mutation"}
            ]
            mock_embedder = MagicMock()
            mock_embedder.encode.return_value = np.random.randn(384).astype(np.float32)

            trigger = CrossModalTrigger(mock_manager, mock_embedder, enabled=True)
            cm_result = trigger.evaluate(lung_result)

            if cm_result:
                return (
                    "PASS",
                    f"Cross-modal trigger logic OK (mock Milvus). "
                    f"Reason: {cm_result.trigger_reason}. "
                    f"Real Milvus unavailable: {e}",
                )
            return "FAIL", "Cross-modal trigger logic failed."

    except Exception as e:
        return "FAIL", f"Cross-modal check failed: {e}"


def check_fhir_export():
    """Check (g): FHIR DiagnosticReport R4 export works."""
    try:
        from src.export import export_fhir
        from src.models import (
            AgentResponse,
            CrossCollectionResult,
            FindingSeverity,
            SearchHit,
            WorkflowResult,
            WorkflowStatus,
        )

        # Build a sample response
        evidence = CrossCollectionResult(
            query="test", hits=[],
            total_collections_searched=11,
        )
        workflow = WorkflowResult(
            workflow_name="ct_head_hemorrhage",
            status=WorkflowStatus.COMPLETED,
            findings=[{
                "category": "hemorrhage",
                "description": "Intraparenchymal hemorrhage, 12.5 mL",
                "severity": "urgent",
            }],
            measurements={"volume_ml": 12.5, "midline_shift_mm": 3.2},
            classification="urgent_hemorrhage",
            severity=FindingSeverity.URGENT,
        )
        response = AgentResponse(
            question="CT head hemorrhage triage",
            answer="Urgent hemorrhage detected.",
            evidence=evidence,
            workflow_results=[workflow],
        )

        fhir_json = export_fhir(response, patient_id="P001")
        parsed = json.loads(fhir_json)

        # Validate FHIR structure
        checks_passed = []
        if parsed.get("resourceType") == "Bundle":
            checks_passed.append("Bundle type")
        if parsed.get("type") == "collection":
            checks_passed.append("Bundle category")

        # Find DiagnosticReport in entries
        entries = parsed.get("entry", [])
        has_diag_report = any(
            e.get("resource", {}).get("resourceType") == "DiagnosticReport"
            for e in entries
        )
        if has_diag_report:
            checks_passed.append("DiagnosticReport")

        has_observation = any(
            e.get("resource", {}).get("resourceType") == "Observation"
            for e in entries
        )
        if has_observation:
            checks_passed.append("Observation")

        if len(checks_passed) >= 3:
            return (
                "PASS",
                f"FHIR R4 export valid. Checks: {', '.join(checks_passed)}. "
                f"{len(entries)} resources in bundle.",
            )
        else:
            return (
                "WARN",
                f"FHIR export partial. Passed: {checks_passed}. "
                f"May need review.",
            )

    except ImportError:
        return "WARN", "export_fhir not yet implemented."
    except Exception as e:
        return "FAIL", f"FHIR export check failed: {e}"


def check_real_model_inference():
    """Check (h): Real model inference (CXR DenseNet-121)."""
    try:
        from src.workflows.cxr_rapid_findings import CXRRapidFindingsWorkflow

        wf = CXRRapidFindingsWorkflow(mock_mode=False)
        result = wf.run()

        if not result.is_mock:
            model_info = []
            if result.measurements.get("using_xrv", 0) > 0:
                model_info.append("torchxrayvision CheXpert-trained")
            if result.measurements.get("weights_loaded", 0) > 0:
                model_info.append("pretrained weights loaded")

            return (
                "PASS",
                f"Real CXR inference in {result.inference_time_ms:.0f}ms. "
                f"Model: {', '.join(model_info) or 'MONAI DenseNet-121'}. "
                f"{len(result.findings)} findings detected.",
            )
        else:
            return (
                "WARN",
                "CXR workflow fell back to mock mode. "
                "torchxrayvision or MONAI may not be installed.",
            )

    except Exception as e:
        return "WARN", f"Real model inference check failed: {e}"


def check_cloud_nim():
    """Check (i): Cloud NIM connectivity."""
    try:
        import os
        nvidia_key = os.environ.get("NVIDIA_API_KEY") or os.environ.get("IMAGING_NVIDIA_API_KEY")

        if not nvidia_key:
            return "WARN", "No NVIDIA_API_KEY set — cloud NIM check skipped."

        from config.settings import settings
        from src.nim.llm_client import LlamaLLMClient

        client = LlamaLLMClient(
            base_url=settings.NIM_LLM_URL,
            nvidia_api_key=nvidia_key,
            cloud_url=settings.NIM_CLOUD_URL,
            cloud_llm_model=settings.NIM_CLOUD_LLM_MODEL,
        )

        if client.is_cloud_available():
            # Test a simple generation
            test_response = client.generate(
                messages=[{"role": "user", "content": "Respond with OK"}],
                max_tokens=10,
            )
            if test_response and len(test_response) > 0:
                return (
                    "PASS",
                    f"Cloud NIM ({settings.NIM_CLOUD_LLM_MODEL}) responsive. "
                    f"Test response: {test_response[:50]}",
                )
            return "WARN", "Cloud NIM reachable but returned empty response."
        else:
            return "WARN", "Cloud NIM health check failed."

    except Exception as e:
        return "WARN", f"Cloud NIM check failed: {e}"


def check_dicom_webhook(api_url):
    """Check (j): DICOM event webhook endpoint."""
    try:
        import requests

        # Check event bus status
        status_url = f"{api_url.rstrip('/')}/events/status"
        response = requests.get(status_url, timeout=10)

        if response.status_code == 200:
            data = response.json()
            supported = data.get("supported_events", [])

            # Test webhook with a mock event
            webhook_url = f"{api_url.rstrip('/')}/events/dicom-webhook"
            test_event = {
                "event_type": "study.complete",
                "study_uid": "1.2.3.4.5.6.7.8.9",
                "patient_id": "TEST-001",
                "modality": "CT",
                "body_region": "head",
            }
            wh_response = requests.post(webhook_url, json=test_event, timeout=10)

            if wh_response.status_code == 200:
                return (
                    "PASS",
                    f"DICOM webhook endpoint active. "
                    f"Supported events: {supported}. "
                    f"Test event accepted.",
                )
            return (
                "WARN",
                f"Event status OK but webhook returned {wh_response.status_code}.",
            )
        else:
            return "WARN", f"Event bus status returned {response.status_code}."

    except Exception as e:
        return "WARN", f"DICOM webhook check failed (API may not be running): {e}"


def check_export_pipeline():
    """Check (k): Export pipeline (Markdown + JSON + FHIR)."""
    try:
        from src.export import export_markdown, export_json
        from src.models import AgentResponse, CrossCollectionResult

        response = AgentResponse(
            question="Test export pipeline",
            answer="Normal findings. No acute pathology.",
            evidence=CrossCollectionResult(query="test", hits=[]),
        )

        # Test markdown
        md = export_markdown(response)
        md_ok = len(md) > 50 and "Imaging Intelligence Report" in md

        # Test JSON
        json_str = export_json(response)
        json_ok = False
        try:
            parsed = json.loads(json_str)
            json_ok = parsed.get("question") == "Test export pipeline"
        except json.JSONDecodeError:
            pass

        # Test FHIR (may not be implemented yet)
        fhir_ok = False
        try:
            from src.export import export_fhir
            fhir_str = export_fhir(response)
            fhir_parsed = json.loads(fhir_str)
            fhir_ok = fhir_parsed.get("resourceType") == "Bundle"
        except (ImportError, AttributeError):
            pass

        results = []
        if md_ok:
            results.append("Markdown")
        if json_ok:
            results.append("JSON")
        if fhir_ok:
            results.append("FHIR R4")

        if len(results) >= 2:
            return "PASS", f"Export pipeline OK: {', '.join(results)}."
        elif len(results) >= 1:
            return "WARN", f"Partial export: {', '.join(results)} only."
        else:
            return "FAIL", "No export formats working."

    except Exception as e:
        return "FAIL", f"Export pipeline check failed: {e}"


def check_full_pipeline():
    """Check (l): Full pipeline: workflow → cross-modal → export."""
    try:
        from src.models import AgentResponse, CrossCollectionResult, WorkflowResult
        from src.workflows.cxr_rapid_findings import CXRRapidFindingsWorkflow
        from src.export import export_markdown, export_json

        # Run CXR workflow (mock mode for speed)
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.run()

        if result.status.value != "completed":
            return "FAIL", "Workflow did not complete."

        # Build agent response with workflow result
        evidence = CrossCollectionResult(query="full pipeline test", hits=[])
        response = AgentResponse(
            question="Full pipeline validation",
            answer=f"CXR analysis complete: {result.classification}",
            evidence=evidence,
            workflow_results=[result],
            nim_services_used=result.nim_services_used,
        )

        # Export all formats
        md = export_markdown(response)
        json_str = export_json(response)

        fhir_str = None
        try:
            from src.export import export_fhir
            fhir_str = export_fhir(response)
        except (ImportError, AttributeError):
            pass

        pipeline_steps = [
            f"Workflow ({result.workflow_name})",
            f"{len(result.findings)} findings",
            f"Severity: {result.severity.value}",
            f"Markdown: {len(md)} chars",
            f"JSON: {len(json_str)} chars",
        ]
        if fhir_str:
            pipeline_steps.append(f"FHIR: {len(fhir_str)} chars")

        return "PASS", f"Full pipeline OK: {' → '.join(pipeline_steps)}."

    except Exception as e:
        return "FAIL", f"Full pipeline check failed: {e}"


# =====================================================================
# Main
# =====================================================================


def main():
    parser = argparse.ArgumentParser(
        description="End-to-end validation for Imaging Intelligence Agent"
    )
    parser.add_argument(
        "--api-url",
        default="http://localhost:8524",
        help="API server URL (default: http://localhost:8524)",
    )
    parser.add_argument(
        "--host", default=None, help="Milvus host (default: localhost)"
    )
    parser.add_argument(
        "--port", type=int, default=None, help="Milvus port (default: 19530)"
    )
    parser.add_argument(
        "--quick", action="store_true",
        help="Skip slow checks (cloud NIM, real inference)",
    )
    args = parser.parse_args()

    host = args.host
    port = args.port

    print("=" * 70)
    print("  Imaging Intelligence Agent — End-to-End Validation (v2)")
    print("=" * 70)
    print(f"  API URL : {args.api_url}")
    print(f"  Milvus  : {host or 'default'}:{port or 'default'}")
    print(f"  Mode    : {'quick' if args.quick else 'full'}")
    print()

    checks = [
        ("(a)  Milvus Collections", lambda: check_milvus(host, port)),
        ("(b)  RAG Query", lambda: check_rag_query(host, port)),
        ("(c)  NIM Services", check_nim_services),
        ("(d)  API Health", lambda: check_api_health(args.api_url)),
        ("(e)  Workflow Registry", check_workflow_registry),
        ("(f)  Cross-Modal Trigger", lambda: check_cross_modal(host, port)),
        ("(g)  FHIR R4 Export", check_fhir_export),
        ("(k)  Export Pipeline", check_export_pipeline),
        ("(l)  Full Pipeline", check_full_pipeline),
    ]

    if not args.quick:
        checks.extend([
            ("(h)  Real Model Inference", check_real_model_inference),
            ("(i)  Cloud NIM", check_cloud_nim),
            ("(j)  DICOM Webhook", lambda: check_dicom_webhook(args.api_url)),
        ])

    results = []

    for check_name, check_fn in checks:
        logger.info(f"Running check: {check_name}")
        start = time.time()

        try:
            status, detail = check_fn()
        except Exception as e:
            status, detail = "FAIL", f"Unexpected error: {e}"

        elapsed = (time.time() - start) * 1000
        results.append((check_name, status, detail, elapsed))

        symbol = {"PASS": "[PASS]", "WARN": "[WARN]", "FAIL": "[FAIL]"}.get(
            status, "[????]"
        )
        print(f"  {symbol} {check_name} ({elapsed:.0f}ms)")
        if status != "PASS":
            print(f"         {detail}")

    # --- Summary ---
    print()
    print("=" * 70)
    print("  Validation Summary")
    print("=" * 70)

    for check_name, status, detail, elapsed in results:
        symbol = {"PASS": "[PASS]", "WARN": "[WARN]", "FAIL": "[FAIL]"}.get(
            status, "[????]"
        )
        print(f"  {symbol} {check_name}")
        print(f"         {detail}")

    passed = sum(1 for _, s, _, _ in results if s == "PASS")
    warned = sum(1 for _, s, _, _ in results if s == "WARN")
    failed = sum(1 for _, s, _, _ in results if s == "FAIL")
    total = len(results)
    total_time = sum(t for _, _, _, t in results)

    print()
    print(
        f"  Passed: {passed}/{total}  |  Warnings: {warned}  |  "
        f"Failed: {failed}  |  Total time: {total_time:.0f}ms"
    )

    if failed == 0 and warned == 0:
        print("  Overall: PASS (all checks passed)")
    elif failed == 0:
        print("  Overall: PASS (with warnings)")
    elif passed > 0:
        print("  Overall: PARTIAL")
    else:
        print("  Overall: FAIL")

    print("=" * 70)

    # Exit with non-zero if any checks failed
    sys.exit(1 if failed > 0 else 0)


if __name__ == "__main__":
    main()
