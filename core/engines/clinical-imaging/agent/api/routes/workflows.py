"""Imaging workflow endpoints for Imaging Intelligence Agent.

Provides endpoints to list, inspect, and run reference imaging
analysis workflows (CT head hemorrhage, lung nodule, CXR triage,
MRI MS lesions).

Author: Adam Jones
Date: February 2026
"""

import time
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException
from loguru import logger
from pydantic import BaseModel, Field

from src.cross_modal import CrossModalTrigger
from src.workflows import WORKFLOW_REGISTRY

# =====================================================================
# Response Models
# =====================================================================


class WorkflowInfo(BaseModel):
    """Metadata for a single imaging workflow."""
    name: str
    modality: str
    body_region: str
    target_latency_sec: float
    models_used: List[str]
    mock_mode: bool


class WorkflowListResponse(BaseModel):
    """List of all available imaging workflows."""
    workflows: List[WorkflowInfo]
    total: int


class FindingDetail(BaseModel):
    """A single clinical finding from workflow execution."""
    category: str = ""
    description: str = ""
    severity: str = ""
    recommendation: str = ""
    extra: Dict[str, Any] = {}


class CrossModalResponse(BaseModel):
    """Cross-modal genomics enrichment attached to a workflow result."""
    trigger_reason: str = ""
    genomic_context: List[str] = []
    genomic_hit_count: int = 0
    query_count: int = 0
    enrichment_summary: str = ""


class WorkflowRunResponse(BaseModel):
    """Result from executing an imaging workflow."""
    workflow_name: str
    status: str
    findings: List[FindingDetail]
    measurements: Dict[str, float]
    classification: str
    severity: str
    inference_time_ms: float
    nim_services_used: List[str]
    is_mock: bool
    cross_modal: Optional[CrossModalResponse] = None


class WorkflowRunRequest(BaseModel):
    """Request to run an imaging workflow."""
    input_path: str = Field("", description="Path to input DICOM/NIfTI (empty for mock)")
    mock_mode: bool = Field(True, description="Run in mock mode with synthetic data")


# =====================================================================
# Router
# =====================================================================

router = APIRouter()


# =====================================================================
# Endpoints
# =====================================================================


@router.get("/workflows", response_model=WorkflowListResponse)
async def list_workflows():
    """List all available imaging analysis workflows.

    Returns metadata for each registered workflow including modality,
    body region, target latency, and models used.
    """
    workflows = []
    for name, wf_class in WORKFLOW_REGISTRY.items():
        try:
            instance = wf_class(mock_mode=True)
            info = instance.get_workflow_info()
            workflows.append(WorkflowInfo(
                name=info["name"],
                modality=info["modality"],
                body_region=info["body_region"],
                target_latency_sec=info["target_latency_sec"],
                models_used=info["models_used"],
                mock_mode=info["mock_mode"],
            ))
        except Exception as e:
            logger.warning(f"Failed to get info for workflow {name}: {e}")
            workflows.append(WorkflowInfo(
                name=name,
                modality="unknown",
                body_region="unknown",
                target_latency_sec=0.0,
                models_used=[],
                mock_mode=True,
            ))

    return WorkflowListResponse(
        workflows=workflows,
        total=len(workflows),
    )


@router.post("/workflow/{name}/run", response_model=WorkflowRunResponse)
async def run_workflow(name: str, request: WorkflowRunRequest):
    """Run a named imaging analysis workflow.

    Executes the specified workflow in mock or live mode.
    Mock mode returns clinically realistic synthetic results
    suitable for demonstrations and testing.

    When cross-modal integration is enabled and the workflow result
    meets a severity threshold (e.g., Lung-RADS 4A+), the response
    includes genomic context retrieved from the genomic_evidence
    collection.

    Available workflows:
      - ct_head_hemorrhage: Emergency intracranial hemorrhage triage
      - ct_chest_lung_nodule: Lung cancer screening nodule analysis
      - cxr_rapid_findings: Chest X-ray multi-label triage
      - mri_brain_ms_lesion: MS lesion quantification and activity
    """
    if name not in WORKFLOW_REGISTRY:
        available = list(WORKFLOW_REGISTRY.keys())
        raise HTTPException(
            status_code=404,
            detail=f"Workflow '{name}' not found. Available: {available}",
        )

    wf_class = WORKFLOW_REGISTRY[name]

    try:
        workflow = wf_class(mock_mode=request.mock_mode)
        result = workflow.run(input_path=request.input_path)
    except NotImplementedError as e:
        raise HTTPException(
            status_code=501,
            detail=f"Workflow requires live NIM deployment: {e}",
        )
    except Exception as e:
        logger.error(f"Workflow {name} execution failed: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"Workflow execution failed: {e}",
        )

    # Convert findings to FindingDetail objects
    finding_details = []
    for f in result.findings:
        extra = {k: v for k, v in f.items() if k not in ("category", "description", "severity", "recommendation")}
        finding_details.append(FindingDetail(
            category=f.get("category", ""),
            description=f.get("description", ""),
            severity=f.get("severity", ""),
            recommendation=f.get("recommendation", ""),
            extra=extra,
        ))

    # Evaluate cross-modal trigger
    cross_modal_response = None
    cross_modal_trigger = _get_cross_modal_trigger()
    if cross_modal_trigger:
        cm_result = cross_modal_trigger.evaluate(result)
        if cm_result:
            cross_modal_response = CrossModalResponse(
                trigger_reason=cm_result.trigger_reason,
                genomic_context=cm_result.genomic_context,
                genomic_hit_count=cm_result.genomic_hit_count,
                query_count=cm_result.query_count,
                enrichment_summary=cm_result.enrichment_summary,
            )
            logger.info(
                f"Cross-modal enrichment for {name}: "
                f"{cm_result.genomic_hit_count} genomic hits"
            )

    return WorkflowRunResponse(
        workflow_name=result.workflow_name,
        status=result.status.value,
        findings=finding_details,
        measurements=result.measurements,
        classification=result.classification,
        severity=result.severity.value,
        inference_time_ms=result.inference_time_ms,
        nim_services_used=result.nim_services_used,
        is_mock=result.is_mock,
        cross_modal=cross_modal_response,
    )


def _get_cross_modal_trigger() -> Optional[CrossModalTrigger]:
    """Retrieve the CrossModalTrigger from application state, if available.

    Returns None if the trigger has not been initialized (e.g., during
    unit tests or when cross-modal is disabled).
    """
    try:
        from api.main import _state
        return _state.get("cross_modal_trigger")
    except Exception:
        return None


@router.get("/workflow/{name}/info", response_model=WorkflowInfo)
async def workflow_info(name: str):
    """Get metadata for a specific imaging workflow.

    Returns modality, body region, target latency, models used,
    and current mode (mock or live).
    """
    if name not in WORKFLOW_REGISTRY:
        available = list(WORKFLOW_REGISTRY.keys())
        raise HTTPException(
            status_code=404,
            detail=f"Workflow '{name}' not found. Available: {available}",
        )

    wf_class = WORKFLOW_REGISTRY[name]

    try:
        instance = wf_class(mock_mode=True)
        info = instance.get_workflow_info()
        return WorkflowInfo(
            name=info["name"],
            modality=info["modality"],
            body_region=info["body_region"],
            target_latency_sec=info["target_latency_sec"],
            models_used=info["models_used"],
            mock_mode=info["mock_mode"],
        )
    except Exception as e:
        logger.error(f"Failed to get workflow info for {name}: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"Failed to retrieve workflow info: {e}",
        )
