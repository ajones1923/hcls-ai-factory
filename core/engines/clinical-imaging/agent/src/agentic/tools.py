"""Tool definitions for the Imaging Intelligence agentic reasoning loop.

Each tool wraps a capability of the Imaging Intelligence Agent (RAG search,
workflow execution, cross-modal genomics, cross-agent queries, radiomics,
report parsing) behind a uniform ImagingTool interface.

Tools can be registered with NVIDIA AIQ Toolkit when available, or used
standalone by the built-in reasoning loop.

Open source (AIQ Toolkit).  Apache 2.0 compatible.
"""

from __future__ import annotations

import time
from typing import Any, Callable, Dict, List, Optional

from loguru import logger
from pydantic import BaseModel, Field


# ═══════════════════════════════════════════════════════════════════════════════
# Tool model
# ═══════════════════════════════════════════════════════════════════════════════


class ImagingTool(BaseModel):
    """Descriptor for a tool available to the agentic reasoning loop."""

    name: str = Field(..., description="Unique tool identifier")
    description: str = Field(..., description="Human-readable tool description")
    function: Callable = Field(..., description="Callable implementing the tool")

    model_config = {"arbitrary_types_allowed": True}


# ═══════════════════════════════════════════════════════════════════════════════
# Tool implementations
# ═══════════════════════════════════════════════════════════════════════════════


async def search_evidence(
    query: str,
    rag_engine: Any,
    *,
    collections_filter: Optional[List[str]] = None,
    top_k_per_collection: int = 5,
    modality_filter: Optional[str] = None,
    body_region_filter: Optional[str] = None,
) -> Dict[str, Any]:
    """Search imaging evidence collections with weighted multi-collection RAG.

    Wraps ``ImagingRAGEngine.retrieve()`` and returns a dict summary
    suitable for the reasoning chain.

    Args:
        query: Natural-language search query.
        rag_engine: An ``ImagingRAGEngine`` instance.
        collections_filter: Optional list of collection names to restrict search.
        top_k_per_collection: Results per collection (default 5).
        modality_filter: Restrict to a single modality (e.g. ``"ct"``).
        body_region_filter: Restrict to a single body region (e.g. ``"head"``).

    Returns:
        Dict with ``hit_count``, ``top_scores``, ``collections_searched``,
        ``search_time_ms``, and the raw ``evidence`` object.
    """
    start = time.time()
    kwargs: Dict[str, Any] = {"top_k_per_collection": top_k_per_collection}
    if collections_filter:
        kwargs["collections_filter"] = collections_filter
    if modality_filter:
        kwargs["modality_filter"] = modality_filter
    if body_region_filter:
        kwargs["body_region_filter"] = body_region_filter

    evidence = rag_engine.retrieve(query, **kwargs)
    elapsed = (time.time() - start) * 1000

    top_scores = [h.score for h in evidence.hits[:5]]
    logger.debug(
        f"search_evidence: {evidence.hit_count} hits, "
        f"top scores {top_scores}, {elapsed:.0f}ms"
    )

    return {
        "hit_count": evidence.hit_count,
        "top_scores": top_scores,
        "collections_searched": evidence.total_collections_searched,
        "search_time_ms": elapsed,
        "knowledge_context": evidence.knowledge_context,
        "evidence": evidence,
    }


async def run_workflow(
    workflow_name: str,
    workflow_registry: Dict[str, Any],
    *,
    input_path: str = "",
    cross_modal_trigger: Any = None,
) -> Dict[str, Any]:
    """Execute a clinical imaging workflow (CT head, lung, CXR, etc.).

    Instantiates the workflow from the registry, runs it, and optionally
    evaluates the cross-modal trigger.

    Args:
        workflow_name: Registry key (e.g. ``"ct_head_hemorrhage"``).
        workflow_registry: Dict mapping workflow names to workflow classes.
        input_path: Path to input volume or image.
        cross_modal_trigger: Optional ``CrossModalTrigger`` instance.

    Returns:
        Dict with ``workflow_name``, ``status``, ``severity``,
        ``classification``, ``cross_modal_triggered``, and raw results.
    """
    if workflow_name not in workflow_registry:
        logger.warning(f"run_workflow: unknown workflow '{workflow_name}'")
        return {
            "workflow_name": workflow_name,
            "status": "not_found",
            "severity": None,
            "classification": None,
            "cross_modal_triggered": False,
            "result": None,
            "cross_modal_result": None,
        }

    start = time.time()
    workflow_cls = workflow_registry[workflow_name]
    workflow = workflow_cls(mock_mode=True, nim_clients={})
    result = workflow.run(input_path)

    cross_modal_result = None
    if result and cross_modal_trigger:
        cross_modal_result = cross_modal_trigger.evaluate(result)

    elapsed = (time.time() - start) * 1000
    logger.debug(
        f"run_workflow: {workflow_name} -> {result.status.value if result else 'None'} "
        f"in {elapsed:.0f}ms"
    )

    return {
        "workflow_name": workflow_name,
        "status": result.status.value if result else "failed",
        "severity": result.severity.value if result else None,
        "classification": result.classification if result else None,
        "cross_modal_triggered": cross_modal_result is not None,
        "result": result,
        "cross_modal_result": cross_modal_result,
        "duration_ms": elapsed,
    }


async def query_genomics(
    query: str,
    cross_modal_trigger: Any,
) -> Dict[str, Any]:
    """Cross-modal query to genomic evidence collection for precision medicine.

    Directly queries the genomic_evidence Milvus collection via the
    ``CrossModalTrigger._query_genomics()`` method.

    Args:
        query: Genomic-domain query string.
        cross_modal_trigger: A ``CrossModalTrigger`` instance.

    Returns:
        Dict with ``hit_count``, ``enrichment_summary``, and raw result.
    """
    start = time.time()
    result = cross_modal_trigger._query_genomics(
        queries=[query],
        trigger_reason=f"Agentic refinement: {query[:80]}",
    )
    elapsed = (time.time() - start) * 1000

    logger.debug(
        f"query_genomics: {result.genomic_hit_count} hits in {elapsed:.0f}ms"
    )

    return {
        "hit_count": result.genomic_hit_count,
        "enrichment_summary": result.enrichment_summary,
        "genomic_context": result.genomic_context,
        "result": result,
        "duration_ms": elapsed,
    }


async def query_agent(
    agent_name: str,
    query_params: Dict[str, Any],
) -> Dict[str, Any]:
    """Query another HCLS AI Factory intelligence agent.

    Dispatches to the appropriate cross-agent query function based on
    ``agent_name``.

    Args:
        agent_name: One of ``"oncology"``, ``"clinical_trial"``,
            ``"cardiology"``, ``"neurology"``.
        query_params: Parameters passed to the cross-agent query function.

    Returns:
        Dict with ``status``, ``agent``, ``recommendations``, and raw data.
    """
    from src.cross_agent import (
        query_cardiology_agent,
        query_neurology_agent,
        query_oncology_agent,
        query_trial_agent,
    )

    dispatch = {
        "oncology": lambda p: query_oncology_agent(
            cancer_type=p.get("cancer_type", ""),
            patient_context=p.get("patient_context"),
        ),
        "clinical_trial": lambda p: query_trial_agent(
            patient_profile=p.get("patient_profile", {}),
        ),
        "cardiology": lambda p: query_cardiology_agent(
            patient_id=p.get("patient_id", ""),
            imaging_context=p.get("imaging_context"),
        ),
        "neurology": lambda p: query_neurology_agent(
            patient_id=p.get("patient_id", ""),
            neurological_context=p.get("neurological_context"),
        ),
    }

    start = time.time()
    func = dispatch.get(agent_name)
    if func is None:
        logger.warning(f"query_agent: unknown agent '{agent_name}'")
        return {"status": "unknown_agent", "agent": agent_name, "recommendations": []}

    result = func(query_params)
    elapsed = (time.time() - start) * 1000

    logger.debug(
        f"query_agent: {agent_name} -> {result.get('status', '?')} in {elapsed:.0f}ms"
    )

    return result


async def extract_radiomics(
    input_path: str,
    mask_path: Optional[str] = None,
    feature_classes: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Extract quantitative radiomics features from segmented imaging volumes.

    Placeholder for PyRadiomics / MONAI radiomics integration.  Returns
    a mock result when the radiomics engine is not available.

    Args:
        input_path: Path to NIfTI volume.
        mask_path: Path to segmentation mask.
        feature_classes: List of feature classes to extract
            (e.g. ``["firstorder", "shape", "glcm"]``).

    Returns:
        Dict with ``feature_count``, ``features``, ``duration_ms``.
    """
    start = time.time()
    classes = feature_classes or ["firstorder", "shape", "glcm"]

    # Placeholder — real implementation would call radiomics_engine
    features: Dict[str, float] = {
        f"{cls}_placeholder": 0.0 for cls in classes
    }

    elapsed = (time.time() - start) * 1000
    logger.debug(f"extract_radiomics: {len(features)} features in {elapsed:.0f}ms")

    return {
        "feature_count": len(features),
        "features": features,
        "feature_classes": classes,
        "duration_ms": elapsed,
    }


async def parse_report(
    report_text: str,
) -> Dict[str, Any]:
    """Parse a radiology report into structured sections, findings, and measurements.

    Uses the ``ReportParser`` from ``src.report_parser`` when available,
    otherwise returns a minimal parsed structure.

    Args:
        report_text: Raw radiology report text.

    Returns:
        Dict with ``sections``, ``findings_count``, ``measurements_count``,
        ``critical_finding``, and raw ``parsed_report``.
    """
    start = time.time()

    try:
        from src.report_parser import ReportParser

        parser = ReportParser()
        parsed = parser.parse(report_text)
        elapsed = (time.time() - start) * 1000

        return {
            "sections": list(parsed.sections.keys()),
            "findings_count": len(parsed.entities),
            "measurements_count": len(parsed.measurements),
            "critical_finding": parsed.critical_finding,
            "parsed_report": parsed,
            "duration_ms": elapsed,
        }
    except ImportError:
        elapsed = (time.time() - start) * 1000
        logger.debug("ReportParser not available, returning minimal parse")
        return {
            "sections": [],
            "findings_count": 0,
            "measurements_count": 0,
            "critical_finding": False,
            "parsed_report": None,
            "duration_ms": elapsed,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# Tool registry
# ═══════════════════════════════════════════════════════════════════════════════


IMAGING_TOOLS: List[ImagingTool] = [
    ImagingTool(
        name="search_imaging_evidence",
        description=(
            "Search imaging evidence collections with weighted "
            "multi-collection RAG"
        ),
        function=search_evidence,
    ),
    ImagingTool(
        name="run_clinical_workflow",
        description=(
            "Execute a clinical imaging workflow "
            "(CT head, lung, CXR, etc.)"
        ),
        function=run_workflow,
    ),
    ImagingTool(
        name="query_genomic_evidence",
        description=(
            "Cross-modal query to genomic evidence collection "
            "for precision medicine"
        ),
        function=query_genomics,
    ),
    ImagingTool(
        name="query_cross_agent",
        description=(
            "Query another HCLS AI Factory intelligence agent "
            "(oncology, cardiology, etc.)"
        ),
        function=query_agent,
    ),
    ImagingTool(
        name="extract_radiomics",
        description=(
            "Extract quantitative radiomics features from "
            "segmented imaging volumes"
        ),
        function=extract_radiomics,
    ),
    ImagingTool(
        name="parse_radiology_report",
        description=(
            "Parse a radiology report into structured sections, "
            "findings, and measurements"
        ),
        function=parse_report,
    ),
]
