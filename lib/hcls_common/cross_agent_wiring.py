"""
Cross-Agent Event Wiring for the HCLS AI Factory.

Defines the standard event subscriptions and handlers that enable
intelligence flow between agents. Each agent registers its handlers
at startup to participate in the cross-agent intelligence fabric.

Event Flow Architecture:
    Genomics -> VCF_READY -> RAG/Chat -> TARGETS_IDENTIFIED -> Drug Discovery
                                   |
                            VARIANT_CLASSIFIED
                                   |
    Biomarker Agent <- biomarkers + genotypes -> PGX_RESULT_READY
         |                                          |
    BIOMARKER_ALERT                          PGX_DRUG_FILTER
         |                                          |
    Oncology Agent                          Drug Discovery (re-rank)
         |
    THERAPY_RANKED -> CAR-T Agent
                  -> Imaging Agent (CONCORDANCE_CHECK)
                  -> Autoimmune Agent (if autoimmune pathway)
"""

from __future__ import annotations

from typing import Any, Callable

from loguru import logger


# Handler registry type
HandlerSpec = dict[str, Any]  # {"event_types": [...], "handler": callable, "name": str}


def biomarker_agent_handlers() -> list[HandlerSpec]:
    """Event handlers that the Biomarker Agent should register.

    Listens for:
    - VCF_READY: Extract genotypes for PGx analysis
    - IMAGING_FINDING: Correlate imaging findings with biomarker data
    - ONCOLOGY_CASE_CREATED: Provide biomarker context for oncology cases
    """
    return [
        {
            "event_types": ["vcf_ready"],
            "name": "biomarker_genotype_extractor",
            "description": "Extracts PGx-relevant genotypes from VCF for pharmacogenomic analysis",
        },
        {
            "event_types": ["imaging_finding"],
            "name": "biomarker_imaging_correlator",
            "description": "Correlates imaging findings with biomarker trajectories for concordance scoring",
        },
        {
            "event_types": ["oncology_case_created"],
            "name": "biomarker_oncology_context",
            "description": "Provides inflammatory and metabolic biomarker context to oncology cases",
        },
    ]


def oncology_agent_handlers() -> list[HandlerSpec]:
    """Event handlers that the Oncology Agent should register.

    Listens for:
    - BIOMARKER_ALERT: Receive biomarker alerts for tumor monitoring
    - PGX_RESULT_READY: Filter therapy recommendations by PGx
    - IMAGING_FINDING: Incorporate imaging into tumor assessment
    - TARGETS_IDENTIFIED: Receive genomic targets for molecular profiling
    """
    return [
        {
            "event_types": ["biomarker_alert"],
            "name": "oncology_biomarker_monitor",
            "description": "Monitors biomarker alerts for tumor-related signals",
        },
        {
            "event_types": ["pgx_result_ready"],
            "name": "oncology_pgx_therapy_filter",
            "description": "Filters therapy recommendations based on PGx metabolizer status",
        },
        {
            "event_types": ["imaging_finding"],
            "name": "oncology_imaging_correlator",
            "description": "Incorporates imaging findings into tumor staging and response assessment",
        },
        {
            "event_types": ["targets_identified"],
            "name": "oncology_target_profiler",
            "description": "Receives genomic targets for molecular tumor profiling",
        },
    ]


def imaging_agent_handlers() -> list[HandlerSpec]:
    """Event handlers that the Imaging Agent should register.

    Listens for:
    - BIOMARKER_ALERT: Correlate biomarker alerts with imaging findings
    - ONCOLOGY_CASE_CREATED: Provide imaging context for oncology cases
    - DISEASE_TRAJECTORY_ALERT: Correlate disease trajectory with imaging progression
    """
    return [
        {
            "event_types": ["biomarker_alert"],
            "name": "imaging_biomarker_correlator",
            "description": "Correlates biomarker alerts with imaging findings for concordance",
        },
        {
            "event_types": ["oncology_case_created"],
            "name": "imaging_oncology_context",
            "description": "Provides imaging context for oncology tumor board cases",
        },
        {
            "event_types": ["disease_trajectory_alert"],
            "name": "imaging_disease_progression",
            "description": "Correlates disease trajectory predictions with imaging progression",
        },
    ]


def cart_agent_handlers() -> list[HandlerSpec]:
    """Event handlers that the CAR-T Agent should register.

    Listens for:
    - THERAPY_RANKED: Evaluate CAR-T therapy ranking among options
    - BIOMARKER_ALERT: Assess T-cell fitness from blood biomarkers
    - IMAGING_FINDING: Evaluate tumor burden from imaging
    """
    return [
        {
            "event_types": ["therapy_ranked"],
            "name": "cart_therapy_evaluator",
            "description": "Evaluates CAR-T as a therapy option when ranked by oncology agent",
        },
        {
            "event_types": ["biomarker_alert"],
            "name": "cart_tcell_fitness",
            "description": "Assesses T-cell fitness from biomarker data for manufacturing feasibility",
        },
        {
            "event_types": ["imaging_finding"],
            "name": "cart_tumor_burden",
            "description": "Evaluates tumor burden from imaging for CAR-T eligibility",
        },
    ]


def autoimmune_agent_handlers() -> list[HandlerSpec]:
    """Event handlers that the Autoimmune Agent should register.

    Listens for:
    - BIOMARKER_ALERT: Monitor inflammatory markers for flare prediction
    - PGX_RESULT_READY: Filter biologic therapy by PGx profile
    - IMAGING_FINDING: Correlate joint/organ imaging with disease activity
    """
    return [
        {
            "event_types": ["biomarker_alert"],
            "name": "autoimmune_inflammation_monitor",
            "description": "Monitors inflammatory biomarker alerts for flare prediction",
        },
        {
            "event_types": ["pgx_result_ready"],
            "name": "autoimmune_biologic_filter",
            "description": "Filters biologic therapy recommendations based on PGx profile",
        },
        {
            "event_types": ["imaging_finding"],
            "name": "autoimmune_imaging_correlator",
            "description": "Correlates joint/organ imaging with disease activity scores",
        },
    ]


def drug_discovery_handlers() -> list[HandlerSpec]:
    """Event handlers for the Drug Discovery Pipeline.

    Listens for:
    - PGX_DRUG_FILTER: Re-rank candidates based on PGx profile
    - TARGETS_IDENTIFIED: Start drug discovery for new targets
    """
    return [
        {
            "event_types": ["pgx_drug_filter"],
            "name": "drug_discovery_pgx_rerank",
            "description": "Re-ranks drug candidates based on patient PGx metabolizer profile",
        },
        {
            "event_types": ["targets_identified"],
            "name": "drug_discovery_target_receiver",
            "description": "Initiates drug discovery pipeline for newly identified targets",
        },
    ]


def get_all_handler_specs() -> dict[str, list[HandlerSpec]]:
    """Return all handler specifications grouped by agent.

    This is the master registry of cross-agent event wiring.
    Used by the orchestrator to verify all handlers are registered.
    """
    return {
        "precision_biomarker_agent": biomarker_agent_handlers(),
        "precision_oncology_agent": oncology_agent_handlers(),
        "imaging_intelligence_agent": imaging_agent_handlers(),
        "cart_intelligence_agent": cart_agent_handlers(),
        "precision_autoimmune_agent": autoimmune_agent_handlers(),
        "drug_discovery_pipeline": drug_discovery_handlers(),
    }


def print_wiring_diagram() -> str:
    """Generate a text-based wiring diagram of cross-agent events.

    Returns a human-readable diagram showing how events flow
    between agents in the HCLS AI Factory.
    """
    diagram = """
+========================================================================+
|              HCLS AI Factory -- Cross-Agent Intelligence Fabric         |
+========================================================================+

  +--------------+    VCF_READY     +------------------+
  |  Genomics    |----------------->|    RAG / Chat    |
  |  Pipeline    |                  |    Pipeline      |
  +--------------+                  +--------+---------+
                                             | TARGETS_IDENTIFIED
                                             v
  +--------------+  BIOMARKER_ALERT  +------------------+
  |  Precision   |<-----------------|  Precision        |
  |  Biomarker   |  PGX_RESULT_READY|  Oncology        |
  |  Agent       |----------------->|  Agent            |
  +------+-------+                  +--------+---------+
         |                                   |
         | PGX_DRUG_FILTER                   | THERAPY_RANKED
         v                                   v
  +--------------+                  +------------------+
  |    Drug      |                  |    CAR-T         |
  |  Discovery   |                  |  Intelligence    |
  |  Pipeline    |                  |    Agent         |
  +--------------+                  +------------------+
         ^                                   ^
         |                                   |
  +------+-------+  IMAGING_FINDING  +-------+----------+
  |  Imaging     |----------------->|  Precision        |
  | Intelligence |  CONCORDANCE     |  Autoimmune       |
  |    Agent     |<-----------------|    Agent          |
  +--------------+                  +------------------+
"""
    return diagram
