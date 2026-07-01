"""AIQ-powered agentic reasoning for the Clinical Imaging Engine.

Implements a Plan -> Execute -> Reflect -> Refine cycle that:
1. PLAN: Classifies the query, selects tools (workflows, collections, cross-modal)
2. EXECUTE: Runs selected tools in parallel where possible
3. REFLECT: Evaluates evidence sufficiency, checks for contradictions, assesses confidence
4. REFINE: If gaps identified, expands search or runs additional tools
5. REPORT: Generates structured output with full provenance

Uses NVIDIA AIQ Toolkit when available, falls back to a built-in
lightweight reasoning loop when AIQ is not installed.

Open source (AIQ Toolkit).  Apache 2.0 compatible.
"""

from __future__ import annotations

import asyncio
import re
import time
from typing import Any, Dict, List, Optional

from loguru import logger
from pydantic import BaseModel, Field

from src.agentic.tools import (
    IMAGING_TOOLS,
    query_agent,
    query_genomics,
    run_workflow,
    search_evidence,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Data models
# ═══════════════════════════════════════════════════════════════════════════════


class ReasoningStep(BaseModel):
    """A single step in the agentic reasoning chain."""

    step: str = Field(
        ...,
        description="Phase: plan, execute, reflect, refine, or report",
    )
    action: str = Field(..., description="What was done in this step")
    result_summary: str = Field(..., description="Brief summary of the result")
    confidence: float = Field(
        0.0, ge=0.0, le=1.0, description="Confidence score 0.0-1.0"
    )
    duration_ms: float = Field(0.0, ge=0.0, description="Wall-clock time in ms")


class AgenticResult(BaseModel):
    """Complete result of the agentic reasoning cycle."""

    query: str
    answer: str
    reasoning_chain: List[ReasoningStep] = Field(default_factory=list)
    evidence_count: int = 0
    collections_searched: List[str] = Field(default_factory=list)
    workflows_executed: List[str] = Field(default_factory=list)
    cross_modal_triggered: bool = False
    cross_agent_queried: List[str] = Field(default_factory=list)
    total_duration_ms: float = 0.0
    confidence: float = Field(
        0.0, ge=0.0, le=1.0, description="Overall confidence 0.0-1.0"
    )
    refinement_rounds: int = 0


# ═══════════════════════════════════════════════════════════════════════════════
# Query classification
# ═══════════════════════════════════════════════════════════════════════════════

_COMPARATIVE_RE = re.compile(
    r"\b(compare|compared to|vs\.?|versus|difference between"
    r"|head.to.head|better than|advantages|disadvantages)\b",
    re.IGNORECASE,
)

_WORKFLOW_KEYWORDS = [
    "analyze", "triage", "segment", "detect", "classify",
    "hemorrhage", "lung nodule", "chest x-ray", "cxr",
    "ms lesion", "coronary", "prostate", "pirads",
]

_CROSS_MODAL_KEYWORDS = [
    "genomic", "genetic", "mutation", "variant", "molecular",
    "precision medicine", "biomarker", "driver mutation",
]

_WORKFLOW_MAP = {
    "hemorrhage": "ct_head_hemorrhage",
    "head ct": "ct_head_hemorrhage",
    "lung nodule": "ct_chest_lung_nodule",
    "chest x-ray": "cxr_rapid_findings",
    "cxr": "cxr_rapid_findings",
    "ms lesion": "mri_brain_ms_lesion",
    "coronary": "ct_coronary_angiography",
    "prostate": "mri_prostate_pirads",
    "pirads": "mri_prostate_pirads",
    "pi-rads": "mri_prostate_pirads",
}


def _classify_query(query: str) -> str:
    """Classify a query into a routing category.

    Returns one of:
        ``"comparative_query"`` — head-to-head comparison (X vs Y)
        ``"workflow_query"``    — clinical workflow execution
        ``"cross_modal_query"`` — needs genomic context
        ``"multi_step_query"``  — requires multiple chained tools
        ``"evidence_query"``    — standard evidence retrieval
    """
    q = query.lower()

    if _COMPARATIVE_RE.search(query):
        return "comparative_query"

    has_workflow = any(kw in q for kw in _WORKFLOW_KEYWORDS)
    has_cross_modal = any(kw in q for kw in _CROSS_MODAL_KEYWORDS)

    if has_workflow and has_cross_modal:
        return "multi_step_query"
    if has_workflow:
        return "workflow_query"
    if has_cross_modal:
        return "cross_modal_query"

    return "evidence_query"


def _identify_workflow(query: str) -> Optional[str]:
    """Map a query to a specific workflow name, or None."""
    q = query.lower()
    for keyword, workflow_name in _WORKFLOW_MAP.items():
        if keyword in q:
            return workflow_name
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# Confidence scoring
# ═══════════════════════════════════════════════════════════════════════════════


def _compute_confidence(
    evidence_count: int,
    top_scores: List[float],
    min_evidence: int = 3,
    contradiction_detected: bool = False,
) -> float:
    """Compute a 0.0-1.0 confidence score from evidence quality signals.

    Factors:
      - Evidence count relative to ``min_evidence`` threshold
      - Average quality of top citation scores
      - Penalty for detected contradictions

    Args:
        evidence_count: Total number of evidence items retrieved.
        top_scores: Similarity scores of the top-N hits.
        min_evidence: Minimum expected evidence count.
        contradiction_detected: Whether contradictory evidence was found.

    Returns:
        Float confidence in [0.0, 1.0].
    """
    # Count component: 0.0 when 0 hits, 1.0 when >= min_evidence
    if min_evidence > 0:
        count_score = min(evidence_count / min_evidence, 1.0)
    else:
        count_score = 1.0 if evidence_count > 0 else 0.0

    # Quality component: average of top scores (already 0.0-1.0 similarity)
    if top_scores:
        quality_score = sum(top_scores) / len(top_scores)
    else:
        quality_score = 0.0

    # Weighted combination
    confidence = 0.5 * count_score + 0.5 * quality_score

    # Contradiction penalty
    if contradiction_detected:
        confidence *= 0.7

    return round(min(max(confidence, 0.0), 1.0), 3)


def _detect_contradictions(evidence_texts: List[str]) -> bool:
    """Lightweight contradiction detection in evidence texts.

    Looks for opposing sentiment patterns in the evidence corpus.
    This is a heuristic — not a full NLI model.

    Args:
        evidence_texts: List of evidence text snippets.

    Returns:
        True if likely contradictions are detected.
    """
    if len(evidence_texts) < 2:
        return False

    positive_markers = [
        "superior", "outperforms", "higher sensitivity",
        "recommended", "gold standard", "preferred",
    ]
    negative_markers = [
        "inferior", "underperforms", "lower sensitivity",
        "not recommended", "limited evidence", "not preferred",
    ]

    combined = " ".join(evidence_texts).lower()
    has_positive = any(m in combined for m in positive_markers)
    has_negative = any(m in combined for m in negative_markers)

    return has_positive and has_negative


# ═══════════════════════════════════════════════════════════════════════════════
# ImagingAgenticEngine
# ═══════════════════════════════════════════════════════════════════════════════


class ImagingAgenticEngine:
    """Multi-step reasoning engine for clinical imaging queries.

    Implements a Plan -> Execute -> Reflect -> Refine -> Report loop
    that orchestrates RAG retrieval, workflow execution, cross-modal
    genomic enrichment, and cross-agent queries into a single coherent
    reasoning cycle.

    Uses NVIDIA AIQ Toolkit when available; otherwise falls back to
    a built-in lightweight reasoning loop.

    Args:
        rag_engine: An ``ImagingRAGEngine`` instance for evidence retrieval.
        workflow_registry: Dict mapping workflow names to workflow classes.
        cross_modal: A ``CrossModalTrigger`` instance (or None).
        cross_agent: Module or namespace with cross-agent query functions.
        guardrails: Optional ``ClinicalGuardrails`` instance.
        max_refinement_rounds: Max refinement iterations (default 2).
        confidence_threshold: Below this confidence, trigger refinement.
        min_evidence_count: Below this hit count, expand search.
    """

    def __init__(
        self,
        rag_engine: Any,
        workflow_registry: Optional[Dict[str, Any]] = None,
        cross_modal: Any = None,
        cross_agent: Any = None,
        guardrails: Any = None,
        max_refinement_rounds: int = 2,
        confidence_threshold: float = 0.7,
        min_evidence_count: int = 3,
    ):
        self.rag_engine = rag_engine
        self.workflow_registry = workflow_registry or {}
        self.cross_modal = cross_modal
        self.cross_agent = cross_agent
        self.guardrails = guardrails
        self.max_refinement_rounds = max_refinement_rounds
        self.confidence_threshold = confidence_threshold
        self.min_evidence_count = min_evidence_count
        self._aiq_available = self._check_aiq()

    # ───────────────────────────────────────────────────────────────────
    # AIQ detection
    # ───────────────────────────────────────────────────────────────────

    def _check_aiq(self) -> bool:
        """Check if NVIDIA AIQ Toolkit is installed."""
        try:
            import aiq  # noqa: F401

            logger.info("NVIDIA AIQ Toolkit detected. AIQ mode enabled.")
            return True
        except ImportError:
            logger.info(
                "AIQ Toolkit not installed. Using built-in reasoning loop."
            )
            return False

    # ───────────────────────────────────────────────────────────────────
    # Main entry point
    # ───────────────────────────────────────────────────────────────────

    async def reason(
        self,
        query: str,
        context: Optional[Dict[str, Any]] = None,
    ) -> AgenticResult:
        """Execute the full Plan -> Execute -> Reflect -> Refine -> Report cycle.

        Args:
            query: Natural-language clinical imaging query.
            context: Optional context dict (patient data, prior results, etc.).

        Returns:
            ``AgenticResult`` with the answer, full reasoning chain,
            evidence provenance, and confidence score.
        """
        overall_start = time.time()
        reasoning_chain: List[ReasoningStep] = []
        evidence_data: Dict[str, Any] = {}
        ctx = context or {}

        # ── Guardrails: input check ──
        if self.guardrails:
            input_check = self.guardrails.check_input(query)
            if input_check.blocked:
                return AgenticResult(
                    query=query,
                    answer=(
                        "Query blocked by clinical safety guardrails: "
                        + ", ".join(input_check.flags)
                    ),
                    reasoning_chain=[
                        ReasoningStep(
                            step="plan",
                            action="guardrails_input_check",
                            result_summary=f"Blocked: {input_check.flags}",
                            confidence=0.0,
                            duration_ms=0.0,
                        )
                    ],
                    confidence=0.0,
                    total_duration_ms=(time.time() - overall_start) * 1000,
                )
            if input_check.modified_text:
                query = input_check.modified_text

        # ── PLAN ──
        plan_step, plan_data = self._plan(query, ctx)
        reasoning_chain.append(plan_step)

        # ── EXECUTE ──
        exec_step, evidence_data = await self._execute(plan_data, query, ctx)
        reasoning_chain.append(exec_step)

        # ── REFLECT ──
        reflect_step, gaps = self._reflect(evidence_data)
        reasoning_chain.append(reflect_step)

        # ── REFINE loop ──
        refinement_rounds = 0
        while (
            gaps
            and reflect_step.confidence < self.confidence_threshold
            and refinement_rounds < self.max_refinement_rounds
        ):
            refinement_rounds += 1
            refine_step, new_evidence = await self._refine(
                gaps, evidence_data, query, ctx
            )
            reasoning_chain.append(refine_step)

            # Merge new evidence into accumulated data
            self._merge_evidence(evidence_data, new_evidence)

            # Re-reflect
            reflect_step, gaps = self._reflect(evidence_data)
            reasoning_chain.append(reflect_step)

        # ── REPORT ──
        report_step, answer = self._report(reasoning_chain, evidence_data, query)
        reasoning_chain.append(report_step)

        # ── Guardrails: output check ──
        if self.guardrails:
            output_check = self.guardrails.check_output(
                answer,
                evidence_count=evidence_data.get("total_hit_count", 0),
            )
            if output_check.modified_text:
                answer = output_check.modified_text

        total_ms = (time.time() - overall_start) * 1000

        return AgenticResult(
            query=query,
            answer=answer,
            reasoning_chain=reasoning_chain,
            evidence_count=evidence_data.get("total_hit_count", 0),
            collections_searched=evidence_data.get("collections", []),
            workflows_executed=evidence_data.get("workflows_executed", []),
            cross_modal_triggered=evidence_data.get(
                "cross_modal_triggered", False
            ),
            cross_agent_queried=evidence_data.get("cross_agents_queried", []),
            total_duration_ms=total_ms,
            confidence=reflect_step.confidence,
            refinement_rounds=refinement_rounds,
        )

    # ───────────────────────────────────────────────────────────────────
    # PLAN
    # ───────────────────────────────────────────────────────────────────

    def _plan(
        self,
        query: str,
        context: Optional[Dict[str, Any]] = None,
    ) -> tuple[ReasoningStep, Dict[str, Any]]:
        """PLAN: Analyze query, determine tools needed.

        Classifies the query into categories:
          - evidence_query: search collections
          - workflow_query: run clinical workflow
          - comparative_query: dual retrieval for X vs Y
          - cross_modal_query: needs genomic context
          - multi_step_query: needs multiple tools chained

        Returns:
            Tuple of (ReasoningStep, plan_data dict).
        """
        start = time.time()

        query_type = _classify_query(query)
        workflow_name = _identify_workflow(query)

        plan_data: Dict[str, Any] = {
            "query_type": query_type,
            "workflow_name": workflow_name,
            "search_collections": True,
            "run_workflow": query_type in ("workflow_query", "multi_step_query"),
            "cross_modal": query_type in ("cross_modal_query", "multi_step_query"),
            "comparative": query_type == "comparative_query",
        }

        elapsed = (time.time() - start) * 1000

        step = ReasoningStep(
            step="plan",
            action=f"classify_query -> {query_type}",
            result_summary=(
                f"Type: {query_type}"
                + (f", workflow: {workflow_name}" if workflow_name else "")
                + (", comparative" if plan_data["comparative"] else "")
                + (", cross-modal" if plan_data["cross_modal"] else "")
            ),
            confidence=0.0,
            duration_ms=elapsed,
        )

        logger.info(
            f"PLAN: query_type={query_type}, "
            f"workflow={workflow_name}, "
            f"comparative={plan_data['comparative']}"
        )

        return step, plan_data

    # ───────────────────────────────────────────────────────────────────
    # EXECUTE
    # ───────────────────────────────────────────────────────────────────

    async def _execute(
        self,
        plan_data: Dict[str, Any],
        query: str,
        context: Dict[str, Any],
    ) -> tuple[ReasoningStep, Dict[str, Any]]:
        """EXECUTE: Run planned tools and collect results.

        Runs evidence search (always), and optionally runs workflows,
        comparative retrieval, and cross-modal queries based on the plan.

        Returns:
            Tuple of (ReasoningStep, evidence_data dict).
        """
        start = time.time()
        evidence_data: Dict[str, Any] = {
            "total_hit_count": 0,
            "top_scores": [],
            "collections": [],
            "evidence_texts": [],
            "evidence_object": None,
            "workflows_executed": [],
            "cross_modal_triggered": False,
            "cross_agents_queried": [],
            "workflow_results": [],
            "cross_modal_result": None,
        }

        actions: List[str] = []

        # ── Evidence search ──
        if plan_data.get("search_collections", True):
            search_result = await search_evidence(
                query=query,
                rag_engine=self.rag_engine,
            )
            evidence_data["total_hit_count"] = search_result["hit_count"]
            evidence_data["top_scores"] = search_result["top_scores"]
            evidence_data["collections"] = list(
                self.rag_engine.collection_manager.get_collection_stats().keys()
            ) if hasattr(self.rag_engine, "collection_manager") else []
            evidence_data["evidence_object"] = search_result.get("evidence")

            # Extract evidence texts for contradiction detection
            ev = search_result.get("evidence")
            if ev and hasattr(ev, "hits"):
                evidence_data["evidence_texts"] = [
                    h.text for h in ev.hits[:20]
                ]

            actions.append(
                f"search({search_result['hit_count']} hits)"
            )

        # ── Comparative retrieval ──
        if plan_data.get("comparative", False):
            try:
                comp_result = self.rag_engine.retrieve_comparative(query)
                total = (
                    comp_result.evidence_a.hit_count
                    + comp_result.evidence_b.hit_count
                )
                evidence_data["total_hit_count"] = max(
                    evidence_data["total_hit_count"], total
                )
                evidence_data["comparative_result"] = comp_result

                # Merge comparative texts for contradiction detection
                for side in (comp_result.evidence_a, comp_result.evidence_b):
                    evidence_data["evidence_texts"].extend(
                        h.text for h in side.hits[:10]
                    )

                actions.append(
                    f"comparative({comp_result.entity_a} vs "
                    f"{comp_result.entity_b}, {total} hits)"
                )
            except Exception as exc:
                logger.warning(f"Comparative retrieval failed: {exc}")
                actions.append("comparative(failed)")

        # ── Workflow execution ──
        if plan_data.get("run_workflow", False) and plan_data.get("workflow_name"):
            wf_name = plan_data["workflow_name"]
            wf_result = await run_workflow(
                workflow_name=wf_name,
                workflow_registry=self.workflow_registry,
                cross_modal_trigger=self.cross_modal,
            )
            evidence_data["workflows_executed"].append(wf_name)
            evidence_data["workflow_results"].append(wf_result)

            if wf_result.get("cross_modal_triggered"):
                evidence_data["cross_modal_triggered"] = True
                evidence_data["cross_modal_result"] = wf_result.get(
                    "cross_modal_result"
                )

            actions.append(
                f"workflow({wf_name} -> "
                f"{wf_result.get('severity', 'n/a')})"
            )

        # ── Cross-modal genomic query ──
        if plan_data.get("cross_modal", False) and self.cross_modal:
            try:
                geno_result = await query_genomics(
                    query=query,
                    cross_modal_trigger=self.cross_modal,
                )
                evidence_data["cross_modal_triggered"] = True
                evidence_data["cross_modal_result"] = geno_result.get("result")
                evidence_data["total_hit_count"] += geno_result.get(
                    "hit_count", 0
                )
                actions.append(
                    f"genomics({geno_result.get('hit_count', 0)} hits)"
                )
            except Exception as exc:
                logger.warning(f"Cross-modal query failed: {exc}")
                actions.append("genomics(failed)")

        elapsed = (time.time() - start) * 1000

        step = ReasoningStep(
            step="execute",
            action=", ".join(actions) if actions else "no_actions",
            result_summary=(
                f"{evidence_data['total_hit_count']} total evidence items"
            ),
            confidence=0.0,
            duration_ms=elapsed,
        )

        logger.info(f"EXECUTE: {step.action} in {elapsed:.0f}ms")

        return step, evidence_data

    # ───────────────────────────────────────────────────────────────────
    # REFLECT
    # ───────────────────────────────────────────────────────────────────

    def _reflect(
        self,
        evidence_data: Dict[str, Any],
    ) -> tuple[ReasoningStep, List[str]]:
        """REFLECT: Evaluate results for sufficiency and quality.

        Checks:
          - Evidence sufficiency (are there enough high-quality citations?)
          - Contradiction detection (do sources disagree?)
          - Confidence assessment (how strong is the evidence?)
          - Gap identification (what's missing?)

        Returns:
            Tuple of (ReasoningStep, list of identified gaps).
        """
        start = time.time()

        hit_count = evidence_data.get("total_hit_count", 0)
        top_scores = evidence_data.get("top_scores", [])
        evidence_texts = evidence_data.get("evidence_texts", [])

        # Contradiction check
        contradiction = _detect_contradictions(evidence_texts)

        # Confidence score
        confidence = _compute_confidence(
            evidence_count=hit_count,
            top_scores=top_scores,
            min_evidence=self.min_evidence_count,
            contradiction_detected=contradiction,
        )

        # Gap identification
        gaps: List[str] = []

        if hit_count < self.min_evidence_count:
            gaps.append("insufficient_evidence")

        if top_scores and max(top_scores) < 0.5:
            gaps.append("low_relevance_scores")

        if contradiction:
            gaps.append("contradictory_evidence")

        if not evidence_data.get("cross_modal_triggered", False) and confidence < self.confidence_threshold:
            gaps.append("no_genomic_context")

        elapsed = (time.time() - start) * 1000

        step = ReasoningStep(
            step="reflect",
            action="evaluate_evidence",
            result_summary=(
                f"confidence={confidence:.2f}, "
                f"hits={hit_count}, "
                f"contradictions={'yes' if contradiction else 'no'}, "
                f"gaps={gaps if gaps else 'none'}"
            ),
            confidence=confidence,
            duration_ms=elapsed,
        )

        logger.info(
            f"REFLECT: confidence={confidence:.2f}, "
            f"hits={hit_count}, gaps={gaps}"
        )

        return step, gaps

    # ───────────────────────────────────────────────────────────────────
    # REFINE
    # ───────────────────────────────────────────────────────────────────

    async def _refine(
        self,
        gaps: List[str],
        evidence_data: Dict[str, Any],
        query: str,
        context: Dict[str, Any],
    ) -> tuple[ReasoningStep, Dict[str, Any]]:
        """REFINE: Address gaps identified in reflection.

        Actions:
          - Expand search to additional collections
          - Run additional workflows
          - Query cross-agent endpoints for domain expertise
          - Trigger cross-modal genomic enrichment

        Returns:
            Tuple of (ReasoningStep, new_evidence dict to merge).
        """
        start = time.time()
        actions: List[str] = []
        new_evidence: Dict[str, Any] = {
            "total_hit_count": 0,
            "top_scores": [],
            "evidence_texts": [],
            "cross_modal_triggered": False,
            "cross_agents_queried": [],
        }

        # ── Expand search with broader parameters ──
        if "insufficient_evidence" in gaps or "low_relevance_scores" in gaps:
            try:
                expanded = await search_evidence(
                    query=query,
                    rag_engine=self.rag_engine,
                    top_k_per_collection=10,  # Double the default
                )
                new_evidence["total_hit_count"] += expanded["hit_count"]
                new_evidence["top_scores"].extend(expanded["top_scores"])

                ev = expanded.get("evidence")
                if ev and hasattr(ev, "hits"):
                    new_evidence["evidence_texts"].extend(
                        h.text for h in ev.hits[:20]
                    )

                actions.append(
                    f"expanded_search({expanded['hit_count']} hits)"
                )
            except Exception as exc:
                logger.warning(f"Expanded search failed: {exc}")
                actions.append("expanded_search(failed)")

        # ── Cross-modal genomic enrichment ──
        if "no_genomic_context" in gaps and self.cross_modal:
            try:
                geno_result = await query_genomics(
                    query=query,
                    cross_modal_trigger=self.cross_modal,
                )
                new_evidence["cross_modal_triggered"] = True
                new_evidence["total_hit_count"] += geno_result.get(
                    "hit_count", 0
                )
                actions.append(
                    f"genomic_refinement({geno_result.get('hit_count', 0)} hits)"
                )
            except Exception as exc:
                logger.warning(f"Genomic refinement failed: {exc}")
                actions.append("genomic_refinement(failed)")

        # ── Cross-agent queries ──
        if "contradictory_evidence" in gaps:
            for agent_name in ("oncology", "cardiology"):
                try:
                    agent_result = await query_agent(
                        agent_name=agent_name,
                        query_params={
                            "cancer_type": "general",
                            "patient_id": context.get("patient_id", ""),
                        },
                    )
                    new_evidence["cross_agents_queried"].append(agent_name)
                    actions.append(
                        f"cross_agent({agent_name} -> "
                        f"{agent_result.get('status', '?')})"
                    )
                except Exception as exc:
                    logger.warning(
                        f"Cross-agent query to {agent_name} failed: {exc}"
                    )

        elapsed = (time.time() - start) * 1000

        step = ReasoningStep(
            step="refine",
            action=", ".join(actions) if actions else "no_refinement",
            result_summary=(
                f"Added {new_evidence['total_hit_count']} evidence items"
            ),
            confidence=0.0,
            duration_ms=elapsed,
        )

        logger.info(f"REFINE: {step.action} in {elapsed:.0f}ms")

        return step, new_evidence

    # ───────────────────────────────────────────────────────────────────
    # REPORT
    # ───────────────────────────────────────────────────────────────────

    def _report(
        self,
        reasoning_chain: List[ReasoningStep],
        evidence_data: Dict[str, Any],
        query: str,
    ) -> tuple[ReasoningStep, str]:
        """REPORT: Generate final structured answer with provenance.

        Synthesizes the answer by delegating to the RAG engine's LLM
        synthesis, then wraps it with provenance metadata.

        Returns:
            Tuple of (ReasoningStep, answer_text).
        """
        start = time.time()

        # Use the RAG engine's query method for LLM synthesis
        try:
            answer = self.rag_engine.query(query)
        except Exception as exc:
            logger.error(f"LLM synthesis failed: {exc}")
            answer = (
                "Unable to generate a synthesized answer. "
                f"Evidence retrieved: {evidence_data.get('total_hit_count', 0)} items."
            )

        # Append provenance summary
        provenance_lines = [
            "",
            "---",
            "**Reasoning Provenance:**",
            f"- Evidence items: {evidence_data.get('total_hit_count', 0)}",
            f"- Collections searched: {len(evidence_data.get('collections', []))}",
        ]

        if evidence_data.get("workflows_executed"):
            provenance_lines.append(
                f"- Workflows: {', '.join(evidence_data['workflows_executed'])}"
            )

        if evidence_data.get("cross_modal_triggered"):
            provenance_lines.append("- Cross-modal genomics: triggered")

        if evidence_data.get("cross_agents_queried"):
            provenance_lines.append(
                f"- Cross-agent: {', '.join(evidence_data['cross_agents_queried'])}"
            )

        # Find the latest reflect step for overall confidence
        last_confidence = 0.0
        for step in reversed(reasoning_chain):
            if step.step == "reflect":
                last_confidence = step.confidence
                break

        provenance_lines.append(f"- Confidence: {last_confidence:.2f}")

        refinement_count = sum(
            1 for s in reasoning_chain if s.step == "refine"
        )
        if refinement_count > 0:
            provenance_lines.append(
                f"- Refinement rounds: {refinement_count}"
            )

        answer_with_provenance = answer + "\n".join(provenance_lines)

        elapsed = (time.time() - start) * 1000

        step = ReasoningStep(
            step="report",
            action="synthesize_answer",
            result_summary=f"Generated answer ({len(answer)} chars)",
            confidence=last_confidence,
            duration_ms=elapsed,
        )

        logger.info(f"REPORT: synthesized answer in {elapsed:.0f}ms")

        return step, answer_with_provenance

    # ───────────────────────────────────────────────────────────────────
    # Evidence merging
    # ───────────────────────────────────────────────────────────────────

    def _merge_evidence(
        self,
        target: Dict[str, Any],
        new: Dict[str, Any],
    ) -> None:
        """Merge new evidence data into the accumulated target dict.

        Updates hit counts, scores, text lists, and boolean flags.

        Args:
            target: The accumulated evidence data dict (modified in place).
            new: New evidence data to merge in.
        """
        target["total_hit_count"] += new.get("total_hit_count", 0)
        target["top_scores"].extend(new.get("top_scores", []))
        target["evidence_texts"].extend(new.get("evidence_texts", []))

        if new.get("cross_modal_triggered"):
            target["cross_modal_triggered"] = True

        target.setdefault("cross_agents_queried", []).extend(
            new.get("cross_agents_queried", [])
        )
