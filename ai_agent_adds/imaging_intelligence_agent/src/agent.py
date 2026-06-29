"""Imaging Intelligence Agent — top-level orchestrator.

Implements plan->search->NIM->evaluate->synthesize pipeline.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from loguru import logger

from src.cross_modal import CrossModalTrigger
from src.models import AgentQuery, AgentResponse, CrossModalResult, WorkflowResult
from src.rag_engine import ImagingRAGEngine


@dataclass
class SearchPlan:
    question: str
    identified_topics: List[str] = field(default_factory=list)
    modalities: List[str] = field(default_factory=list)
    body_regions: List[str] = field(default_factory=list)
    search_strategy: str = "broad"  # broad, targeted, comparative
    sub_questions: List[str] = field(default_factory=list)
    recommended_nims: List[str] = field(default_factory=list)


class ImagingIntelligenceAgent:

    MODALITY_KEYWORDS = {
        "ct": ["ct", "computed tomography", "cat scan"],
        "mri": ["mri", "magnetic resonance", "mr imaging"],
        "xray": ["x-ray", "xray", "radiograph"],
        "cxr": ["chest x-ray", "cxr", "chest film"],
        "ultrasound": ["ultrasound", "sonography", "us "],
        "pet": ["pet", "positron emission"],
        "mammography": ["mammography", "mammogram", "breast imaging"],
    }

    REGION_KEYWORDS = {
        "head": ["head", "brain", "cranial", "intracranial"],
        "chest": ["chest", "lung", "pulmonary", "thorax", "thoracic"],
        "abdomen": ["abdomen", "liver", "kidney", "renal", "hepatic", "pancreas"],
        "spine": ["spine", "spinal", "vertebral", "lumbar", "cervical"],
        "cardiac": ["heart", "cardiac", "coronary"],
        "breast": ["breast", "mammary"],
    }

    NIM_KEYWORDS = {
        "vista3d": ["segment", "segmentation", "organ", "anatomy", "vista", "3d"],
        "maisi": ["synthetic", "generate", "augment", "training data"],
        "vilam3": ["interpret", "describe", "report", "what does this show", "vlm"],
    }

    def __init__(
        self,
        rag_engine: ImagingRAGEngine,
        workflow_registry=None,
        nim_manager=None,
        cross_modal_trigger: Optional[CrossModalTrigger] = None,
    ):
        self.rag = rag_engine
        self.workflows = workflow_registry or {}
        self.nim_manager = nim_manager
        self.cross_modal_trigger = cross_modal_trigger

    def search_plan(self, question: str) -> SearchPlan:
        """Analyze question to create a search plan."""
        q = question.lower()
        plan = SearchPlan(question=question)

        # Identify modalities
        for mod, keywords in self.MODALITY_KEYWORDS.items():
            if any(kw in q for kw in keywords):
                plan.modalities.append(mod)

        # Identify body regions
        for region, keywords in self.REGION_KEYWORDS.items():
            if any(kw in q for kw in keywords):
                plan.body_regions.append(region)

        # Determine search strategy
        if self.rag._is_comparative(question):
            plan.search_strategy = "comparative"
        elif plan.modalities or plan.body_regions:
            plan.search_strategy = "targeted"

        # Recommend NIMs
        for nim, keywords in self.NIM_KEYWORDS.items():
            if any(kw in q for kw in keywords):
                plan.recommended_nims.append(nim)

        return plan

    def evaluate_evidence(self, evidence) -> str:
        """Evaluate if evidence is sufficient."""
        if evidence.hit_count >= 10:
            return "sufficient"
        elif evidence.hit_count >= 3:
            return "partial"
        return "insufficient"

    def invoke_workflow(
        self, workflow_name: str, input_path: str = ""
    ) -> tuple[Optional[WorkflowResult], Optional[CrossModalResult]]:
        """Run a reference imaging workflow and evaluate cross-modal triggers.

        Returns:
            Tuple of (WorkflowResult, CrossModalResult).  The cross-modal
            result is None when the trigger is disabled or the workflow
            result does not meet any severity threshold.
        """
        if workflow_name not in self.workflows:
            logger.warning(f"Unknown workflow: {workflow_name}")
            return None, None
        workflow_cls = self.workflows[workflow_name]
        workflow = workflow_cls(mock_mode=True, nim_clients={})
        if self.nim_manager:
            workflow.nim_clients = {
                "vista3d": self.nim_manager.vista3d,
                "vilam3": self.nim_manager.vilam3,
            }
        result = workflow.run(input_path)

        # Evaluate cross-modal trigger
        cross_modal_result = None
        if result and self.cross_modal_trigger:
            cross_modal_result = self.cross_modal_trigger.evaluate(result)
            if cross_modal_result:
                logger.info(
                    f"Cross-modal enrichment attached: "
                    f"{cross_modal_result.genomic_hit_count} genomic hits "
                    f"({cross_modal_result.trigger_reason})"
                )

        return result, cross_modal_result

    def run(self, query: AgentQuery) -> AgentResponse:
        """Full agent pipeline: plan -> search -> NIM -> evaluate -> synthesize."""
        plan = self.search_plan(query.question)
        logger.info(f"Search plan: strategy={plan.search_strategy}, modalities={plan.modalities}")

        # Search
        kwargs = {}
        if plan.modalities and len(plan.modalities) == 1:
            kwargs["modality_filter"] = plan.modalities[0]
        if plan.body_regions and len(plan.body_regions) == 1:
            kwargs["body_region_filter"] = plan.body_regions[0]

        if plan.search_strategy == "comparative":
            comp = self.rag.retrieve_comparative(query.question, **kwargs)
            evidence = comp.evidence_a  # Use combined for answer
        else:
            evidence = self.rag.retrieve(query.question, **kwargs)

        # NIM invocation (if requested and available)
        workflow_results = []
        nim_services_used = []
        cross_modal_result = None

        if query.include_nim and plan.recommended_nims and self.nim_manager:
            for nim_name in plan.recommended_nims:
                try:
                    client = self.nim_manager.get_client(nim_name)
                    if client.is_available() or client.mock_enabled:
                        nim_services_used.append(nim_name)
                except Exception as e:
                    logger.warning(f"NIM {nim_name} invocation skipped: {e}")

        # Evaluate cross-modal triggers on any workflow results
        if workflow_results and self.cross_modal_trigger:
            for wr in workflow_results:
                cm = self.cross_modal_trigger.evaluate(wr)
                if cm is not None:
                    cross_modal_result = cm
                    break  # Use the first trigger match

        # Generate answer
        answer = self.rag.query(query.question)

        return AgentResponse(
            question=query.question,
            answer=answer,
            evidence=evidence,
            workflow_results=workflow_results,
            nim_services_used=nim_services_used,
            knowledge_used=[k for k in ["pathology", "modality", "anatomy"] if k in (evidence.knowledge_context or "").lower()],
            cross_modal=cross_modal_result,
        )

    def generate_report(self, response: AgentResponse) -> str:
        """Generate a markdown report from an agent response."""
        lines = [
            "# Imaging Intelligence Report",
            "",
            f"**Query:** {response.question}",
            f"**Timestamp:** {response.timestamp}",
            "",
            "## Answer",
            "",
            response.answer,
            "",
            "## Evidence Summary",
            "",
            f"- **Total evidence items:** {response.evidence.hit_count}",
            f"- **Collections searched:** {response.evidence.total_collections_searched}",
            f"- **Search time:** {response.evidence.search_time_ms:.0f} ms",
        ]

        if response.nim_services_used:
            lines.extend([
                "",
                "## NIM Services Used",
                "",
            ])
            for nim in response.nim_services_used:
                lines.append(f"- {nim}")

        if response.workflow_results:
            lines.extend([
                "",
                "## Workflow Results",
                "",
            ])
            for wr in response.workflow_results:
                lines.append(f"### {wr.workflow_name}")
                lines.append(f"- Status: {wr.status.value}")
                lines.append(f"- Severity: {wr.severity.value}")
                if wr.measurements:
                    for k, v in wr.measurements.items():
                        lines.append(f"- {k}: {v}")

        if response.cross_modal:
            cm = response.cross_modal
            lines.extend([
                "",
                "## Cross-Modal Genomics Enrichment",
                "",
                f"**Trigger:** {cm.trigger_reason}",
                f"**Genomic hits:** {cm.genomic_hit_count}",
                f"**Queries executed:** {cm.query_count}",
                "",
                cm.enrichment_summary,
            ])
            if cm.genomic_context:
                lines.extend([
                    "",
                    "### Genomic Evidence",
                    "",
                ])
                for ctx in cm.genomic_context[:10]:  # Limit to top 10
                    lines.append(f"- {ctx}")

        lines.extend([
            "",
            "---",
            "*Research use only. All findings require clinician review.*",
        ])

        return "\n".join(lines)
