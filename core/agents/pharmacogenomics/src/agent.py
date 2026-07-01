"""Pharmacogenomics Intelligence Agent -- autonomous reasoning across data silos.

Implements the plan -> search -> evaluate -> synthesize -> report pattern from the
intelligent agent. The agent can:

1. Parse complex multi-part questions about pharmacogenomics
2. Plan a search strategy across relevant collections
3. Execute multi-collection retrieval via the RAG engine
4. Evaluate evidence quality and completeness
5. Synthesize cross-functional insights with clinical alerts
6. Generate structured reports with PGx-specific formatting

Mapping to the AI platform:
  - agent engine entry point: PGxIntelligenceAgent.run()
  - Plan -> search_plan()
  - Execute -> rag_engine.retrieve()
  - Reflect -> evaluate_evidence()
  - Report -> generate_report()

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from typing import Dict, List

from .knowledge import (
    PHARMACOGENES,
    DRUG_CATEGORIES,
    HLA_DRUG_ASSOCIATIONS,
)
from .models import (
    AgentQuery,
    AlertLevel,
    CrossCollectionResult,
    PGxAlert,
    PGxResponse,
    PGxWorkflowType,
)


# =====================================================================
# Build HLA allele lookup from knowledge.py's drug-keyed associations
# =====================================================================

def _build_hla_alleles() -> dict:
    """Convert knowledge.py HLA_DRUG_ASSOCIATIONS to allele-keyed lookup."""
    alleles = {}
    for key, data in HLA_DRUG_ASSOCIATIONS.items():
        hla = data.get("hla_allele", "")
        if not hla:
            continue
        if hla not in alleles:
            alleles[hla] = {
                "drugs": [],
                "reaction_types": [],
                "populations_at_risk": [],
            }
        drug = data.get("drug", "")
        if drug and drug not in alleles[hla]["drugs"]:
            alleles[hla]["drugs"].append(drug)
        reaction = data.get("reaction_type", "")
        if reaction and reaction not in alleles[hla]["reaction_types"]:
            alleles[hla]["reaction_types"].append(reaction)
        prev = data.get("prevalence_by_population", {})
        for pop, freq in prev.items():
            note = f"{pop}: {freq}"
            if note not in alleles[hla]["populations_at_risk"]:
                alleles[hla]["populations_at_risk"].append(note)
    return alleles

HLA_ALLELES = _build_hla_alleles()


# PGx-specific system prompt for LLM synthesis
PGX_SYSTEM_PROMPT = """You are a pharmacogenomics intelligence agent with deep expertise in:

1. **Gene-Drug Interactions** -- CYP enzymes, transporters, HLA alleles, and their impact on drug metabolism
2. **Phenotype Interpretation** -- metabolizer status (UM/RM/NM/IM/PM), clinical implications
3. **Dosing Guidance** -- genotype-guided dose adjustments per CPIC/DPWG guidelines
4. **Drug Safety** -- HLA-mediated hypersensitivity, enzyme deficiency toxicity risks
5. **Phenoconversion** -- drug-induced inhibition converting genotypic to phenotypic metabolizer status
6. **Population Pharmacogenomics** -- allele frequency variation across ethnicities, health equity

You synthesize evidence from multiple knowledge bases:
- CPIC and DPWG clinical guidelines
- FDA pharmacogenomic drug labels
- PharmGKB clinical annotations
- PharmVar allele definitions
- PubMed clinical evidence
- ClinicalTrials.gov PGx studies
- Population allele frequency databases

Always provide:
- Specific gene-drug pairs with CPIC evidence levels
- Actionable clinical recommendations with alternatives
- Alert levels (CRITICAL/WARNING/INFO) for safety-relevant findings
- Population-specific considerations where relevant
- PMID citations when available

IMPORTANT: Flag CRITICAL alerts for contraindicated gene-drug combinations,
WARNING for dose adjustments needed, INFO for general pharmacogenomic context."""


# Workflow-specific collection boost weights
WORKFLOW_COLLECTION_BOOST: Dict[PGxWorkflowType, Dict[str, float]] = {
    PGxWorkflowType.GENE_QUERY: {
        "pgx_gene_reference": 2.0,
        "pgx_drug_guidelines": 1.5,
        "pgx_population_data": 1.3,
        "pgx_clinical_evidence": 1.2,
    },
    PGxWorkflowType.DRUG_QUERY: {
        "pgx_drug_guidelines": 2.0,
        "pgx_drug_interactions": 1.8,
        "pgx_fda_labels": 1.5,
        "pgx_drug_alternatives": 1.3,
        "pgx_dosing_algorithms": 1.2,
    },
    PGxWorkflowType.PROFILE_QUERY: {
        "pgx_patient_profiles": 2.0,
        "pgx_gene_reference": 1.5,
        "pgx_drug_guidelines": 1.5,
        "pgx_hla_hypersensitivity": 1.3,
    },
    PGxWorkflowType.INTERACTION_QUERY: {
        "pgx_drug_interactions": 2.0,
        "pgx_phenoconversion": 1.8,
        "pgx_drug_guidelines": 1.3,
        "pgx_fda_labels": 1.2,
    },
    PGxWorkflowType.DOSING_QUERY: {
        "pgx_dosing_algorithms": 2.0,
        "pgx_drug_guidelines": 1.8,
        "pgx_clinical_evidence": 1.5,
        "pgx_fda_labels": 1.3,
    },
    PGxWorkflowType.HLA_SCREEN: {
        "pgx_hla_hypersensitivity": 2.5,
        "pgx_population_data": 1.5,
        "pgx_drug_alternatives": 1.3,
        "pgx_drug_guidelines": 1.2,
    },
}


@dataclass
class SearchPlan:
    """Agent's plan for answering a pharmacogenomics question."""
    question: str
    identified_topics: List[str] = field(default_factory=list)
    genes: List[str] = field(default_factory=list)
    drugs: List[str] = field(default_factory=list)
    relevant_workflows: List[PGxWorkflowType] = field(default_factory=list)
    search_strategy: str = "broad"  # broad, targeted, comparative, clinical
    sub_questions: List[str] = field(default_factory=list)


class PGxIntelligenceAgent:
    """Autonomous Pharmacogenomics Intelligence Agent.

    Wraps the multi-collection RAG engine with planning and reasoning
    capabilities. Designed to answer complex cross-functional questions
    about pharmacogenomics, drug-gene interactions, and genotype-guided
    prescribing.

    Example queries this agent handles:
    - "Is codeine safe for a CYP2D6 ultra-rapid metabolizer?"
    - "What dose of warfarin for a CYP2C9 *1/*3, VKORC1 -1639 A/A patient?"
    - "What drug interactions should I watch for with fluoxetine and codeine?"
    - "What does CYP2D6 *4/*4 mean for my patient?"
    - "Should I screen for HLA-B*57:01 before starting abacavir?"
    - "Give me a full PGx profile for this patient"
    - "Compare clopidogrel vs ticagrelor in CYP2C19 poor metabolizers"
    - "How do CYP2D6 allele frequencies differ across populations?"

    Usage:
        agent = PGxIntelligenceAgent(rag_engine)
        response = agent.run("Is codeine safe for a CYP2D6 UM?")
    """

    def __init__(self, rag_engine):
        """Initialize agent with a configured RAG engine.

        Args:
            rag_engine: PGxRAGEngine instance with all collections connected
        """
        self.rag = rag_engine

    def run(self, question: str, **kwargs) -> PGxResponse:
        """Execute the full agent pipeline: plan -> search -> evaluate -> synthesize -> report.

        Args:
            question: Natural language question about pharmacogenomics
            **kwargs: Additional query parameters (patient_id, medication_list,
                      gene_filter, drug_filter)

        Returns:
            PGxResponse with answer, evidence, alerts, and metadata
        """
        # Phase 1: Plan
        plan = self.search_plan(question)

        # Phase 2: Search (via RAG engine) with workflow-boosted collections
        query = AgentQuery(
            question=question,
            gene=kwargs.get("gene_filter", plan.genes[0] if plan.genes else None),
            drug=kwargs.get("drug_filter", plan.drugs[0] if plan.drugs else None),
        )

        evidence = self.rag.retrieve(
            query, workflows=plan.relevant_workflows,
        )

        # Phase 3: Evaluate evidence quality
        quality = self.evaluate_evidence(evidence)

        # Phase 4: If evidence is thin, try sub-questions
        if quality == "insufficient" and plan.sub_questions:
            for sub_q in plan.sub_questions[:2]:
                sub_query = AgentQuery(question=sub_q)
                sub_evidence = self.rag.retrieve(sub_query)
                evidence.hits.extend(sub_evidence.hits)

        # Phase 5: Generate answer (reuse already-retrieved evidence)
        prompt = self.rag._build_prompt(question, evidence)
        from .rag_engine import PGX_SYSTEM_PROMPT as RAG_SYSTEM_PROMPT
        answer = self.rag.llm.generate(
            prompt=prompt,
            system_prompt=RAG_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

        # Phase 6: Extract knowledge used and generate alerts
        knowledge_used = []
        if evidence.knowledge_context:
            for line in evidence.knowledge_context.split("\n"):
                if line.strip().startswith("## "):
                    knowledge_used.append(line.strip().lstrip("# "))

        alerts = self._generate_alerts(plan, evidence)

        return PGxResponse(
            question=question,
            answer=answer,
            evidence=evidence,
            knowledge_used=knowledge_used,
            alerts=alerts,
        )

    def search_plan(self, question: str) -> SearchPlan:
        """Analyze a question and plan the search strategy.

        Identifies genes, drugs, relevant PGx workflow types,
        and decomposes complex questions into sub-queries.

        Args:
            question: The user's question

        Returns:
            SearchPlan with identified topics, genes, drugs, and strategy
        """
        plan = SearchPlan(question=question)
        q_upper = question.upper()

        # Identify pharmacogenes mentioned
        plan.genes = [g for g in PHARMACOGENES if g.upper() in q_upper]

        # Also check HLA alleles
        for hla in HLA_ALLELES:
            hla_check = hla.upper().replace("*", "").replace(":", "")
            q_check = q_upper.replace("*", "").replace(":", "")
            if hla_check in q_check:
                plan.genes.append(hla)

        # Identify drugs mentioned
        for category, cat_info in DRUG_CATEGORIES.items():
            for drug in cat_info.get("drugs", []):
                if drug.upper() in q_upper:
                    plan.drugs.append(drug)
                    if category not in plan.identified_topics:
                        plan.identified_topics.append(category)

        # Also check gene substrates for implicit drug mentions
        for gene_name, gene_info in PHARMACOGENES.items():
            if gene_name.upper() in q_upper:
                plan.identified_topics.append(gene_name)

        # Determine PGx workflow types
        workflow_keywords = {
            PGxWorkflowType.GENE_QUERY: [
                "GENE", "ALLELE", "STAR ALLELE", "DIPLOTYPE", "GENOTYPE",
                "VARIANT", "POLYMORPHISM", "PHENOTYPE", "METABOLIZER",
                "ACTIVITY SCORE",
            ],
            PGxWorkflowType.DRUG_QUERY: [
                "DRUG", "MEDICATION", "PRESCRI", "SAFE", "CONTRAINDIC",
                "LABEL", "FDA",
            ],
            PGxWorkflowType.PROFILE_QUERY: [
                "PROFILE", "COMPREHENSIVE", "ALL GENES", "FULL",
                "PASSPORT", "PANEL",
            ],
            PGxWorkflowType.INTERACTION_QUERY: [
                "INTERACTION", "DDI", "DRUG-DRUG", "INHIBIT",
                "PHENOCONVERSION", "CONCOMITANT", "CO-PRESCRI",
            ],
            PGxWorkflowType.DOSING_QUERY: [
                "DOSE", "DOSING", "DOSAGE", "TITRAT", "ALGORITHM",
                "ADJUST", "REDUCE", "INCREASE", "MG",
            ],
            PGxWorkflowType.HLA_SCREEN: [
                "HLA", "HYPERSENSITIV", "SJS", "TEN", "DRESS",
                "SCREEN", "ALLERGY", "ADVERSE REACTION",
            ],
        }

        for workflow, keywords in workflow_keywords.items():
            if any(kw in q_upper for kw in keywords):
                plan.relevant_workflows.append(workflow)

        # Default workflow if none detected
        if not plan.relevant_workflows:
            if plan.genes:
                plan.relevant_workflows.append(PGxWorkflowType.GENE_QUERY)
            elif plan.drugs:
                plan.relevant_workflows.append(PGxWorkflowType.DRUG_QUERY)
            else:
                plan.relevant_workflows.append(PGxWorkflowType.GENE_QUERY)

        # Determine search strategy
        if "COMPARE" in q_upper or " VS " in q_upper or "VERSUS" in q_upper:
            plan.search_strategy = "comparative"
        elif "POPULATION" in q_upper or "ETHNIC" in q_upper or "FREQUENC" in q_upper:
            plan.search_strategy = "clinical"
        elif (plan.genes and len(plan.genes) <= 1) and (plan.drugs and len(plan.drugs) <= 1):
            plan.search_strategy = "targeted"
        elif plan.genes or plan.drugs:
            plan.search_strategy = "targeted"
        else:
            plan.search_strategy = "broad"

        # Decompose complex questions into sub-queries (8 patterns)
        gene_label = plan.genes[0] if plan.genes else "pharmacogene"
        drug_label = plan.drugs[0] if plan.drugs else "the medication"

        # Pattern 1: DRUG_SAFETY
        if ("SAFE" in q_upper or "CONTRAINDIC" in q_upper or "AVOID" in q_upper
                or "RISK" in q_upper):
            plan.sub_questions = [
                f"What is the {gene_label} metabolizer phenotype impact on {drug_label}?",
                f"Are there phenoconversion risks with {drug_label} and concomitant medications?",
                f"What are alternative drugs to {drug_label} for patients with {gene_label} variants?",
            ]

        # Pattern 2: DOSING
        elif ("DOSE" in q_upper or "DOSING" in q_upper or "DOSAGE" in q_upper
              or "ADJUST" in q_upper):
            plan.sub_questions = [
                f"What are the CPIC dosing recommendations for {drug_label} based on {gene_label}?",
                f"What genotype-guided dosing algorithms exist for {drug_label}?",
                f"What monitoring parameters should be used when adjusting {drug_label} dose?",
            ]

        # Pattern 3: INTERACTION
        elif ("INTERACTION" in q_upper or "DDI" in q_upper
              or "PHENOCONVERSION" in q_upper):
            plan.sub_questions = [
                f"What drug-drug interactions affect {gene_label} metabolism?",
                f"What is the phenoconversion risk with {drug_label} and CYP inhibitors?",
                f"What is the clinical cascade of {gene_label} phenoconversion?",
            ]

        # Pattern 4: GENE_INTERPRETATION
        elif ("WHAT DOES" in q_upper or "MEAN" in q_upper or "INTERPRET" in q_upper
              or "DIPLOTYPE" in q_upper):
            plan.sub_questions = [
                f"What phenotype does {gene_label} diplotype confer?",
                f"What drugs are affected by {gene_label} variant status?",
                f"What are the clinical implications of {gene_label} variant phenotype?",
            ]

        # Pattern 5: HLA_SCREENING
        elif "HLA" in q_upper or "HYPERSENSITIV" in q_upper or "SJS" in q_upper:
            plan.sub_questions = [
                f"What is the {gene_label} allele frequency in relevant populations?",
                f"What hypersensitivity reactions are associated with {drug_label}?",
                f"What alternative medications are available if {gene_label} is positive?",
            ]

        # Pattern 6: PROFILE
        elif "PROFILE" in q_upper or "COMPREHENSIVE" in q_upper or "PASSPORT" in q_upper:
            plan.sub_questions = [
                "What is the metabolizer status for each pharmacogene?",
                "What critical gene-drug alerts should be generated?",
                "What drugs require dose adjustment or avoidance based on the profile?",
            ]

        # Pattern 7: COMPARISON
        elif "COMPARE" in q_upper or " VS " in q_upper or "VERSUS" in q_upper:
            items = plan.identified_topics[:2] if len(plan.identified_topics) >= 2 else [drug_label, "alternative"]
            plan.sub_questions = [
                f"What are the metabolic pathways of {items[0]}?",
                f"How does genotype affect efficacy of {items[0]} vs {items[1] if len(items) > 1 else 'alternative'}?",
            ]

        # Pattern 8: POPULATION
        elif "POPULATION" in q_upper or "ETHNIC" in q_upper or "FREQUENC" in q_upper:
            plan.sub_questions = [
                f"What is the allele frequency of {gene_label} variants across populations?",
                f"How do population differences in {gene_label} affect drug response?",
                f"What health equity considerations exist for {gene_label} pharmacogenomics?",
            ]

        return plan

    def evaluate_evidence(self, evidence: CrossCollectionResult) -> str:
        """Evaluate the quality and coverage of retrieved evidence.

        Uses the same thresholds as the CAR-T agent:
          - sufficient: >= 3 collections with >= 10 total hits
          - partial: >= 2 collections with >= 5 total hits
          - insufficient: everything else

        Returns:
            "sufficient", "partial", or "insufficient"
        """
        if evidence.hit_count == 0:
            return "insufficient"

        by_coll = evidence.hits_by_collection()
        collections_with_hits = len(by_coll)

        if collections_with_hits >= 3 and evidence.hit_count >= 10:
            return "sufficient"
        elif collections_with_hits >= 2 and evidence.hit_count >= 5:
            return "partial"
        else:
            return "insufficient"

    def generate_report(self, response: PGxResponse) -> str:
        """Generate a structured pharmacogenomics analysis report.

        Args:
            response: PGxResponse from a completed query

        Returns:
            Formatted markdown report
        """
        by_coll = response.evidence.hits_by_collection()

        report_lines = [
            "# Pharmacogenomics Intelligence Report",
            f"**Query:** {response.question}",
            f"**Generated:** {response.timestamp}",
            f"**Collections Searched:** {response.evidence.total_collections_searched}",
            f"**Evidence Items:** {response.evidence.hit_count}",
            f"**Search Time:** {response.evidence.search_time_ms:.0f} ms",
            "",
        ]

        # Alerts section
        if response.alerts:
            report_lines.extend([
                "---",
                "",
                "## Clinical Alerts",
                "",
            ])
            for alert in response.alerts:
                {"critical": "!!!", "warning": "!!", "info": "i"}.get(
                    alert.alert_level.value, "i"
                )
                report_lines.append(
                    f"- **[{alert.alert_level.value.upper()}]** "
                    f"{alert.gene} / {alert.drug}: {alert.recommendation}"
                )
            report_lines.append("")

        report_lines.extend([
            "---",
            "",
            "## Analysis",
            "",
            response.answer,
            "",
            "---",
            "",
            "## Evidence Sources",
            "",
        ])

        for coll_name, hits in by_coll.items():
            report_lines.append(f"### {coll_name} ({len(hits)} results)")
            for hit in hits[:5]:
                report_lines.append(
                    f"- **{hit.id}** (score: {hit.score:.3f}): "
                    f"{hit.text[:200]}..."
                )
            report_lines.append("")

        if response.knowledge_used:
            report_lines.extend([
                "## Knowledge Graph",
                "",
            ])
            for k in response.knowledge_used:
                report_lines.append(f"- {k}")

        return "\n".join(report_lines)

    # ── Private helpers ──────────────────────────────────────────────

    def _generate_alerts(
        self, plan: SearchPlan, evidence: CrossCollectionResult
    ) -> List[PGxAlert]:
        """Generate clinical decision support alerts based on identified
        gene-drug pairs and evidence.

        Args:
            plan: The search plan with identified genes and drugs
            evidence: Retrieved evidence from collections

        Returns:
            List of PGxAlert instances
        """
        alerts: List[PGxAlert] = []

        # Check for known critical HLA-drug pairs
        for hla_allele, hla_info in HLA_ALLELES.items():
            if hla_allele in plan.genes:
                for drug in hla_info.get("drugs", []):
                    if drug in plan.drugs or not plan.drugs:
                        reaction = ", ".join(hla_info.get("reaction_types", ["adverse reaction"]))
                        alerts.append(PGxAlert(
                            alert_level=AlertLevel.CRITICAL,
                            gene=hla_allele,
                            drug=drug,
                            phenotype=f"{hla_allele} carrier",
                            recommendation=(
                                f"Screen for {hla_allele} before prescribing {drug}. "
                                f"Risk of {reaction}."
                            ),
                        ))

        # Check for critical enzyme-drug combinations
        critical_combos = {
            ("DPYD", "fluorouracil"): "DPYD deficiency can cause fatal fluoropyrimidine toxicity. Pre-treatment genotyping required.",
            ("DPYD", "capecitabine"): "DPYD deficiency can cause fatal fluoropyrimidine toxicity. Pre-treatment genotyping required.",
            ("G6PD", "rasburicase"): "G6PD deficiency: rasburicase is CONTRAINDICATED. Risk of severe hemolytic anemia.",
            ("TPMT", "azathioprine"): "TPMT deficiency requires significant dose reduction of azathioprine to prevent myelosuppression.",
            ("TPMT", "mercaptopurine"): "TPMT deficiency requires significant dose reduction of mercaptopurine to prevent myelosuppression.",
        }

        for (gene, drug), recommendation in critical_combos.items():
            if gene in plan.genes and drug in plan.drugs:
                alerts.append(PGxAlert(
                    alert_level=AlertLevel.CRITICAL,
                    gene=gene,
                    drug=drug,
                    phenotype="deficient/poor metabolizer",
                    recommendation=recommendation,
                ))

        # Check for warning-level dose adjustment pairs
        dose_adjust_combos = {
            ("CYP2D6", "codeine"): "CYP2D6 ultra-rapid metabolizers: increased morphine formation, respiratory depression risk. CYP2D6 poor metabolizers: lack of analgesic effect.",
            ("CYP2C19", "clopidogrel"): "CYP2C19 poor/intermediate metabolizers: reduced clopidogrel activation, increased cardiovascular event risk. Consider alternative antiplatelet.",
            ("SLCO1B1", "simvastatin"): "SLCO1B1 decreased function: increased simvastatin exposure, myopathy risk 5-17x. Consider dose reduction or alternative statin.",
            ("CYP2C9", "warfarin"): "CYP2C9 poor metabolizer: reduced warfarin clearance. Use genotype-guided dosing algorithm.",
        }

        for (gene, drug), recommendation in dose_adjust_combos.items():
            if gene in plan.genes and drug in plan.drugs:
                alerts.append(PGxAlert(
                    alert_level=AlertLevel.WARNING,
                    gene=gene,
                    drug=drug,
                    phenotype="variant metabolizer",
                    recommendation=recommendation,
                ))

        return alerts
