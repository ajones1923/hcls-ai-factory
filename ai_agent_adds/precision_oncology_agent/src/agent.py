"""
Oncology Intelligence Agent.

Plan-search-synthesize pattern for multi-collection oncology evidence retrieval.
Analyzes queries for gene targets, cancer types, and therapeutic contexts,
then orchestrates cross-collection searches with adaptive retry logic.

Author: Adam Jones
Date: February 2026
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Dict, List, Set

from .models import AgentQuery, AgentResponse, CrossCollectionResult

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Known oncology gene panels and cancer type vocabularies
# ---------------------------------------------------------------------------

KNOWN_GENES: Set[str] = {
    "BRAF", "EGFR", "ALK", "ROS1", "KRAS", "HER2",
    "NTRK", "RET", "MET", "FGFR", "PIK3CA", "IDH1",
    "IDH2", "BRCA", "BRCA1", "BRCA2", "TP53", "PTEN", "CDKN2A",
    "STK11", "ESR1", "ERBB2", "NRAS", "APC", "VHL",
    "KIT", "PDGFRA", "FLT3", "NPM1", "DNMT3A",
}

KNOWN_CANCER_TYPES: Set[str] = {
    "NSCLC", "BREAST", "MELANOMA", "COLORECTAL", "PANCREATIC",
    "OVARIAN", "PROSTATE", "GLIOMA", "GLIOBLASTOMA", "AML",
    "CML", "CLL", "DLBCL", "BLADDER", "RENAL", "HEPATOCELLULAR",
    "GASTRIC", "ESOPHAGEAL", "THYROID", "ENDOMETRIAL", "CERVICAL",
    "HEAD_AND_NECK", "SARCOMA", "CHOLANGIOCARCINOMA", "MESOTHELIOMA",
}

# Aliases so natural-language mentions resolve to canonical names
_CANCER_ALIASES: Dict[str, str] = {
    "lung": "NSCLC",
    "lung cancer": "NSCLC",
    "non-small cell lung": "NSCLC",
    "non small cell lung": "NSCLC",
    "nsclc": "NSCLC",
    "lung adenocarcinoma": "NSCLC",
    "small cell lung": "SCLC",
    "sclc": "SCLC",
    "breast cancer": "BREAST",
    "breast": "BREAST",
    "triple negative breast": "BREAST",
    "tnbc": "BREAST",
    "colon": "COLORECTAL",
    "colon cancer": "COLORECTAL",
    "colorectal cancer": "COLORECTAL",
    "crc": "COLORECTAL",
    "rectal cancer": "COLORECTAL",
    "melanoma": "MELANOMA",
    "skin cancer": "MELANOMA",
    "cutaneous melanoma": "MELANOMA",
    "pancreatic": "PANCREATIC",
    "pancreatic cancer": "PANCREATIC",
    "pdac": "PANCREATIC",
    "ovarian": "OVARIAN",
    "ovarian cancer": "OVARIAN",
    "prostate": "PROSTATE",
    "prostate cancer": "PROSTATE",
    "crpc": "PROSTATE",
    "glioma": "GLIOMA",
    "gbm": "GLIOBLASTOMA",
    "glioblastoma": "GLIOBLASTOMA",
    "aml": "AML",
    "acute myeloid": "AML",
    "acute myeloid leukemia": "AML",
    "cml": "CML",
    "chronic myeloid": "CML",
    "chronic myeloid leukemia": "CML",
    "cll": "CLL",
    "chronic lymphocytic leukemia": "CLL",
    "bladder": "BLADDER",
    "bladder cancer": "BLADDER",
    "urothelial": "BLADDER",
    "renal": "RENAL",
    "kidney": "RENAL",
    "kidney cancer": "RENAL",
    "rcc": "RENAL",
    "liver": "HEPATOCELLULAR",
    "hcc": "HEPATOCELLULAR",
    "liver cancer": "HEPATOCELLULAR",
    "hepatocellular carcinoma": "HEPATOCELLULAR",
    "gastric": "GASTRIC",
    "stomach": "GASTRIC",
    "stomach cancer": "GASTRIC",
    "gastric cancer": "GASTRIC",
    "thyroid": "THYROID",
    "thyroid cancer": "THYROID",
    "endometrial": "ENDOMETRIAL",
    "uterine": "ENDOMETRIAL",
    "endometrial cancer": "ENDOMETRIAL",
    "sarcoma": "SARCOMA",
    "esophageal": "ESOPHAGEAL",
    "esophageal cancer": "ESOPHAGEAL",
    "cholangiocarcinoma": "CHOLANGIOCARCINOMA",
    "bile duct cancer": "CHOLANGIOCARCINOMA",
    "head and neck": "HEAD_AND_NECK",
    "hnscc": "HEAD_AND_NECK",
    "cervical": "CERVICAL",
    "cervical cancer": "CERVICAL",
    "mesothelioma": "MESOTHELIOMA",
}


# ---------------------------------------------------------------------------
# SearchPlan dataclass
# ---------------------------------------------------------------------------

@dataclass
class SearchPlan:
    """Structured plan produced by the agent before executing searches."""

    question: str
    identified_topics: List[str] = field(default_factory=list)
    target_genes: List[str] = field(default_factory=list)
    relevant_cancer_types: List[str] = field(default_factory=list)
    search_strategy: str = "broad"  # "broad" | "targeted" | "comparative"
    sub_questions: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# OncoIntelligenceAgent
# ---------------------------------------------------------------------------

class OncoIntelligenceAgent:
    """Oncology Intelligence Agent.

    Implements the plan -> search -> evaluate -> (retry) -> synthesize loop
    to answer complex oncology questions by aggregating evidence across
    multiple Milvus collections (literature, clinical trials, CIViC, etc.).
    """

    # Retry configuration
    MAX_RETRIES: int = 2
    MIN_SUFFICIENT_HITS: int = 3
    MIN_COLLECTIONS_FOR_SUFFICIENT: int = 2

    def __init__(self, rag_engine) -> None:
        """Initialise the agent.

        Parameters
        ----------
        rag_engine:
            Backend responsible for embedding queries, searching Milvus
            collections, and calling the LLM for synthesis.
        """
        self.rag_engine = rag_engine

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def run(self, question: str, **kwargs) -> AgentResponse:
        """Execute the full plan-search-synthesize pipeline.

        Parameters
        ----------
        question:
            Free-text oncology question from the user or an upstream system.
        **kwargs:
            Forwarded to the RAG engine search calls (e.g. ``top_k``,
            ``collection_filter``).

        Returns
        -------
        AgentResponse
            Consolidated response with evidence, citations, and a
            structured markdown report.
        """
        start = time.time()
        logger.info("OncoIntelligenceAgent.run  question=%s", question[:120])

        # 1. Plan -----------------------------------------------------------
        plan = self.search_plan(question)
        logger.info(
            "Search plan: strategy=%s  genes=%s  cancers=%s  sub_questions=%d",
            plan.search_strategy,
            plan.target_genes,
            plan.relevant_cancer_types,
            len(plan.sub_questions),
        )

        # 2. Search (with adaptive retry) -----------------------------------
        all_evidence: List[CrossCollectionResult] = []
        queries_to_run = [plan.question] + plan.sub_questions

        for attempt in range(1, self.MAX_RETRIES + 2):  # 1-based
            for q in queries_to_run:
                query = AgentQuery(question=q)
                results = self.rag_engine.cross_collection_search(query)
                if results:
                    all_evidence.extend(results)

            # 3. Evaluate ----------------------------------------------------
            verdict = self.evaluate_evidence(all_evidence)
            logger.info(
                "Evidence evaluation (attempt %d/%d): %s  hits=%d",
                attempt,
                self.MAX_RETRIES + 1,
                verdict,
                len(all_evidence),
            )

            if verdict == "sufficient" or attempt > self.MAX_RETRIES:
                break

            # Broaden the search for the next attempt
            logger.info("Insufficient evidence – broadening search for retry")
            if plan.search_strategy == "targeted":
                plan.search_strategy = "broad"
            queries_to_run = self._generate_fallback_queries(plan)

        # 4. Synthesize -----------------------------------------------------
        response: AgentResponse = self.rag_engine.synthesize(
            question=question,
            evidence=all_evidence,
            plan=plan,
        )

        # Attach plan and formatted report
        response.plan = plan
        response.report = self.generate_report(response)

        elapsed = time.time() - start
        logger.info("OncoIntelligenceAgent.run completed in %.2fs", elapsed)
        return response

    # ------------------------------------------------------------------
    # Planning
    # ------------------------------------------------------------------

    def search_plan(self, question: str) -> SearchPlan:
        """Analyse the question and produce a structured search plan.

        Identifies gene targets, cancer types, and selects an appropriate
        search strategy. Complex questions are decomposed into focused
        sub-queries.
        """
        q_upper = question.upper()
        q_lower = question.lower()

        # --- Identify genes ------------------------------------------------
        target_genes: List[str] = [
            gene for gene in KNOWN_GENES if gene in q_upper
        ]

        # --- Identify cancer types -----------------------------------------
        relevant_cancer_types: List[str] = []
        # Check canonical names
        for ct in KNOWN_CANCER_TYPES:
            if ct in q_upper:
                relevant_cancer_types.append(ct)
        # Check aliases
        for alias, canonical in _CANCER_ALIASES.items():
            if alias in q_lower and canonical not in relevant_cancer_types:
                relevant_cancer_types.append(canonical)

        # --- Identify topics -----------------------------------------------
        identified_topics: List[str] = []
        topic_keywords = {
            "resistance": "therapeutic resistance",
            "biomarker": "biomarker identification",
            "prognosis": "prognostic significance",
            "survival": "survival outcomes",
            "immunotherapy": "immunotherapy response",
            "targeted therapy": "targeted therapy",
            "clinical trial": "clinical trials",
            "combination": "combination therapy",
            "mutation": "mutation landscape",
            "variant": "variant interpretation",
            "expression": "gene expression",
            "amplification": "gene amplification",
            "fusion": "gene fusion",
            "methylation": "epigenetic regulation",
            "liquid biopsy": "liquid biopsy / ctDNA",
            "ctdna": "liquid biopsy / ctDNA",
            "pdl1": "PD-L1 / immune checkpoint",
            "pd-l1": "PD-L1 / immune checkpoint",
            "checkpoint": "PD-L1 / immune checkpoint",
            "tumor mutational burden": "TMB",
            "tmb": "TMB",
            "microsatellite": "MSI / microsatellite instability",
            "msi": "MSI / microsatellite instability",
        }
        for kw, topic in topic_keywords.items():
            if kw in q_lower and topic not in identified_topics:
                identified_topics.append(topic)

        # --- Select strategy -----------------------------------------------
        comparative_signals = {"compare", "vs", "versus", "difference between", "head to head"}
        if any(sig in q_lower for sig in comparative_signals):
            search_strategy = "comparative"
        elif target_genes and relevant_cancer_types:
            search_strategy = "targeted"
        else:
            search_strategy = "broad"

        # --- Decompose into sub-questions ----------------------------------
        sub_questions = self._decompose_question(
            question, target_genes, relevant_cancer_types, identified_topics,
        )

        return SearchPlan(
            question=question,
            identified_topics=identified_topics,
            target_genes=target_genes,
            relevant_cancer_types=relevant_cancer_types,
            search_strategy=search_strategy,
            sub_questions=sub_questions,
        )

    # ------------------------------------------------------------------
    # Evidence evaluation
    # ------------------------------------------------------------------

    # Minimum similarity score for evidence to count as relevant
    MIN_SIMILARITY_SCORE: float = 0.30

    def evaluate_evidence(self, evidence: List[CrossCollectionResult]) -> str:
        """Rate the adequacy of retrieved evidence.

        Considers hit count, collection diversity, and similarity scores.
        Evidence items with scores below MIN_SIMILARITY_SCORE are treated
        as low-quality and discounted.

        Returns
        -------
        str
            ``"sufficient"`` if we have enough evidence to synthesize a
            confident answer, ``"partial"`` if there is some relevant data
            but gaps remain, and ``"insufficient"`` if there is too little
            to work with.
        """
        if not evidence:
            return "insufficient"

        # Filter to items with meaningful similarity scores
        scored_evidence = []
        for e in evidence:
            score = getattr(e, "score", None)
            if isinstance(score, (int, float)) and score < self.MIN_SIMILARITY_SCORE:
                continue
            scored_evidence.append(e)

        hit_count = len(scored_evidence)
        collections_represented: Set[str] = {
            e.collection for e in scored_evidence if hasattr(e, "collection")
        }

        # Calculate average score for quality assessment
        scores = [
            getattr(e, "score", 0) for e in scored_evidence
            if isinstance(getattr(e, "score", None), (int, float))
        ]
        avg_score = sum(scores) / len(scores) if scores else 0.0

        if (
            hit_count >= self.MIN_SUFFICIENT_HITS
            and len(collections_represented) >= self.MIN_COLLECTIONS_FOR_SUFFICIENT
            and avg_score >= 0.50
        ):
            return "sufficient"

        if (
            hit_count >= self.MIN_SUFFICIENT_HITS
            and len(collections_represented) >= self.MIN_COLLECTIONS_FOR_SUFFICIENT
        ):
            # Enough hits and diversity but lower scores — still usable
            return "sufficient"

        if hit_count > 0:
            return "partial"

        return "insufficient"

    # ------------------------------------------------------------------
    # Report generation
    # ------------------------------------------------------------------

    def generate_report(self, response: AgentResponse) -> str:
        """Produce a formatted Markdown report from the agent response.

        Sections:
        - Query & Analysis
        - Evidence Sources
        - Knowledge Graph (if available)
        - Synthesis / Answer
        """
        sections: List[str] = []

        # Header
        sections.append("# Oncology Intelligence Report")
        sections.append("")

        # Query section
        sections.append("## Query")
        sections.append(f"**Question:** {response.question}")
        sections.append("")

        # Analysis section
        if response.plan:
            plan = response.plan
            sections.append("## Analysis")
            sections.append(f"- **Strategy:** {plan.search_strategy}")
            if plan.target_genes:
                sections.append(f"- **Target genes:** {', '.join(plan.target_genes)}")
            if plan.relevant_cancer_types:
                sections.append(f"- **Cancer types:** {', '.join(plan.relevant_cancer_types)}")
            if plan.identified_topics:
                sections.append(f"- **Topics:** {', '.join(plan.identified_topics)}")
            if plan.sub_questions:
                sections.append("- **Sub-questions:**")
                for sq in plan.sub_questions:
                    sections.append(f"  - {sq}")
            sections.append("")

        # Evidence sources
        if response.evidence:
            sections.append("## Evidence Sources")
            sections.append(f"Total evidence items: **{response.evidence.hit_count}**")
            sections.append("")

            # Group by collection
            by_collection: Dict[str, List] = {}
            for item in response.evidence.hits:
                col = getattr(item, "collection", "unknown")
                by_collection.setdefault(col, []).append(item)

            for collection, items in sorted(by_collection.items()):
                sections.append(f"### {collection} ({len(items)} hits)")
                for idx, item in enumerate(items[:10], 1):  # cap display
                    title = getattr(item, "title", None) or getattr(item, "id", f"item-{idx}")
                    score = getattr(item, "score", None)
                    score_str = f" (score: {score:.3f})" if score is not None else ""
                    sections.append(f"{idx}. {title}{score_str}")
                if len(items) > 10:
                    sections.append(f"   _...and {len(items) - 10} more_")
                sections.append("")

        # Knowledge graph section
        if hasattr(response, "knowledge_graph") and response.knowledge_graph:
            sections.append("## Knowledge Graph")
            kg = response.knowledge_graph
            if isinstance(kg, dict):
                for key, value in kg.items():
                    sections.append(f"- **{key}:** {value}")
            else:
                sections.append(str(kg))
            sections.append("")

        # Synthesis
        sections.append("## Synthesis")
        answer_text = getattr(response, "answer", None) or getattr(response, "synthesis", "")
        sections.append(str(answer_text))
        sections.append("")

        # Footer
        sections.append("---")
        sections.append("*Generated by Oncology Intelligence Agent*")

        return "\n".join(sections)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _decompose_question(
        self,
        question: str,
        genes: List[str],
        cancer_types: List[str],
        topics: List[str],
    ) -> List[str]:
        """Break a complex question into focused sub-queries."""
        sub_questions: List[str] = []

        # If multiple genes, create per-gene queries
        if len(genes) > 1:
            for gene in genes:
                cancer_ctx = f" in {cancer_types[0]}" if cancer_types else ""
                sub_questions.append(
                    f"What is the role of {gene}{cancer_ctx}?"
                )

        # If multiple cancer types, create per-cancer queries
        if len(cancer_types) > 1:
            gene_ctx = f"{genes[0]} " if genes else ""
            for ct in cancer_types:
                sub_questions.append(
                    f"{gene_ctx}therapeutic landscape in {ct}"
                )

        # Topic-driven sub-questions
        if "therapeutic resistance" in topics and genes:
            sub_questions.append(
                f"Mechanisms of resistance to {genes[0]} inhibitors"
            )
        if "clinical trials" in topics:
            gene_ctx = f" targeting {genes[0]}" if genes else ""
            cancer_ctx = f" in {cancer_types[0]}" if cancer_types else ""
            sub_questions.append(
                f"Active clinical trials{gene_ctx}{cancer_ctx}"
            )
        if "biomarker identification" in topics:
            cancer_ctx = f" for {cancer_types[0]}" if cancer_types else ""
            sub_questions.append(
                f"Predictive biomarkers{cancer_ctx}"
            )
        if "combination therapy" in topics and genes:
            sub_questions.append(
                f"Combination strategies with {genes[0]} inhibitors"
            )

        # Deduplicate while preserving order
        seen: Set[str] = set()
        unique: List[str] = []
        for sq in sub_questions:
            if sq not in seen:
                seen.add(sq)
                unique.append(sq)

        return unique

    def _generate_fallback_queries(self, plan: SearchPlan) -> List[str]:
        """Generate broader fallback queries when initial search is insufficient."""
        fallbacks: List[str] = []

        if plan.target_genes:
            for gene in plan.target_genes:
                fallbacks.append(f"{gene} oncology therapeutic implications")
                fallbacks.append(f"{gene} mutation clinical significance")

        if plan.relevant_cancer_types:
            for ct in plan.relevant_cancer_types:
                fallbacks.append(f"{ct} precision medicine current landscape")

        if not fallbacks:
            fallbacks.append(f"{plan.question} precision oncology")

        return fallbacks
