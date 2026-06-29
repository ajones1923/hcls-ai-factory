"""Biomarker Intelligence Agent -- autonomous reasoning across biomarker data silos.

Implements the plan -> analyze -> search -> synthesize -> report pattern from the
VAST AI OS AgentEngine model. The agent integrates four analysis modules
(biological age, disease trajectory, pharmacogenomics, genotype adjustment)
alongside the multi-collection RAG engine.

Key differences from other Biomarker Intelligence Agent designs:
- Integrates 4 analysis modules alongside RAG (not RAG-only)
- Patient profile drives module execution (biological age, trajectories, PGx)
- Critical alerts extracted across all analysis results
- Cross-modal triggers generated for other HCLS AI Factory agents

Mapping to VAST AI OS:
  - AgentEngine entry point: PrecisionBiomarkerAgent.run()
  - Plan -> search_plan()
  - Analyze -> analyze_patient()
  - Execute -> rag_engine.retrieve()
  - Reflect -> evaluate_evidence()
  - Report -> report via ReportGenerator (separate module)

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from typing import List, Optional

from loguru import logger

from .biological_age import BiologicalAgeCalculator
from .critical_values import CriticalValueEngine
from .discordance_detector import DiscordanceDetector
from .disease_trajectory import DiseaseTrajectoryAnalyzer
from .genotype_adjustment import GenotypeAdjuster
from .lab_range_interpreter import LabRangeInterpreter
from .models import (
    AgentQuery,
    AgentResponse,
    AnalysisResult,
    BiologicalAgeResult,
    CrossCollectionResult,
    DiseaseCategory,
    DiseaseTrajectoryResult,
    GenotypeAdjustmentResult,
    MetabolizerPhenotype,
    PatientProfile,
    PGxResult,
    RiskLevel,
)
from .pharmacogenomics import PharmacogenomicMapper
from .rag_engine import BIOMARKER_SYSTEM_PROMPT


@dataclass
class SearchPlan:
    """Agent's plan for answering a question."""
    question: str
    identified_topics: List[str] = field(default_factory=list)
    disease_areas: List[str] = field(default_factory=list)
    relevant_modules: List[str] = field(default_factory=list)
    search_strategy: str = "broad"  # broad, targeted, domain-specific
    sub_questions: List[str] = field(default_factory=list)


class PrecisionBiomarkerAgent:
    """Autonomous Biomarker Intelligence Agent.

    Wraps the multi-collection RAG engine with planning, analysis modules,
    and reasoning capabilities. Designed to answer complex cross-functional
    questions about precision biomarker interpretation, pharmacogenomics,
    biological aging, and disease trajectory analysis.

    Analysis Modules:
    - BiologicalAgeCalculator: PhenoAge and GrimAge surrogate estimation
    - DiseaseTrajectoryAnalyzer: Pre-symptomatic disease trajectory analysis
    - PharmacogenomicMapper: CPIC-guided star allele to phenotype mapping
    - GenotypeAdjuster: Genotype-adjusted reference range calculation

    Example queries this agent handles:
    - "What does my HbA1c of 5.8% mean given my TCF7L2 CT genotype?"
    - "Calculate my biological age from these blood markers"
    - "What are the PGx implications of CYP2D6 *4/*4 for pain management?"
    - "Assess my cardiovascular risk given elevated Lp(a) and APOE e3/e4"
    - "How do my MTHFR and FADS1 genotypes affect my supplement protocol?"

    Usage:
        agent = PrecisionBiomarkerAgent(rag_engine)
        response = agent.run("What does my HbA1c of 5.8% mean?")
        # Or with patient profile:
        response = agent.run(question, patient_profile=profile)
    """

    def __init__(self, rag_engine, bio_age_calc=None, trajectory_analyzer=None,
                 pgx_mapper=None, genotype_adjuster=None):
        """Initialize agent with a configured RAG engine and optional analysis modules.

        Args:
            rag_engine: BiomarkerRAGEngine instance with all collections connected.
            bio_age_calc: BiologicalAgeCalculator instance (creates default if None).
            trajectory_analyzer: DiseaseTrajectoryAnalyzer instance.
            pgx_mapper: PharmacogenomicMapper instance.
            genotype_adjuster: GenotypeAdjuster instance.
        """
        self.rag = rag_engine
        self.bio_age = bio_age_calc or BiologicalAgeCalculator()
        self.trajectory = trajectory_analyzer or DiseaseTrajectoryAnalyzer()
        self.pgx = pgx_mapper or PharmacogenomicMapper()
        self.adjuster = genotype_adjuster or GenotypeAdjuster()
        self.critical_values = CriticalValueEngine()
        self.discordance = DiscordanceDetector()
        self.lab_ranges = LabRangeInterpreter()

    def run(self, question: str, patient_profile: Optional[PatientProfile] = None,
            **kwargs) -> AgentResponse:
        """Execute the full agent pipeline: plan -> analyze -> search -> synthesize.

        Args:
            question: Natural language question about biomarkers, PGx, etc.
            patient_profile: Optional patient data for personalized analysis.
            **kwargs: Additional query parameters (collections_filter, year_min, etc.)

        Returns:
            AgentResponse with answer, evidence, and analysis results.
        """
        logger.info(f"Agent run started: question='{question[:80]}...', has_profile={patient_profile is not None}")

        # Phase 1: If patient profile provided, run analysis modules
        analysis = None
        if patient_profile:
            analysis = self.analyze_patient(patient_profile)

        # Phase 2: Plan search strategy
        plan = self.search_plan(question)

        # Phase 3: Search via RAG engine
        query = AgentQuery(
            question=question,
            patient_profile=patient_profile,
            include_genomic=kwargs.get("include_genomic", True),
        )
        # NOTE: rag_engine.retrieve() accepts conversation_context for multi-turn
        # support. Currently unused -- wire through when AgentQuery gains a
        # conversation_context field or run() accepts conversation history.
        evidence = self.rag.retrieve(
            query,
            collections_filter=kwargs.get("collections_filter"),
            year_min=kwargs.get("year_min"),
            year_max=kwargs.get("year_max"),
        )

        # Phase 4: Evaluate evidence quality
        quality = self.evaluate_evidence(evidence)

        # Phase 5: If evidence is thin, try sub-questions
        if quality == "insufficient" and plan.sub_questions:
            for sub_q in plan.sub_questions[:2]:
                sub_query = AgentQuery(question=sub_q, include_genomic=False)
                sub_evidence = self.rag.retrieve(sub_query)
                evidence.hits.extend(sub_evidence.hits)

        # Phase 6: Build enhanced prompt with analysis results
        prompt = self._build_enhanced_prompt(question, evidence, analysis)

        # Phase 7: Generate answer via LLM
        answer = self.rag.llm.generate(
            prompt=prompt,
            system_prompt=BIOMARKER_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

        # Build response
        return AgentResponse(
            question=question,
            answer=answer,
            evidence=evidence,
            biological_age=analysis.biological_age if analysis else None,
            disease_trajectories=analysis.disease_trajectories if analysis else None,
            pgx_results=analysis.pgx_results if analysis else None,
            genotype_adjustments=analysis.genotype_adjustments if analysis else None,
            critical_alerts=analysis.critical_alerts if analysis else [],
        )

    def analyze_patient(self, profile: PatientProfile) -> AnalysisResult:
        """Run all analysis modules on a patient profile.

        Executes biological age calculation, disease trajectory analysis,
        pharmacogenomic mapping, and genotype-adjusted reference ranges.

        Args:
            profile: PatientProfile with biomarkers, genotypes, and star alleles.

        Returns:
            AnalysisResult combining all sub-analyses with critical alerts.
        """
        logger.info(f"Analyzing patient {profile.patient_id}: age={profile.age}, sex={profile.sex}, "
                     f"biomarkers={len(profile.biomarkers)}, genotypes={len(profile.genotypes)}, "
                     f"star_alleles={len(profile.star_alleles)}")

        # 1. Biological age calculation
        bio_age_raw = self.bio_age.calculate(profile.age, profile.biomarkers)
        bio_age = BiologicalAgeResult(
            chronological_age=profile.age,
            biological_age=bio_age_raw.get("biological_age", profile.age),
            age_acceleration=bio_age_raw.get("age_acceleration", 0.0),
            phenoage_score=bio_age_raw.get("biological_age", profile.age),
            grimage_score=(
                bio_age_raw["grimage"]["grimage_score"]
                if bio_age_raw.get("grimage") else None
            ),
            grimage_data=bio_age_raw.get("grimage"),
            mortality_risk=bio_age_raw.get("phenoage", {}).get("mortality_score", 0.0),
            aging_drivers=bio_age_raw.get("phenoage", {}).get("top_aging_drivers", []),
            confidence_interval=bio_age_raw.get("phenoage", {}).get("confidence_interval"),
            risk_confidence=bio_age_raw.get("phenoage", {}).get("risk_confidence"),
        )

        # 2. Disease trajectory analysis (returns List[Dict], convert to model objects)
        trajectory_dicts = self.trajectory.analyze_all(
            profile.biomarkers, profile.genotypes, profile.age, profile.sex,
        )
        # Map disease_trajectory.py disease strings to DiseaseCategory enum values
        _disease_str_map = {
            "type2_diabetes": "diabetes",
            "cardiovascular": "cardiovascular",
            "liver": "liver",
            "thyroid": "thyroid",
            "iron": "iron",
            "nutritional": "nutritional",
            "kidney": "kidney",
            "bone_health": "bone_health",
            "cognitive": "cognitive",
        }
        trajectories: List[DiseaseTrajectoryResult] = []
        for td in trajectory_dicts:
            try:
                disease_str = td.get("disease", "")
                disease_str = _disease_str_map.get(disease_str, disease_str)
                disease = DiseaseCategory(disease_str) if disease_str in [e.value for e in DiseaseCategory] else None
                if disease is None:
                    continue
                risk_str = td.get("risk_level", "LOW").lower()
                risk_level = RiskLevel(risk_str) if risk_str in [e.value for e in RiskLevel] else RiskLevel.NORMAL
                # Extract genetic risk factor labels from list of dicts
                genetic_factors = [
                    f"{grf.get('gene', '')} {grf.get('genotype', '')}"
                    for grf in td.get("genetic_risk_factors", [])
                ]
                trajectories.append(DiseaseTrajectoryResult(
                    disease=disease,
                    risk_level=risk_level,
                    current_markers=td.get("current_markers", {}),
                    genetic_risk_factors=genetic_factors,
                    years_to_onset_estimate=td.get("years_to_onset_estimate"),
                    intervention_recommendations=td.get("recommendations", []),
                ))
            except Exception as e:
                logger.warning(f"Skipping malformed trajectory result: {e}")

        # 3. Pharmacogenomic mapping (returns Dict, convert to PGxResult model objects)
        pgx_raw = self.pgx.map_all(profile.star_alleles, profile.genotypes)
        pgx_results: List[PGxResult] = []
        for gr in pgx_raw.get("gene_results", []):
            try:
                phenotype_str = gr.get("phenotype", "normal_metabolizer")
                if phenotype_str is None:
                    continue
                # Normalize phenotype string to enum value
                phenotype_val = phenotype_str.lower().replace(" ", "_").replace("-", "_")
                phenotype = MetabolizerPhenotype(phenotype_val) if phenotype_val in [e.value for e in MetabolizerPhenotype] else MetabolizerPhenotype.NORMAL
                pgx_results.append(PGxResult(
                    gene=gr.get("gene", ""),
                    star_alleles=gr.get("star_alleles", "") or "",
                    phenotype=phenotype,
                    drugs_affected=gr.get("affected_drugs", []),
                ))
            except Exception as e:
                logger.warning(f"Skipping malformed PGx result: {e}")

        # 4. Genotype-adjusted reference ranges (returns Dict, convert to model objects)
        adj_raw = self.adjuster.adjust_all(
            profile.biomarkers, profile.genotypes,
        )
        adjustments: List[GenotypeAdjustmentResult] = []
        for adj in adj_raw.get("adjustments", []):
            try:
                std_range = adj.get("standard_range", {})
                adj_range = adj.get("adjusted_range", {})
                unit = adj.get("unit", "")
                adjustments.append(GenotypeAdjustmentResult(
                    biomarker=adj.get("biomarker", ""),
                    standard_range=f"{std_range.get('lower', '')}-{std_range.get('upper', '')} {unit}".strip(),
                    adjusted_range=f"{adj_range.get('lower', '')}-{adj_range.get('upper', '')} {unit}".strip(),
                    genotype=adj.get("genotype_value", ""),
                    gene=adj.get("gene_display_name", ""),
                    rationale=adj.get("rationale", ""),
                ))
            except Exception as e:
                logger.warning(f"Skipping malformed adjustment result: {e}")

        # 5. Extract critical alerts
        critical_alerts = self._extract_critical_alerts(
            bio_age, trajectories, pgx_results,
        )

        # 6. Critical value threshold checking
        cv_alerts = self.critical_values.check(profile.biomarkers)
        for cv in cv_alerts:
            critical_alerts.append(cv.to_alert_string())

        # 7. Cross-biomarker discordance detection
        discordances = self.discordance.check(profile.biomarkers)
        for disc in discordances:
            critical_alerts.append(disc.to_alert_string())

        # 8. Lab range optimization (standard vs optimal)
        lab_comparisons = self.lab_ranges.get_discrepancies(
            profile.biomarkers, sex=profile.sex,
        )
        if lab_comparisons:
            for comp in lab_comparisons:
                critical_alerts.append(
                    f"OPTIMIZATION: {comp.to_interpretation()}"
                )

        # 9. Age-stratified reference range adjustments
        age_adjustments = self.adjuster.apply_age_adjustments(
            profile.biomarkers, profile.age, profile.sex,
        )
        for aa in age_adjustments:
            if aa["flag"] != "NORMAL":
                critical_alerts.append(
                    f"AGE-ADJUSTED: {aa['biomarker'].upper()} {aa['value']} is {aa['flag']} "
                    f"for age {profile.age} (age-adjusted range: "
                    f"{aa['age_adjusted_range']['low']}-{aa['age_adjusted_range']['high']}). "
                    f"{aa['note']}"
                )

        return AnalysisResult(
            patient_profile=profile,
            biological_age=bio_age,
            disease_trajectories=trajectories,
            pgx_results=pgx_results,
            genotype_adjustments=adjustments,
            critical_alerts=critical_alerts,
        )

    def search_plan(self, question: str) -> SearchPlan:
        """Analyze a question and determine search strategy.

        Identifies disease areas, relevant analysis modules, and
        decomposes complex questions into sub-queries.

        Args:
            question: The user's question.

        Returns:
            SearchPlan with identified topics and strategy.
        """
        plan = SearchPlan(question=question)
        q_lower = question.lower()
        q_upper = question.upper()

        # Identify disease areas
        domain_keywords = {
            "diabetes": ["diabetes", "hba1c", "glucose", "insulin", "homa", "metabolic"],
            "cardiovascular": ["cardiovascular", "cardiac", "heart", "lipid", "cholesterol",
                               "lp(a)", "apob", "ldl", "cvd", "atherosclerosis"],
            "liver": ["liver", "hepatic", "nafld", "nash", "masld", "alt", "ast",
                       "ggt", "fib-4", "pnpla3", "fibrosis"],
            "thyroid": ["thyroid", "tsh", "t3", "t4", "hashimoto"],
            "iron": ["iron", "ferritin", "hemochromatosis", "hfe", "transferrin"],
            "nutritional": ["vitamin", "folate", "b12", "omega", "nutrient", "mthfr", "supplement"],
        }
        for area, keywords in domain_keywords.items():
            if any(kw in q_lower for kw in keywords):
                plan.disease_areas.append(area)
                plan.identified_topics.append(area)

        # Identify relevant analysis modules
        module_keywords = {
            "biological_age": ["biological age", "phenoage", "grimage", "aging",
                               "epigenetic clock", "age acceleration"],
            "disease_trajectory": ["trajectory", "risk", "progression", "pre-diabetes",
                                   "pre-symptomatic", "disease risk"],
            "pharmacogenomics": ["pgx", "pharmacogenomic", "cyp", "metabolizer",
                                 "star allele", "cpic", "drug interaction", "dpyd",
                                 "vkorc1", "slco1b1"],
            "genotype_adjustment": ["genotype adjust", "reference range", "personalized range",
                                    "genetic modifier"],
        }
        for module, keywords in module_keywords.items():
            if any(kw in q_lower for kw in keywords):
                plan.relevant_modules.append(module)
                plan.identified_topics.append(module)

        # Determine search strategy
        if len(plan.disease_areas) == 1 and len(plan.relevant_modules) <= 1:
            plan.search_strategy = "domain-specific"
        elif plan.relevant_modules:
            plan.search_strategy = "targeted"
        else:
            plan.search_strategy = "broad"

        # Decompose complex questions
        if "WHY" in q_upper and ("ELEVATED" in q_upper or "HIGH" in q_upper):
            plan.sub_questions = [
                "What genetic variants cause elevated biomarker levels?",
                "What lifestyle factors contribute to elevated biomarker levels?",
                "What medications affect biomarker levels?",
            ]
        elif "COMPARE" in q_upper or " VS " in q_upper:
            plan.sub_questions = [
                "What are the differences in clinical interpretation?",
                "What are the genotype-specific considerations?",
            ]
        elif "SUPPLEMENT" in q_upper or "TREATMENT" in q_upper:
            plan.sub_questions = [
                "What are the evidence-based interventions for this condition?",
                "What genetic factors affect treatment response?",
            ]

        return plan

    def evaluate_evidence(self, evidence: CrossCollectionResult) -> str:
        """Evaluate the quality and coverage of retrieved evidence.

        Returns:
            'sufficient', 'partial', or 'insufficient'.
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

    def _extract_critical_alerts(
        self,
        bio_age: BiologicalAgeResult,
        trajectories: List[DiseaseTrajectoryResult],
        pgx_results: List[PGxResult],
    ) -> List[str]:
        """Extract critical findings requiring immediate attention.

        Checks for:
        - Severe age acceleration (>5 years)
        - Critical/high disease risk trajectories
        - PGx findings with safety implications (DPYD poor, CYP2D6 ultra-rapid on opioids)
        - Cross-modal triggers

        Args:
            bio_age: Biological age calculation result.
            trajectories: Disease trajectory analysis results.
            pgx_results: Pharmacogenomic mapping results.

        Returns:
            List of critical alert strings.
        """
        logger.info(f"Extracting critical alerts: trajectories={len(trajectories)}, pgx_results={len(pgx_results)}")

        alerts = []

        # Biological age alerts
        if bio_age.age_acceleration > 5:
            alerts.append(
                f"CRITICAL: Biological age acceleration of {bio_age.age_acceleration:.1f} years "
                f"(PhenoAge {bio_age.biological_age:.1f} vs chronological {bio_age.chronological_age}). "
                f"Indicates significantly accelerated aging."
            )

        # Disease trajectory alerts
        for traj in trajectories:
            if traj.risk_level in (RiskLevel.CRITICAL, RiskLevel.HIGH):
                alerts.append(
                    f"HIGH RISK: {traj.disease.value} trajectory at {traj.risk_level.value} level. "
                    + (f"Estimated {traj.years_to_onset_estimate:.0f} years to clinical onset. "
                       if traj.years_to_onset_estimate else "")
                    + "Immediate intervention recommended."
                )

        # PGx safety alerts
        for pgx in pgx_results:
            gene_upper = pgx.gene.upper()
            if gene_upper == "DPYD" and pgx.phenotype.value in ("poor", "intermediate"):
                alerts.append(
                    f"CRITICAL PGx: {pgx.gene} {pgx.star_alleles} -- {pgx.phenotype.value} metabolizer. "
                    f"Fluoropyrimidine (5-FU/capecitabine) toxicity risk. Contraindicated or dose reduction required."
                )
            elif gene_upper == "CYP2D6" and pgx.phenotype.value == "ultra_rapid":
                alerts.append(
                    f"PGx ALERT: {pgx.gene} {pgx.star_alleles} -- ultra-rapid metabolizer. "
                    f"Avoid codeine/tramadol (respiratory depression risk). Monitor opioid therapy."
                )
            elif gene_upper == "CYP2C19" and pgx.phenotype.value in ("poor", "intermediate"):
                # Check if clopidogrel is relevant
                alerts.append(
                    f"PGx ALERT: {pgx.gene} {pgx.star_alleles} -- {pgx.phenotype.value} metabolizer. "
                    f"Clopidogrel may be ineffective; consider prasugrel or ticagrelor."
                )

        return alerts

    def _build_enhanced_prompt(
        self,
        question: str,
        evidence: CrossCollectionResult,
        analysis: Optional[AnalysisResult],
    ) -> str:
        """Build enhanced LLM prompt incorporating RAG evidence and analysis results.

        Args:
            question: The user's question.
            evidence: Retrieved evidence from multi-collection search.
            analysis: Optional AnalysisResult from patient analysis modules.

        Returns:
            Complete prompt string for LLM generation.
        """
        # Start with standard RAG prompt
        prompt = self.rag._build_prompt(
            question, evidence,
            analysis.patient_profile if analysis else None,
        )

        # Append analysis results if available
        if analysis:
            analysis_sections = ["\n\n## Analysis Module Results\n"]

            # Biological age
            ba = analysis.biological_age
            analysis_sections.append(
                f"### Biological Age Assessment\n"
                f"- Chronological Age: {ba.chronological_age}\n"
                f"- Biological Age (PhenoAge): {ba.biological_age:.1f}\n"
                f"- Age Acceleration: {ba.age_acceleration:+.1f} years\n"
                f"- Mortality Risk Score: {ba.mortality_risk:.4f}\n"
            )
            if ba.aging_drivers:
                top_drivers = ba.aging_drivers[:3]
                driver_strs = [
                    f"{d.get('biomarker', 'unknown')}: {d.get('direction', '')} "
                    f"(contribution: {d.get('contribution', 0):.3f})"
                    for d in top_drivers
                ]
                analysis_sections.append(
                    f"- Top Aging Drivers: {'; '.join(driver_strs)}\n"
                )

            # Disease trajectories
            if analysis.disease_trajectories:
                analysis_sections.append("### Disease Trajectory Analysis\n")
                for traj in analysis.disease_trajectories:
                    analysis_sections.append(
                        f"- **{traj.disease.value}**: Risk Level = {traj.risk_level.value}"
                    )
                    if traj.years_to_onset_estimate:
                        analysis_sections.append(
                            f"  (est. {traj.years_to_onset_estimate:.0f} years to onset)"
                        )
                    if traj.genetic_risk_factors:
                        analysis_sections.append(
                            f"  Genetic factors: {', '.join(traj.genetic_risk_factors)}"
                        )
                analysis_sections.append("")

            # PGx results
            if analysis.pgx_results:
                analysis_sections.append("### Pharmacogenomic Profile\n")
                for pgx in analysis.pgx_results:
                    analysis_sections.append(
                        f"- **{pgx.gene}** {pgx.star_alleles}: "
                        f"{pgx.phenotype.value} metabolizer"
                    )
                    if pgx.drugs_affected:
                        for drug_info in pgx.drugs_affected[:3]:
                            analysis_sections.append(
                                f"  - {drug_info.get('drug', '')}: "
                                f"{drug_info.get('recommendation', '')}"
                            )
                analysis_sections.append("")

            # Genotype adjustments
            if analysis.genotype_adjustments:
                analysis_sections.append("### Genotype-Adjusted Reference Ranges\n")
                for adj in analysis.genotype_adjustments:
                    analysis_sections.append(
                        f"- **{adj.biomarker}** ({adj.gene} {adj.genotype}): "
                        f"Standard {adj.standard_range} -> Adjusted {adj.adjusted_range}"
                    )
                analysis_sections.append("")

            # Critical alerts
            if analysis.critical_alerts:
                analysis_sections.append("### Critical Alerts\n")
                for alert in analysis.critical_alerts:
                    analysis_sections.append(f"- {alert}")
                analysis_sections.append("")

            prompt += "\n".join(analysis_sections)

        return prompt
