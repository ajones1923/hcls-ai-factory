"""
Autoimmune Intelligence Agent -- Main orchestrator.

Coordinates autoimmune disease analysis:
1. Autoantibody panel interpretation
2. HLA-disease association analysis
3. Disease activity scoring
4. Flare prediction
5. Biologic therapy recommendation with PGx context
6. Cross-agent integration (Biomarker -> inflammation, Imaging -> joints)

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

from typing import Any, Dict, List

from loguru import logger

from .knowledge import (
    AUTOANTIBODY_DISEASE_MAP,
    BIOLOGIC_THERAPIES,
    DISEASE_ACTIVITY_THRESHOLDS,
    FLARE_BIOMARKER_PATTERNS,
    HLA_DISEASE_ASSOCIATIONS,
)
from .models import (
    AutoantibodyPanel,
    AutoimmuneAnalysisResult,
    AutoimmunePatientProfile,
    BiologicTherapy,
    DiseaseActivityLevel,
    DiseaseActivityScore,
    FlarePredictor,
    FlareRisk,
    HLAProfile,
)


class AutoimmuneAgent:
    """Autoimmune Intelligence Agent.

    Analyzes autoimmune disease risk, activity, and treatment options
    using autoantibody panels, HLA typing, biomarkers, and genomics.
    Integrates with the Biomarker Agent for inflammation monitoring
    and the Imaging Agent for joint/organ assessment.
    """

    def __init__(self, settings=None):
        self._settings = settings
        logger.info("Autoimmune Intelligence Agent initialized")

    def analyze_patient(self, profile: AutoimmunePatientProfile) -> AutoimmuneAnalysisResult:
        """Run complete autoimmune analysis pipeline.

        Steps:
        1. Interpret autoantibody panel -> disease associations
        2. Analyze HLA profile -> genetic susceptibility
        3. Calculate disease activity scores
        4. Predict flare risk
        5. Recommend biologic therapies with PGx filtering
        """
        logger.info(f"Analyzing autoimmune profile for patient {profile.patient_id}")

        result = AutoimmuneAnalysisResult(patient_id=profile.patient_id)

        # 1. Autoantibody interpretation
        if profile.autoantibody_panel:
            ab_findings = self.interpret_autoantibodies(profile.autoantibody_panel)
            result.cross_agent_findings.extend(ab_findings)

        # 2. HLA association analysis
        if profile.hla_profile:
            hla_assocs = self.analyze_hla_associations(profile.hla_profile)
            result.hla_associations = hla_assocs

        # 3. Disease activity scoring
        if profile.biomarkers and profile.diagnosed_conditions:
            scores = self.calculate_disease_activity(
                profile.biomarkers, profile.diagnosed_conditions
            )
            result.disease_activity_scores = scores

        # 4. Flare prediction
        if profile.biomarkers and profile.diagnosed_conditions:
            predictions = self.predict_flares(
                profile.biomarkers, profile.diagnosed_conditions
            )
            result.flare_predictions = predictions

        # 5. Biologic therapy recommendations
        if profile.diagnosed_conditions:
            recommendations = self.recommend_biologics(
                profile.diagnosed_conditions, profile.genotypes
            )
            result.biologic_recommendations = recommendations

        # 6. Generate critical alerts
        result.critical_alerts = self._generate_alerts(result)

        logger.info(
            f"Analysis complete: {len(result.disease_activity_scores)} activity scores, "
            f"{len(result.flare_predictions)} flare predictions, "
            f"{len(result.biologic_recommendations)} therapy recommendations, "
            f"{len(result.critical_alerts)} critical alerts"
        )

        return result

    def interpret_autoantibodies(self, panel: AutoantibodyPanel) -> List[Dict[str, Any]]:
        """Interpret autoantibody panel against disease associations."""
        findings = []

        for result in panel.results:
            if not result.positive:
                continue

            associations = AUTOANTIBODY_DISEASE_MAP.get(result.antibody, [])
            for assoc in associations:
                findings.append({
                    "antibody": result.antibody,
                    "disease": assoc["disease"],
                    "sensitivity": assoc.get("sensitivity", 0),
                    "specificity": assoc.get("specificity", 0),
                    "value": result.value,
                    "titer": result.titer,
                    "pattern": result.pattern,
                    "note": assoc.get("note", ""),
                })

        return findings

    def analyze_hla_associations(self, hla_profile: HLAProfile) -> List[Dict[str, Any]]:
        """Analyze HLA alleles for disease susceptibility."""
        associations = []

        for allele in hla_profile.all_alleles:
            # Check exact match and partial match (e.g., B*27:05 matches B*27)
            matches = HLA_DISEASE_ASSOCIATIONS.get(f"HLA-{allele}", [])

            # Also check broader allele groups
            if not matches and ":" in allele:
                broad = allele.split(":")[0]
                for key, assocs in HLA_DISEASE_ASSOCIATIONS.items():
                    if broad in key:
                        matches.extend(assocs)

            for match in matches:
                associations.append({
                    "allele": allele,
                    "disease": match["disease"],
                    "odds_ratio": match["odds_ratio"],
                    "pmid": match.get("pmid", ""),
                    "note": match.get("note", ""),
                })

        # Sort by odds ratio (highest risk first)
        associations.sort(key=lambda x: x["odds_ratio"], reverse=True)

        return associations

    def calculate_disease_activity(
        self,
        biomarkers: Dict[str, float],
        conditions: list,
    ) -> List[DiseaseActivityScore]:
        """Calculate disease activity scores for diagnosed conditions."""
        scores = []

        for condition in conditions:
            condition_str = condition.value if hasattr(condition, "value") else str(condition)

            # Find applicable scoring systems
            for score_name, score_info in DISEASE_ACTIVITY_THRESHOLDS.items():
                if score_info["disease"] != condition_str:
                    continue

                # Check if we have enough biomarker data
                # Simple CRP-based scoring as fallback
                crp = biomarkers.get("CRP", biomarkers.get("crp", None))
                esr = biomarkers.get("ESR", biomarkers.get("esr", None))

                if crp is not None or esr is not None:
                    thresholds = score_info["thresholds"]

                    # Simplified activity estimation from inflammatory markers
                    marker_value = crp if crp is not None else (esr / 10 if esr else 0)

                    if marker_value < thresholds.get("remission", 0):
                        level = DiseaseActivityLevel.REMISSION
                    elif marker_value < thresholds.get("low", 0):
                        level = DiseaseActivityLevel.LOW
                    elif marker_value < thresholds.get("moderate", 0):
                        level = DiseaseActivityLevel.MODERATE
                    else:
                        level = DiseaseActivityLevel.HIGH

                    scores.append(DiseaseActivityScore(
                        disease=condition,
                        score_name=score_name,
                        score_value=round(marker_value, 2),
                        level=level,
                        components={"CRP": crp} if crp else {"ESR": esr},
                        thresholds=thresholds,
                    ))

        return scores

    def predict_flares(
        self,
        biomarkers: Dict[str, float],
        conditions: list,
    ) -> List[FlarePredictor]:
        """Predict flare risk based on biomarker patterns."""
        predictions = []

        for condition in conditions:
            condition_str = condition.value if hasattr(condition, "value") else str(condition)
            pattern = FLARE_BIOMARKER_PATTERNS.get(condition_str)

            if not pattern:
                continue

            contributing = []
            protective = []
            risk_score = 0.3  # Base risk

            # Check early warning biomarkers
            for marker in pattern.get("early_warning_biomarkers", []):
                value = biomarkers.get(marker, biomarkers.get(marker.lower()))
                if value is not None:
                    # Elevated inflammatory markers increase risk
                    if marker in ("CRP", "ESR", "IL-6", "calprotectin") and value > 5:
                        contributing.append(f"Elevated {marker}: {value}")
                        risk_score += 0.15
                    elif marker in ("complement_C3", "complement_C4") and value < 80:
                        contributing.append(f"Low {marker}: {value}")
                        risk_score += 0.15
                    elif marker == "albumin" and value < 3.5:
                        contributing.append(f"Low albumin: {value}")
                        risk_score += 0.1
                    else:
                        protective.append(f"Normal {marker}: {value}")

            # Clamp risk score
            risk_score = min(max(risk_score, 0.0), 1.0)

            # Determine risk level from configurable thresholds
            imminent_thresh = getattr(self._settings, 'FLARE_RISK_IMMINENT', 0.8) if self._settings else 0.8
            high_thresh = getattr(self._settings, 'FLARE_RISK_HIGH', 0.6) if self._settings else 0.6
            moderate_thresh = getattr(self._settings, 'FLARE_RISK_MODERATE', 0.4) if self._settings else 0.4

            if risk_score >= imminent_thresh:
                risk_level = FlareRisk.IMMINENT
            elif risk_score >= high_thresh:
                risk_level = FlareRisk.HIGH
            elif risk_score >= moderate_thresh:
                risk_level = FlareRisk.MODERATE
            else:
                risk_level = FlareRisk.LOW

            monitoring = [f"Repeat {m} in 2-4 weeks" for m in pattern["early_warning_biomarkers"][:3]]

            predictions.append(FlarePredictor(
                disease=condition,
                current_activity=DiseaseActivityLevel.MODERATE,  # Simplified
                predicted_risk=risk_level,
                risk_score=round(risk_score, 3),
                contributing_factors=contributing,
                protective_factors=protective,
                recommended_monitoring=monitoring,
            ))

        return predictions

    def recommend_biologics(
        self,
        conditions: list,
        genotypes: Dict[str, str] = None,
    ) -> List[BiologicTherapy]:
        """Recommend biologic therapies filtered by PGx considerations."""
        recommendations = []
        condition_strs = [
            c.value if hasattr(c, "value") else str(c) for c in conditions
        ]

        for therapy_data in BIOLOGIC_THERAPIES:
            indicated = therapy_data.get("indicated_diseases", [])
            if not any(c in indicated for c in condition_strs):
                continue

            therapy = BiologicTherapy(
                drug_name=therapy_data["drug_name"],
                drug_class=therapy_data["drug_class"],
                mechanism=therapy_data.get("mechanism", ""),
                indicated_diseases=conditions,
                pgx_considerations=therapy_data.get("pgx_considerations", []),
                contraindications=therapy_data.get("contraindications", []),
                monitoring_requirements=therapy_data.get("monitoring_requirements", []),
            )
            recommendations.append(therapy)

        return recommendations

    def _generate_alerts(self, result: AutoimmuneAnalysisResult) -> List[str]:
        """Generate critical clinical alerts from analysis results."""
        alerts = []

        # High disease activity
        for score in result.disease_activity_scores:
            if score.level in (DiseaseActivityLevel.HIGH, DiseaseActivityLevel.VERY_HIGH):
                alerts.append(
                    f"HIGH DISEASE ACTIVITY: {score.disease.value} "
                    f"({score.score_name} = {score.score_value})"
                )

        # Imminent flare risk
        for pred in result.flare_predictions:
            if pred.predicted_risk in (FlareRisk.HIGH, FlareRisk.IMMINENT):
                alerts.append(
                    f"FLARE RISK {pred.predicted_risk.value.upper()}: "
                    f"{pred.disease.value} (score: {pred.risk_score:.2f})"
                )

        # Strong HLA associations
        for assoc in result.hla_associations:
            if assoc.get("odds_ratio", 0) > 5:
                alerts.append(
                    f"STRONG HLA ASSOCIATION: {assoc['allele']} -> "
                    f"{assoc['disease']} (OR={assoc['odds_ratio']})"
                )

        return alerts

    # ── Cross-agent integration ────────────────────────────────────────

    def request_biomarker_context(
        self,
        patient_id: str,
        biomarker_names: List[str],
    ) -> Dict[str, Any]:
        """Request inflammation context from the Biomarker Agent.

        Integration point: calls Biomarker Agent API for longitudinal
        biomarker trends and PGx context.
        """
        # Stub: returns structured placeholder until Biomarker Agent API is available
        logger.info(f"Cross-agent: requesting biomarker context for {patient_id}: {biomarker_names}")
        return {
            "source": "biomarker_agent",
            "patient_id": patient_id,
            "status": "stub",
            "requested_markers": biomarker_names,
            "trends": {marker: {"trend": "unavailable", "note": "Biomarker Agent integration pending"} for marker in biomarker_names},
        }

    def request_imaging_context(
        self,
        patient_id: str,
        body_regions: List[str],
    ) -> Dict[str, Any]:
        """Request imaging findings from the Imaging Agent.

        Integration point: correlates HLA/antibodies with joint damage
        scoring and organ involvement on imaging.
        """
        logger.info(f"Cross-agent: requesting imaging context for {patient_id}: {body_regions}")
        return {
            "source": "imaging_agent",
            "patient_id": patient_id,
            "status": "stub",
            "requested_regions": body_regions,
            "findings": {region: {"status": "unavailable", "note": "Imaging Agent integration pending"} for region in body_regions},
        }

    def publish_diagnosis_event(
        self,
        patient_id: str,
        disease: str,
        confidence: float,
        supporting_evidence: List[str],
    ) -> Dict[str, Any]:
        """Publish a diagnosis event for other agents to consume.

        Event-driven architecture: when a new autoimmune diagnosis is made,
        this event notifies Biomarker, Imaging, and Oncology agents to
        update their analysis context.
        """
        event = {
            "event_type": "autoimmune_diagnosis",
            "source_agent": "precision_autoimmune",
            "patient_id": patient_id,
            "disease": disease,
            "confidence": confidence,
            "supporting_evidence": supporting_evidence,
            "timestamp": None,  # Set by event bus
        }
        logger.info(f"Cross-agent: publishing diagnosis event for {patient_id}: {disease} (confidence={confidence:.2f})")
        # Stub: in production, this would publish to an event bus (Redis, Kafka, etc.)
        return {"status": "published_stub", "event": event}

    def enrich_analysis_with_cross_agent(
        self,
        result: AutoimmuneAnalysisResult,
        patient_id: str,
    ) -> AutoimmuneAnalysisResult:
        """Enrich analysis results with cross-agent data.

        Calls Biomarker and Imaging agents to add context to the analysis.
        """
        # Request biomarker trends for contributing factors
        biomarker_names = []
        for pred in result.flare_predictions:
            biomarker_names.extend([f.split(":")[0].replace("Elevated ", "").replace("Low ", "").strip()
                                     for f in pred.contributing_factors])

        if biomarker_names:
            bio_context = self.request_biomarker_context(patient_id, list(set(biomarker_names)))
            result.cross_agent_findings.append(bio_context)

        # Request imaging for diseases with joint/organ involvement
        imaging_diseases = {"rheumatoid_arthritis": ["hands", "feet", "knees"],
                           "ankylosing_spondylitis": ["spine", "sacroiliac_joints"],
                           "systemic_sclerosis": ["lungs", "heart"],
                           "systemic_lupus_erythematosus": ["kidneys", "joints"]}

        for score in result.disease_activity_scores:
            disease_str = score.disease.value if hasattr(score.disease, "value") else str(score.disease)
            regions = imaging_diseases.get(disease_str, [])
            if regions:
                img_context = self.request_imaging_context(patient_id, regions)
                result.cross_agent_findings.append(img_context)
                break  # Only request once

        return result
