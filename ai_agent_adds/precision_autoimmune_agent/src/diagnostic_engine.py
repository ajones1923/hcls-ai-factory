"""
Precision Autoimmune Agent — Diagnostic Engine

Implements clinical decision-support logic:
- ACR/EULAR classification criteria evaluation
- Diagnostic odyssey timeline analysis
- Differential diagnosis scoring
- Cross-disease overlap detection
- Autoantibody pattern recognition

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from .knowledge import (
    AUTOANTIBODY_DISEASE_MAP,
    HLA_DISEASE_ASSOCIATIONS,
)
from .models import (
    AutoimmuneDisease,
)

# ── ACR/EULAR Classification Criteria ────────────────────────────────────

# Simplified criteria point systems for supported diseases
CLASSIFICATION_CRITERIA = {
    AutoimmuneDisease.RHEUMATOID_ARTHRITIS: {
        "name": "2010 ACR/EULAR RA Classification",
        "threshold": 6,  # ≥6 points = definite RA
        "criteria": {
            "joint_involvement": {
                "1_large": 0, "2-10_large": 1, "1-3_small": 2,
                "4-10_small": 3, ">10_at_least_1_small": 5,
            },
            "serology": {
                "negative_RF_negative_CCP": 0, "low_positive_RF_or_CCP": 2,
                "high_positive_RF_or_CCP": 3,
            },
            "acute_phase": {
                "normal_CRP_normal_ESR": 0, "abnormal_CRP_or_ESR": 1,
            },
            "duration": {
                "<6_weeks": 0, ">=6_weeks": 1,
            },
        },
    },
    AutoimmuneDisease.SYSTEMIC_LUPUS: {
        "name": "2019 ACR/EULAR SLE Classification",
        "threshold": 10,  # ≥10 points + positive ANA
        "entry_criterion": "ANA ≥ 1:80",
        "criteria": {
            "constitutional": {"fever": 2},
            "hematologic": {
                "leukopenia": 3, "thrombocytopenia": 4, "autoimmune_hemolysis": 4,
            },
            "neuropsychiatric": {"delirium": 2, "psychosis": 3, "seizure": 5},
            "mucocutaneous": {
                "alopecia": 2, "oral_ulcers": 2, "subacute_cutaneous": 4,
                "acute_cutaneous": 6, "discoid": 4,
            },
            "serosal": {"pleural_or_pericardial": 5, "acute_pericarditis": 6},
            "musculoskeletal": {"joint_involvement": 6},
            "renal": {
                "proteinuria_>0.5g": 4, "class_II_V_nephritis": 8, "class_III_IV_nephritis": 10,
            },
            "immunology": {
                "anti_dsDNA_or_anti_Smith": 6, "low_complement": 3, "antiphospholipid": 2,
            },
        },
    },
    AutoimmuneDisease.ANKYLOSING_SPONDYLITIS: {
        "name": "ASAS Axial SpA Classification",
        "threshold": 1,  # imaging arm OR clinical arm
        "entry_criterion": "Back pain ≥3 months, onset <45 years",
        "criteria": {
            "imaging_arm": {"sacroiliitis_on_imaging_plus_1_SpA_feature": 1},
            "clinical_arm": {"HLA_B27_plus_2_SpA_features": 1},
            "spa_features": [
                "inflammatory_back_pain", "arthritis", "enthesitis", "uveitis",
                "dactylitis", "psoriasis", "IBD", "good_NSAID_response",
                "family_history_SpA", "HLA_B27", "elevated_CRP",
            ],
        },
    },
    AutoimmuneDisease.SYSTEMIC_SCLEROSIS: {
        "name": "2013 ACR/EULAR Systemic Sclerosis Classification",
        "threshold": 9,
        "criteria": {
            "skin_thickening": {
                "skin_thickening_fingers_extending_proximal_to_MCP": 9,
                "puffy_fingers": 2,
                "skin_thickening_fingers_distal_to_MCP": 4,
            },
            "fingertip_lesions": {
                "digital_tip_ulcers": 2,
                "fingertip_pitting_scars": 3,
            },
            "telangiectasia": {"telangiectasia": 2},
            "abnormal_nailfold_capillaries": {"abnormal_nailfold_capillaries": 2},
            "pulmonary": {
                "PAH": 2,
                "ILD": 2,
            },
            "raynauds": {"raynauds_phenomenon": 3},
            "ssc_autoantibodies": {
                "anti_centromere": 3,
                "anti_scl_70": 3,
                "anti_rna_polymerase_iii": 3,
            },
        },
    },
    AutoimmuneDisease.SJOGRENS_SYNDROME: {
        "name": "2016 ACR/EULAR Sjögren's Syndrome Classification",
        "threshold": 4,
        "entry_criterion": "At least one symptom of ocular or oral dryness",
        "criteria": {
            "histopathology": {
                "focal_lymphocytic_sialadenitis_focus_score_ge_1": 3,
            },
            "serology": {
                "anti_SSA_Ro_positive": 3,
            },
            "ocular": {
                "ocular_staining_score_ge_5": 1,
                "schirmer_test_le_5mm_5min": 1,
            },
            "salivary": {
                "unstimulated_salivary_flow_le_0_1_mL_min": 1,
            },
        },
    },
    AutoimmuneDisease.MULTIPLE_SCLEROSIS: {
        "name": "2017 McDonald Criteria for MS",
        "threshold": 2,
        "entry_criterion": "At least 1 clinical attack (CIS) or progressive neurological dysfunction",
        "criteria": {
            "dissemination_in_space": {
                "ge_1_T2_lesion_in_ge_2_MS_typical_CNS_regions": 1,
            },
            "dissemination_in_time": {
                "new_T2_or_gadolinium_lesion_on_followup_MRI": 1,
                "simultaneous_gadolinium_enhancing_and_nonenhancing_lesions": 1,
                "CSF_oligoclonal_bands": 1,
            },
        },
    },
    AutoimmuneDisease.MYASTHENIA_GRAVIS: {
        "name": "Clinical Diagnostic Criteria for Myasthenia Gravis",
        "threshold": 3,
        "criteria": {
            "clinical": {
                "fluctuating_weakness_ocular_bulbar_limb": 1,
                "fatigable_weakness_on_examination": 1,
            },
            "serological": {
                "AChR_antibody_positive": 2,
                "MuSK_antibody_positive": 2,
            },
            "electrophysiology": {
                "decremental_response_on_RNS": 1,
                "increased_jitter_on_SFEMG": 1,
            },
            "pharmacological": {
                "positive_edrophonium_test": 1,
                "positive_ice_pack_test": 1,
            },
        },
    },
    AutoimmuneDisease.CELIAC_DISEASE: {
        "name": "ESPGHAN Celiac Disease Diagnostic Criteria",
        "threshold": 3,
        "criteria": {
            "serology": {
                "anti_tTG_IgA_positive": 2,
                "anti_tTG_IgA_ge_10x_ULN": 3,
                "anti_endomysial_antibody_positive": 2,
            },
            "histology": {
                "villous_atrophy_marsh_3": 2,
                "crypt_hyperplasia_marsh_2": 1,
            },
            "genetics": {
                "HLA_DQ2_or_DQ8_positive": 1,
            },
            "clinical_response": {
                "symptomatic_improvement_on_GFD": 1,
            },
        },
    },
    AutoimmuneDisease.INFLAMMATORY_BOWEL: {
        "name": "Montreal Classification for IBD",
        "threshold": 3,
        "criteria": {
            "clinical": {
                "chronic_diarrhea_ge_4_weeks": 1,
                "bloody_stools": 1,
                "abdominal_pain": 1,
            },
            "endoscopic": {
                "mucosal_inflammation_on_colonoscopy": 2,
                "skip_lesions_crohns": 1,
                "continuous_inflammation_UC": 1,
            },
            "histological": {
                "chronic_inflammatory_infiltrate": 1,
                "granulomas_crohns": 2,
                "crypt_architecture_distortion": 1,
            },
            "imaging": {
                "bowel_wall_thickening_on_CT_MRI": 1,
                "strictures_or_fistulae": 1,
            },
        },
    },
    AutoimmuneDisease.PSORIASIS: {
        "name": "Clinical Psoriasis Diagnostic Criteria",
        "threshold": 3,
        "criteria": {
            "clinical_features": {
                "well_demarcated_erythematous_plaques": 2,
                "silvery_white_scale": 2,
                "auspitz_sign_positive": 1,
                "koebner_phenomenon": 1,
            },
            "distribution": {
                "extensor_surfaces_scalp_nails": 1,
                "nail_pitting_or_onycholysis": 1,
            },
            "histology": {
                "epidermal_hyperplasia_parakeratosis": 1,
                "munro_microabscesses": 1,
            },
        },
    },
}

# Cross-disease overlap syndromes
OVERLAP_SYNDROMES = {
    "mixed_connective_tissue_disease": {
        "required": ["anti-RNP"],
        "features_from": [
            AutoimmuneDisease.SYSTEMIC_LUPUS,
            AutoimmuneDisease.SYSTEMIC_SCLEROSIS,
            AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
        ],
        "min_features": 2,
    },
    "overlap_syndrome_sle_ra": {
        "diseases": [AutoimmuneDisease.SYSTEMIC_LUPUS, AutoimmuneDisease.RHEUMATOID_ARTHRITIS],
        "shared_markers": ["anti-CCP", "ANA", "RF"],
    },
    "sjogrens_sle_overlap": {
        "diseases": [AutoimmuneDisease.SJOGRENS_SYNDROME, AutoimmuneDisease.SYSTEMIC_LUPUS],
        "shared_markers": ["anti-SSA/Ro", "ANA", "anti-dsDNA"],
    },
    "pots_eds_mcas_triad": {
        "components": ["POTS", "hEDS", "MCAS"],
        "diagnostic_markers": ["tilt_table_positive", "beighton_score_>=5", "tryptase_elevated"],
    },
    "ssc_myositis_overlap": {
        "required": ["anti-Pm-Scl"],
        "features_from": [
            AutoimmuneDisease.SYSTEMIC_SCLEROSIS,
        ],
        "additional_features": ["myositis", "ILD", "mechanic_hands"],
        "min_features": 2,
    },
    "t1d_celiac_overlap": {
        "diseases": [AutoimmuneDisease.TYPE_1_DIABETES, AutoimmuneDisease.CELIAC_DISEASE],
        "shared_markers": ["anti-tTG IgA", "HLA-DQ2", "HLA-DQ8"],
        "prevalence": "5-10% of T1D patients have celiac disease",
    },
    "autoimmune_thyroid_t1d_overlap": {
        "diseases": [AutoimmuneDisease.HASHIMOTO_THYROIDITIS, AutoimmuneDisease.TYPE_1_DIABETES],
        "shared_markers": ["anti-TPO", "anti-GAD"],
        "prevalence": "15-30% of T1D patients develop autoimmune thyroiditis (APS-2)",
    },
    "ra_sjogrens_overlap": {
        "diseases": [AutoimmuneDisease.RHEUMATOID_ARTHRITIS, AutoimmuneDisease.SJOGRENS_SYNDROME],
        "shared_markers": ["RF", "anti-SSA/Ro", "ANA"],
        "prevalence": "Up to 31% of RA patients have secondary Sjogrens",
    },
    "lupus_aps_overlap": {
        "diseases": [AutoimmuneDisease.SYSTEMIC_LUPUS],
        "shared_markers": ["anti-cardiolipin", "lupus anticoagulant", "anti-beta2-glycoprotein"],
        "prevalence": "30-40% of SLE patients have antiphospholipid antibodies",
        "additional_risk": "Thrombosis, pregnancy morbidity",
    },
}


class DiagnosticEngine:
    """Clinical diagnostic reasoning engine for autoimmune diseases."""

    def __init__(self, agent=None, rag_engine=None):
        self.agent = agent
        self.rag = rag_engine

    # ── Classification criteria evaluation ────────────────────────────
    def evaluate_classification_criteria(
        self,
        disease: AutoimmuneDisease,
        clinical_data: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Evaluate ACR/EULAR classification criteria for a disease."""
        criteria_set = CLASSIFICATION_CRITERIA.get(disease)
        if not criteria_set:
            return {"disease": disease.value, "supported": False, "message": "No criteria defined"}

        total_points = 0
        met_criteria = []
        unmet_criteria = []

        criteria = criteria_set["criteria"]
        for category, items in criteria.items():
            if isinstance(items, dict):
                for item, points in items.items():
                    if clinical_data.get(item):
                        total_points += points
                        met_criteria.append(f"{category}: {item} (+{points})")
                    else:
                        unmet_criteria.append(f"{category}: {item} ({points} pts)")

        threshold = criteria_set["threshold"]
        meets_criteria = total_points >= threshold

        return {
            "disease": disease.value,
            "criteria_set": criteria_set["name"],
            "total_points": total_points,
            "threshold": threshold,
            "meets_criteria": meets_criteria,
            "met_criteria": met_criteria,
            "unmet_criteria": unmet_criteria,
            "entry_criterion": criteria_set.get("entry_criterion"),
        }

    # ── Differential diagnosis ────────────────────────────────────────
    def generate_differential(
        self,
        positive_antibodies: List[str],
        hla_alleles: Optional[List[str]] = None,
        symptoms: Optional[List[str]] = None,
        biomarkers: Optional[Dict[str, float]] = None,
    ) -> List[Dict[str, Any]]:
        """Generate ranked differential diagnosis from available data."""
        disease_scores: Dict[str, float] = {}
        evidence: Dict[str, List[str]] = {}

        # Score from autoantibodies
        for ab in positive_antibodies:
            associations = AUTOANTIBODY_DISEASE_MAP.get(ab, [])
            for assoc in associations:
                disease = assoc["disease"]
                # Weight by specificity (high specificity = more diagnostic value)
                score = assoc.get("specificity", 0.5) * 2.0
                disease_scores[disease] = disease_scores.get(disease, 0) + score
                evidence.setdefault(disease, []).append(
                    f"Positive {ab} (spec={assoc.get('specificity', '?')})"
                )

        # Score from HLA alleles
        if hla_alleles:
            for allele in hla_alleles:
                key = f"HLA-{allele}"
                assocs = HLA_DISEASE_ASSOCIATIONS.get(key, [])
                for assoc in assocs:
                    disease = assoc["disease"]
                    # Log-scale odds ratio contribution
                    import math
                    or_score = math.log2(max(assoc["odds_ratio"], 1.0)) * 0.5
                    disease_scores[disease] = disease_scores.get(disease, 0) + or_score
                    evidence.setdefault(disease, []).append(
                        f"HLA {allele} (OR={assoc['odds_ratio']})"
                    )

        # Sort by score
        ranked = sorted(disease_scores.items(), key=lambda x: x[1], reverse=True)

        return [
            {
                "disease": disease,
                "score": round(score, 2),
                "evidence": evidence.get(disease, []),
                "rank": i + 1,
            }
            for i, (disease, score) in enumerate(ranked)
        ]

    # ── Diagnostic odyssey analysis ───────────────────────────────────
    def analyze_diagnostic_odyssey(
        self,
        timeline_events: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Analyze a patient's diagnostic odyssey from timeline events.

        Identifies:
        - Time from first symptom to diagnosis
        - Number of specialists seen
        - Misdiagnoses
        - Key turning points
        - Missed opportunities
        """
        if not timeline_events:
            return {"status": "no_data"}

        # Sort by date
        events = sorted(timeline_events, key=lambda e: e.get("event_date", ""))

        first_symptom = None
        diagnosis_date = None
        specialists = set()
        misdiagnoses = []
        key_tests = []

        for event in events:
            event_type = event.get("event_type", "")
            specialty = event.get("specialty", "")
            description = event.get("description", "")

            if specialty:
                specialists.add(specialty)

            if event_type == "symptom_onset" and first_symptom is None:
                first_symptom = event.get("event_date")

            if event_type == "diagnosis":
                diagnosis_date = event.get("event_date")

            if event_type == "misdiagnosis":
                misdiagnoses.append({
                    "date": event.get("event_date"),
                    "wrong_diagnosis": description,
                    "provider": event.get("provider"),
                })

            if event_type in ("lab_result", "imaging", "biopsy", "genetic_test"):
                key_tests.append({
                    "date": event.get("event_date"),
                    "test": description,
                    "type": event_type,
                })

        # Calculate diagnostic delay
        delay_info = None
        if first_symptom and diagnosis_date:
            try:
                from datetime import datetime
                fmt = "%Y-%m-%d"
                d1 = datetime.strptime(first_symptom[:10], fmt)
                d2 = datetime.strptime(diagnosis_date[:10], fmt)
                delay_days = (d2 - d1).days
                delay_info = {
                    "days": delay_days,
                    "months": round(delay_days / 30.44, 1),
                    "years": round(delay_days / 365.25, 1),
                }
            except (ValueError, TypeError):
                pass

        return {
            "first_symptom_date": first_symptom,
            "diagnosis_date": diagnosis_date,
            "diagnostic_delay": delay_info,
            "num_specialists_seen": len(specialists),
            "specialists": sorted(specialists),
            "num_misdiagnoses": len(misdiagnoses),
            "misdiagnoses": misdiagnoses,
            "key_diagnostic_tests": key_tests,
            "total_events": len(events),
        }

    # ── Cross-disease overlap detection ───────────────────────────────
    def detect_overlap_syndromes(
        self,
        positive_antibodies: List[str],
        diagnosed_conditions: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """Detect potential overlap syndromes from antibody profile."""
        detected = []

        for syndrome_name, criteria in OVERLAP_SYNDROMES.items():

            # Check required antibodies
            required = criteria.get("required", [])
            if required and not all(ab in positive_antibodies for ab in required):
                continue

            # Check shared markers
            shared = criteria.get("shared_markers", [])
            marker_matches = [m for m in shared if m in positive_antibodies]

            # Check diagnostic markers
            criteria.get("diagnostic_markers", [])

            if marker_matches or required:
                detected.append({
                    "syndrome": syndrome_name,
                    "matched_markers": marker_matches or required,
                    "diseases_involved": [
                        d.value if hasattr(d, "value") else str(d)
                        for d in criteria.get("diseases", criteria.get("features_from", []))
                    ],
                    "confidence": "high" if len(marker_matches) >= 2 else "moderate",
                })

        return detected
