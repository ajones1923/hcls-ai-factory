"""Validated clinical scale calculators for the Neurology Intelligence Agent.

Implements 10 evidence-based neurological scoring instruments with
automated interpretation, threshold identification, and clinical
recommendations.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional

from src.models import ClinicalScaleType, ScaleResult

logger = logging.getLogger(__name__)


# ===================================================================
# 1. NIHSS Calculator  (NIH Stroke Scale)
# ===================================================================


class NIHSSCalculator:
    """NIH Stroke Scale — 15 items, score 0-42.

    Quantifies stroke severity for tPA eligibility and prognosis.

    Items
    -----
    1a_loc          : Level of consciousness (0-3)
    1b_loc_questions: LOC questions (0-2)
    1c_loc_commands : LOC commands (0-2)
    2_gaze          : Best gaze (0-2)
    3_visual        : Visual fields (0-3)
    4_facial        : Facial palsy (0-3)
    5a_left_arm     : Left arm motor (0-4)
    5b_right_arm    : Right arm motor (0-4)
    6a_left_leg     : Left leg motor (0-4)
    6b_right_leg    : Right leg motor (0-4)
    7_ataxia        : Limb ataxia (0-2)
    8_sensory       : Sensory (0-2)
    9_language      : Best language (0-3)
    10_dysarthria   : Dysarthria (0-2)
    11_extinction   : Extinction/inattention (0-2)
    """

    ITEMS: Dict[str, int] = {
        "1a_loc": 3,
        "1b_loc_questions": 2,
        "1c_loc_commands": 2,
        "2_gaze": 2,
        "3_visual": 3,
        "4_facial": 3,
        "5a_left_arm": 4,
        "5b_right_arm": 4,
        "6a_left_leg": 4,
        "6b_right_leg": 4,
        "7_ataxia": 2,
        "8_sensory": 2,
        "9_language": 3,
        "10_dysarthria": 2,
        "11_extinction": 2,
    }

    MAX_SCORE = 42

    @classmethod
    def calculate(cls, scores: Dict[str, int]) -> ScaleResult:
        """Calculate NIHSS total from item scores.

        Parameters
        ----------
        scores : dict
            Mapping of item names to integer scores.
        """
        total = 0
        for item, max_val in cls.ITEMS.items():
            val = scores.get(item, 0)
            val = max(0, min(val, max_val))
            total += val

        if total == 0:
            category = "No stroke symptoms"
            interpretation = "No measurable neurological deficit."
        elif total <= 4:
            category = "Minor stroke"
            interpretation = "Minor neurological deficit. Consider IV tPA if within window."
        elif total <= 15:
            category = "Moderate stroke"
            interpretation = (
                "Moderate deficit. IV tPA strongly recommended if within 4.5-hour "
                "window and no contraindications."
            )
        elif total <= 20:
            category = "Moderate-to-severe stroke"
            interpretation = (
                "Moderate-to-severe deficit. IV tPA indicated. Consider mechanical "
                "thrombectomy if large vessel occlusion."
            )
        elif total <= 25:
            category = "Severe stroke"
            interpretation = (
                "Severe deficit. IV tPA indicated. Strong consideration for "
                "endovascular thrombectomy (DAWN/DEFUSE-3 criteria)."
            )
        else:
            category = "Very severe stroke"
            interpretation = (
                "Very severe deficit. Evaluate for thrombectomy; high risk of "
                "hemorrhagic transformation with tPA. ICU-level care required."
            )

        recommendations = []
        if 1 <= total <= 25:
            recommendations.append(
                "Administer IV alteplase (0.9 mg/kg, max 90 mg) if within 4.5 hours "
                "of symptom onset and no contraindications"
            )
        if total >= 6:
            recommendations.append(
                "Obtain CTA head/neck to evaluate for large vessel occlusion"
            )
        if total >= 6:
            recommendations.append(
                "Consider mechanical thrombectomy if LVO confirmed and within "
                "treatment window (up to 24h per DAWN/DEFUSE-3)"
            )
        if total > 25:
            recommendations.append(
                "ICU admission with close neurological monitoring; "
                "neurosurgical consultation for possible decompressive craniectomy"
            )

        return ScaleResult(
            scale_type=ClinicalScaleType.NIHSS,
            score=total,
            max_score=cls.MAX_SCORE,
            interpretation=interpretation,
            severity_category=category,
            thresholds={
                "tpa_consideration": 1.0,
                "moderate_stroke": 5.0,
                "severe_stroke": 21.0,
                "thrombectomy_consideration": 6.0,
            },
            recommendations=recommendations,
        )


# ===================================================================
# 2. GCS Calculator  (Glasgow Coma Scale)
# ===================================================================


class GCSCalculator:
    """Glasgow Coma Scale — Eye + Verbal + Motor = 3-15.

    Parameters
    ----------
    eye : int (1-4)
        1=None, 2=To pressure, 3=To voice, 4=Spontaneous
    verbal : int (1-5)
        1=None, 2=Sounds, 3=Words, 4=Confused, 5=Oriented
    motor : int (1-6)
        1=None, 2=Extension, 3=Abnormal flexion, 4=Withdrawal,
        5=Localising, 6=Obeys commands
    """

    MAX_SCORE = 15

    @classmethod
    def calculate(cls, eye: int, verbal: int, motor: int) -> ScaleResult:
        eye = max(1, min(eye, 4))
        verbal = max(1, min(verbal, 5))
        motor = max(1, min(motor, 6))
        total = eye + verbal + motor

        if total <= 8:
            category = "Severe"
            interpretation = (
                "Severe brain injury. Intubation and airway protection likely "
                "required. Neurosurgical consultation recommended."
            )
        elif total <= 12:
            category = "Moderate"
            interpretation = (
                "Moderate brain injury. Close monitoring required. "
                "Consider CT head and neurosurgical evaluation."
            )
        else:
            category = "Mild"
            interpretation = (
                "Mild brain injury or normal. Monitor for deterioration."
            )

        recommendations = []
        if total <= 8:
            recommendations.extend([
                "Secure airway — consider intubation (GCS <= 8)",
                "Emergent CT head without contrast",
                "Neurosurgical consultation",
                "ICU admission with neuro-checks every 1-2 hours",
            ])
        elif total <= 12:
            recommendations.extend([
                "CT head without contrast",
                "Neurological observation with serial GCS assessments",
                "Neurosurgery consult if structural lesion identified",
            ])
        else:
            recommendations.append(
                "Serial neurological assessments; CT head if mechanism warrants"
            )

        return ScaleResult(
            scale_type=ClinicalScaleType.GCS,
            score=total,
            max_score=cls.MAX_SCORE,
            interpretation=interpretation,
            severity_category=category,
            thresholds={
                "intubation_threshold": 8.0,
                "moderate_injury": 9.0,
                "mild_injury": 13.0,
            },
            recommendations=recommendations,
        )


# ===================================================================
# 3. MoCA Calculator  (Montreal Cognitive Assessment)
# ===================================================================


class MoCACalculator:
    """Montreal Cognitive Assessment — 8 domains, score 0-30.

    Domains
    -------
    visuospatial   : Visuospatial/executive (0-5)
    naming         : Naming (0-3)
    attention      : Attention (0-6)
    language       : Language (0-3)
    abstraction    : Abstraction (0-2)
    delayed_recall : Delayed recall (0-5)
    orientation    : Orientation (0-6)
    education_adj  : +1 if <= 12 years education (0-1)
    """

    DOMAINS: Dict[str, int] = {
        "visuospatial": 5,
        "naming": 3,
        "attention": 6,
        "language": 3,
        "abstraction": 2,
        "delayed_recall": 5,
        "orientation": 6,
    }

    MAX_SCORE = 30

    @classmethod
    def calculate(
        cls, domain_scores: Dict[str, int], education_years: int = 13
    ) -> ScaleResult:
        total = 0
        for domain, max_val in cls.DOMAINS.items():
            val = domain_scores.get(domain, 0)
            val = max(0, min(val, max_val))
            total += val

        # Education adjustment: +1 point if <= 12 years of education
        if education_years <= 12:
            total = min(total + 1, cls.MAX_SCORE)

        if total >= 26:
            category = "Normal"
            interpretation = "Cognitive performance within normal limits."
        elif total >= 18:
            category = "Mild cognitive impairment"
            interpretation = (
                "Score suggests mild cognitive impairment (MCI). "
                "Further neuropsychological testing and biomarker evaluation recommended."
            )
        elif total >= 10:
            category = "Moderate cognitive impairment"
            interpretation = (
                "Score suggests moderate cognitive impairment consistent with dementia. "
                "Comprehensive evaluation with neuroimaging and biomarkers indicated."
            )
        else:
            category = "Severe cognitive impairment"
            interpretation = (
                "Score suggests severe cognitive impairment. "
                "Full dementia workup including structural MRI, FDG-PET, "
                "and CSF biomarkers indicated."
            )

        recommendations = []
        if total < 26:
            recommendations.append(
                "Refer for comprehensive neuropsychological testing"
            )
            recommendations.append(
                "Obtain structural MRI brain with volumetric analysis"
            )
        if total < 22:
            recommendations.append(
                "Consider amyloid PET or CSF AD biomarkers (Abeta42, p-tau, t-tau)"
            )
            recommendations.append(
                "Evaluate for reversible causes: B12, TSH, RPR, HIV, metabolic panel"
            )
        if total < 18:
            recommendations.append(
                "Assess functional status (ADLs/IADLs) and caregiver support needs"
            )
        if 18 <= total <= 25:
            recommendations.append(
                "Consider anti-amyloid therapy eligibility (lecanemab/donanemab) "
                "if amyloid-positive on PET or CSF"
            )

        # Flag impaired domains
        impaired = []
        for domain, max_val in cls.DOMAINS.items():
            val = domain_scores.get(domain, 0)
            if val < max_val * 0.5:
                impaired.append(domain)
        if impaired:
            recommendations.append(
                f"Notably impaired domains: {', '.join(impaired)}"
            )

        return ScaleResult(
            scale_type=ClinicalScaleType.MOCA,
            score=total,
            max_score=cls.MAX_SCORE,
            interpretation=interpretation,
            severity_category=category,
            thresholds={
                "normal_cutoff": 26.0,
                "mci_range_low": 18.0,
                "moderate_impairment": 10.0,
                "anti_amyloid_consideration": 18.0,
            },
            recommendations=recommendations,
        )


# ===================================================================
# 4. UPDRS Part III Calculator (Unified Parkinson's Disease Rating Scale)
# ===================================================================


class UPDRSCalculator:
    """MDS-UPDRS Part III Motor Examination — 33 scores, 0-132.

    18 items (some bilateral) rated 0-4 each. Total 33 sub-scores.
    """

    ITEMS: Dict[str, int] = {
        "3.1_speech": 4,
        "3.2_facial_expression": 4,
        "3.3a_rigidity_neck": 4,
        "3.3b_rigidity_rue": 4,
        "3.3c_rigidity_lue": 4,
        "3.3d_rigidity_rle": 4,
        "3.3e_rigidity_lle": 4,
        "3.4a_finger_tapping_r": 4,
        "3.4b_finger_tapping_l": 4,
        "3.5a_hand_movements_r": 4,
        "3.5b_hand_movements_l": 4,
        "3.6a_pronation_supination_r": 4,
        "3.6b_pronation_supination_l": 4,
        "3.7a_toe_tapping_r": 4,
        "3.7b_toe_tapping_l": 4,
        "3.8a_leg_agility_r": 4,
        "3.8b_leg_agility_l": 4,
        "3.9_arising_from_chair": 4,
        "3.10_gait": 4,
        "3.11_freezing_of_gait": 4,
        "3.12_postural_stability": 4,
        "3.13_posture": 4,
        "3.14_body_bradykinesia": 4,
        "3.15a_postural_tremor_r": 4,
        "3.15b_postural_tremor_l": 4,
        "3.16a_kinetic_tremor_r": 4,
        "3.16b_kinetic_tremor_l": 4,
        "3.17a_rest_tremor_rue": 4,
        "3.17b_rest_tremor_lue": 4,
        "3.17c_rest_tremor_rle": 4,
        "3.17d_rest_tremor_lle": 4,
        "3.17e_rest_tremor_jaw": 4,
        "3.18_constancy_of_rest_tremor": 4,
    }

    MAX_SCORE = 132

    @classmethod
    def calculate(cls, scores: Dict[str, int]) -> ScaleResult:
        total = 0
        for item, max_val in cls.ITEMS.items():
            val = scores.get(item, 0)
            val = max(0, min(val, max_val))
            total += val

        if total <= 10:
            category = "Minimal"
            interpretation = "Minimal motor signs. May not meet clinical threshold."
        elif total <= 32:
            category = "Mild"
            interpretation = "Mild motor impairment. Dopaminergic therapy may be considered."
        elif total <= 58:
            category = "Moderate"
            interpretation = (
                "Moderate motor impairment. Optimise dopaminergic therapy. "
                "Consider adjunctive agents."
            )
        elif total <= 80:
            category = "Severe"
            interpretation = (
                "Severe motor impairment despite medication. "
                "Evaluate for DBS candidacy (CAPSIT-PD criteria)."
            )
        else:
            category = "Very severe"
            interpretation = (
                "Very severe motor impairment. Consider advanced therapies: "
                "DBS, Duopa, apomorphine infusion."
            )

        recommendations = []
        if total >= 40:
            recommendations.append(
                "Review and optimise levodopa dosing and timing"
            )
        if total >= 59:
            recommendations.append(
                "Evaluate DBS candidacy: levodopa-responsive, no significant "
                "cognitive impairment (MoCA >= 26), age consideration"
            )
        if total >= 80:
            recommendations.append(
                "Consider advanced therapies: DBS (STN or GPi), "
                "LCIG (Duopa), subcutaneous apomorphine infusion"
            )

        return ScaleResult(
            scale_type=ClinicalScaleType.UPDRS,
            score=total,
            max_score=cls.MAX_SCORE,
            interpretation=interpretation,
            severity_category=category,
            thresholds={
                "mild_threshold": 11.0,
                "moderate_threshold": 33.0,
                "dbs_candidacy_threshold": 59.0,
                "advanced_therapy_threshold": 80.0,
            },
            recommendations=recommendations,
        )


# ===================================================================
# 5. EDSS Calculator  (Expanded Disability Status Scale)
# ===================================================================


class EDSSCalculator:
    """Kurtzke EDSS — 7 functional systems, score 0.0-10.0.

    Functional Systems
    ------------------
    visual      : Visual / optic (0-6)
    brainstem   : Brainstem (0-5)
    pyramidal   : Pyramidal / motor (0-6)
    cerebellar  : Cerebellar (0-5)
    sensory     : Sensory (0-6)
    bowel_bladder: Bowel & bladder (0-6)
    cerebral    : Cerebral / mental (0-5)
    ambulation  : Ambulatory status descriptor (string or float)
    """

    FS_MAXIMA: Dict[str, int] = {
        "visual": 6,
        "brainstem": 5,
        "pyramidal": 6,
        "cerebellar": 5,
        "sensory": 6,
        "bowel_bladder": 6,
        "cerebral": 5,
    }

    MAX_SCORE = 10.0

    @classmethod
    def calculate(
        cls, fs_scores: Dict[str, int], edss_step: float
    ) -> ScaleResult:
        """Calculate EDSS from functional system scores and overall EDSS step.

        The EDSS step (0.0-10.0 in 0.5 increments) is determined by the
        clinician based on FS scores plus ambulation. This calculator
        validates and interprets the score.
        """
        edss = max(0.0, min(edss_step, 10.0))
        # Round to nearest 0.5
        edss = round(edss * 2) / 2

        if edss <= 1.5:
            category = "Minimal disability"
            interpretation = "No or minimal disability. Fully ambulatory."
        elif edss <= 3.5:
            category = "Mild disability"
            interpretation = (
                "Mild-to-moderate disability. Fully ambulatory but with "
                "some functional system impairment."
            )
        elif edss <= 5.5:
            category = "Moderate disability"
            interpretation = (
                "Moderate disability with ambulatory limitation. "
                "May require assistive device for longer distances."
            )
        elif edss <= 6.5:
            category = "Walking aid required"
            interpretation = (
                "Requires unilateral or bilateral walking aid for ambulation. "
                "Consider high-efficacy DMT and rehabilitation."
            )
        elif edss <= 7.5:
            category = "Wheelchair dependent"
            interpretation = (
                "Unable to walk more than a few steps. Wheelchair is principal "
                "mode of mobility. Evaluate for advanced supportive care."
            )
        elif edss <= 9.0:
            category = "Restricted to bed"
            interpretation = (
                "Restricted to bed or chair. Dependent for most functions. "
                "Focus on symptom management and quality of life."
            )
        else:
            category = "Death due to MS"
            interpretation = "Death attributable to MS."

        recommendations = []
        if edss <= 6.5:
            recommendations.append(
                "Continue or escalate disease-modifying therapy"
            )
        if edss >= 4.0:
            recommendations.append(
                "Refer to neuro-rehabilitation and physical therapy"
            )
        if edss >= 6.0:
            recommendations.append(
                "Assess for symptomatic management: spasticity, bladder, fatigue"
            )

        # NEDA-3 progression: >= 1.0 point increase if baseline EDSS <= 5.5,
        # or >= 0.5 point increase if baseline EDSS >= 6.0
        return ScaleResult(
            scale_type=ClinicalScaleType.EDSS,
            score=edss,
            max_score=cls.MAX_SCORE,
            interpretation=interpretation,
            severity_category=category,
            thresholds={
                "neda3_progression_low_baseline": 1.0,
                "neda3_progression_high_baseline": 0.5,
                "walking_aid_threshold": 6.0,
                "wheelchair_threshold": 7.0,
            },
            recommendations=recommendations,
        )


# ===================================================================
# 6. mRS Calculator  (Modified Rankin Scale)
# ===================================================================


class mRSCalculator:
    """Modified Rankin Scale — global disability after stroke, 0-6.

    0 = No symptoms
    1 = No significant disability despite symptoms
    2 = Slight disability (unable to carry out all previous activities)
    3 = Moderate disability (requires some help, walks without assistance)
    4 = Moderately severe (unable to walk/attend to needs without assistance)
    5 = Severe disability (bedridden, incontinent, requires constant care)
    6 = Dead
    """

    MAX_SCORE = 6

    _DESCRIPTIONS = {
        0: ("No symptoms", "No residual symptoms at all."),
        1: (
            "No significant disability",
            "Despite symptoms, able to carry out all usual duties and activities.",
        ),
        2: (
            "Slight disability",
            "Unable to carry out all previous activities but able to look "
            "after own affairs without assistance.",
        ),
        3: (
            "Moderate disability",
            "Requires some help but able to walk without assistance.",
        ),
        4: (
            "Moderately severe disability",
            "Unable to walk without assistance and unable to attend to "
            "own bodily needs without assistance.",
        ),
        5: (
            "Severe disability",
            "Bedridden, incontinent, requires constant nursing care and attention.",
        ),
        6: ("Dead", "Death."),
    }

    @classmethod
    def calculate(cls, score: int) -> ScaleResult:
        score = max(0, min(score, 6))
        category, interpretation = cls._DESCRIPTIONS[score]

        recommendations = []
        if score <= 2:
            recommendations.append(
                "Good functional outcome. Focus on secondary stroke prevention."
            )
        elif score == 3:
            recommendations.append(
                "Moderate disability. Outpatient rehabilitation and "
                "secondary prevention indicated."
            )
        elif score >= 4 and score <= 5:
            recommendations.append(
                "Significant disability. Inpatient rehabilitation, "
                "caregiver support, and long-term care planning."
            )

        return ScaleResult(
            scale_type=ClinicalScaleType.MRS,
            score=score,
            max_score=cls.MAX_SCORE,
            interpretation=interpretation,
            severity_category=category,
            thresholds={
                "good_outcome": 2.0,
                "functional_independence": 3.0,
                "poor_outcome": 4.0,
            },
            recommendations=recommendations,
        )


# ===================================================================
# 7. HIT-6 Calculator  (Headache Impact Test)
# ===================================================================


class HIT6Calculator:
    """Headache Impact Test-6 — 6 items, score 36-78.

    Each item scored: Never=6, Rarely=8, Sometimes=10, Very Often=11, Always=13

    Items assess: pain severity, social functioning, role functioning,
    cognitive functioning, psychological distress, vitality.
    """

    ITEM_VALUES = {"never": 6, "rarely": 8, "sometimes": 10, "very_often": 11, "always": 13}
    MIN_SCORE = 36
    MAX_SCORE = 78

    @classmethod
    def calculate(cls, responses: List[str]) -> ScaleResult:
        """Calculate HIT-6 from 6 categorical responses.

        Parameters
        ----------
        responses : list of str
            Six responses, each one of: 'never', 'rarely', 'sometimes',
            'very_often', 'always'.
        """
        if len(responses) < 6:
            responses = responses + ["never"] * (6 - len(responses))

        total = sum(cls.ITEM_VALUES.get(r.lower(), 6) for r in responses[:6])
        total = max(cls.MIN_SCORE, min(total, cls.MAX_SCORE))

        if total <= 49:
            category = "Little or no impact"
            interpretation = "Headaches have little or no impact on daily life."
        elif total <= 55:
            category = "Some impact"
            interpretation = (
                "Headaches have some impact on daily life. "
                "Discuss acute treatment optimization."
            )
        elif total <= 59:
            category = "Substantial impact"
            interpretation = (
                "Headaches have a substantial impact on daily life. "
                "Preventive therapy should be strongly considered."
            )
        else:
            category = "Severe impact"
            interpretation = (
                "Headaches have a severe and disabling impact on daily life. "
                "Preventive therapy is indicated."
            )

        recommendations = []
        if total >= 50:
            recommendations.append(
                "Optimise acute treatment: consider triptans if not already prescribed"
            )
        if total >= 56:
            recommendations.append(
                "Initiate preventive therapy: first-line options include "
                "topiramate, propranolol, amitriptyline, or valproate"
            )
        if total >= 60:
            recommendations.append(
                "Consider CGRP monoclonal antibody therapy "
                "(erenumab, fremanezumab, galcanezumab) if >= 2 preventives have failed"
            )
            recommendations.append(
                "Evaluate for medication-overuse headache"
            )

        return ScaleResult(
            scale_type=ClinicalScaleType.HIT6,
            score=total,
            max_score=cls.MAX_SCORE,
            interpretation=interpretation,
            severity_category=category,
            thresholds={
                "some_impact": 50.0,
                "substantial_impact": 56.0,
                "severe_impact": 60.0,
                "cgrp_consideration": 60.0,
            },
            recommendations=recommendations,
        )


# ===================================================================
# 8. ALSFRS-R Calculator  (ALS Functional Rating Scale - Revised)
# ===================================================================


class ALSFRSCalculator:
    """ALS Functional Rating Scale - Revised — 12 items, score 0-48.

    Higher scores indicate better function. Rate of decline is prognostic.

    Items (each 0-4, where 4 = normal function):
    1. Speech
    2. Salivation
    3. Swallowing
    4. Handwriting
    5. Cutting food (with/without gastrostomy)
    6. Dressing and hygiene
    7. Turning in bed
    8. Walking
    9. Climbing stairs
    10. Dyspnea
    11. Orthopnea
    12. Respiratory insufficiency
    """

    ITEMS = [
        "speech", "salivation", "swallowing", "handwriting",
        "cutting_food", "dressing_hygiene", "turning_in_bed", "walking",
        "climbing_stairs", "dyspnea", "orthopnea", "respiratory_insufficiency",
    ]

    MAX_SCORE = 48

    @classmethod
    def calculate(
        cls,
        scores: Dict[str, int],
        months_since_onset: Optional[float] = None,
    ) -> ScaleResult:
        total = 0
        for item in cls.ITEMS:
            val = scores.get(item, 4)
            val = max(0, min(val, 4))
            total += val

        # Calculate decline rate if onset info available
        decline_rate = None
        if months_since_onset and months_since_onset > 0:
            decline_rate = (cls.MAX_SCORE - total) / months_since_onset

        if total >= 40:
            category = "Mild functional impairment"
            interpretation = "Mild functional impairment. Early disease stage."
        elif total >= 30:
            category = "Moderate functional impairment"
            interpretation = (
                "Moderate functional impairment. Consider multidisciplinary "
                "ALS clinic engagement."
            )
        elif total >= 20:
            category = "Significant functional impairment"
            interpretation = (
                "Significant functional decline. Evaluate for assistive devices, "
                "PEG placement, and respiratory support."
            )
        else:
            category = "Severe functional impairment"
            interpretation = (
                "Severe functional impairment. Hospice or palliative care "
                "discussion appropriate. Non-invasive ventilation evaluation."
            )

        recommendations = []
        # Respiratory sub-scores
        resp_items = ["dyspnea", "orthopnea", "respiratory_insufficiency"]
        resp_total = sum(max(0, min(scores.get(i, 4), 4)) for i in resp_items)
        if resp_total <= 8:
            recommendations.append(
                "Obtain pulmonary function tests (FVC); consider non-invasive "
                "ventilation (BiPAP) if FVC < 50%"
            )

        bulbar_items = ["speech", "salivation", "swallowing"]
        bulbar_total = sum(max(0, min(scores.get(i, 4), 4)) for i in bulbar_items)
        if bulbar_total <= 8:
            recommendations.append(
                "Refer to speech-language pathology; evaluate PEG tube placement"
            )

        if decline_rate is not None:
            if decline_rate > 1.0:
                recommendations.append(
                    f"Rapid decline rate ({decline_rate:.2f} points/month). "
                    "Consider aggressive supportive interventions and clinical trials."
                )
            interpretation += f" Decline rate: {decline_rate:.2f} points/month."

        if total <= 30:
            recommendations.append(
                "Multidisciplinary ALS clinic: neurology, PT/OT, respiratory, "
                "nutrition, social work, palliative care"
            )

        return ScaleResult(
            scale_type=ClinicalScaleType.ALSFRS,
            score=total,
            max_score=cls.MAX_SCORE,
            interpretation=interpretation,
            severity_category=category,
            thresholds={
                "mild_impairment": 40.0,
                "moderate_impairment": 30.0,
                "significant_impairment": 20.0,
                "rapid_decline_rate": 1.0,
            },
            recommendations=recommendations,
        )


# ===================================================================
# 9. ASPECTS Calculator  (Alberta Stroke Program Early CT Score)
# ===================================================================


class ASPECTSCalculator:
    """ASPECTS — 10 MCA territory CT regions, score 0-10.

    Score starts at 10; subtract 1 for each region with early ischemic
    change on non-contrast CT.

    Regions
    -------
    C  : Caudate
    L  : Lentiform nucleus
    IC : Internal capsule
    I  : Insular ribbon
    M1 : Anterior MCA cortex
    M2 : MCA cortex lateral to insular ribbon
    M3 : Posterior MCA cortex
    M4 : Anterior MCA territory above M1
    M5 : Lateral MCA territory above M2
    M6 : Posterior MCA territory above M3
    """

    REGIONS = ["C", "L", "IC", "I", "M1", "M2", "M3", "M4", "M5", "M6"]
    MAX_SCORE = 10

    @classmethod
    def calculate(cls, affected_regions: List[str]) -> ScaleResult:
        """Calculate ASPECTS from list of affected regions.

        Parameters
        ----------
        affected_regions : list of str
            Region codes showing early ischemic changes.
        """
        valid_affected = [r.upper() for r in affected_regions if r.upper() in cls.REGIONS]
        unique_affected = list(set(valid_affected))
        score = cls.MAX_SCORE - len(unique_affected)
        score = max(0, score)

        if score >= 8:
            category = "Favorable for intervention"
            interpretation = (
                f"ASPECTS {score}/10 — small ischemic core. Favorable for "
                "IV thrombolysis and mechanical thrombectomy."
            )
        elif score >= 6:
            category = "Intermediate"
            interpretation = (
                f"ASPECTS {score}/10 — moderate ischemic core. IV thrombolysis "
                "may still be considered. Thrombectomy decision should consider "
                "clinical-imaging mismatch."
            )
        else:
            category = "Unfavorable for intervention"
            interpretation = (
                f"ASPECTS {score}/10 — large established infarct core. "
                "Higher risk of hemorrhagic transformation with reperfusion. "
                "Thrombectomy benefit less certain (excluded from MR CLEAN, ESCAPE)."
            )

        recommendations = []
        if score >= 6:
            recommendations.append(
                "Eligible for IV thrombolysis if within time window and no contraindications"
            )
        if score >= 6:
            recommendations.append(
                "Mechanical thrombectomy indicated if LVO confirmed (ICA/M1)"
            )
        if score < 6:
            recommendations.append(
                "Discuss risks/benefits of intervention carefully — "
                "large core infarct with limited salvageable tissue"
            )
            recommendations.append(
                "Consider CTP (CT perfusion) to assess mismatch before thrombectomy"
            )
        if unique_affected:
            recommendations.append(
                f"Affected regions: {', '.join(sorted(unique_affected))}"
            )

        return ScaleResult(
            scale_type=ClinicalScaleType.ASPECTS,
            score=score,
            max_score=cls.MAX_SCORE,
            interpretation=interpretation,
            severity_category=category,
            thresholds={
                "thrombectomy_favorable": 6.0,
                "small_core": 8.0,
                "large_core_cutoff": 5.0,
            },
            recommendations=recommendations,
        )


# ===================================================================
# 10. Hoehn & Yahr Calculator
# ===================================================================


class HoehnYahrCalculator:
    """Hoehn & Yahr Scale — Parkinson's disease staging, 1-5.

    Stage 1  : Unilateral involvement only
    Stage 1.5: Unilateral and axial involvement
    Stage 2  : Bilateral involvement without balance impairment
    Stage 2.5: Mild bilateral disease with recovery on pull test
    Stage 3  : Mild-to-moderate bilateral disease; postural instability;
               physically independent
    Stage 4  : Severe disability; still able to walk or stand unassisted
    Stage 5  : Wheelchair-bound or bedridden unless aided
    """

    MAX_SCORE = 5.0

    _STAGES = {
        1.0: (
            "Stage 1 — Unilateral",
            "Unilateral involvement only. Minimal or no functional disability.",
        ),
        1.5: (
            "Stage 1.5 — Unilateral + axial",
            "Unilateral and axial involvement. Minimal functional disability.",
        ),
        2.0: (
            "Stage 2 — Bilateral, no balance impairment",
            "Bilateral involvement without impairment of balance. "
            "Mild functional disability.",
        ),
        2.5: (
            "Stage 2.5 — Mild bilateral with recovery on pull test",
            "Mild bilateral disease with recovery on retropulsion (pull) test.",
        ),
        3.0: (
            "Stage 3 — Bilateral with postural instability",
            "Mild-to-moderate bilateral disease with postural instability; "
            "physically independent.",
        ),
        4.0: (
            "Stage 4 — Severe, can walk/stand",
            "Severe disability; still able to walk or stand unassisted "
            "but markedly incapacitated.",
        ),
        5.0: (
            "Stage 5 — Wheelchair/bedridden",
            "Wheelchair-bound or bedridden unless aided. "
            "Requires full-time caregiver assistance.",
        ),
    }

    @classmethod
    def calculate(cls, stage: float) -> ScaleResult:
        # Snap to valid stages
        valid_stages = sorted(cls._STAGES.keys())
        stage = max(1.0, min(stage, 5.0))
        # Find nearest valid stage
        closest = min(valid_stages, key=lambda s: abs(s - stage))
        stage = closest

        category, interpretation = cls._STAGES[stage]

        recommendations = []
        if stage <= 2.0:
            recommendations.append(
                "Consider dopaminergic therapy if symptoms are bothersome"
            )
            recommendations.append(
                "Begin regular exercise program (shown to slow progression)"
            )
        if stage == 3.0:
            recommendations.append(
                "Optimise dopaminergic therapy; add MAO-B or COMT inhibitor if needed"
            )
            recommendations.append(
                "Physical therapy for balance and fall prevention"
            )
        if stage >= 3.0:
            recommendations.append(
                "Evaluate for DBS candidacy if levodopa-responsive with "
                "motor fluctuations or dyskinesias"
            )
        if stage >= 4.0:
            recommendations.append(
                "Multidisciplinary care: PT, OT, speech therapy, nutrition"
            )
            recommendations.append(
                "Assess cognitive function (MoCA) and screen for PD dementia"
            )
        if stage >= 5.0:
            recommendations.append(
                "Caregiver support and palliative care consultation"
            )

        return ScaleResult(
            scale_type=ClinicalScaleType.HOEHN_YAHR,
            score=stage,
            max_score=cls.MAX_SCORE,
            interpretation=interpretation,
            severity_category=category,
            thresholds={
                "bilateral_onset": 2.0,
                "postural_instability": 3.0,
                "severe_disability": 4.0,
                "wheelchair_bound": 5.0,
            },
            recommendations=recommendations,
        )
