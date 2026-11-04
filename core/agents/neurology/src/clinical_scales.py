"""Validated clinical scale calculators for the Neurology Intelligence Agent.

Implements 10 evidence-based neurological scoring instruments with
automated interpretation, threshold identification, and clinical
recommendations.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional, Tuple

from src.models import ClinicalScaleType, ScaleResult

logger = logging.getLogger(__name__)


class ScaleDataError(ValueError):
    """The submitted items cannot produce a score for this instrument.

    Raised instead of substituting a value for an item nobody recorded. Callers surface it
    as an explicit insufficient-data answer naming the missing or unusable input; the HTTP
    layer maps it to 422.
    """


def _coerce_int(value) -> Optional[int]:
    """Return *value* as an int, or None if it is not an integer-valued scalar.

    `str(v).isdigit()` was the old test for "is this a score", and it is false for -1, for
    "2.0", and for True. Anything this returns None for is reported as unusable rather than
    silently skipped.
    """
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value) if float(value).is_integer() else None
    if isinstance(value, str):
        text = value.strip()
        try:
            as_float = float(text)
        except ValueError:
            return None
        return int(as_float) if as_float.is_integer() else None
    return None


def _normalise_key(key: str) -> str:
    """Lower-case, strip punctuation runs, collapse separators.

    The Streamlit forms, the route payloads and the calculators were each written with their
    own spelling of the same item (`dressing` vs `dressing_hygiene`, `3.15a` vs
    `3.15a_postural_tremor_r`). A mismatch used to be absorbed by a default, which is how an
    unrecorded ALSFRS-R item became 4 = normal function.
    """
    return re.sub(r"[^a-z0-9]+", "_", str(key).strip().lower()).strip("_")


def _resolve_items(
    raw: Dict[str, object],
    canonical: Dict[str, int],
    aliases: Dict[str, str],
) -> Tuple[Dict[str, int], List[str]]:
    """Map a submitted item dict onto *canonical* item names.

    Returns (resolved, problems). `resolved` holds only items that were actually submitted,
    clamped to each item's maximum; `problems` names every key that could not be matched and
    every value that was out of range. Nothing is defaulted here.
    """
    lookup = {_normalise_key(name): name for name in canonical}
    lookup.update({_normalise_key(a): target for a, target in aliases.items()})

    resolved: Dict[str, int] = {}
    problems: List[str] = []
    for key, value in (raw or {}).items():
        norm = _normalise_key(key)
        target = lookup.get(norm)
        if target is None:
            # Try a prefix match: "3.15a" -> "3.15a_postural_tremor_r".
            matches = [name for name in canonical if _normalise_key(name).startswith(norm + "_")]
            if len(matches) == 1:
                target = matches[0]
        if target is None:
            problems.append(f"unrecognized item '{key}' ignored")
            continue
        as_int = _coerce_int(value)
        if as_int is None:
            problems.append(f"item '{key}' value {value!r} is not an integer score -- ignored")
            continue
        max_val = canonical[target]
        if as_int < 0 or as_int > max_val:
            problems.append(
                f"item '{target}' value {as_int} outside the valid range 0-{max_val} -- clamped"
            )
            as_int = max(0, min(as_int, max_val))
        if target in resolved and resolved[target] != as_int:
            problems.append(f"item '{target}' submitted twice with different values")
        resolved[target] = as_int
    return resolved, problems


def _missing_note(missing: List[str], total_items: int) -> str:
    shown = ", ".join(missing[:6])
    if len(missing) > 6:
        shown += f" (+{len(missing) - 6} more)"
    return (
        f"INCOMPLETE: {len(missing)} of {total_items} items were not recorded ({shown}). "
        "The total below is computed from the recorded items only and is not the patient's score."
    )


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

    # Spellings the UI and the route payloads use for the same items.
    ALIASES: Dict[str, str] = {
        "4_facial_palsy": "4_facial",
        "5a_motor_left_arm": "5a_left_arm",
        "5b_motor_right_arm": "5b_right_arm",
        "6a_motor_left_leg": "6a_left_leg",
        "6b_motor_right_leg": "6b_right_leg",
        "7_limb_ataxia": "7_ataxia",
        "11_extinction_inattention": "11_extinction",
    }

    @classmethod
    def calculate(cls, scores: Dict[str, int]) -> ScaleResult:
        """Calculate NIHSS total from item scores.

        Parameters
        ----------
        scores : dict
            Mapping of item names to integer scores. Items that were not examined are
            reported as not examined; they are NOT scored 0.

        Notes
        -----
        Two defects were fixed here.

        *Unexamined items were scored 0 = no deficit* (N20). `scores.get(item, 0)` meant a
        two-item exam produced a confident total: `{"1a_loc": 2, "9_language": 3}` returned
        5/42 "Moderate stroke" with nothing saying 13 items had never been examined. NIHSS >= 4
        gates alteplase and >= 6 gates thrombectomy in this agent, so a partial exam
        systematically under-triages. The total is now accompanied by `items_scored`,
        `complete` and a warning, and the callers refuse to run the treatment thresholds off an
        incomplete exam.

        *The interpretation carried treatment advice* (N4). Every band above "minor" read
        "IV tPA indicated" / "Consider mechanical thrombectomy", regardless of stroke type --
        and `AcuteStrokeTriageWorkflow` runs NIHSS for hemorrhagic stroke too, with
        `export_stroke_report` printing each scale's interpretation verbatim. An exported ICH
        report therefore endorsed alteplase three lines above "Intracerebral haemorrhage
        identified". Alteplase in ICH is lethal. NIHSS measures severity and knows nothing
        about stroke type, so it now describes severity only; all thrombolysis and thrombectomy
        language belongs to the stroke workflow, which does know the type.
        """
        resolved, problems = _resolve_items(scores or {}, cls.ITEMS, cls.ALIASES)
        total = sum(resolved.values())
        missing = [item for item in cls.ITEMS if item not in resolved]
        complete = not missing
        warnings = list(problems)
        if missing:
            warnings.append(_missing_note(missing, len(cls.ITEMS)))

        if total == 0:
            category = "No stroke symptoms"
            interpretation = "No measurable neurological deficit on the items recorded."
        elif total <= 4:
            category = "Minor stroke"
            interpretation = "Minor neurological deficit."
        elif total <= 15:
            category = "Moderate stroke"
            interpretation = "Moderate neurological deficit."
        elif total <= 20:
            category = "Moderate-to-severe stroke"
            interpretation = "Moderate-to-severe neurological deficit."
        elif total <= 25:
            category = "Severe stroke"
            interpretation = "Severe neurological deficit."
        else:
            category = "Very severe stroke"
            interpretation = "Very severe neurological deficit; ICU-level care required."

        if not complete:
            category = f"Indeterminate -- incomplete NIHSS ({len(resolved)}/{len(cls.ITEMS)} items)"
            interpretation = (
                f"NIHSS cannot be scored: {len(resolved)} of {len(cls.ITEMS)} items were "
                f"recorded. The sum of the recorded items is {total}, which is a lower bound "
                f"on the true total, not the patient's NIHSS. Complete the examination before "
                f"using any NIHSS threshold."
            )

        recommendations: List[str] = []
        # No treatment recommendation is issued from this calculator. NIHSS is a severity
        # instrument; thrombolysis and thrombectomy eligibility depend on stroke type, onset
        # time, imaging and contraindications, none of which are inputs here.
        if not complete:
            recommendations.append(
                "Complete the remaining NIHSS items before applying the alteplase (>= 4) or "
                "thrombectomy (>= 6) thresholds"
            )
        elif total > 25:
            recommendations.append(
                "ICU admission with close neurological monitoring"
            )
        recommendations.append(
            "Treatment eligibility is determined by the stroke-type-specific workflow, not by "
            "this score alone"
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
            items_scored=len(resolved),
            items_expected=len(cls.ITEMS),
            complete=complete,
            warnings=warnings,
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

    ALIASES: Dict[str, str] = {
        "visuospatial_executive": "visuospatial",
        "executive": "visuospatial",
        "recall": "delayed_recall",
    }

    @classmethod
    def calculate(
        cls, domain_scores: Dict[str, int], education_years: Optional[int] = None
    ) -> ScaleResult:
        """Score MoCA from the seven domain sub-scores.

        Parameters
        ----------
        domain_scores : dict
            The seven MoCA domains. `education_years` must NOT be passed inside this dict --
            it is a demographic, not a domain, and it is rejected here if it appears.
        education_years : int | None
            Years of formal education. The published MoCA adjustment adds 1 point when
            education is 12 years or fewer, applied once, and only when the raw total is
            below 30. With education unknown no adjustment is applied and the result says so;
            an unknown demographic is not grounds for giving the patient a free point.

        Notes
        -----
        N1: the live route twin summed *every* numeric value in the submitted items, and the
        UI puts `education_years` (default 16) in that same dict. A true 12/30 exam therefore
        reported 28/30 "Normal cognition" -- and because the whole workup block is gated on
        total < 26, the response carried no recommendations at all: no neuropsychological
        testing, no MRI, no B12/TSH/RPR screen for reversible causes. Moderate dementia read
        as a normal cognitive screen. Only the seven known domains are summed now, and the
        demographic is a separate parameter that cannot be added to a score.
        """
        raw = dict(domain_scores or {})
        stray_education = None
        for key in list(raw):
            if _normalise_key(key) in ("education_years", "education", "years_of_education"):
                stray_education = _coerce_int(raw.pop(key))

        resolved, problems = _resolve_items(raw, cls.DOMAINS, cls.ALIASES)
        total = sum(resolved.values())
        missing = [d for d in cls.DOMAINS if d not in resolved]
        complete = not missing
        warnings = list(problems)
        if missing:
            warnings.append(_missing_note(missing, len(cls.DOMAINS)))

        if education_years is None and stray_education is not None:
            # Tolerated so an existing client is not broken, but never summed into the total.
            education_years = stray_education
            warnings.append(
                "education_years was submitted among the MoCA domain scores; it is a "
                "demographic and was NOT added to the total"
            )

        education_adjustment = 0
        if education_years is None:
            warnings.append(
                "Years of education not provided -- the published +1 adjustment for <= 12 "
                "years of education was not applied"
            )
        elif education_years <= 12 and total < cls.MAX_SCORE:
            # Published MoCA scoring instruction: add one point if the examinee has 12 years
            # or fewer of formal education, provided the total is less than 30.
            education_adjustment = 1
            total += 1

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

        # Flag impaired domains. Only domains that were actually scored are considered:
        # an unrecorded domain used to read as 0 and therefore as "notably impaired".
        # The < 50%-of-maximum rule is this repo's own convention, not a published MoCA
        # rule, and is labelled as such.
        impaired = [
            domain for domain, max_val in cls.DOMAINS.items()
            if domain in resolved and resolved[domain] < max_val * 0.5
        ]
        if impaired:
            recommendations.append(
                "Domains scoring below half of their maximum (local flag, not a published "
                f"MoCA rule): {', '.join(impaired)}"
            )

        if not complete:
            # An incomplete MoCA must never be reported as normal: the unrecorded domains can
            # only move the total down.
            category = f"Indeterminate -- incomplete MoCA ({len(resolved)}/{len(cls.DOMAINS)} domains)"
            interpretation = (
                f"MoCA cannot be scored: {len(resolved)} of {len(cls.DOMAINS)} domains were "
                f"recorded. The recorded domains sum to {total}/{cls.MAX_SCORE}, an upper "
                f"bound on the true total. Complete the instrument before interpreting it."
            )
            recommendations.insert(
                0, "Administer the remaining MoCA domains before interpreting the total"
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
                "education_adjustment_applied": float(education_adjustment),
            },
            recommendations=recommendations,
            items_scored=len(resolved),
            items_expected=len(cls.DOMAINS),
            complete=complete,
            warnings=warnings,
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

    # The Streamlit form posts bare item names and bare numbers ("speech", "3.10"); the
    # canonical keys carry both ("3.1_speech", "3.10_gait"). `_resolve_items` prefix-matches
    # "3.10" on its own; these are the name-only spellings.
    ALIASES: Dict[str, str] = {
        "speech": "3.1_speech",
        "facial_expression": "3.2_facial_expression",
        "rigidity_neck": "3.3a_rigidity_neck",
        "rigidity_rue": "3.3b_rigidity_rue",
        "rigidity_lue": "3.3c_rigidity_lue",
        "rigidity_rle": "3.3d_rigidity_rle",
        "rigidity_lle": "3.3e_rigidity_lle",
        "finger_tapping_r": "3.4a_finger_tapping_r",
        "finger_tapping_l": "3.4b_finger_tapping_l",
        "hand_movements_r": "3.5a_hand_movements_r",
        "hand_movements_l": "3.5b_hand_movements_l",
        "pronation_supination_r": "3.6a_pronation_supination_r",
        "pronation_supination_l": "3.6b_pronation_supination_l",
        "toe_tapping_r": "3.7a_toe_tapping_r",
        "toe_tapping_l": "3.7b_toe_tapping_l",
        "leg_agility_r": "3.8a_leg_agility_r",
        "leg_agility_l": "3.8b_leg_agility_l",
        "arising_from_chair": "3.9_arising_from_chair",
        "gait": "3.10_gait",
        "freezing_of_gait": "3.11_freezing_of_gait",
        "postural_stability": "3.12_postural_stability",
        "posture": "3.13_posture",
        "body_bradykinesia": "3.14_body_bradykinesia",
        "postural_tremor_r": "3.15a_postural_tremor_r",
        "postural_tremor_l": "3.15b_postural_tremor_l",
        "kinetic_tremor_r": "3.16a_kinetic_tremor_r",
        "kinetic_tremor_l": "3.16b_kinetic_tremor_l",
        "rest_tremor_rue": "3.17a_rest_tremor_rue",
        "rest_tremor_lue": "3.17b_rest_tremor_lue",
        "rest_tremor_rle": "3.17c_rest_tremor_rle",
        "rest_tremor_lle": "3.17d_rest_tremor_lle",
        "rest_tremor_jaw": "3.17e_rest_tremor_jaw",
        "constancy_of_rest_tremor": "3.18_constancy_of_rest_tremor",
    }

    # Item groups for tremor-dominant / PIGD subtyping, named explicitly rather than matched
    # by substring. N13: `"postural" in key` also matched 3.15a/3.15b postural *tremor*, so a
    # pure-tremor patient contributed tremor points to the PIGD sum.
    TREMOR_ITEMS: List[str] = [
        "3.15a_postural_tremor_r", "3.15b_postural_tremor_l",
        "3.16a_kinetic_tremor_r", "3.16b_kinetic_tremor_l",
        "3.17a_rest_tremor_rue", "3.17b_rest_tremor_lue",
        "3.17c_rest_tremor_rle", "3.17d_rest_tremor_lle",
        "3.17e_rest_tremor_jaw", "3.18_constancy_of_rest_tremor",
    ]
    PIGD_ITEMS: List[str] = [
        "3.9_arising_from_chair", "3.10_gait", "3.11_freezing_of_gait",
        "3.12_postural_stability", "3.13_posture",
    ]

    @classmethod
    def calculate(cls, scores: Dict[str, int]) -> ScaleResult:
        resolved, problems = _resolve_items(scores or {}, cls.ITEMS, cls.ALIASES)
        total = sum(resolved.values())
        missing = [item for item in cls.ITEMS if item not in resolved]
        complete = not missing
        warnings = list(problems)
        if missing:
            warnings.append(_missing_note(missing, len(cls.ITEMS)))

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

        if not complete:
            category = (
                f"Indeterminate -- incomplete MDS-UPDRS III "
                f"({len(resolved)}/{len(cls.ITEMS)} items)"
            )
            interpretation = (
                f"MDS-UPDRS Part III cannot be scored: {len(resolved)} of {len(cls.ITEMS)} "
                f"items were recorded, summing to {total}. That sum is a lower bound, not the "
                f"patient's motor score; the severity bands and the DBS threshold do not apply "
                f"to it."
            )
            recommendations = [
                "Complete the remaining MDS-UPDRS Part III items before interpreting the total"
            ]

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
            items_scored=len(resolved),
            items_expected=len(cls.ITEMS),
            complete=complete,
            warnings=warnings,
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
        cls, fs_scores: Dict[str, int], edss_step: Optional[float]
    ) -> ScaleResult:
        """Interpret a clinician-assigned EDSS step, cross-checked against the FS scores.

        The EDSS step (0.0-10.0 in 0.5 increments) is assigned by the clinician from the
        functional system scores plus ambulation. This method does not derive it.

        Parameters
        ----------
        fs_scores : dict
            Functional system grades. Used to cross-check the step, and reported.
        edss_step : float | None
            The clinician's EDSS step. None means it was not recorded, and this method then
            returns an explicit insufficient-data result rather than a number.

        Notes
        -----
        N8, two parts.

        *The method ignored `fs_scores` entirely* while its docstring claimed it "validates"
        the score, so `calculate({"pyramidal": 6, "cerebellar": 5}, 0.0)` returned 0.0 "No or
        minimal disability. Fully ambulatory." The FS scores are now reported and an
        impossible pairing is flagged. They are still not used to *derive* the step: no
        published FS-to-step mapping is in hand here, and the route twin that did invent one
        (`edss = min(float(max_fs), 9.5)`) assigned EDSS 6.0 -- "requires a walking aid" -- to
        a blind but fully ambulatory patient on a visual FS of 6.

        *The top band was off by one step.* `elif edss <= 9.0 ... else "Death due to MS"` made
        9.5 read as "Death attributable to MS." EDSS 9.5 is a living, totally helpless patient;
        10.0 is death. This repo's own reference table agrees -- `src/knowledge.py`
        `CLINICAL_SCALES["edss"]` maps "7.0-9.5" to "Wheelchair / bed-bound" and "10.0" to
        "Death due to MS" -- so the band boundary contradicted the agent's own documented
        scale. Corrected to `<= 9.5` for bed-bound and 10.0 for death.
        """
        fs_resolved, warnings = _resolve_items(fs_scores or {}, cls.FS_MAXIMA, {})

        if edss_step is None:
            warnings.insert(
                0,
                "EDSS step not recorded. EDSS is assigned by the examiner from the functional "
                "system grades plus ambulation; it cannot be derived from the FS grades alone, "
                "so no score is reported.",
            )
            return ScaleResult(
                scale_type=ClinicalScaleType.EDSS,
                score=0.0,
                max_score=cls.MAX_SCORE,
                interpretation=(
                    "INSUFFICIENT DATA: no EDSS step was recorded. Functional system grades "
                    "received: "
                    + (", ".join(f"{k}={v}" for k, v in sorted(fs_resolved.items())) or "none")
                    + ". Record the examiner's EDSS step (0.0-10.0 in 0.5 increments)."
                ),
                severity_category="Indeterminate -- EDSS step not recorded",
                thresholds={},
                recommendations=[
                    "Record the examiner-assigned EDSS step before using any EDSS threshold",
                ],
                items_scored=len(fs_resolved),
                items_expected=len(cls.FS_MAXIMA),
                complete=False,
                warnings=warnings,
            )

        edss = max(0.0, min(edss_step, 10.0))
        # Round to nearest 0.5
        edss = round(edss * 2) / 2

        if fs_resolved:
            worst_fs = max(fs_resolved.values())
            if edss == 0.0 and worst_fs >= 1:
                # EDSS 0 is by definition a normal neurological examination (all FS grade 0).
                warnings.append(
                    f"EDSS step 0.0 was submitted with functional system grades up to "
                    f"{worst_fs}. EDSS 0 means a normal neurological examination, so the step "
                    f"and the FS grades are inconsistent -- verify before acting on either."
                )

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
        elif edss <= 9.5:
            # 7.0-9.5 is bed/chair-bound in this repo's own CLINICAL_SCALES table; 10.0 alone
            # is death. The previous boundary (<= 9.0) labelled a living EDSS 9.5 patient dead.
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

        if fs_resolved:
            recommendations.append(
                "Functional system grades recorded: "
                + ", ".join(f"{k}={v}" for k, v in sorted(fs_resolved.items()))
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
            items_scored=len(fs_resolved),
            items_expected=len(cls.FS_MAXIMA),
            complete=True,
            warnings=warnings,
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
    ITEM_COUNT = 6
    MIN_SCORE = 36
    MAX_SCORE = 78

    @classmethod
    def calculate(cls, responses: List[str]) -> ScaleResult:
        """Calculate HIT-6 from exactly 6 categorical responses.

        Parameters
        ----------
        responses : list
            Six responses, each one of 'never', 'rarely', 'sometimes', 'very_often',
            'always' -- or the corresponding published weight (6, 8, 10, 11, 13), which is
            what the Streamlit form posts.

        Raises
        ------
        ScaleDataError
            If fewer or more than six responses are supplied, or if any response is not a
            recognised category or weight.

        Notes
        -----
        N9: short lists were padded with "never" and `ITEM_VALUES.get(r.lower(), 6)` mapped
        anything unrecognised to "never" as well -- 6 points, the best possible answer on
        every item. `calculate(["always", "always"])` returned 50, "Some impact", which is
        below the >= 56 preventive-therapy trigger this agent uses: a patient who answered
        both available items at maximum severity was denied preventive therapy. A typo such as
        "very often" (space instead of underscore) scored 6 instead of 11. HIT-6 has no
        partial-credit scoring and its floor is 36 from six answered items, so an incomplete
        instrument has no score; it is refused rather than imputed.
        """
        items = list(responses or [])
        if len(items) != cls.ITEM_COUNT:
            raise ScaleDataError(
                f"HIT-6 requires exactly {cls.ITEM_COUNT} responses; {len(items)} were "
                f"supplied. There is no valid partial HIT-6 score -- an unanswered item "
                f"cannot be imputed."
            )

        weights: List[int] = []
        valid_weights = set(cls.ITEM_VALUES.values())
        for index, response in enumerate(items, 1):
            as_int = _coerce_int(response)
            if as_int is not None:
                if as_int not in valid_weights:
                    raise ScaleDataError(
                        f"HIT-6 item {index}: {as_int} is not a HIT-6 item weight. Valid "
                        f"weights are {sorted(valid_weights)} "
                        f"(never/rarely/sometimes/very_often/always)."
                    )
                weights.append(as_int)
                continue
            key = _normalise_key(response)
            if key not in cls.ITEM_VALUES:
                raise ScaleDataError(
                    f"HIT-6 item {index}: {response!r} is not a recognised response. Valid "
                    f"responses are {sorted(cls.ITEM_VALUES)}."
                )
            weights.append(cls.ITEM_VALUES[key])

        total = sum(weights)

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
            items_scored=cls.ITEM_COUNT,
            items_expected=cls.ITEM_COUNT,
            complete=True,
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

    ALIASES: Dict[str, str] = {
        # The Streamlit form's key for item 6.
        "dressing": "dressing_hygiene",
        "hygiene": "dressing_hygiene",
        "cutting_food_handling_utensils": "cutting_food",
        "turning_in_bed_adjusting_bed_clothes": "turning_in_bed",
        "dyspnoea": "dyspnea",
        "orthopnoea": "orthopnea",
        "respiratory": "respiratory_insufficiency",
    }

    MAX_SCORE = 48

    @classmethod
    def calculate(
        cls,
        scores: Dict[str, int],
        months_since_onset: Optional[float] = None,
    ) -> ScaleResult:
        """Score ALSFRS-R from the twelve item scores (each 0-4, 4 = normal function).

        Notes
        -----
        N5: `scores.get(item, 4)` scored every unrecorded item as 4 = normal function. An ALS
        patient assessed on two items -- `{"walking": 1, "climbing_stairs": 0}` -- returned
        41/48, "Mild functional impairment. Early disease stage.", with zero recommendations:
        the respiratory sub-score computed 12/12 out of defaults, so the PFT/BiPAP
        recommendation never fired, and the bulbar sub-score likewise. The decline rate,
        computed as (48 - total) / months, was understated by the same mechanism. The default
        also absorbed key-name mismatches: the UI posts `dressing` where the calculator
        expected `dressing_hygiene`, so that item silently scored 4.

        Only recorded items are summed now. An incomplete instrument is reported as
        indeterminate, the sub-score recommendations fire only when their items were actually
        recorded, and the decline rate is computed only from a complete instrument.
        """
        item_maxima = {item: 4 for item in cls.ITEMS}
        resolved, problems = _resolve_items(scores or {}, item_maxima, cls.ALIASES)
        total = sum(resolved.values())
        missing = [item for item in cls.ITEMS if item not in resolved]
        complete = not missing
        warnings = list(problems)
        if missing:
            warnings.append(_missing_note(missing, len(cls.ITEMS)))

        # Decline rate requires a complete instrument: (48 - total) / months read every
        # unrecorded item as full function, which understates decline.
        decline_rate = None
        if complete and months_since_onset and months_since_onset > 0:
            decline_rate = (cls.MAX_SCORE - total) / months_since_onset
        elif not complete and months_since_onset:
            warnings.append(
                "Decline rate not computed: it requires all 12 ALSFRS-R items"
            )

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

        # Sub-scores are evaluated only when every item in the sub-score was recorded.
        # A sub-score built out of defaults is what kept the respiratory recommendation from
        # ever firing for a patient whose respiratory items were never asked about.
        def _subscore(label: str, items: List[str], advice: str) -> None:
            present = [i for i in items if i in resolved]
            if len(present) == len(items):
                if sum(resolved[i] for i in items) <= 8:
                    recommendations.append(advice)
            else:
                warnings.append(
                    f"{label} sub-score not evaluated: "
                    f"{', '.join(i for i in items if i not in resolved)} not recorded"
                )

        _subscore(
            "Respiratory",
            ["dyspnea", "orthopnea", "respiratory_insufficiency"],
            "Obtain pulmonary function tests (FVC); consider non-invasive "
            "ventilation (BiPAP) if FVC < 50%",
        )
        _subscore(
            "Bulbar",
            ["speech", "salivation", "swallowing"],
            "Refer to speech-language pathology; evaluate PEG tube placement",
        )

        if decline_rate is not None:
            if decline_rate > 1.0:
                recommendations.append(
                    f"Rapid decline rate ({decline_rate:.2f} points/month). "
                    "Consider aggressive supportive interventions and clinical trials."
                )
            interpretation += f" Decline rate: {decline_rate:.2f} points/month."

        if complete and total <= 30:
            recommendations.append(
                "Multidisciplinary ALS clinic: neurology, PT/OT, respiratory, "
                "nutrition, social work, palliative care"
            )

        if not complete:
            category = (
                f"Indeterminate -- incomplete ALSFRS-R "
                f"({len(resolved)}/{len(cls.ITEMS)} items)"
            )
            interpretation = (
                f"ALSFRS-R cannot be scored: {len(resolved)} of {len(cls.ITEMS)} items were "
                f"recorded, summing to {total}. That sum is a lower bound on the patient's "
                f"total, not their score, and the functional-impairment bands do not apply to "
                f"it. Unrecorded items are not normal function."
            )
            recommendations.insert(
                0, "Record the remaining ALSFRS-R items before interpreting the total"
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
            items_scored=len(resolved),
            items_expected=len(cls.ITEMS),
            complete=complete,
            warnings=warnings,
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

    # Long names a radiologist may write instead of the short codes.
    REGION_ALIASES: Dict[str, str] = {
        "caudate": "C",
        "caudate_nucleus": "C",
        "caudate_head": "C",
        "l": "L",
        "lentiform": "L",
        "lentiform_nucleus": "L",
        "putamen": "L",
        "globus_pallidus": "L",
        "ic": "IC",
        "internal_capsule": "IC",
        "i": "I",
        "insula": "I",
        "insular_ribbon": "I",
        "insular_cortex": "I",
        "m1_cortex": "M1",
        "m2_cortex": "M2",
        "m3_cortex": "M3",
        "m4_cortex": "M4",
        "m5_cortex": "M5",
        "m6_cortex": "M6",
        "anterior_mca_cortex": "M1",
        "posterior_mca_cortex": "M3",
    }

    @classmethod
    def calculate(cls, affected_regions: List[str]) -> ScaleResult:
        """Calculate ASPECTS from the list of regions showing early ischemic change.

        Parameters
        ----------
        affected_regions : list of str
            Region codes (C, L, IC, I, M1-M6) or their long names.

        Raises
        ------
        ScaleDataError
            If any token cannot be resolved to one of the ten ASPECTS regions. The caller
            cannot score a CT whose regions it does not understand.

        Notes
        -----
        N17: `[r.upper() for r in affected_regions if r.upper() in cls.REGIONS]` discarded
        every token it did not recognise, in silence. `calculate(["Caudate", "Insula",
        "M1 cortex"])` returned 10/10, "Favorable for intervention", "small ischemic core ...
        Favorable for IV thrombolysis and mechanical thrombectomy" -- a radiologist writing
        words instead of codes turned a 7/10 CT into a perfect one, and this score feeds the
        thrombectomy gate in the stroke workflow. Long names are now accepted as aliases, and
        a token that still cannot be resolved stops the calculation instead of improving the
        score.
        """
        unresolved: List[str] = []
        resolved_regions: List[str] = []
        for raw in affected_regions or []:
            key = _normalise_key(raw)
            code = raw.strip().upper() if raw.strip().upper() in cls.REGIONS else None
            if code is None:
                code = cls.REGION_ALIASES.get(key)
            if code is None and key.upper() in cls.REGIONS:
                code = key.upper()
            if code is None:
                unresolved.append(str(raw))
            else:
                resolved_regions.append(code)

        if unresolved:
            raise ScaleDataError(
                "ASPECTS region(s) not recognised: "
                + ", ".join(repr(u) for u in unresolved)
                + f". Valid regions are {cls.REGIONS} (long names such as 'caudate', "
                "'lentiform', 'internal capsule', 'insular ribbon' are also accepted). "
                "Unrecognised regions are not scored as normal."
            )

        unique_affected = sorted(set(resolved_regions), key=cls.REGIONS.index)
        score = cls.MAX_SCORE - len(unique_affected)
        score = max(0, min(score, cls.MAX_SCORE))

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
                f"Affected regions: {', '.join(unique_affected)}"
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
            items_scored=len(cls.REGIONS),
            items_expected=len(cls.REGIONS),
            complete=True,
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


# ===================================================================
# Single dispatch entry point
# ===================================================================
#
# Neurology shipped the same clinical logic three times: these calculators, a divergent
# inline set in `api/routes/neuro_clinical.py` (`_calculate_nihss`, `_calculate_moca`, ...),
# and a third set of severity bands in `src/knowledge.py:CLINICAL_SCALES`. The shipped
# Streamlit UI called the route-local set, so the duplicates were what a user actually saw,
# and they disagreed with these on bandings, on defaults and on tPA eligibility -- which is
# why a tPA fix that landed in `clinical_workflows.py` did not change what the UI showed.
#
# `calculate_scale` is now the one implementation. The route-local `_calculate_*` functions
# are thin adapters over it, so there is one answer per scale per input.


def _single_value(items: Dict[str, object], names: List[str]) -> Optional[float]:
    """First usable numeric value among *names* in *items*, else None."""
    normalised = {_normalise_key(k): v for k, v in (items or {}).items()}
    for name in names:
        if name in normalised:
            value = normalised[name]
            if isinstance(value, bool):
                continue
            try:
                return float(value)
            except (TypeError, ValueError):
                continue
    return None


def _scale_nihss(items: Dict[str, object]) -> ScaleResult:
    return NIHSSCalculator.calculate(items)


def _scale_gcs(items: Dict[str, object]) -> ScaleResult:
    """GCS from eye/verbal/motor.

    A missing component used to default to 1 -- the worst possible value for that component,
    invented and then printed as a total. All three are required.
    """
    components = {}
    for name, aliases in (
        ("eye", ["eye", "e", "eye_opening"]),
        ("verbal", ["verbal", "v", "verbal_response"]),
        ("motor", ["motor", "m", "motor_response"]),
    ):
        value = _single_value(items, aliases)
        if value is None:
            raise ScaleDataError(
                f"GCS requires the {name} component; it was not provided. GCS is the sum of "
                f"three observations and cannot be scored from a subset."
            )
        components[name] = int(value)
    return GCSCalculator.calculate(**components)


def _scale_moca(items: Dict[str, object]) -> ScaleResult:
    domains = {
        k: v for k, v in (items or {}).items()
        if _normalise_key(k) not in ("education_years", "education", "years_of_education")
    }
    education = _single_value(items, ["education_years", "education", "years_of_education"])
    return MoCACalculator.calculate(
        domains, education_years=int(education) if education is not None else None
    )


def _scale_updrs(items: Dict[str, object]) -> ScaleResult:
    return UPDRSCalculator.calculate(items)


def _scale_edss(items: Dict[str, object]) -> ScaleResult:
    """EDSS from the examiner's step plus functional system grades.

    The route twin derived a step from `max(fs_scores)` -- `edss = min(float(max_fs), 9.5)` --
    which assigned EDSS 6.0, "requires a walking aid", to a blind but fully ambulatory patient
    whose visual FS was 6, and labelled its own 6.0-6.5 band "restricted to wheelchair". No
    published FS-to-step mapping is in hand, so nothing is derived: the examiner's step is
    required, and without it the result is explicitly indeterminate.
    """
    fs = {
        k: v for k, v in (items or {}).items()
        if _normalise_key(k) not in ("ambulation", "edss_step", "score", "edss")
    }
    step = _single_value(items, ["edss_step", "edss", "score"])
    if step is None:
        # The Streamlit form posts `ambulation` with 0.0 documented as "use FS only", so 0.0
        # there is a sentinel for "not recorded" rather than a genuine EDSS of 0.
        ambulation = _single_value(items, ["ambulation"])
        if ambulation is not None and ambulation > 0:
            step = ambulation
    return EDSSCalculator.calculate(fs, step)


def _scale_mrs(items: Dict[str, object]) -> ScaleResult:
    score = _single_value(items, ["score", "mrs", "mrs_score", "rankin"])
    if score is None:
        raise ScaleDataError(
            "mRS requires the clinician's 0-6 grade; none was provided. An absent grade is "
            "not 0 (no symptoms)."
        )
    return mRSCalculator.calculate(int(score))


def _scale_hit6(items: Dict[str, object]) -> ScaleResult:
    """HIT-6 from the six item responses.

    Accepts either the six published weights (what the Streamlit form posts) or the six
    category names. The route twin summed whatever numbers it was sent with no check that they
    were HIT-6 weights at all, so six 1-5 Likert values totalled ~6-30 and came back as
    "No significant impact".
    """
    if isinstance(items, (list, tuple)):
        responses = list(items)
    else:
        responses = list((items or {}).values())
    return HIT6Calculator.calculate(responses)


def _scale_alsfrs(items: Dict[str, object]) -> ScaleResult:
    months = _single_value(items, ["months_since_onset", "months_from_onset"])
    scores = {
        k: v for k, v in (items or {}).items()
        if _normalise_key(k) not in ("months_since_onset", "months_from_onset")
    }
    return ALSFRSCalculator.calculate(scores, months_since_onset=months)


def _scale_aspects(items: Dict[str, object]) -> ScaleResult:
    """ASPECTS from either a list of affected regions or a per-region dict.

    The dict form accepts 1/0 (1 = normal, 0 = early ischemic change, as the UI posts) and
    booleans (True = affected). The route twin counted only values equal to the *string* "0",
    so `{"M1": true}` -- a region explicitly marked affected -- yielded 10/10.
    """
    if isinstance(items, (list, tuple, set)):
        return ASPECTSCalculator.calculate(list(items))

    affected: List[str] = []
    unusable: List[str] = []
    for key, value in (items or {}).items():
        if _normalise_key(key) in ("affected_regions", "regions"):
            if isinstance(value, (list, tuple, set)):
                affected.extend(str(v) for v in value)
                continue
            unusable.append(f"{key}={value!r}")
            continue
        if isinstance(value, bool):
            if value:
                affected.append(key)
            continue
        as_int = _coerce_int(value)
        if as_int is None:
            unusable.append(f"{key}={value!r}")
            continue
        if as_int == 0:
            affected.append(key)
        elif as_int != 1:
            unusable.append(f"{key}={value!r}")

    if unusable:
        raise ScaleDataError(
            "ASPECTS per-region values must be 1 (normal), 0 (early ischemic change) or a "
            "boolean. Unusable: " + ", ".join(unusable)
        )
    return ASPECTSCalculator.calculate(affected)


def _scale_hoehn_yahr(items: Dict[str, object]) -> ScaleResult:
    stage = _single_value(items, ["stage", "score", "hoehn_yahr", "hy_stage"])
    if stage is None:
        raise ScaleDataError(
            "Hoehn & Yahr requires the clinician's stage; none was provided."
        )
    return HoehnYahrCalculator.calculate(stage)


_CALCULATORS = {
    "nihss": _scale_nihss,
    "gcs": _scale_gcs,
    "moca": _scale_moca,
    "updrs": _scale_updrs,
    "updrs_part_iii": _scale_updrs,
    "updrs_iii": _scale_updrs,
    "edss": _scale_edss,
    "mrs": _scale_mrs,
    "hit6": _scale_hit6,
    "hit_6": _scale_hit6,
    "alsfrs": _scale_alsfrs,
    "alsfrs_r": _scale_alsfrs,
    "aspects": _scale_aspects,
    "hoehn_yahr": _scale_hoehn_yahr,
}

SUPPORTED_SCALES = sorted(_CALCULATORS)


def calculate_scale(scale_name: str, items) -> ScaleResult:
    """Score any supported neurological scale from a submitted item mapping.

    Parameters
    ----------
    scale_name : str
        One of SUPPORTED_SCALES (case- and hyphen-insensitive).
    items : dict | list
        The submitted items. Shapes the UI and the API actually post are accepted; nothing
        is defaulted on the patient's behalf.

    Raises
    ------
    KeyError
        Unknown scale name.
    ScaleDataError
        The submitted items cannot produce a score for this instrument.
    """
    key = _normalise_key(scale_name)
    calculator = _CALCULATORS.get(key)
    if calculator is None:
        raise KeyError(
            f"Unknown scale '{scale_name}'. Supported scales: {SUPPORTED_SCALES}"
        )
    return calculator(items)
