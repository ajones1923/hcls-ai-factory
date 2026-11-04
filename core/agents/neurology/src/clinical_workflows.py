"""Clinical workflows for the Neurology Intelligence Agent.

Author: Adam Jones
Date: March 2026

Implements eight evidence-based clinical workflows that integrate imaging,
electrophysiology, biomarker, and genomic data to produce actionable
neurology assessments.  Each workflow follows the BaseNeuroWorkflow
contract (preprocess -> execute -> postprocess) and is registered in the
WorkflowEngine for unified dispatch.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from typing import Dict, List, Optional

from src.clinical_scales import (
    ScaleDataError,
    ALSFRSCalculator,
    ASPECTSCalculator,
    EDSSCalculator,
    HIT6Calculator,
    HoehnYahrCalculator,
    MoCACalculator,
    NIHSSCalculator,
    UPDRSCalculator,
)
from src.models import (
    DementiaSubtype,
    HeadacheType,
    NeuroDomain,
    NeuroWorkflowType,
    NMJPattern,
    ParkinsonsSubtype,
    ScaleResult,
    SeverityLevel,
    TumorGrade,
    TumorMolecularMarker,
    WorkflowResult,
)

logger = logging.getLogger(__name__)


# ===================================================================
# HELPERS
# ===================================================================

_SEVERITY_ORDER: List[SeverityLevel] = [
    SeverityLevel.INFORMATIONAL,
    SeverityLevel.LOW,
    SeverityLevel.MODERATE,
    SeverityLevel.HIGH,
    SeverityLevel.CRITICAL,
]


def _max_severity(*levels: SeverityLevel) -> SeverityLevel:
    """Return the highest severity among the given levels."""
    return max(levels, key=lambda s: _SEVERITY_ORDER.index(s))


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _trigger_string(trigger_type: str, genes: List[str], reason: str) -> str:
    """Build a human-readable cross-modal trigger string."""
    gene_str = ", ".join(genes[:8])
    if len(genes) > 8:
        gene_str += f" (+{len(genes) - 8} more)"
    return f"[{trigger_type}] Genes: {gene_str} -- {reason}"


# ===================================================================
# BASE CLASS
# ===================================================================


class BaseNeuroWorkflow(ABC):
    """Abstract base for all neurology clinical workflows."""

    workflow_type: NeuroWorkflowType
    domain: Optional[NeuroDomain] = None

    # -- template-method orchestrator ----------------------------------
    def run(self, inputs: dict) -> WorkflowResult:
        """Orchestrate preprocess -> execute -> postprocess."""
        logger.info("Running workflow %s", self.workflow_type.value)
        processed_inputs = self.preprocess(inputs)
        result = self.execute(processed_inputs)
        result = self.postprocess(result)
        # Inject any validation warnings as findings
        warnings = processed_inputs.get("_validation_warnings", [])
        if warnings:
            result.findings = [
                f"[INPUT WARNING] {w}" for w in warnings
            ] + result.findings
        return result

    def preprocess(self, inputs: dict) -> dict:
        """Validate and normalise raw inputs.  Override for workflow-specific logic."""
        return dict(inputs)

    @abstractmethod
    def execute(self, inputs: dict) -> WorkflowResult:
        """Core clinical logic.  Must be implemented by each workflow."""
        ...

    def postprocess(self, result: WorkflowResult) -> WorkflowResult:
        """Shared enrichment after execution.  Override to add workflow-specific post-steps."""
        try:
            from api.routes.events import publish_event
            publish_event("workflow_complete", {
                "workflow": result.workflow_type.value if hasattr(result.workflow_type, 'value') else str(result.workflow_type),
                "severity": result.severity.value if hasattr(result.severity, 'value') else str(result.severity),
                "findings_count": len(result.findings),
            })
        except Exception:
            pass  # Don't break workflow for event publishing failure
        return result

    @staticmethod
    def _init_warnings(inp: dict) -> list:
        """Initialise and return the validation warnings list on *inp*."""
        warnings: list = inp.setdefault("_validation_warnings", [])
        return warnings


# ===================================================================
# WORKFLOW 1 -- Acute Stroke Triage
# ===================================================================


class AcuteStrokeTriageWorkflow(BaseNeuroWorkflow):
    """Acute stroke triage: NIHSS, ASPECTS, tPA eligibility, thrombectomy
    criteria (DAWN/DEFUSE-3), and stroke type classification.

    Inputs
    ------
    nihss_scores : dict
        Per-item NIHSS scores (see NIHSSCalculator).
    onset_hours : float
        Hours since symptom onset.  Use 0 for unknown / wake-up stroke.
    ct_affected_regions : list[str]
        ASPECTS regions with early ischemic change on NCCT.
    age : int
        Patient age in years.
    blood_pressure_systolic : int
        Current systolic BP in mmHg.
    blood_glucose : float
        Blood glucose in mg/dL.
    on_anticoagulation : bool
        Currently on anticoagulation therapy.
    inr : float | None
        International normalised ratio (if on warfarin).
    platelet_count : int | None
        Platelet count (x10^3/uL).
    stroke_type_hint : str | None
        One of 'ischemic', 'hemorrhagic', 'tia', 'sah' if already known.
    large_vessel_occlusion : bool | None
        Whether CTA shows LVO (ICA, M1, basilar).
    core_volume_ml : float | None
        Perfusion CT core volume (mL).
    mismatch_ratio : float | None
        Perfusion CT penumbra:core mismatch ratio.
    """

    workflow_type = NeuroWorkflowType.ACUTE_STROKE
    domain = NeuroDomain.CEREBROVASCULAR

    # tPA contraindications (absolute)
    _TPA_CONTRAINDICATIONS = {
        "active_internal_bleeding",
        "recent_intracranial_surgery_3mo",
        "recent_stroke_3mo",
        "intracranial_hemorrhage",
        "intracranial_neoplasm",
        "arteriovenous_malformation",
        "aortic_dissection",
        "bleeding_diathesis",
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        # Default onset to unknown
        if "onset_hours" not in inp:
            inp["onset_hours"] = 0.0
            warnings.append("Onset time not specified; treating as unknown/wake-up stroke")
        if "nihss_scores" not in inp:
            inp["nihss_scores"] = {}
            warnings.append("No NIHSS item scores provided; NIHSS cannot be computed")
        elif inp["nihss_scores"]:
            # The warning above fired only when the dict was entirely absent, so a PARTIAL
            # NIHSS passed through silently and the missing items scored 0 = no deficit (N20).
            scored = [k for k in NIHSSCalculator.ITEMS if k in inp["nihss_scores"]]
            aliased = [
                k for k in inp["nihss_scores"]
                if k in NIHSSCalculator.ALIASES
            ]
            if len(scored) + len(aliased) < len(NIHSSCalculator.ITEMS):
                warnings.append(
                    f"NIHSS examination incomplete: "
                    f"{len(scored) + len(aliased)}/{len(NIHSSCalculator.ITEMS)} items "
                    f"recorded. Unexamined items are NOT scored 0; the treatment thresholds "
                    f"(>= 4 alteplase, >= 6 thrombectomy) are not applied to a partial exam."
                )
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        guidelines: List[str] = []
        scale_results: List[ScaleResult] = []
        severity = SeverityLevel.HIGH

        # -- NIHSS --
        nihss_result = NIHSSCalculator.calculate(inputs.get("nihss_scores", {}))
        scale_results.append(nihss_result)
        nihss = nihss_result.score
        findings.append(f"NIHSS total: {int(nihss)} ({nihss_result.severity_category})")

        nihss_complete = nihss_result.complete

        # -- ASPECTS --
        ct_regions = inputs.get("ct_affected_regions", [])
        aspects_result = None
        aspects_scored = False
        if ct_regions:
            try:
                aspects_result = ASPECTSCalculator.calculate(ct_regions)
            except ScaleDataError as exc:
                # A region name the calculator does not recognise used to be dropped in
                # silence, turning a 7/10 CT into 10/10 "Favorable for intervention" (N17).
                findings.append(f"ASPECTS NOT SCORED: {exc}")
            else:
                aspects_scored = True
                scale_results.append(aspects_result)
                findings.append(
                    f"ASPECTS: {int(aspects_result.score)}/10 "
                    f"({aspects_result.severity_category})"
                )
        elif inputs.get("ct_performed"):
            findings.append(
                "ASPECTS not reported although a CT was performed -- an unscored CT is not a "
                "favourable CT"
            )

        # -- Stroke type classification --
        stroke_hint = inputs.get("stroke_type_hint")
        if stroke_hint:
            stroke_type_str = stroke_hint
        else:
            # Default to ischemic if not specified
            stroke_type_str = "ischemic"
        findings.append(f"Stroke type: {stroke_type_str}")

        onset = inputs.get("onset_hours", 0.0) or 0.0
        # 0.0 is this workflow's standing marker for "onset time unknown" (see the branches
        # below), and a MISSING onset_hours defaulted to it. That is the right marker and the
        # wrong consequence: an unknown onset used to leave tPA eligibility untouched, so the
        # verdict came out "tPA ELIGIBLE -- administer IV alteplase urgently" for a patient
        # whose time of onset nobody knew. Thrombolysis cannot be cleared against a window that
        # has not been established.
        onset_known = onset > 0.0
        age = inputs.get("age")
        # N3: these defaulted to 0 / 0.0 / False, so an unmeasured blood pressure ran no check
        # and emitted no warning, and an unmeasured glucose printed the fabricated value
        # "glucose 0.0 mg/dL" as a finding. None now means "not measured", which blocks
        # thrombolysis pending the result rather than clearing it.
        sbp = inputs.get("blood_pressure_systolic")
        dbp = inputs.get("blood_pressure_diastolic")
        glucose = inputs.get("blood_glucose")
        on_anticoag = inputs.get("on_anticoagulation")
        inr = inputs.get("inr")
        platelets = inputs.get("platelet_count")
        contraindications = set(inputs.get("contraindications", []))
        lvo = inputs.get("large_vessel_occlusion")
        prestroke_mrs = inputs.get("prestroke_mrs", inputs.get("prior_mrs"))
        core_vol = inputs.get("core_volume_ml")
        mismatch = inputs.get("mismatch_ratio")
        mismatch_vol = inputs.get("mismatch_volume_ml")

        # -- tPA eligibility (ischemic only) --
        if stroke_type_str in ("ischemic", "tia"):
            tpa_eligible = True
            tpa_reasons: List[str] = []

            # Time window
            if onset > 4.5:
                tpa_eligible = False
                tpa_reasons.append(f"Onset {onset:.1f}h exceeds 4.5-hour window")
            elif not onset_known:
                # NOT eligible, rather than eligible-with-a-note. AHA/ASA 2019 requires a
                # established onset within the window; for wake-up or unwitnessed stroke the
                # route to thrombolysis is DWI-FLAIR mismatch (WAKE-UP, NEJM 2018), which is a
                # finding this workflow does not have.
                tpa_eligible = False
                findings.append(
                    "Unknown onset: consider MRI DWI-FLAIR mismatch for tPA eligibility"
                )
                tpa_reasons.append(
                    "Onset time not established -- tPA requires a known onset within 4.5h, or "
                    "DWI-FLAIR mismatch on MRI (WAKE-UP); neither is available here"
                )

            # Age exclusion for 3-4.5h window
            if 3.0 < onset <= 4.5:
                if age is None:
                    tpa_eligible = False
                    tpa_reasons.append(
                        "INSUFFICIENT DATA: age not recorded, and the extended (3-4.5h) "
                        "window excludes patients over 80"
                    )
                elif age > 80:
                    tpa_eligible = False
                    tpa_reasons.append("Age >80 in extended (3-4.5h) window")

            # NIHSS completeness. An incomplete exam understates the score, and the score
            # gates alteplase at >= 4 (N20).
            if not nihss_complete:
                tpa_eligible = False
                tpa_reasons.append(
                    f"INSUFFICIENT DATA: NIHSS incomplete "
                    f"({nihss_result.items_scored}/{nihss_result.items_expected} items) -- "
                    f"the alteplase threshold cannot be applied to a partial examination"
                )

            # BP thresholds. AHA/ASA 2019 Early Management of Acute Ischemic Stroke requires
            # BP below 185/110 mmHg before IV alteplase. Above that it is an exclusion, which
            # is what this block failed to do: it appended a reason, left `tpa_eligible` True,
            # and the report then carried both "Reduce BP to <185/110 before IV alteplase" and
            # "tPA ELIGIBLE -- administer IV alteplase urgently" (N3).
            if sbp is None:
                tpa_eligible = False
                tpa_reasons.append(
                    "INSUFFICIENT DATA: blood pressure not measured -- tPA requires "
                    "BP <185/110 mmHg"
                )
            elif sbp > 185:
                tpa_eligible = False
                tpa_reasons.append(
                    f"SBP {sbp} mmHg exceeds the 185 mmHg threshold -- NOT eligible until BP "
                    f"is lowered; re-triage after treatment"
                )
                recommendations.append(
                    "Reduce BP to <185/110 mmHg, then re-assess alteplase eligibility"
                )
            if dbp is not None and dbp > 110:
                tpa_eligible = False
                tpa_reasons.append(
                    f"DBP {dbp} mmHg exceeds the 110 mmHg threshold -- NOT eligible until BP "
                    f"is lowered"
                )

            # Glucose
            if glucose is None:
                tpa_eligible = False
                tpa_reasons.append(
                    "INSUFFICIENT DATA: blood glucose not measured -- hypoglycaemia and "
                    "hyperglycaemia both mimic stroke and must be excluded before tPA"
                )
            elif glucose < 50:
                tpa_eligible = False
                tpa_reasons.append(f"Hypoglycaemia (glucose {glucose} mg/dL)")
            elif glucose > 400:
                tpa_eligible = False
                tpa_reasons.append(
                    f"Hyperglycaemia (glucose {glucose} mg/dL) -- correct and reassess "
                    f"before tPA"
                )

            # Anticoagulation. "On anticoagulation -- verify INR/DOAC levels" was emitted
            # while eligibility stayed True, i.e. the patient was cleared for a thrombolytic
            # on the strength of a verification that had not happened (N3).
            if on_anticoag is None:
                tpa_eligible = False
                tpa_reasons.append(
                    "INSUFFICIENT DATA: anticoagulation status not recorded"
                )
            elif on_anticoag:
                if inr is None:
                    tpa_eligible = False
                    tpa_reasons.append(
                        "On anticoagulation with no INR/DOAC level available -- NOT eligible "
                        "pending coagulation results"
                    )
                elif inr > 1.7:
                    tpa_eligible = False
                    # NOTE (uncited): INR > 1.7 is consistent with AHA/ASA 2019 as commonly
                    # published; the primary table was not retrieved, so the value is
                    # unchanged.
                    tpa_reasons.append(f"INR {inr} > 1.7 while on anticoagulation")

            # Platelets
            if platelets is None:
                tpa_eligible = False
                tpa_reasons.append(
                    "INSUFFICIENT DATA: platelet count not available -- tPA requires "
                    "platelets >= 100K"
                )
            elif platelets < 100:
                # NOTE (uncited): the < 100K threshold is unchanged for the same reason as
                # the INR value above.
                tpa_eligible = False
                tpa_reasons.append(f"Platelet count {platelets}K < 100K")

            # Absolute contraindications. The set intersection silently discarded any
            # contraindication code not in the hard-coded list, so `recent_gi_bleed_21d`,
            # `recent_major_surgery_14d` and `infective_endocarditis` vanished without comment
            # and the patient was cleared (N3).
            matched_contra = contraindications & self._TPA_CONTRAINDICATIONS
            if matched_contra:
                tpa_eligible = False
                tpa_reasons.append(
                    f"Absolute contraindications: {', '.join(sorted(matched_contra))}"
                )
            unrecognised_contra = contraindications - self._TPA_CONTRAINDICATIONS
            if unrecognised_contra:
                tpa_eligible = False
                tpa_reasons.append(
                    f"UNRECOGNISED CONTRAINDICATION(S) reported and therefore NOT cleared: "
                    f"{', '.join(sorted(unrecognised_contra))}. A contraindication this "
                    f"workflow does not know is not an absent contraindication -- review "
                    f"each by hand."
                )

            # Minor stroke exclusion
            if nihss == 0 and nihss_complete:
                tpa_eligible = False
                tpa_reasons.append("NIHSS 0 -- no measurable deficit")

            if tpa_eligible:
                findings.append("tPA ELIGIBLE -- administer IV alteplase urgently")
                recommendations.append(
                    "Administer IV alteplase 0.9 mg/kg (max 90 mg): 10% bolus, "
                    "remainder over 60 min"
                )
                severity = SeverityLevel.CRITICAL
            else:
                findings.append(f"tPA NOT eligible: {'; '.join(tpa_reasons)}")

            guidelines.append(
                "AHA/ASA 2019: Guidelines for Early Management of Acute Ischemic Stroke"
            )

            # -- Thrombectomy eligibility --
            #
            # N6: entry was `if lvo or nihss >= 6`, and the standard-window branch then
            # re-checked only `nihss >= 6` -- a demonstrated occlusion was never required.
            # `{onset_hours: 3.0, NIHSS 10, age: 60}` with no CTA and no CT regions yielded
            # "THROMBECTOMY ELIGIBLE (standard window <6h, NIHSS >= 6)" and "Activate
            # neurointerventional team". Every 2015 trial and the AHA/ASA 2019 Class I
            # recommendation require a demonstrated ICA/M1 occlusion plus pre-stroke
            # functional independence. ASPECTS was permissive when absent
            # (`aspects_result is None or score >= 6`), so an unscored CT read as a
            # favourable one.
            thrombectomy_eligible = False
            thrombectomy_blockers: List[str] = []
            if lvo is None:
                thrombectomy_blockers.append(
                    "no vessel imaging result -- thrombectomy requires a demonstrated ICA/M1 "
                    "or basilar occlusion (obtain CTA/MRA)"
                )
            elif not lvo:
                thrombectomy_blockers.append("no large vessel occlusion demonstrated")
            if not nihss_complete:
                thrombectomy_blockers.append(
                    f"NIHSS incomplete ({nihss_result.items_scored}/"
                    f"{nihss_result.items_expected} items) -- the NIHSS >= 6 threshold cannot "
                    f"be applied to a partial examination"
                )
            if not aspects_scored:
                thrombectomy_blockers.append(
                    "ASPECTS not scored -- an unscored CT is not a favourable CT"
                )
            elif aspects_result.score < 6:
                thrombectomy_blockers.append(
                    f"ASPECTS {int(aspects_result.score)} < 6 -- large established infarct"
                )
            if prestroke_mrs is None:
                thrombectomy_blockers.append(
                    "pre-stroke mRS not recorded -- every thrombectomy trial enrolled on "
                    "pre-stroke functional independence"
                )
            elif prestroke_mrs >= 3:
                # Outside the enrollment range of every trial behind the recommendation:
                # the 2015 trials enrolled pre-stroke mRS 0-1, DAWN and DEFUSE-3 enrolled 0-2.
                thrombectomy_blockers.append(
                    f"pre-stroke mRS {prestroke_mrs} -- outside the enrollment range of the "
                    f"trials behind this recommendation (0-1 in the 2015 trials, 0-2 in "
                    f"DAWN/DEFUSE-3)"
                )
            elif prestroke_mrs == 2:
                findings.append(
                    "Pre-stroke mRS 2: within the DAWN/DEFUSE-3 range but outside the 2015 "
                    "trials' (0-1) -- confirm with the treating neurointerventionalist"
                )

            if thrombectomy_blockers:
                findings.append(
                    "THROMBECTOMY NOT CLEARED by this workflow: "
                    + "; ".join(thrombectomy_blockers)
                )
                recommendations.append(
                    "Resolve before thrombectomy: " + "; ".join(thrombectomy_blockers)
                )
            elif lvo or nihss >= 6:
                # Standard window: 0-6 hours. `onset_known` matters here too: with onset 0.0
                # standing for "unknown", `onset <= 6.0` was true for every unknown-onset
                # patient, so they were cleared on the standard window and the extended-window
                # branch below -- which says `or onset == 0` and routes exactly those patients
                # through DAWN/DEFUSE-3 imaging criteria -- could never be reached for anyone
                # with NIHSS >= 6.
                if onset_known and onset <= 6.0 and nihss >= 6:
                    # `lvo` is True and ASPECTS >= 6 by the blocker checks above.
                    thrombectomy_eligible = True
                    findings.append(
                        "THROMBECTOMY ELIGIBLE (standard window <6h, NIHSS >= 6, "
                        "LVO demonstrated, ASPECTS >= 6, pre-stroke mRS within trial range)"
                    )
                    recommendations.append(
                        "Activate neurointerventional team for mechanical thrombectomy"
                    )
                    guidelines.extend([
                        "MR CLEAN (2015): Endovascular treatment for acute ischemic stroke",
                        "ESCAPE (2015): Endovascular treatment guidelines",
                    ])

                # Extended window: 6-24 hours (DAWN/DEFUSE-3)
                elif (onset_known and 6.0 < onset <= 24.0) or not onset_known:
                    dawn_eligible = False
                    defuse3_eligible = False

                    # DAWN criteria. The three age/NIHSS/core strata are transcribed
                    # correctly; `age` is now required rather than defaulting to 0, which
                    # used to put every age-unknown patient into the `age < 80` strata.
                    if core_vol is not None and age is not None:
                        if age >= 80 and nihss >= 10 and core_vol < 21:
                            dawn_eligible = True
                        elif age < 80 and nihss >= 10 and core_vol < 31:
                            dawn_eligible = True
                        elif age < 80 and nihss >= 20 and core_vol < 51:
                            dawn_eligible = True
                    elif core_vol is not None and age is None:
                        findings.append(
                            "DAWN criteria not evaluated: age not recorded, and all three "
                            "DAWN strata are age-dependent"
                        )

                    # DEFUSE-3 criteria.
                    #
                    # N10: this branch was reached for onset up to 24 h and tested
                    # `core_vol < 70 and mismatch >= 1.8` only, with no NIHSS re-check after
                    # the entry gate. `{onset_hours: 20.0, NIHSS 3, LVO, core 40, mismatch
                    # 2.0}` returned "THROMBECTOMY ELIGIBLE (extended window, DEFUSE-3
                    # criteria met)" -- contradicting the guideline string printed two lines
                    # below it, "DEFUSE-3 (2018): Thrombectomy 6-16h with perfusion mismatch".
                    #
                    # Primary source: DEFUSE-3 protocol / statistical analysis plan,
                    # NCT02586415
                    # (https://cdn.clinicaltrials.gov/large-docs/15/NCT02586415/SAP_000.pdf).
                    # Enrollment window 6-16 h; target mismatch = core < 70 mL AND ratio
                    # >= 1.8 AND mismatch volume >= 15 mL by RAPID; baseline NIHSS >= 6 and
                    # still >= 6 at randomization; ICA or M1 occlusion; pre-stroke mRS 0-2.
                    if onset_known and onset > 16.0:
                        findings.append(
                            "DEFUSE-3 not applicable: enrollment window is 6-16h "
                            "(NCT02586415) and onset is "
                            f"{onset:.1f}h -- DAWN (6-24h) is the extended-window criterion "
                            "beyond 16h"
                        )
                    elif core_vol is not None and mismatch is not None and nihss >= 6:
                        if mismatch_vol is None:
                            findings.append(
                                "DEFUSE-3 not evaluated: mismatch volume (mL) not reported, "
                                "and the trial's target-mismatch profile requires mismatch "
                                "volume >= 15 mL alongside core < 70 mL and ratio >= 1.8"
                            )
                        elif core_vol < 70 and mismatch >= 1.8 and mismatch_vol >= 15:
                            defuse3_eligible = True
                    elif core_vol is not None and mismatch is not None and nihss < 6:
                        findings.append(
                            f"DEFUSE-3 not met: baseline NIHSS {int(nihss)} is below the "
                            f"trial's NIHSS >= 6 inclusion criterion"
                        )

                    if dawn_eligible or defuse3_eligible:
                        thrombectomy_eligible = True
                        criteria = []
                        if dawn_eligible:
                            criteria.append("DAWN")
                        if defuse3_eligible:
                            criteria.append("DEFUSE-3")
                        findings.append(
                            f"THROMBECTOMY ELIGIBLE (extended window, "
                            f"{'/'.join(criteria)} criteria met)"
                        )
                        recommendations.append(
                            "Activate neurointerventional team -- extended window "
                            "thrombectomy per DAWN/DEFUSE-3"
                        )
                        if dawn_eligible:
                            guidelines.append(
                                "DAWN (2018): Thrombectomy 6-24h with clinical-imaging mismatch"
                            )
                        if defuse3_eligible:
                            guidelines.append(
                                "DEFUSE-3 (2018, NCT02586415): Thrombectomy 6-16h with "
                                "perfusion mismatch"
                            )
                    else:
                        findings.append(
                            "Extended window thrombectomy: DAWN/DEFUSE-3 criteria NOT met "
                            "(or perfusion data unavailable)"
                        )

            if not thrombectomy_eligible and nihss >= 6:
                recommendations.append(
                    "Obtain CTA head/neck urgently to evaluate for LVO"
                )

        elif stroke_type_str == "hemorrhagic":
            severity = SeverityLevel.CRITICAL
            findings.append("Intracerebral haemorrhage identified")
            recommendations.extend([
                "Reverse any anticoagulation urgently",
                "BP target <140 mmHg systolic (INTERACT2/ATACH-2)",
                "Neurosurgical consultation for possible evacuation",
                "CT angiography to evaluate for spot sign (active extravasation)",
            ])
            guidelines.append(
                "AHA/ASA 2022: Guidelines for Management of Spontaneous ICH"
            )

        elif stroke_type_str == "subarachnoid_hemorrhage" or stroke_type_str == "sah":
            severity = SeverityLevel.CRITICAL
            findings.append("Subarachnoid haemorrhage suspected/confirmed")
            recommendations.extend([
                "Emergent CT angiography to identify aneurysm",
                "Neurosurgical/neurointerventional consultation",
                "Nimodipine 60 mg q4h for 21 days (vasospasm prophylaxis)",
                "External ventricular drain if hydrocephalus present",
                "ICU admission with close neurological monitoring",
            ])
            guidelines.append(
                "AHA/ASA 2023: Guidelines for Management of Aneurysmal SAH"
            )

        # Hemorrhagic stroke or high NIHSS -> critical
        if nihss >= 21:
            severity = SeverityLevel.CRITICAL
        elif nihss >= 5:
            severity = _max_severity(severity, SeverityLevel.HIGH)

        return WorkflowResult(
            workflow_type=self.workflow_type,
            domain=self.domain,
            findings=findings,
            scale_results=scale_results,
            recommendations=recommendations,
            guideline_references=guidelines,
            severity=severity,
        )


# ===================================================================
# WORKFLOW 2 -- Dementia Evaluation
# ===================================================================


class DementiaEvaluationWorkflow(BaseNeuroWorkflow):
    """Dementia evaluation: atrophy pattern -> differential (AD vs FTD vs DLB
    vs vascular vs PSP vs MSA), ATN biomarker staging, MoCA/MMSE scoring,
    anti-amyloid eligibility (lecanemab/donanemab).

    Inputs
    ------
    moca_domain_scores : dict
        Per-domain MoCA scores.
    education_years : int
        Years of education for MoCA adjustment.
    mmse_score : int | None
        MMSE total if available.
    atrophy_pattern : str | None
        MRI atrophy pattern: 'hippocampal', 'frontotemporal', 'parietal',
        'diffuse', 'periventricular_wmh', 'caudate', 'midbrain'.
    amyloid_status : str | None
        'positive', 'negative', or None.
    tau_status : str | None
        'positive', 'negative', or None.
    neurodegeneration_status : str | None
        'positive', 'negative', or None.
    visual_hallucinations : bool
    parkinsonism : bool
    rem_sleep_disorder : bool
    cognitive_fluctuations : bool
    behavioral_changes : bool
        Personality change, disinhibition, apathy.
    language_variant : str | None
        'semantic', 'nonfluent', 'logopenic' for PPA.
    vascular_risk_factors : list[str]
    step_wise_decline : bool
    age : int
    apoe_genotype : str | None
        E.g. 'e3/e4', 'e4/e4'.
    """

    workflow_type = NeuroWorkflowType.DEMENTIA_EVALUATION
    domain = NeuroDomain.NEURODEGENERATIVE

    # Atrophy pattern to suspected subtype mapping
    _ATROPHY_MAP: Dict[str, DementiaSubtype] = {
        "hippocampal": DementiaSubtype.ALZHEIMERS,
        "medial_temporal": DementiaSubtype.ALZHEIMERS,
        "frontotemporal": DementiaSubtype.FRONTOTEMPORAL,
        "frontal": DementiaSubtype.FRONTOTEMPORAL,
        "temporal_anterior": DementiaSubtype.FRONTOTEMPORAL,
        "parietal": DementiaSubtype.ALZHEIMERS,
        "posterior_cortical": DementiaSubtype.ALZHEIMERS,
        "periventricular_wmh": DementiaSubtype.VASCULAR,
        "lacunar": DementiaSubtype.VASCULAR,
        "caudate": DementiaSubtype.FRONTOTEMPORAL,
        "midbrain": DementiaSubtype.PSP,
        "hummingbird_sign": DementiaSubtype.PSP,
        "hot_cross_bun": DementiaSubtype.MSA,
        "pontocerebellar": DementiaSubtype.MSA,
        "diffuse": DementiaSubtype.MIXED,
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        if "moca_domain_scores" not in inp and "mmse_score" not in inp:
            warnings.append("No cognitive test scores provided")
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        guidelines: List[str] = []
        scale_results: List[ScaleResult] = []
        severity = SeverityLevel.MODERATE

        # -- MoCA scoring --
        moca_domains = inputs.get("moca_domain_scores")
        if moca_domains:
            # Education defaulted to 13 -- a made-up demographic, chosen so that the
            # adjustment did not fire. None now means "not recorded", which the calculator
            # reports rather than assuming.
            edu = inputs.get("education_years")
            moca_result = MoCACalculator.calculate(moca_domains, edu)
            scale_results.append(moca_result)
            findings.append(
                f"MoCA: {int(moca_result.score)}/30 ({moca_result.severity_category})"
            )
            for warning in moca_result.warnings:
                findings.append(f"MoCA data: {warning}")

        mmse = inputs.get("mmse_score")
        if mmse is not None:
            findings.append(f"MMSE: {mmse}/30")
            if mmse < 24:
                findings.append("MMSE below dementia threshold (<24)")

        # -- Atrophy pattern differential --
        atrophy = inputs.get("atrophy_pattern")
        primary_subtype: Optional[DementiaSubtype] = None
        differential: List[DementiaSubtype] = []

        if atrophy:
            primary_subtype = self._ATROPHY_MAP.get(atrophy.lower())
            if primary_subtype:
                findings.append(
                    f"Atrophy pattern '{atrophy}' suggests: {primary_subtype.value}"
                )

        # -- Clinical features for DLB --
        vis_halluc = inputs.get("visual_hallucinations", False)
        parkinsonism = inputs.get("parkinsonism", False)
        rem_sleep = inputs.get("rem_sleep_disorder", False)
        cog_fluctuations = inputs.get("cognitive_fluctuations", False)

        dlb_core_features = sum([vis_halluc, parkinsonism, rem_sleep, cog_fluctuations])
        if dlb_core_features >= 2:
            findings.append(
                f"Probable DLB: {dlb_core_features}/4 core features present "
                "(visual hallucinations, parkinsonism, REM sleep behavior disorder, "
                "cognitive fluctuations)"
            )
            differential.append(DementiaSubtype.LEWY_BODY)
            if not primary_subtype:
                primary_subtype = DementiaSubtype.LEWY_BODY
        elif dlb_core_features == 1:
            findings.append("Possible DLB: 1 core feature present")
            differential.append(DementiaSubtype.LEWY_BODY)

        # -- FTD features --
        behavioral = inputs.get("behavioral_changes", False)
        lang_variant = inputs.get("language_variant")
        if behavioral:
            findings.append(
                "Behavioral changes (disinhibition/apathy) -- consider behavioral variant FTD"
            )
            differential.append(DementiaSubtype.FRONTOTEMPORAL)
        if lang_variant:
            findings.append(f"Language variant: {lang_variant} -- primary progressive aphasia")
            differential.append(DementiaSubtype.FRONTOTEMPORAL)

        # -- Vascular features --
        vasc_rf = inputs.get("vascular_risk_factors", [])
        step_wise = inputs.get("step_wise_decline", False)
        if step_wise or len(vasc_rf) >= 3:
            findings.append(
                "Step-wise decline and/or multiple vascular risk factors -- "
                "consider vascular cognitive impairment"
            )
            differential.append(DementiaSubtype.VASCULAR)

        # -- ATN biomarker staging --
        amyloid = inputs.get("amyloid_status")
        tau = inputs.get("tau_status")
        neurodegeneration = inputs.get("neurodegeneration_status")

        if amyloid or tau or neurodegeneration:
            # N11: `a/t/n = "+" if status == "positive" else "-"` printed an analyte that was
            # never ordered as NEGATIVE, and the code then asserted "A+T- = Alzheimer's
            # pathological change (amyloid-positive, pre-tangle stage)" -- a staging
            # conclusion drawn from a tau test that was never run. `?` marks unmeasured, which
            # is what the route twin already did.
            def _atn(status) -> str:
                if status in ("positive", "+", "abnormal"):
                    return "+"
                if status in ("negative", "-", "normal"):
                    return "-"
                return "?"

            a, t, n = _atn(amyloid), _atn(tau), _atn(neurodegeneration)
            atn_str = f"A{a}T{t}N{n}"
            unmeasured = [
                label for label, mark in (("amyloid", a), ("tau", t), ("neurodegeneration", n))
                if mark == "?"
            ]
            findings.append(
                f"ATN biomarker profile: {atn_str}"
                + (f" ('?' = not measured: {', '.join(unmeasured)})" if unmeasured else "")
            )

            if a == "+":
                if t == "+" and n == "+":
                    findings.append(
                        "A+T+N+ = Alzheimer's disease with neurodegeneration "
                        "(full AD pathological cascade)"
                    )
                elif t == "+":
                    findings.append(
                        "A+T+ = Alzheimer's disease (amyloid and tau pathology confirmed)"
                    )
                elif t == "-":
                    findings.append(
                        "A+T- = Alzheimer's pathological change "
                        "(amyloid-positive, pre-tangle stage)"
                    )
                else:
                    # The A/T/N narrative is a staging statement about two analytes. With tau
                    # unmeasured there is no stage to state.
                    findings.append(
                        "Amyloid-positive with tau NOT MEASURED -- the ATN stage cannot be "
                        "assigned. Order CSF p-tau or tau PET to distinguish Alzheimer's "
                        "pathological change (A+T-) from Alzheimer's disease (A+T+)."
                    )
                if not primary_subtype:
                    primary_subtype = DementiaSubtype.ALZHEIMERS

            guidelines.append(
                "NIA-AA 2018: ATN Research Framework for Alzheimer's disease"
            )

        # -- Anti-amyloid eligibility --
        age = inputs.get("age")
        # An incomplete MoCA is not a cognitive stage: its total is an upper bound only, so it
        # must not be the number an anti-amyloid decision is cleared against.
        if moca_domains and moca_result.complete:
            moca_score = moca_result.score
        elif mmse is not None:
            moca_score = mmse
        else:
            moca_score = None

        if amyloid == "positive":
            # N11: `eligible` started True and was only ever cleared by a cognitive score
            # below 18. With no MoCA and no MMSE supplied, `moca_score is None` and the check
            # was skipped entirely, so `{amyloid_status: "positive", age: 72}` produced
            # "ANTI-AMYLOID THERAPY: Patient may be eligible for lecanemab (Leqembi) or
            # donanemab (Kisunla)" from an amyloid result alone -- with "[INPUT WARNING] No
            # cognitive test scores provided" sitting in the same findings list. Age >90 and
            # APOE e4/e4 were demoted to caveats that did not affect eligibility.
            #
            # Anti-amyloid therapy now requires the inputs its own safety monitoring turns on:
            # a cognitive stage, anticoagulation status, and prior ICH / microhaemorrhage
            # count. Anything absent makes eligibility indeterminate, not met.
            eligible = True
            reasons: List[str] = []
            missing: List[str] = []

            if moca_score is None:
                missing.append("a cognitive stage (MoCA or MMSE)")
            # NOTE (uncited, unchanged): the MoCA < 18 floor has no citation in this repo. The
            # Lecanemab Appropriate Use Recommendations and the trial inclusion criteria are
            # specified on MMSE 22-30 and CDR-GS 0.5-1, not on MoCA; no primary source was
            # retrieved, so the value is left as it stands.
            elif moca_score < 18:
                eligible = False
                reasons.append("Cognitive score too low (moderate-to-severe dementia)")

            if age is None:
                missing.append("age")
            elif age > 90:
                # NOTE (uncited): the >90 cut-off is unchanged; it remains a caveat.
                reasons.append("Age >90 -- limited trial data")

            anticoag = inputs.get("on_anticoagulation")
            if anticoag is None:
                missing.append("anticoagulation status")
            elif anticoag:
                eligible = False
                reasons.append(
                    "On anticoagulation -- substantially increased intracerebral haemorrhage "
                    "risk with anti-amyloid therapy"
                )

            prior_ich = inputs.get("prior_ich")
            microhem = inputs.get("microhemorrhage_count")
            if prior_ich is None:
                missing.append("prior intracerebral haemorrhage (yes/no)")
            elif prior_ich:
                eligible = False
                reasons.append("Prior intracerebral haemorrhage")
            if microhem is None:
                missing.append("cerebral microhaemorrhage count on MRI")

            if missing:
                eligible = False
                findings.append(
                    "ANTI-AMYLOID THERAPY: eligibility INDETERMINATE -- not assessable "
                    "without " + "; ".join(missing) + ". An unmeasured input is not a "
                    "reassuring one."
                )
                recommendations.append(
                    "Before considering lecanemab or donanemab, obtain: "
                    + "; ".join(missing)
                )

            apoe = inputs.get("apoe_genotype", "")
            if "e4/e4" in str(apoe).lower():
                reasons.append(
                    "APOE e4/e4 homozygote -- higher ARIA risk with anti-amyloid therapy"
                )

            if eligible:
                findings.append(
                    "ANTI-AMYLOID THERAPY: Patient may be eligible for "
                    "lecanemab (Leqembi) or donanemab (Kisunla)"
                )
                recommendations.extend([
                    "Obtain baseline brain MRI (ARIA surveillance protocol)",
                    "Discuss risks/benefits of anti-amyloid therapy with patient/family",
                    "If initiating lecanemab: 10 mg/kg IV q2weeks with MRI monitoring "
                    "at weeks 7, 14, 52",
                ])
                if reasons:
                    findings.append(f"Eligibility caveats: {'; '.join(reasons)}")
            elif reasons:
                findings.append(
                    f"Anti-amyloid therapy not recommended: {'; '.join(reasons)}"
                )

            guidelines.append(
                "AAN 2023: Practice Guideline on Aducanumab, Lecanemab, and Donanemab"
            )

        # -- General recommendations --
        if not primary_subtype and differential:
            primary_subtype = differential[0]

        if primary_subtype:
            findings.append(f"Leading diagnosis: {primary_subtype.value}")

        recommendations.extend([
            "Obtain structural MRI brain with volumetric analysis if not done",
            "Complete reversible dementia workup: TSH, B12, folate, RPR, HIV, CMP",
        ])

        if primary_subtype == DementiaSubtype.ALZHEIMERS:
            recommendations.append(
                "Consider cholinesterase inhibitor (donepezil 5-10 mg) for mild-moderate AD"
            )
        elif primary_subtype == DementiaSubtype.LEWY_BODY:
            recommendations.extend([
                "AVOID antipsychotics (severe neuroleptic sensitivity in DLB)",
                "Consider cholinesterase inhibitor (rivastigmine preferred in DLB)",
            ])
            severity = _max_severity(severity, SeverityLevel.HIGH)
        elif primary_subtype == DementiaSubtype.FRONTOTEMPORAL:
            recommendations.append(
                "Consider SSRI for behavioral symptoms; refer to speech therapy for PPA"
            )

        guidelines.extend([
            "McKhann et al. 2011: NIA-AA criteria for Alzheimer's disease",
            "McKeith et al. 2017: DLB Consortium diagnostic criteria",
            "Rascovsky et al. 2011: International criteria for bvFTD",
        ])

        return WorkflowResult(
            workflow_type=self.workflow_type,
            domain=self.domain,
            findings=findings,
            scale_results=scale_results,
            recommendations=recommendations,
            guideline_references=guidelines,
            severity=severity,
        )


# ===================================================================
# WORKFLOW 3 -- Epilepsy Focus Localization
# ===================================================================


class EpilepsyFocusWorkflow(BaseNeuroWorkflow):
    """Epilepsy focus localization: ILAE seizure classification, syndrome
    identification, concordance matrix, ASM recommendations, Dravet SCN1A warning.

    Inputs
    ------
    seizure_semiology : str
        Description of seizure semiology.
    seizure_type : str
        ILAE seizure type (maps to SeizureType enum values).
    eeg_findings : dict
        Keys: 'interictal_focus' (str), 'ictal_onset' (str),
        'generalized_discharges' (bool), 'frequency_hz' (float | None).
    mri_findings : dict
        Keys: 'lesion_location' (str | None), 'lesion_type' (str | None),
        'mesial_temporal_sclerosis' (bool), 'fcd' (bool),
        'tumor' (bool), 'normal' (bool).
    pet_findings : dict | None
        Keys: 'hypometabolism_location' (str | None).
    neuropsych_lateralization : str | None
        'left', 'right', or None.
    age_at_onset : int
    seizure_frequency : str
        E.g. '2-3/month', 'daily'.
    current_asms : list[str]
        Current anti-seizure medications.
    num_failed_asms : int
        Number of adequately tried ASMs that failed.
    family_history : list[str]
        Relevant family history.
    genetic_testing : dict | None
        Keys: gene name -> variant status.
    """

    workflow_type = NeuroWorkflowType.EPILEPSY_FOCUS
    domain = NeuroDomain.EPILEPSY

    # ASM recommendations by seizure type
    _ASM_BY_SEIZURE_TYPE: Dict[str, List[str]] = {
        "focal_aware": [
            "Carbamazepine", "Lamotrigine", "Levetiracetam",
            "Oxcarbazepine", "Lacosamide",
        ],
        "focal_impaired_awareness": [
            "Carbamazepine", "Lamotrigine", "Levetiracetam",
            "Oxcarbazepine", "Lacosamide",
        ],
        "focal_to_bilateral_tonic_clonic": [
            "Carbamazepine", "Lamotrigine", "Levetiracetam",
            "Oxcarbazepine", "Lacosamide", "Brivaracetam",
        ],
        "generalized_tonic_clonic": [
            "Valproate", "Lamotrigine", "Levetiracetam",
            "Perampanel",
        ],
        "generalized_absence": [
            "Ethosuximide", "Valproate", "Lamotrigine",
        ],
        "generalized_myoclonic": [
            "Valproate", "Levetiracetam", "Clonazepam",
        ],
        "generalized_atonic": [
            "Valproate", "Lamotrigine", "Rufinamide", "Clobazam",
        ],
        "unknown_onset": [
            "Levetiracetam", "Lamotrigine", "Valproate",
        ],
    }

    # Dangerous ASMs for specific syndromes
    _CONTRAINDICATED_ASMS: Dict[str, List[str]] = {
        "dravet": [
            "Carbamazepine", "Oxcarbazepine", "Phenytoin", "Lamotrigine",
        ],
        "juvenile_myoclonic": [
            "Carbamazepine", "Oxcarbazepine", "Phenytoin",
        ],
        "childhood_absence": [
            "Carbamazepine", "Oxcarbazepine", "Phenytoin",
        ],
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        if "seizure_type" not in inp:
            inp["seizure_type"] = "unknown_onset"
            warnings.append("Seizure type not specified; defaulting to unknown onset")
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        guidelines: List[str] = []
        severity = SeverityLevel.MODERATE
        cross_triggers: List[str] = []

        seizure_type_val = inputs.get("seizure_type", "unknown_onset")
        findings.append(f"ILAE seizure classification: {seizure_type_val}")

        age_onset = inputs.get("age_at_onset")
        if age_onset is None:
            findings.append(
                "Age at seizure onset: NOT RECORDED -- several syndrome patterns below are "
                "age-gated and cannot be evaluated"
            )
        else:
            findings.append(f"Age at seizure onset: {age_onset} years")

        # -- Syndrome identification --
        eeg = inputs.get("eeg_findings", {})
        mri = inputs.get("mri_findings", {})
        gen_discharges = eeg.get("generalized_discharges", False)
        freq_hz = eeg.get("frequency_hz")

        syndrome: Optional[str] = None

        # A syndrome the clinician has already diagnosed.
        #
        # N19(a): `syndrome = "dravet"` was reachable ONLY through a pathogenic SCN1A result,
        # and the contraindication scan further down is gated on
        # `syndrome in _CONTRAINDICATED_ASMS`. A clinically diagnosed Dravet patient taking
        # carbamazepine, with no genetic testing in the payload, therefore got no
        # "DANGER: ... CONTRAINDICATED" finding and no severity escalation. Sodium-channel
        # blockers worsen Dravet seizures.
        clinical_dx = str(inputs.get("clinical_syndrome_diagnosis", "") or "").lower()
        for _name, _key in (
            ("dravet", "dravet"),
            ("juvenile myoclonic", "juvenile_myoclonic"),
            ("childhood absence", "childhood_absence"),
            ("lennox", "lennox_gastaut"),
            ("west", "west"),
            ("temporal lobe", "temporal_lobe"),
        ):
            if _name in clinical_dx:
                syndrome = _key
                findings.append(
                    f"Clinically diagnosed syndrome supplied: {_key} -- syndrome-specific "
                    f"drug contraindications apply regardless of genetic testing status"
                )
                if _key in ("dravet", "west", "lennox_gastaut"):
                    severity = SeverityLevel.CRITICAL
                break

        # Dravet syndrome check
        genetic = inputs.get("genetic_testing", {})
        scn1a = genetic.get("SCN1A", genetic.get("scn1a"))
        if scn1a and "pathogenic" in str(scn1a).lower():
            syndrome = "dravet"
            findings.append(
                "WARNING: SCN1A pathogenic variant detected -- Dravet syndrome"
            )
            severity = SeverityLevel.CRITICAL
            cross_triggers.append(
                _trigger_string("PHARMACOGENOMIC", ["SCN1A"],
                                "Dravet syndrome -- sodium channel blockers CONTRAINDICATED")
            )

        # JME pattern.
        # N19(b): the code matched `3 <= age_onset <= 25` while the finding it printed said
        # "onset 8-25y", and the route twin used 12-18 -- three ranges, two of them in the
        # same sentence. OPEN: the ILAE 2022 syndrome classification paper is the issuing
        # body's own definition and was not retrieved, so the matched range is left as it is
        # and the printed text is corrected to state the range actually applied.
        if (seizure_type_val in ("generalized_myoclonic", "generalized_tonic_clonic")
                and gen_discharges and age_onset is not None and 3 <= age_onset <= 25):
            if freq_hz and 3.5 <= freq_hz <= 6.0:
                syndrome = "juvenile_myoclonic"
                findings.append(
                    "Pattern consistent with Juvenile Myoclonic Epilepsy (JME): "
                    "myoclonic seizures, generalized discharges, onset age within the "
                    "3-25y range this workflow screens on (confirm against the ILAE 2022 "
                    "syndrome definition)"
                )

        # Childhood absence
        if (seizure_type_val == "generalized_absence"
                and gen_discharges and age_onset is not None and 4 <= age_onset <= 10):
            if freq_hz and 2.5 <= freq_hz <= 3.5:
                syndrome = "childhood_absence"
                findings.append(
                    "Pattern consistent with Childhood Absence Epilepsy: "
                    "3 Hz spike-and-wave, onset 4-10y"
                )

        # Temporal lobe epilepsy
        mts = mri.get("mesial_temporal_sclerosis", False)
        interictal = eeg.get("interictal_focus", "")
        if mts or "temporal" in str(interictal).lower():
            if not syndrome:
                syndrome = "temporal_lobe"
                findings.append("Temporal lobe epilepsy pattern identified")

        # West syndrome
        if (age_onset is not None and age_onset < 2
                and seizure_type_val in ("generalized_tonic", "unknown_onset")):
            hyps = "hypsarrhythmia" in str(eeg).lower()
            if hyps:
                syndrome = "west"
                findings.append("West syndrome (infantile spasms with hypsarrhythmia)")
                severity = SeverityLevel.CRITICAL

        # Lennox-Gastaut
        if (seizure_type_val in ("generalized_tonic", "generalized_atonic")
                and age_onset is not None and age_onset < 8 and gen_discharges):
            if freq_hz and freq_hz < 2.5:
                syndrome = "lennox_gastaut"
                findings.append("Lennox-Gastaut syndrome pattern: slow spike-wave <2.5 Hz")
                severity = SeverityLevel.HIGH

        if syndrome:
            findings.append(f"Epilepsy syndrome: {syndrome}")

        # -- Concordance matrix --
        concordance_zones: Dict[str, str] = {}
        if eeg.get("interictal_focus"):
            concordance_zones["EEG_interictal"] = eeg["interictal_focus"]
        if eeg.get("ictal_onset"):
            concordance_zones["EEG_ictal"] = eeg["ictal_onset"]
        if mri.get("lesion_location"):
            concordance_zones["MRI_lesion"] = mri["lesion_location"]
        pet = inputs.get("pet_findings", {})
        if pet and pet.get("hypometabolism_location"):
            concordance_zones["PET_hypometabolism"] = pet["hypometabolism_location"]
        neuropsych_lat = inputs.get("neuropsych_lateralization")
        if neuropsych_lat:
            concordance_zones["Neuropsych"] = neuropsych_lat

        if concordance_zones:
            findings.append(f"Concordance matrix: {concordance_zones}")
            # Check if all zones agree
            set(
                v.lower().replace("left ", "").replace("right ", "")
                for v in concordance_zones.values()
            )
            lateralities = set()
            for v in concordance_zones.values():
                vl = v.lower()
                if "left" in vl:
                    lateralities.add("left")
                elif "right" in vl:
                    lateralities.add("right")

            if len(concordance_zones) >= 3:
                if len(lateralities) == 1:
                    findings.append(
                        f"Concordant lateralization: {lateralities.pop()} hemisphere"
                    )
                elif len(lateralities) > 1:
                    findings.append(
                        "DISCORDANT lateralization across modalities -- "
                        "consider Phase II (intracranial EEG) monitoring"
                    )

        # -- ASM recommendations --
        current_asms = [a.lower() for a in inputs.get("current_asms", [])]
        num_failed = inputs.get("num_failed_asms", 0)

        recommended_asms = self._ASM_BY_SEIZURE_TYPE.get(
            seizure_type_val,
            self._ASM_BY_SEIZURE_TYPE["unknown_onset"],
        )

        # Filter out contraindicated ASMs
        contraindicated = []
        if syndrome and syndrome in self._CONTRAINDICATED_ASMS:
            contraindicated = self._CONTRAINDICATED_ASMS[syndrome]
            safe_asms = [a for a in recommended_asms if a not in contraindicated]
            # Warn if patient is on a contraindicated ASM
            for asm in current_asms:
                for contra in contraindicated:
                    if contra.lower() in asm:
                        findings.append(
                            f"DANGER: {contra} is CONTRAINDICATED in {syndrome} -- "
                            "may worsen seizures"
                        )
                        severity = _max_severity(severity, SeverityLevel.CRITICAL)
        else:
            safe_asms = recommended_asms

        # Dravet-specific ASMs
        if syndrome == "dravet":
            safe_asms = [
                "Stiripentol + Clobazam + Valproate (Dravet first-line)",
                "Cannabidiol (Epidiolex)",
                "Fenfluramine (Fintepla)",
            ]

        new_asms = [a for a in safe_asms if a.lower() not in current_asms]
        if new_asms:
            recommendations.append(
                f"Consider ASM options: {', '.join(new_asms[:4])}"
            )

        if contraindicated:
            recommendations.append(
                f"AVOID in this syndrome: {', '.join(contraindicated)}"
            )

        # Drug-resistant epilepsy
        if num_failed >= 2:
            findings.append(
                f"Drug-resistant epilepsy: {num_failed} adequate ASM trials failed "
                "(ILAE definition: >= 2)"
            )
            recommendations.append(
                "Refer to comprehensive epilepsy center for surgical evaluation"
            )
            if concordance_zones and len(concordance_zones) >= 2:
                recommendations.append(
                    "Consider Phase I presurgical evaluation: video-EEG, "
                    "3T epilepsy-protocol MRI, neuropsych, PET, +/- MEG"
                )

        guidelines.extend([
            "ILAE 2017: Operational Classification of Seizure Types",
            "AAN/AES 2018: Practice guideline on when to initiate ASM treatment",
            "Kwan & Brodie 2000: Definition of drug-resistant epilepsy",
        ])

        return WorkflowResult(
            workflow_type=self.workflow_type,
            domain=self.domain,
            findings=findings,
            recommendations=recommendations,
            guideline_references=guidelines,
            severity=severity,
            cross_modal_triggers=cross_triggers,
        )


# ===================================================================
# WORKFLOW 4 -- Brain Tumor Grading
# ===================================================================


class BrainTumorGradingWorkflow(BaseNeuroWorkflow):
    """Brain tumor grading: WHO 2021 classification, molecular markers
    (IDH, MGMT, 1p/19q, H3K27M), treatment protocol (Stupp for GBM),
    RANO criteria.

    Inputs
    ------
    histology : str
        Histological diagnosis (e.g. 'glioblastoma', 'oligodendroglioma',
        'astrocytoma', 'meningioma', 'medulloblastoma').
    molecular_markers : list[str]
        List of marker identifiers (maps to TumorMolecularMarker values).
    location : str
        Anatomical tumor location.
    size_cm : float
        Maximum dimension in cm.
    age : int
        Patient age.
    kps : int
        Karnofsky Performance Status (0-100).
    prior_treatments : list[str]
        Previous treatments.
    enhancing : bool
        Contrast-enhancing on MRI.
    recurrence : bool
        Whether this is recurrent disease.
    """

    workflow_type = NeuroWorkflowType.BRAIN_TUMOR
    domain = NeuroDomain.NEURO_ONCOLOGY

    # WHO 2021 integrated diagnosis rules (simplified)
    _CLASSIFICATION_RULES = {
        "glioblastoma": {
            "default_grade": TumorGrade.GRADE_4,
            "required_markers": [TumorMolecularMarker.IDH_WILDTYPE],
            "who_name": "Glioblastoma, IDH-wildtype (WHO Grade 4)",
        },
        "astrocytoma": {
            "default_grade": TumorGrade.GRADE_2,
            "preferred_markers": [TumorMolecularMarker.IDH_MUTANT],
            "who_name": "Astrocytoma, IDH-mutant",
        },
        "oligodendroglioma": {
            "default_grade": TumorGrade.GRADE_2,
            "required_markers": [
                TumorMolecularMarker.IDH_MUTANT,
                TumorMolecularMarker.CODEL_1P19Q,
            ],
            "who_name": "Oligodendroglioma, IDH-mutant and 1p/19q-codeleted",
        },
        "meningioma": {
            "default_grade": TumorGrade.GRADE_1,
            "who_name": "Meningioma",
        },
        "medulloblastoma": {
            "default_grade": TumorGrade.GRADE_4,
            "who_name": "Medulloblastoma",
        },
        "ependymoma": {
            "default_grade": TumorGrade.GRADE_2,
            "who_name": "Ependymoma",
        },
        "dmg": {
            "default_grade": TumorGrade.GRADE_4,
            "required_markers": [TumorMolecularMarker.H3K27M_MUTANT],
            "who_name": "Diffuse midline glioma, H3 K27-altered (WHO Grade 4)",
        },
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        if "histology" not in inp:
            warnings.append("No histological diagnosis provided")
            inp["histology"] = "unknown"
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        guidelines: List[str] = []
        severity = SeverityLevel.HIGH

        histology = inputs.get("histology", "unknown").lower()
        marker_strs = inputs.get("molecular_markers", [])
        location = inputs.get("location", "not specified")
        # N12: `size = inputs.get("size_cm", 0.0)` and `kps = inputs.get("kps", 100)` were
        # both PRINTED AS FINDINGS, so a GBM case with neither supplied reported
        # "Location: not specified, Size: 0.0 cm" and "KPS: 100". The fabricated 100 also kept
        # the `kps < 60` palliative-care recommendation from ever firing, and for meningioma
        # `if size < 3.0 and not enhancing` gave an unmeasured tumour "Observation with serial
        # MRI ... for small, asymptomatic meningioma".
        size = inputs.get("size_cm")
        age = inputs.get("age")
        kps = inputs.get("kps")
        enhancing = inputs.get("enhancing", False)
        recurrence = inputs.get("recurrence", False)
        prior_tx = inputs.get("prior_treatments", [])

        # Parse molecular markers
        markers = set()
        for m in marker_strs:
            try:
                markers.add(TumorMolecularMarker(m))
            except ValueError:
                pass

        # -- WHO 2021 classification --
        rule = self._CLASSIFICATION_RULES.get(histology)
        if rule:
            who_name = rule["who_name"]
            grade = rule["default_grade"]

            # Grade adjustment for astrocytoma based on markers
            if histology == "astrocytoma":
                if TumorMolecularMarker.IDH_MUTANT in markers:
                    # Check for higher grade features.
                    # OPEN (flagged, unchanged): neither the `enhancing or size > 5.0` upgrade
                    # here nor the `atrx_loss` upgrade in the route twin appears in WHO CNS5,
                    # which uses histology plus homozygous CDKN2A/B deletion for grade 4. The
                    # WHO Classification of Tumours of the CNS, 5th edition was not
                    # accessible, so the rule is left as it is rather than replaced from
                    # memory. `size` being None no longer makes the comparison crash or
                    # silently read as 0.
                    if enhancing or (size is not None and size > 5.0):
                        grade = TumorGrade.GRADE_3
                        who_name = "Astrocytoma, IDH-mutant, WHO Grade 3"
                    else:
                        who_name = "Astrocytoma, IDH-mutant, WHO Grade 2"
                elif TumorMolecularMarker.IDH_WILDTYPE in markers:
                    # IDH-wildtype astrocytoma -> classify as GBM per WHO 2021
                    if (TumorMolecularMarker.TERT_MUTANT in markers
                            or TumorMolecularMarker.EGFR_AMPLIFIED in markers):
                        grade = TumorGrade.GRADE_4
                        who_name = (
                            "Glioblastoma, IDH-wildtype (reclassified from astrocytoma per "
                            "WHO 2021 molecular criteria: TERT/EGFR)"
                        )

            findings.append(f"WHO 2021 integrated diagnosis: {who_name}")
            findings.append(f"WHO Grade: {grade.value}")
        else:
            grade = TumorGrade.GRADE_2  # default
            findings.append(f"Histology: {histology} (not in WHO 2021 lookup)")

        findings.append(
            f"Location: {location}, Size: "
            + (f"{size:.1f} cm" if size is not None else "NOT MEASURED")
        )
        if kps is None:
            findings.append(
                "KPS: NOT RECORDED -- performance status drives treatment intensity and the "
                "palliative-care decision; record it"
            )
            recommendations.append(
                "Record Karnofsky Performance Status: treatment intensity cannot be "
                "recommended without it"
            )
        else:
            findings.append(f"KPS: {kps}")
        if size is None:
            recommendations.append(
                "Record the maximum tumour diameter: surveillance-versus-resection decisions "
                "depend on it"
            )

        # -- Molecular marker reporting --
        # MGMT status is a temozolomide-response marker, so this applies to the diffuse
        # gliomas, not to every histology.
        if histology in ("glioblastoma", "astrocytoma", "oligodendroglioma", "dmg") and not (
            {TumorMolecularMarker.MGMT_METHYLATED,
             TumorMolecularMarker.MGMT_UNMETHYLATED} & markers
        ):
            findings.append(
                "MGMT promoter methylation: NOT TESTED -- absence of a result is not an "
                "unmethylated result, and temozolomide benefit cannot be predicted without it"
            )
            recommendations.append("Order MGMT promoter methylation testing")
        for marker in markers:
            if marker == TumorMolecularMarker.MGMT_METHYLATED:
                findings.append(
                    "MGMT promoter: METHYLATED -- predicts better response to temozolomide"
                )
            elif marker == TumorMolecularMarker.MGMT_UNMETHYLATED:
                findings.append(
                    "MGMT promoter: UNMETHYLATED -- reduced temozolomide benefit"
                )
            elif marker == TumorMolecularMarker.IDH_MUTANT:
                findings.append("IDH: MUTANT -- favorable prognostic factor")
            elif marker == TumorMolecularMarker.IDH_WILDTYPE:
                findings.append("IDH: WILDTYPE -- aggressive biology")
            elif marker == TumorMolecularMarker.CODEL_1P19Q:
                findings.append(
                    "1p/19q: CO-DELETED -- oligodendroglial lineage, "
                    "predicts chemosensitivity (PCV)"
                )
            elif marker == TumorMolecularMarker.H3K27M_MUTANT:
                findings.append(
                    "H3 K27M: MUTANT -- diffuse midline glioma, universally Grade 4, "
                    "poor prognosis"
                )
                severity = SeverityLevel.CRITICAL

        # -- Treatment protocol --
        if grade == TumorGrade.GRADE_4 and histology in ("glioblastoma", "astrocytoma", "dmg"):
            severity = SeverityLevel.CRITICAL
            if not recurrence:
                recommendations.extend([
                    "Stupp protocol: maximal safe resection followed by concurrent "
                    "temozolomide + radiation (60 Gy/30 fractions) then adjuvant "
                    "temozolomide x6 cycles",
                    "Consider tumor-treating fields (TTFields/Optune) during adjuvant phase",
                ])
                if (TumorMolecularMarker.MGMT_UNMETHYLATED in markers
                        and age is not None and age > 65):
                    recommendations.append(
                        "Consider hypofractionated RT (40 Gy/15 fx) +/- temozolomide "
                        "given age >65 and MGMT unmethylated status"
                    )
            else:
                recommendations.extend([
                    "Recurrent GBM: consider bevacizumab (Avastin) or lomustine (CCNU)",
                    "Evaluate clinical trial eligibility (immunotherapy, targeted therapy)",
                    "Repeat molecular profiling on recurrent tissue if available",
                ])

        elif histology == "oligodendroglioma":
            recommendations.extend([
                "PCV chemotherapy (procarbazine, CCNU, vincristine) is standard "
                "for 1p/19q-codeleted oligodendroglioma",
                "Radiation therapy 54 Gy/30 fractions for Grade 3 or symptomatic Grade 2",
            ])

        elif histology == "meningioma":
            if grade == TumorGrade.GRADE_1:
                if size is None:
                    recommendations.append(
                        "Meningioma management cannot be recommended: the tumour was not "
                        "measured. Observation applies to small, asymptomatic lesions, and an "
                        "unmeasured lesion is not a small one."
                    )
                elif size < 3.0 and not enhancing:
                    recommendations.append(
                        "Observation with serial MRI (6 months, then annually) "
                        "for small, asymptomatic meningioma"
                    )
                else:
                    recommendations.append(
                        "Consider surgical resection (Simpson Grade I-III if feasible)"
                    )
            severity = SeverityLevel.MODERATE

        # -- RANO criteria for treatment response --
        if recurrence or prior_tx:
            findings.append("RANO criteria applicable for treatment response assessment")
            recommendations.append(
                "Follow-up MRI with RANO-based response assessment: "
                "measure enhancing and non-enhancing (T2/FLAIR) tumor components"
            )

        # General recommendations
        recommendations.extend([
            "Multidisciplinary tumor board review (neurosurgery, neuro-oncology, "
            "radiation oncology, neuropathology)",
            "Seizure prophylaxis with levetiracetam if presenting with seizures",
        ])

        # NOTE: this workflow uses KPS < 60 while api/routes/neuro_clinical.py uses KPS < 70
        # for the same decision. Neither cut-point is cited in this repo and no primary source
        # was retrieved, so both values are left as they are and the divergence is recorded
        # rather than silently resolved in one direction.
        if kps is not None and kps < 60:
            recommendations.append(
                "KPS <60: consider palliative/supportive care focus; "
                "limited benefit from aggressive treatment"
            )

        guidelines.extend([
            "WHO 2021: Classification of Tumors of the Central Nervous System (5th ed.)",
            "Stupp et al. 2005/2017: Temozolomide + RT for newly diagnosed GBM",
            "RANO Working Group 2010: Updated response criteria for high-grade gliomas",
            "NCCN CNS Cancers Guidelines v2.2025",
        ])

        return WorkflowResult(
            workflow_type=self.workflow_type,
            domain=self.domain,
            findings=findings,
            recommendations=recommendations,
            guideline_references=guidelines,
            severity=severity,
        )


# ===================================================================
# WORKFLOW 5 -- MS Monitoring
# ===================================================================


class MSMonitoringWorkflow(BaseNeuroWorkflow):
    """MS monitoring: McDonald 2017 criteria, new T2/enhancing lesions,
    NEDA-3 status, DMT escalation triggers, PML risk (JCV index),
    brain volume change.

    Inputs
    ------
    phenotype : str
        MS phenotype: 'rrms', 'spms', 'ppms', 'cis'.
    edss_current : float
        Current EDSS score.
    edss_prior : float | None
        EDSS at last assessment.
    edss_fs_scores : dict | None
        Functional system scores for EDSS.
    new_t2_lesions : int
        Count of new/enlarging T2 lesions on MRI.
    gad_enhancing : int
        Count of gadolinium-enhancing lesions.
    relapses_last_year : int
        Number of relapses in the past 12 months.
    relapses_last_2yr : int
        Number of relapses in the past 24 months.
    current_dmt : str | None
        Current disease-modifying therapy name.
    dmt_category : str | None
        'platform', 'moderate_efficacy', 'high_efficacy'.
    jcv_antibody_index : float | None
        JCV antibody index value.
    brain_volume_change_pct : float | None
        Annualised brain volume change (%, negative = atrophy).
    nfl_level : float | None
        Serum neurofilament light chain (pg/mL).
    disease_duration_years : float
    age : int
    """

    workflow_type = NeuroWorkflowType.MS_MONITORING
    domain = NeuroDomain.MS

    # DMT escalation hierarchy
    _PLATFORM_DMTS = {
        "interferon_beta", "glatiramer_acetate", "teriflunomide",
        "dimethyl_fumarate", "diroximel_fumarate",
    }
    _MODERATE_DMTS = {
        "fingolimod", "siponimod", "ozanimod", "ponesimod", "cladribine",
    }
    _HIGH_EFFICACY_DMTS = {
        "natalizumab", "ocrelizumab", "ofatumumab", "rituximab",
        "alemtuzumab", "ublituximab",
    }

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        guidelines: List[str] = []
        scale_results: List[ScaleResult] = []
        severity = SeverityLevel.MODERATE

        phenotype = inputs.get("phenotype", "rrms")
        findings.append(f"MS phenotype: {phenotype.upper()}")

        # N7(a): these all defaulted to 0 / 0.0 and `neda3` started True, so
        # `MSMonitoringWorkflow().run({})` returned "EDSS: 0.0", "NEDA-3 STATUS: ACHIEVED --
        # no disease activity detected" and severity LOW: a patient with no MRI, no relapse
        # history and no recorded EDSS was certified disease-free. NEDA-3 means no EVIDENCE of
        # disease activity, which requires evidence.
        edss_current = inputs.get("edss_current")
        edss_prior = inputs.get("edss_prior")
        new_t2 = inputs.get("new_t2_lesions")
        gad_enh = inputs.get("gad_enhancing")
        relapses_1yr = inputs.get("relapses_last_year")
        relapses_2yr = inputs.get("relapses_last_2yr", 0)
        current_dmt = inputs.get("current_dmt")
        dmt_cat = inputs.get("dmt_category")
        jcv_index = inputs.get("jcv_antibody_index")
        bv_change = inputs.get("brain_volume_change_pct")
        nfl = inputs.get("nfl_level")

        # -- EDSS calculation --
        fs_scores = inputs.get("edss_fs_scores")
        if fs_scores or edss_current is not None:
            edss_result = EDSSCalculator.calculate(fs_scores or {}, edss_current)
            scale_results.append(edss_result)
        if edss_current is None:
            findings.append("EDSS: NOT RECORDED")
        else:
            findings.append(f"EDSS: {edss_current}")

        # -- NEDA-3 assessment --
        # NEDA-3: No Evidence of Disease Activity
        # 1. No relapses
        # 2. No new/enlarging T2 or Gd+ lesions
        # 3. No EDSS progression
        #
        # Each limb is now three-valued: True (clear), False (activity), None (not assessed).
        # NEDA-3 is achieved only when all three are True, and is "not assessable" -- not
        # achieved -- when any is None.
        neda_failures: List[str] = []
        neda_unassessed: List[str] = []

        if relapses_1yr is None:
            neda_unassessed.append("relapse count for the past year")
            relapse_limb = None
        elif relapses_1yr > 0:
            relapse_limb = False
            neda_failures.append(f"{relapses_1yr} relapse(s) in past year")
        else:
            relapse_limb = True

        if new_t2 is None and gad_enh is None:
            neda_unassessed.append("MRI lesion activity (new T2 / Gd-enhancing)")
            mri_limb = None
        elif (new_t2 or 0) > 0 or (gad_enh or 0) > 0:
            mri_limb = False
            neda_failures.append(
                f"{new_t2 or 0} new T2 lesion(s), {gad_enh or 0} Gd-enhancing lesion(s)"
            )
        else:
            mri_limb = True

        edss_progression = False
        if edss_prior is None or edss_current is None:
            neda_unassessed.append("EDSS progression (requires a prior and a current EDSS)")
            progression_limb = None
        else:
            delta = edss_current - edss_prior
            if edss_prior <= 5.5 and delta >= 1.0:
                edss_progression = True
            elif edss_prior >= 6.0 and delta >= 0.5:
                edss_progression = True
            progression_limb = not edss_progression
            if edss_progression:
                neda_failures.append(
                    f"EDSS progression: {edss_prior} -> {edss_current} (+{delta:.1f})"
                )

        limbs = (relapse_limb, mri_limb, progression_limb)
        if any(limb is False for limb in limbs):
            neda3 = False
            findings.append(
                f"NEDA-3 STATUS: NOT MET -- {'; '.join(neda_failures)}"
            )
            severity = SeverityLevel.HIGH
        elif all(limb is True for limb in limbs):
            neda3 = True
            findings.append("NEDA-3 STATUS: ACHIEVED -- no disease activity detected")
            severity = SeverityLevel.LOW
        else:
            neda3 = None
            findings.append(
                "NEDA-3 STATUS: NOT ASSESSABLE -- missing "
                + "; ".join(neda_unassessed)
                + ". Absent data is not absence of disease activity."
            )
            recommendations.append(
                "Record, before reporting disease-activity status: "
                + "; ".join(neda_unassessed)
            )

        # -- Brain volume --
        if bv_change is not None:
            findings.append(
                f"Annualised brain volume change: {bv_change:+.2f}%"
            )
            if bv_change < -0.4:
                findings.append(
                    "Accelerated brain atrophy (>0.4%/year exceeds normal ageing)"
                )
                severity = _max_severity(severity, SeverityLevel.HIGH)

        # -- NFL --
        if nfl is not None:
            findings.append(f"Serum neurofilament light: {nfl:.1f} pg/mL")
            if nfl > 16.0:
                findings.append("Elevated sNfL -- suggests ongoing neuroaxonal damage")

        # -- DMT escalation triggers --
        if current_dmt:
            findings.append(f"Current DMT: {current_dmt}")
            escalate = False
            escalation_reasons: List[str] = []

            # `neda3` is now three-valued; `not None` would have made "not assessed" read as
            # "activity present".
            if neda3 is False:
                escalate = True
                escalation_reasons.append("Disease activity despite current DMT")

            if relapses_2yr >= 2:
                escalate = True
                escalation_reasons.append(f"{relapses_2yr} relapses in 2 years")

            if new_t2 is not None and new_t2 >= 2:
                escalate = True
                escalation_reasons.append(f"{new_t2} new T2 lesions")

            if gad_enh is not None and gad_enh >= 1:
                escalate = True
                escalation_reasons.append(f"{gad_enh} enhancing lesion(s)")

            if escalate:
                findings.append(
                    f"DMT ESCALATION TRIGGER: {'; '.join(escalation_reasons)}"
                )
                # Recommend next tier
                dmt_cat_lower = (dmt_cat or "").lower()
                if dmt_cat_lower == "platform":
                    recommendations.append(
                        "Escalate to high-efficacy DMT: consider ocrelizumab, "
                        "natalizumab (if JCV-), or ofatumumab"
                    )
                elif dmt_cat_lower == "moderate_efficacy":
                    recommendations.append(
                        "Escalate to high-efficacy DMT: ocrelizumab, natalizumab, "
                        "or alemtuzumab"
                    )
                elif dmt_cat_lower == "high_efficacy":
                    recommendations.append(
                        "Already on high-efficacy DMT with breakthrough activity. "
                        "Consider: switch within class, add rituximab, "
                        "or evaluate for HSCT in eligible patients"
                    )
            else:
                findings.append("No DMT escalation trigger identified")
        else:
            recommendations.append(
                "No DMT documented. Initiate disease-modifying therapy -- "
                "early treatment improves long-term outcomes"
            )

        # -- PML risk (JCV) --
        if jcv_index is not None:
            findings.append(f"JCV antibody index: {jcv_index:.2f}")
            if jcv_index > 1.5:
                findings.append(
                    "HIGH PML RISK: JCV index >1.5 -- natalizumab is CONTRAINDICATED"
                )
                severity = _max_severity(severity, SeverityLevel.HIGH)
                recommendations.append(
                    "Do NOT use natalizumab. If currently on natalizumab, "
                    "switch to alternative high-efficacy DMT (ocrelizumab, ofatumumab)"
                )
            elif jcv_index > 0.9:
                findings.append(
                    "Moderate PML risk: JCV index 0.9-1.5 -- use natalizumab with caution"
                )
                recommendations.append(
                    "If on natalizumab: consider extended-interval dosing (q6 weeks) "
                    "or switch to anti-CD20"
                )
            elif jcv_index <= 0.9:
                findings.append("Low PML risk: JCV index <= 0.9")

        # -- PPMS specific --
        if phenotype == "ppms":
            recommendations.append(
                "Ocrelizumab is the only FDA-approved DMT for PPMS -- "
                "consider if not already on therapy"
            )
            guidelines.append(
                "Montalban et al. 2017 (ORATORIO): Ocrelizumab in PPMS"
            )

        # -- Standard monitoring --
        recommendations.extend([
            "MRI brain (with/without Gd) annually or sooner if clinical change",
            "Annual EDSS assessment and neurological examination",
            "Monitor for DMT-specific adverse effects",
        ])

        guidelines.extend([
            "Thompson et al. 2018: McDonald 2017 criteria for MS diagnosis",
            "AAN 2018: Practice guideline on DMT for MS",
            "Kappos et al. 2021: NEDA-3 and NEDA-4 as treatment targets",
        ])

        return WorkflowResult(
            workflow_type=self.workflow_type,
            domain=self.domain,
            findings=findings,
            scale_results=scale_results,
            recommendations=recommendations,
            guideline_references=guidelines,
            severity=severity,
        )


# ===================================================================
# WORKFLOW 6 -- Parkinson's Assessment
# ===================================================================


class ParkinsonsAssessmentWorkflow(BaseNeuroWorkflow):
    """Parkinson's disease assessment: MDS diagnostic criteria, UPDRS Part III,
    Hoehn & Yahr, tremor-dominant vs PIGD, DaT-SPECT, DBS candidacy
    (CAPSIT-PD), atypical red flags.

    Inputs
    ------
    updrs_scores : dict
        Per-item UPDRS Part III scores.
    hoehn_yahr : float
        Hoehn & Yahr stage (1-5).
    symptom_duration_years : float
    age_at_onset : int
    age : int
    dat_scan_result : str | None
        'abnormal', 'normal', or None.
    levodopa_response : str | None
        'excellent' (>30% UPDRS improvement), 'moderate', 'poor', 'not_tested'.
    motor_fluctuations : bool
        Wearing off, dyskinesias.
    dyskinesias : bool
    cognitive_status : str | None
        'normal', 'mci', 'dementia'.
    moca_domain_scores : dict | None
    autonomic_symptoms : list[str]
        E.g. 'orthostatic_hypotension', 'urinary', 'constipation'.
    red_flags : list[str]
        Atypical features: 'early_falls', 'rapid_progression',
        'cerebellar_signs', 'vertical_gaze_palsy', 'early_dementia',
        'early_autonomic_failure', 'poor_levodopa_response',
        'symmetric_onset', 'inspiratory_stridor', 'anterocollis'.
    """

    workflow_type = NeuroWorkflowType.PARKINSONS_ASSESSMENT
    domain = NeuroDomain.MOVEMENT_DISORDERS

    # MDS criteria: absolute exclusion criteria for PD
    _EXCLUSION_RED_FLAGS = {
        "cerebellar_signs": "Cerebellar ataxia -- consider MSA-C",
        "vertical_gaze_palsy": "Supranuclear vertical gaze palsy -- consider PSP",
        "early_dementia": "Dementia preceding or within 1 year of motor onset -- consider DLB",
        "cortical_sensory_loss": "Cortical sensory loss -- consider CBS",
        "early_autonomic_failure": "Severe autonomic failure within 5 years -- consider MSA",
        "inspiratory_stridor": "Inspiratory stridor -- consider MSA",
        "poor_levodopa_response": "Absence of levodopa response at adequate doses",
    }

    # Supportive criteria (need >= 2 for clinically established PD)
    _SUPPORTIVE_CRITERIA = {
        "rest_tremor", "excellent_levodopa_response", "levodopa_dyskinesias",
        "olfactory_loss", "cardiac_sympathetic_denervation",
    }

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        guidelines: List[str] = []
        scale_results: List[ScaleResult] = []
        severity = SeverityLevel.MODERATE
        cross_triggers: List[str] = []

        # These were printed as findings with their invented defaults: a request with none of
        # them reported "Age: 0, onset age: 0, duration: 0.0 years" and then used duration 0
        # as grounds for "Disease duration <5 years".
        age = inputs.get("age")
        age_onset = inputs.get("age_at_onset")
        duration = inputs.get(
            "symptom_duration_years", inputs.get("disease_duration_years")
        )
        findings.append(
            "Age: {}, onset age: {}, duration: {}".format(
                age if age is not None else "not recorded",
                age_onset if age_onset is not None else "not recorded",
                f"{duration:.1f} years" if duration is not None else "not recorded",
            )
        )

        # -- UPDRS Part III --
        updrs_scores = inputs.get("updrs_scores")
        if updrs_scores:
            updrs_result = UPDRSCalculator.calculate(updrs_scores)
            scale_results.append(updrs_result)
            findings.append(
                f"UPDRS Part III: {int(updrs_result.score)}/132 "
                f"({updrs_result.severity_category})"
            )

            # Tremor-dominant vs PIGD subtyping
            tremor_items = [
                "3.15a_postural_tremor_r", "3.15b_postural_tremor_l",
                "3.16a_kinetic_tremor_r", "3.16b_kinetic_tremor_l",
                "3.17a_rest_tremor_rue", "3.17b_rest_tremor_lue",
                "3.17c_rest_tremor_rle", "3.17d_rest_tremor_lle",
                "3.17e_rest_tremor_jaw", "3.18_constancy_of_rest_tremor",
            ]
            pigd_items = [
                "3.10_gait", "3.11_freezing_of_gait",
                "3.12_postural_stability", "3.9_arising_from_chair",
            ]

            tremor_score = sum(updrs_scores.get(i, 0) for i in tremor_items)
            pigd_score = sum(updrs_scores.get(i, 0) for i in pigd_items)

            tremor_mean = tremor_score / len(tremor_items) if tremor_items else 0
            pigd_mean = pigd_score / len(pigd_items) if pigd_items else 0

            if tremor_mean > 0 or pigd_mean > 0:
                ratio = tremor_mean / pigd_mean if pigd_mean > 0 else float('inf')
                if ratio >= 1.5:
                    subtype = ParkinsonsSubtype.TREMOR_DOMINANT
                elif ratio <= 0.67:
                    subtype = ParkinsonsSubtype.PIGD
                else:
                    subtype = ParkinsonsSubtype.INDETERMINATE
                findings.append(
                    f"Motor subtype: {subtype.value} "
                    f"(tremor mean {tremor_mean:.2f}, PIGD mean {pigd_mean:.2f})"
                )
                if subtype == ParkinsonsSubtype.PIGD:
                    findings.append(
                        "PIGD subtype carries higher risk of cognitive decline, "
                        "falls, and disability progression"
                    )

        # -- Hoehn & Yahr --
        hy = inputs.get("hoehn_yahr")
        if hy is not None:
            hy_result = HoehnYahrCalculator.calculate(hy)
            scale_results.append(hy_result)
            findings.append(f"Hoehn & Yahr: {hy_result.severity_category}")

        # -- DaT-SPECT --
        dat = inputs.get("dat_scan_result")
        if dat:
            if dat == "abnormal":
                findings.append(
                    "DaT-SPECT: ABNORMAL (reduced dopamine transporter uptake) -- "
                    "supports dopaminergic deficit"
                )
            elif dat == "normal":
                findings.append(
                    "DaT-SPECT: NORMAL -- diagnosis of degenerative parkinsonism "
                    "less likely; consider essential tremor, drug-induced, "
                    "functional parkinsonism"
                )
                severity = SeverityLevel.LOW

        # -- Levodopa response --
        levo_resp = inputs.get("levodopa_response")
        if levo_resp:
            findings.append(f"Levodopa response: {levo_resp}")
            if levo_resp == "excellent":
                findings.append("Excellent levodopa response (>30% UPDRS improvement) -- supportive of PD")
            elif levo_resp == "poor":
                findings.append(
                    "Poor levodopa response -- consider atypical parkinsonism "
                    "(MSA, PSP, CBS)"
                )
                severity = _max_severity(severity, SeverityLevel.HIGH)

        # -- Red flags for atypical parkinsonism --
        red_flags = inputs.get("red_flags", [])
        atypical_findings: List[str] = []
        for flag in red_flags:
            desc = self._EXCLUSION_RED_FLAGS.get(flag)
            if desc:
                atypical_findings.append(desc)

        if atypical_findings:
            findings.append(
                f"ATYPICAL RED FLAGS ({len(atypical_findings)}): "
                + "; ".join(atypical_findings)
            )
            severity = _max_severity(severity, SeverityLevel.HIGH)
            recommendations.append(
                "Atypical features identified -- consider alternative diagnoses: "
                "MSA, PSP, CBS, DLB"
            )
            recommendations.append(
                "Consider brain MRI to evaluate for atypical parkinsonism patterns "
                "(hot cross bun, hummingbird sign, cortical atrophy)"
            )

        # -- DBS candidacy (CAPSIT-PD) --
        motor_fluct = inputs.get("motor_fluctuations", False)
        dyskinesias = inputs.get("dyskinesias", False)
        # N13(a): this defaulted to "normal", so a patient whose cognition was never assessed
        # yielded "DBS CANDIDATE: ... no cognitive contraindication" -- an assertion about a
        # domain nobody had examined. A MoCA, when supplied, was computed and printed but
        # never read by the DBS decision.
        cog_status = inputs.get("cognitive_status")

        # MoCA for cognitive screening
        moca_domains = inputs.get("moca_domain_scores")
        moca_result = None
        if moca_domains:
            moca_result = MoCACalculator.calculate(
                moca_domains, inputs.get("education_years")
            )
            scale_results.append(moca_result)
            findings.append(f"MoCA: {int(moca_result.score)}/30 ({moca_result.severity_category})")

        dbs_candidate = True
        dbs_reasons: List[str] = []

        if not motor_fluct and not dyskinesias:
            dbs_candidate = False
            dbs_reasons.append("No motor fluctuations or dyskinesias")
        if cog_status == "dementia":
            dbs_candidate = False
            dbs_reasons.append("Cognitive impairment (dementia) is a contraindication")
        elif cog_status is None:
            if moca_result is not None and moca_result.complete:
                # Use the MoCA that was already computed rather than asserting a status.
                if moca_result.score < 26:
                    dbs_candidate = False
                    dbs_reasons.append(
                        f"MoCA {int(moca_result.score)}/30 is below the cognitive threshold "
                        f"this agent applies to DBS candidacy (MoCA >= 26) -- formal "
                        f"neuropsychological assessment required"
                    )
            else:
                dbs_candidate = False
                dbs_reasons.append(
                    "Cognitive status NOT ASSESSED -- DBS candidacy cannot be cleared on an "
                    "unexamined cognitive domain; obtain neuropsychological testing"
                )
        if levo_resp == "poor":
            dbs_candidate = False
            dbs_reasons.append("Poor levodopa response")
        if red_flags:
            dbs_candidate = False
            dbs_reasons.append("Atypical features present")
        if duration is None:
            dbs_reasons.append(
                "Disease duration not recorded -- the diagnosis should be well-established "
                "before DBS"
            )
        elif duration < 5:
            # NOTE: this workflow uses < 5 years while api/routes/neuro_clinical.py uses
            # >= 4. Neither is cited; both left unchanged, divergence recorded.
            dbs_reasons.append("Disease duration <5 years -- diagnosis should be well-established")
        if age is not None and age > 75:
            dbs_reasons.append("Age >75 -- relative contraindication; consider carefully")

        if dbs_candidate and updrs_scores:
            updrs_total = sum(
                max(0, min(updrs_scores.get(i, 0), 4))
                for i in UPDRSCalculator.ITEMS
            )
            # OPEN (uncited, unchanged): the >= 59 gate is cited nowhere and is not a
            # CAPSIT-PD criterion; and if `updrs_scores` is absent this whole block prints
            # nothing at all, which is why the "not currently indicated" branch below now
            # also fires when the DBS assessment could not be completed.
            if updrs_total >= 59:
                findings.append(
                    "DBS CANDIDATE: motor fluctuations/dyskinesias present, "
                    "levodopa-responsive, and no cognitive contraindication on the "
                    "assessment supplied"
                )
                recommendations.append(
                    "Refer for DBS evaluation (CAPSIT-PD protocol): "
                    "neuropsych testing, ON/OFF levodopa challenge, brain MRI"
                )
                recommendations.append(
                    "DBS targets: STN (bilateral) preferred for medication reduction, "
                    "GPi for dyskinesia control"
                )
        elif dbs_reasons:
            findings.append(
                f"DBS not currently indicated: {'; '.join(dbs_reasons)}"
            )

        # -- General recommendations --
        if not levo_resp or levo_resp == "not_tested":
            recommendations.append(
                "Initiate levodopa trial (carbidopa/levodopa 25/100 TID) "
                "to assess response"
            )
        if motor_fluct:
            recommendations.extend([
                "Optimise levodopa dosing: consider controlled-release, "
                "add COMT inhibitor (entacapone), or MAO-B inhibitor (rasagiline)",
                "Consider LCIG (Duopa) for refractory motor fluctuations",
            ])

        guidelines.extend([
            "MDS 2015: Diagnostic Criteria for Parkinson's Disease",
            "AAN 2006: Practice parameter on neuroprotective strategies in PD",
            "CAPSIT-PD: Criteria for DBS candidate selection",
            "MDS Evidence-Based Review: Treatments for motor symptoms of PD",
        ])

        return WorkflowResult(
            workflow_type=self.workflow_type,
            domain=self.domain,
            findings=findings,
            scale_results=scale_results,
            recommendations=recommendations,
            guideline_references=guidelines,
            severity=severity,
            cross_modal_triggers=cross_triggers,
        )


# ===================================================================
# WORKFLOW 7 -- Headache Classification
# ===================================================================


class HeadacheClassificationWorkflow(BaseNeuroWorkflow):
    """Headache classification: ICHD-3, SNOOP red flags, acute treatment,
    preventive candidacy, CGRP therapy eligibility, HIT-6/MIDAS scoring.

    Inputs
    ------
    headache_features : dict
        Keys: 'location' (str), 'quality' (str: throbbing/pressing/stabbing),
        'intensity' (str: mild/moderate/severe),
        'unilateral' (bool), 'duration_hours' (float),
        'frequency_per_month' (int), 'nausea' (bool), 'vomiting' (bool),
        'photophobia' (bool), 'phonophobia' (bool),
        'aura_present' (bool), 'aura_type' (str | None),
        'autonomic_features' (list[str]),  # lacrimation, rhinorrhea, ptosis
        'worsened_by_activity' (bool).
    hit6_responses : list[str] | None
        Six HIT-6 categorical responses.
    midas_score : int | None
        MIDAS disability score (0-270+).
    red_flags : list[str]
        SNOOP red flags present.
    age : int
    medication_overuse : bool
        Acute medication use >= 10-15 days/month.
    current_acute_meds : list[str]
    current_preventives : list[str]
    num_failed_preventives : int
    """

    workflow_type = NeuroWorkflowType.HEADACHE_CLASSIFICATION
    domain = NeuroDomain.HEADACHE

    # SNOOP4 red flag descriptors
    _SNOOP_FLAGS: Dict[str, str] = {
        "systemic_symptoms": "Systemic symptoms (fever, weight loss, cancer hx) -- S",
        "neurologic_signs": "Neurologic signs or symptoms -- N",
        "onset_sudden": "Sudden/thunderclap onset -- O (SAH until proven otherwise)",
        "onset_after_40": "New headache onset after age 40 -- O",
        "pattern_change": "Progressive pattern or change from prior headaches -- P",
        "positional": "Positional headache (worse upright or supine) -- P",
        "papilledema": "Papilledema -- P",
        "precipitated_by_valsalva": "Precipitated by Valsalva or exertion -- P",
        "pregnancy": "Pregnancy or postpartum -- P",
        "painful_eye": "Painful eye with autonomic features -- consider TAC",
        "posttraumatic": "Post-traumatic onset -- P",
        "immunosuppressed": "Immunosuppressed patient -- consider CNS infection",
    }

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        guidelines: List[str] = []
        scale_results: List[ScaleResult] = []
        severity = SeverityLevel.MODERATE

        hf = inputs.get("headache_features", {})
        inputs.get("age", 0)

        # -- SNOOP red flags --
        red_flags = inputs.get("red_flags", [])
        if red_flags:
            severity = SeverityLevel.HIGH
            findings.append(f"RED FLAGS PRESENT ({len(red_flags)}):")
            for flag in red_flags:
                desc = self._SNOOP_FLAGS.get(flag, flag)
                findings.append(f"  - {desc}")

            if "onset_sudden" in red_flags:
                severity = SeverityLevel.CRITICAL
                recommendations.extend([
                    "EMERGENT: Rule out subarachnoid haemorrhage -- CT head STAT, "
                    "if negative consider LP or CTA",
                    "Thunderclap headache protocol: CT -> CTA -> LP if needed",
                ])

            recommendations.append(
                "Secondary headache workup indicated: MRI brain with/without contrast, "
                "consider MRA/MRV"
            )

            if "papilledema" in red_flags:
                recommendations.append(
                    "Urgent ophthalmology consultation; consider LP with opening "
                    "pressure for idiopathic intracranial hypertension"
                )

        # -- ICHD-3 classification --
        unilateral = hf.get("unilateral", False)
        throbbing = "throb" in str(hf.get("quality", "")).lower()
        nausea = hf.get("nausea", False)
        photophobia = hf.get("photophobia", False)
        phonophobia = hf.get("phonophobia", False)
        aura = hf.get("aura_present", False)
        duration = hf.get("duration_hours", 0)
        frequency = hf.get("frequency_per_month", 0)
        worsened_activity = hf.get("worsened_by_activity", False)
        autonomic = hf.get("autonomic_features", [])
        quality = hf.get("quality", "")
        intensity = hf.get("intensity", "")

        headache_dx: Optional[HeadacheType] = None
        med_overuse = inputs.get("medication_overuse", False)

        # Migraine criteria (ICHD-3)
        migraine_features = sum([
            unilateral, throbbing, intensity in ("moderate", "severe"),
            worsened_activity,
        ])
        migraine_assoc = sum([nausea, photophobia and phonophobia])

        if migraine_features >= 2 and migraine_assoc >= 1 and 4 <= duration <= 72:
            if frequency >= 15 and duration >= 4:
                headache_dx = HeadacheType.CHRONIC_MIGRAINE
                # N16: the finding asserted ">= 15 headache days/month, >= 8 with migraine
                # features, >= 3 months" while only the >= 15 was ever tested -- two criteria
                # were stated as verified and never checked. They are reported as outstanding
                # unless the caller supplies them.
                migraine_days = hf.get("migraine_days_per_month")
                months = hf.get("duration_months")
                verified = [f"{frequency} headache days/month (>= 15)"]
                outstanding = []
                if migraine_days is None:
                    outstanding.append(">= 8 days/month with migraine features (not supplied)")
                elif migraine_days >= 8:
                    verified.append(f"{migraine_days} migraine days/month (>= 8)")
                else:
                    outstanding.append(
                        f"only {migraine_days} migraine days/month (ICHD-3 requires >= 8)"
                    )
                if months is None:
                    outstanding.append("duration > 3 months (not supplied)")
                elif months > 3:
                    verified.append(f"{months} months (> 3)")
                else:
                    outstanding.append(
                        f"only {months} months (ICHD-3 requires > 3)"
                    )
                findings.append(
                    "ICHD-3 chronic migraine -- verified: " + "; ".join(verified)
                    + (
                        ". NOT VERIFIED: " + "; ".join(outstanding)
                        if outstanding else ". All ICHD-3 1.3 criteria verified."
                    )
                )
            elif aura:
                headache_dx = HeadacheType.MIGRAINE_WITH_AURA
                findings.append("ICHD-3: Migraine with aura")
                aura_type = hf.get("aura_type", "visual")
                findings.append(f"  Aura type: {aura_type}")
            else:
                headache_dx = HeadacheType.MIGRAINE_WITHOUT_AURA
                findings.append("ICHD-3: Migraine without aura")

        # Tension-type headache
        elif ("press" in str(quality).lower() or "tight" in str(quality).lower()):
            if not unilateral and not nausea and not throbbing:
                if frequency >= 15:
                    headache_dx = HeadacheType.TENSION_TYPE_CHRONIC
                    findings.append("ICHD-3: Chronic tension-type headache")
                else:
                    headache_dx = HeadacheType.TENSION_TYPE_EPISODIC
                    findings.append("ICHD-3: Episodic tension-type headache")

        # Cluster / TAC
        elif unilateral and autonomic and intensity == "severe":
            # N16: `duration <= 3` had no lower bound, so the `duration = 0` default -- i.e.
            # a duration nobody recorded -- satisfied a criterion whose own printed text says
            # "duration 15-180 min". ICHD-3 3.1 specifies attacks of 15-180 minutes, which is
            # 0.25-3 hours in this field's units.
            if 0.25 <= duration <= 3 and len(autonomic) >= 1:
                headache_dx = HeadacheType.CLUSTER
                findings.append(
                    "ICHD-3: Cluster headache (severe unilateral with "
                    "autonomic features, duration 15-180 min)"
                )
                severity = SeverityLevel.HIGH

        # Medication overuse
        if med_overuse:
            findings.append(
                "Medication overuse headache likely: acute medication use >= 10-15 days/month"
            )
            if not headache_dx:
                headache_dx = HeadacheType.MEDICATION_OVERUSE
            recommendations.append(
                "Discontinue overused acute medication with bridging strategy "
                "(naproxen, nerve blocks, or short steroid taper)"
            )

        if not headache_dx and not red_flags:
            headache_dx = HeadacheType.MIGRAINE_WITHOUT_AURA  # most common
            findings.append(
                "Classification uncertain -- defaulting to migraine workup "
                "(most common primary headache)"
            )

        # -- HIT-6 --
        # HIT6Calculator now refuses a short or unrecognised response set rather than padding
        # it with "never" (N9), so the refusal is reported as a finding instead of producing a
        # best-case score.
        hit6_resp = inputs.get("hit6_responses")
        hit6_result = None
        if hit6_resp:
            try:
                HIT6Calculator.calculate(hit6_resp)
            except ScaleDataError as exc:
                findings.append(f"HIT-6 NOT SCORED: {exc}")
                hit6_resp = None
        if hit6_resp:
            hit6_result = HIT6Calculator.calculate(hit6_resp)
            scale_results.append(hit6_result)
            findings.append(
                f"HIT-6: {int(hit6_result.score)} ({hit6_result.severity_category})"
            )

        # -- MIDAS --
        midas = inputs.get("midas_score")
        if midas is not None:
            findings.append(f"MIDAS disability score: {midas}")
            if midas <= 5:
                findings.append("MIDAS Grade I: Little or no disability")
            elif midas <= 10:
                findings.append("MIDAS Grade II: Mild disability")
            elif midas <= 20:
                findings.append("MIDAS Grade III: Moderate disability")
            else:
                findings.append("MIDAS Grade IV: Severe disability")

        # -- Acute treatment plan --
        if headache_dx in (
            HeadacheType.MIGRAINE_WITHOUT_AURA,
            HeadacheType.MIGRAINE_WITH_AURA,
            HeadacheType.CHRONIC_MIGRAINE,
        ):
            current_acute = inputs.get("current_acute_meds", [])
            if not current_acute:
                recommendations.extend([
                    "Acute migraine treatment: triptan (sumatriptan 50-100 mg PO, "
                    "or rizatriptan 10 mg) for moderate-severe attacks",
                    "NSAIDs (ibuprofen 400-600 mg, naproxen 500 mg) for mild attacks",
                    "Anti-emetic (metoclopramide 10 mg) if nausea/vomiting prominent",
                ])
            recommendations.append(
                "Limit acute medication use to <10 days/month to prevent MOH"
            )

        elif headache_dx == HeadacheType.CLUSTER:
            recommendations.extend([
                "Acute cluster treatment: sumatriptan 6 mg SC or "
                "high-flow oxygen (12-15 L/min x 15 min via non-rebreather)",
                "Transitional prevention: prednisone taper or greater occipital "
                "nerve block",
                "Maintenance prevention: verapamil (start 240 mg/day, "
                "titrate up with ECG monitoring)",
            ])

        # -- Preventive candidacy --
        inputs.get("current_preventives", [])
        num_failed = inputs.get("num_failed_preventives", 0)

        preventive_indicated = (
            frequency >= 4
            or (midas is not None and midas > 10)
            or (hit6_result is not None and hit6_result.score >= 56)
            or headache_dx == HeadacheType.CHRONIC_MIGRAINE
        )

        if preventive_indicated:
            findings.append("Preventive therapy INDICATED")
            if num_failed == 0:
                recommendations.append(
                    "First-line preventives: topiramate (25-100 mg), propranolol "
                    "(40-240 mg), amitriptyline (10-75 mg), or valproate (500-1500 mg)"
                )
            elif num_failed == 1:
                recommendations.append(
                    "Try alternative first-line preventive from different class"
                )

            # CGRP therapy eligibility (after >= 2 preventive failures)
            if num_failed >= 2:
                findings.append(
                    "CGRP THERAPY ELIGIBLE: >= 2 prior preventive failures"
                )
                recommendations.extend([
                    "CGRP monoclonal antibody: erenumab (70-140 mg SC monthly), "
                    "fremanezumab (225 mg monthly or 675 mg quarterly), "
                    "or galcanezumab (240 mg loading, then 120 mg monthly)",
                    "Alternative: CGRP receptor antagonist (gepant class) -- "
                    "atogepant (60 mg daily) or rimegepant (75 mg every other day)",
                ])
            if num_failed >= 3:
                recommendations.append(
                    "Consider onabotulinum toxin A (Botox) 155-195 units "
                    "for chronic migraine (>= 15 days/month)"
                )
                recommendations.append(
                    "Refer to headache specialist or comprehensive headache center"
                )

        guidelines.extend([
            "ICHD-3 2018: International Classification of Headache Disorders (3rd edition)",
            "AHS 2021: Consensus Statement on CGRP-targeting therapies",
            "AAN 2012: Evidence-based guideline on episodic migraine prevention",
        ])

        return WorkflowResult(
            workflow_type=self.workflow_type,
            domain=self.domain,
            findings=findings,
            scale_results=scale_results,
            recommendations=recommendations,
            guideline_references=guidelines,
            severity=severity,
        )


# ===================================================================
# WORKFLOW 8 -- Neuromuscular Evaluation
# ===================================================================


class NeuromuscularEvaluationWorkflow(BaseNeuroWorkflow):
    """Neuromuscular evaluation: EMG/NCS pattern classification (axonal vs
    demyelinating vs NMJ vs myopathic), differential diagnosis, antibody
    panel recommendations, genetic testing.

    Inputs
    ------
    emg_ncs : dict
        Keys: 'motor_cv' (list[float], m/s), 'sensory_cv' (list[float]),
        'distal_latency_prolonged' (bool), 'f_wave_prolonged' (bool),
        'conduction_block' (bool), 'temporal_dispersion' (bool),
        'fibrillations' (bool), 'fasciculations' (bool),
        'positive_sharp_waves' (bool), 'myopathic_units' (bool),
        'decremental_response' (bool), 'incremental_response' (bool),
        'motor_amplitude_reduced' (bool), 'sensory_amplitude_reduced' (bool),
        'jitter_increased' (bool).
    weakness_pattern : str
        'proximal', 'distal', 'bulbar', 'respiratory', 'generalized',
        'asymmetric_limb', 'ocular'.
    distribution : str
        'symmetric', 'asymmetric'.
    progression : str
        'acute' (<4 weeks), 'subacute' (4-8 weeks), 'chronic' (>8 weeks).
    sensory_involvement : bool
    ck_level : float | None
    age : int
    family_history_neuromuscular : bool
    existing_antibodies : dict
        Already-obtained antibody results.
    alsfrs_scores : dict | None
        ALSFRS-R item scores if ALS suspected.
    months_since_onset : float | None
    """

    workflow_type = NeuroWorkflowType.NEUROMUSCULAR_EVALUATION
    domain = NeuroDomain.NEUROMUSCULAR

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        guidelines: List[str] = []
        scale_results: List[ScaleResult] = []
        severity = SeverityLevel.MODERATE
        cross_triggers: List[str] = []

        emg = inputs.get("emg_ncs", {})
        weakness = inputs.get("weakness_pattern", "generalized")
        distribution = inputs.get("distribution", "symmetric")
        # N14: this defaulted to "chronic", and "chronic" is printed as a finding and routes
        # the demyelinating differential away from Guillain-Barre -- which carries the only
        # urgent actions in this file: "URGENT: Monitor respiratory function (FVC q4h) --
        # intubate if FVC < 20 mL/kg", IVIg/plasmapheresis, and severity CRITICAL. A
        # demyelinating EMG with `progression` absent reported "Progression: chronic" (a value
        # nobody supplied), a CIDP/MMN/anti-MAG/CMT1 differential, severity MODERATE, and an
        # antibody panel plus a CMT gene panel -- with no airway monitoring for a
        # demyelinating polyneuropathy of unknown tempo. GBS respiratory failure is the
        # non-recoverable error in this domain.
        progression = inputs.get("progression")
        sensory = inputs.get("sensory_involvement", False)
        ck = inputs.get("ck_level")
        age = inputs.get("age")
        fam_hx = inputs.get("family_history_neuromuscular", False)

        findings.append(f"Weakness pattern: {weakness}, distribution: {distribution}")
        if progression is None:
            findings.append(
                "Progression: NOT RECORDED -- the tempo separates an acute (GBS) from a "
                "chronic (CIDP) presentation and determines whether respiratory monitoring "
                "is required"
            )
        else:
            findings.append(f"Progression: {progression}")
        if sensory:
            findings.append("Sensory involvement present")
        if ck is not None:
            findings.append(f"CK level: {ck:.0f} U/L")
            if ck > 1000:
                findings.append("Markedly elevated CK -- suggests myopathic process")
            elif ck > 200:
                findings.append("Mildly elevated CK")

        # -- EMG/NCS pattern classification --
        pattern: Optional[NMJPattern] = None
        differentials: List[str] = []

        motor_cv = emg.get("motor_cv", [])
        emg.get("sensory_cv", [])
        fibs = emg.get("fibrillations", False)
        fascs = emg.get("fasciculations", False)
        psw = emg.get("positive_sharp_waves", False)
        myopathic_units = emg.get("myopathic_units", False)
        cond_block = emg.get("conduction_block", False)
        temp_disp = emg.get("temporal_dispersion", False)
        decremental = emg.get("decremental_response", False)
        incremental = emg.get("incremental_response", False)
        motor_amp_low = emg.get("motor_amplitude_reduced", False)
        sensory_amp_low = emg.get("sensory_amplitude_reduced", False)
        dl_prolonged = emg.get("distal_latency_prolonged", False)
        f_wave_prolonged = emg.get("f_wave_prolonged", False)
        jitter = emg.get("jitter_increased", False)

        # Check for demyelinating pattern
        slow_motor = any(cv < 35 for cv in motor_cv) if motor_cv else False
        is_demyelinating = (
            slow_motor or cond_block or temp_disp
            or dl_prolonged or f_wave_prolonged
        )

        # Check for NMJ pattern
        is_nmj = decremental or incremental or jitter

        if is_nmj:
            if decremental:
                pattern = NMJPattern.NMJ_POSTSYNAPTIC
                findings.append(
                    "EMG pattern: NEUROMUSCULAR JUNCTION -- POSTSYNAPTIC "
                    "(decremental response on repetitive nerve stimulation)"
                )
                differentials = [
                    "Myasthenia gravis (autoimmune)",
                    "Congenital myasthenic syndrome",
                ]
            elif incremental:
                pattern = NMJPattern.NMJ_PRESYNAPTIC
                findings.append(
                    "EMG pattern: NEUROMUSCULAR JUNCTION -- PRESYNAPTIC "
                    "(incremental response on high-frequency RNS)"
                )
                differentials = [
                    "Lambert-Eaton myasthenic syndrome (LEMS)",
                    "Botulism",
                ]
            if jitter:
                findings.append("Increased jitter on single-fibre EMG confirms NMJ disorder")

        elif myopathic_units:
            pattern = NMJPattern.MYOPATHIC
            findings.append(
                "EMG pattern: MYOPATHIC (short-duration, low-amplitude, "
                "polyphasic motor units with early recruitment)"
            )
            differentials = [
                "Inflammatory myopathy (dermatomyositis, polymyositis, IBM)",
                "Muscular dystrophy",
                "Metabolic myopathy (acid maltase deficiency, mitochondrial)",
                "Toxic/drug-induced myopathy (statins, steroids)",
            ]
            if ck and ck > 5000:
                differentials.insert(0, "Necrotizing autoimmune myopathy")

        elif is_demyelinating:
            pattern = NMJPattern.DEMYELINATING
            findings.append(
                "EMG pattern: DEMYELINATING (slow conduction velocities, "
                "conduction block, temporal dispersion, prolonged distal latencies)"
            )
            if progression == "acute":
                differentials = [
                    "Guillain-Barre syndrome (AIDP)",
                    "AMAN (axonal GBS variant)",
                ]
                severity = SeverityLevel.CRITICAL
                recommendations.extend([
                    "URGENT: Monitor respiratory function (FVC q4h) -- "
                    "intubate if FVC < 20 mL/kg or declining rapidly",
                    "Start IVIg (0.4 g/kg/day x 5 days) or plasmapheresis (5 exchanges)",
                ])
            elif progression is None:
                # Tempo not established. GBS cannot be excluded, so it stays on the
                # differential and keeps its airway monitoring.
                differentials = [
                    "Guillain-Barre syndrome (AIDP) -- CANNOT BE EXCLUDED: tempo not recorded",
                    "CIDP (chronic inflammatory demyelinating polyneuropathy)",
                    "MMN (multifocal motor neuropathy with conduction block)",
                    "Anti-MAG neuropathy",
                    "Charcot-Marie-Tooth (hereditary demyelinating, CMT1)",
                ]
                severity = SeverityLevel.CRITICAL
                recommendations.extend([
                    "Establish the tempo of progression (symptom onset date, nadir, whether "
                    "still worsening) -- it decides between GBS and CIDP",
                    "URGENT until GBS is excluded: monitor respiratory function (FVC q4h, "
                    "negative inspiratory force, single-breath count) -- intubate if "
                    "FVC < 20 mL/kg or declining rapidly; monitor for autonomic instability",
                ])
            else:
                differentials = [
                    "CIDP (chronic inflammatory demyelinating polyneuropathy)",
                    "MMN (multifocal motor neuropathy with conduction block)",
                    "Anti-MAG neuropathy",
                    "Charcot-Marie-Tooth (hereditary demyelinating, CMT1)",
                ]

        elif (fibs or fascs or psw) and motor_amp_low and not sensory_amp_low:
            pattern = NMJPattern.AXONAL_MOTOR
            findings.append(
                "EMG pattern: AXONAL MOTOR (fibrillations/fasciculations, "
                "reduced motor amplitudes, preserved sensory)"
            )
            differentials = [
                "Amyotrophic lateral sclerosis (ALS)",
                "Progressive muscular atrophy",
                "Spinal muscular atrophy",
                "Multifocal motor neuropathy",
            ]
            if fascs and fibs and weakness in ("asymmetric_limb", "generalized", "bulbar"):
                findings.append(
                    "Pattern highly suggestive of ALS (Awaji criteria): "
                    "fasciculations + fibrillations + UMN signs"
                )
                severity = SeverityLevel.HIGH

        elif motor_amp_low and sensory_amp_low:
            pattern = NMJPattern.AXONAL_SENSORIMOTOR
            findings.append(
                "EMG pattern: AXONAL SENSORIMOTOR (reduced motor and sensory amplitudes)"
            )
            differentials = [
                "Diabetic polyneuropathy",
                "Toxic neuropathy (alcohol, chemotherapy)",
                "Vasculitic neuropathy",
                "Paraneoplastic neuropathy",
                "B12 deficiency",
                "Hereditary (CMT2)",
            ]

        elif sensory_amp_low and not motor_amp_low:
            pattern = NMJPattern.AXONAL_SENSORY
            findings.append(
                "EMG pattern: AXONAL SENSORY (isolated sensory amplitude reduction)"
            )
            differentials = [
                "Small fibre neuropathy",
                "Sensory ganglionopathy (paraneoplastic, Sjogren's)",
                "Vitamin B12/B6 deficiency",
                "Hereditary sensory neuropathy",
            ]

        if differentials:
            findings.append(f"Differential diagnosis: {', '.join(differentials)}")

        # -- ALSFRS-R if ALS suspected --
        alsfrs_scores = inputs.get("alsfrs_scores")
        if alsfrs_scores:
            months = inputs.get("months_since_onset")
            alsfrs_result = ALSFRSCalculator.calculate(alsfrs_scores, months)
            scale_results.append(alsfrs_result)
            findings.append(
                f"ALSFRS-R: {int(alsfrs_result.score)}/48 "
                f"({alsfrs_result.severity_category})"
            )

        # -- Antibody panel recommendations --
        existing_abs = inputs.get("existing_antibodies", {})
        recommended_abs: List[str] = []

        if pattern == NMJPattern.NMJ_POSTSYNAPTIC:
            recommended_abs = [
                "AChR binding antibodies",
                "AChR blocking antibodies",
                "AChR modulating antibodies",
                "Anti-MuSK antibodies",
                "Anti-LRP4 antibodies",
            ]
            recommendations.append(
                "If seropositive MG: start pyridostigmine (60 mg TID), "
                "consider immunosuppression (prednisone + steroid-sparing agent)"
            )
        elif pattern == NMJPattern.NMJ_PRESYNAPTIC:
            recommended_abs = [
                "P/Q-type VGCC antibodies (Lambert-Eaton)",
                "N-type VGCC antibodies",
                "SOX1 antibodies (paraneoplastic marker)",
            ]
            recommendations.append(
                "Screen for underlying malignancy (CT chest for SCLC) -- "
                "LEMS is paraneoplastic in ~60% of cases"
            )
        elif pattern == NMJPattern.MYOPATHIC:
            recommended_abs = [
                "ANA, anti-Jo-1 (myositis-specific)",
                "Anti-Mi-2, anti-MDA5, anti-TIF1-gamma",
                "Anti-SRP, anti-HMGCR (necrotizing myopathy)",
                "Anti-cN1A (inclusion body myositis)",
            ]
            recommendations.append(
                "Consider muscle biopsy for definitive diagnosis if inflammatory "
                "myopathy suspected"
            )
        elif pattern == NMJPattern.DEMYELINATING and progression != "acute":
            recommended_abs = [
                "Anti-MAG antibodies",
                "Anti-GM1 antibodies (MMN)",
                "SPEP/UPEP with immunofixation (paraprotein)",
                "Anti-neurofascin (NF155, NF186) for CIDP subtype",
            ]
        elif pattern in (NMJPattern.AXONAL_MOTOR,) and fascs:
            recommended_abs = [
                "Anti-GM1 antibodies (rule out MMN)",
            ]

        # Remove already-tested antibodies
        new_abs = [a for a in recommended_abs if a.split("(")[0].strip().lower()
                   not in {k.lower() for k in existing_abs}]
        if new_abs:
            recommendations.append(
                f"Recommended antibody panel: {', '.join(new_abs)}"
            )

        # -- Genetic testing recommendations --
        # A CMT gene panel is a chronic-neuropathy work-up. It used to be recommended for an
        # unknown-tempo patient because `progression` defaulted to "chronic", and `age < 40`
        # was true for the `age` default of 0.
        if fam_hx or (pattern == NMJPattern.DEMYELINATING and progression == "chronic"
                       and age is not None and age < 40):
            recommendations.append(
                "Consider genetic testing: CMT gene panel (PMP22, MPZ, GJB1, MFN2) "
                "for hereditary neuropathy"
            )
            if pattern == NMJPattern.MYOPATHIC and fam_hx:
                recommendations.append(
                    "Consider dystrophin gene testing and/or muscular dystrophy "
                    "gene panel"
                )

        if pattern == NMJPattern.AXONAL_MOTOR and age < 50 and fam_hx:
            recommendations.append(
                "Consider SMA genetic testing (SMN1/SMN2) if proximal weakness"
            )
            cross_triggers.append(
                _trigger_string("GENETIC", ["SMN1", "SMN2"],
                                "Hereditary motor neuropathy suspected")
            )

        # General
        if not recommendations:
            recommendations.append(
                "Comprehensive metabolic workup: glucose, HbA1c, B12, folate, "
                "TSH, SPEP, CMP, ESR, CRP"
            )

        guidelines.extend([
            "AAN 2009: Practice parameter for electrodiagnostic studies in NMD",
            "AAN 2021: Practice guideline on management of GBS",
            "Gilhus et al. 2019: Myasthenia gravis -- pathogenesis, diagnosis, treatment",
            "AAN 2012: Evidence-based guideline for CIDP",
        ])

        return WorkflowResult(
            workflow_type=self.workflow_type,
            domain=self.domain,
            findings=findings,
            scale_results=scale_results,
            recommendations=recommendations,
            guideline_references=guidelines,
            severity=severity,
            cross_modal_triggers=cross_triggers,
        )


# ===================================================================
# WORKFLOW ENGINE
# ===================================================================


class WorkflowEngine:
    """Central dispatcher that maps NeuroWorkflowType to the appropriate
    workflow implementation and handles query-based workflow detection."""

    _KEYWORD_MAP: Dict[str, NeuroWorkflowType] = {
        # Acute Stroke
        "stroke": NeuroWorkflowType.ACUTE_STROKE,
        "nihss": NeuroWorkflowType.ACUTE_STROKE,
        "tpa": NeuroWorkflowType.ACUTE_STROKE,
        "alteplase": NeuroWorkflowType.ACUTE_STROKE,
        "thrombectomy": NeuroWorkflowType.ACUTE_STROKE,
        "aspects": NeuroWorkflowType.ACUTE_STROKE,
        "large vessel occlusion": NeuroWorkflowType.ACUTE_STROKE,
        "lvo": NeuroWorkflowType.ACUTE_STROKE,
        "ischemic stroke": NeuroWorkflowType.ACUTE_STROKE,
        "hemorrhagic stroke": NeuroWorkflowType.ACUTE_STROKE,
        "cerebral hemorrhage": NeuroWorkflowType.ACUTE_STROKE,
        "subarachnoid": NeuroWorkflowType.ACUTE_STROKE,
        "tia": NeuroWorkflowType.ACUTE_STROKE,
        "dawn criteria": NeuroWorkflowType.ACUTE_STROKE,
        "defuse": NeuroWorkflowType.ACUTE_STROKE,
        # Dementia Evaluation
        "dementia": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "alzheimer": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "cognitive decline": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "cognitive impairment": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "moca": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "mmse": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "frontotemporal": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "lewy body": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "atn biomarker": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "amyloid": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "lecanemab": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "donanemab": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "apoe": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "hippocampal atrophy": NeuroWorkflowType.DEMENTIA_EVALUATION,
        "memory loss": NeuroWorkflowType.DEMENTIA_EVALUATION,
        # Epilepsy Focus
        "epilepsy": NeuroWorkflowType.EPILEPSY_FOCUS,
        "seizure": NeuroWorkflowType.EPILEPSY_FOCUS,
        "epileptic": NeuroWorkflowType.EPILEPSY_FOCUS,
        "eeg": NeuroWorkflowType.EPILEPSY_FOCUS,
        "antiseizure": NeuroWorkflowType.EPILEPSY_FOCUS,
        "anti-seizure": NeuroWorkflowType.EPILEPSY_FOCUS,
        "asm": NeuroWorkflowType.EPILEPSY_FOCUS,
        "aed": NeuroWorkflowType.EPILEPSY_FOCUS,
        "dravet": NeuroWorkflowType.EPILEPSY_FOCUS,
        "scn1a": NeuroWorkflowType.EPILEPSY_FOCUS,
        "temporal lobe epilepsy": NeuroWorkflowType.EPILEPSY_FOCUS,
        "drug-resistant epilepsy": NeuroWorkflowType.EPILEPSY_FOCUS,
        "epilepsy surgery": NeuroWorkflowType.EPILEPSY_FOCUS,
        "lennox-gastaut": NeuroWorkflowType.EPILEPSY_FOCUS,
        "absence seizure": NeuroWorkflowType.EPILEPSY_FOCUS,
        # Brain Tumor Grading
        "brain tumor": NeuroWorkflowType.BRAIN_TUMOR,
        "glioblastoma": NeuroWorkflowType.BRAIN_TUMOR,
        "gbm": NeuroWorkflowType.BRAIN_TUMOR,
        "glioma": NeuroWorkflowType.BRAIN_TUMOR,
        "oligodendroglioma": NeuroWorkflowType.BRAIN_TUMOR,
        "astrocytoma": NeuroWorkflowType.BRAIN_TUMOR,
        "meningioma": NeuroWorkflowType.BRAIN_TUMOR,
        "idh": NeuroWorkflowType.BRAIN_TUMOR,
        "mgmt": NeuroWorkflowType.BRAIN_TUMOR,
        "1p/19q": NeuroWorkflowType.BRAIN_TUMOR,
        "1p19q": NeuroWorkflowType.BRAIN_TUMOR,
        "temozolomide": NeuroWorkflowType.BRAIN_TUMOR,
        "stupp": NeuroWorkflowType.BRAIN_TUMOR,
        "rano": NeuroWorkflowType.BRAIN_TUMOR,
        "who grade": NeuroWorkflowType.BRAIN_TUMOR,
        "h3k27m": NeuroWorkflowType.BRAIN_TUMOR,
        "medulloblastoma": NeuroWorkflowType.BRAIN_TUMOR,
        # MS Monitoring
        "multiple sclerosis": NeuroWorkflowType.MS_MONITORING,
        " ms ": NeuroWorkflowType.MS_MONITORING,
        "relapsing remitting": NeuroWorkflowType.MS_MONITORING,
        "rrms": NeuroWorkflowType.MS_MONITORING,
        "spms": NeuroWorkflowType.MS_MONITORING,
        "ppms": NeuroWorkflowType.MS_MONITORING,
        "edss": NeuroWorkflowType.MS_MONITORING,
        "neda": NeuroWorkflowType.MS_MONITORING,
        "disease modifying therapy": NeuroWorkflowType.MS_MONITORING,
        "dmt": NeuroWorkflowType.MS_MONITORING,
        "ocrelizumab": NeuroWorkflowType.MS_MONITORING,
        "natalizumab": NeuroWorkflowType.MS_MONITORING,
        "pml": NeuroWorkflowType.MS_MONITORING,
        "jcv": NeuroWorkflowType.MS_MONITORING,
        "mcdonald criteria": NeuroWorkflowType.MS_MONITORING,
        "t2 lesion": NeuroWorkflowType.MS_MONITORING,
        "gadolinium enhancing": NeuroWorkflowType.MS_MONITORING,
        "neurofilament": NeuroWorkflowType.MS_MONITORING,
        "ofatumumab": NeuroWorkflowType.MS_MONITORING,
        # Parkinson's Assessment
        "parkinson": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "parkinsons": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "updrs": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "hoehn and yahr": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "tremor": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "bradykinesia": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "rigidity": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "levodopa": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "dbs": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "deep brain stimulation": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "dat scan": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "dopamine transporter": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "pigd": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "dyskinesia": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "motor fluctuation": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        "carbidopa": NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        # Headache Classification
        "headache": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "migraine": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "cluster headache": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "tension headache": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "cgrp": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "triptan": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "aura": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "hit-6": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "midas": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "thunderclap": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "medication overuse": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "erenumab": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "fremanezumab": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "galcanezumab": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        "preventive headache": NeuroWorkflowType.HEADACHE_CLASSIFICATION,
        # Neuromuscular Evaluation
        "neuromuscular": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "emg": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "nerve conduction": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "ncs": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "myasthenia": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "als": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "amyotrophic lateral": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "guilain-barre": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "guillain-barre": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "gbs": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "cidp": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "neuropathy": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "myopathy": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "lambert-eaton": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "muscular dystrophy": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "axonal": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "demyelinating neuropathy": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "fasciculation": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
        "alsfrs": NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
    }

    def __init__(self) -> None:
        workflow_instances: List[BaseNeuroWorkflow] = [
            AcuteStrokeTriageWorkflow(),
            DementiaEvaluationWorkflow(),
            EpilepsyFocusWorkflow(),
            BrainTumorGradingWorkflow(),
            MSMonitoringWorkflow(),
            ParkinsonsAssessmentWorkflow(),
            HeadacheClassificationWorkflow(),
            NeuromuscularEvaluationWorkflow(),
        ]
        self._workflows: Dict[NeuroWorkflowType, BaseNeuroWorkflow] = {
            wf.workflow_type: wf for wf in workflow_instances
        }

    # -- public API ----------------------------------------------------

    def run_workflow(
        self, workflow_type: NeuroWorkflowType, inputs: dict
    ) -> WorkflowResult:
        """Execute a specific workflow by type."""
        wf = self._workflows.get(workflow_type)
        if wf is None:
            raise ValueError(f"Unknown workflow type: {workflow_type}")
        return wf.run(inputs)

    def detect_workflow(self, query: str) -> Optional[NeuroWorkflowType]:
        """Detect the most relevant workflow from a free-text query.

        Returns the workflow type with the most keyword matches, or None
        if no keywords match.
        """
        query_lower = query.lower()
        scores: Dict[NeuroWorkflowType, int] = {}
        for keyword, wf_type in self._KEYWORD_MAP.items():
            if keyword in query_lower:
                scores[wf_type] = scores.get(wf_type, 0) + 1
        if not scores:
            return None
        return max(scores, key=scores.get)  # type: ignore[arg-type]

    def get_available_workflows(self) -> List[NeuroWorkflowType]:
        """Return all registered workflow types."""
        return list(self._workflows.keys())
