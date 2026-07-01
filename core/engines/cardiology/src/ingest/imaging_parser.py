"""Cardiac imaging data ingest parser for Cardiology Intelligence Agent.

Provides curated reference data for four cardiac imaging modalities:
  - Echocardiography (transthoracic, transesophageal, stress echo)
  - Cardiac CT (coronary CTA, calcium scoring, structural CT)
  - Cardiac MRI (cine, LGE, T1/T2 mapping, stress perfusion)
  - Nuclear cardiology (SPECT MPI, PET, pyrophosphate scan)

Each seed method produces IngestRecord instances with imaging-specific
metadata (modality, protocol, normal ranges, abnormality criteria,
clinical significance) targeting the ``cardio_imaging`` Milvus collection.

Author: Adam Jones
Date: March 2026
"""

from typing import List

from .base import BaseIngestParser, IngestRecord


# ═══════════════════════════════════════════════════════════════════════
# IMAGING PARSER
# ═══════════════════════════════════════════════════════════════════════


class ImagingParser(BaseIngestParser):
    """Ingest parser for cardiac imaging protocols and reference data.

    Seeds four imaging modality sub-collections into the unified
    ``cardio_imaging`` Milvus collection.  Each modality has its own
    seed method, or all can be loaded at once via :meth:`run`.

    Usage::

        parser = ImagingParser()
        all_records = parser.run()

        # Or individual modalities:
        echo = parser.seed_echo_measurements()
        ct = parser.seed_ct_protocols()
        cmr = parser.seed_cmr_protocols()
        nuc = parser.seed_nuclear_protocols()
    """

    def __init__(self):
        """Initialize the imaging parser."""
        super().__init__("cardio_imaging")

    # ── Fetch / Parse interface ───────────────────────────────────────

    def fetch(self, **kwargs) -> List[dict]:
        """Return all imaging seed data as raw dictionaries.

        Aggregates records from all four modality seed methods.

        Returns:
            List of imaging record dictionaries.
        """
        all_records = []
        all_records.extend(self._echo_data())
        all_records.extend(self._ct_data())
        all_records.extend(self._cmr_data())
        all_records.extend(self._nuclear_data())
        self.logger.info(f"Loaded {len(all_records)} imaging reference entries")
        return all_records

    def parse(self, raw_data: List[dict]) -> List[IngestRecord]:
        """Convert imaging dictionaries into IngestRecord instances.

        Args:
            raw_data: List of imaging dictionaries from :meth:`fetch`.

        Returns:
            List of :class:`IngestRecord` instances.
        """
        records: List[IngestRecord] = []
        for entry in raw_data:
            text = (
                f"{entry.get('modality', '')} — "
                f"{entry.get('finding', entry.get('measurement_name', ''))}. "
                f"Protocol: {entry.get('protocol', '')}. "
                f"Normal: {entry.get('normal_range', '')}. "
                f"Abnormal criteria: {entry.get('abnormal_criteria', '')}. "
                f"Significance: {entry.get('clinical_significance', '')}."
            )

            metadata = {
                "modality": self.truncate(entry.get("modality", ""), 62),
                "protocol": self.truncate(entry.get("protocol", ""), 254),
                "finding": self.truncate(entry.get("finding", ""), 510),
                "measurement_name": self.truncate(
                    entry.get("measurement_name", ""), 126
                ),
                "normal_range": self.truncate(
                    entry.get("normal_range", ""), 126
                ),
                "abnormal_criteria": self.truncate(
                    entry.get("abnormal_criteria", ""), 510
                ),
                "clinical_significance": self.truncate(
                    entry.get("clinical_significance", ""), 1022
                ),
                "guideline_society": self.truncate(
                    entry.get("guideline_society", ""), 62
                ),
                "cross_modal_trigger": entry.get("cross_modal_trigger", False),
            }

            records.append(
                IngestRecord(
                    text=self.truncate(text, 2048),
                    metadata=metadata,
                    collection=self.collection,
                    source=entry.get("guideline_society", "Imaging Reference"),
                )
            )
        return records

    # ── Convenience seed methods ──────────────────────────────────────

    def seed_imaging_protocols(self) -> List[IngestRecord]:
        """Seed all imaging modalities in one call.

        Returns:
            List of valid :class:`IngestRecord` instances.
        """
        return self.run()

    def seed_echo_measurements(self) -> List[IngestRecord]:
        """Seed echocardiography measurements and criteria.

        Returns:
            List of :class:`IngestRecord` for echo data.
        """
        return self.filter_valid(self.parse(self._echo_data()))

    def seed_ct_protocols(self) -> List[IngestRecord]:
        """Seed cardiac CT protocols and scoring.

        Returns:
            List of :class:`IngestRecord` for cardiac CT data.
        """
        return self.filter_valid(self.parse(self._ct_data()))

    def seed_cmr_protocols(self) -> List[IngestRecord]:
        """Seed cardiac MRI protocols and tissue characterization.

        Returns:
            List of :class:`IngestRecord` for CMR data.
        """
        return self.filter_valid(self.parse(self._cmr_data()))

    def seed_nuclear_protocols(self) -> List[IngestRecord]:
        """Seed nuclear cardiology protocols.

        Returns:
            List of :class:`IngestRecord` for nuclear data.
        """
        return self.filter_valid(self.parse(self._nuclear_data()))

    # ── Private data methods ──────────────────────────────────────────

    @staticmethod
    def _echo_data() -> List[dict]:
        """Echocardiography reference data (ASE guidelines).

        Returns:
            List of echo measurement dictionaries.
        """
        return [
            {
                "modality": "echo",
                "protocol": "TTE — 2D and Doppler assessment",
                "measurement_name": "Left Ventricular Ejection Fraction (LVEF)",
                "finding": "LV systolic function assessment",
                "normal_range": "52-72% (men), 54-74% (women)",
                "abnormal_criteria": (
                    "Mildly reduced: 41-51%; Moderately reduced: 30-40%; "
                    "Severely reduced: <30%"
                ),
                "clinical_significance": (
                    "LVEF is the primary determinant for HF classification "
                    "(HFrEF/HFmrEF/HFpEF) and guides GDMT, device therapy "
                    "(ICD, CRT), and transplant evaluation."
                ),
                "guideline_society": "ASE",
                "cross_modal_trigger": True,
            },
            {
                "modality": "echo",
                "protocol": "TTE — Speckle-tracking echocardiography",
                "measurement_name": "Global Longitudinal Strain (GLS)",
                "finding": "LV myocardial deformation assessment",
                "normal_range": "-18% to -22% (more negative = better)",
                "abnormal_criteria": (
                    "Abnormal: > -16% (less negative). Subclinical dysfunction "
                    "may be detected before LVEF decline."
                ),
                "clinical_significance": (
                    "GLS detects subclinical LV dysfunction earlier than LVEF. "
                    "Critical for cardio-oncology surveillance (>15% relative "
                    "decrease from baseline triggers intervention). Prognostic "
                    "in HFpEF and valvular disease."
                ),
                "guideline_society": "ASE/EACVI",
                "cross_modal_trigger": True,
            },
            {
                "modality": "echo",
                "protocol": "TTE — Doppler / tissue Doppler",
                "measurement_name": "E/e' ratio",
                "finding": "LV diastolic function — filling pressure estimate",
                "normal_range": "< 8 (normal filling pressures)",
                "abnormal_criteria": (
                    "Indeterminate: 8-14; Elevated: > 14. Consistent with "
                    "elevated LV filling pressures when >=15."
                ),
                "clinical_significance": (
                    "E/e' is a key non-invasive surrogate for LV filling "
                    "pressure. Elevated E/e' supports diagnosis of HFpEF, "
                    "diastolic dysfunction grading, and pre-procedural "
                    "hemodynamic assessment."
                ),
                "guideline_society": "ASE/EACVI",
                "cross_modal_trigger": False,
            },
            {
                "modality": "echo",
                "protocol": "TTE — Right heart assessment",
                "measurement_name": "TAPSE",
                "finding": "RV systolic function — tricuspid annular excursion",
                "normal_range": ">= 17 mm",
                "abnormal_criteria": (
                    "Abnormal: < 17 mm. Suggests RV systolic dysfunction."
                ),
                "clinical_significance": (
                    "TAPSE is a simple, reproducible measure of RV longitudinal "
                    "function. Reduced TAPSE is prognostic in HF, pulmonary "
                    "hypertension, and post-cardiac surgery."
                ),
                "guideline_society": "ASE",
                "cross_modal_trigger": False,
            },
            {
                "modality": "echo",
                "protocol": "TTE — Aortic valve assessment",
                "measurement_name": "Aortic Valve Area (AVA) / Mean Gradient",
                "finding": "Aortic stenosis severity grading",
                "normal_range": "AVA 3.0-4.0 cm2; Mean gradient < 5 mmHg",
                "abnormal_criteria": (
                    "Mild: AVA 1.5-2.0 cm2, gradient 20-39 mmHg; "
                    "Moderate: AVA 1.0-1.5 cm2; "
                    "Severe: AVA < 1.0 cm2, gradient >= 40 mmHg, Vmax >= 4.0 m/s"
                ),
                "clinical_significance": (
                    "Severe AS with symptoms (angina, syncope, HF) or LVEF <50% "
                    "warrants AVR (SAVR or TAVR). Low-flow low-gradient AS "
                    "requires dobutamine stress echo for confirmation."
                ),
                "guideline_society": "ASE/ACC/AHA",
                "cross_modal_trigger": True,
            },
            {
                "modality": "echo",
                "protocol": "TEE — Transesophageal echocardiography",
                "measurement_name": "LAA thrombus assessment",
                "finding": "Left atrial appendage evaluation pre-cardioversion",
                "normal_range": "No thrombus; LAA velocity > 40 cm/s",
                "abnormal_criteria": (
                    "Thrombus present; spontaneous echo contrast (SEC); "
                    "LAA velocity < 20 cm/s (high thrombus risk)"
                ),
                "clinical_significance": (
                    "TEE is required before cardioversion in AF > 48 hours "
                    "without adequate anticoagulation. LAA thrombus is an "
                    "absolute contraindication to cardioversion."
                ),
                "guideline_society": "ASE/ACC",
                "cross_modal_trigger": False,
            },
            {
                "modality": "echo",
                "protocol": "TTE — Right heart assessment",
                "measurement_name": "RV Fractional Area Change (FAC)",
                "finding": "RV systolic function — area-based assessment",
                "normal_range": ">35%",
                "abnormal_criteria": (
                    "Borderline: 25-35%; Abnormal: <25%. Measured in apical "
                    "4-chamber view: (end-diastolic area - end-systolic area) "
                    "/ end-diastolic area x 100."
                ),
                "clinical_significance": (
                    "RV FAC is a simple 2D measure of RV systolic function. "
                    "Better correlates with CMR RVEF than TAPSE. Prognostic "
                    "in pulmonary hypertension, HF, and post-MI RV involvement."
                ),
                "guideline_society": "ASE",
                "cross_modal_trigger": False,
            },
            {
                "modality": "echo",
                "protocol": "TTE — Left atrial assessment",
                "measurement_name": "LA Volume Index (LAVI)",
                "finding": "Left atrial size — diastolic dysfunction marker",
                "normal_range": "<34 mL/m2",
                "abnormal_criteria": (
                    "Mildly enlarged: 34-41 mL/m2; Moderately enlarged: "
                    "42-48 mL/m2; Severely enlarged: >48 mL/m2. "
                    "Indexed to body surface area."
                ),
                "clinical_significance": (
                    "LAVI is a marker of chronic diastolic dysfunction and "
                    "reflects cumulative LV filling pressure burden. Enlarged "
                    "LAVI predicts AF, stroke, HF, and cardiovascular death. "
                    "Key parameter in ASE diastolic function grading algorithm."
                ),
                "guideline_society": "ASE/EACVI",
                "cross_modal_trigger": False,
            },
            {
                "modality": "echo",
                "protocol": "TTE — Doppler / TR jet velocity",
                "measurement_name": "Pulmonary Artery Systolic Pressure (PASP)",
                "finding": "Pulmonary pressure estimation via TR jet",
                "normal_range": "<35 mmHg",
                "abnormal_criteria": (
                    "Mild PH: 36-50 mmHg; Moderate PH: 51-60 mmHg; "
                    "Severe PH: >60 mmHg. Estimated by modified Bernoulli "
                    "equation: 4 x (TR jet velocity)^2 + RAP."
                ),
                "clinical_significance": (
                    "Elevated PASP is a screening indicator for pulmonary "
                    "hypertension. Right heart catheterization is required "
                    "for definitive diagnosis (mPAP >20 mmHg). Echo may "
                    "over- or underestimate pressures by 10+ mmHg."
                ),
                "guideline_society": "ASE/ESC",
                "cross_modal_trigger": True,
            },
            {
                "modality": "echo",
                "protocol": "TTE — 2D and M-mode",
                "measurement_name": "LV Mass Index",
                "finding": "LV hypertrophy assessment",
                "normal_range": (
                    "Male: <115 g/m2; Female: <95 g/m2"
                ),
                "abnormal_criteria": (
                    "LVH present when LV mass index exceeds sex-specific "
                    "thresholds. Concentric hypertrophy: increased mass with "
                    "normal cavity. Eccentric: increased mass with dilated "
                    "cavity. Concentric remodeling: normal mass, increased RWT."
                ),
                "clinical_significance": (
                    "LVH is an independent predictor of cardiovascular events. "
                    "Geometry pattern (concentric vs eccentric) has prognostic "
                    "implications. Key criterion for hypertension target organ "
                    "damage and HCM assessment."
                ),
                "guideline_society": "ASE",
                "cross_modal_trigger": True,
            },
            {
                "modality": "echo",
                "protocol": "TTE — Speckle-tracking echocardiography",
                "measurement_name": "RV Free Wall Longitudinal Strain (RVFWLS)",
                "finding": "RV systolic function — strain-based assessment",
                "normal_range": "> -20% (more negative = better function)",
                "abnormal_criteria": (
                    "Abnormal: > -20% (less negative). Early RV dysfunction "
                    "can be detected by strain before TAPSE or FAC become "
                    "abnormal."
                ),
                "clinical_significance": (
                    "RVFWLS provides early detection of RV dysfunction before "
                    "conventional parameters decline. Prognostic in pulmonary "
                    "hypertension, HF, and post-PE. Superior to TAPSE for "
                    "detecting subclinical RV dysfunction."
                ),
                "guideline_society": "ASE/EACVI",
                "cross_modal_trigger": True,
            },
            {
                "modality": "echo",
                "protocol": "TTE — Doppler / tissue Doppler",
                "measurement_name": "Diastolic Function Grading",
                "finding": "LV diastolic function — integrated grading",
                "normal_range": "Normal: E/A 0.8-2.0, DT 160-200 ms, e' >=10 cm/s",
                "abnormal_criteria": (
                    "Grade I (impaired relaxation): E/A <0.8, prolonged DT, "
                    "normal filling pressures. Grade II (pseudonormal): E/A "
                    "0.8-1.5, E/e' 10-14, requires tissue Doppler to unmask. "
                    "Grade III (restrictive): E/A >=2, DT <160 ms, E/e' >=15, "
                    "elevated filling pressures."
                ),
                "clinical_significance": (
                    "Diastolic dysfunction grading integrates E/A ratio, DT, "
                    "E/e', LAVI, and TR velocity per ASE/EACVI 2016 algorithm. "
                    "Grade II-III associated with HFpEF, increased mortality, "
                    "and hospitalization."
                ),
                "guideline_society": "ASE/EACVI",
                "cross_modal_trigger": False,
            },
        ]

    @staticmethod
    def _ct_data() -> List[dict]:
        """Cardiac CT reference data (SCCT guidelines).

        Returns:
            List of cardiac CT protocol dictionaries.
        """
        return [
            {
                "modality": "CT",
                "protocol": "Coronary CT Angiography (CCTA)",
                "measurement_name": "CAD-RADS score",
                "finding": "Coronary artery stenosis grading",
                "normal_range": "CAD-RADS 0: No stenosis",
                "abnormal_criteria": (
                    "CAD-RADS 1: 1-24% (minimal); 2: 25-49% (mild); "
                    "3: 50-69% (moderate); 4A: 70-99% single vessel (severe); "
                    "4B: LM >50% or 3-vessel >70%; 5: Total occlusion"
                ),
                "clinical_significance": (
                    "CAD-RADS standardises CCTA reporting. CAD-RADS >=3 "
                    "warrants functional testing or invasive angiography. "
                    "CAD-RADS 4B may proceed directly to catheterisation."
                ),
                "guideline_society": "SCCT",
                "cross_modal_trigger": True,
            },
            {
                "modality": "CT",
                "protocol": "Coronary Artery Calcium (CAC) Scoring",
                "measurement_name": "Agatston Score",
                "finding": "Coronary calcium quantification for risk stratification",
                "normal_range": "CAC = 0 (very low 10-year ASCVD risk)",
                "abnormal_criteria": (
                    "1-99: Mild (mildly increased risk); "
                    "100-299: Moderate (moderately increased risk); "
                    ">=300: Severe (high risk, >75th percentile for age/sex)"
                ),
                "clinical_significance": (
                    "CAC = 0 can down-classify borderline-risk patients, "
                    "deferring statin therapy. CAC >= 100 supports statin "
                    "initiation. CAC is a powerful negative predictor: "
                    "10-year ASCVD event rate < 1% when CAC = 0."
                ),
                "guideline_society": "SCCT/ACC/AHA",
                "cross_modal_trigger": False,
            },
            {
                "modality": "CT",
                "protocol": "Structural CT — TAVR Planning",
                "measurement_name": "Aortic annulus dimensions",
                "finding": "Pre-TAVR anatomical assessment",
                "normal_range": "Annulus area 314-621 mm2 (device sizing dependent)",
                "abnormal_criteria": (
                    "Annulus too small (<314 mm2) or too large (>683 mm2) "
                    "for available prostheses; bicuspid morphology; "
                    "extensive LVOT calcification"
                ),
                "clinical_significance": (
                    "Multidetector CT is mandatory for TAVR planning: "
                    "annulus sizing, access route assessment (iliofemoral "
                    "dimensions >=5.0 mm), coronary height, and calcium "
                    "distribution determine procedural approach and risk."
                ),
                "guideline_society": "SCCT/ACC",
                "cross_modal_trigger": False,
            },
        ]

    @staticmethod
    def _cmr_data() -> List[dict]:
        """Cardiac MRI reference data (SCMR guidelines).

        Returns:
            List of CMR protocol dictionaries.
        """
        return [
            {
                "modality": "CMR",
                "protocol": "Cine SSFP — Volumetric analysis",
                "measurement_name": "CMR LVEF / LV volumes",
                "finding": "Gold-standard LV function assessment",
                "normal_range": (
                    "LVEF 57-77% (men), 61-79% (women); "
                    "LVEDV 106-214 mL (men), 76-156 mL (women)"
                ),
                "abnormal_criteria": (
                    "Mildly reduced LVEF: 46-56% (men), 50-60% (women); "
                    "Dilated LVEDV: >214 mL (men), >156 mL (women)"
                ),
                "clinical_significance": (
                    "CMR is the gold standard for LV volume and EF "
                    "quantification. Particularly valuable when echo "
                    "windows are limited or discrepant values affect "
                    "clinical decisions (ICD implantation, transplant)."
                ),
                "guideline_society": "SCMR",
                "cross_modal_trigger": False,
            },
            {
                "modality": "CMR",
                "protocol": "Late Gadolinium Enhancement (LGE)",
                "measurement_name": "Myocardial scar / fibrosis pattern",
                "finding": "Tissue characterisation — scar detection",
                "normal_range": "No LGE enhancement",
                "abnormal_criteria": (
                    "Subendocardial: ischemic (CAD); "
                    "Mid-wall: non-ischemic (dilated CMP, myocarditis); "
                    "Epicardial: myocarditis, sarcoidosis; "
                    "RV insertion: PH overload; "
                    "Diffuse: amyloidosis, Anderson-Fabry"
                ),
                "clinical_significance": (
                    "LGE pattern differentiates ischemic from non-ischemic "
                    "cardiomyopathy, guides arrhythmia risk stratification, "
                    "and identifies myocarditis. LGE extent >15% of LV mass "
                    "increases SCD risk even with preserved LVEF."
                ),
                "guideline_society": "SCMR",
                "cross_modal_trigger": True,
            },
            {
                "modality": "CMR",
                "protocol": "T1 Mapping (native and post-contrast)",
                "measurement_name": "Native T1 / Extracellular Volume (ECV)",
                "finding": "Diffuse myocardial fibrosis quantification",
                "normal_range": (
                    "Native T1 (1.5T): 950-1050 ms; ECV: 25-30%"
                ),
                "abnormal_criteria": (
                    "Elevated T1 / ECV: diffuse fibrosis (HFpEF, amyloid, "
                    "myocarditis); Low T1: iron overload, Anderson-Fabry; "
                    "ECV > 40%: strongly suggestive of amyloidosis"
                ),
                "clinical_significance": (
                    "T1 mapping detects diffuse fibrosis not visible on LGE. "
                    "ECV > 40% is highly specific for cardiac amyloidosis. "
                    "Native T1 elevation in HFpEF correlates with outcomes."
                ),
                "guideline_society": "SCMR",
                "cross_modal_trigger": True,
            },
            {
                "modality": "CMR",
                "protocol": "Stress Perfusion CMR (Adenosine / Regadenoson)",
                "measurement_name": "Myocardial perfusion defect",
                "finding": "Inducible ischemia assessment",
                "normal_range": "Homogeneous myocardial perfusion, no defects",
                "abnormal_criteria": (
                    "Perfusion defect in >=2 contiguous segments; "
                    "matched defect (scar) vs unmatched (ischemia); "
                    "ischemic burden >=10% of LV myocardium"
                ),
                "clinical_significance": (
                    "Stress perfusion CMR has sensitivity >90% and specificity "
                    ">80% for detecting significant CAD. MR-INFORM showed "
                    "CMR-guided management is non-inferior to FFR-guided "
                    "management for MACE outcomes."
                ),
                "guideline_society": "SCMR/ACC",
                "cross_modal_trigger": True,
            },
            {
                "modality": "CMR",
                "protocol": "T2 Mapping",
                "measurement_name": "Myocardial T2 values",
                "finding": "Myocardial edema detection",
                "normal_range": "Normal T2 <50 ms at 1.5T",
                "abnormal_criteria": (
                    "Elevated T2 (>50 ms at 1.5T): myocardial edema in acute "
                    "myocarditis, Takotsubo cardiomyopathy, acute MI. "
                    "Anderson-Fabry disease shows LOW T1 with normal T2 "
                    "(distinguishing feature)."
                ),
                "clinical_significance": (
                    "T2 mapping detects acute myocardial inflammation and "
                    "edema, distinguishing acute from chronic processes. "
                    "Combined T1/T2 mapping improves diagnostic accuracy for "
                    "myocarditis (Lake Louise criteria). Normal T2 in Fabry "
                    "with low T1 is a key diagnostic pattern."
                ),
                "guideline_society": "SCMR",
                "cross_modal_trigger": True,
            },
            {
                "modality": "CMR",
                "protocol": "4D Flow MRI",
                "measurement_name": "Volumetric flow quantification",
                "finding": "Velocity-encoded volumetric flow assessment",
                "normal_range": "Qp/Qs = 1.0 (no shunt); normal aortic flow patterns",
                "abnormal_criteria": (
                    "Qp/Qs >1.5: significant left-to-right shunt. "
                    "Abnormal regurgitant fraction (>20% mild, >40% severe). "
                    "Elevated aortic wall shear stress associated with "
                    "aortopathy and bicuspid aortic valve."
                ),
                "clinical_significance": (
                    "4D flow MRI provides comprehensive hemodynamic assessment "
                    "in a single acquisition. Quantifies regurgitation volumes, "
                    "shunt ratios (Qp/Qs), and aortic wall shear stress. "
                    "Valuable for congenital heart disease, valvular disease, "
                    "and aortopathy evaluation."
                ),
                "guideline_society": "SCMR",
                "cross_modal_trigger": True,
            },
            {
                "modality": "CMR",
                "protocol": "Fat Imaging — Dixon Technique",
                "measurement_name": "Myocardial fat content",
                "finding": "Myocardial fatty infiltration assessment",
                "normal_range": "No myocardial fat signal on fat-only images",
                "abnormal_criteria": (
                    "Fat signal in RV free wall: ARVC (major criterion). "
                    "Fat in subendocardial regions: lipomatous metaplasia "
                    "post-MI (chronic infarct). Fat in interventricular "
                    "septum: lipomatous hypertrophy."
                ),
                "clinical_significance": (
                    "Dixon fat-water separation confirms fatty infiltration "
                    "in ARVC (2010 Task Force Criteria), identifies "
                    "lipomatous metaplasia in chronic MI scars, and detects "
                    "cardiac lipomas. Complementary to LGE for tissue "
                    "characterization."
                ),
                "guideline_society": "SCMR",
                "cross_modal_trigger": True,
            },
        ]

    @staticmethod
    def _nuclear_data() -> List[dict]:
        """Nuclear cardiology reference data (ASNC guidelines).

        Returns:
            List of nuclear cardiology protocol dictionaries.
        """
        return [
            {
                "modality": "nuclear",
                "protocol": "SPECT MPI (Tc-99m sestamibi / tetrofosmin)",
                "measurement_name": "Summed Stress Score (SSS) / Ischemic Burden",
                "finding": "Myocardial perfusion imaging for CAD detection",
                "normal_range": "SSS 0-3 (normal perfusion)",
                "abnormal_criteria": (
                    "SSS 4-8: mildly abnormal; SSS 9-13: moderately abnormal; "
                    "SSS >=14: severely abnormal. Ischemic burden >10% "
                    "associated with revascularization benefit."
                ),
                "clinical_significance": (
                    "SPECT MPI remains the most widely used non-invasive "
                    "ischemia test. Normal SPECT confers <1%/year MACE risk. "
                    "Ischemic burden >10% may benefit from revascularization "
                    "per ISCHEMIA-informed shared decision-making."
                ),
                "guideline_society": "ASNC",
                "cross_modal_trigger": True,
            },
            {
                "modality": "nuclear",
                "protocol": "Cardiac PET (Rubidium-82 / N-13 ammonia)",
                "measurement_name": "Myocardial Blood Flow (MBF) / Flow Reserve",
                "finding": "Quantitative perfusion and coronary flow reserve",
                "normal_range": (
                    "Rest MBF: 0.6-1.0 mL/g/min; "
                    "Stress MBF: 2.5-3.5 mL/g/min; "
                    "CFR >= 2.0"
                ),
                "abnormal_criteria": (
                    "CFR < 2.0: reduced flow reserve (epicardial CAD or "
                    "microvascular dysfunction); CFR < 1.5: severely impaired"
                ),
                "clinical_significance": (
                    "Cardiac PET provides absolute quantification of "
                    "myocardial blood flow, enabling detection of balanced "
                    "ischemia missed by relative perfusion imaging. "
                    "CFR < 2.0 is an independent predictor of cardiac events."
                ),
                "guideline_society": "ASNC",
                "cross_modal_trigger": True,
            },
            {
                "modality": "nuclear",
                "protocol": "Tc-99m Pyrophosphate (PYP) Scan",
                "measurement_name": "Heart-to-Contralateral (H/CL) ratio",
                "finding": "Cardiac amyloidosis screening — ATTR type",
                "normal_range": "H/CL ratio < 1.0 (no uptake)",
                "abnormal_criteria": (
                    "H/CL >= 1.5: highly suggestive of ATTR cardiac amyloid; "
                    "Visual grade 2-3 at 1 hour: positive. "
                    "Requires serum/urine immunofixation to exclude AL amyloid."
                ),
                "clinical_significance": (
                    "PYP scan is >97% specific for ATTR cardiac amyloidosis "
                    "when AL amyloid is excluded by serum immunofixation. "
                    "Enables non-invasive diagnosis without endomyocardial "
                    "biopsy. Critical given availability of tafamidis therapy."
                ),
                "guideline_society": "ASNC/ACC",
                "cross_modal_trigger": True,
            },
            {
                "modality": "nuclear",
                "protocol": "FDG PET for Cardiac Sarcoidosis",
                "measurement_name": "FDG uptake pattern",
                "finding": "Cardiac sarcoidosis — active inflammation detection",
                "normal_range": (
                    "No focal FDG uptake after adequate myocardial "
                    "suppression (72-hour high-fat/low-carb diet)"
                ),
                "abnormal_criteria": (
                    "Focal FDG uptake = active granulomatous inflammation. "
                    "Focal on diffuse pattern: active inflammation with "
                    "inadequate suppression. Combined FDG + perfusion defect: "
                    "mismatch pattern (inflammation + scar)."
                ),
                "clinical_significance": (
                    "FDG PET is the primary non-invasive tool for detecting "
                    "active cardiac sarcoidosis. Requires strict 72-hour "
                    "high-fat, low-carbohydrate diet preparation to suppress "
                    "normal myocardial glucose uptake. Guides initiation and "
                    "monitoring of immunosuppressive therapy. Follow-up PET "
                    "assesses treatment response."
                ),
                "guideline_society": "ASNC/HRS",
                "cross_modal_trigger": True,
            },
            {
                "modality": "nuclear",
                "protocol": "I-123 MIBG Cardiac Sympathetic Innervation",
                "measurement_name": "Heart-to-Mediastinum (H/M) ratio",
                "finding": "Cardiac sympathetic innervation assessment",
                "normal_range": "H/M ratio >1.6 (late imaging at 4 hours)",
                "abnormal_criteria": (
                    "H/M ratio <=1.6: reduced cardiac sympathetic innervation. "
                    "Lower H/M ratio correlates with worse prognosis. "
                    "ADMIRE-HF trial: H/M <1.6 associated with 2-year cardiac "
                    "event rate of 37% vs 15%."
                ),
                "clinical_significance": (
                    "MIBG imaging quantifies cardiac sympathetic denervation. "
                    "Prognostic in HFrEF: H/M ratio <1.6 identifies patients "
                    "at higher risk of arrhythmic death. Aids ICD "
                    "decision-making in borderline candidates. ADMIRE-HF "
                    "trial validated prognostic value."
                ),
                "guideline_society": "ASNC",
                "cross_modal_trigger": True,
            },
        ]
