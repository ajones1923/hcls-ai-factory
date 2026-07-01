"""Hemodynamic reference data ingest parser for Cardiology Intelligence Agent.

Provides curated reference data for invasive and non-invasive hemodynamic
parameters measured during right heart catheterisation, left heart
catheterisation, and echocardiographic Doppler assessment.

Also includes catheterisation lab protocols for common diagnostic and
interventional procedures.

Targets the ``cardio_hemodynamics`` Milvus collection.

Author: Adam Jones
Date: March 2026
"""

from typing import List

from .base import BaseIngestParser, IngestRecord


# ═══════════════════════════════════════════════════════════════════════
# HEMODYNAMICS PARSER
# ═══════════════════════════════════════════════════════════════════════


class HemodynamicsParser(BaseIngestParser):
    """Ingest parser for hemodynamic parameters and cath lab protocols.

    Seeds two sub-categories into the ``cardio_hemodynamics`` Milvus
    collection:
      1. Hemodynamic parameters (:meth:`seed_hemodynamic_parameters`)
      2. Cath lab protocols (:meth:`seed_cathlab_protocols`)

    Usage::

        parser = HemodynamicsParser()
        all_records = parser.run()

        # Or individually:
        params = parser.seed_hemodynamic_parameters()
        protocols = parser.seed_cathlab_protocols()
    """

    def __init__(self):
        """Initialize the hemodynamics parser."""
        super().__init__("cardio_hemodynamics")

    # ── Fetch / Parse interface ───────────────────────────────────────

    def fetch(self, **kwargs) -> List[dict]:
        """Return all hemodynamic seed data as raw dictionaries.

        Returns:
            Combined list of parameter and protocol data.
        """
        all_data: List[dict] = []
        all_data.extend(self._parameter_data())
        all_data.extend(self._protocol_data())
        self.logger.info(
            f"Loaded {len(all_data)} hemodynamic reference entries"
        )
        return all_data

    def parse(self, raw_data: List[dict]) -> List[IngestRecord]:
        """Convert hemodynamic dictionaries into IngestRecord instances.

        Args:
            raw_data: List of hemodynamic dictionaries from :meth:`fetch`.

        Returns:
            List of :class:`IngestRecord` instances.
        """
        records: List[IngestRecord] = []
        for entry in raw_data:
            text = (
                f"{entry.get('parameter_name', '')}. "
                f"Method: {entry.get('measurement_method', '')}. "
                f"Normal: {entry.get('normal_range', '')}. "
                f"Abnormal: {entry.get('abnormal_criteria', '')}. "
                f"Significance: {entry.get('clinical_significance', '')}. "
                f"Formula: {entry.get('calculation_formula', '')}."
            )

            metadata = {
                "parameter_name": self.truncate(
                    entry.get("parameter_name", ""), 126
                ),
                "measurement_method": self.truncate(
                    entry.get("measurement_method", ""), 254
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
                "calculation_formula": self.truncate(
                    entry.get("calculation_formula", ""), 510
                ),
                "cathlab_context": self.truncate(
                    entry.get("cathlab_context", ""), 62
                ),
            }

            records.append(
                IngestRecord(
                    text=self.truncate(text, 2048),
                    metadata=metadata,
                    collection=self.collection,
                    source="Hemodynamic Reference",
                )
            )
        return records

    # ── Convenience seed methods ──────────────────────────────────────

    def seed_hemodynamic_parameters(self) -> List[IngestRecord]:
        """Seed hemodynamic parameter reference data.

        Returns:
            List of :class:`IngestRecord` for hemodynamic parameters.
        """
        return self.filter_valid(self.parse(self._parameter_data()))

    def seed_cathlab_protocols(self) -> List[IngestRecord]:
        """Seed catheterisation lab protocol data.

        Returns:
            List of :class:`IngestRecord` for cath lab protocols.
        """
        return self.filter_valid(self.parse(self._protocol_data()))

    # ── Private data methods ──────────────────────────────────────────

    @staticmethod
    def _parameter_data() -> List[dict]:
        """Core hemodynamic parameters (Swan-Ganz, Fick, echo-Doppler).

        Returns:
            List of hemodynamic parameter dictionaries.
        """
        return [
            {
                "parameter_name": "Pulmonary Capillary Wedge Pressure (PCWP)",
                "measurement_method": "Swan-Ganz catheter, balloon-tipped in pulmonary artery",
                "normal_range": "6-12 mmHg",
                "abnormal_criteria": (
                    "Elevated >12 mmHg: volume overload; "
                    ">18 mmHg: pulmonary congestion; "
                    ">25 mmHg: pulmonary oedema threshold. "
                    "Low <6 mmHg: hypovolemia."
                ),
                "clinical_significance": (
                    "PCWP is the primary surrogate for left atrial pressure "
                    "and LV filling pressure. Distinguishes pre-capillary "
                    "from post-capillary pulmonary hypertension (PCWP >15 = "
                    "post-capillary). Critical for heart transplant evaluation "
                    "and HF hemodynamic tailoring."
                ),
                "calculation_formula": "Direct measurement (wedge tracing after balloon inflation)",
                "cathlab_context": "right_heart_cath",
            },
            {
                "parameter_name": "Cardiac Output (CO) / Cardiac Index (CI)",
                "measurement_method": "Thermodilution (Swan-Ganz) or Fick method",
                "normal_range": "CO: 4.0-8.0 L/min; CI: 2.5-4.0 L/min/m2",
                "abnormal_criteria": (
                    "CI < 2.2 L/min/m2: cardiogenic shock if with elevated "
                    "PCWP; CI < 1.8 L/min/m2: severe low-output state; "
                    "CI > 4.0 L/min/m2: high-output state (sepsis, AV fistula, "
                    "thyrotoxicosis)."
                ),
                "clinical_significance": (
                    "CI is the primary determinant of organ perfusion. "
                    "Low CI with elevated PCWP = cardiogenic shock (wet-cold). "
                    "Low CI with low PCWP = hypovolemia (dry-cold). "
                    "Guides vasopressor/inotrope selection in shock."
                ),
                "calculation_formula": (
                    "Fick: CO = VO2 / (CaO2 - CvO2). "
                    "Thermodilution: CO = V_injectate * (Tb - Ti) * K / "
                    "integral(delta_T * dt). CI = CO / BSA."
                ),
                "cathlab_context": "right_heart_cath",
            },
            {
                "parameter_name": "Systemic Vascular Resistance (SVR)",
                "measurement_method": "Calculated from MAP, CVP, and CO",
                "normal_range": "800-1200 dynes*s/cm5 (10-16 Wood units)",
                "abnormal_criteria": (
                    "High SVR (>1200): cardiogenic shock, hypovolemia, "
                    "vasoconstrictive states. Low SVR (<800): distributive "
                    "shock (sepsis), cirrhosis, AV fistula."
                ),
                "clinical_significance": (
                    "SVR reflects LV afterload. Elevated SVR in HF increases "
                    "myocardial work and worsens forward flow. SVR-guided "
                    "therapy: vasodilators (nitroprusside) for high SVR, "
                    "vasopressors (norepinephrine) for low SVR."
                ),
                "calculation_formula": "SVR = (MAP - CVP) / CO * 80 dynes*s/cm5",
                "cathlab_context": "right_heart_cath",
            },
            {
                "parameter_name": "Pulmonary Vascular Resistance (PVR)",
                "measurement_method": "Calculated from mPAP, PCWP, and CO",
                "normal_range": "< 3 Wood units (< 240 dynes*s/cm5)",
                "abnormal_criteria": (
                    "Mild PH: 3-5 WU; Moderate: 5-8 WU; Severe: >8 WU. "
                    "PVR > 5 WU with PCWP <=15: pre-capillary PH. "
                    "Fixed PVR > 6 WU: relative contraindication to heart "
                    "transplant."
                ),
                "clinical_significance": (
                    "PVR differentiates pre-capillary (PAH) from post-capillary "
                    "(left heart disease) pulmonary hypertension. "
                    "Elevated PVR in HF indicates secondary PH and impacts "
                    "transplant candidacy. Vasoreactivity testing with iNO "
                    "or prostacyclin guides PAH therapy."
                ),
                "calculation_formula": "PVR = (mPAP - PCWP) / CO in Wood units",
                "cathlab_context": "right_heart_cath",
            },
            {
                "parameter_name": "Mean Pulmonary Artery Pressure (mPAP)",
                "measurement_method": "Swan-Ganz catheter, direct measurement",
                "normal_range": "9-18 mmHg (new definition: normal <= 20 mmHg)",
                "abnormal_criteria": (
                    "PH defined as mPAP > 20 mmHg (ESC/ERS 2022, revised "
                    "from > 25 mmHg). Pre-capillary: mPAP > 20 + PVR > 2 WU "
                    "+ PCWP <= 15. Post-capillary: mPAP > 20 + PCWP > 15."
                ),
                "clinical_significance": (
                    "mPAP is the defining hemodynamic parameter for pulmonary "
                    "hypertension. Combined with PCWP and PVR, classifies PH "
                    "into WHO groups I-V. Prognosis correlates with mPAP "
                    "and response to therapy."
                ),
                "calculation_formula": "mPAP = (sPAP + 2*dPAP) / 3. Estimated by echo: mPAP = 0.61 * sPAP + 2",
                "cathlab_context": "right_heart_cath",
            },
            {
                "parameter_name": "Central Venous Pressure (CVP) / Right Atrial Pressure",
                "measurement_method": "Right heart catheter or central venous line",
                "normal_range": "2-8 mmHg (mean RA pressure)",
                "abnormal_criteria": (
                    "Elevated >12 mmHg: RV failure, fluid overload, "
                    "tricuspid regurgitation, constrictive pericarditis, "
                    "cardiac tamponade. Low <2 mmHg: hypovolemia."
                ),
                "clinical_significance": (
                    "CVP reflects RV preload and right-sided filling pressures. "
                    "Elevated CVP with low PCWP: isolated RV failure. "
                    "Elevated CVP with elevated PCWP: biventricular failure. "
                    "Kussmaul sign (CVP rise with inspiration): constrictive "
                    "pericarditis."
                ),
                "calculation_formula": "Direct measurement (mean of A-wave tracing)",
                "cathlab_context": "right_heart_cath",
            },
            {
                "parameter_name": "Fractional Flow Reserve (FFR)",
                "measurement_method": "Intracoronary pressure wire during maximal hyperemia",
                "normal_range": "FFR = 1.0 (no stenosis)",
                "abnormal_criteria": (
                    "FFR <= 0.80: hemodynamically significant stenosis "
                    "(revascularization benefit). FFR 0.81-0.85: grey zone "
                    "(consider iFR or clinical context). FFR > 0.85: defer "
                    "revascularization."
                ),
                "clinical_significance": (
                    "FFR is the invasive gold standard for physiological "
                    "assessment of coronary stenosis significance. "
                    "FAME trials: FFR-guided PCI reduces MACE vs angiography-guided. "
                    "FFR-guided revascularization is a Class I recommendation."
                ),
                "calculation_formula": "FFR = Pd / Pa (distal coronary pressure / aortic pressure during hyperemia)",
                "cathlab_context": "diagnostic_cath",
            },
            {
                "parameter_name": "Aortic Valve Area (AVA) — Gorlin Formula",
                "measurement_method": "Simultaneous LV and aortic pressure measurement + CO",
                "normal_range": "3.0-4.0 cm2",
                "abnormal_criteria": (
                    "Mild AS: AVA 1.5-2.0 cm2; Moderate: 1.0-1.5 cm2; "
                    "Severe: < 1.0 cm2 (or < 0.6 cm2/m2 indexed)."
                ),
                "clinical_significance": (
                    "Invasive AVA by Gorlin formula is the reference standard "
                    "when non-invasive assessment is inconclusive. "
                    "Low-flow low-gradient AS: dobutamine challenge to "
                    "distinguish true-severe from pseudo-severe AS."
                ),
                "calculation_formula": (
                    "AVA = CO / (44.3 * SEP * HR * sqrt(mean_gradient)). "
                    "SEP = systolic ejection period."
                ),
                "cathlab_context": "diagnostic_cath",
            },
            {
                "parameter_name": "Instantaneous Wave-Free Ratio (iFR)",
                "measurement_method": (
                    "Resting coronary pressure measurement during the "
                    "wave-free period of diastole (no adenosine required)"
                ),
                "normal_range": "iFR > 0.89 (no intervention needed)",
                "abnormal_criteria": (
                    "iFR <= 0.89: hemodynamically significant stenosis, "
                    "revascularization recommended. Validated non-inferior "
                    "to FFR by DEFINE-FLAIR and iFR-SWEDEHEART trials."
                ),
                "clinical_significance": (
                    "Alternative to FFR that does not require adenosine "
                    "administration. Class I recommendation (ACC/AHA 2021). "
                    "DEFINE-FLAIR and iFR-SWEDEHEART trials demonstrated "
                    "non-inferiority to FFR-guided revascularization for "
                    "MACE outcomes at 1 year."
                ),
                "calculation_formula": "iFR = Pd / Pa during the wave-free period of diastole",
                "cathlab_context": "diagnostic_cath",
            },
            {
                "parameter_name": "dP/dt (Rate of LV Pressure Rise)",
                "measurement_method": (
                    "Derived from mitral regurgitation jet CW Doppler "
                    "(1 m/s to 3 m/s interval) or direct LV pressure tracing"
                ),
                "normal_range": "> 1200 mmHg/s",
                "abnormal_criteria": (
                    "dP/dt < 1200 mmHg/s: LV systolic dysfunction; "
                    "dP/dt < 800 mmHg/s: severe LV systolic dysfunction. "
                    "Values correlate with invasive contractility indices."
                ),
                "clinical_significance": (
                    "Relatively load-independent measure of LV contractility. "
                    "Predicts response to inotropic therapy. Useful for "
                    "serial assessment of LV function in patients with MR. "
                    "Non-invasive estimation: dP/dt = 32 mmHg / time(1m/s to 3m/s)."
                ),
                "calculation_formula": (
                    "dP/dt = 32 mmHg / delta_t (where delta_t is the time "
                    "interval for MR jet velocity to rise from 1 m/s to 3 m/s)"
                ),
                "cathlab_context": "diagnostic_cath",
            },
            {
                "parameter_name": "Transpulmonary Gradient (TPG)",
                "measurement_method": "Calculated from mPAP and PCWP (Swan-Ganz catheter)",
                "normal_range": "< 12 mmHg",
                "abnormal_criteria": (
                    "TPG > 12 mmHg: pre-capillary component to pulmonary "
                    "hypertension present. TPG > 15 mmHg: significant "
                    "pre-capillary PH component; relative contraindication "
                    "to heart transplantation without vasoreactivity testing."
                ),
                "clinical_significance": (
                    "Differentiates reactive (reversible) from fixed pulmonary "
                    "hypertension. Critical transplant eligibility criterion. "
                    "Combined with DPG (diastolic pulmonary gradient) for "
                    "comprehensive PH classification per ESC/ERS 2022."
                ),
                "calculation_formula": "TPG = mPAP - PCWP",
                "cathlab_context": "right_heart_cath",
            },
            {
                "parameter_name": "Pulmonary Artery Pulsatility Index (PAPi)",
                "measurement_method": (
                    "Calculated from systolic PA pressure, diastolic PA "
                    "pressure, and CVP (Swan-Ganz catheter)"
                ),
                "normal_range": "> 1.85",
                "abnormal_criteria": (
                    "PAPi < 1.85: RV dysfunction; PAPi < 1.0: severe RV "
                    "failure with high mortality risk. PAPi < 0.9 post-MI: "
                    "predictor of need for RV mechanical support."
                ),
                "clinical_significance": (
                    "Predicts RV failure post-LVAD implantation and post-MI. "
                    "Guides mechanical circulatory support candidacy. "
                    "Low PAPi identifies patients who may need biventricular "
                    "support rather than isolated LVAD."
                ),
                "calculation_formula": "PAPi = (sPAP - dPAP) / CVP",
                "cathlab_context": "right_heart_cath",
            },
            {
                "parameter_name": "Fick Cardiac Output",
                "measurement_method": (
                    "Assumed or directly measured VO2, arterial and mixed "
                    "venous (PA) oxygen content sampling"
                ),
                "normal_range": "CO: 4.0-8.0 L/min; CI: 2.5-4.0 L/min/m2",
                "abnormal_criteria": (
                    "CI < 2.2 L/min/m2: low output state. "
                    "CI < 1.8 L/min/m2: cardiogenic shock. "
                    "Assumed VO2 (125 mL/min/m2) may overestimate CO in "
                    "obese or critically ill patients; direct VO2 preferred."
                ),
                "clinical_significance": (
                    "Gold standard cardiac output measurement. More accurate "
                    "than thermodilution in low output states, severe "
                    "tricuspid regurgitation, and intracardiac shunts. "
                    "Required when shunt quantification (Qp/Qs) is needed."
                ),
                "calculation_formula": (
                    "CO = VO2 / (CaO2 - CvO2). "
                    "CaO2 = (Hb * 1.36 * SaO2) + (0.003 * PaO2). "
                    "CI = CO / BSA."
                ),
                "cathlab_context": "right_heart_cath",
            },
            {
                "parameter_name": "Qp/Qs (Shunt Ratio)",
                "measurement_method": (
                    "Oxygen saturation sampling from SVC, IVC, PA, and "
                    "pulmonary veins (or assumed values) during right and "
                    "left heart catheterisation"
                ),
                "normal_range": "Qp/Qs = 1.0 (no shunt)",
                "abnormal_criteria": (
                    "Qp/Qs > 1.5: significant left-to-right shunt, consider "
                    "closure (ASD, VSD, PFO). Qp/Qs > 2.0: large shunt with "
                    "RV volume overload. Qp/Qs < 1.0: right-to-left shunt "
                    "(Eisenmenger physiology if < 0.7)."
                ),
                "clinical_significance": (
                    "Quantifies intracardiac shunt severity. ASD/VSD closure "
                    "threshold is Qp/Qs >= 1.5 with symptoms or RV dilation. "
                    "Eisenmenger syndrome (Qp/Qs < 1.0 with cyanosis) is a "
                    "contraindication to shunt closure."
                ),
                "calculation_formula": (
                    "Qp/Qs = (SaO2 - SvO2) / (SpvO2 - SpaO2). "
                    "SvO2 = mixed venous O2 sat (PA). SpvO2 = pulmonary "
                    "venous O2 sat (assumed 98% if not sampled). "
                    "SpaO2 = pulmonary artery O2 sat."
                ),
                "cathlab_context": "diagnostic_cath",
            },
        ]

    @staticmethod
    def _protocol_data() -> List[dict]:
        """Catheterisation lab protocol reference data.

        Returns:
            List of cath lab protocol dictionaries.
        """
        return [
            {
                "parameter_name": "Right Heart Catheterisation Protocol",
                "measurement_method": "Swan-Ganz catheter via IJ/femoral/subclavian access",
                "normal_range": "N/A — procedural protocol",
                "abnormal_criteria": (
                    "Sequential measurement: RA, RV, PA, PCWP. "
                    "Perform thermodilution CO (or Fick if shunt suspected). "
                    "Calculate PVR, SVR, CI, transpulmonary gradient."
                ),
                "clinical_significance": (
                    "RHC is required for: (1) Heart transplant evaluation "
                    "(assess PVR, TPG, DPG), (2) WHO Group 1 PAH diagnosis, "
                    "(3) Hemodynamic-guided HF management, (4) Constrictive "
                    "vs restrictive physiology differentiation."
                ),
                "calculation_formula": (
                    "TPG = mPAP - PCWP; DPG = dPAP - PCWP; "
                    "PVR = TPG / CO; SVR = (MAP - CVP) / CO * 80"
                ),
                "cathlab_context": "right_heart_cath",
            },
            {
                "parameter_name": "Coronary Angiography Protocol",
                "measurement_method": "Selective catheterisation via radial or femoral access",
                "normal_range": "N/A — procedural protocol",
                "abnormal_criteria": (
                    "Stenosis grading: <50% mild, 50-69% moderate, "
                    "70-99% severe (>50% for left main). "
                    "TIMI flow grading: 0 (no flow) to 3 (normal). "
                    "Coronary dominance: right (85%), left (8%), co-dominant (7%)."
                ),
                "clinical_significance": (
                    "Coronary angiography remains the gold standard for "
                    "coronary anatomy definition. Guides revascularization "
                    "strategy: PCI vs CABG based on SYNTAX score, number "
                    "of vessels, LM involvement, and LV function."
                ),
                "calculation_formula": (
                    "SYNTAX score: quantitative coronary angiographic scoring "
                    "system. Low (0-22): PCI reasonable. Intermediate (23-32): "
                    "Heart team. High (>=33): CABG preferred."
                ),
                "cathlab_context": "diagnostic_cath",
            },
            {
                "parameter_name": "Endomyocardial Biopsy Protocol",
                "measurement_method": "RV septum biopsy via IJ or femoral venous access",
                "normal_range": "N/A — procedural protocol",
                "abnormal_criteria": (
                    "Indications: unexplained new-onset HF < 2 weeks with "
                    "hemodynamic compromise; unexplained HF 2 weeks to 3 months "
                    "with dilated LV unresponsive to therapy; suspected "
                    "cardiac amyloidosis (if non-invasive testing inconclusive); "
                    "post-transplant rejection surveillance."
                ),
                "clinical_significance": (
                    "EMB provides histological diagnosis for infiltrative "
                    "cardiomyopathies (amyloid, sarcoid, iron), myocarditis "
                    "(Dallas criteria, immunohistochemistry), and transplant "
                    "rejection (ISHLT grading). Complication rate ~1% "
                    "(perforation, arrhythmia)."
                ),
                "calculation_formula": "N/A — histopathological analysis",
                "cathlab_context": "diagnostic_cath",
            },
            {
                "parameter_name": "Pre-Transplant Hemodynamic Evaluation",
                "measurement_method": "Comprehensive RHC with vasoreactivity testing",
                "normal_range": "N/A — transplant candidacy assessment",
                "abnormal_criteria": (
                    "Transplant contraindication thresholds: "
                    "PVR > 5 WU (or > 3 WU after vasodilator challenge); "
                    "TPG > 15 mmHg; DPG > 7 mmHg. "
                    "Reversibility with iNO/milrinone/nitroprusside "
                    "required to qualify for transplant listing."
                ),
                "clinical_significance": (
                    "Pre-transplant RHC is mandatory (ISHLT listing criteria). "
                    "Elevated, fixed PVR causes acute RV failure of the donor "
                    "heart post-transplant. If PVR is elevated but reversible, "
                    "LVAD as bridge-to-transplant may lower PVR over time."
                ),
                "calculation_formula": (
                    "PVR = (mPAP - PCWP) / CO. "
                    "Vasoreactivity: repeat after iNO 40 ppm or IV "
                    "nitroprusside. Target PVR < 3 WU."
                ),
                "cathlab_context": "pre_transplant",
            },
            {
                "parameter_name": "PCI Technique Details",
                "measurement_method": "Interventional catheterisation via radial or femoral access",
                "normal_range": "N/A — procedural protocol",
                "abnormal_criteria": (
                    "Atherectomy indications: severely calcified lesions "
                    "(rotational, orbital, or laser atherectomy). "
                    "Thrombectomy: large thrombus burden in STEMI (aspiration "
                    "or rheolytic). Bifurcation strategies: provisional "
                    "stenting (default), two-stent (DK-crush, culotte, "
                    "TAP) for true bifurcations with large side branch. "
                    "CTO approach: antegrade wire escalation, antegrade "
                    "dissection re-entry (ADR), retrograde techniques."
                ),
                "clinical_significance": (
                    "Technique selection impacts procedural success and "
                    "long-term outcomes. Intravascular imaging (IVUS/OCT) "
                    "recommended for left main PCI, CTO, and stent "
                    "optimisation. Physiology-guided complete "
                    "revascularization (COMPLETE trial) in STEMI with "
                    "multivessel disease."
                ),
                "calculation_formula": (
                    "Lesion preparation: balloon-to-artery ratio 1:1. "
                    "Stent sizing: IVUS MLA-guided or OCT lumen-based. "
                    "SYNTAX score for revascularization strategy."
                ),
                "cathlab_context": "interventional",
            },
            {
                "parameter_name": "Impella / Mechanical Circulatory Support Management",
                "measurement_method": (
                    "Percutaneous or surgical insertion of Impella "
                    "(2.5, CP, 5.0, 5.5, RP) or other MCS devices"
                ),
                "normal_range": "N/A — procedural and management protocol",
                "abnormal_criteria": (
                    "Insertion: femoral arterial access (14Fr for CP, "
                    "21Fr for 5.0), positioned across aortic valve into LV. "
                    "P-level settings: P2 (low support ~1.5 L/min) to "
                    "P8 (max support ~3.5-5.0 L/min depending on device). "
                    "Weaning protocol: reduce P-level by 2 every 6-12 hours "
                    "if MAP >60, CI >2.2, PCWP <18, lactate normalising. "
                    "Hemolysis monitoring: plasma-free Hgb q6h (target <40 mg/dL), "
                    "LDH, haptoglobin, urine colour."
                ),
                "clinical_significance": (
                    "Impella provides active LV unloading, reducing myocardial "
                    "oxygen demand and PCWP. Indicated for cardiogenic shock "
                    "(SCAI stage C-E), high-risk PCI support, and acute MI "
                    "complicated by shock. Complications: hemolysis, limb "
                    "ischemia, device migration, aortic regurgitation. "
                    "Impella RP for isolated RV failure post-LVAD or post-MI."
                ),
                "calculation_formula": (
                    "Cardiac power output (CPO) = MAP * CO / 451. "
                    "CPO < 0.6 W: cardiogenic shock. "
                    "Monitor placement: inlet area 3.5 cm below aortic valve."
                ),
                "cathlab_context": "interventional",
            },
        ]
