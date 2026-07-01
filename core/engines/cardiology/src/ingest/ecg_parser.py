"""ECG and electrophysiology data ingest parser for Cardiology Intelligence Agent.

Provides curated reference data for ECG interpretation criteria and
arrhythmia management protocols, targeting the ``cardio_electrophysiology``
Milvus collection.

Covers:
  - Standard 12-lead ECG diagnostic criteria (LVH, STEMI, bundle branch
    blocks, pre-excitation, Brugada, long QT, etc.)
  - Arrhythmia management pathways (AF, VT, SVT, bradycardia, ACLS)

Author: Adam Jones
Date: March 2026
"""

from typing import List

from .base import BaseIngestParser, IngestRecord


# ═══════════════════════════════════════════════════════════════════════
# ECG PARSER
# ═══════════════════════════════════════════════════════════════════════


class ECGParser(BaseIngestParser):
    """Ingest parser for ECG / electrophysiology reference data.

    Seeds two sub-categories into the ``cardio_electrophysiology``
    Milvus collection:
      1. ECG diagnostic criteria (:meth:`seed_ecg_criteria`)
      2. Arrhythmia management protocols (:meth:`seed_arrhythmia_management`)

    Usage::

        parser = ECGParser()
        all_records = parser.run()

        # Or individually:
        ecg = parser.seed_ecg_criteria()
        arrhythmia = parser.seed_arrhythmia_management()
    """

    def __init__(self):
        """Initialize the ECG parser."""
        super().__init__("cardio_electrophysiology")

    # ── Fetch / Parse interface ───────────────────────────────────────

    def fetch(self, **kwargs) -> List[dict]:
        """Return all ECG/EP seed data as raw dictionaries.

        Returns:
            Combined list of ECG criteria and arrhythmia management data.
        """
        all_data: List[dict] = []
        all_data.extend(self._ecg_criteria_data())
        all_data.extend(self._arrhythmia_data())
        self.logger.info(f"Loaded {len(all_data)} ECG/EP reference entries")
        return all_data

    def parse(self, raw_data: List[dict]) -> List[IngestRecord]:
        """Convert ECG/EP dictionaries into IngestRecord instances.

        Args:
            raw_data: List of ECG/EP dictionaries from :meth:`fetch`.

        Returns:
            List of :class:`IngestRecord` instances.
        """
        records: List[IngestRecord] = []
        for entry in raw_data:
            text = (
                f"{entry.get('category', '')} — "
                f"{entry.get('finding', '')}. "
                f"Criteria: {entry.get('criteria', '')}. "
                f"Significance: {entry.get('clinical_significance', '')}. "
                f"Management: {entry.get('management', '')}."
            )

            metadata = {
                "category": self.truncate(entry.get("category", ""), 62),
                "finding": self.truncate(entry.get("finding", ""), 254),
                "criteria": self.truncate(entry.get("criteria", ""), 1022),
                "clinical_significance": self.truncate(
                    entry.get("clinical_significance", ""), 1022
                ),
                "urgency": self.truncate(entry.get("urgency", "routine"), 30),
                "associated_conditions": self.truncate(
                    entry.get("associated_conditions", ""), 510
                ),
                "management": self.truncate(
                    entry.get("management", ""), 1022
                ),
            }

            records.append(
                IngestRecord(
                    text=self.truncate(text, 2048),
                    metadata=metadata,
                    collection=self.collection,
                    source="ECG/EP Reference",
                )
            )
        return records

    # ── Convenience seed methods ──────────────────────────────────────

    def seed_ecg_criteria(self) -> List[IngestRecord]:
        """Seed ECG diagnostic criteria.

        Returns:
            List of :class:`IngestRecord` for ECG criteria.
        """
        return self.filter_valid(self.parse(self._ecg_criteria_data()))

    def seed_arrhythmia_management(self) -> List[IngestRecord]:
        """Seed arrhythmia management protocols.

        Returns:
            List of :class:`IngestRecord` for arrhythmia management.
        """
        return self.filter_valid(self.parse(self._arrhythmia_data()))

    # ── Private data methods ──────────────────────────────────────────

    @staticmethod
    def _ecg_criteria_data() -> List[dict]:
        """Standard 12-lead ECG diagnostic criteria.

        Returns:
            List of ECG criteria dictionaries.
        """
        return [
            {
                "category": "ECG",
                "finding": "ST-Elevation Myocardial Infarction (STEMI)",
                "criteria": (
                    "New ST elevation at the J point in >=2 contiguous leads: "
                    ">=1 mm in all leads except V2-V3 where >=2 mm in men >=40y, "
                    ">=2.5 mm in men <40y, >=1.5 mm in women. "
                    "New LBBB with clinical suspicion of MI (Sgarbossa criteria)."
                ),
                "clinical_significance": (
                    "STEMI requires emergent reperfusion (primary PCI within "
                    "90 minutes door-to-balloon or fibrinolysis within 30 "
                    "minutes if PCI unavailable within 120 minutes)."
                ),
                "urgency": "emergent",
                "associated_conditions": (
                    "Acute coronary syndrome, coronary artery occlusion"
                ),
                "management": (
                    "Activate catheterization lab; aspirin 325 mg + P2Y12 "
                    "inhibitor; heparin; emergent PCI; consider GPIIb/IIIa "
                    "inhibitor if high thrombus burden."
                ),
            },
            {
                "category": "ECG",
                "finding": "Left Ventricular Hypertrophy (LVH)",
                "criteria": (
                    "Sokolow-Lyon: S(V1) + R(V5 or V6) >= 35 mm. "
                    "Cornell voltage: R(aVL) + S(V3) >= 28 mm (men), "
                    ">= 20 mm (women). Cornell product: (voltage) x QRS "
                    "duration >= 2440 mm*ms."
                ),
                "clinical_significance": (
                    "LVH on ECG indicates increased LV mass, commonly due "
                    "to hypertension, aortic stenosis, or hypertrophic "
                    "cardiomyopathy. Associated with increased risk of HF, "
                    "arrhythmias, and cardiovascular death."
                ),
                "urgency": "routine",
                "associated_conditions": (
                    "Hypertension, aortic stenosis, HCM, athlete's heart"
                ),
                "management": (
                    "Echocardiography for structural confirmation; treat "
                    "underlying cause (BP control, valve intervention); "
                    "consider CMR if HCM suspected."
                ),
            },
            {
                "category": "ECG",
                "finding": "Left Bundle Branch Block (LBBB)",
                "criteria": (
                    "QRS >= 120 ms; broad monophasic R in I, aVL, V5-V6; "
                    "absent Q waves in I, V5-V6; rS or QS in V1-V2; "
                    "ST-T changes discordant to QRS."
                ),
                "clinical_significance": (
                    "New LBBB may indicate acute MI (use Sgarbossa criteria). "
                    "Chronic LBBB with HFrEF is a CRT indication (QRS >=150 ms). "
                    "LBBB can cause dyssynchrony worsening HF."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "CAD, dilated cardiomyopathy, hypertension, "
                    "conduction system disease"
                ),
                "management": (
                    "If new: evaluate for acute MI. If chronic with HFrEF: "
                    "evaluate for CRT. Echocardiography for structural "
                    "assessment."
                ),
            },
            {
                "category": "ECG",
                "finding": "Brugada Pattern (Type 1)",
                "criteria": (
                    "Coved ST elevation >= 2 mm in >=1 right precordial lead "
                    "(V1-V2) followed by negative T wave. Type 1 pattern is "
                    "the only diagnostic pattern (saddleback = Type 2, "
                    "suggestive only)."
                ),
                "clinical_significance": (
                    "Brugada syndrome is a channelopathy causing SCD from "
                    "polymorphic VT/VF. Type 1 pattern + syncope or "
                    "documented arrhythmia = high risk."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "SCN5A mutation, familial SCD, unexplained syncope"
                ),
                "management": (
                    "ICD for secondary prevention (survived arrest) or high-risk "
                    "primary prevention. Avoid Class IA/IC antiarrhythmics, "
                    "fever (treat aggressively). Genetic testing for SCN5A. "
                    "EP study for risk stratification."
                ),
            },
            {
                "category": "ECG",
                "finding": "Long QT Syndrome",
                "criteria": (
                    "QTc > 480 ms (Bazett correction) on repeated ECG. "
                    "QTc > 500 ms significantly increases TdP risk. "
                    "Schwartz score >= 3.5 = high probability."
                ),
                "clinical_significance": (
                    "Congenital LQTS (LQT1, LQT2, LQT3) or acquired "
                    "(drug-induced) QT prolongation increases risk of torsades "
                    "de pointes and SCD. LQT1: exercise-triggered; "
                    "LQT2: auditory/emotional; LQT3: sleep."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "KCNQ1/KCNH2/SCN5A mutations, drug-induced "
                    "(sotalol, dofetilide, QT-prolonging drugs), "
                    "electrolyte imbalance"
                ),
                "management": (
                    "Beta-blockers (nadolol preferred) for congenital LQTS. "
                    "Discontinue offending drugs. IV magnesium for TdP. "
                    "ICD for high risk or cardiac arrest survivors. "
                    "Genetic testing and family screening."
                ),
            },
            {
                "category": "ECG",
                "finding": "Wolff-Parkinson-White (WPW) Pattern",
                "criteria": (
                    "Short PR interval (< 120 ms); delta wave (slurred QRS "
                    "upstroke); wide QRS (> 120 ms). PR shortening + delta "
                    "wave = ventricular pre-excitation."
                ),
                "clinical_significance": (
                    "WPW indicates accessory pathway (bypass tract) capable "
                    "of re-entrant tachycardia (AVRT). Risk of rapid "
                    "pre-excited AF degenerating to VF (avoid AV nodal "
                    "blockers in pre-excited AF)."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "AVRT, pre-excited atrial fibrillation, Ebstein anomaly"
                ),
                "management": (
                    "Catheter ablation of accessory pathway (curative, >95% "
                    "success). Procainamide or ibutilide for pre-excited AF "
                    "(NEVER adenosine, verapamil, or digoxin). EP study for "
                    "risk stratification."
                ),
            },
            {
                "category": "ECG",
                "finding": "Complete Heart Block (Third-Degree AV Block)",
                "criteria": (
                    "Complete AV dissociation: no relationship between P waves "
                    "and QRS complexes. Atrial rate > ventricular rate. "
                    "Escape rhythm: narrow = junctional (40-60 bpm); "
                    "wide = ventricular (20-40 bpm)."
                ),
                "clinical_significance": (
                    "Complete heart block causes hemodynamic compromise, "
                    "syncope, and risk of asystole. Requires urgent "
                    "temporary pacing if symptomatic, followed by "
                    "permanent pacemaker."
                ),
                "urgency": "emergent",
                "associated_conditions": (
                    "Inferior MI, conduction system disease, post-TAVR, "
                    "Lyme disease, drug toxicity (digoxin, beta-blockers)"
                ),
                "management": (
                    "Atropine 0.5 mg IV (may help if junctional escape). "
                    "Temporary transcutaneous or transvenous pacing. "
                    "Permanent pacemaker implantation. Treat reversible "
                    "causes (Lyme, drug toxicity)."
                ),
            },
            {
                "category": "ECG",
                "finding": "Early Repolarization Pattern",
                "criteria": (
                    "J-point elevation >=1 mm in >=2 contiguous inferior "
                    "(II, III, aVF) or lateral (I, aVL, V4-V6) leads. "
                    "J-point elevation may manifest as QRS slurring or notching."
                ),
                "clinical_significance": (
                    "Early repolarization in inferior leads is associated with "
                    "idiopathic ventricular fibrillation and increased risk of "
                    "arrhythmic death. Previously considered benign; Haissaguerre "
                    "et al. demonstrated association with VF."
                ),
                "urgency": "routine",
                "associated_conditions": (
                    "Idiopathic VF, J-wave syndromes, athlete's heart"
                ),
                "management": (
                    "Risk stratify based on symptoms and family history. "
                    "ICD for survivors of cardiac arrest. Isoproterenol for "
                    "VF storm. Quinidine as adjunctive therapy. Avoid vagotonic "
                    "states and excessive bradycardia."
                ),
            },
            {
                "category": "ECG",
                "finding": "Right Ventricular Hypertrophy (RVH)",
                "criteria": (
                    "R in V1 >7 mm; R/S ratio >1 in V1; right axis deviation "
                    ">110 degrees; persistent S waves in V5-V6; may have "
                    "right atrial enlargement (P pulmonale)."
                ),
                "clinical_significance": (
                    "RVH on ECG indicates increased RV mass due to pressure "
                    "or volume overload. Common in pulmonary hypertension, "
                    "pulmonary embolism, congenital heart disease, and COPD."
                ),
                "urgency": "routine",
                "associated_conditions": (
                    "Pulmonary hypertension, pulmonary embolism, congenital "
                    "heart disease, COPD, mitral stenosis"
                ),
                "management": (
                    "Echocardiography for RV assessment and PASP estimation. "
                    "Evaluate for underlying cause: CTPA if PE suspected, "
                    "PFTs for lung disease, right heart catheterization for "
                    "PH confirmation. Treat underlying etiology."
                ),
            },
            {
                "category": "ECG",
                "finding": "Epsilon Waves",
                "criteria": (
                    "Small positive deflections (low amplitude signals) at the "
                    "end of QRS complex in leads V1-V3, best seen with "
                    "increased gain and reduced sweep speed. Represents delayed "
                    "RV activation through fibrofatty tissue."
                ),
                "clinical_significance": (
                    "Epsilon waves are pathognomonic for arrhythmogenic right "
                    "ventricular cardiomyopathy (ARVC). A major diagnostic "
                    "criterion in the 2010 Task Force Criteria for ARVC."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "ARVC, desmosomal gene mutations (PKP2, DSP, DSG2), "
                    "familial sudden cardiac death"
                ),
                "management": (
                    "CMR for RV assessment and fibrofatty infiltration. "
                    "Genetic testing for desmosomal mutations. ICD for SCD "
                    "prevention if high risk. Exercise restriction. Family "
                    "screening with ECG, echo, and CMR."
                ),
            },
            {
                "category": "ECG",
                "finding": "Fragmented QRS (fQRS)",
                "criteria": (
                    "Additional R wave (R') or notching of R or S wave in "
                    ">=2 contiguous leads corresponding to a myocardial "
                    "territory. Narrow fQRS (<120 ms) or wide fQRS (>=120 ms). "
                    "Not due to bundle branch block."
                ),
                "clinical_significance": (
                    "Fragmented QRS is associated with myocardial scar from "
                    "prior MI, non-ischemic cardiomyopathy, or infiltrative "
                    "disease. Predicts sudden cardiac death, arrhythmic events, "
                    "and all-cause mortality in multiple populations."
                ),
                "urgency": "routine",
                "associated_conditions": (
                    "Prior MI, myocardial scar, non-ischemic cardiomyopathy, "
                    "sarcoidosis, ARVC, Brugada syndrome"
                ),
                "management": (
                    "Echocardiography and CMR for scar assessment. Evaluate "
                    "for CAD if not previously known. Consider ICD evaluation "
                    "if associated with reduced LVEF. Monitor for arrhythmias."
                ),
            },
            {
                "category": "ECG",
                "finding": "T-Wave Alternans",
                "criteria": (
                    "Beat-to-beat alternation of T-wave amplitude, morphology, "
                    "or polarity. Can be macroscopic (visible on surface ECG) "
                    "or microvolt (detected by specialized analysis). "
                    "Typically assessed during exercise or pacing at HR 100-110."
                ),
                "clinical_significance": (
                    "T-wave alternans reflects repolarization heterogeneity "
                    "and predicts ventricular tachycardia and ventricular "
                    "fibrillation. Microvolt TWA testing aids SCD risk "
                    "stratification in patients with reduced LVEF."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "Ischemic cardiomyopathy, non-ischemic cardiomyopathy, "
                    "long QT syndrome, heart failure, SCD risk"
                ),
                "management": (
                    "Evaluate LVEF and underlying substrate. Consider ICD "
                    "implantation if positive TWA with reduced LVEF. "
                    "Optimise HF therapy (beta-blockers reduce TWA). "
                    "Correct electrolyte abnormalities."
                ),
            },
            {
                "category": "ECG",
                "finding": "Wellens Syndrome",
                "criteria": (
                    "Type A (25%): biphasic T waves (positive-negative) in "
                    "V2-V3. Type B (75%): deeply and symmetrically inverted "
                    "T waves in V2-V3, may extend V1-V6. Occurs in pain-free "
                    "interval. Minimal or no ST elevation. Normal or minimally "
                    "elevated troponin. Preserved R-wave progression."
                ),
                "clinical_significance": (
                    "Wellens syndrome indicates critical proximal LAD stenosis "
                    "(typically >90%). High risk of anterior wall MI if not "
                    "recognized and treated. Stress testing is contraindicated."
                ),
                "urgency": "emergent",
                "associated_conditions": (
                    "Critical LAD stenosis, acute coronary syndrome, "
                    "anterior wall MI"
                ),
                "management": (
                    "Admit to monitored bed. Anticoagulation and antiplatelet "
                    "therapy. Urgent coronary angiography with intent to "
                    "revascularize (PCI or CABG). Stress testing is "
                    "CONTRAINDICATED — risk of precipitating anterior STEMI."
                ),
            },
            {
                "category": "ECG",
                "finding": "de Winter Pattern",
                "criteria": (
                    "Upsloping ST-segment depression >1 mm at the J-point in "
                    "leads V1-V6 with tall, prominent, symmetric T waves. "
                    "Absence of ST elevation in precordial leads. May have "
                    "1-3 mm ST elevation in aVR."
                ),
                "clinical_significance": (
                    "The de Winter pattern is a STEMI equivalent representing "
                    "acute LAD occlusion. Present in approximately 2% of "
                    "anterior MI presentations. Often not recognized by "
                    "automated ECG algorithms."
                ),
                "urgency": "emergent",
                "associated_conditions": (
                    "Acute LAD occlusion, acute coronary syndrome, "
                    "STEMI equivalent"
                ),
                "management": (
                    "Treat as STEMI equivalent — activate catheterization lab "
                    "for emergent PCI. Aspirin 325 mg + P2Y12 inhibitor; "
                    "anticoagulation with heparin. Do not wait for serial "
                    "ECGs to show ST elevation."
                ),
            },
            {
                "category": "ECG",
                "finding": "Sgarbossa Criteria",
                "criteria": (
                    "STEMI diagnosis in setting of LBBB: (1) Concordant ST "
                    "elevation >1 mm in leads with positive QRS (5 points); "
                    "(2) Concordant ST depression >1 mm in V1-V3 (3 points); "
                    "(3) Discordant ST elevation >5 mm in leads with negative "
                    "QRS (2 points). Score >=3 = STEMI. Modified Sgarbossa: "
                    "discordant ST/S ratio >0.25 replaces criterion 3."
                ),
                "clinical_significance": (
                    "Sgarbossa criteria enable STEMI diagnosis in the "
                    "presence of LBBB, which normally obscures MI diagnosis. "
                    "Modified Sgarbossa (Smith criteria) improves sensitivity "
                    "from 36% to 91% while maintaining high specificity."
                ),
                "urgency": "emergent",
                "associated_conditions": (
                    "STEMI with LBBB, acute coronary syndrome, "
                    "ventricular paced rhythm"
                ),
                "management": (
                    "If Sgarbossa or modified Sgarbossa positive: activate "
                    "catheterization lab for emergent PCI. Standard STEMI "
                    "protocol: aspirin, P2Y12 inhibitor, heparin. "
                    "Applies also to ventricular paced rhythms."
                ),
            },
        ]

    @staticmethod
    def _arrhythmia_data() -> List[dict]:
        """Arrhythmia management protocols.

        Returns:
            List of arrhythmia management dictionaries.
        """
        return [
            {
                "category": "EP",
                "finding": "Atrial Fibrillation — Rate vs Rhythm Control",
                "criteria": (
                    "Diagnosis: irregularly irregular rhythm, absent P waves. "
                    "Paroxysmal (<7 days), persistent (>7 days), long-standing "
                    "persistent (>12 months), permanent."
                ),
                "clinical_significance": (
                    "AF increases stroke risk 5-fold. CHA2DS2-VASc score "
                    "guides anticoagulation. EAST-AFNET 4 supports early "
                    "rhythm control in recently diagnosed AF."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "Stroke/TIA, heart failure, hypertension, valvular disease, "
                    "thyroid disease, sleep apnea"
                ),
                "management": (
                    "Anticoagulation per CHA2DS2-VASc (DOACs preferred). "
                    "Rate control: beta-blocker or CCB. Rhythm control: "
                    "flecainide, amiodarone, or catheter ablation (PVI). "
                    "Manage modifiable risk factors (weight, BP, OSA, alcohol)."
                ),
            },
            {
                "category": "EP",
                "finding": "Ventricular Tachycardia (VT) — Sustained",
                "criteria": (
                    "Wide complex tachycardia >= 120 ms, rate >= 100 bpm, "
                    "sustained >= 30 seconds or requiring intervention. "
                    "AV dissociation, fusion/capture beats, concordance "
                    "= VT. Brugada criteria for differentiation from SVT."
                ),
                "clinical_significance": (
                    "Sustained VT is a medical emergency with risk of "
                    "degeneration to VF and cardiac arrest. Monomorphic VT "
                    "in structural heart disease is often scar-mediated."
                ),
                "urgency": "emergent",
                "associated_conditions": (
                    "Ischemic cardiomyopathy, non-ischemic cardiomyopathy, "
                    "ARVC, sarcoidosis, channelopathies"
                ),
                "management": (
                    "Hemodynamically unstable: synchronised cardioversion. "
                    "Stable: IV amiodarone or procainamide. Long-term: ICD "
                    "implantation, catheter ablation for recurrent VT, "
                    "optimise HF therapy. Evaluate for ischemia."
                ),
            },
            {
                "category": "EP",
                "finding": "Supraventricular Tachycardia (SVT)",
                "criteria": (
                    "Narrow complex tachycardia (QRS < 120 ms), regular, "
                    "rate typically 140-250 bpm. Types: AVNRT (most common), "
                    "AVRT (accessory pathway), atrial tachycardia."
                ),
                "clinical_significance": (
                    "SVT is usually not life-threatening but causes "
                    "palpitations, dizziness, presyncope. Recurrent SVT "
                    "may be curable by catheter ablation."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "WPW syndrome, structural heart disease (atrial "
                    "tachycardia), thyroid disease, caffeine/stimulants"
                ),
                "management": (
                    "Acute: vagal manoeuvres (modified Valsalva), adenosine "
                    "6-12 mg IV rapid push. If refractory: IV verapamil or "
                    "beta-blocker. Catheter ablation for recurrent SVT "
                    "(>95% cure rate for AVNRT)."
                ),
            },
            {
                "category": "EP",
                "finding": "Torsades de Pointes (TdP)",
                "criteria": (
                    "Polymorphic VT with twisting QRS axis around baseline, "
                    "occurring in setting of prolonged QTc (>500 ms). "
                    "'Short-long-short' initiation sequence."
                ),
                "clinical_significance": (
                    "TdP may degenerate to VF and cardiac arrest. Usually "
                    "drug-induced or associated with electrolyte abnormalities. "
                    "Distinct management from other VT forms."
                ),
                "urgency": "emergent",
                "associated_conditions": (
                    "Drug-induced QT prolongation, hypokalemia, "
                    "hypomagnesemia, congenital long QT syndrome, bradycardia"
                ),
                "management": (
                    "IV magnesium sulfate 2g bolus (first line). "
                    "Isoproterenol or temporary overdrive pacing to increase "
                    "heart rate and shorten QT. Discontinue all QT-prolonging "
                    "drugs. Correct K+ to >4.0 mEq/L. Defibrillation if "
                    "degenerates to VF."
                ),
            },
            {
                "category": "EP",
                "finding": "Atrial Flutter — Typical and Atypical",
                "criteria": (
                    "Typical (isthmus-dependent): sawtooth flutter waves in "
                    "II, III, aVF at ~300 bpm with 2:1 block (ventricular "
                    "rate ~150 bpm). Counterclockwise (common) or clockwise. "
                    "Atypical: non-isthmus-dependent, variable morphology."
                ),
                "clinical_significance": (
                    "Typical atrial flutter has high success rate with "
                    "catheter ablation. Often coexists with AF. Stroke risk "
                    "and anticoagulation approach same as AF per CHA2DS2-VASc."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "Atrial fibrillation, post-cardiac surgery, COPD, "
                    "pulmonary embolism, congenital heart disease"
                ),
                "management": (
                    "Typical flutter: catheter ablation of cavotricuspid "
                    "isthmus is first-line (Class I, >95% success). Atypical: "
                    "detailed electroanatomic mapping required. "
                    "Anticoagulation same as AF per CHA2DS2-VASc score. "
                    "Rate control with beta-blockers or CCBs."
                ),
            },
            {
                "category": "EP",
                "finding": "Bradycardia — Sick Sinus Syndrome",
                "criteria": (
                    "Sinus node dysfunction: sinus bradycardia <40 bpm, "
                    "sinus pauses >3 seconds, sinoatrial exit block, "
                    "chronotropic incompetence (failure to achieve 80% "
                    "predicted max HR with exercise), tachy-brady syndrome."
                ),
                "clinical_significance": (
                    "Sick sinus syndrome is the most common indication for "
                    "permanent pacemaker implantation. Symptoms include "
                    "syncope, presyncope, fatigue, exercise intolerance, "
                    "and heart failure."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "Aging/fibrosis of sinus node, post-cardiac surgery, "
                    "drug-induced (beta-blockers, CCBs, digoxin), "
                    "infiltrative disease, sleep apnea"
                ),
                "management": (
                    "Symptomatic sinus node dysfunction: permanent pacemaker "
                    "(dual-chamber preferred). Chronotropic incompetence: "
                    "rate-adaptive pacing. Review and reduce/discontinue "
                    "offending medications (beta-blockers, CCBs, digoxin). "
                    "Atropine or isoproterenol for acute symptomatic episodes."
                ),
            },
            {
                "category": "EP",
                "finding": "AV Block Management — Graded Approach",
                "criteria": (
                    "First degree: PR >200 ms, all P waves conducted. "
                    "Second degree Mobitz I (Wenckebach): progressive PR "
                    "prolongation then dropped QRS. Mobitz II: constant PR "
                    "with sudden dropped QRS. Third degree: complete AV "
                    "dissociation."
                ),
                "clinical_significance": (
                    "First degree and Mobitz I are usually benign (AV nodal "
                    "level). Mobitz II (infra-nodal) has high risk of "
                    "progression to complete heart block. Complete heart "
                    "block requires emergent intervention."
                ),
                "urgency": "urgent",
                "associated_conditions": (
                    "Conduction system disease, inferior MI, anterior MI, "
                    "post-TAVR, Lyme disease, cardiac surgery, medications"
                ),
                "management": (
                    "First degree: observe, no treatment needed. Mobitz I "
                    "(Wenckebach): usually benign, monitor; pacemaker if "
                    "symptomatic. Mobitz II: high risk of progression — "
                    "pacemaker indicated even if asymptomatic. Complete "
                    "heart block: emergent temporary pacing followed by "
                    "permanent pacemaker."
                ),
            },
            {
                "category": "EP",
                "finding": "Premature Ventricular Complexes (PVCs)",
                "criteria": (
                    "Wide QRS (>120 ms) without preceding P wave. PVC burden "
                    "quantified on 24-hour Holter monitor. Morphology "
                    "indicates origin: LBBB pattern = RV origin, RBBB "
                    "pattern = LV origin."
                ),
                "clinical_significance": (
                    "PVC burden >10-15% is associated with PVC-induced "
                    "cardiomyopathy (reversible with successful ablation). "
                    "Frequent PVCs may impair cardiac output and cause "
                    "symptoms of heart failure."
                ),
                "urgency": "routine",
                "associated_conditions": (
                    "Idiopathic (RVOT origin most common), ischemic heart "
                    "disease, non-ischemic cardiomyopathy, electrolyte "
                    "abnormalities, stimulants"
                ),
                "management": (
                    "PVC burden >10-15%: evaluate for PVC-induced "
                    "cardiomyopathy with echocardiography for LVEF. "
                    "Catheter ablation if symptomatic or high burden with "
                    "LV dysfunction (>90% success for RVOT PVCs). Beta-blockers "
                    "or CCBs for symptom relief. Serial echo for LVEF monitoring."
                ),
            },
            {
                "category": "EP",
                "finding": "Ventricular Fibrillation (VF) — Acute Management",
                "criteria": (
                    "Chaotic, irregular ventricular activity with no "
                    "identifiable QRS complexes. No cardiac output. "
                    "Hemodynamic collapse and cardiac arrest."
                ),
                "clinical_significance": (
                    "VF is the most common initial rhythm in out-of-hospital "
                    "cardiac arrest. Survival depends on time to "
                    "defibrillation — each minute delay reduces survival by "
                    "7-10%."
                ),
                "urgency": "emergent",
                "associated_conditions": (
                    "Acute MI, ischemic cardiomyopathy, channelopathies "
                    "(Brugada, LQTS), WPW with pre-excited AF, commotio cordis, "
                    "electrolyte abnormalities"
                ),
                "management": (
                    "Immediate defibrillation (200J biphasic). CPR with "
                    "minimal interruptions. Epinephrine 1 mg IV q3-5 min. "
                    "Amiodarone 300 mg IV bolus for refractory VF. Post-arrest: "
                    "targeted temperature management 32-36 degrees C for 24h. "
                    "ICD for secondary prevention. Genetic testing if no "
                    "structural heart disease identified."
                ),
            },
            {
                "category": "EP",
                "finding": "Inappropriate Sinus Tachycardia (IST)",
                "criteria": (
                    "Resting sinus rate >100 bpm or mean 24-hour heart rate "
                    ">90 bpm without identifiable physiologic cause. Normal "
                    "P-wave morphology. Exaggerated heart rate response to "
                    "minimal exertion."
                ),
                "clinical_significance": (
                    "IST is a diagnosis of exclusion causing significant "
                    "symptoms including palpitations, dyspnea, dizziness, and "
                    "exercise intolerance. More common in young women. "
                    "Mechanism involves enhanced sinus node automaticity or "
                    "autonomic dysregulation."
                ),
                "urgency": "routine",
                "associated_conditions": (
                    "Autonomic dysfunction, postural orthostatic tachycardia "
                    "syndrome (POTS), deconditioning, anxiety"
                ),
                "management": (
                    "Exclude secondary causes: anemia, thyroid dysfunction, "
                    "pulmonary embolism, infection, hypovolemia. Ivabradine "
                    "(sinus node If channel blocker) is first-line therapy. "
                    "Beta-blockers if tolerated. Exercise reconditioning. "
                    "Sinus node modification/ablation as last resort."
                ),
            },
        ]
