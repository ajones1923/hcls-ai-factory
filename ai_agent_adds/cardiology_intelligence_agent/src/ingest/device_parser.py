"""FDA cardiac device and implantable device ingest parser.

Provides curated reference data for two device categories:
  1. FDA-cleared AI/ML cardiac diagnostic devices (software as medical
     device — SaMD) for ECG, echo, CT, and remote monitoring.
  2. Implantable cardiac devices (pacemakers, ICDs, CRT, LVADs, LAA
     occluders) with indications and evidence.

Targets the ``cardio_devices`` Milvus collection.

Author: Adam Jones
Date: March 2026
"""

from typing import Dict, List

from .base import BaseIngestParser, IngestRecord


# ═══════════════════════════════════════════════════════════════════════
# DEVICE PARSER
# ═══════════════════════════════════════════════════════════════════════


class DeviceParser(BaseIngestParser):
    """Ingest parser for cardiac devices (AI/ML and implantable).

    Seeds two sub-categories into the ``cardio_devices`` Milvus collection:
      1. FDA-cleared AI/ML cardiac devices (:meth:`seed_fda_devices`)
      2. Implantable cardiac devices (:meth:`seed_implantable_devices`)

    Usage::

        parser = DeviceParser()
        all_records = parser.run()

        # Or individually:
        fda = parser.seed_fda_devices()
        implantable = parser.seed_implantable_devices()
    """

    FDA_AI_CARDIAC_DEVICES: List[Dict[str, str]] = [
        {
            "name": "Viz.ai ECG",
            "type": "AI_diagnostic",
            "manufacturer": "Viz.ai",
            "application": "STEMI detection and notification",
            "fda_status": "510(k) cleared",
            "evidence": (
                "AI algorithm analyzes 12-lead ECG for STEMI patterns and "
                "sends direct notification to interventional cardiology team. "
                "Reduces door-to-balloon time by average 12 minutes."
            ),
            "integration": "EHR-integrated, automated ECG routing",
        },
        {
            "name": "Caption Guidance (Caption Health)",
            "type": "AI_diagnostic",
            "manufacturer": "Caption Health (GE HealthCare)",
            "application": "AI-guided echocardiography acquisition",
            "fda_status": "De Novo (DEN200028)",
            "evidence": (
                "AI provides real-time guidance to non-expert users for "
                "acquiring diagnostic-quality echocardiographic images. "
                "Enables point-of-care cardiac ultrasound by nurses. "
                "FDA-cleared for LVEF assessment."
            ),
            "integration": "Ultrasound device integration, cloud-based analysis",
        },
        {
            "name": "Eko AI (Eko Health)",
            "type": "AI_diagnostic",
            "manufacturer": "Eko Health",
            "application": "Heart murmur and AF detection via digital stethoscope",
            "fda_status": "510(k) cleared",
            "evidence": (
                "AI algorithm detects heart murmurs with sensitivity >87% "
                "and atrial fibrillation with sensitivity >99%. "
                "Deployed in >500 healthcare systems."
            ),
            "integration": "Digital stethoscope hardware + mobile app + EHR",
        },
        {
            "name": "AliveCor KardiaMobile 6L",
            "type": "AI_diagnostic",
            "manufacturer": "AliveCor",
            "application": "Personal 6-lead ECG with AF and QTc detection",
            "fda_status": "510(k) cleared",
            "evidence": (
                "6-lead personal ECG detects AF with sensitivity >97%. "
                "FDA-cleared KardiaAI detects AF, bradycardia, tachycardia, "
                "PVCs, and wide QRS. Used in Apple Heart Study."
            ),
            "integration": "Consumer device, physician portal, Apple Health",
        },
        {
            "name": "HeartFlow FFRct",
            "type": "AI_diagnostic",
            "manufacturer": "HeartFlow",
            "application": "Non-invasive FFR from coronary CTA",
            "fda_status": "De Novo cleared",
            "evidence": (
                "Computes fractional flow reserve from standard CCTA using "
                "computational fluid dynamics. PLATFORM trial: reduced "
                "unnecessary catheterization by 61%. NXT trial: per-vessel "
                "diagnostic accuracy of 86%."
            ),
            "integration": "Cloud-based analysis, CCTA DICOM upload",
        },
        {
            "name": "Cleerly (AI Coronary Analysis)",
            "type": "AI_diagnostic",
            "manufacturer": "Cleerly",
            "application": "AI quantification of coronary atherosclerosis from CCTA",
            "fda_status": "510(k) cleared",
            "evidence": (
                "Quantifies total plaque volume, low-density non-calcified "
                "plaque, stenosis severity, and pericoronary fat attenuation "
                "from standard CCTA. SCOT-HEART sub-analysis: identified "
                "high-risk plaque predicting events."
            ),
            "integration": "Cloud-based, DICOM integration, structured reporting",
        },
        {
            "name": "Ultromics EchoGo",
            "type": "AI_diagnostic",
            "manufacturer": "Ultromics",
            "application": "Automated echocardiography analysis (LVEF, strain, diastology)",
            "fda_status": "510(k) cleared",
            "evidence": (
                "Automated LVEF, GLS, and diastolic function assessment "
                "from 2D echo. Reduces echo reporting time by 50%. "
                "Validated against expert cardiologist interpretation."
            ),
            "integration": "PACS integration, vendor-neutral, structured reporting",
        },
        {
            "name": "Tempus ACS (Philips)",
            "type": "AI_diagnostic",
            "manufacturer": "Philips",
            "application": "Pre-hospital STEMI detection in ambulance setting",
            "fda_status": "510(k) cleared",
            "evidence": (
                "AI-powered 12-lead ECG analysis for pre-hospital STEMI "
                "detection with automated transmission to receiving hospital. "
                "Reduces first medical contact to device time."
            ),
            "integration": "Ambulance ECG monitors, hospital receiving systems",
        },
        {
            "name": "Biofourmis (RhythmAnalytics)",
            "type": "AI_diagnostic",
            "manufacturer": "Biofourmis",
            "application": "Continuous remote cardiac monitoring and arrhythmia detection",
            "fda_status": "510(k) cleared",
            "evidence": (
                "AI-powered RPM platform for continuous multi-parameter "
                "cardiac monitoring. Detects AF, VT, and HF decompensation "
                "signals. Used in post-discharge HF monitoring programs."
            ),
            "integration": "Wearable biosensors, cloud platform, EHR alerts",
        },
        {
            "name": "ADAS 3D (Galgo Medical)",
            "type": "AI_diagnostic",
            "manufacturer": "Galgo Medical",
            "application": "AI-powered left atrial fibrosis mapping for AF ablation",
            "fda_status": "CE marked / 510(k) cleared",
            "evidence": (
                "Automated LGE-CMR analysis to quantify left atrial fibrosis "
                "and generate patient-specific ablation targets. Utah staging "
                "of atrial fibrosis predicts ablation outcomes."
            ),
            "integration": "CMR DICOM import, EP lab mapping system export",
        },
        {
            "name": "Apple Watch ECG",
            "type": "AI_diagnostic",
            "manufacturer": "Apple",
            "application": "Consumer single-lead ECG with AF notification",
            "fda_status": "De Novo cleared",
            "evidence": (
                "Single-lead ECG app detects AF with sensitivity >98% and "
                "specificity >99%. Apple Heart Study (N=419,297): 0.52% "
                "received irregular pulse notification; PPV 84% for AF."
            ),
            "integration": "Consumer wearable, Apple Health, physician export",
        },
        {
            "name": "BioTelemetry (Philips) MCOT",
            "type": "AI_diagnostic",
            "manufacturer": "BioTelemetry (Philips)",
            "application": "AI-enhanced mobile cardiac outpatient telemetry",
            "fda_status": "510(k) cleared",
            "evidence": (
                "Continuous ECG monitoring with AI-triggered alert for "
                "arrhythmias. 30-day MCOT detects AF in 16.7% of "
                "cryptogenic stroke patients vs 3.2% with 24h Holter."
            ),
            "integration": "Wearable patch, cloud monitoring center, EHR",
        },
        {
            "name": "Anumana ECG-AI (nCode)",
            "type": "AI_diagnostic",
            "manufacturer": "Anumana (Mayo Clinic)",
            "application": "AI detection of low LVEF from standard 12-lead ECG",
            "fda_status": "510(k) cleared",
            "evidence": (
                "Deep learning algorithm identifies LVEF <=35% from standard "
                "12-lead ECG with AUC 0.93. Enables screening for "
                "asymptomatic LV dysfunction without echocardiography. "
                "Validated across multiple health systems."
            ),
            "integration": "ECG system integration, automated screening alerts",
        },
        {
            "name": "RocketAF",
            "type": "AI_diagnostic",
            "manufacturer": "Rockethealth (wearable partners)",
            "application": "AI-powered AF detection on Apple Watch and wearables",
            "fda_status": "510(k) cleared (2023)",
            "evidence": (
                "AI algorithm continuously monitors wearable PPG and "
                "accelerometer data for AF episodes. FDA-cleared for "
                "passive AF detection with high sensitivity on consumer "
                "wearables including Apple Watch."
            ),
            "integration": "Apple Watch, consumer wearables, mobile app",
        },
        {
            "name": "Cardiologs",
            "type": "AI_diagnostic",
            "manufacturer": "Cardiologs (Philips)",
            "application": "AI ECG analysis for 30+ arrhythmias",
            "fda_status": "510(k) cleared, CE marked",
            "evidence": (
                "Deep learning platform analyzes ECG recordings for over "
                "30 arrhythmia types including AF, VT, SVT, and conduction "
                "abnormalities. Validated across multi-center studies with "
                "cardiologist-level accuracy. Reduces ECG interpretation "
                "turnaround time significantly."
            ),
            "integration": "Cloud-based, Holter/event monitor integration, EHR",
        },
        {
            "name": "Us2.ai",
            "type": "AI_diagnostic",
            "manufacturer": "Us2.ai",
            "application": (
                "Automated echocardiography analysis: LVEF, GLS, "
                "diastology, valve assessment"
            ),
            "fda_status": "510(k) cleared (2023)",
            "evidence": (
                "Fully automated echocardiographic analysis providing "
                "LVEF, global longitudinal strain, diastolic function "
                "grading, and valve assessment from standard 2D echo. "
                "Validated against expert cardiologist reads with high "
                "concordance. Reduces reporting time by up to 80%."
            ),
            "integration": "PACS integration, vendor-neutral, structured reporting",
        },
        {
            "name": "Arterys Cardio AI",
            "type": "AI_diagnostic",
            "manufacturer": "Arterys (Tempus)",
            "application": "Automated CMR analysis: LV/RV volumes, function, flow",
            "fda_status": "510(k) cleared",
            "evidence": (
                "AI-powered cardiac MRI analysis providing automated "
                "LV and RV volumetric quantification, ejection fraction, "
                "and 4D flow analysis. Reduces CMR post-processing time "
                "from 30+ minutes to under 5 minutes. Validated against "
                "expert manual contouring."
            ),
            "integration": "Cloud-based, DICOM integration, PACS viewer",
        },
        {
            "name": "Viz Aorto",
            "type": "AI_diagnostic",
            "manufacturer": "Viz.ai",
            "application": "AI detection of aortic emergencies on CT",
            "fda_status": "510(k) cleared (2024)",
            "evidence": (
                "AI algorithm automatically detects aortic dissection, "
                "aneurysm, and other aortic emergencies from CT "
                "angiography. Sends real-time notification to vascular "
                "and cardiothoracic surgery teams. Reduces time to "
                "specialist notification for aortic emergencies."
            ),
            "integration": "CT scanner DICOM, automated routing, mobile notification",
        },
    ]
    """FDA-cleared AI/ML cardiac diagnostic devices."""

    def __init__(self):
        """Initialize the device parser."""
        super().__init__("cardio_devices")

    # ── Fetch / Parse interface ───────────────────────────────────────

    def fetch(self, **kwargs) -> List[dict]:
        """Return all device seed data as raw dictionaries.

        Returns:
            Combined list of FDA AI devices and implantable devices.
        """
        all_data: List[dict] = []
        all_data.extend(self._fda_device_data())
        all_data.extend(self._implantable_device_data())
        self.logger.info(f"Loaded {len(all_data)} cardiac device entries")
        return all_data

    def parse(self, raw_data: List[dict]) -> List[IngestRecord]:
        """Convert device dictionaries into IngestRecord instances.

        Args:
            raw_data: List of device dictionaries from :meth:`fetch`.

        Returns:
            List of :class:`IngestRecord` instances.
        """
        records: List[IngestRecord] = []
        for entry in raw_data:
            text = (
                f"{entry.get('name', '')} ({entry.get('manufacturer', '')}). "
                f"Type: {entry.get('type', '')}. "
                f"Application: {entry.get('application', '')}. "
                f"FDA status: {entry.get('fda_status', '')}. "
                f"Evidence: {entry.get('evidence', '')}."
            )

            metadata = {
                "device_name": self.truncate(entry.get("name", ""), 254),
                "device_type": self.truncate(entry.get("type", ""), 62),
                "manufacturer": self.truncate(
                    entry.get("manufacturer", ""), 126
                ),
                "fda_status": self.truncate(
                    entry.get("fda_status", ""), 62
                ),
                "clinical_application": self.truncate(
                    entry.get("application", ""), 510
                ),
                "evidence": self.truncate(
                    entry.get("evidence", ""), 1022
                ),
                "integration_notes": self.truncate(
                    entry.get("integration", ""), 510
                ),
            }

            records.append(
                IngestRecord(
                    text=self.truncate(text, 2048),
                    metadata=metadata,
                    collection=self.collection,
                    source="Cardiac Device Registry",
                )
            )
        return records

    # ── Convenience seed methods ──────────────────────────────────────

    def seed_fda_devices(self) -> List[IngestRecord]:
        """Seed FDA-cleared AI/ML cardiac devices.

        Returns:
            List of :class:`IngestRecord` for FDA AI devices.
        """
        return self.filter_valid(self.parse(self._fda_device_data()))

    def seed_implantable_devices(self) -> List[IngestRecord]:
        """Seed implantable cardiac devices.

        Returns:
            List of :class:`IngestRecord` for implantable devices.
        """
        return self.filter_valid(self.parse(self._implantable_device_data()))

    # ── Private data methods ──────────────────────────────────────────

    def _fda_device_data(self) -> List[dict]:
        """Return FDA AI cardiac device data.

        Returns:
            List of device dictionaries from :attr:`FDA_AI_CARDIAC_DEVICES`.
        """
        return list(self.FDA_AI_CARDIAC_DEVICES)

    @staticmethod
    def _implantable_device_data() -> List[dict]:
        """Implantable cardiac device reference data.

        Returns:
            List of implantable device dictionaries.
        """
        return [
            {
                "name": "Implantable Cardioverter-Defibrillator (ICD)",
                "type": "implantable",
                "manufacturer": "Medtronic, Abbott, Boston Scientific",
                "application": (
                    "Primary and secondary prevention of sudden cardiac death "
                    "in patients with reduced LVEF or prior cardiac arrest"
                ),
                "fda_status": "PMA approved",
                "evidence": (
                    "MADIT-II: ICD reduced all-cause mortality by 31% in "
                    "post-MI patients with LVEF <=30%. SCD-HeFT: ICD reduced "
                    "mortality by 23% in NYHA II-III HF with LVEF <=35%. "
                    "DANISH: No mortality benefit for ICD in non-ischemic CMP "
                    "(but benefit in patients <68 years)."
                ),
                "integration": "Remote monitoring (CareLink, Merlin, LATITUDE)",
            },
            {
                "name": "Cardiac Resynchronization Therapy (CRT-D / CRT-P)",
                "type": "implantable",
                "manufacturer": "Medtronic, Abbott, Boston Scientific",
                "application": (
                    "Biventricular pacing for HFrEF with electrical "
                    "dyssynchrony (LBBB, wide QRS)"
                ),
                "fda_status": "PMA approved",
                "evidence": (
                    "COMPANION: CRT-D reduced all-cause mortality by 36% in "
                    "NYHA III-IV HF with QRS >=120 ms. MADIT-CRT: CRT-D reduced "
                    "HF events by 34% in NYHA I-II HFrEF with LBBB >=130 ms. "
                    "Greatest benefit: LBBB >=150 ms."
                ),
                "integration": "Remote monitoring, device-based HF diagnostics",
            },
            {
                "name": "Permanent Pacemaker (PPM)",
                "type": "implantable",
                "manufacturer": "Medtronic, Abbott, Boston Scientific, Biotronik",
                "application": (
                    "Bradycardia therapy: sinus node dysfunction, AV block, "
                    "post-TAVR conduction disturbance"
                ),
                "fda_status": "PMA approved",
                "evidence": (
                    "Guideline-indicated for symptomatic sinus node dysfunction, "
                    "Mobitz type II and third-degree AV block. Leadless pacemakers "
                    "(Micra, AVEIR) reduce lead-related complications. "
                    "Conduction system pacing (His, LBBP) emerging as physiologic "
                    "alternative."
                ),
                "integration": "Remote monitoring, MRI-conditional designs",
            },
            {
                "name": "Left Ventricular Assist Device (LVAD) — HeartMate 3",
                "type": "implantable",
                "manufacturer": "Abbott",
                "application": (
                    "Mechanical circulatory support for advanced HF (bridge "
                    "to transplant or destination therapy)"
                ),
                "fda_status": "PMA approved",
                "evidence": (
                    "MOMENTUM 3: HeartMate 3 fully magnetically levitated "
                    "centrifugal pump was superior to HeartMate II for "
                    "survival free of disabling stroke or reoperation at "
                    "2 years (77.9% vs 56.4%). Reduced pump thrombosis "
                    "rate to near zero."
                ),
                "integration": "Patient controllers, remote monitoring, INR management",
            },
            {
                "name": "WATCHMAN FLX (LAA Occluder)",
                "type": "implantable",
                "manufacturer": "Boston Scientific",
                "application": (
                    "Left atrial appendage occlusion for stroke prevention in "
                    "non-valvular AF with contraindication to long-term "
                    "anticoagulation"
                ),
                "fda_status": "PMA approved",
                "evidence": (
                    "PROTECT AF and PREVAIL: WATCHMAN non-inferior to warfarin "
                    "for stroke prevention with significant reduction in "
                    "hemorrhagic stroke and mortality at 5 years. "
                    "CHAMPION-AF (pending): WATCHMAN FLX vs DOACs."
                ),
                "integration": "TEE/CT guidance for sizing, short-term DAPT post-implant",
            },
            {
                "name": "Implantable Loop Recorder (ILR) — LINQ II",
                "type": "implantable",
                "manufacturer": "Medtronic",
                "application": (
                    "Continuous cardiac monitoring for cryptogenic stroke, "
                    "unexplained syncope, AF detection"
                ),
                "fda_status": "510(k) cleared",
                "evidence": (
                    "CRYSTAL AF: ILR detected AF in 30% of cryptogenic stroke "
                    "patients at 3 years vs 3% with standard monitoring. "
                    "Battery life 4+ years. AI-enhanced algorithms reduce "
                    "false-positive alerts."
                ),
                "integration": "CareLink remote monitoring, auto-transmitted alerts",
            },
            {
                "name": "Subcutaneous ICD (S-ICD) — EMBLEM",
                "type": "implantable",
                "manufacturer": "Boston Scientific",
                "application": (
                    "Subcutaneous defibrillation without transvenous leads; "
                    "primary and secondary SCD prevention"
                ),
                "fda_status": "PMA approved",
                "evidence": (
                    "EMBLEM S-ICD eliminates transvenous lead complications "
                    "(endocarditis, lead fracture, venous occlusion). Requires "
                    "screening ECG to confirm appropriate sensing. No "
                    "anti-bradycardia or antitachycardia pacing capability. "
                    "Appropriate for young patients, those with limited venous "
                    "access, and prior device infections. UNTOUCHED trial: "
                    "99.2% freedom from inappropriate shocks at 18 months."
                ),
                "integration": "LATITUDE remote monitoring, MRI-conditional",
            },
            {
                "name": "Leadless Pacemaker — Micra AV",
                "type": "implantable",
                "manufacturer": "Medtronic",
                "application": (
                    "Leadless single-chamber pacemaker with AV synchrony "
                    "(VVI/VDD modes) for bradycardia"
                ),
                "fda_status": "PMA approved",
                "evidence": (
                    "Self-contained capsule implanted directly in the RV via "
                    "femoral vein catheter. No leads, no subcutaneous pocket. "
                    "Micra AV provides AV synchronous pacing (VDD) using an "
                    "accelerometer to sense atrial contraction. 63% lower rate "
                    "of major complications vs transvenous pacemakers. "
                    "Lower infection risk. Battery life ~12 years."
                ),
                "integration": "CareLink remote monitoring, MRI-conditional",
            },
            {
                "name": "Left Bundle Branch Area Pacing (LBBAP) System",
                "type": "implantable",
                "manufacturer": "Medtronic, Abbott, Boston Scientific",
                "application": (
                    "Physiologic conduction system pacing as alternative "
                    "to CRT or RV pacing"
                ),
                "fda_status": "Utilized with PMA-approved pulse generators",
                "evidence": (
                    "LBBAP captures the left bundle branch area to achieve "
                    "physiologic ventricular activation with narrow QRS. "
                    "Emerging as alternative to traditional CRT with higher "
                    "success rates than His bundle pacing. Rising adoption "
                    "for CRT-eligible patients and pacing-dependent patients. "
                    "Observational data shows comparable or superior "
                    "outcomes to BiV-CRT."
                ),
                "integration": "Standard pulse generators, remote monitoring platforms",
            },
            {
                "name": "Impella Heart Pump",
                "type": "implantable",
                "manufacturer": "Abiomed (Johnson & Johnson)",
                "application": (
                    "Percutaneous mechanical circulatory support for "
                    "cardiogenic shock and high-risk PCI"
                ),
                "fda_status": "PMA approved",
                "evidence": (
                    "Percutaneous axial-flow pump providing temporary MCS. "
                    "Models include CP (3.5 L/min), 5.0 (5.0 L/min), "
                    "5.5 (6.2 L/min), and RP (right-sided support). "
                    "Used in cardiogenic shock and protected high-risk PCI. "
                    "PROTECT III trial evaluating outcomes in AMI "
                    "cardiogenic shock. Provides direct LV unloading."
                ),
                "integration": "Automated Impella Controller, SmartAssist platform",
            },
            {
                "name": "WATCHMAN FLX Pro (LAA Occluder)",
                "type": "implantable",
                "manufacturer": "Boston Scientific",
                "application": (
                    "Next-generation left atrial appendage occlusion for "
                    "stroke prevention in non-valvular AF"
                ),
                "fda_status": "PMA approved",
                "evidence": (
                    "Next-generation LAA occluder with improved design for "
                    "higher single-procedure complete closure rates. "
                    "Enhanced conformability and dual-seal mechanism for "
                    "superior peridevice leak reduction compared to FLX. "
                    "Builds on PROTECT AF, PREVAIL, and FLX registry data. "
                    "Designed for single-procedure seal without need for "
                    "reintervention."
                ),
                "integration": "TEE/CT/ICE guidance, short-term DAPT post-implant",
            },
        ]
