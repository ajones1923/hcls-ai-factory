"""Imaging Intelligence Agent — Domain Knowledge Graph.

Extends the Clinker pattern from rag-chat-pipeline/src/knowledge.py,
adapted for the medical imaging domain. Contains:

1. IMAGING_PATHOLOGIES: ~15 pathology entries with imaging characteristics
2. IMAGING_MODALITIES: ~8 modality entries with physics and protocols
3. IMAGING_ANATOMY: ~15 anatomy entries with VISTA-3D labels and structures

Helper functions mirror the CAR-T agent pattern:
    - get_pathology_context / get_modality_context / get_anatomy_context
    - get_nim_recommendation
    - resolve_comparison_entity / get_comparison_context
    - get_knowledge_stats

Author: Adam Jones
Date: February 2026
"""

from typing import Any, Dict, List, Optional


# ═══════════════════════════════════════════════════════════════════════════════
# 1. IMAGING_PATHOLOGIES — Pathology knowledge graph (~15 entries)
# ═══════════════════════════════════════════════════════════════════════════════

IMAGING_PATHOLOGIES: Dict[str, Dict[str, Any]] = {
    "intracranial_hemorrhage": {
        "icd10": "I62.9",
        "display_name": "Intracranial Hemorrhage",
        "modalities": ["ct", "mri"],
        "body_region": "head",
        "subtypes": [
            "epidural", "subdural", "subarachnoid",
            "intraparenchymal", "intraventricular",
        ],
        "ct_characteristics": (
            "Hyperdense (acute 50-70 HU), hypodense (chronic). "
            "Blood window 0-80 HU."
        ),
        "mri_characteristics": (
            "T1 hyperintense (subacute), T2 variable, "
            "SWI blooming artifact"
        ),
        "severity_criteria": {
            "critical": "Volume > 30 mL OR midline shift > 5 mm OR thickness > 10 mm",
            "urgent": "Volume > 5 mL",
            "routine": "Volume < 5 mL, no shift",
        },
        "key_measurements": ["volume_ml", "midline_shift_mm", "max_thickness_mm"],
        "clinical_guidelines": [
            "Brain Trauma Foundation",
            "AHA/ASA Stroke Guidelines",
        ],
        "ai_models": ["3D U-Net (MONAI)", "VISTA-3D"],
        "nim_workflow": "ct_head_hemorrhage",
    },
    "lung_nodule": {
        "icd10": "R91.1",
        "display_name": "Pulmonary Nodule",
        "modalities": ["ct"],
        "body_region": "chest",
        "subtypes": ["solid", "ground_glass", "part_solid"],
        "ct_characteristics": (
            "Solid: soft tissue density. GGO: hazy opacity. "
            "Part-solid: mixed."
        ),
        "mri_characteristics": (
            "Limited role; DWI may show restricted diffusion "
            "in malignant nodules"
        ),
        "severity_criteria": {
            "critical": "Lung-RADS 4B or 4X, size > 15 mm",
            "urgent": "Lung-RADS 4A, size 8-15 mm or growing",
            "routine": "Lung-RADS 1-3, size < 6 mm",
        },
        "key_measurements": [
            "long_axis_mm", "short_axis_mm",
            "volume_mm3", "doubling_time_days",
        ],
        "clinical_guidelines": [
            "ACR Lung-RADS v2022",
            "Fleischner Society 2017",
        ],
        "ai_models": ["RetinaNet (MONAI)", "SegResNet (MONAI)", "VISTA-3D"],
        "nim_workflow": "ct_chest_lung_nodule",
    },
    "pneumonia": {
        "icd10": "J18.9",
        "display_name": "Pneumonia",
        "modalities": ["cxr", "ct"],
        "body_region": "chest",
        "subtypes": ["lobar", "bronchopneumonia", "interstitial", "aspiration"],
        "ct_characteristics": (
            "Consolidation with air bronchograms, "
            "ground-glass opacities"
        ),
        "mri_characteristics": (
            "Not typically used; T2 hyperintense consolidation "
            "if performed"
        ),
        "severity_criteria": {
            "critical": "Bilateral extensive consolidation, respiratory failure",
            "urgent": "Multilobar involvement",
            "routine": "Focal consolidation, single lobe",
        },
        "key_measurements": ["affected_lobes", "consolidation_volume_percent"],
        "clinical_guidelines": ["ATS/IDSA CAP Guidelines", "CURB-65"],
        "ai_models": ["DenseNet-121 (MONAI)"],
        "nim_workflow": "cxr_rapid_findings",
    },
    "pulmonary_embolism": {
        "icd10": "I26.99",
        "display_name": "Pulmonary Embolism",
        "modalities": ["ct"],
        "body_region": "chest",
        "subtypes": ["acute", "chronic", "saddle"],
        "ct_characteristics": (
            "Filling defect in pulmonary arteries on CTPA, "
            "RV/LV ratio > 1.0 = strain"
        ),
        "mri_characteristics": "MRA can detect PE but CT is standard",
        "severity_criteria": {
            "critical": "Saddle PE or RV strain (RV/LV > 1.0)",
            "urgent": "Segmental PE without strain",
            "routine": "Subsegmental PE",
        },
        "key_measurements": ["rv_lv_ratio", "clot_burden_score"],
        "clinical_guidelines": ["ESC 2019 PE Guidelines", "PIOPED II"],
        "ai_models": ["VISTA-3D"],
        "nim_workflow": None,
    },
    "stroke_ischemic": {
        "icd10": "I63.9",
        "display_name": "Ischemic Stroke",
        "modalities": ["ct", "mri"],
        "body_region": "brain",
        "subtypes": [
            "large_vessel_occlusion", "lacunar",
            "watershed", "embolic",
        ],
        "ct_characteristics": (
            "Early: loss of gray-white differentiation, hyperdense vessel. "
            "Late: hypodense territory"
        ),
        "mri_characteristics": (
            "DWI restriction (bright) within minutes, "
            "FLAIR hyperintense (hours)"
        ),
        "severity_criteria": {
            "critical": "Large vessel occlusion, ASPECTS < 6",
            "urgent": "ASPECTS 6-8, partial territory",
            "routine": "Lacunar infarct, small territory",
        },
        "key_measurements": [
            "aspects_score", "infarct_volume_ml", "penumbra_volume_ml",
        ],
        "clinical_guidelines": ["AHA/ASA 2019 Stroke Guidelines"],
        "ai_models": ["VISTA-3D"],
        "nim_workflow": None,
    },
    "brain_tumor": {
        "icd10": "C71.9",
        "display_name": "Brain Tumor",
        "modalities": ["mri", "ct"],
        "body_region": "brain",
        "subtypes": [
            "glioblastoma", "meningioma", "metastasis",
            "astrocytoma", "oligodendroglioma",
        ],
        "ct_characteristics": (
            "Variable density mass, surrounding edema, "
            "possible calcification"
        ),
        "mri_characteristics": (
            "T1 hypo/iso, T2/FLAIR hyperintense, variable enhancement, "
            "DWI restricted in high-grade"
        ),
        "severity_criteria": {
            "critical": "Mass effect with herniation risk, midline shift > 5 mm",
            "urgent": "Growing tumor, new enhancement",
            "routine": "Stable known tumor, surveillance",
        },
        "key_measurements": [
            "tumor_volume_ml", "enhancing_volume_ml", "midline_shift_mm",
        ],
        "clinical_guidelines": ["NCCN CNS Cancers", "RANO Criteria"],
        "ai_models": ["3D U-Net (MONAI)", "VISTA-3D", "VILA-M3"],
        "nim_workflow": None,
    },
    "ms_lesion": {
        "icd10": "G35",
        "display_name": "Multiple Sclerosis Lesion",
        "modalities": ["mri"],
        "body_region": "brain",
        "subtypes": ["new", "enlarging", "stable", "black_hole"],
        "ct_characteristics": (
            "Poorly sensitive; may show hypodense lesions "
            "in severe cases"
        ),
        "mri_characteristics": (
            "FLAIR/T2 hyperintense periventricular, juxtacortical, "
            "infratentorial, spinal. Enhancing = active."
        ),
        "severity_criteria": {
            "critical": "Tumefactive MS (>2 cm) or acute presentation",
            "urgent": "New or enlarging lesions on follow-up",
            "routine": "Stable lesion burden",
        },
        "key_measurements": [
            "lesion_count", "total_lesion_volume_ml",
            "new_lesion_count", "t1_black_hole_count",
        ],
        "clinical_guidelines": [
            "McDonald Criteria 2017",
            "MAGNIMS Guidelines",
        ],
        "ai_models": ["3D U-Net (MONAI)"],
        "nim_workflow": "mri_brain_ms_lesion",
    },
    "pneumothorax": {
        "icd10": "J93.9",
        "display_name": "Pneumothorax",
        "modalities": ["cxr", "ct"],
        "body_region": "chest",
        "subtypes": ["simple", "tension", "traumatic", "iatrogenic"],
        "ct_characteristics": (
            "Air in pleural space, lung collapse, "
            "mediastinal shift in tension"
        ),
        "mri_characteristics": "Not used for diagnosis",
        "severity_criteria": {
            "critical": "Tension pneumothorax (mediastinal shift)",
            "urgent": "Large (>2 cm at apex), symptomatic",
            "routine": "Small (<2 cm), asymptomatic",
        },
        "key_measurements": ["size_cm_at_apex", "lung_collapse_percent"],
        "clinical_guidelines": ["BTS Pleural Disease Guidelines 2010"],
        "ai_models": ["DenseNet-121 (MONAI)"],
        "nim_workflow": "cxr_rapid_findings",
    },
    "pleural_effusion": {
        "icd10": "J91.8",
        "display_name": "Pleural Effusion",
        "modalities": ["cxr", "ct", "ultrasound"],
        "body_region": "chest",
        "subtypes": ["transudative", "exudative", "hemorrhagic", "empyema"],
        "ct_characteristics": (
            "Dependent fluid density, meniscus sign, "
            "may show enhancement in empyema"
        ),
        "mri_characteristics": (
            "T2 hyperintense fluid, T1 signal depends on content"
        ),
        "severity_criteria": {
            "critical": "Massive effusion with respiratory compromise",
            "urgent": "Large effusion (>1/3 hemithorax)",
            "routine": "Small to moderate effusion",
        },
        "key_measurements": ["estimated_volume_ml", "max_depth_mm"],
        "clinical_guidelines": [
            "BTS Pleural Disease Guidelines",
            "Light's Criteria",
        ],
        "ai_models": ["DenseNet-121 (MONAI)"],
        "nim_workflow": "cxr_rapid_findings",
    },
    "fracture": {
        "icd10": "T14.8",
        "display_name": "Fracture",
        "modalities": ["xray", "ct"],
        "body_region": "extremity",
        "subtypes": [
            "simple", "comminuted", "compression",
            "stress", "pathologic",
        ],
        "ct_characteristics": (
            "Cortical disruption, fragment displacement, "
            "associated soft tissue"
        ),
        "mri_characteristics": (
            "Bone marrow edema on STIR, fracture line on T1"
        ),
        "severity_criteria": {
            "critical": (
                "Open fracture, neurovascular compromise, "
                "spinal with cord compression"
            ),
            "urgent": "Displaced fracture requiring reduction",
            "routine": "Non-displaced, stable pattern",
        },
        "key_measurements": ["displacement_mm", "angulation_degrees"],
        "clinical_guidelines": [
            "AO/OTA Classification",
            "Ottawa Ankle/Knee Rules",
        ],
        "ai_models": ["DenseNet-121 (MONAI)"],
        "nim_workflow": "cxr_rapid_findings",
    },
    "cardiomegaly": {
        "icd10": "I51.7",
        "display_name": "Cardiomegaly",
        "modalities": ["cxr"],
        "body_region": "cardiac",
        "subtypes": [
            "global", "left_ventricular",
            "right_ventricular", "biventricular",
        ],
        "ct_characteristics": (
            "Enlarged cardiac silhouette, chamber measurements on CT"
        ),
        "mri_characteristics": (
            "Cardiac MRI gold standard for volumes and function"
        ),
        "severity_criteria": {
            "critical": "Acute decompensation with pulmonary edema",
            "urgent": "New cardiomegaly, CTR > 0.6",
            "routine": "Stable known cardiomegaly",
        },
        "key_measurements": ["cardiothoracic_ratio"],
        "clinical_guidelines": ["ACC/AHA Heart Failure Guidelines"],
        "ai_models": ["DenseNet-121 (MONAI)"],
        "nim_workflow": "cxr_rapid_findings",
    },
    "aortic_dissection": {
        "icd10": "I71.0",
        "display_name": "Aortic Dissection",
        "modalities": ["ct"],
        "body_region": "chest",
        "subtypes": ["stanford_a", "stanford_b"],
        "ct_characteristics": (
            "Intimal flap, true and false lumen, "
            "may see intramural hematoma"
        ),
        "mri_characteristics": (
            "T1 dark-blood shows intimal flap, no contrast needed"
        ),
        "severity_criteria": {
            "critical": "Stanford Type A (ascending aorta involvement)",
            "urgent": "Stanford Type B with complications",
            "routine": "Uncomplicated Type B, chronic dissection",
        },
        "key_measurements": ["aortic_diameter_mm", "false_lumen_diameter_mm"],
        "clinical_guidelines": [
            "ESC 2014 Aortic Disease Guidelines",
            "Stanford Classification",
        ],
        "ai_models": ["VISTA-3D"],
        "nim_workflow": None,
    },
    "liver_lesion": {
        "icd10": "K76.9",
        "display_name": "Liver Lesion",
        "modalities": ["ct", "mri", "ultrasound"],
        "body_region": "abdomen",
        "subtypes": [
            "hemangioma", "hcc", "metastasis",
            "fnh", "adenoma", "cyst",
        ],
        "ct_characteristics": (
            "Enhancement pattern: arterial hyperenhancement (HCC), "
            "peripheral nodular (hemangioma)"
        ),
        "mri_characteristics": (
            "Hepatobiliary agents (Eovist/Primovist), "
            "DWI restriction in malignant"
        ),
        "severity_criteria": {
            "critical": "Ruptured HCC, portal vein invasion",
            "urgent": "LI-RADS 5, growing lesion",
            "routine": "LI-RADS 1-2, typical benign",
        },
        "key_measurements": ["diameter_mm", "enhancement_pattern"],
        "clinical_guidelines": [
            "ACR LI-RADS",
            "AASLD HCC Guidelines",
        ],
        "ai_models": ["VISTA-3D"],
        "nim_workflow": None,
    },
    "renal_mass": {
        "icd10": "N28.1",
        "display_name": "Renal Mass",
        "modalities": ["ct", "mri", "ultrasound"],
        "body_region": "abdomen",
        "subtypes": [
            "clear_cell_rcc", "papillary_rcc", "chromophobe",
            "oncocytoma", "angiomyolipoma", "cyst",
        ],
        "ct_characteristics": (
            "Enhancement > 20 HU = solid. Fat (-20 to -80 HU) = AML. "
            "Bosniak cyst classification."
        ),
        "mri_characteristics": (
            "Clear cell RCC: T2 hyperintense, signal drop "
            "on opposed-phase in some"
        ),
        "severity_criteria": {
            "critical": "Large mass with vena cava invasion",
            "urgent": "Growing solid enhancing mass > 3 cm",
            "routine": "Bosniak I-II cyst, small AML",
        },
        "key_measurements": ["diameter_mm", "enhancement_hu"],
        "clinical_guidelines": [
            "Bosniak Classification 2019",
            "AUA Renal Mass Guidelines",
        ],
        "ai_models": ["VISTA-3D"],
        "nim_workflow": None,
    },
    "pulmonary_fibrosis": {
        "icd10": "J84.10",
        "display_name": "Pulmonary Fibrosis",
        "modalities": ["ct"],
        "body_region": "chest",
        "subtypes": [
            "uip", "nsip",
            "organizing_pneumonia", "hypersensitivity_pneumonitis",
        ],
        "ct_characteristics": (
            "Honeycombing, traction bronchiectasis, "
            "ground-glass opacity, reticular pattern"
        ),
        "mri_characteristics": (
            "Limited role; research applications with UTE sequences"
        ),
        "severity_criteria": {
            "critical": "Acute exacerbation with new ground-glass",
            "urgent": "Progressive fibrosis on follow-up",
            "routine": "Stable UIP pattern",
        },
        "key_measurements": ["fibrosis_extent_percent", "honeycombing_present"],
        "clinical_guidelines": [
            "ATS/ERS/JRS/ALAT IPF Guidelines 2022",
        ],
        "ai_models": ["VISTA-3D"],
        "nim_workflow": None,
    },
}

# Future Phase 2 pathologies: emphysema, lymphadenopathy, breast_mass,
# thyroid_nodule, pancreatic_mass, bowel_obstruction, appendicitis,
# coronary_calcification, osteoporosis, degenerative_disc.


# ═══════════════════════════════════════════════════════════════════════════════
# 2. IMAGING_MODALITIES — Modality knowledge graph (~8 entries)
# ═══════════════════════════════════════════════════════════════════════════════

IMAGING_MODALITIES: Dict[str, Dict[str, Any]] = {
    "ct": {
        "full_name": "Computed Tomography",
        "dicom_modality_code": "CT",
        "physics": (
            "X-ray tube rotates around patient; attenuation measured by "
            "detector array; filtered back-projection or iterative "
            "reconstruction produces cross-sectional images in Hounsfield "
            "Units (HU). Air = -1000 HU, water = 0 HU, bone = +1000 HU."
        ),
        "strengths": [
            "Fast acquisition (sub-second rotation)",
            "Excellent spatial resolution (0.5-0.625 mm isotropic)",
            "Superior bone and calcification detail",
            "Widely available, 24/7 emergency use",
            "CTA/CTV for vascular mapping",
        ],
        "limitations": [
            "Ionizing radiation exposure",
            "Limited soft tissue contrast vs MRI",
            "Iodinated contrast risks (allergy, nephrotoxicity)",
            "Beam hardening and metal artifacts",
        ],
        "common_protocols": [
            "Non-contrast head CT",
            "CT angiography (CTA head/neck, chest, abdomen)",
            "High-resolution CT (HRCT) chest",
            "CT abdomen/pelvis with contrast",
            "Low-dose lung cancer screening CT",
            "CT perfusion (stroke)",
        ],
        "typical_dose_msv": {
            "head": 2.0,
            "chest": 7.0,
            "abdomen_pelvis": 10.0,
            "low_dose_lung_screening": 1.5,
            "ct_angiography": 5.0,
        },
        "nim_models": [
            "VISTA-3D (universal segmentation)",
            "3D U-Net (hemorrhage, tumor)",
            "RetinaNet (nodule detection)",
            "SegResNet (lung segmentation)",
        ],
    },
    "mri": {
        "full_name": "Magnetic Resonance Imaging",
        "dicom_modality_code": "MR",
        "physics": (
            "Strong static magnetic field (1.5T/3T) aligns hydrogen protons. "
            "RF pulses perturb alignment; relaxation signals (T1, T2, T2*) "
            "are spatially encoded via gradient coils. No ionizing radiation."
        ),
        "strengths": [
            "Superior soft tissue contrast",
            "No ionizing radiation",
            "Multiplanar and multiparametric imaging",
            "Functional imaging (fMRI, DWI, perfusion, spectroscopy)",
            "Cardiac MRI for volumes and function",
        ],
        "limitations": [
            "Long acquisition times (20-60 min per study)",
            "Motion-sensitive",
            "MRI-unsafe implants contraindicated",
            "Claustrophobia / patient tolerance",
            "Higher cost than CT",
            "Gadolinium risks (NSF in renal failure)",
        ],
        "common_protocols": [
            "Brain MRI with and without contrast",
            "MRI spine (cervical, thoracic, lumbar)",
            "Cardiac MRI (function, late gadolinium enhancement)",
            "MR angiography (MRA)",
            "Breast MRI",
            "Prostate mpMRI (PI-RADS)",
            "Liver MRI with hepatobiliary agent (Eovist)",
        ],
        "typical_dose_msv": None,  # Non-ionizing
        "nim_models": [
            "VISTA-3D (universal segmentation)",
            "3D U-Net (brain tumor, MS lesion)",
            "VILA-M3 (multimodal report generation)",
        ],
    },
    "xray": {
        "full_name": "Radiography (X-ray)",
        "dicom_modality_code": "DX",
        "physics": (
            "Single X-ray exposure through body onto digital detector. "
            "2D projection image based on differential attenuation. "
            "Bone appears white (high attenuation), air appears black."
        ),
        "strengths": [
            "Fast (seconds), low cost",
            "Widely available at all care levels",
            "Excellent for fracture detection",
            "Low radiation dose",
            "Portable units available for bedside",
        ],
        "limitations": [
            "2D projection only (superimposition of structures)",
            "Limited soft tissue detail",
            "Low sensitivity for subtle pathology",
            "Ionizing radiation (though low dose)",
        ],
        "common_protocols": [
            "PA and lateral chest",
            "Extremity (AP/lateral/oblique)",
            "Abdomen series (KUB, upright, decubitus)",
            "Spine series",
            "Pelvis AP",
        ],
        "typical_dose_msv": {
            "chest_pa": 0.02,
            "extremity": 0.001,
            "abdomen": 0.7,
            "spine_lumbar": 1.5,
            "pelvis": 0.6,
        },
        "nim_models": [
            "DenseNet-121 (multi-label classification)",
        ],
    },
    "cxr": {
        "full_name": "Chest X-Ray",
        "dicom_modality_code": "DX",
        "physics": (
            "PA and/or lateral projection radiograph of the chest. "
            "Standard 72-inch SID, upright positioning preferred. "
            "Digital detector (DR) or computed radiography (CR)."
        ),
        "strengths": [
            "First-line chest imaging",
            "Fast screening for pneumonia, effusion, pneumothorax",
            "Low dose (~0.02 mSv PA)",
            "Portable for ICU/ED patients",
            "High throughput for AI triage",
        ],
        "limitations": [
            "Low sensitivity vs CT (misses small nodules <1 cm)",
            "2D projection limits localization",
            "Reader variability in interpretation",
            "Limited mediastinal detail",
        ],
        "common_protocols": [
            "PA and lateral (standard)",
            "AP portable (ICU/ED)",
            "Lateral decubitus (effusion layering)",
            "Expiratory view (pneumothorax)",
        ],
        "typical_dose_msv": {
            "pa": 0.02,
            "lateral": 0.04,
            "ap_portable": 0.02,
        },
        "nim_models": [
            "DenseNet-121 (14-finding classification)",
            "VILA-M3 (report generation)",
        ],
    },
    "pet_ct": {
        "full_name": "Positron Emission Tomography / CT",
        "dicom_modality_code": "PT",
        "physics": (
            "Radiotracer (e.g., 18F-FDG) injected IV; positron emission "
            "produces annihilation photons detected by ring detectors. "
            "Co-registered with CT for anatomic localization. SUV "
            "(Standardized Uptake Value) quantifies metabolic activity."
        ),
        "strengths": [
            "Functional/metabolic imaging",
            "Whole-body staging in oncology",
            "Treatment response assessment (PERCIST)",
            "Detection of occult metastases",
            "SUV provides quantitative metric",
        ],
        "limitations": [
            "Combined radiation dose (PET + CT)",
            "Low spatial resolution (~4-5 mm)",
            "False positives (inflammation, infection)",
            "FDG not specific to malignancy",
            "Expensive, limited availability",
            "Patient must fast 4-6 hours",
        ],
        "common_protocols": [
            "FDG PET/CT whole body (skull base to mid-thigh)",
            "FDG PET/CT brain",
            "PSMA PET/CT (prostate cancer)",
            "Dotatate PET/CT (neuroendocrine tumors)",
            "FLT PET (proliferation imaging)",
        ],
        "typical_dose_msv": {
            "fdg_whole_body": 14.0,  # ~7 mSv PET + ~7 mSv low-dose CT
            "fdg_brain": 7.0,
            "psma_whole_body": 12.0,
        },
        "nim_models": [
            "VISTA-3D (CT component segmentation)",
        ],
    },
    "ultrasound": {
        "full_name": "Ultrasound (Sonography)",
        "dicom_modality_code": "US",
        "physics": (
            "Piezoelectric transducer emits high-frequency sound waves "
            "(2-18 MHz). Reflected echoes create real-time grayscale "
            "images. Doppler mode measures blood flow velocity. "
            "No ionizing radiation."
        ),
        "strengths": [
            "Real-time imaging",
            "No ionizing radiation (safe in pregnancy)",
            "Portable and bedside capable",
            "Low cost",
            "Doppler for vascular flow assessment",
            "Image-guided procedures (biopsy, drainage)",
        ],
        "limitations": [
            "Operator-dependent",
            "Limited by body habitus (obese patients)",
            "Air and bone attenuate sound (poor for lungs, brain)",
            "Limited field of view",
            "Reduced reproducibility",
        ],
        "common_protocols": [
            "Abdominal (liver, gallbladder, kidneys, pancreas)",
            "Pelvic (transabdominal, transvaginal)",
            "Thyroid",
            "Carotid Doppler",
            "Renal Doppler",
            "Obstetric (1st/2nd/3rd trimester)",
            "Echocardiography (TTE, TEE)",
            "Point-of-care (FAST exam)",
        ],
        "typical_dose_msv": None,  # Non-ionizing
        "nim_models": [],  # No current NIM models for US
    },
    "mammography": {
        "full_name": "Mammography",
        "dicom_modality_code": "MG",
        "physics": (
            "Low-energy X-ray (25-35 kVp) with breast compression. "
            "Digital mammography (DM) or digital breast tomosynthesis "
            "(DBT/3D). Anode/filter combinations (Mo/Mo, Rh/Rh, W/Rh) "
            "optimize contrast for breast tissue."
        ),
        "strengths": [
            "Proven screening reduces breast cancer mortality",
            "High sensitivity in fatty breasts",
            "DBT reduces recall rates vs 2D",
            "BI-RADS standardized reporting",
            "Cost-effective screening tool",
        ],
        "limitations": [
            "Reduced sensitivity in dense breasts",
            "Ionizing radiation (though low dose)",
            "Compression can be uncomfortable",
            "Limited for implant evaluation",
            "False positives lead to biopsies",
        ],
        "common_protocols": [
            "Screening mammogram (CC and MLO bilateral)",
            "Diagnostic mammogram (additional views, spot compression)",
            "Digital breast tomosynthesis (3D)",
            "Contrast-enhanced mammography (CEM)",
        ],
        "typical_dose_msv": {
            "bilateral_screening": 0.4,
            "diagnostic": 0.6,
            "tomosynthesis": 0.5,
        },
        "nim_models": [],  # Future mammography NIM planned
    },
    "fluoroscopy": {
        "full_name": "Fluoroscopy",
        "dicom_modality_code": "RF",
        "physics": (
            "Continuous or pulsed X-ray beam produces real-time moving "
            "images on a fluorescent screen or digital detector. "
            "Used for dynamic studies and image-guided procedures. "
            "Dose management via pulsed mode and collimation."
        ),
        "strengths": [
            "Real-time dynamic imaging",
            "Essential for GI studies (barium swallow, enema)",
            "Interventional procedure guidance",
            "Assessment of joint motion",
            "Catheter and device placement",
        ],
        "limitations": [
            "Significant cumulative radiation dose",
            "Limited soft tissue contrast",
            "Operator-dependent technique",
            "2D projection only",
            "Patient and staff radiation exposure",
        ],
        "common_protocols": [
            "Upper GI series (barium swallow)",
            "Barium enema",
            "Voiding cystourethrogram (VCUG)",
            "Hysterosalpingography (HSG)",
            "Arthrography",
            "Catheter angiography",
        ],
        "typical_dose_msv": {
            "upper_gi": 6.0,
            "barium_enema": 8.0,
            "vcug": 1.0,
            "angiography": 10.0,
        },
        "nim_models": [],  # No current NIM models for fluoroscopy
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3. IMAGING_ANATOMY — Anatomy knowledge graph (~15 entries)
# ═══════════════════════════════════════════════════════════════════════════════

IMAGING_ANATOMY: Dict[str, Dict[str, Any]] = {
    "brain": {
        "display_name": "Brain",
        "structures": [
            "cerebral cortex", "white matter", "basal ganglia",
            "thalamus", "brainstem", "cerebellum",
            "ventricles", "corpus callosum", "hippocampus",
        ],
        "vista3d_labels": [1, 2, 3, 4, 5, 6, 7, 8, 9],
        "common_pathologies": [
            "intracranial_hemorrhage", "stroke_ischemic",
            "brain_tumor", "ms_lesion",
        ],
        "preferred_modality": "mri",
        "key_sequences": [
            "T1 (anatomy)", "T2 (edema/fluid)", "FLAIR (periventricular)",
            "DWI/ADC (ischemia)", "SWI (hemorrhage/calcification)",
            "T1+Gd (enhancement/BBB breakdown)",
        ],
    },
    "lungs": {
        "display_name": "Lungs",
        "structures": [
            "right upper lobe", "right middle lobe", "right lower lobe",
            "left upper lobe", "lingula", "left lower lobe",
            "main bronchi", "segmental bronchi", "pulmonary vasculature",
        ],
        "vista3d_labels": [10, 11, 12, 13, 14, 15],
        "common_pathologies": [
            "lung_nodule", "pneumonia", "pneumothorax",
            "pulmonary_embolism", "pulmonary_fibrosis",
        ],
        "preferred_modality": "ct",
        "protocols": [
            "Standard chest CT (soft tissue + lung windows)",
            "HRCT (1 mm slices, prone if needed)",
            "Low-dose lung cancer screening",
            "CT pulmonary angiography (CTPA)",
        ],
    },
    "heart": {
        "display_name": "Heart",
        "structures": [
            "left ventricle", "right ventricle",
            "left atrium", "right atrium",
            "interventricular septum", "aortic valve",
            "mitral valve", "tricuspid valve", "pulmonic valve",
            "pericardium", "coronary arteries",
        ],
        "vista3d_labels": [20, 21, 22, 23, 24, 25],
        "common_pathologies": ["cardiomegaly"],
        "preferred_modality": "mri",
        "key_sequences": [
            "Cine SSFP (function, EF)",
            "Late gadolinium enhancement (fibrosis/scar)",
            "T1/T2 mapping (edema, infiltration)",
            "Phase contrast (flow quantification)",
        ],
    },
    "liver": {
        "display_name": "Liver",
        "structures": [
            "right lobe (segments V-VIII)",
            "left lobe (segments II-IV)",
            "caudate lobe (segment I)",
            "hepatic veins", "portal vein",
            "hepatic artery", "gallbladder",
            "common bile duct",
        ],
        "vista3d_labels": [30, 31, 32, 33, 34],
        "common_pathologies": ["liver_lesion"],
        "preferred_modality": "mri",
        "key_sequences": [
            "T1 in/opposed-phase (fat/iron)",
            "T2 (lesion characterization)",
            "DWI (malignancy detection)",
            "Multiphase contrast (arterial, portal venous, delayed)",
            "Hepatobiliary phase (Eovist 20 min)",
        ],
    },
    "kidneys": {
        "display_name": "Kidneys",
        "structures": [
            "renal cortex", "renal medulla",
            "renal pelvis", "calyces",
            "renal artery", "renal vein",
            "ureter (proximal)",
        ],
        "vista3d_labels": [40, 41, 42, 43],
        "common_pathologies": ["renal_mass"],
        "preferred_modality": "ct",
        "protocols": [
            "CT urogram (non-contrast, nephrographic, delayed)",
            "Renal mass protocol (multiphase)",
            "MRI with subtraction (for small lesions)",
        ],
    },
    "spine": {
        "display_name": "Spine",
        "structures": [
            "cervical vertebrae (C1-C7)",
            "thoracic vertebrae (T1-T12)",
            "lumbar vertebrae (L1-L5)",
            "sacrum", "coccyx",
            "intervertebral discs",
            "spinal cord", "nerve roots",
            "spinal canal", "facet joints",
        ],
        "vista3d_labels": [50, 51, 52, 53, 54, 55],
        "common_pathologies": ["fracture"],
        "preferred_modality": "mri",
        "key_sequences": [
            "T1 sagittal (marrow signal, anatomy)",
            "T2 sagittal (cord, disc, fluid)",
            "STIR (edema, acute fracture)",
            "T2 axial (disc herniation, stenosis)",
            "T1+Gd (infection, tumor enhancement)",
        ],
    },
    "aorta": {
        "display_name": "Aorta",
        "structures": [
            "aortic root", "ascending aorta",
            "aortic arch", "descending thoracic aorta",
            "abdominal aorta", "aortic bifurcation",
            "iliac arteries",
        ],
        "vista3d_labels": [60, 61, 62, 63],
        "common_pathologies": ["aortic_dissection"],
        "preferred_modality": "ct",
        "protocols": [
            "CT angiography aorta (ECG-gated for root)",
            "CTA runoff (aorto-iliac to pedal)",
            "MRA aorta (non-contrast or gadolinium-enhanced)",
        ],
    },
    "mediastinum": {
        "display_name": "Mediastinum",
        "structures": [
            "anterior mediastinum (thymus, lymph nodes)",
            "middle mediastinum (heart, great vessels, trachea)",
            "posterior mediastinum (esophagus, descending aorta, spine)",
            "hilar lymph nodes", "paratracheal lymph nodes",
            "subcarinal lymph nodes",
        ],
        "vista3d_labels": [70, 71, 72],
        "common_pathologies": [
            "lung_nodule",  # mediastinal lymphadenopathy
        ],
        "preferred_modality": "ct",
        "protocols": [
            "Contrast-enhanced chest CT",
            "PET/CT for lymph node staging",
            "MRI for posterior mediastinal masses",
        ],
    },
    "pelvis": {
        "display_name": "Pelvis",
        "structures": [
            "bladder", "rectum", "prostate (male)",
            "uterus (female)", "ovaries (female)",
            "pelvic sidewall musculature",
            "iliac vessels", "pelvic lymph nodes",
            "sacroiliac joints",
        ],
        "vista3d_labels": [80, 81, 82, 83, 84],
        "common_pathologies": [],  # Phase 2: cervical/prostate/rectal cancers
        "preferred_modality": "mri",
        "key_sequences": [
            "T2 high-resolution (anatomy, staging)",
            "DWI (malignancy detection)",
            "T1+Gd (enhancement, lymph nodes)",
            "Dynamic contrast (cervical cancer)",
        ],
    },
    "breast": {
        "display_name": "Breast",
        "structures": [
            "fibroglandular tissue", "adipose tissue",
            "ducts", "lobules",
            "Cooper ligaments", "pectoralis muscle",
            "axillary lymph nodes",
        ],
        "vista3d_labels": [90, 91],
        "common_pathologies": [],  # Phase 2: breast_mass
        "preferred_modality": "mammography",
        "protocols": [
            "Screening mammogram (CC + MLO)",
            "Diagnostic mammogram (additional views)",
            "Breast MRI (dynamic contrast-enhanced)",
            "Breast ultrasound (targeted, whole-breast)",
        ],
    },
    "thyroid": {
        "display_name": "Thyroid",
        "structures": [
            "right lobe", "left lobe", "isthmus",
            "pyramidal lobe (variant)",
            "cervical lymph nodes (levels II-VI)",
        ],
        "vista3d_labels": [100, 101],
        "common_pathologies": [],  # Phase 2: thyroid_nodule
        "preferred_modality": "ultrasound",
        "protocols": [
            "Thyroid ultrasound with Doppler",
            "US-guided FNA biopsy",
            "CT neck with contrast (for large goiter/invasion)",
            "Nuclear medicine thyroid scan (I-123/Tc-99m)",
        ],
    },
    "pancreas": {
        "display_name": "Pancreas",
        "structures": [
            "head", "uncinate process", "neck",
            "body", "tail",
            "pancreatic duct (Wirsung)",
            "common bile duct (intrapancreatic portion)",
        ],
        "vista3d_labels": [110, 111],
        "common_pathologies": [],  # Phase 2: pancreatic_mass
        "preferred_modality": "ct",
        "protocols": [
            "Pancreas protocol CT (arterial + portal venous phase)",
            "MRCP (ductal evaluation)",
            "Endoscopic ultrasound (EUS) for small lesions",
        ],
    },
    "bowel": {
        "display_name": "Bowel",
        "structures": [
            "duodenum", "jejunum", "ileum",
            "cecum", "ascending colon", "transverse colon",
            "descending colon", "sigmoid colon", "rectum",
            "appendix",
        ],
        "vista3d_labels": [120, 121, 122, 123, 124],
        "common_pathologies": [],  # Phase 2: bowel_obstruction, appendicitis
        "preferred_modality": "ct",
        "protocols": [
            "CT abdomen/pelvis with IV + oral contrast",
            "CT enterography (small bowel Crohn's)",
            "CT colonography (virtual colonoscopy)",
            "MR enterography (radiation-free alternative)",
        ],
    },
    "ribs": {
        "display_name": "Ribs",
        "structures": [
            "true ribs (1-7)", "false ribs (8-10)",
            "floating ribs (11-12)",
            "costal cartilage", "costochondral junction",
            "intercostal spaces",
        ],
        "vista3d_labels": [130, 131, 132],
        "common_pathologies": ["fracture"],
        "preferred_modality": "ct",
        "protocols": [
            "Dedicated rib CT (thin-slice, bone windows)",
            "Chest X-ray (limited sensitivity for non-displaced)",
            "Bone scan (occult fractures, metastases)",
        ],
    },
    "skull": {
        "display_name": "Skull",
        "structures": [
            "frontal bone", "parietal bones", "temporal bones",
            "occipital bone", "sphenoid bone", "ethmoid bone",
            "maxilla", "mandible",
            "orbits", "paranasal sinuses",
            "skull base foramina",
        ],
        "vista3d_labels": [140, 141, 142, 143],
        "common_pathologies": ["fracture", "intracranial_hemorrhage"],
        "preferred_modality": "ct",
        "protocols": [
            "Non-contrast head CT (trauma, hemorrhage)",
            "CT facial bones (thin-slice, multiplanar reformat)",
            "CT temporal bones (0.6 mm, bone algorithm)",
            "CT orbits with contrast",
        ],
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 4. HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def get_pathology_context(key: str) -> str:
    """Return formatted knowledge context for an imaging pathology.

    Args:
        key: Pathology key (e.g. 'intracranial_hemorrhage', 'lung_nodule').
             Case-insensitive, spaces converted to underscores.

    Returns:
        Formatted markdown string with pathology knowledge, or empty string
        if not found.
    """
    normalized = key.strip().lower().replace(" ", "_").replace("-", "_")
    data = IMAGING_PATHOLOGIES.get(normalized)
    if not data:
        # Partial/substring match
        for k in IMAGING_PATHOLOGIES:
            if normalized in k or k in normalized:
                data = IMAGING_PATHOLOGIES[k]
                normalized = k
                break
    if not data:
        return ""

    lines = [
        f"## Pathology: {data['display_name']}",
        f"- **ICD-10:** {data['icd10']}",
        f"- **Body Region:** {data['body_region']}",
        f"- **Modalities:** {', '.join(data['modalities'])}",
        f"- **Subtypes:** {', '.join(data['subtypes'])}",
        f"- **CT Characteristics:** {data['ct_characteristics']}",
        f"- **MRI Characteristics:** {data['mri_characteristics']}",
    ]

    if data.get("severity_criteria"):
        lines.append("- **Severity Criteria:**")
        for level, criteria in data["severity_criteria"].items():
            lines.append(f"  - {level.title()}: {criteria}")

    if data.get("key_measurements"):
        lines.append(
            f"- **Key Measurements:** {', '.join(data['key_measurements'])}"
        )
    if data.get("clinical_guidelines"):
        lines.append(
            f"- **Clinical Guidelines:** {', '.join(data['clinical_guidelines'])}"
        )
    if data.get("ai_models"):
        lines.append(f"- **AI Models:** {', '.join(data['ai_models'])}")
    if data.get("nim_workflow"):
        lines.append(f"- **NIM Workflow:** {data['nim_workflow']}")

    return "\n".join(lines)


def get_modality_context(key: str) -> str:
    """Return formatted knowledge context for an imaging modality.

    Args:
        key: Modality key (e.g. 'ct', 'mri', 'cxr').
             Case-insensitive.

    Returns:
        Formatted markdown string with modality knowledge, or empty string
        if not found.
    """
    normalized = key.strip().lower().replace(" ", "_").replace("-", "_")
    data = IMAGING_MODALITIES.get(normalized)
    if not data:
        # Try matching by full name
        for k, v in IMAGING_MODALITIES.items():
            if normalized in v["full_name"].lower().replace(" ", "_"):
                data = v
                normalized = k
                break
    if not data:
        return ""

    lines = [
        f"## Modality: {data['full_name']}",
        f"- **DICOM Code:** {data['dicom_modality_code']}",
        f"- **Physics:** {data['physics']}",
    ]

    if data.get("strengths"):
        lines.append("- **Strengths:**")
        for s in data["strengths"]:
            lines.append(f"  - {s}")

    if data.get("limitations"):
        lines.append("- **Limitations:**")
        for lim in data["limitations"]:
            lines.append(f"  - {lim}")

    if data.get("common_protocols"):
        lines.append("- **Common Protocols:**")
        for p in data["common_protocols"]:
            lines.append(f"  - {p}")

    if data.get("typical_dose_msv") is not None:
        lines.append("- **Typical Dose (mSv):**")
        for exam, dose in data["typical_dose_msv"].items():
            label = exam.replace("_", " ").title()
            lines.append(f"  - {label}: {dose} mSv")
    elif data.get("typical_dose_msv") is None:
        lines.append("- **Radiation:** Non-ionizing (no radiation dose)")

    if data.get("nim_models"):
        lines.append(f"- **NIM Models:** {', '.join(data['nim_models'])}")

    return "\n".join(lines)


def get_anatomy_context(key: str) -> str:
    """Return formatted knowledge context for an anatomical region.

    Args:
        key: Anatomy key (e.g. 'brain', 'lungs', 'liver').
             Case-insensitive, spaces converted to underscores.

    Returns:
        Formatted markdown string with anatomy knowledge, or empty string
        if not found.
    """
    normalized = key.strip().lower().replace(" ", "_").replace("-", "_")
    data = IMAGING_ANATOMY.get(normalized)
    if not data:
        for k in IMAGING_ANATOMY:
            if normalized in k or k in normalized:
                data = IMAGING_ANATOMY[k]
                normalized = k
                break
    if not data:
        return ""

    lines = [
        f"## Anatomy: {data['display_name']}",
        f"- **Structures:** {', '.join(data['structures'])}",
        f"- **VISTA-3D Labels:** {data['vista3d_labels']}",
        f"- **Preferred Modality:** {data['preferred_modality']}",
    ]

    if data.get("common_pathologies"):
        lines.append(
            f"- **Common Pathologies:** "
            f"{', '.join(data['common_pathologies'])}"
        )

    # MRI-primary regions have key_sequences; CT-primary have protocols
    if data.get("key_sequences"):
        lines.append("- **Key MRI Sequences:**")
        for seq in data["key_sequences"]:
            lines.append(f"  - {seq}")
    if data.get("protocols"):
        lines.append("- **Protocols:**")
        for proto in data["protocols"]:
            lines.append(f"  - {proto}")

    return "\n".join(lines)


def get_nim_recommendation(pathology_key: str) -> Optional[str]:
    """Return the recommended NIM workflow name for a pathology.

    Args:
        pathology_key: Pathology key (e.g. 'lung_nodule').

    Returns:
        NIM workflow string (e.g. 'ct_chest_lung_nodule'), or None if
        no NIM workflow is configured for that pathology.
    """
    normalized = pathology_key.strip().lower().replace(" ", "_").replace("-", "_")
    data = IMAGING_PATHOLOGIES.get(normalized)
    if not data:
        for k in IMAGING_PATHOLOGIES:
            if normalized in k or k in normalized:
                data = IMAGING_PATHOLOGIES[k]
                break
    if not data:
        return None
    return data.get("nim_workflow")


def resolve_comparison_entity(raw: str) -> Optional[Dict[str, str]]:
    """Resolve a raw text string to a known imaging entity for comparison.

    Searches all three knowledge dictionaries (pathologies, modalities,
    anatomy) in priority order.

    Args:
        raw: Free-text entity name (e.g. 'CT', 'brain tumor', 'lungs').

    Returns:
        Dict with 'type' and 'canonical' keys, or None if not recognized.
    """
    cleaned = raw.strip().lower().replace(" ", "_").replace("-", "_")

    # 1. Pathologies (exact match)
    if cleaned in IMAGING_PATHOLOGIES:
        return {"type": "pathology", "canonical": cleaned}

    # 2. Modalities (exact match)
    if cleaned in IMAGING_MODALITIES:
        return {"type": "modality", "canonical": cleaned}

    # 3. Anatomy (exact match)
    if cleaned in IMAGING_ANATOMY:
        return {"type": "anatomy", "canonical": cleaned}

    # 4. Pathologies (partial/substring match)
    for k in IMAGING_PATHOLOGIES:
        if cleaned in k or k in cleaned:
            return {"type": "pathology", "canonical": k}

    # 5. Modalities (partial match by full_name)
    for k, v in IMAGING_MODALITIES.items():
        full = v["full_name"].lower().replace(" ", "_")
        if cleaned in full or full in cleaned or cleaned in k:
            return {"type": "modality", "canonical": k}

    # 6. Anatomy (partial/substring match)
    for k in IMAGING_ANATOMY:
        if cleaned in k or k in cleaned:
            return {"type": "anatomy", "canonical": k}

    return None


def get_comparison_context(
    entity_a: Dict[str, str],
    entity_b: Dict[str, str],
) -> str:
    """Build side-by-side knowledge graph context for two imaging entities.

    Reuses get_pathology_context / get_modality_context / get_anatomy_context
    depending on entity type.

    Args:
        entity_a: Dict with 'type' and 'canonical' keys (from
                  resolve_comparison_entity).
        entity_b: Dict with 'type' and 'canonical' keys.

    Returns:
        Formatted comparison context string with both entities' data,
        separated by a horizontal rule.
    """

    def _get_entity_context(entity: Dict[str, str]) -> str:
        etype = entity["type"]
        canonical = entity["canonical"]
        if etype == "pathology":
            return get_pathology_context(canonical)
        elif etype == "modality":
            return get_modality_context(canonical)
        elif etype == "anatomy":
            return get_anatomy_context(canonical)
        return ""

    sections: List[str] = []
    ctx_a = _get_entity_context(entity_a)
    ctx_b = _get_entity_context(entity_b)

    if ctx_a:
        sections.append(f"### {entity_a['canonical']}\n{ctx_a}")
    if ctx_b:
        sections.append(f"### {entity_b['canonical']}\n{ctx_b}")

    return "\n\n---\n\n".join(sections)


def get_knowledge_stats() -> Dict[str, int]:
    """Return statistics about the imaging knowledge graph.

    Returns:
        Dict with counts for each knowledge dictionary plus derived
        aggregate metrics.
    """
    pathologies_with_nim = sum(
        1 for p in IMAGING_PATHOLOGIES.values()
        if p.get("nim_workflow") is not None
    )
    modalities_ionizing = sum(
        1 for m in IMAGING_MODALITIES.values()
        if m.get("typical_dose_msv") is not None
    )
    anatomy_with_pathologies = sum(
        1 for a in IMAGING_ANATOMY.values()
        if a.get("common_pathologies")
    )

    return {
        "pathologies": len(IMAGING_PATHOLOGIES),
        "pathologies_with_nim_workflow": pathologies_with_nim,
        "modalities": len(IMAGING_MODALITIES),
        "modalities_ionizing": modalities_ionizing,
        "anatomy_regions": len(IMAGING_ANATOMY),
        "anatomy_with_pathologies": anatomy_with_pathologies,
    }
