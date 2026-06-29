"""Neuroimaging protocol and findings parser for the Neurology Intelligence Agent.

Parses neuroimaging protocols and seeds 70 key neuroimaging protocols and
findings covering MRI sequences, CT protocols, PET tracers, SPECT studies,
and characteristic imaging findings for major neurological conditions.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA: 70 NEUROIMAGING PROTOCOLS AND FINDINGS
# ===================================================================

NEUROIMAGING_PROTOCOLS: List[Dict[str, Any]] = [
    # --- MRI Stroke Protocols ---
    {
        "protocol_id": "NI-001",
        "name": "Acute Stroke MRI Protocol",
        "modality": "MRI",
        "domain": "stroke",
        "sequences": ["DWI", "ADC", "FLAIR", "GRE/SWI", "MRA", "PWI"],
        "indication": "Acute ischemic stroke evaluation within therapeutic window",
        "key_findings": "DWI restriction in ischemic core; PWI-DWI mismatch defines penumbra; GRE/SWI for hemorrhagic transformation",
        "clinical_relevance": "Guides thrombolysis and thrombectomy decisions; ASPECTS scoring on DWI",
    },
    {
        "protocol_id": "NI-002",
        "name": "CT Perfusion for Stroke",
        "modality": "CT",
        "domain": "stroke",
        "sequences": ["NCCT", "CTA", "CTP"],
        "indication": "Extended window stroke assessment for thrombectomy candidacy",
        "key_findings": "CBF<30% of contralateral = ischemic core; Tmax>6s = hypoperfused tissue; mismatch ratio for DAWN/DEFUSE eligibility",
        "clinical_relevance": "Identifies patients eligible for thrombectomy 6-24h post-onset per DAWN/DEFUSE criteria",
    },
    {
        "protocol_id": "NI-003",
        "name": "Carotid MRA Protocol",
        "modality": "MRI",
        "domain": "cerebrovascular",
        "sequences": ["3D TOF MRA", "CE-MRA", "vessel wall imaging"],
        "indication": "Extracranial and intracranial arterial stenosis evaluation",
        "key_findings": "NASCET-criteria stenosis measurement; plaque characterization on vessel wall imaging",
        "clinical_relevance": "Guides carotid endarterectomy vs stenting decisions; identifies vulnerable plaque",
    },
    # --- Dementia Imaging ---
    {
        "protocol_id": "NI-004",
        "name": "Alzheimer's Disease MRI Protocol",
        "modality": "MRI",
        "domain": "dementia",
        "sequences": ["3D T1 MPRAGE", "FLAIR", "T2*", "ASL perfusion"],
        "indication": "Evaluation of suspected Alzheimer's disease",
        "key_findings": "Medial temporal lobe atrophy (Scheltens scale); posterior cortical atrophy; global cortical atrophy (GCA scale)",
        "clinical_relevance": "Supports ATN framework staging; correlates with Braak staging",
    },
    {
        "protocol_id": "NI-005",
        "name": "Amyloid PET Imaging",
        "modality": "PET",
        "domain": "dementia",
        "sequences": ["18F-florbetapir", "18F-florbetaben", "18F-flutemetamol"],
        "indication": "In vivo amyloid-beta plaque detection for Alzheimer's diagnosis",
        "key_findings": "Positive: diffuse cortical retention; Negative: white matter only; supports ATN A+ classification",
        "clinical_relevance": "Appropriate use criteria per AUC 2.0; required for anti-amyloid therapy eligibility",
    },
    {
        "protocol_id": "NI-006",
        "name": "Tau PET Imaging",
        "modality": "PET",
        "domain": "dementia",
        "sequences": ["18F-flortaucipir (AV-1451)", "18F-MK6240", "18F-PI-2620"],
        "indication": "In vivo tau neurofibrillary tangle detection",
        "key_findings": "Temporal lobe uptake in AD; distribution correlates with Braak stages; frontal predominance in FTD variants",
        "clinical_relevance": "Supports ATN T+ classification; correlates with cognitive decline better than amyloid",
    },
    {
        "protocol_id": "NI-007",
        "name": "FDG-PET for Dementia Differential",
        "modality": "PET",
        "domain": "dementia",
        "sequences": ["18F-FDG"],
        "indication": "Metabolic pattern analysis for dementia differential diagnosis",
        "key_findings": "AD: temporoparietal hypometabolism; FTD: frontal/anterior temporal; DLB: occipital hypometabolism; cingulate island sign in DLB",
        "clinical_relevance": "Distinguishes AD, FTD, DLB with >85% sensitivity; normal in depression-related cognitive complaints",
    },
    {
        "protocol_id": "NI-008",
        "name": "Frontotemporal Dementia MRI Pattern",
        "modality": "MRI",
        "domain": "dementia",
        "sequences": ["3D T1 MPRAGE", "FLAIR"],
        "indication": "Evaluation of suspected frontotemporal dementia",
        "key_findings": "Frontal and anterior temporal atrophy (knife-edge pattern); behavioral variant: orbitofrontal/dorsolateral; semantic: anterior temporal pole",
        "clinical_relevance": "Distinguishes behavioral variant from semantic and nonfluent/agrammatic subtypes",
    },
    # --- Epilepsy Imaging ---
    {
        "protocol_id": "NI-009",
        "name": "Epilepsy MRI Protocol (3T)",
        "modality": "MRI",
        "domain": "epilepsy",
        "sequences": ["3D T1 MPRAGE", "coronal T2 FLAIR", "coronal T2", "SWI", "3D T1 post-contrast"],
        "indication": "Structural evaluation for seizure focus localization",
        "key_findings": "Hippocampal sclerosis (T2 hyperintensity + volume loss); cortical dysplasia; cavernomas; tumors",
        "clinical_relevance": "Identifies surgical targets; MRI-negative epilepsy requires advanced post-processing",
    },
    {
        "protocol_id": "NI-010",
        "name": "SISCOM (Subtraction Ictal SPECT)",
        "modality": "SPECT",
        "domain": "epilepsy",
        "sequences": ["ictal 99mTc-HMPAO SPECT", "interictal SPECT", "co-registered MRI"],
        "indication": "Seizure focus localization in MRI-negative epilepsy",
        "key_findings": "Hyperperfusion at seizure onset zone on subtraction images co-registered with MRI",
        "clinical_relevance": "Adjunctive localization tool for epilepsy surgery planning when MRI is non-lesional",
    },
    {
        "protocol_id": "NI-011",
        "name": "Epilepsy PET (FDG)",
        "modality": "PET",
        "domain": "epilepsy",
        "sequences": ["18F-FDG interictal"],
        "indication": "Interictal metabolic assessment for seizure focus",
        "key_findings": "Focal hypometabolism in epileptogenic zone; temporal lobe sensitivity >80%",
        "clinical_relevance": "Concordance with EEG and MRI increases surgical success rate",
    },
    {
        "protocol_id": "NI-012",
        "name": "MEG for Epilepsy Mapping",
        "modality": "MEG",
        "domain": "epilepsy",
        "sequences": ["whole-head MEG", "equivalent current dipole modeling"],
        "indication": "Non-invasive source localization of epileptiform discharges",
        "key_findings": "Dipole clustering at seizure onset zone; concordance with intracranial EEG",
        "clinical_relevance": "Guides intracranial electrode placement; especially useful for neocortical epilepsy",
    },
    # --- MS Imaging ---
    {
        "protocol_id": "NI-013",
        "name": "MS Brain MRI Protocol",
        "modality": "MRI",
        "domain": "ms",
        "sequences": ["3D FLAIR", "3D T1 pre/post-Gd", "DIR", "SWI"],
        "indication": "Multiple sclerosis diagnosis and monitoring",
        "key_findings": "Periventricular (Dawson fingers), juxtacortical, infratentorial, cortical lesions; Gd-enhancing = active",
        "clinical_relevance": "McDonald 2017 criteria: DIS requires 2+ of 4 locations; DIT by simultaneous Gd+ and non-Gd+ lesions",
    },
    {
        "protocol_id": "NI-014",
        "name": "MS Spinal Cord MRI Protocol",
        "modality": "MRI",
        "domain": "ms",
        "sequences": ["sagittal T2", "sagittal STIR", "axial T2", "post-Gd T1"],
        "indication": "Spinal cord lesion detection in MS",
        "key_findings": "Short-segment (<2 vertebral segments) peripheral lesions; dorsal column involvement; mild cord swelling if acute",
        "clinical_relevance": "Spinal cord lesions count toward DIS; distinguish from NMOSD (longitudinally extensive >3 segments)",
    },
    {
        "protocol_id": "NI-015",
        "name": "Central Vein Sign Assessment",
        "modality": "MRI",
        "domain": "ms",
        "sequences": ["T2*-weighted", "FLAIR*", "3T or 7T susceptibility"],
        "indication": "Differentiating MS from mimics using central vein sign",
        "key_findings": "Central vein visible in >40% of white matter lesions suggests MS; absent in migraine, SVD, NMOSD lesions",
        "clinical_relevance": "Proposed biomarker with >90% specificity for MS; may replace DIT requirement in future criteria",
    },
    # --- Parkinson's / Movement Disorders ---
    {
        "protocol_id": "NI-016",
        "name": "DaTscan (123I-Ioflupane SPECT)",
        "modality": "SPECT",
        "domain": "parkinson",
        "sequences": ["123I-ioflupane SPECT"],
        "indication": "Presynaptic dopaminergic deficit evaluation",
        "key_findings": "Reduced putaminal uptake (comma to period sign); asymmetric loss; normal in essential tremor",
        "clinical_relevance": "Distinguishes degenerative parkinsonism from essential tremor, drug-induced, psychogenic tremor",
    },
    {
        "protocol_id": "NI-017",
        "name": "SWI for Nigrosome-1 (Swallow Tail Sign)",
        "modality": "MRI",
        "domain": "parkinson",
        "sequences": ["SWI", "QSM (quantitative susceptibility mapping)"],
        "indication": "Substantia nigra evaluation in parkinsonism",
        "key_findings": "Loss of nigrosome-1 hyperintensity (swallow tail sign) on SWI; increased iron deposition on QSM",
        "clinical_relevance": "High sensitivity (>90%) for PD at 3T; may reduce need for DaTscan in some patients",
    },
    {
        "protocol_id": "NI-018",
        "name": "MIBG Cardiac Scintigraphy",
        "modality": "SPECT",
        "domain": "parkinson",
        "sequences": ["123I-MIBG cardiac"],
        "indication": "Cardiac sympathetic denervation in Lewy body disorders",
        "key_findings": "Low heart-to-mediastinum ratio in PD, DLB; preserved in MSA, PSP",
        "clinical_relevance": "Distinguishes PD/DLB from MSA/PSP with high specificity; reflects peripheral autonomic degeneration",
    },
    # --- Brain Tumors ---
    {
        "protocol_id": "NI-019",
        "name": "Brain Tumor MRI Protocol",
        "modality": "MRI",
        "domain": "neuro_oncology",
        "sequences": ["T1 pre/post-Gd", "T2", "FLAIR", "DWI", "SWI", "perfusion (DSC)", "MR spectroscopy"],
        "indication": "Brain tumor characterization and grading",
        "key_findings": "Enhancement pattern, perfusion rCBV, NAA/Cho/Cr ratios on MRS; DWI restriction in high-grade tumors",
        "clinical_relevance": "RANO criteria for response assessment; rCBV>1.75 suggests high-grade; lipid/lactate peaks in necrosis",
    },
    {
        "protocol_id": "NI-020",
        "name": "MR Spectroscopy for Glioma Grading",
        "modality": "MRI",
        "domain": "neuro_oncology",
        "sequences": ["single-voxel SVS", "multivoxel CSI at TE 135ms and 30ms"],
        "indication": "Non-invasive metabolic characterization of brain tumors",
        "key_findings": "Elevated Cho/NAA and Cho/Cr in tumor; 2-HG peak in IDH-mutant gliomas; lipid-lactate in high-grade",
        "clinical_relevance": "2-HG detection has >90% sensitivity for IDH mutation; helps distinguish recurrence from radiation necrosis",
    },
    {
        "protocol_id": "NI-021",
        "name": "Amino Acid PET for Brain Tumors",
        "modality": "PET",
        "domain": "neuro_oncology",
        "sequences": ["18F-FET PET", "11C-MET PET"],
        "indication": "Brain tumor delineation and recurrence detection",
        "key_findings": "Increased uptake in viable tumor; TBR>1.6 suggests high-grade; static and dynamic analysis for grading",
        "clinical_relevance": "Superior to MRI for tumor delineation; distinguishes true progression from pseudoprogression",
    },
    # --- Headache ---
    {
        "protocol_id": "NI-022",
        "name": "Headache Red Flag MRI Protocol",
        "modality": "MRI",
        "domain": "headache",
        "sequences": ["T1 pre/post-Gd", "FLAIR", "DWI", "MRA", "MRV"],
        "indication": "Evaluation of secondary headache (thunderclap, new-onset, atypical features)",
        "key_findings": "SAH on FLAIR; venous sinus thrombosis on MRV; pituitary apoplexy; RCVS vasoconstriction on MRA",
        "clinical_relevance": "Must exclude SAH, CVT, mass lesion, RCVS, IIH before diagnosing primary headache",
    },
    {
        "protocol_id": "NI-023",
        "name": "MRI for Spontaneous Intracranial Hypotension",
        "modality": "MRI",
        "domain": "headache",
        "sequences": ["brain MRI with Gd", "spine MRI with Gd", "CT myelography"],
        "indication": "Orthostatic headache with suspected CSF leak",
        "key_findings": "Diffuse pachymeningeal enhancement; brain sagging; subdural collections; spinal epidural fluid collections",
        "clinical_relevance": "Guides targeted epidural blood patch; dynamic CT myelography for fast leak localization",
    },
    # --- Neuromuscular ---
    {
        "protocol_id": "NI-024",
        "name": "Muscle MRI Protocol",
        "modality": "MRI",
        "domain": "neuromuscular",
        "sequences": ["axial T1", "axial STIR", "coronal STIR"],
        "indication": "Evaluation of myopathy pattern and muscle inflammation",
        "key_findings": "T1: fatty infiltration pattern; STIR: active inflammation/edema; distribution guides differential (proximal vs distal, symmetric vs asymmetric)",
        "clinical_relevance": "Guides muscle biopsy site selection; distinguishes inflammatory from dystrophic myopathies",
    },
    {
        "protocol_id": "NI-025",
        "name": "Brachial Plexus MRI Protocol",
        "modality": "MRI",
        "domain": "neuromuscular",
        "sequences": ["coronal STIR", "axial T1", "3D SPACE", "post-Gd T1"],
        "indication": "Evaluation of brachial plexopathy",
        "key_findings": "Nerve thickening and T2 hyperintensity; mass lesions; root avulsion; denervation edema in muscles",
        "clinical_relevance": "Distinguishes Parsonage-Turner from infiltrative, traumatic, and radiation-induced plexopathy",
    },
    # --- Additional MRI Protocols ---
    {
        "protocol_id": "NI-026",
        "name": "Neuromyelitis Optica Spectrum MRI",
        "modality": "MRI",
        "domain": "ms",
        "sequences": ["brain MRI", "full spine MRI", "orbital MRI"],
        "indication": "Evaluation of suspected NMOSD",
        "key_findings": "LETM (>3 segments) with central cord involvement; area postrema lesion; bilateral optic neuritis; periependymal lesions",
        "clinical_relevance": "Distinguishes NMOSD from MS; AQP4-Ab+ defines seronegative/seropositive; guides treatment selection",
    },
    {
        "protocol_id": "NI-027",
        "name": "Normal Pressure Hydrocephalus CT/MRI",
        "modality": "MRI",
        "domain": "dementia",
        "sequences": ["T1 MPRAGE", "FLAIR", "phase-contrast CSF flow study"],
        "indication": "Evaluation of suspected NPH (triad: gait, dementia, incontinence)",
        "key_findings": "Ventriculomegaly disproportionate to sulcal atrophy (DESH); Evans index >0.3; callosal angle <90 degrees; tight high-convexity sulci",
        "clinical_relevance": "Predicts response to shunt surgery; CSF tap test and extended lumbar drainage for confirmation",
    },
    {
        "protocol_id": "NI-028",
        "name": "Susceptibility-Weighted Imaging for Microbleeds",
        "modality": "MRI",
        "domain": "cerebrovascular",
        "sequences": ["SWI", "T2*-GRE"],
        "indication": "Detection of cerebral microbleeds for hemorrhagic risk stratification",
        "key_findings": "Lobar microbleeds: cerebral amyloid angiopathy; deep/infratentorial: hypertensive microangiopathy; superficial siderosis",
        "clinical_relevance": "Informs anticoagulation risk; Boston criteria for CAA; microbleed burden predicts ICH risk",
    },
    {
        "protocol_id": "NI-029",
        "name": "DTI Tractography",
        "modality": "MRI",
        "domain": "neuro_oncology",
        "sequences": ["diffusion tensor imaging", "multi-shell DWI", "probabilistic tractography"],
        "indication": "White matter tract mapping for surgical planning",
        "key_findings": "Corticospinal tract relationship to tumor; arcuate fasciculus for language mapping; optic radiations",
        "clinical_relevance": "Reduces postoperative deficits; informs extent of resection decisions; combines with fMRI",
    },
    {
        "protocol_id": "NI-030",
        "name": "Functional MRI Motor and Language Mapping",
        "modality": "MRI",
        "domain": "neuro_oncology",
        "sequences": ["BOLD fMRI", "task-based paradigms", "resting-state fMRI"],
        "indication": "Pre-surgical eloquent cortex localization",
        "key_findings": "Motor cortex activation; Broca/Wernicke area lateralization; hemispheric dominance index",
        "clinical_relevance": "Wada test replacement for language lateralization; guides craniotomy planning; identifies reorganized networks",
    },
    # --- Additional Specialized Protocols ---
    {
        "protocol_id": "NI-031",
        "name": "CJD MRI Pattern",
        "modality": "MRI",
        "domain": "dementia",
        "sequences": ["DWI", "FLAIR", "T2"],
        "indication": "Evaluation of rapidly progressive dementia",
        "key_findings": "Cortical ribboning on DWI; caudate and putamen hyperintensity; pulvinar sign (variant CJD)",
        "clinical_relevance": "DWI has >90% sensitivity for sporadic CJD; RT-QuIC assay in CSF is confirmatory",
    },
    {
        "protocol_id": "NI-032",
        "name": "MRI Patterns in Autoimmune Encephalitis",
        "modality": "MRI",
        "domain": "epilepsy",
        "sequences": ["FLAIR", "T2", "post-Gd T1"],
        "indication": "Suspected autoimmune encephalitis (anti-NMDAR, anti-LGI1, anti-CASPR2)",
        "key_findings": "Medial temporal FLAIR hyperintensity; bilateral in LGI1; normal or subtle in anti-NMDAR; basal ganglia involvement in anti-CASPR2",
        "clinical_relevance": "Early immunotherapy improves outcomes; MRI may be normal in 50% of anti-NMDAR encephalitis",
    },
    {
        "protocol_id": "NI-033",
        "name": "Cervical Spine MRI for Myelopathy",
        "modality": "MRI",
        "domain": "neuromuscular",
        "sequences": ["sagittal T2", "sagittal T1", "axial T2", "post-Gd if needed"],
        "indication": "Cervical spondylotic myelopathy evaluation",
        "key_findings": "Cord compression; T2 intramedullary hyperintensity (myelomalacia); T1 hypointensity (poor prognosis); snake-eye appearance",
        "clinical_relevance": "Surgical decompression timing; modified Japanese Orthopaedic Association (mJOA) score correlation",
    },
    {
        "protocol_id": "NI-034",
        "name": "Whole-Body MRI for Neurological Metastases",
        "modality": "MRI",
        "domain": "neuro_oncology",
        "sequences": ["brain MRI with Gd", "whole-spine MRI", "diffusion whole-body"],
        "indication": "Staging and surveillance of CNS metastases",
        "key_findings": "Enhancing parenchymal lesions; leptomeningeal enhancement; vertebral bone marrow metastases",
        "clinical_relevance": "Guides SRS vs WBRT decisions; RANO-BM criteria for response; identifies drop metastases",
    },
    {
        "protocol_id": "NI-035",
        "name": "Transcranial Doppler for Vasospasm",
        "modality": "Ultrasound",
        "domain": "cerebrovascular",
        "sequences": ["TCD bilateral MCA, ACA, PCA", "Lindegaard ratio"],
        "indication": "Post-SAH vasospasm monitoring",
        "key_findings": "MCA velocity >120 cm/s concerning, >200 cm/s severe vasospasm; Lindegaard ratio >3 suggests vasospasm",
        "clinical_relevance": "Daily monitoring days 3-14 post-SAH; guides induced hypertension and angioplasty timing",
    },
    {
        "protocol_id": "NI-036",
        "name": "High-Resolution Vessel Wall MRI",
        "modality": "MRI",
        "domain": "cerebrovascular",
        "sequences": ["3D T1 SPACE/CUBE pre/post-Gd", "T2 SPACE", "PD SPACE"],
        "indication": "Intracranial vasculopathy characterization",
        "key_findings": "Eccentric wall enhancement: atherosclerosis; concentric: vasculitis; smooth enhancement: reversible vasoconstriction",
        "clinical_relevance": "Distinguishes atherosclerotic from inflammatory and non-inflammatory vasculopathies; guides treatment",
    },
    {
        "protocol_id": "NI-037",
        "name": "MRI Neurography for Peripheral Nerves",
        "modality": "MRI",
        "domain": "neuromuscular",
        "sequences": ["3D STIR SPACE", "3D T2 SPACE", "DTI of nerves"],
        "indication": "Peripheral nerve pathology evaluation",
        "key_findings": "Nerve caliber changes; T2 hyperintensity; fascicular pattern; mass lesions (schwannoma vs neurofibroma)",
        "clinical_relevance": "Localizes entrapment; distinguishes inflammatory from compressive neuropathy; guides surgical planning",
    },
    {
        "protocol_id": "NI-038",
        "name": "Iron-Sensitive MRI for Neurodegeneration",
        "modality": "MRI",
        "domain": "parkinson",
        "sequences": ["QSM", "R2* mapping", "SWI"],
        "indication": "Iron deposition quantification in neurodegenerative diseases",
        "key_findings": "Increased nigral iron in PD; putaminal iron in MSA-P; dentate/red nucleus in PSP; globus pallidus in NBIA",
        "clinical_relevance": "Helps differentiate parkinsonian syndromes; NBIA phenotyping; potential disease progression biomarker",
    },
    {
        "protocol_id": "NI-039",
        "name": "Perfusion MRI for Tumefactive MS vs Tumor",
        "modality": "MRI",
        "domain": "ms",
        "sequences": ["DSC perfusion", "conventional MRI", "MR spectroscopy"],
        "indication": "Distinguishing tumefactive demyelination from neoplasm",
        "key_findings": "Low rCBV in demyelination (vs elevated in high-grade tumor); incomplete ring enhancement; leading edge of restricted diffusion",
        "clinical_relevance": "Avoids unnecessary biopsy; MRS shows elevated Cho/Cr in both but lipid/lactate pattern differs",
    },
    {
        "protocol_id": "NI-040",
        "name": "CT Angiography for Intracranial Aneurysm",
        "modality": "CT",
        "domain": "cerebrovascular",
        "sequences": ["CTA with 3D reconstruction", "bone subtraction"],
        "indication": "Aneurysm detection in SAH or screening",
        "key_findings": "Aneurysm location, size, morphology; PHASES score components; parent vessel relationship",
        "clinical_relevance": "Sensitivity >95% for aneurysms >3mm; guides coiling vs clipping decision; 3D printing for planning",
    },
    # --- Additional clinical imaging findings ---
    {
        "protocol_id": "NI-041",
        "name": "MRI Findings in Status Epilepticus",
        "modality": "MRI",
        "domain": "epilepsy",
        "sequences": ["DWI", "FLAIR", "T2", "post-Gd T1"],
        "indication": "Post-status epilepticus brain injury evaluation",
        "key_findings": "Cortical DWI restriction (cytotoxic edema); hippocampal swelling/signal change; pulvinar sign; crossed cerebellar diaschisis",
        "clinical_relevance": "Predicts outcome; permanent hippocampal sclerosis if prolonged; guides ongoing seizure management",
    },
    {
        "protocol_id": "NI-042",
        "name": "MRI for Moyamoya Disease",
        "modality": "MRI",
        "domain": "cerebrovascular",
        "sequences": ["TOF MRA", "FLAIR", "DWI", "ASL perfusion", "SWI"],
        "indication": "Evaluation of suspected moyamoya disease or syndrome",
        "key_findings": "Bilateral ICA terminal stenosis/occlusion; basal moyamoya vessels; ivy sign on FLAIR; ASL territorial perfusion asymmetry",
        "clinical_relevance": "Suzuki staging; revascularization planning (STA-MCA bypass vs EDAS); perfusion reserve assessment",
    },
    {
        "protocol_id": "NI-043",
        "name": "MRI Spectroscopy for Metabolic Leukodystrophy",
        "modality": "MRI",
        "domain": "neuromuscular",
        "sequences": ["single-voxel MRS (short and long TE)", "conventional MRI"],
        "indication": "Metabolic characterization of white matter disease",
        "key_findings": "NAA reduction (axonal loss); elevated myo-inositol (gliosis); galactitol peak (Krabbe); NAA peak in Canavan disease",
        "clinical_relevance": "Non-invasive metabolic phenotyping; guides genetic testing; monitors treatment response in leukodystrophies",
    },
    {
        "protocol_id": "NI-044",
        "name": "Trigeminal Neuralgia MRI Protocol",
        "modality": "MRI",
        "domain": "headache",
        "sequences": ["3D CISS/FIESTA", "3D TOF MRA", "post-Gd T1"],
        "indication": "Evaluation of trigeminal neuralgia for neurovascular conflict",
        "key_findings": "Vascular loop compressing trigeminal nerve root entry zone; demyelinating plaque (MS-related); CPA mass lesion",
        "clinical_relevance": "Identifies surgical candidates for microvascular decompression; distinguishes classic from secondary TN",
    },
    {
        "protocol_id": "NI-045",
        "name": "MRI for Wilson Disease",
        "modality": "MRI",
        "domain": "movement",
        "sequences": ["T1", "T2", "FLAIR", "SWI"],
        "indication": "CNS involvement in Wilson disease evaluation",
        "key_findings": "T2/FLAIR hyperintensity in putamen, caudate, thalamus; face of giant panda sign in midbrain; T1 pallidal hyperintensity",
        "clinical_relevance": "Distinguishes hepatolenticular degeneration from other movement disorders; monitors treatment response to chelation",
    },
    {
        "protocol_id": "NI-046",
        "name": "Gadolinium-Enhanced MRI for Meningitis",
        "modality": "MRI",
        "domain": "cerebrovascular",
        "sequences": ["post-Gd T1", "FLAIR", "DWI"],
        "indication": "Suspected meningitis/encephalitis complications",
        "key_findings": "Leptomeningeal enhancement; subdural empyema; cerebritis with DWI restriction; ventriculitis; hydrocephalus",
        "clinical_relevance": "Identifies complications requiring neurosurgical intervention; abscess DWI restriction distinguishes from tumor",
    },
    {
        "protocol_id": "NI-047",
        "name": "Cavernous Sinus MRI Protocol",
        "modality": "MRI",
        "domain": "headache",
        "sequences": ["coronal T2", "coronal post-Gd T1 fat-sat", "MRA"],
        "indication": "Evaluation of cavernous sinus pathology with cranial neuropathy",
        "key_findings": "Tolosa-Hunt granulomatous inflammation; cavernous sinus thrombosis; meningioma; pituitary apoplexy",
        "clinical_relevance": "Guides treatment of painful ophthalmoplegia; distinguishes inflammatory from neoplastic causes",
    },
    {
        "protocol_id": "NI-048",
        "name": "DaTscan Patterns in Atypical Parkinsonism",
        "modality": "SPECT",
        "domain": "parkinson",
        "sequences": ["123I-ioflupane SPECT with semi-quantitative analysis"],
        "indication": "Differential diagnosis of parkinsonian syndromes",
        "key_findings": "PD: asymmetric putaminal loss (eagle wing); MSA-P: symmetric loss; PSP: symmetric with caudate involvement; DLB: similar to PD",
        "clinical_relevance": "Cannot distinguish PD from MSA/PSP/DLB; useful to differentiate degenerative from non-degenerative tremor",
    },
    {
        "protocol_id": "NI-049",
        "name": "MRI for Chiari Malformation",
        "modality": "MRI",
        "domain": "headache",
        "sequences": ["sagittal T1", "sagittal T2", "cine phase-contrast CSF flow"],
        "indication": "Evaluation of Chiari malformation with headache and syringomyelia",
        "key_findings": "Tonsillar herniation >5mm below foramen magnum; syringomyelia; absent posterior CSF flow on cine study",
        "clinical_relevance": "Abnormal CSF flow dynamics predict surgical benefit from posterior fossa decompression",
    },
    {
        "protocol_id": "NI-050",
        "name": "ASL Perfusion for Seizure Lateralization",
        "modality": "MRI",
        "domain": "epilepsy",
        "sequences": ["pulsed or pseudo-continuous ASL", "co-registered with structural"],
        "indication": "Non-invasive perfusion assessment for epilepsy lateralization",
        "key_findings": "Interictal hypoperfusion in seizure onset zone; postictal hyperperfusion if recent seizure",
        "clinical_relevance": "Non-contrast alternative to SPECT; can be added to routine epilepsy MRI protocol; supports lateralization",
    },
    {
        "protocol_id": "NI-051",
        "name": "Whole-Brain Volumetry for Neurodegenerative Monitoring",
        "modality": "MRI",
        "domain": "dementia",
        "sequences": ["3D T1 MPRAGE", "automated volumetric analysis"],
        "indication": "Longitudinal brain volume monitoring in neurodegeneration",
        "key_findings": "Hippocampal volume percentile; annualized whole-brain atrophy rate; regional cortical thickness maps",
        "clinical_relevance": "ARIA monitoring in anti-amyloid therapy; atrophy rate >1.5%/year concerning; NeuroQuant/BrainSuite quantification",
    },
    {
        "protocol_id": "NI-052",
        "name": "MRI for Spinal Dural AV Fistula",
        "modality": "MRI",
        "domain": "neuromuscular",
        "sequences": ["sagittal T2", "post-Gd T1", "MR angiography of spine"],
        "indication": "Suspected spinal dural AV fistula with progressive myelopathy",
        "key_findings": "Flow voids along spinal cord surface; cord T2 hyperintensity; centromedullary enhancement; dilated perimedullary veins",
        "clinical_relevance": "Often misdiagnosed as MS or spondylotic myelopathy; early treatment (embolization/surgery) prevents irreversible cord damage",
    },
    {
        "protocol_id": "NI-053",
        "name": "MRI FLAIR Hyperintense Vessels in Stroke",
        "modality": "MRI",
        "domain": "stroke",
        "sequences": ["FLAIR", "DWI", "MRA"],
        "indication": "Assessment of collateral flow in acute ischemic stroke",
        "key_findings": "Hyperintense vessel sign on FLAIR distal to arterial occlusion; indicates slow flow in leptomeningeal collaterals",
        "clinical_relevance": "Positive sign correlates with salvageable tissue; prognostic marker for reperfusion therapy benefit",
    },
    {
        "protocol_id": "NI-054",
        "name": "Posterior Reversible Encephalopathy Syndrome (PRES) MRI",
        "modality": "MRI",
        "domain": "cerebrovascular",
        "sequences": ["FLAIR", "DWI", "ADC", "post-Gd T1", "SWI"],
        "indication": "Suspected PRES in hypertensive emergency, eclampsia, immunosuppression",
        "key_findings": "Bilateral posterior white matter FLAIR/T2 hyperintensity; vasogenic edema (elevated ADC); hemorrhage on SWI; rare DWI restriction (poor prognosis)",
        "clinical_relevance": "Reversible with blood pressure control and offending agent removal; cytotoxic edema predicts permanent injury",
    },
    {
        "protocol_id": "NI-055",
        "name": "Neonatal Brain MRI Protocol",
        "modality": "MRI",
        "domain": "epilepsy",
        "sequences": ["T1 MPRAGE", "T2", "DWI", "SWI", "MRS"],
        "indication": "Neonatal seizures, HIE, congenital brain malformations",
        "key_findings": "HIE patterns on DWI (watershed vs basal ganglia/thalamus); cortical malformations; myelination assessment; MRS lactate/NAA ratio",
        "clinical_relevance": "Timing of MRI 3-5 days post-HIE for prognostication; guides therapeutic hypothermia response assessment",
    },
    # --- Expanded Neuroimaging Protocols (15 new entries) ---
    {
        "protocol_id": "NI-056",
        "name": "DWI Acute Stroke Pattern Analysis",
        "modality": "MRI",
        "domain": "stroke",
        "sequences": ["DWI", "ADC map", "trace-weighted imaging"],
        "indication": "Acute ischemic stroke diffusion restriction pattern characterization",
        "key_findings": "Territorial infarct patterns (MCA, PCA, ACA, vertebrobasilar); lacunar DWI patterns; watershed distribution; scattered embolic pattern; cortical vs subcortical restriction",
        "clinical_relevance": "DWI pattern determines stroke mechanism (large vessel, small vessel, embolic, hemodynamic); guides acute treatment and secondary prevention strategy",
    },
    {
        "protocol_id": "NI-057",
        "name": "FLAIR MS Lesion Characterization",
        "modality": "MRI",
        "domain": "ms",
        "sequences": ["3D FLAIR", "sagittal FLAIR", "coronal FLAIR"],
        "indication": "MS lesion detection and classification by anatomical distribution",
        "key_findings": "Periventricular lesions (Dawson fingers perpendicular to ventricles); juxtacortical/cortical lesions; infratentorial lesions; corpus callosum involvement (sagittal FLAIR)",
        "clinical_relevance": "McDonald 2017 DIS requires lesions in 2+ of 4 locations; FLAIR is most sensitive for periventricular and juxtacortical lesions; lesion volume correlates with disability",
    },
    {
        "protocol_id": "NI-058",
        "name": "SWI Microbleed Mapping",
        "modality": "MRI",
        "domain": "cerebrovascular",
        "sequences": ["SWI", "T2*-GRE", "QSM"],
        "indication": "Cerebral microbleed detection and distribution mapping for hemorrhagic risk stratification",
        "key_findings": "Lobar distribution (cerebral amyloid angiopathy); deep/infratentorial distribution (hypertensive microangiopathy); mixed pattern; cortical superficial siderosis; microbleed count and burden score",
        "clinical_relevance": "Boston criteria v2.0 for CAA diagnosis; microbleed burden informs anticoagulation risk-benefit; >5 lobar microbleeds associated with increased ICH risk on anticoagulation",
    },
    {
        "protocol_id": "NI-059",
        "name": "DTI White Matter Tractography",
        "modality": "MRI",
        "domain": "neuromuscular",
        "sequences": ["multi-direction DTI (30+ directions)", "fractional anisotropy maps", "mean diffusivity maps", "probabilistic tractography"],
        "indication": "White matter tract integrity assessment in neurodegeneration and demyelination",
        "key_findings": "Reduced FA in corticospinal tract (ALS); reduced FA in corpus callosum (MS); altered connectivity in default mode network (AD); tract-specific injury patterns",
        "clinical_relevance": "Quantitative biomarker for white matter degeneration; early detection of tract-specific injury before conventional MRI changes; research biomarker for disease progression",
    },
    {
        "protocol_id": "NI-060",
        "name": "ASL Perfusion Without Contrast",
        "modality": "MRI",
        "domain": "cerebrovascular",
        "sequences": ["pseudo-continuous ASL (pCASL)", "3D background-suppressed acquisition"],
        "indication": "Non-invasive cerebral perfusion assessment without gadolinium contrast",
        "key_findings": "Territorial perfusion asymmetry in stroke; hypoperfusion in epilepsy focus (interictal); hyperperfusion in luxury perfusion post-stroke; global hypoperfusion patterns in dementia",
        "clinical_relevance": "Contrast-free perfusion alternative for patients with renal insufficiency; repeatable for longitudinal monitoring; added to routine epilepsy and dementia protocols",
    },
    {
        "protocol_id": "NI-061",
        "name": "fMRI Presurgical Functional Mapping",
        "modality": "MRI",
        "domain": "neuro_oncology",
        "sequences": ["BOLD fMRI motor paradigm", "BOLD fMRI language paradigm", "resting-state fMRI"],
        "indication": "Presurgical functional cortex localization for tumor and epilepsy surgery planning",
        "key_findings": "Primary motor cortex activation relative to tumor; Broca and Wernicke area lateralization; hemispheric dominance index; functional reorganization in slow-growing tumors",
        "clinical_relevance": "Replaces Wada test for language lateralization; guides extent of resection to preserve eloquent cortex; combined with DTI tractography for comprehensive surgical planning",
    },
    {
        "protocol_id": "NI-062",
        "name": "CT Angiography Circle of Willis for LVO Detection",
        "modality": "CT",
        "domain": "stroke",
        "sequences": ["CTA with bolus tracking", "thin-slice multiplanar reformats", "MIP reconstructions"],
        "indication": "Large vessel occlusion detection in acute ischemic stroke for thrombectomy triage",
        "key_findings": "ICA terminus occlusion; M1/M2 MCA occlusion; basilar artery occlusion; tandem lesions; collateral grade assessment (Tan score, CS-ASPECTS)",
        "clinical_relevance": "Sensitivity >95% for proximal LVO detection; collateral status predicts thrombectomy outcome; rapid acquisition enables field-to-cath-lab triage",
    },
    {
        "protocol_id": "NI-063",
        "name": "CT Perfusion Core-Penumbra Mismatch",
        "modality": "CT",
        "domain": "stroke",
        "sequences": ["CTP with automated RAPID/Viz.ai processing", "CBF maps", "CBV maps", "Tmax maps"],
        "indication": "Ischemic core and penumbra quantification for extended-window thrombectomy eligibility",
        "key_findings": "Ischemic core (CBF <30% or rCBF <30%); critically hypoperfused tissue (Tmax >6s); mismatch ratio (penumbra/core); mismatch volume; DAWN and DEFUSE-3 eligibility criteria",
        "clinical_relevance": "Automated processing enables 6-24h thrombectomy selection; mismatch ratio >1.8 and core <70 mL (DEFUSE-3) or clinical-core mismatch (DAWN) guide intervention",
    },
    {
        "protocol_id": "NI-064",
        "name": "MR Spectroscopy for Tumor Grading",
        "modality": "MRI",
        "domain": "neuro_oncology",
        "sequences": ["single-voxel MRS (TE 135ms)", "single-voxel MRS (TE 30ms)", "multivoxel CSI"],
        "indication": "Non-invasive metabolic characterization and grading of brain tumors",
        "key_findings": "NAA/Cho ratio reduction in tumor (neuronal loss); elevated Cho/Cr (membrane turnover); 2-HG peak at 2.25 ppm in IDH-mutant gliomas; lipid-lactate peaks in high-grade necrosis; myo-inositol elevation in low-grade gliomas",
        "clinical_relevance": "2-HG detection >90% sensitive for IDH mutation without biopsy; Cho/NAA ratio correlates with WHO grade; distinguishes recurrence from radiation necrosis; guides biopsy targeting",
    },
    {
        "protocol_id": "NI-065",
        "name": "MIBG Cardiac Scintigraphy for Autonomic Neuropathy",
        "modality": "SPECT",
        "domain": "parkinson",
        "sequences": ["123I-MIBG planar imaging", "early (15 min) and delayed (4 h) heart-to-mediastinum ratio"],
        "indication": "Cardiac sympathetic denervation assessment in Lewy body disorders vs other parkinsonisms",
        "key_findings": "Low H/M ratio (<1.6 delayed) in PD and DLB (cardiac sympathetic denervation); preserved H/M ratio in MSA, PSP, CBD; washout rate >43% abnormal",
        "clinical_relevance": "Distinguishes PD/DLB from MSA/PSP with high specificity; reflects peripheral autonomic degeneration; supportive diagnostic criterion for DLB; not affected by dopaminergic medications",
    },
    {
        "protocol_id": "NI-066",
        "name": "Amyloid PET with Florbetapir",
        "modality": "PET",
        "domain": "dementia",
        "sequences": ["18F-florbetapir (Amyvid)", "18F-florbetaben (Neuraceq)", "11C-PiB (Pittsburgh Compound B)"],
        "indication": "In vivo amyloid-beta plaque quantification for Alzheimer's diagnosis and treatment eligibility",
        "key_findings": "Positive scan: diffuse cortical retention in frontal, parietal, temporal, precuneus; negative scan: white matter only uptake; Centiloid scale for quantification; SUVR >1.11 positive",
        "clinical_relevance": "Required for anti-amyloid therapy eligibility (lecanemab, donanemab); AUC 2.0 appropriate use criteria; negative scan essentially excludes AD pathology; serial imaging tracks amyloid clearance",
    },
    {
        "protocol_id": "NI-067",
        "name": "Tau PET with Flortaucipir",
        "modality": "PET",
        "domain": "dementia",
        "sequences": ["18F-flortaucipir (Tauvid)", "18F-MK-6240", "18F-PI-2620"],
        "indication": "In vivo tau neurofibrillary tangle distribution mapping",
        "key_findings": "Braak stage-concordant tau distribution; entorhinal/hippocampal uptake in early AD; widespread neocortical uptake in advanced AD; frontal predominance in behavioral variant FTD; off-target binding in basal ganglia with flortaucipir",
        "clinical_relevance": "Correlates with cognitive decline better than amyloid PET; supports ATN T+ classification; staging in vivo corresponds to Braak stages; emerging endpoint in clinical trials",
    },
    {
        "protocol_id": "NI-068",
        "name": "FDG-PET Brain Hypometabolism Patterns",
        "modality": "PET",
        "domain": "dementia",
        "sequences": ["18F-FDG PET brain"],
        "indication": "Metabolic pattern analysis for dementia differential diagnosis and neurodegeneration assessment",
        "key_findings": "AD pattern: bilateral temporoparietal and posterior cingulate hypometabolism; FTD pattern: frontal and anterior temporal hypometabolism; DLB: occipital hypometabolism with cingulate island sign; CBD: asymmetric frontoparietal",
        "clinical_relevance": "Sensitivity >85% for AD vs FTD vs DLB differential; normal metabolism argues against neurodegeneration; hypometabolism pattern guides clinical diagnosis when biomarkers are unavailable",
    },
    {
        "protocol_id": "NI-069",
        "name": "SPECT Ictal vs Interictal for Epilepsy Localization",
        "modality": "SPECT",
        "domain": "epilepsy",
        "sequences": ["ictal 99mTc-HMPAO SPECT", "interictal 99mTc-HMPAO SPECT", "SISCOM subtraction co-registered to MRI"],
        "indication": "Seizure focus localization in MRI-negative or discordant epilepsy surgery evaluation",
        "key_findings": "Ictal hyperperfusion at seizure onset zone; interictal hypoperfusion at epileptogenic region; SISCOM subtraction highlights ictal-interictal difference co-registered to structural MRI",
        "clinical_relevance": "Injection timing critical (<30 seconds from seizure onset for optimal localization); concordance with EEG and MRI improves surgical outcome; SISCOM increases sensitivity of SPECT localization by 20-30%",
    },
    {
        "protocol_id": "NI-070",
        "name": "Whole Spine MRI Myelopathy Screening Protocol",
        "modality": "MRI",
        "domain": "neuromuscular",
        "sequences": ["sagittal T2 whole spine", "sagittal T1 whole spine", "sagittal STIR", "axial T2 at levels of interest", "post-Gd T1 if indicated"],
        "indication": "Comprehensive myelopathy screening for spinal cord compression, demyelination, or intrinsic cord lesions",
        "key_findings": "Cervical spondylotic myelopathy (cord compression with T2 signal change); longitudinally extensive transverse myelitis (>3 segments in NMOSD); syringomyelia; spinal cord tumors; epidural abscess/hematoma",
        "clinical_relevance": "Distinguishes compressive from non-compressive myelopathy; NMOSD vs MS myelitis differentiation by lesion length; urgent surgical evaluation for progressive cord compression with myelopathy signs",
    },
]


# ===================================================================
# NEUROIMAGING PARSER IMPLEMENTATION
# ===================================================================


class NeuroimagingParser(BaseIngestParser):
    """Parse neuroimaging protocols and findings for the Neurology Intelligence Agent.

    In seed mode, returns the curated NEUROIMAGING_PROTOCOLS list.

    Usage::

        parser = NeuroimagingParser()
        records, stats = parser.run()
    """

    def __init__(
        self,
        collection_manager: Any = None,
        embedder: Any = None,
    ) -> None:
        super().__init__(
            source_name="neuroimaging",
            collection_manager=collection_manager,
            embedder=embedder,
        )

    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch neuroimaging protocol data.

        Returns:
            List of raw neuroimaging protocol dictionaries.
        """
        self.logger.info(
            "Using curated neuroimaging seed data (%d protocols)",
            len(NEUROIMAGING_PROTOCOLS),
        )
        return list(NEUROIMAGING_PROTOCOLS)

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw neuroimaging data into IngestRecord objects.

        Args:
            raw_data: List of neuroimaging protocol dictionaries.

        Returns:
            List of IngestRecord objects.
        """
        records: List[IngestRecord] = []

        for entry in raw_data:
            protocol_id = entry.get("protocol_id", "")
            name = entry.get("name", "")
            modality = entry.get("modality", "")
            domain = entry.get("domain", "")
            sequences = entry.get("sequences", [])
            indication = entry.get("indication", "")
            key_findings = entry.get("key_findings", "")
            clinical_relevance = entry.get("clinical_relevance", "")

            sequences_str = ", ".join(sequences) if sequences else "not specified"
            text = (
                f"Neuroimaging Protocol: {name} ({modality}). "
                f"Domain: {domain}. "
                f"Sequences: {sequences_str}. "
                f"Indication: {indication}. "
                f"Key findings: {key_findings}. "
                f"Clinical relevance: {clinical_relevance}."
            )

            record = IngestRecord(
                text=text,
                metadata={
                    "protocol_id": protocol_id,
                    "name": name,
                    "modality": modality,
                    "domain": domain,
                    "sequences": sequences,
                    "indication": indication,
                    "source_db": "neuroimaging_seed",
                },
                collection_name="neuro_imaging",
                record_id=protocol_id,
                source="neuroimaging",
            )
            records.append(record)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate a neuroimaging IngestRecord.

        Requirements:
            - text must be non-empty
            - must have protocol_id in metadata
            - must have name in metadata
            - must have modality in metadata

        Args:
            record: The record to validate.

        Returns:
            True if the record passes all validation checks.
        """
        if not record.text or not record.text.strip():
            return False

        meta = record.metadata
        if not meta.get("protocol_id"):
            return False
        if not meta.get("name"):
            return False
        if not meta.get("modality"):
            return False

        return True


def get_neuroimaging_protocol_count() -> int:
    """Return the number of curated neuroimaging protocols."""
    return len(NEUROIMAGING_PROTOCOLS)


def get_imaging_modalities() -> List[str]:
    """Return a deduplicated sorted list of imaging modalities."""
    modalities = list({p["modality"] for p in NEUROIMAGING_PROTOCOLS})
    modalities.sort()
    return modalities


def get_imaging_domains() -> List[str]:
    """Return a deduplicated sorted list of imaging domains."""
    domains = list({p["domain"] for p in NEUROIMAGING_PROTOCOLS})
    domains.sort()
    return domains
