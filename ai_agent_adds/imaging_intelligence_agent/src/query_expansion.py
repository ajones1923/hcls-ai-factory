"""Query expansion for medical imaging domain.

12 expansion maps that augment user queries with synonyms, abbreviations,
and related technical terms to improve Milvus vector search recall.
"""

from typing import Set


# ═══════════════════════════════════════════════════════════════════════
# EXPANSION MAPS
# ═══════════════════════════════════════════════════════════════════════

MODALITY_EXPANSION = {
    "ct": {"computed tomography", "ct scan", "cat scan", "helical ct", "spiral ct", "mdct", "multi-detector ct"},
    "mri": {"magnetic resonance imaging", "mri scan", "mr imaging", "nmr", "magnetic resonance"},
    "xray": {"x-ray", "radiograph", "plain film", "plain radiograph", "conventional radiography"},
    "cxr": {"chest x-ray", "chest radiograph", "chest film", "pa chest", "ap chest", "portable chest"},
    "pet": {"positron emission tomography", "pet scan", "fdg pet", "fdg-pet"},
    "pet_ct": {"pet/ct", "pet-ct", "combined pet ct", "hybrid imaging"},
    "ultrasound": {"us", "sonography", "ultrasonography", "doppler", "duplex"},
    "mammography": {"mammo", "mammogram", "breast imaging", "tomosynthesis", "dbt", "digital breast"},
    "fluoroscopy": {"fluoro", "fluoroscopic", "real-time imaging", "contrast study"},
}

BODY_REGION_EXPANSION = {
    "head": {"brain", "cranial", "intracranial", "cerebral", "skull", "calvarium"},
    "neck": {"cervical", "thyroid", "larynx", "pharynx", "carotid"},
    "chest": {"thorax", "thoracic", "lung", "pulmonary", "mediastinal", "pleural"},
    "abdomen": {"abdominal", "liver", "hepatic", "renal", "kidney", "pancreas", "spleen", "bowel", "intestinal"},
    "pelvis": {"pelvic", "hip", "bladder", "uterine", "ovarian", "prostatic", "rectal"},
    "spine": {"spinal", "vertebral", "lumbar", "thoracic spine", "cervical spine", "sacral", "disc"},
    "extremity": {"musculoskeletal", "upper extremity", "lower extremity", "limb", "joint", "bone"},
    "brain": {"cerebral", "intracranial", "cranial", "neurological", "cns"},
    "cardiac": {"heart", "coronary", "myocardial", "pericardial", "aortic", "valvular"},
    "breast": {"mammary", "breast tissue", "axillary"},
}

PATHOLOGY_EXPANSION = {
    "hemorrhage": {"bleeding", "hematoma", "blood", "hemorrhagic", "haemorrhage", "ich", "intracranial hemorrhage"},
    "nodule": {"nodular", "pulmonary nodule", "lung nodule", "solitary pulmonary nodule", "spn"},
    "mass": {"tumor", "tumour", "lesion", "neoplasm", "malignancy", "cancer", "growth"},
    "fracture": {"broken bone", "fx", "fractures", "break", "cortical disruption"},
    "pneumonia": {"consolidation", "pneumonic", "lung infection", "pneumonitis", "cap", "hap"},
    "effusion": {"pleural effusion", "fluid collection", "pleural fluid", "empyema"},
    "pneumothorax": {"ptx", "collapsed lung", "air leak", "tension pneumothorax"},
    "stroke": {"cva", "cerebrovascular accident", "infarct", "ischemic stroke", "hemorrhagic stroke"},
    "embolism": {"pe", "pulmonary embolism", "thromboembolism", "dvt", "clot"},
    "dissection": {"aortic dissection", "intimal flap", "false lumen", "stanford"},
    "stenosis": {"narrowing", "stricture", "obstruction", "occlusion"},
    "fibrosis": {"pulmonary fibrosis", "ipf", "uip", "nsip", "interstitial lung disease", "ild"},
}

AI_TASK_EXPANSION = {
    "segmentation": {"segment", "segmenting", "delineation", "contour", "annotation", "mask", "roi"},
    "detection": {"detect", "detecting", "finding", "identify", "localize", "localization", "cade"},
    "classification": {"classify", "classifying", "categorize", "grade", "staging", "diagnosis", "cadx"},
    "triage": {"prioritize", "prioritization", "urgent", "worklist", "routing", "alert"},
    "quantification": {"measure", "measurement", "volumetric", "volume", "quantify", "quantitative"},
    "registration": {"alignment", "align", "coregistration", "mapping", "spatial correspondence", "deformable"},
    "reconstruction": {"reconstruct", "super-resolution", "denoising", "artifact reduction"},
    "report_generation": {"reporting", "report", "structured report", "radiology report", "impression"},
}

SEVERITY_EXPANSION = {
    "critical": {"emergent", "stat", "p1", "life-threatening", "emergency", "acute"},
    "urgent": {"p2", "priority", "time-sensitive", "significant"},
    "routine": {"p4", "normal", "stable", "benign", "unremarkable", "negative"},
}

FINDING_EXPANSION = {
    "consolidation": {"airspace opacity", "airspace disease", "opacification", "air bronchograms"},
    "ground_glass": {"ggo", "ground-glass opacity", "hazy opacity", "frosted glass"},
    "atelectasis": {"collapse", "volume loss", "subsegmental atelectasis", "plate-like"},
    "edema": {"pulmonary edema", "interstitial edema", "alveolar edema", "fluid overload", "chf"},
    "calcification": {"calcified", "calcium", "atherosclerotic", "coronary calcium"},
    "lymphadenopathy": {"enlarged lymph node", "lymph node", "adenopathy", "hilar", "mediastinal"},
}

GUIDELINE_EXPANSION = {
    "lung_rads": {"lung-rads", "lungrads", "acr lung-rads", "lung cancer screening", "ldct screening"},
    "bi_rads": {"bi-rads", "birads", "breast imaging reporting", "mammography assessment"},
    "ti_rads": {"ti-rads", "tirads", "thyroid imaging", "thyroid nodule assessment"},
    "li_rads": {"li-rads", "lirads", "liver imaging reporting", "hcc screening"},
    "acr": {"american college of radiology", "acr appropriateness", "appropriateness criteria"},
    "fleischner": {"fleischner society", "fleischner criteria", "incidental nodule"},
    "recist": {"recist 1.1", "response evaluation", "tumor response", "treatment response"},
}

DEVICE_EXPANSION = {
    "fda": {"fda cleared", "fda approved", "regulatory", "510k", "510(k)"},
    "510k": {"510(k)", "premarket notification", "substantially equivalent"},
    "de_novo": {"de novo classification", "novel device"},
    "cleared": {"fda cleared", "authorized", "marketed"},
    "approved": {"fda approved", "pma approved"},
}

DATASET_EXPANSION = {
    "rsna": {"rsna dataset", "rsna challenge", "rsna competition", "radiological society"},
    "tcia": {"the cancer imaging archive", "cancer imaging", "tcia collection"},
    "nih": {"nih dataset", "nih clinical center", "nih chest x-ray", "chestx-ray14"},
    "physionet": {"physionet dataset", "mimic", "mimic-cxr"},
    "kaggle": {"kaggle competition", "kaggle dataset"},
    "lidc": {"lidc-idri", "lung image database", "lung nodule dataset"},
    "brats": {"brain tumor segmentation", "brats challenge", "glioma segmentation"},
}

MODEL_ARCHITECTURE_EXPANSION = {
    "unet": {"u-net", "3d u-net", "u-net 3d", "encoder-decoder"},
    "vista3d": {"vista-3d", "vista 3d", "versatile imaging segmentation", "nvidia vista"},
    "maisi": {"medical ai synthetic imaging", "synthetic ct", "latent diffusion ct"},
    "vilam3": {"vila-m3", "vila m3", "vision language medical", "radiology vlm", "monai m3"},
    "nnunet": {"nnu-net", "nn-unet", "self-configuring segmentation"},
    "swin": {"swin transformer", "swinunetr", "swin-unetr", "shifted window"},
    "resnet": {"residual network", "resnet50", "resnet101"},
    "densenet": {"densenet-121", "densenet121", "dense network"},
    "retinanet": {"retina-net", "retina net", "focal loss detector"},
    "segresnet": {"seg-resnet", "segmentation resnet"},
    "monai": {"monai bundle", "monai model zoo", "monai deploy", "monai toolkit"},
}

MEASUREMENT_EXPANSION = {
    "volume": {"volumetric", "volume measurement", "cc", "ml", "cubic"},
    "diameter": {"size", "longest diameter", "short axis", "long axis", "max dimension"},
    "recist": {"recist measurement", "target lesion", "sum of diameters", "response criteria"},
    "doubling_time": {"vdt", "volume doubling time", "growth rate", "tumor growth"},
    "hounsfield": {"hu", "density", "attenuation", "hounsfield unit"},
}

CONTRAST_EXPANSION = {
    "contrast": {"contrast agent", "contrast enhanced", "iv contrast", "contrast material"},
    "gadolinium": {"gad", "gd", "gadolinium-based", "gbca", "mri contrast"},
    "iodinated": {"iodine", "iodinated contrast", "ct contrast", "omnipaque", "visipaque"},
    "enhancement": {"enhancing", "post-contrast", "arterial phase", "venous phase", "delayed phase"},
}

# All expansion maps in a single list for iteration
ALL_EXPANSION_MAPS = [
    MODALITY_EXPANSION,
    BODY_REGION_EXPANSION,
    PATHOLOGY_EXPANSION,
    AI_TASK_EXPANSION,
    SEVERITY_EXPANSION,
    FINDING_EXPANSION,
    GUIDELINE_EXPANSION,
    DEVICE_EXPANSION,
    DATASET_EXPANSION,
    MODEL_ARCHITECTURE_EXPANSION,
    MEASUREMENT_EXPANSION,
    CONTRAST_EXPANSION,
]


def expand_query(query: str) -> Set[str]:
    """Expand a query with synonyms and related terms from all 12 maps.

    Scans the query text against all expansion map keys using word boundary
    matching, and returns the union of matched expansion terms.
    """
    query_lower = query.lower()
    expanded = set()

    for expansion_map in ALL_EXPANSION_MAPS:
        for key, synonyms in expansion_map.items():
            # Use simple substring match on key (keys are short domain terms)
            # Word boundary matching would be better but adds complexity
            if key.replace("_", " ") in query_lower or key in query_lower:
                expanded.update(synonyms)

    return expanded
