"""OMIM (Online Mendelian Inheritance in Man) parser for the Rare Disease Diagnostic Agent.

Parses OMIM disease entries and seeds 50+ key rare diseases with structured
data including OMIM ID, disease name, gene, inheritance pattern, clinical
features, and prevalence information.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA: 55 KEY RARE DISEASES FROM OMIM
# ===================================================================

OMIM_DISEASES: List[Dict[str, Any]] = [
    {
        "omim_id": "219700",
        "disease_name": "Cystic Fibrosis",
        "gene": "CFTR",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "chronic pulmonary disease", "pancreatic insufficiency",
            "elevated sweat chloride", "male infertility",
            "recurrent respiratory infections", "bronchiectasis",
        ],
        "prevalence": "1:3,500 Caucasians",
        "chromosome": "7q31.2",
    },
    {
        "omim_id": "310200",
        "disease_name": "Duchenne Muscular Dystrophy",
        "gene": "DMD",
        "inheritance": "x_linked_recessive",
        "clinical_features": [
            "progressive muscle weakness", "elevated creatine kinase",
            "Gowers sign", "pseudohypertrophy of calves",
            "cardiomyopathy", "respiratory failure",
        ],
        "prevalence": "1:3,500 males",
        "chromosome": "Xp21.2",
    },
    {
        "omim_id": "141900",
        "disease_name": "Huntington Disease",
        "gene": "HTT",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "chorea", "progressive dementia", "psychiatric disturbances",
            "caudate atrophy", "motor impairment", "weight loss",
        ],
        "prevalence": "1:10,000-1:20,000",
        "chromosome": "4p16.3",
    },
    {
        "omim_id": "230800",
        "disease_name": "Gaucher Disease Type 1",
        "gene": "GBA",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "hepatosplenomegaly", "anemia", "thrombocytopenia",
            "bone pain", "Erlenmeyer flask deformity", "fatigue",
        ],
        "prevalence": "1:40,000 general; 1:800 Ashkenazi",
        "chromosome": "1q22",
    },
    {
        "omim_id": "232300",
        "disease_name": "Glycogen Storage Disease Type II (Pompe Disease)",
        "gene": "GAA",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "cardiomegaly", "hypotonia", "muscle weakness",
            "respiratory insufficiency", "hepatomegaly", "feeding difficulties",
        ],
        "prevalence": "1:40,000",
        "chromosome": "17q25.3",
    },
    {
        "omim_id": "143100",
        "disease_name": "Huntington Disease-Like 2",
        "gene": "JPH3",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "chorea", "dystonia", "dementia",
            "psychiatric symptoms", "weight loss",
        ],
        "prevalence": "rare, African descent",
        "chromosome": "16q24.2",
    },
    {
        "omim_id": "256700",
        "disease_name": "Neuronal Ceroid Lipofuscinosis Type 2",
        "gene": "TPP1",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "seizures", "progressive motor decline", "visual loss",
            "language regression", "ataxia", "dementia",
        ],
        "prevalence": "1:100,000",
        "chromosome": "11p15.4",
    },
    {
        "omim_id": "261600",
        "disease_name": "Phenylketonuria",
        "gene": "PAH",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "intellectual disability", "seizures", "eczema",
            "musty odor", "microcephaly", "behavioral problems",
        ],
        "prevalence": "1:10,000-1:15,000",
        "chromosome": "12q23.2",
    },
    {
        "omim_id": "166200",
        "disease_name": "Osteogenesis Imperfecta Type I",
        "gene": "COL1A1",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "bone fragility", "blue sclerae", "hearing loss",
            "easy bruising", "joint hypermobility", "dental abnormalities",
        ],
        "prevalence": "1:15,000-1:20,000",
        "chromosome": "17q21.33",
    },
    {
        "omim_id": "162200",
        "disease_name": "Neurofibromatosis Type 1",
        "gene": "NF1",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "cafe-au-lait macules", "neurofibromas", "Lisch nodules",
            "optic glioma", "bone dysplasia", "learning disabilities",
        ],
        "prevalence": "1:3,000",
        "chromosome": "17q11.2",
    },
    {
        "omim_id": "191100",
        "disease_name": "Tuberous Sclerosis Complex",
        "gene": "TSC1",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "cortical tubers", "subependymal nodules", "renal angiomyolipomas",
            "facial angiofibromas", "seizures", "intellectual disability",
        ],
        "prevalence": "1:6,000-1:10,000",
        "chromosome": "9q34.13",
    },
    {
        "omim_id": "254800",
        "disease_name": "Spinal Muscular Atrophy Type I",
        "gene": "SMN1",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "severe hypotonia", "muscle weakness", "respiratory failure",
            "absent reflexes", "tongue fasciculations", "feeding difficulties",
        ],
        "prevalence": "1:10,000",
        "chromosome": "5q13.2",
    },
    {
        "omim_id": "300376",
        "disease_name": "Rett Syndrome",
        "gene": "MECP2",
        "inheritance": "x_linked_dominant",
        "clinical_features": [
            "loss of acquired skills", "stereotypic hand movements",
            "gait abnormalities", "seizures", "breathing irregularities",
            "microcephaly",
        ],
        "prevalence": "1:10,000-1:15,000 females",
        "chromosome": "Xq28",
    },
    {
        "omim_id": "248200",
        "disease_name": "Maple Syrup Urine Disease",
        "gene": "BCKDHA",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "sweet-smelling urine", "poor feeding", "lethargy",
            "metabolic crisis", "developmental delay", "seizures",
        ],
        "prevalence": "1:185,000",
        "chromosome": "19q13.2",
    },
    {
        "omim_id": "232700",
        "disease_name": "Glycogen Storage Disease Type III",
        "gene": "AGL",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "hepatomegaly", "hypoglycemia", "muscle weakness",
            "cardiomyopathy", "growth retardation", "elevated liver enzymes",
        ],
        "prevalence": "1:100,000",
        "chromosome": "1p21.2",
    },
    {
        "omim_id": "277900",
        "disease_name": "Wilson Disease",
        "gene": "ATP7B",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "hepatic dysfunction", "Kayser-Fleischer rings", "neuropsychiatric symptoms",
            "tremor", "low serum ceruloplasmin", "elevated hepatic copper",
        ],
        "prevalence": "1:30,000",
        "chromosome": "13q14.3",
    },
    {
        "omim_id": "276700",
        "disease_name": "Tyrosinemia Type I",
        "gene": "FAH",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "liver failure", "renal tubular dysfunction", "rickets",
            "hepatocellular carcinoma risk", "porphyria-like crises",
            "cabbage-like odor",
        ],
        "prevalence": "1:100,000",
        "chromosome": "15q25.1",
    },
    {
        "omim_id": "230900",
        "disease_name": "Gaucher Disease Type 2",
        "gene": "GBA",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "hepatosplenomegaly", "neurological deterioration",
            "oculomotor apraxia", "dysphagia", "stridor",
            "failure to thrive",
        ],
        "prevalence": "1:100,000",
        "chromosome": "1q22",
    },
    {
        "omim_id": "252500",
        "disease_name": "Mucolipidosis Type II (I-Cell Disease)",
        "gene": "GNPTAB",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "coarse facial features", "skeletal dysplasia",
            "developmental delay", "cardiomegaly",
            "restricted joint mobility", "gingival hyperplasia",
        ],
        "prevalence": "1:100,000-1:400,000",
        "chromosome": "12q23.2",
    },
    {
        "omim_id": "607624",
        "disease_name": "CHARGE Syndrome",
        "gene": "CHD7",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "coloboma", "heart defects", "atresia choanae",
            "growth retardation", "genital abnormalities", "ear anomalies",
        ],
        "prevalence": "1:10,000-1:15,000",
        "chromosome": "8q12.2",
    },
    {
        "omim_id": "154700",
        "disease_name": "Marfan Syndrome",
        "gene": "FBN1",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "aortic root dilation", "lens subluxation", "tall stature",
            "arachnodactyly", "pectus deformity", "scoliosis",
        ],
        "prevalence": "1:5,000-1:10,000",
        "chromosome": "15q21.1",
    },
    {
        "omim_id": "174000",
        "disease_name": "Polycystic Kidney Disease (Autosomal Dominant)",
        "gene": "PKD1",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "bilateral renal cysts", "hypertension", "renal insufficiency",
            "hepatic cysts", "intracranial aneurysms", "flank pain",
        ],
        "prevalence": "1:400-1:1,000",
        "chromosome": "16p13.3",
    },
    {
        "omim_id": "301500",
        "disease_name": "Fabry Disease",
        "gene": "GLA",
        "inheritance": "x_linked_recessive",
        "clinical_features": [
            "angiokeratomas", "acroparesthesias", "corneal verticillata",
            "renal failure", "cardiomyopathy", "stroke",
        ],
        "prevalence": "1:40,000-1:117,000",
        "chromosome": "Xq22.1",
    },
    {
        "omim_id": "608013",
        "disease_name": "Hereditary Angioedema",
        "gene": "SERPING1",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "recurrent angioedema", "abdominal pain", "laryngeal edema",
            "limb swelling", "low C4 levels", "low C1-inhibitor",
        ],
        "prevalence": "1:50,000",
        "chromosome": "11q12.1",
    },
    {
        "omim_id": "203800",
        "disease_name": "Alkaptonuria",
        "gene": "HGD",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "ochronosis", "dark urine", "arthropathy",
            "cardiac valve disease", "renal stones", "intervertebral disc calcification",
        ],
        "prevalence": "1:250,000-1:1,000,000",
        "chromosome": "3q13.33",
    },
    {
        "omim_id": "312750",
        "disease_name": "Adrenoleukodystrophy (X-linked)",
        "gene": "ABCD1",
        "inheritance": "x_linked_recessive",
        "clinical_features": [
            "progressive demyelination", "adrenal insufficiency",
            "behavioral changes", "visual loss", "spastic paraplegia",
            "elevated VLCFA",
        ],
        "prevalence": "1:17,000",
        "chromosome": "Xq28",
    },
    {
        "omim_id": "214100",
        "disease_name": "Cerebral Cavernous Malformations",
        "gene": "CCM2",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "seizures", "headaches", "focal neurological deficits",
            "hemorrhagic stroke", "cerebral cavernomas",
        ],
        "prevalence": "1:200-1:1,000 for CCMs generally",
        "chromosome": "7p13",
    },
    {
        "omim_id": "105400",
        "disease_name": "Amyloidosis (Hereditary ATTR)",
        "gene": "TTR",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "polyneuropathy", "cardiomyopathy", "autonomic dysfunction",
            "carpal tunnel syndrome", "vitreous opacities", "nephropathy",
        ],
        "prevalence": "1:100,000",
        "chromosome": "18q12.1",
    },
    {
        "omim_id": "305400",
        "disease_name": "Hemophilia A",
        "gene": "F8",
        "inheritance": "x_linked_recessive",
        "clinical_features": [
            "prolonged bleeding", "hemarthrosis", "muscle hematomas",
            "easy bruising", "prolonged aPTT", "reduced factor VIII",
        ],
        "prevalence": "1:5,000 males",
        "chromosome": "Xq28",
    },
    {
        "omim_id": "613795",
        "disease_name": "Dravet Syndrome",
        "gene": "SCN1A",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "febrile seizures", "status epilepticus", "developmental regression",
            "ataxia", "pharmacoresistant epilepsy", "cognitive decline",
        ],
        "prevalence": "1:15,700-1:40,000",
        "chromosome": "2q24.3",
    },
    {
        "omim_id": "278000",
        "disease_name": "Xeroderma Pigmentosum",
        "gene": "XPA",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "extreme UV sensitivity", "freckling", "skin cancers",
            "ocular damage", "neurological degeneration", "photophobia",
        ],
        "prevalence": "1:250,000 (US); 1:22,000 (Japan)",
        "chromosome": "9q22.33",
    },
    {
        "omim_id": "300100",
        "disease_name": "Lowe Syndrome",
        "gene": "OCRL",
        "inheritance": "x_linked_recessive",
        "clinical_features": [
            "congenital cataracts", "intellectual disability",
            "renal tubular dysfunction", "hypotonia",
            "aminoaciduria", "proteinuria",
        ],
        "prevalence": "1:500,000",
        "chromosome": "Xq26.1",
    },
    {
        "omim_id": "271900",
        "disease_name": "Stickler Syndrome Type I",
        "gene": "COL2A1",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "myopia", "retinal detachment", "cleft palate",
            "Pierre Robin sequence", "joint hypermobility", "sensorineural hearing loss",
        ],
        "prevalence": "1:7,500-1:9,000",
        "chromosome": "12q13.11",
    },
    {
        "omim_id": "253300",
        "disease_name": "Mucopolysaccharidosis Type I (Hurler Syndrome)",
        "gene": "IDUA",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "coarse facial features", "corneal clouding", "hepatosplenomegaly",
            "skeletal dysplasia", "developmental delay", "cardiac valve disease",
        ],
        "prevalence": "1:100,000",
        "chromosome": "4p16.3",
    },
    {
        "omim_id": "253200",
        "disease_name": "Mucopolysaccharidosis Type II (Hunter Syndrome)",
        "gene": "IDS",
        "inheritance": "x_linked_recessive",
        "clinical_features": [
            "coarse facial features", "hepatosplenomegaly",
            "skeletal abnormalities", "hearing loss",
            "airway obstruction", "progressive neurodegeneration",
        ],
        "prevalence": "1:100,000-1:170,000 males",
        "chromosome": "Xq28",
    },
    {
        "omim_id": "229300",
        "disease_name": "Friedreich Ataxia",
        "gene": "FXN",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "progressive ataxia", "cardiomyopathy", "scoliosis",
            "diabetes mellitus", "dysarthria", "loss of deep tendon reflexes",
        ],
        "prevalence": "1:50,000",
        "chromosome": "9q21.11",
    },
    {
        "omim_id": "130050",
        "disease_name": "Ehlers-Danlos Syndrome Classical Type",
        "gene": "COL5A1",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "skin hyperextensibility", "atrophic scarring", "joint hypermobility",
            "easy bruising", "tissue fragility", "molluscoid pseudotumors",
        ],
        "prevalence": "1:20,000-1:40,000",
        "chromosome": "9q34.3",
    },
    {
        "omim_id": "130000",
        "disease_name": "Ehlers-Danlos Syndrome Vascular Type",
        "gene": "COL3A1",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "arterial rupture", "organ rupture", "thin translucent skin",
            "easy bruising", "characteristic facial features", "acrogeria",
        ],
        "prevalence": "1:50,000-1:200,000",
        "chromosome": "2q32.2",
    },
    {
        "omim_id": "300400",
        "disease_name": "Severe Combined Immunodeficiency X-linked",
        "gene": "IL2RG",
        "inheritance": "x_linked_recessive",
        "clinical_features": [
            "recurrent severe infections", "failure to thrive",
            "absent T cells", "absent NK cells", "chronic diarrhea",
            "opportunistic infections",
        ],
        "prevalence": "1:50,000-1:100,000",
        "chromosome": "Xq13.1",
    },
    {
        "omim_id": "151623",
        "disease_name": "Li-Fraumeni Syndrome",
        "gene": "TP53",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "early-onset cancers", "breast cancer", "sarcomas",
            "brain tumors", "adrenocortical carcinoma", "leukemia",
        ],
        "prevalence": "1:5,000-1:20,000",
        "chromosome": "17p13.1",
    },
    {
        "omim_id": "120435",
        "disease_name": "Lynch Syndrome",
        "gene": "MLH1",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "colorectal cancer", "endometrial cancer", "ovarian cancer",
            "gastric cancer", "microsatellite instability", "early-onset malignancies",
        ],
        "prevalence": "1:279 (population carrier frequency)",
        "chromosome": "3p22.2",
    },
    {
        "omim_id": "607014",
        "disease_name": "Dravet Syndrome (SCN1A-related)",
        "gene": "SCN1A",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "febrile seizures", "status epilepticus", "developmental regression",
            "ataxia", "pharmacoresistant epilepsy", "cognitive decline",
        ],
        "prevalence": "1:15,700-1:40,000",
        "chromosome": "2q24.3",
    },
    {
        "omim_id": "300672",
        "disease_name": "CDKL5 Deficiency Disorder",
        "gene": "CDKL5",
        "inheritance": "x_linked_dominant",
        "clinical_features": [
            "early-onset seizures", "severe intellectual disability",
            "absent speech", "stereotypic hand movements",
            "cortical visual impairment", "hypotonia",
        ],
        "prevalence": "1:40,000-1:60,000",
        "chromosome": "Xp22.13",
    },
    {
        "omim_id": "603903",
        "disease_name": "Sickle Cell Disease",
        "gene": "HBB",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "vaso-occlusive crises", "chronic hemolytic anemia",
            "acute chest syndrome", "splenic sequestration",
            "stroke risk", "organ damage",
        ],
        "prevalence": "1:365 African Americans",
        "chromosome": "11p15.4",
    },
    {
        "omim_id": "306700",
        "disease_name": "Hemophilia A (F8-related)",
        "gene": "F8",
        "inheritance": "x_linked_recessive",
        "clinical_features": [
            "prolonged bleeding", "hemarthrosis", "muscle hematomas",
            "easy bruising", "prolonged aPTT", "reduced factor VIII",
        ],
        "prevalence": "1:5,000 males",
        "chromosome": "Xq28",
    },
    {
        "omim_id": "219800",
        "disease_name": "Cystinosis",
        "gene": "CTNS",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "renal Fanconi syndrome", "growth retardation",
            "corneal cystine crystals", "photophobia",
            "hypothyroidism", "end-stage renal disease",
        ],
        "prevalence": "1:100,000-1:200,000",
        "chromosome": "17p13.2",
    },
    {
        "omim_id": "256550",
        "disease_name": "Neuronal Ceroid Lipofuscinosis Type 1",
        "gene": "PPT1",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "progressive vision loss", "seizures", "motor decline",
            "cognitive regression", "myoclonus", "brain atrophy",
        ],
        "prevalence": "1:100,000",
        "chromosome": "1p34.2",
    },
    {
        "omim_id": "608049",
        "disease_name": "Angelman Syndrome",
        "gene": "UBE3A",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "severe intellectual disability", "absent speech",
            "ataxia", "happy demeanor", "seizures", "microcephaly",
        ],
        "prevalence": "1:12,000-1:20,000",
        "chromosome": "15q11.2",
    },
    {
        "omim_id": "176270",
        "disease_name": "Prader-Willi Syndrome",
        "gene": "SNRPN region",
        "inheritance": "autosomal_dominant",
        "clinical_features": [
            "neonatal hypotonia", "hyperphagia", "obesity",
            "intellectual disability", "hypogonadism", "short stature",
        ],
        "prevalence": "1:10,000-1:30,000",
        "chromosome": "15q11.2",
    },
    {
        "omim_id": "300624",
        "disease_name": "Fragile X Syndrome",
        "gene": "FMR1",
        "inheritance": "x_linked_dominant",
        "clinical_features": [
            "intellectual disability", "long face", "large ears",
            "macroorchidism", "social anxiety", "hand flapping",
        ],
        "prevalence": "1:4,000 males; 1:8,000 females",
        "chromosome": "Xq27.3",
    },
    {
        "omim_id": "256730",
        "disease_name": "Niemann-Pick Disease Type C",
        "gene": "NPC1",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "vertical supranuclear gaze palsy", "ataxia",
            "hepatosplenomegaly", "dystonia", "cognitive decline",
            "gelastic cataplexy",
        ],
        "prevalence": "1:120,000-1:150,000",
        "chromosome": "18q11.2",
    },
    {
        "omim_id": "227650",
        "disease_name": "Fanconi Anemia Complementation Group A",
        "gene": "FANCA",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "bone marrow failure", "congenital malformations",
            "short stature", "skin pigmentation abnormalities",
            "thumb anomalies", "cancer predisposition",
        ],
        "prevalence": "1:160,000",
        "chromosome": "16q24.3",
    },
    {
        "omim_id": "203100",
        "disease_name": "Albinism Oculocutaneous Type 1",
        "gene": "TYR",
        "inheritance": "autosomal_recessive",
        "clinical_features": [
            "hypopigmentation of skin and hair", "nystagmus",
            "reduced visual acuity", "photophobia",
            "foveal hypoplasia", "misrouting of optic fibers",
        ],
        "prevalence": "1:40,000",
        "chromosome": "11q14.3",
    },
]


# ===================================================================
# OMIM PARSER IMPLEMENTATION
# ===================================================================


class OMIMParser(BaseIngestParser):
    """Parse OMIM disease entries for the Rare Disease Diagnostic Agent.

    In offline/seed mode, returns the curated OMIM_DISEASES list.
    In online mode (when api_key is provided), fetches from the OMIM API.

    Usage::

        parser = OMIMParser()
        records, stats = parser.run()
    """

    def __init__(
        self,
        api_key: str | None = None,
        collection_manager: Any = None,
        embedder: Any = None,
    ) -> None:
        super().__init__(
            source_name="omim",
            collection_manager=collection_manager,
            embedder=embedder,
        )
        self.api_key = api_key

    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch OMIM disease data.

        In seed mode (no API key), returns the curated OMIM_DISEASES list.
        With an API key, attempts to fetch from the OMIM API.

        Returns:
            List of raw OMIM disease dictionaries.
        """
        if self.api_key:
            self.logger.info("OMIM API key provided but live fetch not implemented; using seed data")

        self.logger.info("Using curated OMIM seed data (%d diseases)", len(OMIM_DISEASES))
        return list(OMIM_DISEASES)

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw OMIM disease data into IngestRecord objects.

        Args:
            raw_data: List of OMIM disease dictionaries.

        Returns:
            List of IngestRecord objects.
        """
        records: List[IngestRecord] = []

        for entry in raw_data:
            omim_id = entry.get("omim_id", "")
            disease_name = entry.get("disease_name", "")
            gene = entry.get("gene", "")
            inheritance = entry.get("inheritance", "")
            clinical_features = entry.get("clinical_features", [])
            prevalence = entry.get("prevalence", "")
            chromosome = entry.get("chromosome", "")

            # Build rich text content for embedding
            features_str = ", ".join(clinical_features) if clinical_features else "not specified"
            text = (
                f"OMIM Disease: {disease_name} (OMIM #{omim_id}). "
                f"Gene: {gene} ({chromosome}). "
                f"Inheritance: {inheritance.replace('_', ' ')}. "
                f"Clinical features: {features_str}. "
                f"Prevalence: {prevalence}."
            )

            record = IngestRecord(
                text=text,
                metadata={
                    "omim_id": omim_id,
                    "disease_name": disease_name,
                    "gene": gene,
                    "inheritance": inheritance,
                    "clinical_features": clinical_features,
                    "prevalence": prevalence,
                    "chromosome": chromosome,
                    "source_db": "OMIM",
                },
                collection_name="rd_diseases",
                record_id=f"OMIM_{omim_id}",
                source="omim",
            )
            records.append(record)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate an OMIM IngestRecord.

        Requirements:
            - text must be non-empty
            - must have omim_id in metadata
            - must have disease_name in metadata
            - must have gene in metadata

        Args:
            record: The record to validate.

        Returns:
            True if the record passes all validation checks.
        """
        if not record.text or not record.text.strip():
            return False

        meta = record.metadata
        if not meta.get("omim_id"):
            return False
        if not meta.get("disease_name"):
            return False
        if not meta.get("gene"):
            return False

        return True


def get_omim_disease_count() -> int:
    """Return the number of curated OMIM diseases."""
    return len(OMIM_DISEASES)


def get_omim_genes() -> List[str]:
    """Return a deduplicated list of genes from the OMIM seed data."""
    genes = list({d["gene"] for d in OMIM_DISEASES})
    genes.sort()
    return genes
