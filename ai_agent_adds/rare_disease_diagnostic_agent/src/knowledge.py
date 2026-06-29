"""Rare Disease Diagnostic Agent — Domain Knowledge Base.

Comprehensive rare-disease knowledge covering 13 disease categories,
metabolic/neurological/connective-tissue/hematologic/immunologic/cancer-
predisposition/endocrine/skeletal/renal/pulmonary/dermatologic/ophthalmologic
diseases, approved gene therapies, diagnostic algorithms, ACMG variant
classification criteria, and HPO top-level terms.

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict


# ═══════════════════════════════════════════════════════════════════════════════
# KNOWLEDGE BASE VERSION
# ═══════════════════════════════════════════════════════════════════════════════

KNOWLEDGE_VERSION: Dict[str, Any] = {
    "version": "1.0.0",
    "last_updated": "2026-03-22",
    "revision_notes": "Expanded release — 13 disease categories, 28 metabolic diseases, "
                      "23 neurological diseases, 10 connective-tissue diseases, "
                      "15 hematologic diseases, 13 immunologic diseases, 8 cancer "
                      "predisposition syndromes, 12 approved/recent gene therapies, "
                      "9 diagnostic algorithms, 28 ACMG criteria, 23 HPO top-level terms.",
    "sources": [
        "OMIM (Online Mendelian Inheritance in Man)",
        "Orphanet Rare Disease Database",
        "GeneReviews (NCBI)",
        "ClinGen / ClinVar",
        "Human Phenotype Ontology (HPO)",
        "ACMG/AMP Standards and Guidelines (Richards et al. 2015)",
        "NIH Genetic and Rare Diseases Information Center (GARD)",
        "European Reference Networks (ERNs)",
        "FDA Approved Cellular and Gene Therapy Products",
        "Newborn Screening ACTion (ACT) Sheets — ACMG",
    ],
    "counts": {
        "disease_categories": 13,
        "metabolic_diseases": 28,
        "neurological_diseases": 23,
        "connective_tissue_diseases": 10,
        "hematologic_diseases": 15,
        "immunologic_diseases": 13,
        "cancer_predisposition": 8,
        "gene_therapies_approved": 12,
        "diagnostic_algorithms": 9,
        "acmg_criteria": 28,
        "hpo_top_level_terms": 23,
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. RARE DISEASE CATEGORIES
# ═══════════════════════════════════════════════════════════════════════════════

RARE_DISEASE_CATEGORIES: Dict[str, Dict[str, Any]] = {
    "metabolic": {
        "description": "Inborn errors of metabolism affecting enzyme or transporter "
                       "function in intermediary metabolism, lysosomal storage, "
                       "peroxisomal, or mitochondrial pathways.",
        "example_diseases": [
            "Phenylketonuria (PKU)",
            "Gaucher disease",
            "Fabry disease",
            "Pompe disease",
            "Mucopolysaccharidosis type I (Hurler syndrome)",
            "Maple syrup urine disease (MSUD)",
            "Medium-chain acyl-CoA dehydrogenase deficiency (MCADD)",
            "Galactosemia",
        ],
        "key_genes": [
            "PAH", "GBA", "GLA", "GAA", "IDUA", "BCKDHA", "BCKDHB",
            "ACADM", "GALT", "HEXA", "ASM",
        ],
        "diagnostic_approach": "Newborn screening, tandem mass spectrometry, enzyme "
                               "assays, urine organic acids, plasma amino acids, "
                               "acylcarnitine profile, molecular confirmation.",
        "prevalence_range": "1:10,000 – 1:200,000 individually; ~1:800 collectively",
    },
    "neurological": {
        "description": "Genetic disorders primarily affecting the central or peripheral "
                       "nervous system, including neuromuscular, neurodegenerative, and "
                       "neurodevelopmental conditions.",
        "example_diseases": [
            "Spinal muscular atrophy (SMA)",
            "Duchenne muscular dystrophy (DMD)",
            "Rett syndrome",
            "Huntington disease",
            "Friedreich ataxia",
            "Charcot-Marie-Tooth disease",
            "Tuberous sclerosis complex",
            "Dravet syndrome",
        ],
        "key_genes": [
            "SMN1", "SMN2", "DMD", "MECP2", "HTT", "FXN", "PMP22",
            "TSC1", "TSC2", "SCN1A", "NF1",
        ],
        "diagnostic_approach": "Clinical examination, EMG/NCV, brain MRI, genetic "
                               "panel testing, whole-exome/genome sequencing, repeat "
                               "expansion analysis, muscle biopsy (select cases).",
        "prevalence_range": "1:6,000 – 1:100,000",
    },
    "hematologic": {
        "description": "Inherited disorders of hemoglobin, coagulation factors, "
                       "platelets, or bone marrow failure syndromes.",
        "example_diseases": [
            "Sickle cell disease",
            "Beta-thalassemia",
            "Hemophilia A",
            "Hemophilia B",
            "Von Willebrand disease",
            "Diamond-Blackfan anemia",
            "Fanconi anemia",
            "Hereditary spherocytosis",
        ],
        "key_genes": [
            "HBB", "HBA1", "HBA2", "F8", "F9", "VWF", "RPS19",
            "FANCA", "FANCC", "ANK1",
        ],
        "diagnostic_approach": "Complete blood count, hemoglobin electrophoresis, "
                               "coagulation studies, bone marrow biopsy, chromosomal "
                               "breakage testing (Fanconi), targeted gene panels.",
        "prevalence_range": "1:500 (SCD in African Americans) – 1:100,000",
    },
    "connective_tissue": {
        "description": "Heritable disorders of connective tissue affecting collagen, "
                       "fibrillin, or extracellular matrix components, impacting "
                       "skin, joints, vessels, and skeleton.",
        "example_diseases": [
            "Marfan syndrome",
            "Ehlers-Danlos syndrome (classical)",
            "Ehlers-Danlos syndrome (vascular)",
            "Osteogenesis imperfecta",
            "Loeys-Dietz syndrome",
            "Stickler syndrome",
        ],
        "key_genes": [
            "FBN1", "COL5A1", "COL5A2", "COL3A1", "COL1A1", "COL1A2",
            "TGFBR1", "TGFBR2", "COL2A1", "COL11A1",
        ],
        "diagnostic_approach": "Clinical criteria (Ghent for Marfan, Villefranche/2017 "
                               "for EDS), echocardiography, skin biopsy with electron "
                               "microscopy, collagen biochemistry, molecular testing.",
        "prevalence_range": "1:5,000 – 1:250,000",
    },
    "immunologic": {
        "description": "Primary immunodeficiency disorders (PIDs) affecting innate "
                       "or adaptive immunity, leading to recurrent infections, "
                       "autoimmunity, or lymphoproliferation.",
        "example_diseases": [
            "Severe combined immunodeficiency (SCID)",
            "Chronic granulomatous disease (CGD)",
            "Hyper-IgE syndrome",
            "Common variable immunodeficiency (CVID)",
            "Wiskott-Aldrich syndrome",
            "X-linked agammaglobulinemia (XLA)",
        ],
        "key_genes": [
            "IL2RG", "JAK3", "RAG1", "RAG2", "CYBB", "STAT3",
            "TNFRSF13B", "WAS", "BTK",
        ],
        "diagnostic_approach": "Immunoglobulin levels, lymphocyte subsets, TREC assay "
                               "(newborn screening for SCID), neutrophil oxidative "
                               "burst (DHR for CGD), targeted gene panels, flow cytometry.",
        "prevalence_range": "1:1,200 (CVID) – 1:100,000 (SCID)",
    },
    "cardiac": {
        "description": "Inherited cardiac conditions including cardiomyopathies, "
                       "channelopathies, and structural heart defects with "
                       "Mendelian inheritance.",
        "example_diseases": [
            "Hypertrophic cardiomyopathy (HCM)",
            "Dilated cardiomyopathy (DCM)",
            "Long QT syndrome",
            "Brugada syndrome",
            "Arrhythmogenic right ventricular cardiomyopathy (ARVC)",
            "Familial hypercholesterolemia",
        ],
        "key_genes": [
            "MYH7", "MYBPC3", "TTN", "LMNA", "KCNQ1", "KCNH2",
            "SCN5A", "PKP2", "DSP", "LDLR", "PCSK9",
        ],
        "diagnostic_approach": "ECG, echocardiography, cardiac MRI, exercise stress "
                               "testing, Holter monitoring, cascade family screening, "
                               "targeted gene panels.",
        "prevalence_range": "1:200 (FH) – 1:5,000 (HCM)",
    },
    "cancer_predisposition": {
        "description": "Hereditary cancer predisposition syndromes caused by germline "
                       "pathogenic variants in tumor suppressor genes or DNA repair "
                       "genes, conferring elevated lifetime cancer risk.",
        "example_diseases": [
            "Li-Fraumeni syndrome",
            "Lynch syndrome (HNPCC)",
            "Hereditary breast/ovarian cancer (BRCA1/2)",
            "Familial adenomatous polyposis (FAP)",
            "Multiple endocrine neoplasia type 1 (MEN1)",
            "Multiple endocrine neoplasia type 2 (MEN2)",
            "Von Hippel-Lindau syndrome",
            "Retinoblastoma",
        ],
        "key_genes": [
            "TP53", "MLH1", "MSH2", "MSH6", "PMS2", "BRCA1", "BRCA2",
            "APC", "MEN1", "RET", "VHL", "RB1",
        ],
        "diagnostic_approach": "Family history (3+ generation pedigree), NCCN criteria, "
                               "germline multi-gene panel testing, tumor profiling "
                               "(MSI/IHC for Lynch), predictive testing for at-risk relatives.",
        "prevalence_range": "1:200 (BRCA carriers) – 1:20,000 (Li-Fraumeni)",
    },
    "endocrine": {
        "description": "Rare genetic disorders affecting endocrine glands, hormone "
                       "biosynthesis, or hormone receptor signaling.",
        "example_diseases": [
            "Congenital adrenal hyperplasia (CAH)",
            "Turner syndrome",
            "Klinefelter syndrome",
            "Congenital hypothyroidism",
            "Noonan syndrome",
            "Kallmann syndrome",
        ],
        "key_genes": [
            "CYP21A2", "CYP11B1", "PTPN11", "SOS1", "RAF1",
            "KAL1", "FGFR1", "PROKR2", "TSHR",
        ],
        "diagnostic_approach": "Hormone panels, newborn screening (17-OHP for CAH, TSH), "
                               "karyotype, growth curves, targeted gene panels, "
                               "dynamic endocrine testing.",
        "prevalence_range": "1:2,500 (Turner) – 1:15,000 (CAH)",
    },
    "skeletal": {
        "description": "Genetic skeletal dysplasias and bone disorders affecting "
                       "growth plate, bone density, or skeletal patterning.",
        "example_diseases": [
            "Achondroplasia",
            "Hypophosphatasia",
            "Osteogenesis imperfecta",
            "Thanatophoric dysplasia",
            "Spondyloepiphyseal dysplasia",
            "Cleidocranial dysplasia",
        ],
        "key_genes": [
            "FGFR3", "ALPL", "COL1A1", "COL1A2", "COL2A1",
            "RUNX2", "SHOX", "SOX9",
        ],
        "diagnostic_approach": "Skeletal survey, growth charts, bone density (DXA), "
                               "alkaline phosphatase, calcium/phosphate levels, "
                               "skeletal dysplasia gene panels, prenatal ultrasound.",
        "prevalence_range": "1:15,000 (achondroplasia) – 1:100,000+",
    },
    "renal": {
        "description": "Inherited kidney diseases affecting glomerular, tubular, "
                       "or structural development of the kidneys.",
        "example_diseases": [
            "Autosomal dominant polycystic kidney disease (ADPKD)",
            "Autosomal recessive polycystic kidney disease (ARPKD)",
            "Alport syndrome",
            "Cystinosis",
            "Fabry disease (renal manifestation)",
            "Bartter syndrome",
        ],
        "key_genes": [
            "PKD1", "PKD2", "PKHD1", "COL4A3", "COL4A4", "COL4A5",
            "CTNS", "GLA", "SLC12A1", "CLCNKB",
        ],
        "diagnostic_approach": "Renal ultrasound, urinalysis, serum creatinine/GFR, "
                               "kidney biopsy with electron microscopy, slit lamp exam "
                               "(Alport), white cell cystine levels, targeted gene panels.",
        "prevalence_range": "1:400 (ADPKD) – 1:100,000+",
    },
    "pulmonary": {
        "description": "Rare genetic lung disorders affecting surfactant metabolism, "
                       "airway development, or pulmonary vasculature.",
        "example_diseases": [
            "Cystic fibrosis",
            "Primary ciliary dyskinesia",
            "Alpha-1 antitrypsin deficiency",
            "Pulmonary alveolar proteinosis",
            "Hereditary hemorrhagic telangiectasia",
            "Hermansky-Pudlak syndrome (pulmonary fibrosis)",
        ],
        "key_genes": [
            "CFTR", "DNAI1", "DNAH5", "SERPINA1", "CSF2RA",
            "ENG", "ACVRL1", "HPS1", "SFTPC",
        ],
        "diagnostic_approach": "Sweat chloride test, nasal NO, ciliary electron "
                               "microscopy, alpha-1 level, HRCT, pulmonary function "
                               "tests, genetic testing.",
        "prevalence_range": "1:2,500 (CF in Caucasians) – 1:100,000+",
    },
    "dermatologic": {
        "description": "Inherited skin disorders (genodermatoses) affecting keratin, "
                       "pigmentation, adhesion, or structural skin proteins.",
        "example_diseases": [
            "Epidermolysis bullosa (EB)",
            "Ichthyosis (various types)",
            "Xeroderma pigmentosum",
            "Tuberous sclerosis (skin manifestations)",
            "Neurofibromatosis type 1 (skin manifestations)",
            "Albinism",
        ],
        "key_genes": [
            "KRT5", "KRT14", "COL7A1", "LAMB3", "TGM1",
            "XPA", "XPC", "TYR", "OCA2", "NF1",
        ],
        "diagnostic_approach": "Skin biopsy with immunofluorescence mapping, electron "
                               "microscopy, Wood's lamp, dermoscopy, targeted gene panels, "
                               "prenatal testing in severe forms.",
        "prevalence_range": "1:17,000 (NF1) – 1:1,000,000 (XP)",
    },
    "ophthalmologic": {
        "description": "Inherited eye disorders affecting retinal function, lens, "
                       "optic nerve, or anterior segment structures.",
        "example_diseases": [
            "Retinitis pigmentosa",
            "Leber congenital amaurosis",
            "Stargardt disease",
            "Aniridia",
            "Retinoblastoma",
            "Leber hereditary optic neuropathy (LHON)",
        ],
        "key_genes": [
            "RHO", "RPE65", "ABCA4", "PAX6", "RB1",
            "MT-ND4", "RPGR", "RS1", "CRB1", "CEP290",
        ],
        "diagnostic_approach": "Fundoscopy, OCT, electroretinography (ERG), visual "
                               "field testing, dark adaptometry, inherited retinal "
                               "disease gene panels, RPE65 testing for gene therapy.",
        "prevalence_range": "1:3,000 (retinitis pigmentosa) – 1:40,000 (LCA)",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. METABOLIC DISEASES (detailed)
# ═══════════════════════════════════════════════════════════════════════════════

METABOLIC_DISEASES: Dict[str, Dict[str, Any]] = {
    "pku": {
        "name": "Phenylketonuria",
        "gene": "PAH",
        "enzyme": "Phenylalanine hydroxylase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Low-phenylalanine diet, sapropterin (BH4-responsive), "
                     "pegvaliase (enzyme substitution)",
        "omim_id": "261600",
    },
    "gaucher": {
        "name": "Gaucher disease",
        "gene": "GBA1",
        "enzyme": "Glucocerebrosidase (acid beta-glucosidase)",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "Enzyme replacement therapy (imiglucerase, velaglucerase), "
                     "substrate reduction therapy (eliglustat, miglustat)",
        "omim_id": "230800",
    },
    "fabry": {
        "name": "Fabry disease",
        "gene": "GLA",
        "enzyme": "Alpha-galactosidase A",
        "inheritance": "X-linked",
        "newborn_screening": True,
        "treatment": "Enzyme replacement therapy (agalsidase alfa/beta), "
                     "oral chaperone therapy (migalastat)",
        "omim_id": "301500",
    },
    "pompe": {
        "name": "Pompe disease",
        "gene": "GAA",
        "enzyme": "Acid alpha-glucosidase (acid maltase)",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Enzyme replacement therapy (alglucosidase alfa, avalglucosidase alfa)",
        "omim_id": "232300",
    },
    "mps_i": {
        "name": "Mucopolysaccharidosis type I (Hurler/Scheie)",
        "gene": "IDUA",
        "enzyme": "Alpha-L-iduronidase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Enzyme replacement therapy (laronidase), hematopoietic stem cell "
                     "transplant (severe/Hurler form)",
        "omim_id": "607014",
    },
    "mps_ii": {
        "name": "Mucopolysaccharidosis type II (Hunter syndrome)",
        "gene": "IDS",
        "enzyme": "Iduronate-2-sulfatase",
        "inheritance": "X-linked recessive",
        "newborn_screening": False,
        "treatment": "Enzyme replacement therapy (idursulfase), intrathecal ERT "
                     "for CNS (investigational)",
        "omim_id": "309900",
    },
    "mps_iii": {
        "name": "Mucopolysaccharidosis type III (Sanfilippo syndrome)",
        "gene": "SGSH / NAGLU / HGSNAT / GNS",
        "enzyme": "Heparan sulfate degradation enzymes (subtypes A-D)",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "Supportive care; gene therapy trials underway",
        "omim_id": "252900",
    },
    "mps_iv": {
        "name": "Mucopolysaccharidosis type IV (Morquio syndrome)",
        "gene": "GALNS / GLB1",
        "enzyme": "N-acetylgalactosamine-6-sulfatase (type A) / beta-galactosidase (type B)",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "Enzyme replacement therapy (elosulfase alfa for type A)",
        "omim_id": "253000",
    },
    "mps_vi": {
        "name": "Mucopolysaccharidosis type VI (Maroteaux-Lamy syndrome)",
        "gene": "ARSB",
        "enzyme": "N-acetylgalactosamine 4-sulfatase (arylsulfatase B)",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "Enzyme replacement therapy (galsulfase)",
        "omim_id": "253200",
    },
    "mps_vii": {
        "name": "Mucopolysaccharidosis type VII (Sly syndrome)",
        "gene": "GUSB",
        "enzyme": "Beta-glucuronidase",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "Enzyme replacement therapy (vestronidase alfa)",
        "omim_id": "253220",
    },
    "galactosemia": {
        "name": "Classic galactosemia",
        "gene": "GALT",
        "enzyme": "Galactose-1-phosphate uridylyltransferase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Galactose-restricted diet (lifelong)",
        "omim_id": "230400",
    },
    "msud": {
        "name": "Maple syrup urine disease",
        "gene": "BCKDHA / BCKDHB / DBT",
        "enzyme": "Branched-chain alpha-ketoacid dehydrogenase complex",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Branched-chain amino acid-restricted diet, thiamine trial "
                     "(thiamine-responsive variant), liver transplant in severe cases",
        "omim_id": "248600",
    },
    "mcad": {
        "name": "Medium-chain acyl-CoA dehydrogenase deficiency",
        "gene": "ACADM",
        "enzyme": "Medium-chain acyl-CoA dehydrogenase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Avoidance of fasting, emergency protocol with IV dextrose, "
                     "medium-chain triglyceride avoidance",
        "omim_id": "201450",
    },
    "otc_deficiency": {
        "name": "Ornithine transcarbamylase deficiency",
        "gene": "OTC",
        "enzyme": "Ornithine transcarbamylase",
        "inheritance": "X-linked",
        "newborn_screening": False,
        "treatment": "Low-protein diet, nitrogen scavengers (sodium benzoate, "
                     "sodium phenylbutyrate, glycerol phenylbutyrate), "
                     "liver transplant in severe neonatal onset",
        "omim_id": "311250",
    },
    "propionic_acidemia": {
        "name": "Propionic acidemia",
        "gene": "PCCA / PCCB",
        "enzyme": "Propionyl-CoA carboxylase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Protein-restricted diet, carnitine supplementation, "
                     "metronidazole (reduces propionic acid-producing gut bacteria), "
                     "liver transplant",
        "omim_id": "606054",
    },
    "methylmalonic_acidemia": {
        "name": "Methylmalonic acidemia",
        "gene": "MUT / MMAA / MMAB",
        "enzyme": "Methylmalonyl-CoA mutase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Protein-restricted diet, hydroxocobalamin (B12-responsive forms), "
                     "carnitine, liver or combined liver-kidney transplant",
        "omim_id": "251000",
    },
    "niemann_pick_a_b": {
        "name": "Niemann-Pick disease type A/B",
        "gene": "SMPD1",
        "enzyme": "Acid sphingomyelinase",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "Enzyme replacement therapy (olipudase alfa for type B), "
                     "supportive care (type A is severe infantile, usually fatal)",
        "omim_id": "257200",
    },
    "niemann_pick_c": {
        "name": "Niemann-Pick disease type C",
        "gene": "NPC1 / NPC2",
        "enzyme": "NPC1/NPC2 cholesterol transport proteins",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "Miglustat (substrate reduction), arimoclomol (investigational)",
        "omim_id": "257220",
    },
    "tay_sachs": {
        "name": "Tay-Sachs disease",
        "gene": "HEXA",
        "enzyme": "Hexosaminidase A",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "Supportive care; no approved disease-modifying therapy. "
                     "Gene therapy trials underway.",
        "omim_id": "272800",
    },
    "krabbe": {
        "name": "Krabbe disease (globoid cell leukodystrophy)",
        "gene": "GALC",
        "enzyme": "Galactosylceramidase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Hematopoietic stem cell transplant (pre-symptomatic only), "
                     "supportive care",
        "omim_id": "245200",
    },
    "homocystinuria": {
        "name": "Homocystinuria (classical)",
        "gene": "CBS",
        "enzyme": "Cystathionine beta-synthase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Pyridoxine (B6-responsive forms), betaine, methionine-restricted diet, "
                     "folate and B12 supplementation",
        "omim_id": "236200",
    },
    "biotinidase_deficiency": {
        "name": "Biotinidase deficiency",
        "gene": "BTD",
        "enzyme": "Biotinidase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Biotin supplementation (lifelong, 5-20 mg/day oral biotin)",
        "omim_id": "253260",
    },
    "glutaric_aciduria_type1": {
        "name": "Glutaric aciduria type I",
        "gene": "GCDH",
        "enzyme": "Glutaryl-CoA dehydrogenase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Carnitine supplementation, lysine-restricted diet, emergency protocol "
                     "during illness to prevent encephalopathic crises",
        "omim_id": "231670",
    },
    "tyrosinemia_type1": {
        "name": "Tyrosinemia type I (hepatorenal tyrosinemia)",
        "gene": "FAH",
        "enzyme": "Fumarylacetoacetate hydrolase",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Nitisinone (NTBC) — first-line, phenylalanine- and tyrosine-restricted diet, "
                     "liver transplant if refractory",
        "omim_id": "276700",
    },
    "glycogen_storage_disease_1": {
        "name": "Glycogen storage disease type Ia (von Gierke disease)",
        "gene": "G6PC",
        "enzyme": "Glucose-6-phosphatase",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "Frequent feeds with uncooked cornstarch therapy, hepatomegaly monitoring, "
                     "avoidance of fasting, renal and hepatic surveillance",
        "omim_id": "232200",
    },
    "cystinosis": {
        "name": "Cystinosis (nephropathic)",
        "gene": "CTNS",
        "enzyme": "Cystinosin (lysosomal cystine transporter)",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "Cysteamine therapy (oral and ophthalmic), renal Fanconi syndrome management, "
                     "electrolyte replacement, kidney transplant",
        "omim_id": "219800",
    },
    "phenylketonuria_bh4": {
        "name": "Phenylketonuria (BH4-responsive)",
        "gene": "PAH",
        "enzyme": "Phenylalanine hydroxylase (BH4-responsive variants)",
        "inheritance": "autosomal recessive",
        "newborn_screening": True,
        "treatment": "Sapropterin (Kuvan / BH4) — reduces phenylalanine levels in responsive patients, "
                     "dietary phenylalanine restriction may be relaxed",
        "omim_id": "261600",
    },
    "congenital_disorder_glycosylation": {
        "name": "Congenital disorder of glycosylation type Ia (PMM2-CDG)",
        "gene": "PMM2",
        "enzyme": "Phosphomannomutase 2",
        "inheritance": "autosomal recessive",
        "newborn_screening": False,
        "treatment": "No specific treatment; supportive/multisystem management, physical therapy, "
                     "nutritional support, seizure management, mannose supplementation (investigational)",
        "omim_id": "212065",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3. NEUROLOGICAL DISEASES (detailed)
# ═══════════════════════════════════════════════════════════════════════════════

NEUROLOGICAL_DISEASES: Dict[str, Dict[str, Any]] = {
    "sma": {
        "name": "Spinal muscular atrophy",
        "gene": "SMN1 (survival motor neuron 1)",
        "inheritance": "autosomal recessive",
        "age_onset": "Type I: birth–6 months; Type II: 6–18 months; "
                     "Type III: >18 months; Type IV: adult",
        "key_features": [
            "Progressive proximal muscle weakness",
            "Hypotonia (floppy infant)",
            "Absent or diminished deep tendon reflexes",
            "Tongue fasciculations",
            "Respiratory insufficiency",
        ],
        "treatment": "Nusinersen (antisense oligonucleotide), onasemnogene abeparvovec "
                     "(gene therapy), risdiplam (oral SMN2 splicing modifier)",
    },
    "dmd": {
        "name": "Duchenne muscular dystrophy",
        "gene": "DMD (dystrophin)",
        "inheritance": "X-linked recessive",
        "age_onset": "2–5 years",
        "key_features": [
            "Progressive proximal muscle weakness (legs before arms)",
            "Gowers sign (using hands to rise from floor)",
            "Calf pseudohypertrophy",
            "Elevated serum creatine kinase (10-100x normal)",
            "Loss of ambulation by age 12 typically",
            "Cardiomyopathy",
        ],
        "treatment": "Corticosteroids (deflazacort, prednisone), exon-skipping ASOs "
                     "(eteplirsen, golodirsen, viltolarsen, casimersen), cardiac management",
    },
    "bmd": {
        "name": "Becker muscular dystrophy",
        "gene": "DMD (dystrophin — in-frame deletions)",
        "inheritance": "X-linked recessive",
        "age_onset": "5–15 years (later than DMD)",
        "key_features": [
            "Progressive proximal muscle weakness (milder than DMD)",
            "Preserved ambulation into adulthood",
            "Calf pseudohypertrophy",
            "Cardiomyopathy (can be primary manifestation)",
            "Elevated CK",
        ],
        "treatment": "Supportive care, cardiac monitoring, physical therapy, "
                     "corticosteroids in some cases",
    },
    "rett": {
        "name": "Rett syndrome",
        "gene": "MECP2",
        "inheritance": "X-linked dominant (almost exclusively females)",
        "age_onset": "6–18 months (after period of normal development)",
        "key_features": [
            "Regression of acquired skills",
            "Loss of purposeful hand movements (hand stereotypies)",
            "Loss of spoken language",
            "Gait abnormalities / ataxia",
            "Seizures (60-80%)",
            "Breathing irregularities",
            "Microcephaly (acquired, deceleration of head growth)",
        ],
        "treatment": "Trofinetide (approved 2023), supportive care, seizure management, "
                     "physical/occupational therapy",
    },
    "angelman": {
        "name": "Angelman syndrome",
        "gene": "UBE3A (maternal allele)",
        "inheritance": "Imprinting disorder (maternal deletion/UPD/imprinting defect/UBE3A mutation)",
        "age_onset": "6–12 months",
        "key_features": [
            "Severe intellectual disability",
            "Absent or minimal speech",
            "Happy demeanor with frequent laughing",
            "Movement/balance disorder (ataxia, tremor)",
            "Seizures (80-90%)",
            "Microcephaly",
        ],
        "treatment": "Seizure management, behavioral therapy, ASO therapy (investigational), "
                     "gene therapy (investigational)",
    },
    "prader_willi": {
        "name": "Prader-Willi syndrome",
        "gene": "15q11.2-q13 (paternal deletion/maternal UPD/imprinting defect)",
        "inheritance": "Imprinting disorder (loss of paternal contribution)",
        "age_onset": "Neonatal (hypotonia/feeding difficulties); hyperphagia onset 2-6 years",
        "key_features": [
            "Neonatal hypotonia and poor feeding",
            "Hyperphagia and obesity (onset ~2-6 years)",
            "Short stature",
            "Hypogonadism",
            "Intellectual disability (mild to moderate)",
            "Behavioral issues (tantrums, OCD-like features)",
        ],
        "treatment": "Growth hormone therapy, dietary management, behavioral interventions, "
                     "sex hormone replacement",
    },
    "huntington": {
        "name": "Huntington disease",
        "gene": "HTT (CAG trinucleotide repeat expansion)",
        "inheritance": "autosomal dominant",
        "age_onset": "30–50 years (juvenile form <20 years with >60 repeats)",
        "key_features": [
            "Chorea (involuntary dance-like movements)",
            "Progressive cognitive decline",
            "Psychiatric symptoms (depression, irritability)",
            "Dystonia and rigidity (later stages)",
            "Relentless progression over 15-20 years",
        ],
        "treatment": "Tetrabenazine/deutetrabenazine (chorea), psychiatric medications, "
                     "supportive care, genetic counseling",
    },
    "friedreich_ataxia": {
        "name": "Friedreich ataxia",
        "gene": "FXN (GAA trinucleotide repeat expansion)",
        "inheritance": "autosomal recessive",
        "age_onset": "5–15 years",
        "key_features": [
            "Progressive gait and limb ataxia",
            "Dysarthria",
            "Loss of deep tendon reflexes",
            "Hypertrophic cardiomyopathy",
            "Diabetes mellitus (10-30%)",
            "Scoliosis",
            "Pes cavus (high-arched feet)",
        ],
        "treatment": "Omaveloxolone (Nrf2 activator, approved 2023), cardiac management, "
                     "physical therapy, orthopedic interventions",
    },
    "cmt": {
        "name": "Charcot-Marie-Tooth disease",
        "gene": "PMP22 (CMT1A), MFN2 (CMT2A), GJB1 (CMTX1), and 100+ others",
        "inheritance": "Autosomal dominant (most common), AR, X-linked",
        "age_onset": "First or second decade",
        "key_features": [
            "Slowly progressive distal muscle weakness and atrophy",
            "Foot drop and high-arched feet (pes cavus)",
            "Distal sensory loss",
            "Reduced nerve conduction velocities (demyelinating forms)",
            "Hammer toes",
        ],
        "treatment": "Supportive care (AFOs, physical therapy, surgery for foot deformities), "
                     "PXT3003 and gene therapy (investigational)",
    },
    "tuberous_sclerosis": {
        "name": "Tuberous sclerosis complex",
        "gene": "TSC1 (hamartin) / TSC2 (tuberin)",
        "inheritance": "autosomal dominant (2/3 de novo)",
        "age_onset": "Variable; often prenatal (cardiac rhabdomyomas) to infancy",
        "key_features": [
            "Cortical tubers and subependymal nodules",
            "Seizures (infantile spasms, focal seizures)",
            "Intellectual disability (variable)",
            "Skin findings (hypomelanotic macules, angiofibromas, shagreen patches)",
            "Renal angiomyolipomas",
            "Pulmonary lymphangioleiomyomatosis (LAM)",
        ],
        "treatment": "mTOR inhibitors (everolimus for SEGA, renal AML, LAM), "
                     "vigabatrin (infantile spasms), seizure management",
    },
    "nf1": {
        "name": "Neurofibromatosis type 1",
        "gene": "NF1 (neurofibromin)",
        "inheritance": "autosomal dominant (50% de novo)",
        "age_onset": "Childhood (cafe-au-lait macules from birth)",
        "key_features": [
            "Cafe-au-lait macules (6+)",
            "Neurofibromas (cutaneous, plexiform)",
            "Freckling in axillary/inguinal regions",
            "Lisch nodules (iris hamartomas)",
            "Optic pathway gliomas",
            "Learning disabilities",
            "Increased cancer risk (MPNST)",
        ],
        "treatment": "Selumetinib (MEK inhibitor for plexiform neurofibromas), "
                     "surgical management, monitoring/surveillance",
    },
    "dravet": {
        "name": "Dravet syndrome (severe myoclonic epilepsy of infancy)",
        "gene": "SCN1A",
        "inheritance": "autosomal dominant (usually de novo)",
        "age_onset": "5–8 months",
        "key_features": [
            "Prolonged febrile seizures in infancy",
            "Multiple seizure types (myoclonic, focal, generalized tonic-clonic)",
            "Status epilepticus episodes",
            "Developmental regression after age 2",
            "Gait abnormalities (crouch gait)",
            "Temperature sensitivity",
        ],
        "treatment": "Stiripentol, cannabidiol (Epidiolex), fenfluramine; avoid sodium "
                     "channel blockers (carbamazepine, phenytoin)",
    },
    "sma_type_0": {
        "name": "Spinal muscular atrophy type 0 (prenatal onset)",
        "gene": "SMN1",
        "inheritance": "autosomal recessive",
        "age_onset": "Prenatal (decreased fetal movement)",
        "key_features": [
            "Severe neonatal hypotonia",
            "Arthrogryposis",
            "Respiratory failure at birth",
            "Joint contractures",
        ],
        "treatment": "Palliative care; gene therapy considered in some cases when "
                     "diagnosed prenatally",
    },
    "ataxia_telangiectasia": {
        "name": "Ataxia-telangiectasia",
        "gene": "ATM",
        "inheritance": "autosomal recessive",
        "age_onset": "1–4 years",
        "key_features": [
            "Progressive cerebellar ataxia",
            "Oculocutaneous telangiectasias",
            "Immunodeficiency (low IgA, IgG subclass)",
            "Elevated alpha-fetoprotein",
            "Cancer predisposition (lymphoma, leukemia)",
            "Radiosensitivity",
        ],
        "treatment": "Supportive care, immunoglobulin replacement, infection prophylaxis, "
                     "cancer surveillance",
    },
    "wilson_disease": {
        "name": "Wilson disease (hepatolenticular degeneration)",
        "gene": "ATP7B",
        "inheritance": "autosomal recessive",
        "age_onset": "5–35 years",
        "key_features": [
            "Hepatic dysfunction (hepatitis, cirrhosis, acute liver failure)",
            "Neuropsychiatric symptoms (tremor, dystonia, dysarthria, personality changes)",
            "Kayser-Fleischer rings (corneal copper deposits)",
            "Low serum ceruloplasmin",
            "Elevated urinary copper",
        ],
        "treatment": "Copper chelation (D-penicillamine, trientine), zinc acetate, "
                     "liver transplant in acute liver failure",
    },
    "cdkl5_deficiency": {
        "name": "CDKL5 deficiency disorder",
        "gene": "CDKL5",
        "inheritance": "X-linked dominant",
        "age_onset": "First months of life",
        "key_features": [
            "Early-onset seizures (often infantile spasms)",
            "Severe developmental delay / intellectual disability",
            "Absent or limited speech",
            "Hand stereotypies (similar to Rett but distinct)",
            "Cortical visual impairment",
            "Hypotonia",
        ],
        "treatment": "Ganaxolone (approved 2022 for seizures in CDKL5), antiseizure medications, "
                     "supportive therapies",
    },
    "scn2a_related": {
        "name": "SCN2A-related disorder",
        "gene": "SCN2A",
        "inheritance": "autosomal dominant (usually de novo)",
        "age_onset": "Neonatal to early infantile",
        "key_features": [
            "Neonatal or infantile seizures",
            "Epileptic encephalopathy (severe forms)",
            "Autism spectrum disorder (some variants)",
            "Intellectual disability (variable)",
            "Movement disorders",
        ],
        "treatment": "Sodium channel blockers (carbamazepine, phenytoin — for gain-of-function variants), "
                     "avoid sodium channel blockers in loss-of-function variants",
    },
    "canavan_disease": {
        "name": "Canavan disease",
        "gene": "ASPA",
        "inheritance": "autosomal recessive",
        "age_onset": "3–6 months",
        "key_features": [
            "Macrocephaly (progressive)",
            "Severe developmental delay and regression",
            "Hypotonia progressing to spasticity",
            "Leukodystrophy on MRI (diffuse white matter involvement)",
            "Elevated N-acetylaspartic acid (NAA) in urine",
            "Optic atrophy",
        ],
        "treatment": "Supportive care; experimental gene therapy (AAV-based ASPA delivery) "
                     "in clinical trials",
    },
    "alexander_disease": {
        "name": "Alexander disease",
        "gene": "GFAP",
        "inheritance": "autosomal dominant (usually de novo)",
        "age_onset": "Infantile (most common), juvenile, or adult forms",
        "key_features": [
            "Leukodystrophy with frontal predominance",
            "Macrocephaly (infantile form)",
            "Seizures",
            "Developmental regression",
            "Spasticity",
            "Rosenthal fiber accumulation (astrocytic inclusions)",
        ],
        "treatment": "Supportive care; no disease-modifying therapy; antisense oligonucleotide "
                     "approaches under investigation",
    },
    "pelizaeus_merzbacher": {
        "name": "Pelizaeus-Merzbacher disease",
        "gene": "PLP1",
        "inheritance": "X-linked recessive",
        "age_onset": "Neonatal to early infancy",
        "key_features": [
            "Nystagmus (often first sign)",
            "Hypomyelination on MRI",
            "Progressive spasticity",
            "Ataxia and tremor",
            "Cognitive impairment",
            "Stridor (laryngeal involvement)",
        ],
        "treatment": "Supportive care (physical therapy, spasticity management); "
                     "stem cell transplant and gene therapy investigational",
    },
    "giant_axonal_neuropathy": {
        "name": "Giant axonal neuropathy",
        "gene": "GAN",
        "inheritance": "autosomal recessive",
        "age_onset": "Early childhood (3–5 years)",
        "key_features": [
            "Progressive peripheral neuropathy (sensory and motor)",
            "CNS involvement (ataxia, dysarthria, intellectual disability)",
            "Characteristically kinky/curly hair",
            "Giant axonal swellings on nerve biopsy",
            "Cranial nerve involvement",
        ],
        "treatment": "Supportive care; intrathecal gene therapy (AAV9-GAN) in clinical trials",
    },
    "cockayne_syndrome": {
        "name": "Cockayne syndrome",
        "gene": "ERCC6 (CSB) / ERCC8 (CSA)",
        "inheritance": "autosomal recessive",
        "age_onset": "First 1–2 years of life",
        "key_features": [
            "Photosensitivity (without increased skin cancer risk)",
            "Dwarfism / cachectic appearance",
            "Progressive neurodegeneration (leukodystrophy, cerebellar atrophy)",
            "Sensorineural hearing loss",
            "Retinal dystrophy / cataracts",
            "Microcephaly",
            "Dental caries",
        ],
        "treatment": "Supportive care; sun protection, nutritional support, hearing aids, "
                     "dental care; no disease-modifying therapy",
    },
    "neuronal_ceroid_lipofuscinosis": {
        "name": "Neuronal ceroid lipofuscinosis (Batten disease)",
        "gene": "CLN3 / CLN5 / CLN6 / CLN8 (also PPT1, TPP1 for CLN1/CLN2)",
        "inheritance": "autosomal recessive",
        "age_onset": "Variable — infantile (CLN1), late infantile (CLN2), juvenile (CLN3)",
        "key_features": [
            "Progressive vision loss (retinal dystrophy)",
            "Seizures (often drug-resistant)",
            "Progressive cognitive decline and dementia",
            "Motor deterioration",
            "Behavioral changes",
            "Curvilinear or fingerprint profiles on electron microscopy",
        ],
        "treatment": "Cerliponase alfa (Brineura) — intracerebroventricular ERT for CLN2 disease; "
                     "supportive care for other forms; gene therapy investigational",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 4. CONNECTIVE TISSUE DISEASES (detailed)
# ═══════════════════════════════════════════════════════════════════════════════

CONNECTIVE_TISSUE_DISEASES: Dict[str, Dict[str, Any]] = {
    "marfan": {
        "name": "Marfan syndrome",
        "gene": "FBN1 (fibrillin-1)",
        "inheritance": "autosomal dominant",
        "key_features": [
            "Tall stature with disproportionately long limbs (dolichostenomelia)",
            "Arachnodactyly (long, slender fingers)",
            "Pectus excavatum or carinatum",
            "Ectopia lentis (lens subluxation)",
            "Aortic root dilation / dissection",
            "Mitral valve prolapse",
            "Dural ectasia",
        ],
        "diagnostic_criteria": "Revised Ghent nosology (2010) — aortic root Z-score "
                               "and systemic score with/without FBN1 variant",
        "treatment": "Beta-blockers or ARBs (losartan), aortic root monitoring, "
                     "prophylactic aortic root replacement, lens management",
        "omim_id": "154700",
    },
    "eds_classical": {
        "name": "Ehlers-Danlos syndrome, classical type",
        "gene": "COL5A1 / COL5A2",
        "inheritance": "autosomal dominant",
        "key_features": [
            "Skin hyperextensibility",
            "Atrophic scarring (cigarette-paper or papyraceous scars)",
            "Generalized joint hypermobility",
            "Easy bruising",
            "Smooth, velvety skin texture",
        ],
        "diagnostic_criteria": "2017 International EDS Classification — major and minor "
                               "criteria with molecular confirmation",
        "treatment": "Skin protection, joint stabilization, physical therapy, "
                     "wound care modifications",
        "omim_id": "130000",
    },
    "eds_vascular": {
        "name": "Ehlers-Danlos syndrome, vascular type",
        "gene": "COL3A1",
        "inheritance": "autosomal dominant",
        "key_features": [
            "Thin, translucent skin with visible veins",
            "Arterial/intestinal/uterine fragility and rupture",
            "Characteristic facial features (thin nose/lips, prominent eyes)",
            "Easy bruising",
            "Clubfoot (talipes equinovarus)",
            "Median age of first major event: 20-30 years",
        ],
        "diagnostic_criteria": "Clinical features plus confirmed COL3A1 pathogenic variant",
        "treatment": "Celiprolol (reduces vascular events), avoidance of invasive procedures, "
                     "surveillance imaging, emergency protocols",
        "omim_id": "130050",
    },
    "eds_hypermobile": {
        "name": "Ehlers-Danlos syndrome, hypermobile type",
        "gene": "Unknown (no confirmed gene as of 2026)",
        "inheritance": "Autosomal dominant (presumed)",
        "key_features": [
            "Generalized joint hypermobility (Beighton score >= 5)",
            "Chronic musculoskeletal pain",
            "Recurrent joint dislocations/subluxations",
            "Fatigue",
            "Dysautonomia (POTS-like symptoms)",
            "Functional GI disorders",
        ],
        "diagnostic_criteria": "2017 International EDS Classification clinical criteria "
                               "(no molecular test available)",
        "treatment": "Physical therapy (focus on joint stability), pain management, "
                     "autonomic symptom management, psychological support",
        "omim_id": "130020",
    },
    "oi_type_i": {
        "name": "Osteogenesis imperfecta type I (mild)",
        "gene": "COL1A1 (haploinsufficiency)",
        "inheritance": "autosomal dominant",
        "key_features": [
            "Fractures with mild to moderate trauma",
            "Blue sclerae",
            "Normal or near-normal stature",
            "Hearing loss (adult onset, 50%)",
            "Dentinogenesis imperfecta (some families)",
        ],
        "diagnostic_criteria": "Clinical features, bone density, collagen biochemistry, "
                               "COL1A1/COL1A2 sequencing",
        "treatment": "Bisphosphonates, fracture management, physical therapy, "
                     "hearing aids, orthopedic surgery",
        "omim_id": "166200",
    },
    "oi_type_ii": {
        "name": "Osteogenesis imperfecta type II (perinatal lethal)",
        "gene": "COL1A1 / COL1A2 (dominant-negative structural variants)",
        "inheritance": "autosomal dominant (usually de novo)",
        "key_features": [
            "Multiple fractures in utero",
            "Severe bone deformity",
            "Extremely soft calvarium",
            "Dark blue sclerae",
            "Pulmonary hypoplasia",
            "Usually lethal in perinatal period",
        ],
        "diagnostic_criteria": "Prenatal ultrasound, skeletal radiographs, molecular testing",
        "treatment": "Comfort care; most cases are lethal perinatally",
        "omim_id": "166210",
    },
    "oi_type_iii": {
        "name": "Osteogenesis imperfecta type III (severe, progressive)",
        "gene": "COL1A1 / COL1A2",
        "inheritance": "autosomal dominant (or AR in some forms)",
        "key_features": [
            "Progressive skeletal deformity",
            "Very short stature",
            "Fractures at birth and throughout life",
            "Triangular facies",
            "Variable blue/grey/white sclerae",
            "Dentinogenesis imperfecta (common)",
            "Basilar invagination risk",
        ],
        "diagnostic_criteria": "Clinical and radiographic findings, molecular confirmation",
        "treatment": "IV bisphosphonates (pamidronate/zoledronic acid), rodding surgery, "
                     "physical therapy, anti-sclerostin antibodies (investigational)",
        "omim_id": "259420",
    },
    "oi_type_iv": {
        "name": "Osteogenesis imperfecta type IV (moderate)",
        "gene": "COL1A1 / COL1A2",
        "inheritance": "autosomal dominant",
        "key_features": [
            "Moderate fracture frequency",
            "Mild to moderate short stature",
            "Normal to grey sclerae",
            "Bowing of long bones",
            "Dentinogenesis imperfecta (variable)",
        ],
        "diagnostic_criteria": "Clinical features, molecular testing",
        "treatment": "Bisphosphonates, orthopedic management, physical therapy",
        "omim_id": "166220",
    },
    "loeys_dietz": {
        "name": "Loeys-Dietz syndrome",
        "gene": "TGFBR1 / TGFBR2 / SMAD3 / TGFB2 / TGFB3",
        "inheritance": "autosomal dominant",
        "key_features": [
            "Arterial tortuosity and aneurysms (widespread, not limited to aortic root)",
            "Bifid uvula or cleft palate",
            "Hypertelorism (widely spaced eyes)",
            "Craniosynostosis (some subtypes)",
            "Joint hypermobility or contractures",
            "Aggressive vascular disease (dissection at smaller diameters than Marfan)",
        ],
        "diagnostic_criteria": "Clinical triad (arterial aneurysm/tortuosity + hypertelorism + "
                               "bifid uvula/cleft palate) plus molecular confirmation",
        "treatment": "ARBs (losartan), beta-blockers, early prophylactic vascular surgery "
                     "(lower surgical thresholds than Marfan), head-to-toe vascular imaging",
        "omim_id": "609192",
    },
    "stickler": {
        "name": "Stickler syndrome",
        "gene": "COL2A1 / COL11A1 / COL11A2 / COL9A1 / COL9A2 / COL9A3",
        "inheritance": "Autosomal dominant (most) or autosomal recessive",
        "key_features": [
            "High myopia (early onset, severe)",
            "Vitreoretinal degeneration with retinal detachment risk",
            "Midface hypoplasia (Pierre Robin sequence in some)",
            "Sensorineural and/or conductive hearing loss",
            "Joint hypermobility and early-onset arthropathy",
            "Cleft palate (some subtypes)",
        ],
        "diagnostic_criteria": "Clinical features across ocular, orofacial, auditory, "
                               "and skeletal systems plus molecular confirmation",
        "treatment": "Prophylactic retinal cryotherapy/laser, hearing aids, "
                     "joint physiotherapy, orthopedic management",
        "omim_id": "108300",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 5. HEMATOLOGIC DISEASES (detailed)
# ═══════════════════════════════════════════════════════════════════════════════

HEMATOLOGIC_DISEASES: Dict[str, Dict[str, Any]] = {
    "sickle_cell": {
        "name": "Sickle cell disease",
        "gene": "HBB (p.Glu6Val — HbS)",
        "inheritance": "autosomal recessive",
        "key_features": [
            "Vaso-occlusive pain crises",
            "Chronic hemolytic anemia",
            "Acute chest syndrome",
            "Stroke risk",
            "Splenic sequestration (children)",
            "Organ damage (kidney, retina, bone — avascular necrosis)",
        ],
        "treatment": "Hydroxyurea, L-glutamine, crizanlizumab (anti-P-selectin), "
                     "voxelotor (HbS polymerization inhibitor), blood transfusion, "
                     "hematopoietic stem cell transplant, gene therapy (exagamglogene/Casgevy, "
                     "lovotibeglogene/Lyfgenia)",
        "omim_id": "603903",
    },
    "beta_thalassemia": {
        "name": "Beta-thalassemia",
        "gene": "HBB (300+ pathogenic variants)",
        "inheritance": "autosomal recessive",
        "key_features": [
            "Microcytic hypochromic anemia",
            "Ineffective erythropoiesis",
            "Iron overload (transfusion-dependent)",
            "Hepatosplenomegaly",
            "Skeletal changes (frontal bossing, maxillary hyperplasia)",
            "Growth retardation",
        ],
        "treatment": "Regular transfusions (TDT), iron chelation (deferasirox, deferoxamine, "
                     "deferiprone), luspatercept, gene therapy (betibeglogene/Zynteglo, "
                     "exagamglogene/Casgevy), HSCT",
        "omim_id": "613985",
    },
    "hemophilia_a": {
        "name": "Hemophilia A",
        "gene": "F8 (coagulation factor VIII)",
        "inheritance": "X-linked recessive",
        "key_features": [
            "Recurrent hemarthrosis (joint bleeding)",
            "Deep muscle hematomas",
            "Prolonged bleeding after surgery/trauma",
            "Intracranial hemorrhage risk",
            "Target joint development (chronic arthropathy)",
        ],
        "treatment": "Factor VIII replacement (recombinant or plasma-derived), "
                     "emicizumab (bispecific antibody, factor VIII mimetic), "
                     "fitusiran (anti-antithrombin siRNA), gene therapy "
                     "(valoctocogene roxaparvovec/Roctavian)",
        "omim_id": "306700",
    },
    "hemophilia_b": {
        "name": "Hemophilia B (Christmas disease)",
        "gene": "F9 (coagulation factor IX)",
        "inheritance": "X-linked recessive",
        "key_features": [
            "Clinically indistinguishable from hemophilia A",
            "Recurrent hemarthrosis",
            "Deep tissue bleeding",
            "Prolonged aPTT with normal PT",
        ],
        "treatment": "Factor IX replacement, emicizumab (off-label), gene therapy "
                     "(etranacogene dezaparvovec/Hemgenix, fidanacogene elaparvovec/Beqvez)",
        "omim_id": "306900",
    },
    "von_willebrand": {
        "name": "Von Willebrand disease",
        "gene": "VWF",
        "inheritance": "Autosomal dominant (types 1, 2A, 2B, 2M) or recessive (type 3)",
        "key_features": [
            "Mucocutaneous bleeding (epistaxis, menorrhagia, gingival bleeding)",
            "Easy bruising",
            "Post-surgical bleeding",
            "Type 3: severe — joint and muscle bleeding similar to hemophilia",
        ],
        "treatment": "Desmopressin (DDAVP — types 1 and 2A), VWF-containing concentrates "
                     "(Humate-P, Vonvendi), tranexamic acid",
        "omim_id": "193400",
    },
    "diamond_blackfan": {
        "name": "Diamond-Blackfan anemia",
        "gene": "RPS19 (25%), RPL5, RPL11, RPS26, and other ribosomal protein genes",
        "inheritance": "autosomal dominant (most)",
        "key_features": [
            "Pure red cell aplasia (macrocytic anemia in infancy)",
            "Congenital anomalies (50%) — thumb, cardiac, craniofacial",
            "Short stature",
            "Increased cancer risk (MDS, AML, osteosarcoma)",
        ],
        "treatment": "Corticosteroids (first-line, ~80% respond), chronic transfusion "
                     "(steroid-refractory), hematopoietic stem cell transplant, "
                     "leucine supplementation (investigational)",
        "omim_id": "105650",
    },
    "fanconi_anemia": {
        "name": "Fanconi anemia",
        "gene": "FANCA (60-70%), FANCC, FANCG, and 20+ other FA complementation group genes",
        "inheritance": "autosomal recessive (most); X-linked (FANCB)",
        "key_features": [
            "Progressive bone marrow failure (pancytopenia)",
            "Congenital anomalies (radial ray defects, short stature, skin pigmentation)",
            "Cancer predisposition (AML, MDS, head/neck squamous cell carcinoma)",
            "Chromosomal breakage on DEB/MMC testing",
            "Endocrine abnormalities",
        ],
        "treatment": "Androgens (oxymetholone), HSCT (definitive for marrow failure), "
                     "cancer surveillance, gene therapy (investigational)",
        "omim_id": "227650",
    },
    "hereditary_spherocytosis": {
        "name": "Hereditary spherocytosis",
        "gene": "ANK1 / SLC4A1 / SPTA1 / SPTB / EPB41 / EPB42",
        "inheritance": "Autosomal dominant (75%) or recessive",
        "key_features": [
            "Hemolytic anemia (variable severity)",
            "Jaundice / hyperbilirubinemia",
            "Splenomegaly",
            "Gallstones",
            "Spherocytes on peripheral blood smear",
            "Osmotic fragility increased",
        ],
        "treatment": "Folic acid supplementation, splenectomy (severe cases), "
                     "cholecystectomy for gallstones, transfusion for aplastic crises",
        "omim_id": "182900",
    },
    "hereditary_hemochromatosis": {
        "name": "Hereditary hemochromatosis",
        "gene": "HFE (C282Y, H63D variants)",
        "inheritance": "autosomal recessive",
        "key_features": [
            "Iron overload affecting liver, heart, pancreas, joints, skin",
            "Fatigue and arthralgia (early symptoms)",
            "Elevated transferrin saturation and ferritin",
            "Bronze skin pigmentation",
            "Diabetes mellitus",
            "Cirrhosis and hepatocellular carcinoma risk",
        ],
        "treatment": "Therapeutic phlebotomy (first-line), iron chelation "
                     "(if phlebotomy-intolerant), dietary modification",
        "omim_id": "235200",
    },
    "paroxysmal_nocturnal_hemoglobinuria": {
        "name": "Paroxysmal nocturnal hemoglobinuria",
        "gene": "PIGA (somatic mutation in hematopoietic stem cell)",
        "inheritance": "Acquired (somatic); not inherited",
        "key_features": [
            "Intravascular hemolysis",
            "Hemoglobinuria (dark morning urine)",
            "Thrombosis (venous, often unusual sites — hepatic, cerebral)",
            "Bone marrow failure / aplastic anemia overlap",
            "Flow cytometry: loss of GPI-anchored proteins (CD55, CD59)",
        ],
        "treatment": "Complement inhibitors (eculizumab, ravulizumab — anti-C5; "
                     "iptacopan — factor B inhibitor; danicopan — factor D inhibitor), "
                     "anticoagulation, HSCT",
        "omim_id": "300818",
    },
    "pyruvate_kinase_deficiency": {
        "name": "Pyruvate kinase deficiency",
        "gene": "PKLR",
        "inheritance": "autosomal recessive",
        "key_features": [
            "Chronic hemolytic anemia (variable severity)",
            "Neonatal jaundice (may require exchange transfusion)",
            "Splenomegaly",
            "Iron overload (even without transfusions)",
            "Reticulocytosis",
            "Gallstones",
        ],
        "treatment": "Mitapivat (first-in-class PK activator, approved 2022), "
                     "folic acid, splenectomy (select cases), transfusion support, "
                     "iron chelation",
        "omim_id": "266200",
    },
    "thalassemia_intermedia": {
        "name": "Beta-thalassemia intermedia",
        "gene": "HBB",
        "inheritance": "autosomal recessive (variable severity based on genotype)",
        "key_features": [
            "Moderate anemia (Hb 7-10 g/dL typically)",
            "Ineffective erythropoiesis with extramedullary hematopoiesis",
            "Hepatosplenomegaly",
            "Iron overload (from increased GI absorption)",
            "Skeletal changes (less severe than thalassemia major)",
            "Thrombotic risk",
        ],
        "treatment": "Luspatercept (erythroid maturation agent), intermittent transfusions, "
                     "iron chelation, hydroxyurea (to increase HbF), splenectomy (select cases)",
        "omim_id": "613985",
    },
    "congenital_dyserythropoietic_anemia": {
        "name": "Congenital dyserythropoietic anemia",
        "gene": "SEC23B (type II) / CDAN1 (type I)",
        "inheritance": "autosomal recessive",
        "key_features": [
            "Ineffective erythropoiesis",
            "Morphologically abnormal erythroblasts in bone marrow",
            "Macrocytic or normocytic anemia",
            "Jaundice and splenomegaly",
            "Iron overload",
            "Type II (HEMPAS): positive acidified serum lysis test",
        ],
        "treatment": "Interferon-alpha (type I), splenectomy (type II), iron chelation, "
                     "transfusion support",
        "omim_id": "224120",
    },
    "hereditary_elliptocytosis": {
        "name": "Hereditary elliptocytosis",
        "gene": "SLC4A1 / SPTA1 / EPB41",
        "inheritance": "autosomal dominant",
        "key_features": [
            "Elliptical (oval) red blood cells on peripheral smear",
            "Mild hemolysis (most cases asymptomatic)",
            "Neonatal poikilocytosis may be severe (hereditary pyropoikilocytosis variant)",
            "Splenomegaly (in symptomatic cases)",
        ],
        "treatment": "Observation (most cases), folic acid, splenectomy (severe hemolysis), "
                     "transfusion for aplastic crises",
        "omim_id": "611804",
    },
    "thrombotic_thrombocytopenic_purpura": {
        "name": "Congenital thrombotic thrombocytopenic purpura (Upshaw-Schulman syndrome)",
        "gene": "ADAMTS13",
        "inheritance": "autosomal recessive",
        "key_features": [
            "Thrombocytopenia (severe, episodic)",
            "Microangiopathic hemolytic anemia (schistocytes)",
            "Organ ischemia (renal, neurological)",
            "Severely deficient ADAMTS13 activity (<10%)",
            "Triggered by infection, pregnancy, or surgery",
        ],
        "treatment": "Caplacizumab (anti-VWF nanobody), plasma infusion/exchange, "
                     "rituximab (acquired form), prophylactic plasma infusions (congenital)",
        "omim_id": "274150",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 6. IMMUNOLOGIC DISEASES (detailed)
# ═══════════════════════════════════════════════════════════════════════════════

IMMUNOLOGIC_DISEASES: Dict[str, Dict[str, Any]] = {
    "scid": {
        "name": "Severe combined immunodeficiency",
        "gene": "IL2RG (X-linked, ~50%), JAK3, IL7R, RAG1, RAG2, DCLRE1C (Artemis), ADA",
        "inheritance": "X-linked recessive (most common) or autosomal recessive",
        "key_features": [
            "Onset in first months of life",
            "Severe recurrent infections (bacterial, viral, fungal, opportunistic)",
            "Failure to thrive",
            "Chronic diarrhea",
            "Absent or very low T cells",
            "Absent thymic shadow on chest X-ray",
            "Fatal without treatment (usually by age 1-2)",
        ],
        "treatment": "Hematopoietic stem cell transplant (curative), gene therapy "
                     "(ADA-SCID: Strimvelis; X-SCID: investigational), "
                     "enzyme replacement (pegademase bovine for ADA-SCID)",
        "omim_id": "300400",
    },
    "cgd": {
        "name": "Chronic granulomatous disease",
        "gene": "CYBB (X-linked, 65%), CYBA, NCF1, NCF2, NCF4",
        "inheritance": "X-linked recessive (most common) or autosomal recessive",
        "key_features": [
            "Recurrent severe bacterial and fungal infections",
            "Granuloma formation (skin, liver, lungs, GI tract)",
            "Lymphadenopathy and hepatosplenomegaly",
            "Abscesses (liver, lung, perianal)",
            "Abnormal DHR (dihydrorhodamine) flow cytometry",
            "Catalase-positive organisms (Staphylococcus, Aspergillus, Serratia)",
        ],
        "treatment": "Antimicrobial prophylaxis (TMP-SMX, itraconazole), "
                     "interferon-gamma, HSCT (curative), gene therapy (investigational)",
        "omim_id": "306400",
    },
    "hyper_ige": {
        "name": "Hyper-IgE syndrome (Job syndrome)",
        "gene": "STAT3 (AD form) / DOCK8 (AR form)",
        "inheritance": "Autosomal dominant (STAT3) or autosomal recessive (DOCK8)",
        "key_features": [
            "Markedly elevated serum IgE (>2000 IU/mL)",
            "Recurrent staphylococcal skin abscesses (cold abscesses)",
            "Recurrent pneumonias with pneumatocele formation",
            "Eczema-like dermatitis",
            "Characteristic facies (prominent forehead, wide nasal bridge)",
            "Skeletal anomalies (scoliosis, retained primary teeth)",
            "DOCK8: severe viral skin infections, food allergies, malignancy risk",
        ],
        "treatment": "Antimicrobial prophylaxis (TMP-SMX, antifungals), "
                     "skin care, HSCT for DOCK8 deficiency (curative)",
        "omim_id": "147060",
    },
    "cvid": {
        "name": "Common variable immunodeficiency",
        "gene": "TNFRSF13B (TACI), ICOS, BAFF-R, and others (most cases polygenic/unknown)",
        "inheritance": "Variable (most sporadic; some AD or AR)",
        "key_features": [
            "Recurrent sinopulmonary infections",
            "Low IgG with low IgA and/or IgM",
            "Poor vaccine responses",
            "Autoimmune cytopenias (ITP, AIHA in 20-30%)",
            "Granulomatous disease (lungs, liver, spleen)",
            "Lymphoproliferation / lymphoma risk",
            "Onset typically in second or third decade",
        ],
        "treatment": "Immunoglobulin replacement therapy (IV or subcutaneous), "
                     "antimicrobial prophylaxis, management of autoimmune complications",
        "omim_id": "607594",
    },
    "wiskott_aldrich": {
        "name": "Wiskott-Aldrich syndrome",
        "gene": "WAS (WASP protein)",
        "inheritance": "X-linked recessive",
        "key_features": [
            "Thrombocytopenia with small platelets (pathognomonic)",
            "Eczema",
            "Recurrent infections (bacterial, viral, opportunistic)",
            "Autoimmunity (autoimmune hemolytic anemia, vasculitis)",
            "Lymphoma risk (especially EBV-associated)",
            "Bloody diarrhea in infancy",
        ],
        "treatment": "HSCT (curative — best outcomes with matched sibling donor), "
                     "gene therapy (investigational, promising results), "
                     "splenectomy (improves platelet count but increases infection risk)",
        "omim_id": "301000",
    },
    "xla": {
        "name": "X-linked agammaglobulinemia (Bruton agammaglobulinemia)",
        "gene": "BTK (Bruton tyrosine kinase)",
        "inheritance": "X-linked recessive",
        "key_features": [
            "Onset after 6 months (waning maternal antibodies)",
            "Recurrent sinopulmonary infections (encapsulated bacteria)",
            "Very low or absent immunoglobulins (all classes)",
            "Absent or markedly reduced B cells (<2% CD19+)",
            "Absent/small tonsils and lymph nodes",
            "Susceptibility to enteroviral meningoencephalitis",
        ],
        "treatment": "Lifelong immunoglobulin replacement therapy (IV or subcutaneous), "
                     "prompt antibiotic treatment of infections",
        "omim_id": "300755",
    },
    "ada_scid": {
        "name": "Adenosine deaminase deficiency SCID",
        "gene": "ADA",
        "inheritance": "autosomal recessive",
        "key_features": [
            "Severe combined immunodeficiency (T-B-NK- phenotype)",
            "Skeletal abnormalities (costochondral junction widening, pelvic dysplasia)",
            "Sensorineural hearing loss",
            "Hepatic and pulmonary dysfunction",
            "Elevated deoxyadenosine/dAXP metabolites",
        ],
        "treatment": "PEG-ADA enzyme replacement (pegademase, elapegademase), "
                     "gene therapy (Strimvelis — first approved gene therapy for PID), HSCT",
        "omim_id": "102700",
    },
    "ipex": {
        "name": "IPEX syndrome (immune dysregulation, polyendocrinopathy, "
                "enteropathy, X-linked)",
        "gene": "FOXP3",
        "inheritance": "X-linked recessive",
        "key_features": [
            "Severe intractable diarrhea (autoimmune enteropathy)",
            "Type 1 diabetes mellitus (neonatal onset)",
            "Eczematous dermatitis",
            "Autoimmune cytopenias",
            "Thyroiditis",
            "Absent or dysfunctional regulatory T cells (Tregs)",
        ],
        "treatment": "Immunosuppression (rapamycin/sirolimus, calcineurin inhibitors), "
                     "HSCT (curative)",
        "omim_id": "304790",
    },
    "dock8_deficiency": {
        "name": "DOCK8 deficiency",
        "gene": "DOCK8",
        "inheritance": "autosomal recessive",
        "key_features": [
            "Recurrent severe viral skin infections (HSV, HPV, molluscum contagiosum)",
            "Elevated serum IgE",
            "Eosinophilia",
            "Severe atopic dermatitis",
            "Food allergies and anaphylaxis",
            "Susceptibility to malignancy (squamous cell carcinoma, lymphoma)",
            "T-cell lymphopenia over time",
        ],
        "treatment": "HSCT (curative and recommended early), antimicrobial prophylaxis, "
                     "IgG replacement, careful infection management",
        "omim_id": "243700",
    },
    "stat3_gof": {
        "name": "STAT3 gain-of-function disease",
        "gene": "STAT3",
        "inheritance": "autosomal dominant (gain-of-function)",
        "key_features": [
            "Early-onset multi-organ autoimmunity",
            "Lymphoproliferation",
            "Autoimmune cytopenias (AIHA, ITP)",
            "Type 1 diabetes (early onset)",
            "Autoimmune enteropathy",
            "Short stature",
            "Interstitial lung disease",
        ],
        "treatment": "JAK inhibitors (ruxolitinib, tofacitinib), immunosuppression, "
                     "HSCT in refractory cases",
        "omim_id": "615952",
    },
    "ctla4_haploinsufficiency": {
        "name": "CTLA-4 haploinsufficiency with autoimmune infiltration",
        "gene": "CTLA4",
        "inheritance": "autosomal dominant",
        "key_features": [
            "Autoimmune cytopenias",
            "Lymphocytic infiltration of organs (lungs, brain, gut)",
            "Hypogammaglobulinemia",
            "Lymphadenopathy and splenomegaly",
            "Recurrent respiratory infections",
            "Variable penetrance within families",
        ],
        "treatment": "Abatacept (CTLA4-Ig fusion protein — targeted therapy), "
                     "immunosuppression, immunoglobulin replacement, HSCT in severe cases",
        "omim_id": "616100",
    },
    "lrba_deficiency": {
        "name": "LRBA deficiency",
        "gene": "LRBA",
        "inheritance": "autosomal recessive",
        "key_features": [
            "CTLA4-related immune dysregulation (LRBA regulates CTLA4 trafficking)",
            "Autoimmune cytopenias (AIHA, ITP, autoimmune neutropenia)",
            "Autoimmune enteropathy",
            "Hypogammaglobulinemia",
            "Lymphoproliferation",
            "Organomegaly",
        ],
        "treatment": "Abatacept (targets the CTLA4 pathway deficiency), "
                     "immunosuppression (sirolimus), immunoglobulin replacement, HSCT",
        "omim_id": "614700",
    },
    "activated_pi3k_delta": {
        "name": "Activated PI3K-delta syndrome (APDS)",
        "gene": "PIK3CD (APDS1) / PIK3R1 (APDS2)",
        "inheritance": "autosomal dominant",
        "key_features": [
            "Recurrent sinopulmonary infections",
            "Lymphoproliferation (lymphadenopathy, splenomegaly)",
            "Herpesvirus susceptibility (EBV, CMV)",
            "Autoimmune cytopenias",
            "Bronchiectasis",
            "Increased lymphoma risk",
        ],
        "treatment": "Leniolisib (first approved PI3Kδ inhibitor, 2023), "
                     "immunoglobulin replacement, mTOR inhibitors (rapamycin), HSCT",
        "omim_id": "615513",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 7. CANCER PREDISPOSITION SYNDROMES (detailed)
# ═══════════════════════════════════════════════════════════════════════════════

CANCER_PREDISPOSITION: Dict[str, Dict[str, Any]] = {
    "li_fraumeni": {
        "name": "Li-Fraumeni syndrome",
        "gene": "TP53",
        "inheritance": "autosomal dominant",
        "lifetime_cancer_risk": ">90% by age 70",
        "key_cancers": [
            "Soft tissue sarcoma",
            "Osteosarcoma",
            "Breast cancer (premenopausal)",
            "Brain tumors (glioma, choroid plexus carcinoma)",
            "Adrenocortical carcinoma (childhood)",
            "Leukemia",
        ],
        "surveillance": "Whole-body MRI (annual), breast MRI (age 20+), brain MRI, "
                        "abdominal ultrasound, biochemical markers",
        "omim_id": "151623",
    },
    "lynch": {
        "name": "Lynch syndrome (hereditary non-polyposis colorectal cancer)",
        "gene": "MLH1 / MSH2 / MSH6 / PMS2 / EPCAM",
        "inheritance": "autosomal dominant",
        "lifetime_cancer_risk": "40-80% colorectal; 25-60% endometrial",
        "key_cancers": [
            "Colorectal cancer",
            "Endometrial cancer",
            "Ovarian cancer",
            "Gastric cancer",
            "Urothelial cancer",
            "Hepatobiliary cancer",
            "Sebaceous neoplasms (Muir-Torre variant)",
        ],
        "surveillance": "Colonoscopy every 1-2 years from age 25, endometrial sampling, "
                        "urinalysis, upper GI endoscopy (select genes)",
        "omim_id": "120435",
    },
    "brca": {
        "name": "Hereditary breast and ovarian cancer syndrome",
        "gene": "BRCA1 / BRCA2",
        "inheritance": "autosomal dominant",
        "lifetime_cancer_risk": "45-85% breast; 10-46% ovarian (BRCA1 > BRCA2)",
        "key_cancers": [
            "Breast cancer (female and male)",
            "Ovarian/fallopian tube/peritoneal cancer",
            "Prostate cancer (BRCA2)",
            "Pancreatic cancer (BRCA2)",
        ],
        "surveillance": "Breast MRI and mammography alternating every 6 months from age 25, "
                        "consider risk-reducing mastectomy and/or salpingo-oophorectomy",
        "omim_id": "604370",
    },
    "fap": {
        "name": "Familial adenomatous polyposis",
        "gene": "APC",
        "inheritance": "autosomal dominant",
        "lifetime_cancer_risk": "~100% colorectal without intervention",
        "key_cancers": [
            "Colorectal cancer (from hundreds-thousands of adenomatous polyps)",
            "Duodenal/periampullary cancer",
            "Thyroid cancer (papillary, cribriform-morular variant)",
            "Hepatoblastoma (childhood)",
            "Desmoid tumors",
        ],
        "surveillance": "Sigmoidoscopy/colonoscopy from age 10-12, upper GI endoscopy, "
                        "thyroid ultrasound; prophylactic colectomy when polyp burden high",
        "omim_id": "175100",
    },
    "men1": {
        "name": "Multiple endocrine neoplasia type 1",
        "gene": "MEN1 (menin)",
        "inheritance": "autosomal dominant",
        "lifetime_cancer_risk": ">90% penetrance by age 40",
        "key_cancers": [
            "Parathyroid adenoma/hyperplasia (95%) — primary hyperparathyroidism",
            "Pituitary adenoma (prolactinoma most common)",
            "Pancreatic neuroendocrine tumors (gastrinoma, insulinoma)",
            "Adrenal cortical tumors",
            "Thymic/bronchial carcinoid",
        ],
        "surveillance": "Annual calcium/PTH, pituitary MRI, abdominal CT/MRI, "
                        "fasting gut hormones (gastrin, insulin, glucagon)",
        "omim_id": "131100",
    },
    "men2": {
        "name": "Multiple endocrine neoplasia type 2",
        "gene": "RET (proto-oncogene)",
        "inheritance": "autosomal dominant",
        "lifetime_cancer_risk": ">95% medullary thyroid carcinoma (MEN2A/2B)",
        "key_cancers": [
            "Medullary thyroid carcinoma (virtually 100%)",
            "Pheochromocytoma (50% MEN2A, >50% MEN2B)",
            "Parathyroid adenoma (20-30% MEN2A only)",
        ],
        "surveillance": "Prophylactic thyroidectomy (timing based on RET codon risk), "
                        "annual calcitonin, plasma/urine metanephrines, calcium/PTH",
        "omim_id": "171400",
    },
    "vhl": {
        "name": "Von Hippel-Lindau syndrome",
        "gene": "VHL",
        "inheritance": "autosomal dominant",
        "lifetime_cancer_risk": ">90% penetrance by age 65",
        "key_cancers": [
            "Clear cell renal cell carcinoma (up to 70%)",
            "Hemangioblastoma (CNS: cerebellum, spinal cord; retina)",
            "Pheochromocytoma",
            "Pancreatic neuroendocrine tumors / cysts",
            "Endolymphatic sac tumors",
            "Epididymal/broad ligament cystadenoma",
        ],
        "surveillance": "Annual ophthalmologic exam, abdominal MRI/CT, brain/spine MRI, "
                        "plasma metanephrines, audiologic assessment",
        "omim_id": "193300",
    },
    "rb1": {
        "name": "Hereditary retinoblastoma",
        "gene": "RB1",
        "inheritance": "autosomal dominant (high penetrance)",
        "lifetime_cancer_risk": ">90% retinoblastoma (bilateral); second malignancies ~35%",
        "key_cancers": [
            "Retinoblastoma (bilateral, multifocal)",
            "Osteosarcoma",
            "Soft tissue sarcoma",
            "Melanoma",
            "Pineoblastoma (trilateral retinoblastoma)",
        ],
        "surveillance": "Frequent dilated fundoscopic exams under anesthesia from birth, "
                        "whole-body MRI for second malignancies, genetic counseling",
        "omim_id": "180200",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 8. APPROVED GENE THERAPIES
# ═══════════════════════════════════════════════════════════════════════════════

GENE_THERAPY_APPROVED: Dict[str, Dict[str, Any]] = {
    "nusinersen": {
        "drug_name": "nusinersen",
        "brand_name": "Spinraza",
        "gene_target": "SMN2 (modulates splicing to increase full-length SMN protein)",
        "disease": "Spinal muscular atrophy (all types)",
        "mechanism": "Antisense oligonucleotide — binds SMN2 pre-mRNA intron 7 "
                     "to promote exon 7 inclusion, producing functional SMN protein",
        "approval_year": 2016,
        "route": "Intrathecal injection (loading doses then maintenance every 4 months)",
        "cost_estimate": "~$750,000 first year, ~$375,000/year thereafter",
    },
    "onasemnogene": {
        "drug_name": "onasemnogene abeparvovec",
        "brand_name": "Zolgensma",
        "gene_target": "SMN1 (delivers functional copy via AAV9 vector)",
        "disease": "Spinal muscular atrophy type 1 (children <2 years)",
        "mechanism": "AAV9-based gene replacement — delivers a functional copy of "
                     "the SMN1 gene to motor neurons via a single IV infusion",
        "approval_year": 2019,
        "route": "Single intravenous infusion",
        "cost_estimate": "~$2.1 million (single dose)",
    },
    "risdiplam": {
        "drug_name": "risdiplam",
        "brand_name": "Evrysdi",
        "gene_target": "SMN2 (small molecule splicing modifier)",
        "disease": "Spinal muscular atrophy (all types, ages 2 months and older)",
        "mechanism": "Oral small molecule SMN2 splicing modifier — promotes exon 7 "
                     "inclusion in SMN2 mRNA, increasing functional SMN protein "
                     "production systemically",
        "approval_year": 2020,
        "route": "Oral solution (daily)",
        "cost_estimate": "~$340,000/year (ongoing)",
    },
    "voretigene": {
        "drug_name": "voretigene neparvovec",
        "brand_name": "Luxturna",
        "gene_target": "RPE65 (delivers functional copy via AAV2 vector)",
        "disease": "Biallelic RPE65-mediated inherited retinal dystrophy "
                   "(Leber congenital amaurosis type 2, retinitis pigmentosa type 20)",
        "mechanism": "AAV2-based gene replacement — subretinal injection delivers "
                     "functional RPE65 gene to retinal pigment epithelium cells, "
                     "restoring the visual cycle",
        "approval_year": 2017,
        "route": "Subretinal injection (one injection per eye)",
        "cost_estimate": "~$850,000 (both eyes, one-time)",
    },
    "etranacogene": {
        "drug_name": "etranacogene dezaparvovec",
        "brand_name": "Hemgenix",
        "gene_target": "F9 (delivers Padua gain-of-function variant via AAV5 vector)",
        "disease": "Hemophilia B (adults with factor IX activity <= 2%)",
        "mechanism": "AAV5-based gene therapy — delivers a codon-optimized factor IX "
                     "Padua variant transgene to hepatocytes, enabling endogenous "
                     "factor IX production",
        "approval_year": 2022,
        "route": "Single intravenous infusion",
        "cost_estimate": "~$3.5 million (single dose)",
    },
    "exagamglogene": {
        "drug_name": "exagamglogene autotemcel",
        "brand_name": "Casgevy",
        "gene_target": "BCL11A enhancer (CRISPR-Cas9 editing of autologous HSCs)",
        "disease": "Sickle cell disease and transfusion-dependent beta-thalassemia "
                   "(ages 12 and older)",
        "mechanism": "Ex-vivo CRISPR-Cas9 gene editing — patient's HSCs are edited "
                     "to disrupt the BCL11A erythroid enhancer, reactivating fetal "
                     "hemoglobin (HbF) production; edited cells are reinfused after "
                     "myeloablative conditioning",
        "approval_year": 2023,
        "route": "Autologous HSCT (single infusion after myeloablative conditioning)",
        "cost_estimate": "~$2.2 million (single procedure)",
    },
    "valoctocogene": {
        "drug_name": "valoctocogene roxaparvovec",
        "brand_name": "Roctavian",
        "gene_target": "F8 (delivers functional FVIII via AAV5 vector)",
        "disease": "Hemophilia A (severe, adults without pre-existing AAV5 antibodies)",
        "mechanism": "AAV5-based gene therapy — delivers a B-domain-deleted factor VIII "
                     "transgene to hepatocytes, enabling endogenous FVIII production",
        "approval_year": 2022,
        "route": "Single intravenous infusion",
        "cost_estimate": "~$2.9 million (single dose)",
    },
    "delandistrogene": {
        "drug_name": "delandistrogene moxeparvovec",
        "brand_name": "Elevidys",
        "gene_target": "DMD (delivers micro-dystrophin via AAVrh74 vector)",
        "disease": "Duchenne muscular dystrophy (ambulatory patients aged 4-5 years)",
        "mechanism": "AAVrh74-based gene therapy — delivers a micro-dystrophin transgene "
                     "to skeletal muscle, providing a shortened but functional dystrophin protein",
        "approval_year": 2023,
        "route": "Single intravenous infusion",
        "cost_estimate": "~$3.2 million (single dose)",
    },
    "fidanacogene": {
        "drug_name": "fidanacogene elaparvovec",
        "brand_name": "Beqvez",
        "gene_target": "F9 (delivers high-activity FIX variant via AAV vector)",
        "disease": "Hemophilia B (moderately severe to severe, adults)",
        "mechanism": "AAV-based gene therapy — delivers a codon-optimized factor IX "
                     "Padua variant transgene to hepatocytes for sustained FIX expression",
        "approval_year": 2024,
        "route": "Single intravenous infusion",
        "cost_estimate": "~$3.5 million (single dose)",
    },
    "beremagene": {
        "drug_name": "beremagene geperpavec",
        "brand_name": "Vyjuvek",
        "gene_target": "COL7A1 (delivers functional collagen VII via HSV-1 vector)",
        "disease": "Dystrophic epidermolysis bullosa (DEB, ages 6 months and older)",
        "mechanism": "HSV-1-based gene therapy — topical application delivers functional "
                     "COL7A1 gene to skin keratinocytes and fibroblasts, restoring type VII "
                     "collagen at the dermal-epidermal junction",
        "approval_year": 2023,
        "route": "Topical (applied to wounds, repeated dosing)",
        "cost_estimate": "~$24,250 per mL (repeated application)",
    },
    "lovotibeglogene": {
        "drug_name": "lovotibeglogene autotemcel",
        "brand_name": "Lyfgenia",
        "gene_target": "HBB (lentiviral delivery of anti-sickling beta-globin to HSCs)",
        "disease": "Sickle cell disease (ages 12 and older with recurrent vaso-occlusive events)",
        "mechanism": "Ex-vivo lentiviral gene addition — patient's HSCs are transduced "
                     "with a lentiviral vector encoding an anti-sickling beta-globin "
                     "(βA-T87Q); modified cells reinfused after myeloablative conditioning",
        "approval_year": 2023,
        "route": "Autologous HSCT (single infusion after myeloablative conditioning)",
        "cost_estimate": "~$3.1 million (single procedure)",
    },
    "atidarsagene": {
        "drug_name": "atidarsagene autotemcel",
        "brand_name": "Libmeldy",
        "gene_target": "ARSA (lentiviral delivery of arylsulfatase A to HSCs)",
        "disease": "Metachromatic leukodystrophy (MLD, pre-symptomatic late infantile or "
                   "early juvenile forms)",
        "mechanism": "Ex-vivo lentiviral gene therapy — patient's HSCs are transduced with "
                     "a lentiviral vector encoding ARSA gene; modified cells engraft in bone "
                     "marrow and migrate to CNS as microglia, delivering functional enzyme",
        "approval_year": 2024,
        "route": "Autologous HSCT (single infusion after myeloablative conditioning)",
        "cost_estimate": "~$3.9 million (single procedure)",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 9. DIAGNOSTIC ALGORITHMS
# ═══════════════════════════════════════════════════════════════════════════════

DIAGNOSTIC_ALGORITHMS: Dict[str, Dict[str, Any]] = {
    "neurodevelopmental": {
        "name": "Neurodevelopmental Delay / Intellectual Disability Workup",
        "indication": "Global developmental delay, intellectual disability, autism "
                      "spectrum disorder with regression, unexplained seizures with "
                      "developmental regression",
        "test_sequence": [
            "1. Comprehensive clinical assessment and 3-generation pedigree",
            "2. Chromosomal microarray (CMA) — first-tier test for unexplained ID/ASD",
            "3. Fragile X testing (FMR1 CGG repeat) — if male or family history",
            "4. Metabolic screening (plasma amino acids, urine organic acids, "
            "acylcarnitine profile, urine purines/pyrimidines)",
            "5. Brain MRI with spectroscopy",
            "6. Whole-exome sequencing (WES) — if CMA and FraX negative",
            "7. Whole-genome sequencing (WGS) — if WES non-diagnostic",
            "8. RNA sequencing / transcriptome analysis — research setting",
            "9. Epigenomic testing (methylation array) — for imprinting disorders",
        ],
        "diagnostic_yield": "CMA: 15-20%; WES: 25-40%; WGS: 30-50%",
    },
    "metabolic_crisis": {
        "name": "Acute Metabolic Crisis / Neonatal Metabolic Workup",
        "indication": "Neonatal encephalopathy, recurrent vomiting, hyperammonemia, "
                      "metabolic acidosis with elevated anion gap, hypoglycemia with "
                      "hepatomegaly, abnormal newborn screen",
        "test_sequence": [
            "1. STAT labs: glucose, ammonia, lactate, blood gas (pH, bicarbonate, "
            "anion gap), liver function, electrolytes",
            "2. Plasma amino acids (quantitative)",
            "3. Urine organic acids (GC-MS)",
            "4. Acylcarnitine profile (tandem mass spectrometry)",
            "5. Urine reducing substances, urine ketones",
            "6. Plasma very-long-chain fatty acids (if peroxisomal disorder suspected)",
            "7. Enzyme assay (based on suspected diagnosis)",
            "8. Molecular confirmation (targeted gene or IEM gene panel)",
            "9. Functional studies if variant of uncertain significance",
        ],
        "diagnostic_yield": "Newborn screening + acute biochemistry: 70-85%; "
                            "molecular confirmation: 85-95%",
    },
    "skeletal_dysplasia": {
        "name": "Skeletal Dysplasia Diagnostic Pathway",
        "indication": "Disproportionate short stature, skeletal deformity, multiple "
                      "fractures with minimal trauma, abnormal prenatal ultrasound "
                      "showing shortened limbs",
        "test_sequence": [
            "1. Full skeletal survey (AP and lateral radiographs)",
            "2. Growth parameters and proportionality assessment",
            "3. Bone age",
            "4. Calcium, phosphate, alkaline phosphatase, PTH, vitamin D",
            "5. Classification per Nosology of Genetic Skeletal Disorders (2023 revision)",
            "6. Skeletal dysplasia gene panel (100-400 genes)",
            "7. Whole-exome sequencing if panel non-diagnostic",
            "8. Prenatal: cell-free DNA, CVS, or amniocentesis for molecular testing",
        ],
        "diagnostic_yield": "Skeletal survey + clinical: 60-70% narrow differential; "
                            "molecular: 75-90%",
    },
    "cardiac": {
        "name": "Inherited Cardiac Disease Workup",
        "indication": "Unexplained cardiomyopathy, sudden cardiac death in family, "
                      "early-onset arrhythmia, long QT or Brugada pattern on ECG, "
                      "aortic root dilation in young patient",
        "test_sequence": [
            "1. Detailed 3-generation family history (sudden death, syncope, ICD, pacemaker)",
            "2. 12-lead ECG and signal-averaged ECG",
            "3. Echocardiography (with strain imaging)",
            "4. Cardiac MRI (with late gadolinium enhancement)",
            "5. 24-48 hour Holter monitor",
            "6. Exercise stress test (for LQTS, CPVT)",
            "7. Cardiomyopathy or arrhythmia gene panel (30-200 genes)",
            "8. Cascade family screening (clinical + genotype-specific)",
            "9. Provocative drug testing if indicated (ajmaline for Brugada, "
            "epinephrine for LQTS type 1)",
        ],
        "diagnostic_yield": "HCM panels: 40-60%; DCM panels: 20-35%; "
                            "LQTS panels: 75-80%; Brugada: 20-30%",
    },
    "immunodeficiency": {
        "name": "Primary Immunodeficiency Diagnostic Pathway",
        "indication": "Recurrent severe infections, unusual organisms, family history "
                      "of immunodeficiency, failed newborn screening (TREC), "
                      "autoimmunity with infections",
        "test_sequence": [
            "1. Complete blood count with differential and peripheral smear",
            "2. Quantitative immunoglobulins (IgG, IgA, IgM, IgE)",
            "3. Lymphocyte subsets (CD3, CD4, CD8, CD19, CD16/56 by flow cytometry)",
            "4. Vaccine titers (tetanus, pneumococcal, Hib)",
            "5. DHR (dihydrorhodamine) flow cytometry — for CGD",
            "6. Complement studies (CH50, AH50)",
            "7. TREC/KREC (newborn screening or retrospective from Guthrie card)",
            "8. Targeted PID gene panel (300-400 genes) or WES",
            "9. Functional lymphocyte studies (mitogen stimulation) if indicated",
        ],
        "diagnostic_yield": "Targeted panels: 50-70% for suspected PID; "
                            "WES: 25-40% if panel negative",
    },
    "connective_tissue": {
        "name": "Heritable Connective Tissue Disorder Workup",
        "indication": "Joint hypermobility with complications, aortic root dilation "
                      "in young patient, ectopia lentis, recurrent fractures, "
                      "velvety/hyperextensible skin, arterial dissection in young adult",
        "test_sequence": [
            "1. Beighton score and systemic hypermobility assessment",
            "2. Echocardiography (aortic root Z-score, mitral valve prolapse)",
            "3. Ophthalmologic slit-lamp examination (lens subluxation)",
            "4. Skin biopsy with electron microscopy (EDS subtypes)",
            "5. Ghent criteria assessment (Marfan), 2017 EDS criteria, "
            "Sillence classification (OI)",
            "6. Connective tissue disorder gene panel (FBN1, COL1A1/A2, COL3A1, "
            "COL5A1/A2, TGFBR1/2, SMAD3, etc.)",
            "7. Whole-exome sequencing if panel non-diagnostic",
            "8. Collagen biochemistry (skin fibroblast culture) — select cases",
            "9. Vascular imaging (CTA or MRA — head-to-pelvis for Loeys-Dietz/vEDS)",
        ],
        "diagnostic_yield": "Marfan (FBN1 testing): 90-95%; vEDS (COL3A1): >95%; "
                            "classical EDS: 50-90%; hEDS: clinical diagnosis only",
    },
    "hematologic": {
        "name": "Hematologic / Inherited Anemia Diagnostic Pathway",
        "indication": "Unexplained chronic anemia, hemolysis (elevated LDH, low haptoglobin, "
                      "elevated indirect bilirubin), abnormal hemoglobin electrophoresis, "
                      "recurrent transfusion requirement, family history of blood disorder",
        "test_sequence": [
            "1. Complete blood count (CBC) with red cell indices and reticulocyte count",
            "2. Peripheral blood smear review (morphology: spherocytes, elliptocytes, "
            "sickle cells, target cells, schistocytes)",
            "3. Hemoglobin electrophoresis / HPLC (HbS, HbC, HbF, HbA2 quantification)",
            "4. Iron studies (ferritin, transferrin saturation, TIBC)",
            "5. Hemolysis markers (LDH, haptoglobin, indirect bilirubin, Coombs test)",
            "6. Osmotic fragility test / EMA binding (hereditary spherocytosis)",
            "7. Enzyme assays (G6PD, pyruvate kinase) if indicated",
            "8. Targeted hemoglobinopathy/red cell membrane gene panel",
            "9. Whole-exome sequencing if panel non-diagnostic",
        ],
        "diagnostic_yield": "Hemoglobin electrophoresis + CBC: 70-80% for hemoglobinopathies; "
                            "gene panels: 80-90% for hereditary anemias",
    },
    "ophthalmologic": {
        "name": "Inherited Retinal / Ophthalmic Disease Diagnostic Pathway",
        "indication": "Early-onset vision loss, nyctalopia (night blindness), visual field "
                      "constriction, family history of inherited retinal dystrophy, "
                      "abnormal electroretinography",
        "test_sequence": [
            "1. Comprehensive ophthalmic examination (visual acuity, slit lamp, fundoscopy)",
            "2. Electroretinography (ERG — full-field and multifocal)",
            "3. Optical coherence tomography (OCT — retinal layer assessment)",
            "4. Fundus autofluorescence imaging",
            "5. Visual field testing (Goldmann or Humphrey perimetry)",
            "6. Dark adaptometry (for rod function)",
            "7. Targeted inherited retinal disease gene panel (RPE65, RPGR, RS1, "
            "ABCA4, RHO, and 250+ genes)",
            "8. Whole-exome sequencing if panel non-diagnostic",
            "9. Whole-genome sequencing (for deep intronic and structural variants)",
        ],
        "diagnostic_yield": "Targeted IRD panels: 55-75%; WES: 60-80%; WGS: 70-85%",
    },
    "renal": {
        "name": "Inherited Renal Disease Diagnostic Pathway",
        "indication": "Early-onset chronic kidney disease, bilateral renal cysts, "
                      "persistent hematuria/proteinuria in young patient, renal Fanconi syndrome, "
                      "family history of kidney failure",
        "test_sequence": [
            "1. Renal ultrasound (cyst characterization, kidney size, echogenicity)",
            "2. Urinalysis with microscopy (hematuria, proteinuria, crystals, casts)",
            "3. Serum creatinine, BUN, electrolytes, GFR estimation",
            "4. Urine protein-to-creatinine ratio, albumin-to-creatinine ratio",
            "5. Kidney biopsy with light/electron microscopy and immunofluorescence "
            "(Alport — thin/thick GBM, lamellation)",
            "6. Slit-lamp examination (anterior lenticonus in Alport)",
            "7. Targeted renal disease gene panel (PKD1, PKD2, COL4A3, COL4A4, COL4A5, "
            "CTNS, NPHP1, and 100+ genes)",
            "8. Whole-exome sequencing if panel non-diagnostic",
        ],
        "diagnostic_yield": "Renal gene panels: 50-70% for suspected monogenic kidney disease; "
                            "WES: 25-40% if panel negative",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 10. ACMG/AMP VARIANT CLASSIFICATION CRITERIA
# ═══════════════════════════════════════════════════════════════════════════════

ACMG_CRITERIA: Dict[str, Dict[str, Any]] = {
    # --- Very Strong Pathogenic ---
    "PVS1": {
        "name": "Null variant in a gene with established LOF mechanism",
        "description": "Null variant (nonsense, frameshift, canonical +/-1 or 2 splice "
                       "sites, initiation codon, single or multi-exon deletion) in a gene "
                       "where loss of function is a known mechanism of disease.",
        "strength": "very_strong_pathogenic",
        "examples": ["Nonsense variant in NF1", "Frameshift in BRCA1", "Canonical splice "
                     "site variant in DMD"],
    },
    # --- Strong Pathogenic ---
    "PS1": {
        "name": "Same amino acid change as established pathogenic variant",
        "description": "Same amino acid change as a previously established pathogenic "
                       "variant regardless of nucleotide change.",
        "strength": "strong_pathogenic",
        "examples": ["Known pathogenic p.Arg553Ter in CFTR, novel nucleotide change "
                     "yielding same amino acid change"],
    },
    "PS2": {
        "name": "De novo (maternity and paternity confirmed)",
        "description": "De novo variant (both maternity and paternity confirmed) in a "
                       "patient with the disease and no family history.",
        "strength": "strong_pathogenic",
        "examples": ["De novo MECP2 variant in female with Rett syndrome"],
    },
    "PS3": {
        "name": "Well-established functional studies show damaging effect",
        "description": "Well-established in vitro or in vivo functional studies "
                       "supportive of a damaging effect on the gene or gene product.",
        "strength": "strong_pathogenic",
        "examples": ["Enzyme activity assay showing <5% residual activity"],
    },
    "PS4": {
        "name": "Prevalence significantly increased in affected vs. controls",
        "description": "The prevalence of the variant in affected individuals is "
                       "significantly increased compared with the prevalence in controls.",
        "strength": "strong_pathogenic",
        "examples": ["Variant found in 10 unrelated affected families but absent "
                     "in gnomAD"],
    },
    # --- Moderate Pathogenic ---
    "PM1": {
        "name": "Located in hot spot / functional domain without benign variation",
        "description": "Located in a mutational hot spot and/or critical and "
                       "well-established functional domain without benign variation.",
        "strength": "moderate_pathogenic",
        "examples": ["Variant in the GTP-binding domain of HRAS", "Variant in the "
                     "kinase domain of RET"],
    },
    "PM2": {
        "name": "Absent from controls (or extremely low frequency)",
        "description": "Absent from controls (or at extremely low frequency if recessive) "
                       "in Exome Sequencing Project, 1000 Genomes, ExAC/gnomAD.",
        "strength": "moderate_pathogenic",
        "examples": ["Variant absent from gnomAD (282,000+ alleles)"],
    },
    "PM3": {
        "name": "Detected in trans with a pathogenic variant (recessive)",
        "description": "For recessive disorders, detected in trans with a pathogenic "
                       "variant (phase confirmed).",
        "strength": "moderate_pathogenic",
        "examples": ["Novel missense in CFTR detected in trans with p.Phe508del"],
    },
    "PM4": {
        "name": "Protein length change (in-frame del/ins in non-repeat region)",
        "description": "Protein length changes as a result of in-frame deletions/insertions "
                       "in a non-repeat region or stop-loss variants.",
        "strength": "moderate_pathogenic",
        "examples": ["In-frame deletion of 3 amino acids in COL1A1 Gly-X-Y region"],
    },
    "PM5": {
        "name": "Novel missense at same position as established pathogenic",
        "description": "Novel missense change at an amino acid residue where a different "
                       "missense change determined to be pathogenic has been seen before.",
        "strength": "moderate_pathogenic",
        "examples": ["Novel missense at a TP53 codon where another missense is "
                     "established pathogenic"],
    },
    "PM6": {
        "name": "Assumed de novo (parental testing not performed)",
        "description": "Assumed de novo, but without confirmation of paternity "
                       "and maternity.",
        "strength": "moderate_pathogenic",
        "examples": ["New variant in proband not seen in parents (paternity assumed)"],
    },
    # --- Supporting Pathogenic ---
    "PP1": {
        "name": "Cosegregation with disease in family",
        "description": "Cosegregation with disease in multiple affected family members "
                       "in a gene definitively known to cause the disease.",
        "strength": "supporting_pathogenic",
        "examples": ["Variant segregates with HCM in 5 affected family members, LOD >3"],
    },
    "PP2": {
        "name": "Missense in gene with low rate of benign missense",
        "description": "Missense variant in a gene that has a low rate of benign missense "
                       "variation and in which missense variants are a common mechanism of disease.",
        "strength": "supporting_pathogenic",
        "examples": ["Missense variant in MYH7 (cardiomyopathy gene)"],
    },
    "PP3": {
        "name": "Multiple computational evidence supports deleterious effect",
        "description": "Multiple lines of computational evidence support a deleterious "
                       "effect on the gene or gene product (REVEL, CADD, SIFT, PolyPhen, etc.).",
        "strength": "supporting_pathogenic",
        "examples": ["REVEL score 0.95, CADD score 32, predicted deleterious by "
                     "SIFT and probably damaging by PolyPhen-2"],
    },
    "PP4": {
        "name": "Patient phenotype highly specific for gene",
        "description": "Patient's phenotype or family history is highly specific for a "
                       "disease with a single genetic etiology.",
        "strength": "supporting_pathogenic",
        "examples": ["Classic Marfan phenotype (lens subluxation + aortic dilation + "
                     "systemic score >= 7) in patient with FBN1 variant"],
    },
    "PP5": {
        "name": "Reputable source reports variant as pathogenic",
        "description": "Reputable source recently reports variant as pathogenic, but "
                       "the evidence is not available to the laboratory to perform "
                       "an independent evaluation.",
        "strength": "supporting_pathogenic",
        "examples": ["ClinVar pathogenic assertion from expert panel (ClinGen)"],
    },
    # --- Stand-alone Benign ---
    "BA1": {
        "name": "Allele frequency > 5% in population databases",
        "description": "Allele frequency is >5% in Exome Sequencing Project, "
                       "1000 Genomes Project, or gnomAD.",
        "strength": "stand_alone_benign",
        "examples": ["Variant with 8% allele frequency in gnomAD"],
    },
    # --- Strong Benign ---
    "BS1": {
        "name": "Allele frequency greater than expected for disorder",
        "description": "Allele frequency is greater than expected for disorder "
                       "(using disease-specific thresholds).",
        "strength": "strong_benign",
        "examples": ["Variant at 0.5% frequency for a dominant disorder with "
                     "prevalence of 1:10,000"],
    },
    "BS2": {
        "name": "Observed in healthy adults inconsistent with disease penetrance",
        "description": "Observed in a healthy adult individual for a recessive (homozygous), "
                       "dominant (heterozygous), or X-linked (hemizygous) disorder, "
                       "with full penetrance expected at an early age.",
        "strength": "strong_benign",
        "examples": ["Homozygous variant seen in healthy 60-year-old for an AR "
                     "childhood-onset disorder"],
    },
    "BS3": {
        "name": "Well-established functional studies show no damaging effect",
        "description": "Well-established in vitro or in vivo functional studies "
                       "show no damaging effect on protein function or splicing.",
        "strength": "strong_benign",
        "examples": ["Enzyme assay shows normal activity with the variant protein"],
    },
    "BS4": {
        "name": "Lack of segregation in affected family members",
        "description": "Lack of segregation in affected members of a family "
                       "(with consideration for penetrance).",
        "strength": "strong_benign",
        "examples": ["Variant absent in 3 affected family members with HCM"],
    },
    # --- Supporting Benign ---
    "BP1": {
        "name": "Missense in gene where only truncating cause disease",
        "description": "Missense variant in a gene for which primarily truncating "
                       "variants are known to cause disease.",
        "strength": "supporting_benign",
        "examples": ["Missense in PMS2 where LOF variants are the primary mechanism"],
    },
    "BP2": {
        "name": "Observed in trans/cis with pathogenic variant",
        "description": "Observed in trans with a pathogenic variant for a fully "
                       "penetrant dominant disorder; or observed in cis with a "
                       "pathogenic variant in any inheritance pattern.",
        "strength": "supporting_benign",
        "examples": ["Variant in cis with known pathogenic variant in same allele"],
    },
    "BP3": {
        "name": "In-frame del/ins in repetitive region without known function",
        "description": "In-frame deletions/insertions in a repetitive region without "
                       "a known function.",
        "strength": "supporting_benign",
        "examples": ["In-frame deletion in a polyglutamine tract of unknown function"],
    },
    "BP4": {
        "name": "Multiple computational evidence suggests no impact",
        "description": "Multiple lines of computational evidence suggest no impact on "
                       "gene or gene product (REVEL, CADD, SIFT, PolyPhen, etc.).",
        "strength": "supporting_benign",
        "examples": ["REVEL score 0.05, CADD score 8, predicted tolerated by SIFT, "
                     "benign by PolyPhen-2"],
    },
    "BP5": {
        "name": "Variant found in case with alternate molecular basis",
        "description": "Variant found in a case with an alternate molecular basis "
                       "for disease.",
        "strength": "supporting_benign",
        "examples": ["Second variant found in a gene unrelated to patient's known "
                     "causative variant"],
    },
    "BP6": {
        "name": "Reputable source reports variant as benign",
        "description": "Reputable source recently reports variant as benign, but the "
                       "evidence is not available to the laboratory to perform "
                       "an independent evaluation.",
        "strength": "supporting_benign",
        "examples": ["ClinVar benign assertion from expert panel"],
    },
    "BP7": {
        "name": "Synonymous with no predicted splice impact",
        "description": "A synonymous (silent) variant for which splicing prediction "
                       "algorithms predict no impact to the splice consensus sequence "
                       "nor the creation of a new splice site.",
        "strength": "supporting_benign",
        "examples": ["Synonymous variant in exon middle, no splice predictor flags"],
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 11. HPO TOP-LEVEL TERMS
# ═══════════════════════════════════════════════════════════════════════════════

HPO_TOP_LEVEL_TERMS: Dict[str, Dict[str, str]] = {
    "HP:0000118": {
        "label": "Phenotypic abnormality",
        "description": "Root term for all phenotypic abnormalities in HPO.",
    },
    "HP:0000152": {
        "label": "Abnormality of head or neck",
        "description": "An abnormality of the head or neck region, including cranium, "
                       "face, eyes, ears, nose, mouth, and neck structures.",
    },
    "HP:0000478": {
        "label": "Abnormality of the eye",
        "description": "An abnormality of the eye, including structural, functional, "
                       "and refractive defects.",
    },
    "HP:0000598": {
        "label": "Abnormality of the ear",
        "description": "An abnormality of the ear, including structural, hearing, "
                       "and vestibular abnormalities.",
    },
    "HP:0000707": {
        "label": "Abnormality of the nervous system",
        "description": "An abnormality of the nervous system, including central "
                       "and peripheral nervous system defects.",
    },
    "HP:0000769": {
        "label": "Abnormality of the breast",
        "description": "An abnormality of the breast, including structural and "
                       "developmental defects.",
    },
    "HP:0000818": {
        "label": "Abnormality of the endocrine system",
        "description": "An abnormality of the endocrine system, including pituitary, "
                       "thyroid, parathyroid, adrenal, and gonadal defects.",
    },
    "HP:0000924": {
        "label": "Abnormality of the skeletal system",
        "description": "An abnormality of the skeletal system, including bone density, "
                       "skeletal morphology, and joint abnormalities.",
    },
    "HP:0001197": {
        "label": "Abnormality of prenatal development or birth",
        "description": "An abnormality of prenatal development, including intrauterine "
                       "growth restriction, prematurity, and birth complications.",
    },
    "HP:0001507": {
        "label": "Growth abnormality",
        "description": "An abnormality of growth, including short stature, tall stature, "
                       "and growth patterns.",
    },
    "HP:0001574": {
        "label": "Abnormality of the integument",
        "description": "An abnormality of the skin, hair, nails, or subcutaneous tissue.",
    },
    "HP:0001626": {
        "label": "Abnormality of the cardiovascular system",
        "description": "An abnormality of the cardiovascular system, including "
                       "structural heart defects, cardiomyopathies, and vascular anomalies.",
    },
    "HP:0001871": {
        "label": "Abnormality of blood and blood-forming tissues",
        "description": "An abnormality of the hematopoietic system, including anemias, "
                       "bleeding disorders, and bone marrow failure.",
    },
    "HP:0001939": {
        "label": "Abnormality of metabolism/homeostasis",
        "description": "An abnormality of metabolism, including inborn errors of "
                       "metabolism, electrolyte disturbances, and acid-base disorders.",
    },
    "HP:0002086": {
        "label": "Abnormality of the respiratory system",
        "description": "An abnormality of the respiratory system, including airway, "
                       "lung parenchyma, and pleural defects.",
    },
    "HP:0002664": {
        "label": "Neoplasm",
        "description": "An abnormal tissue growth (benign or malignant) affecting "
                       "any organ system.",
    },
    "HP:0002715": {
        "label": "Abnormality of the immune system",
        "description": "An abnormality of the immune system, including immunodeficiency, "
                       "autoimmunity, and autoinflammatory conditions.",
    },
    "HP:0003011": {
        "label": "Abnormality of the musculature",
        "description": "An abnormality of the muscles, including myopathies, dystrophies, "
                       "and tone abnormalities.",
    },
    "HP:0003549": {
        "label": "Abnormality of connective tissue",
        "description": "An abnormality of connective tissue, including collagen disorders, "
                       "joint hypermobility, and vascular fragility.",
    },
    "HP:0010972": {
        "label": "Abnormality of the genitourinary system",
        "description": "An abnormality of the genitourinary system, including kidney, "
                       "ureter, bladder, and genital abnormalities.",
    },
    "HP:0025031": {
        "label": "Abnormality of the digestive system",
        "description": "An abnormality of the digestive system, including GI tract, "
                       "liver, and pancreatic defects.",
    },
    "HP:0040064": {
        "label": "Abnormality of limbs",
        "description": "An abnormality of the upper or lower limbs, including "
                       "structural, positional, and length discrepancies.",
    },
    "HP:0000119": {
        "label": "Abnormality of the genitourinary system",
        "description": "An abnormality of the genitourinary system, including "
                       "kidney, ureter, bladder, urethra, and reproductive organ "
                       "structural and functional defects.",
    },
}
