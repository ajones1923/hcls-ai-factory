"""Rare Disease Diagnostic Agent — Query Expansion System.

Maps rare-disease terminology to related terms, enabling comprehensive
retrieval across phenotype, gene, disease, and therapy databases.  Each
expansion map covers a domain (diseases, genes, phenotypes, inheritance,
therapies, diagnostics, metabolic pathways) and widens the semantic net
for search queries.

The QueryExpander class orchestrates alias resolution, entity detection,
workflow-aware term boosting, and HPO-style term expansion.

Author: Adam Jones
Date: March 2026
"""

import re
from typing import Dict, List, Optional

from src.models import DiagnosticWorkflowType


# ═══════════════════════════════════════════════════════════════════════
# 1. ENTITY ALIASES (120+ mappings)
# ═══════════════════════════════════════════════════════════════════════

ENTITY_ALIASES: Dict[str, str] = {
    # --- Metabolic disease abbreviations ---
    "PKU": "phenylketonuria",
    "MSUD": "maple syrup urine disease",
    "MCAD": "medium-chain acyl-CoA dehydrogenase deficiency",
    "MCADD": "medium-chain acyl-CoA dehydrogenase deficiency",
    "OTC": "ornithine transcarbamylase deficiency",
    "MMA": "methylmalonic acidemia",
    "PA": "propionic acidemia",
    "MPS": "mucopolysaccharidosis",
    "MPS I": "mucopolysaccharidosis type I",
    "MPS II": "mucopolysaccharidosis type II",
    "MPS III": "mucopolysaccharidosis type III",
    "MPS IV": "mucopolysaccharidosis type IV",
    "MPS VI": "mucopolysaccharidosis type VI",
    "MPS VII": "mucopolysaccharidosis type VII",
    "NPC": "Niemann-Pick disease type C",
    "NPD": "Niemann-Pick disease",
    "GSD": "glycogen storage disease",
    "CDG": "congenital disorders of glycosylation",
    "VLCAD": "very long-chain acyl-CoA dehydrogenase deficiency",
    "IEM": "inborn error of metabolism",
    "LSDs": "lysosomal storage disorders",
    "LSD": "lysosomal storage disease",
    # --- Neurological disease abbreviations ---
    "SMA": "spinal muscular atrophy",
    "DMD": "Duchenne muscular dystrophy",
    "BMD": "Becker muscular dystrophy",
    "CMT": "Charcot-Marie-Tooth disease",
    "TSC": "tuberous sclerosis complex",
    "NF1": "neurofibromatosis type 1",
    "NF2": "neurofibromatosis type 2",
    "HD": "Huntington disease",
    "FRDA": "Friedreich ataxia",
    "AT": "ataxia-telangiectasia",
    "A-T": "ataxia-telangiectasia",
    "RTT": "Rett syndrome",
    "AS": "Angelman syndrome",
    "PWS": "Prader-Willi syndrome",
    "SMEI": "Dravet syndrome",
    "ALS": "amyotrophic lateral sclerosis",
    "HSP": "hereditary spastic paraplegia",
    "LHON": "Leber hereditary optic neuropathy",
    "LCA": "Leber congenital amaurosis",
    "RP": "retinitis pigmentosa",
    # --- Hematologic abbreviations ---
    "SCD": "sickle cell disease",
    "SCA": "sickle cell anemia",
    "HbSS": "sickle cell disease (homozygous)",
    "TDT": "transfusion-dependent thalassemia",
    "DBA": "Diamond-Blackfan anemia",
    "FA": "Fanconi anemia",
    "PNH": "paroxysmal nocturnal hemoglobinuria",
    "VWD": "von Willebrand disease",
    "HS": "hereditary spherocytosis",
    "HH": "hereditary hemochromatosis",
    # --- Connective tissue abbreviations ---
    "EDS": "Ehlers-Danlos syndrome",
    "hEDS": "hypermobile Ehlers-Danlos syndrome",
    "vEDS": "vascular Ehlers-Danlos syndrome",
    "cEDS": "classical Ehlers-Danlos syndrome",
    "OI": "osteogenesis imperfecta",
    "LDS": "Loeys-Dietz syndrome",
    "MFS": "Marfan syndrome",
    # --- Immunologic abbreviations ---
    "SCID": "severe combined immunodeficiency",
    "CGD": "chronic granulomatous disease",
    "CVID": "common variable immunodeficiency",
    "XLA": "X-linked agammaglobulinemia",
    "WAS": "Wiskott-Aldrich syndrome",
    "HIES": "hyper-IgE syndrome",
    "PID": "primary immunodeficiency",
    "PIDD": "primary immunodeficiency disease",
    "IPEX": "immune dysregulation polyendocrinopathy enteropathy X-linked",
    # --- Cancer predisposition abbreviations ---
    "LFS": "Li-Fraumeni syndrome",
    "HNPCC": "Lynch syndrome",
    "FAP": "familial adenomatous polyposis",
    "MEN1": "multiple endocrine neoplasia type 1",
    "MEN2": "multiple endocrine neoplasia type 2",
    "VHL": "von Hippel-Lindau syndrome",
    "RB": "retinoblastoma",
    "HBOC": "hereditary breast and ovarian cancer",
    # --- Cardiac abbreviations ---
    "HCM": "hypertrophic cardiomyopathy",
    "DCM": "dilated cardiomyopathy",
    "ARVC": "arrhythmogenic right ventricular cardiomyopathy",
    "LQTS": "long QT syndrome",
    "BrS": "Brugada syndrome",
    "CPVT": "catecholaminergic polymorphic ventricular tachycardia",
    "FH": "familial hypercholesterolemia",
    # --- Pulmonary abbreviations ---
    "CF": "cystic fibrosis",
    "PCD": "primary ciliary dyskinesia",
    "A1AT": "alpha-1 antitrypsin deficiency",
    "AATD": "alpha-1 antitrypsin deficiency",
    "HHT": "hereditary hemorrhagic telangiectasia",
    # --- Endocrine abbreviations ---
    "CAH": "congenital adrenal hyperplasia",
    # --- Gene names (common) ---
    "SMN1": "survival motor neuron 1",
    "SMN2": "survival motor neuron 2",
    "CFTR": "cystic fibrosis transmembrane conductance regulator",
    "PAH": "phenylalanine hydroxylase",
    "GBA": "glucocerebrosidase",
    "GBA1": "glucocerebrosidase",
    "GLA": "alpha-galactosidase A",
    "GAA": "acid alpha-glucosidase",
    "IDUA": "alpha-L-iduronidase",
    "HEXA": "hexosaminidase A",
    "HBB": "hemoglobin subunit beta",
    "HTT": "huntingtin",
    "FBN1": "fibrillin-1",
    "COL3A1": "collagen type III alpha 1",
    "COL1A1": "collagen type I alpha 1",
    "TP53": "tumor protein p53",
    "BRCA1": "breast cancer 1",
    "BRCA2": "breast cancer 2",
    "RPE65": "retinal pigment epithelium 65",
    "MECP2": "methyl-CpG binding protein 2",
    "DMD": "dystrophin",  # noqa: F601
    "F8": "coagulation factor VIII",
    "F9": "coagulation factor IX",
    # --- Clinical terms ---
    "ERT": "enzyme replacement therapy",
    "SRT": "substrate reduction therapy",
    "HSCT": "hematopoietic stem cell transplant",
    "BMT": "bone marrow transplant",
    "GT": "gene therapy",
    "ASO": "antisense oligonucleotide",
    "AAV": "adeno-associated virus",
    "NBS": "newborn screening",
    "WES": "whole exome sequencing",
    "WGS": "whole genome sequencing",
    "CMA": "chromosomal microarray analysis",
    "NGS": "next-generation sequencing",
    "VUS": "variant of uncertain significance",
    "ACMG": "American College of Medical Genetics and Genomics",
    "HPO": "Human Phenotype Ontology",
    "OMIM": "Online Mendelian Inheritance in Man",
    "LOF": "loss of function",
    "GOF": "gain of function",
    # --- Neuronal ceroid lipofuscinosis aliases ---
    "NCL": "neuronal ceroid lipofuscinosis",
    "CLN2": "CLN2 disease (late infantile neuronal ceroid lipofuscinosis)",
    "CLN3": "CLN3 disease (juvenile neuronal ceroid lipofuscinosis)",
    # --- Additional metabolic abbreviations ---
    "MLD": "metachromatic leukodystrophy",
    # --- Epidermolysis bullosa aliases ---
    "EB": "epidermolysis bullosa",
    "DEB": "dystrophic epidermolysis bullosa",
    # --- Immunodeficiency aliases ---
    "APDS": "activated PI3K delta syndrome",
    "TTP": "thrombotic thrombocytopenic purpura",
    # --- Renal / syndromic aliases ---
    "PKD": "polycystic kidney disease",
    "BBS": "Bardet-Biedl syndrome",
    # --- Glycosylation disorders ---
    "NCF": "Noonan with cardiofaciocutaneous overlap",
    "RASopathy": "RAS-MAPK pathway disorder",
    # --- Approved gene/enzyme therapies (trade names) ---
    "Kuvan": "sapropterin dihydrochloride",
    "Brineura": "cerliponase alfa",
    "Elevidys": "delandistrogene moxeparvovec",
    "Vyjuvek": "beremagene geperpavec",
    "Roctavian": "valoctocogene roxaparvovec",
    "Lyfgenia": "lovotibeglogene autotemcel",
    "XLRP": "X-linked retinitis pigmentosa",
    "AADC": "aromatic L-amino acid decarboxylase deficiency",
}


# ═══════════════════════════════════════════════════════════════════════
# 2. DISEASE SYNONYM MAP (30+ entries)
# ═══════════════════════════════════════════════════════════════════════

DISEASE_SYNONYM_MAP: Dict[str, List[str]] = {
    "phenylketonuria": [
        "PKU", "phenylalanine hydroxylase deficiency", "hyperphenylalaninemia",
        "PAH deficiency", "classic PKU",
    ],
    "spinal_muscular_atrophy": [
        "SMA", "Werdnig-Hoffmann disease", "SMA type 1", "SMA type 2",
        "SMA type 3", "Kugelberg-Welander disease", "SMN1 deletion",
    ],
    "duchenne_muscular_dystrophy": [
        "DMD", "Duchenne", "dystrophinopathy", "DMD-related dystrophy",
        "pseudohypertrophic muscular dystrophy",
    ],
    "cystic_fibrosis": [
        "CF", "mucoviscidosis", "CFTR-related disorder",
        "fibrocystic disease of pancreas",
    ],
    "sickle_cell_disease": [
        "SCD", "sickle cell anemia", "HbSS disease", "drepanocytosis",
        "hemoglobin S disease",
    ],
    "marfan_syndrome": [
        "Marfan", "MFS", "FBN1-related disorder", "fibrillinopathy",
    ],
    "ehlers_danlos_syndrome": [
        "EDS", "Ehlers-Danlos", "cutis hyperelastica",
        "classical EDS", "vascular EDS", "hypermobile EDS",
    ],
    "osteogenesis_imperfecta": [
        "OI", "brittle bone disease", "collagen type I disorder",
        "Lobstein disease", "Vrolik disease",
    ],
    "gaucher_disease": [
        "Gaucher", "glucocerebrosidase deficiency", "GBA1 deficiency",
        "glucosylceramidase deficiency",
    ],
    "fabry_disease": [
        "Fabry", "Anderson-Fabry disease", "alpha-galactosidase A deficiency",
        "GLA deficiency", "angiokeratoma corporis diffusum",
    ],
    "pompe_disease": [
        "Pompe", "glycogen storage disease type II", "GSD II",
        "acid maltase deficiency", "GAA deficiency",
    ],
    "huntington_disease": [
        "Huntington", "HD", "Huntington chorea", "HTT repeat expansion",
    ],
    "rett_syndrome": [
        "Rett", "RTT", "MECP2-related disorder", "MECP2 duplication syndrome",
    ],
    "tay_sachs_disease": [
        "Tay-Sachs", "hexosaminidase A deficiency", "HEXA deficiency",
        "GM2 gangliosidosis",
    ],
    "tuberous_sclerosis": [
        "TSC", "tuberous sclerosis complex", "Bourneville disease",
        "TSC1/TSC2-related disorder",
    ],
    "neurofibromatosis_type_1": [
        "NF1", "von Recklinghausen disease", "neurofibromin deficiency",
    ],
    "hemophilia_a": [
        "hemophilia A", "factor VIII deficiency", "classic hemophilia",
        "F8 deficiency",
    ],
    "hemophilia_b": [
        "hemophilia B", "Christmas disease", "factor IX deficiency",
        "F9 deficiency",
    ],
    "severe_combined_immunodeficiency": [
        "SCID", "bubble boy disease", "combined immunodeficiency",
        "T-B-NK- SCID", "T-B+NK- SCID",
    ],
    "wilson_disease": [
        "Wilson", "hepatolenticular degeneration", "ATP7B deficiency",
        "copper storage disease",
    ],
    "friedreich_ataxia": [
        "FRDA", "Friedreich", "frataxin deficiency", "FXN repeat expansion",
    ],
    "fanconi_anemia": [
        "Fanconi", "FA", "Fanconi pancytopenia", "chromosomal breakage syndrome",
    ],
    "prader_willi_syndrome": [
        "PWS", "Prader-Willi", "15q11 paternal deletion",
    ],
    "angelman_syndrome": [
        "AS", "Angelman", "happy puppet syndrome", "UBE3A deficiency",
    ],
    "niemann_pick_disease": [
        "Niemann-Pick", "NPD", "NPC", "sphingomyelinase deficiency",
        "cholesterol storage disease",
    ],
    "von_hippel_lindau": [
        "VHL", "von Hippel-Lindau", "VHL disease", "hemangioblastomatosis",
    ],
    "lynch_syndrome": [
        "Lynch", "HNPCC", "hereditary nonpolyposis colorectal cancer",
        "mismatch repair deficiency", "dMMR",
    ],
    "li_fraumeni_syndrome": [
        "LFS", "Li-Fraumeni", "TP53 cancer syndrome",
    ],
    "dravet_syndrome": [
        "Dravet", "SMEI", "severe myoclonic epilepsy of infancy",
        "SCN1A-related epilepsy",
    ],
    "krabbe_disease": [
        "Krabbe", "globoid cell leukodystrophy", "galactosylceramidase deficiency",
        "GALC deficiency",
    ],
    "achondroplasia": [
        "FGFR3-related achondroplasia", "short-limbed dwarfism",
    ],
    "retinitis_pigmentosa": [
        "RP", "rod-cone dystrophy", "inherited retinal dystrophy",
        "retinal degeneration",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 3. GENE SYNONYM MAP (30+ entries)
# ═══════════════════════════════════════════════════════════════════════

GENE_SYNONYM_MAP: Dict[str, List[str]] = {
    "SMN1": ["survival motor neuron 1", "SMN", "SMA gene", "5q SMA"],
    "SMN2": ["survival motor neuron 2", "SMN2 copy number", "SMN2 modifier"],
    "DMD": ["dystrophin", "Dp427", "dystrophin gene", "Xp21"],
    "CFTR": ["cystic fibrosis transmembrane conductance regulator", "ABCC7",
             "CF gene", "chloride channel"],
    "PAH": ["phenylalanine hydroxylase", "PKU gene"],
    "GBA1": ["glucocerebrosidase", "GBA", "acid beta-glucosidase", "Gaucher gene"],
    "GLA": ["alpha-galactosidase A", "ceramide trihexosidase", "Fabry gene"],
    "GAA": ["acid alpha-glucosidase", "acid maltase", "Pompe gene"],
    "IDUA": ["alpha-L-iduronidase", "Hurler gene", "MPS I gene"],
    "HEXA": ["hexosaminidase A", "Tay-Sachs gene", "beta-hexosaminidase alpha"],
    "HBB": ["hemoglobin subunit beta", "beta-globin", "sickle gene"],
    "F8": ["coagulation factor VIII", "anti-hemophilic factor", "hemophilia A gene"],
    "F9": ["coagulation factor IX", "Christmas factor", "hemophilia B gene"],
    "FBN1": ["fibrillin-1", "Marfan gene", "fibrillin 1"],
    "COL1A1": ["collagen type I alpha 1 chain", "pro-alpha1(I) collagen", "OI gene"],
    "COL1A2": ["collagen type I alpha 2 chain", "pro-alpha2(I) collagen"],
    "COL3A1": ["collagen type III alpha 1 chain", "vascular EDS gene"],
    "COL5A1": ["collagen type V alpha 1 chain", "classical EDS gene"],
    "HTT": ["huntingtin", "IT15", "Huntington disease gene", "HD gene"],
    "FXN": ["frataxin", "X25", "Friedreich ataxia gene"],
    "MECP2": ["methyl-CpG binding protein 2", "Rett syndrome gene"],
    "SCN1A": ["sodium channel alpha subunit 1", "Nav1.1", "Dravet gene"],
    "NF1": ["neurofibromin", "neurofibromatosis 1 gene"],
    "TSC1": ["hamartin", "tuberous sclerosis 1 gene"],
    "TSC2": ["tuberin", "tuberous sclerosis 2 gene"],
    "TP53": ["tumor protein p53", "p53", "Li-Fraumeni gene"],
    "BRCA1": ["breast cancer 1", "RING finger protein 53"],
    "BRCA2": ["breast cancer 2", "FANCD1", "Fanconi anemia D1"],
    "RPE65": ["retinal pigment epithelium-specific 65 kDa protein",
              "retinoid isomerohydrolase", "Luxturna target"],
    "VWF": ["von Willebrand factor", "factor VIIIR", "VWD gene"],
    "ATP7B": ["ATPase copper transporting beta", "Wilson disease gene",
              "copper-transporting ATPase 2"],
    "IL2RG": ["interleukin 2 receptor subunit gamma", "common gamma chain",
              "X-SCID gene", "gamma-c"],
    "BTK": ["Bruton tyrosine kinase", "agammaglobulinemia tyrosine kinase",
            "XLA gene"],
    "WAS": ["Wiskott-Aldrich syndrome protein", "WASP"],
    "FGFR3": ["fibroblast growth factor receptor 3", "achondroplasia gene"],
}


# ═══════════════════════════════════════════════════════════════════════
# 4. PHENOTYPE MAP (20+ HPO term synonyms)
# ═══════════════════════════════════════════════════════════════════════

PHENOTYPE_MAP: Dict[str, List[str]] = {
    "seizures": [
        "epilepsy", "convulsions", "fits", "epileptic seizure",
        "tonic-clonic seizure", "absence seizure", "myoclonic seizure",
        "focal seizure", "infantile spasms", "West syndrome",
    ],
    "hypotonia": [
        "floppy baby", "low muscle tone", "decreased muscle tone",
        "neonatal hypotonia", "muscular hypotonia", "axial hypotonia",
        "peripheral hypotonia", "central hypotonia",
    ],
    "intellectual_disability": [
        "mental retardation", "cognitive impairment", "learning disability",
        "developmental delay", "global developmental delay",
        "intellectual developmental disorder", "ID",
    ],
    "short_stature": [
        "dwarfism", "growth failure", "growth retardation",
        "short for age", "proportionate short stature",
        "disproportionate short stature",
    ],
    "microcephaly": [
        "small head", "reduced head circumference", "primary microcephaly",
        "secondary microcephaly", "acquired microcephaly",
    ],
    "macrocephaly": [
        "large head", "megalencephaly", "increased head circumference",
    ],
    "hepatomegaly": [
        "enlarged liver", "liver enlargement", "hepatic enlargement",
    ],
    "splenomegaly": [
        "enlarged spleen", "splenic enlargement",
    ],
    "hepatosplenomegaly": [
        "liver and spleen enlargement", "HSM",
    ],
    "ataxia": [
        "cerebellar ataxia", "gait ataxia", "limb ataxia",
        "truncal ataxia", "coordination problems", "unsteady gait",
    ],
    "dystonia": [
        "dystonic posturing", "sustained muscle contraction",
        "involuntary twisting", "generalized dystonia", "focal dystonia",
    ],
    "muscle_weakness": [
        "myopathy", "proximal weakness", "distal weakness",
        "limb-girdle weakness", "progressive weakness",
    ],
    "joint_hypermobility": [
        "hypermobile joints", "joint laxity", "ligamentous laxity",
        "benign joint hypermobility", "generalized hypermobility",
        "Beighton score", "double-jointed",
    ],
    "cardiomyopathy": [
        "heart muscle disease", "hypertrophic cardiomyopathy",
        "dilated cardiomyopathy", "restrictive cardiomyopathy",
    ],
    "aortic_dilation": [
        "aortic root dilation", "aortic aneurysm", "aortic root enlargement",
        "ascending aortic dilation", "aortopathy",
    ],
    "scoliosis": [
        "spinal curvature", "curved spine", "spinal deformity",
    ],
    "hearing_loss": [
        "deafness", "sensorineural hearing loss", "conductive hearing loss",
        "hearing impairment", "prelingual deafness", "SNHL",
    ],
    "visual_impairment": [
        "vision loss", "blindness", "low vision", "visual acuity loss",
        "night blindness", "nyctalopia",
    ],
    "failure_to_thrive": [
        "poor weight gain", "faltering growth", "growth failure",
        "FTT", "poor feeding", "feeding difficulties",
    ],
    "recurrent_infections": [
        "frequent infections", "immunodeficiency", "infection susceptibility",
        "recurrent pneumonia", "recurrent sinusitis", "opportunistic infections",
    ],
    "bleeding": [
        "hemorrhage", "bleeding tendency", "easy bruising",
        "mucocutaneous bleeding", "coagulopathy",
    ],
    "skin_findings": [
        "rash", "dermatitis", "skin lesions", "cafe-au-lait macules",
        "hypopigmented macules", "angiofibromas", "neurofibromas",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 5. INHERITANCE MAP
# ═══════════════════════════════════════════════════════════════════════

INHERITANCE_MAP: Dict[str, List[str]] = {
    "autosomal_dominant": [
        "autosomal dominant", "AD", "dominant inheritance",
        "heterozygous", "50% recurrence risk", "vertical transmission",
        "de novo dominant", "haploinsufficiency",
    ],
    "autosomal_recessive": [
        "autosomal recessive", "AR", "recessive inheritance",
        "homozygous", "compound heterozygous", "25% recurrence risk",
        "carrier", "horizontal transmission", "consanguinity",
    ],
    "x_linked_recessive": [
        "X-linked recessive", "XLR", "X-linked", "hemizygous",
        "carrier female", "affected male", "maternal carrier",
    ],
    "x_linked_dominant": [
        "X-linked dominant", "XLD", "lethal in males",
        "affected females", "de novo X-linked",
    ],
    "mitochondrial": [
        "mitochondrial inheritance", "maternal inheritance",
        "mtDNA mutation", "heteroplasmy", "homoplasmy",
        "mitochondrial DNA", "matrilineal",
    ],
    "imprinting": [
        "genomic imprinting", "imprinting disorder", "uniparental disomy",
        "UPD", "imprinting defect", "parent-of-origin effect",
        "epigenetic", "methylation",
    ],
    "trinucleotide_repeat": [
        "trinucleotide repeat expansion", "repeat expansion",
        "anticipation", "dynamic mutation", "CAG repeat",
        "CGG repeat", "GAA repeat", "CTG repeat",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 6. THERAPY MAP
# ═══════════════════════════════════════════════════════════════════════

THERAPY_MAP: Dict[str, List[str]] = {
    "enzyme_replacement": [
        "enzyme replacement therapy", "ERT", "recombinant enzyme",
        "imiglucerase", "agalsidase", "alglucosidase", "laronidase",
        "idursulfase", "elosulfase", "galsulfase", "vestronidase",
        "olipudase", "infusion therapy",
    ],
    "substrate_reduction": [
        "substrate reduction therapy", "SRT", "miglustat", "eliglustat",
        "substrate inhibition",
    ],
    "gene_therapy": [
        "gene therapy", "gene replacement", "gene transfer",
        "AAV vector", "adeno-associated virus", "lentiviral vector",
        "onasemnogene", "Zolgensma", "voretigene", "Luxturna",
        "etranacogene", "Hemgenix", "valoctocogene", "Roctavian",
        "in vivo gene therapy", "ex vivo gene therapy",
    ],
    "gene_editing": [
        "gene editing", "CRISPR", "CRISPR-Cas9", "base editing",
        "prime editing", "exagamglogene", "Casgevy",
        "genome editing", "zinc finger nuclease", "TALEN",
    ],
    "antisense_oligonucleotide": [
        "antisense oligonucleotide", "ASO", "splice modulation",
        "exon skipping", "nusinersen", "Spinraza", "eteplirsen",
        "golodirsen", "viltolarsen", "casimersen",
    ],
    "small_molecule": [
        "small molecule", "oral therapy", "risdiplam", "Evrysdi",
        "sapropterin", "Kuvan", "migalastat", "Galafold",
        "tafamidis", "ivacaftor", "lumacaftor", "elexacaftor",
        "Trikafta",
    ],
    "stem_cell_transplant": [
        "hematopoietic stem cell transplant", "HSCT", "bone marrow transplant",
        "BMT", "allogeneic transplant", "matched sibling donor",
        "cord blood transplant", "haploidentical transplant",
    ],
    "protein_replacement": [
        "protein replacement", "factor replacement", "factor VIII",
        "factor IX", "immunoglobulin replacement", "IVIG", "SCIG",
        "emicizumab", "Hemlibra",
    ],
    "dietary_management": [
        "dietary therapy", "metabolic diet", "low-protein diet",
        "phenylalanine-restricted diet", "galactose-free diet",
        "leucine-restricted diet", "medical formula",
        "special medical food",
    ],
    "chaperone_therapy": [
        "pharmacological chaperone", "chaperone therapy",
        "migalastat", "ambroxol", "protein folding",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 7. DIAGNOSTIC MAP
# ═══════════════════════════════════════════════════════════════════════

DIAGNOSTIC_MAP: Dict[str, List[str]] = {
    "wes": [
        "whole exome sequencing", "WES", "exome sequencing",
        "clinical exome", "exome", "next-generation sequencing panel",
    ],
    "wgs": [
        "whole genome sequencing", "WGS", "genome sequencing",
        "clinical genome", "short-read sequencing", "long-read sequencing",
    ],
    "cma": [
        "chromosomal microarray analysis", "CMA", "array CGH",
        "SNP array", "copy number analysis", "CNV detection",
        "microarray", "aCGH",
    ],
    "gene_panel": [
        "NGS panel", "gene panel", "targeted panel",
        "multigene panel", "disease-specific panel",
        "comprehensive gene panel",
    ],
    "enzyme_assay": [
        "enzyme assay", "enzyme activity", "leukocyte enzyme assay",
        "dried blood spot enzyme", "fibroblast enzyme assay",
        "plasma enzyme activity",
    ],
    "metabolic_screening": [
        "newborn screening", "NBS", "tandem mass spectrometry",
        "MS/MS", "plasma amino acids", "urine organic acids",
        "acylcarnitine profile", "urine GAGs",
        "urine mucopolysaccharides", "dried blood spot",
    ],
    "karyotype": [
        "karyotype", "chromosome analysis", "G-banding",
        "FISH", "fluorescence in situ hybridization",
    ],
    "single_gene": [
        "single gene testing", "Sanger sequencing",
        "targeted mutation analysis", "deletion/duplication analysis",
        "MLPA", "repeat expansion analysis", "Southern blot",
    ],
    "functional_studies": [
        "functional study", "in vitro assay", "cell-based assay",
        "fibroblast culture", "RNA analysis", "splicing assay",
        "minigene assay", "protein expression study",
    ],
    "imaging": [
        "brain MRI", "skeletal survey", "echocardiography",
        "renal ultrasound", "cardiac MRI", "ophthalmologic exam",
        "fundoscopy", "OCT", "ERG",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 8. METABOLIC MAP
# ═══════════════════════════════════════════════════════════════════════

METABOLIC_MAP: Dict[str, List[str]] = {
    "amino_acid_disorders": [
        "amino acid metabolism", "aminoacidopathy",
        "phenylalanine", "tyrosine", "methionine", "homocysteine",
        "branched-chain amino acids", "leucine", "isoleucine", "valine",
        "maple syrup urine", "phenylketonuria", "homocystinuria",
        "tyrosinemia",
    ],
    "organic_acidemias": [
        "organic acid disorder", "organic acidemia", "organic aciduria",
        "propionic acidemia", "methylmalonic acidemia", "isovaleric acidemia",
        "glutaric aciduria", "3-methylcrotonyl-CoA carboxylase deficiency",
    ],
    "fatty_acid_oxidation": [
        "fatty acid oxidation defect", "FAO defect",
        "MCAD deficiency", "VLCAD deficiency", "LCHAD deficiency",
        "carnitine deficiency", "carnitine palmitoyltransferase",
        "CPT I deficiency", "CPT II deficiency",
        "acylcarnitine", "hypoketotic hypoglycemia",
    ],
    "urea_cycle_defects": [
        "urea cycle disorder", "hyperammonemia",
        "OTC deficiency", "citrullinemia", "argininosuccinic aciduria",
        "arginase deficiency", "CPS1 deficiency", "NAGS deficiency",
        "nitrogen scavenger", "ammonia",
    ],
    "lysosomal_storage": [
        "lysosomal storage disease", "LSD", "lysosomal storage disorder",
        "Gaucher", "Fabry", "Pompe", "Niemann-Pick", "Tay-Sachs",
        "Krabbe", "metachromatic leukodystrophy", "MPS",
        "mucopolysaccharidosis", "mucolipidosis", "sphingolipidosis",
        "glycoproteinosis",
    ],
    "peroxisomal_disorders": [
        "peroxisomal disorder", "peroxisome biogenesis disorder",
        "Zellweger spectrum", "X-linked adrenoleukodystrophy",
        "Refsum disease", "very long-chain fatty acids", "VLCFA",
        "plasmalogen", "phytanic acid",
    ],
    "mitochondrial_disorders": [
        "mitochondrial disease", "respiratory chain defect",
        "MELAS", "MERRF", "Leigh syndrome", "NARP",
        "mitochondrial complex deficiency", "coenzyme Q10 deficiency",
        "pyruvate dehydrogenase deficiency", "lactic acidosis",
        "oxidative phosphorylation",
    ],
    "glycogen_storage": [
        "glycogen storage disease", "GSD", "glycogenosis",
        "von Gierke disease", "GSD type I", "McArdle disease",
        "GSD type V", "Cori disease", "GSD type III",
        "glycogen branching enzyme", "glycogen debranching enzyme",
    ],
    "carbohydrate_disorders": [
        "galactosemia", "fructose intolerance",
        "hereditary fructose intolerance", "galactose-1-phosphate",
        "fructose-1-phosphate", "GALT deficiency",
    ],
    "purine_pyrimidine": [
        "purine metabolism", "pyrimidine metabolism",
        "Lesch-Nyhan syndrome", "HPRT deficiency",
        "adenine phosphoribosyltransferase deficiency",
        "orotic aciduria", "dihydropyrimidine dehydrogenase deficiency",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 9. WORKFLOW TERMS
# ═══════════════════════════════════════════════════════════════════════

_WORKFLOW_TERMS: Dict[DiagnosticWorkflowType, List[str]] = {
    DiagnosticWorkflowType.PHENOTYPE_DRIVEN: [
        "phenotype", "HPO", "clinical features", "signs and symptoms",
        "dysmorphic features", "clinical presentation", "differential diagnosis",
        "phenotype-genotype correlation", "phenotypic overlap",
        "syndromic", "non-syndromic",
    ],
    DiagnosticWorkflowType.WES_WGS_INTERPRETATION: [
        "whole exome sequencing", "whole genome sequencing", "WES", "WGS",
        "variant interpretation", "ACMG classification", "pathogenic variant",
        "likely pathogenic", "variant of uncertain significance", "VUS",
        "loss of function", "missense", "splice site", "copy number variant",
        "structural variant", "trio analysis", "de novo",
    ],
    DiagnosticWorkflowType.METABOLIC_SCREENING: [
        "newborn screening", "tandem mass spectrometry", "metabolic crisis",
        "enzyme assay", "plasma amino acids", "urine organic acids",
        "acylcarnitine profile", "inborn error of metabolism",
        "metabolic acidosis", "hyperammonemia", "hypoglycemia",
        "lactic acidosis",
    ],
    DiagnosticWorkflowType.DYSMORPHOLOGY: [
        "dysmorphic features", "facial gestalt", "congenital anomalies",
        "growth parameters", "physical examination",
        "syndromic diagnosis", "clinical genetics evaluation",
        "Face2Gene", "GestaltMatcher",
    ],
    DiagnosticWorkflowType.NEUROGENETIC: [
        "neurodevelopmental", "seizures", "epilepsy", "ataxia", "dystonia",
        "spasticity", "hypotonia", "regression", "brain MRI",
        "white matter disease", "leukodystrophy", "neuromuscular",
        "EMG", "nerve conduction", "muscle biopsy",
    ],
    DiagnosticWorkflowType.CARDIAC_GENETICS: [
        "cardiomyopathy", "channelopathy", "long QT", "Brugada",
        "sudden cardiac death", "aortic dilation", "echocardiography",
        "cardiac MRI", "cascade screening", "arrhythmia",
        "familial hypercholesterolemia",
    ],
    DiagnosticWorkflowType.CONNECTIVE_TISSUE: [
        "joint hypermobility", "connective tissue disorder", "Marfan",
        "Ehlers-Danlos", "osteogenesis imperfecta", "aortic aneurysm",
        "Beighton score", "collagen", "fibrillin", "skin hyperextensibility",
        "ectopia lentis", "Ghent criteria",
    ],
    DiagnosticWorkflowType.INBORN_ERRORS: [
        "inborn error of metabolism", "IEM", "metabolic disease",
        "enzyme deficiency", "metabolic pathway", "storage disease",
        "lysosomal", "peroxisomal", "mitochondrial",
        "enzyme replacement therapy", "substrate reduction therapy",
    ],
    DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY: [
        "gene therapy", "gene replacement", "AAV vector", "CRISPR",
        "antisense oligonucleotide", "splice modulation",
        "Zolgensma", "Spinraza", "Luxturna", "Hemgenix", "Casgevy",
        "eligibility criteria", "treatment access", "expanded access",
    ],
    DiagnosticWorkflowType.UNDIAGNOSED_DISEASE: [
        "undiagnosed disease", "diagnostic odyssey", "unknown diagnosis",
        "UDN", "Undiagnosed Diseases Network", "unsolved case",
        "phenome analysis", "reanalysis", "research sequencing",
        "novel gene discovery", "functional validation",
        "matchmaker exchange",
    ],
    DiagnosticWorkflowType.GENERAL: [
        "rare disease", "orphan disease", "genetic disease",
        "inherited disorder", "Mendelian disease",
        "clinical genetics", "genetic testing", "genetic counseling",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 10. ENTITY CATEGORIES (for entity detection)
# ═══════════════════════════════════════════════════════════════════════

_ENTITY_CATEGORIES: Dict[str, List[str]] = {
    "disease_categories": [
        "metabolic", "neurological", "hematologic", "connective tissue",
        "immunologic", "cardiac", "cancer predisposition", "endocrine",
        "skeletal", "renal", "pulmonary", "dermatologic", "ophthalmologic",
    ],
    "inheritance_patterns": [
        "autosomal dominant", "autosomal recessive", "X-linked",
        "mitochondrial", "de novo", "imprinting", "trinucleotide repeat",
    ],
    "diagnostic_modalities": [
        "whole exome sequencing", "whole genome sequencing", "gene panel",
        "chromosomal microarray", "karyotype", "FISH", "enzyme assay",
        "newborn screening", "Sanger sequencing", "MLPA",
    ],
    "therapy_types": [
        "enzyme replacement therapy", "gene therapy", "gene editing",
        "CRISPR", "antisense oligonucleotide", "stem cell transplant",
        "substrate reduction", "dietary therapy", "chaperone therapy",
    ],
    "phenotype_domains": [
        "seizures", "hypotonia", "intellectual disability", "ataxia",
        "cardiomyopathy", "hepatosplenomegaly", "skeletal dysplasia",
        "immunodeficiency", "bleeding", "visual impairment",
        "hearing loss", "failure to thrive",
    ],
    "databases": [
        "OMIM", "Orphanet", "ClinVar", "gnomAD", "HPO",
        "GeneReviews", "ClinGen", "HGMD", "DECIPHER", "LOVD",
    ],
    "acmg_classifications": [
        "pathogenic", "likely pathogenic", "VUS",
        "likely benign", "benign",
    ],
    "clinical_settings": [
        "neonatal", "pediatric", "prenatal", "newborn",
        "adult-onset", "late-onset", "early-onset",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 11. NEWBORN SCREENING MAP
# ═══════════════════════════════════════════════════════════════════════

NEWBORN_SCREENING_MAP: Dict[str, List[str]] = {
    "pku_screen": [
        "phenylketonuria", "PKU", "phenylalanine", "PAH",
        "blood phenylalanine level", "Guthrie test",
    ],
    "msud_screen": [
        "maple syrup urine disease", "MSUD", "branched-chain amino acids",
        "leucine", "alloisoleucine", "BCKDHA",
    ],
    "galactosemia_screen": [
        "galactosemia", "GALT", "galactose-1-phosphate",
        "Beutler test", "total galactose",
    ],
    "biotinidase_screen": [
        "biotinidase deficiency", "BTD", "biotinidase activity",
        "biotin", "partial biotinidase deficiency",
    ],
    "cah_screen": [
        "congenital adrenal hyperplasia", "CAH", "17-hydroxyprogesterone",
        "17-OHP", "CYP21A2", "salt-wasting crisis",
    ],
    "cf_screen": [
        "cystic fibrosis", "CF", "immunoreactive trypsinogen", "IRT",
        "CFTR", "sweat chloride",
    ],
    "scd_screen": [
        "sickle cell disease", "SCD", "hemoglobin electrophoresis",
        "HBB", "hemoglobin S", "isoelectric focusing",
    ],
    "scid_screen": [
        "severe combined immunodeficiency", "SCID", "TREC",
        "T-cell receptor excision circles", "lymphopenia",
    ],
    "sma_screen": [
        "spinal muscular atrophy", "SMA", "SMN1",
        "SMN1 deletion", "SMN1 copy number",
    ],
    "mcad_screen": [
        "MCAD deficiency", "MCADD", "ACADM", "C8 acylcarnitine",
        "octanoylcarnitine", "medium-chain acyl-CoA dehydrogenase",
    ],
    "vlcad_screen": [
        "VLCAD deficiency", "ACADVL", "C14:1 acylcarnitine",
        "tetradecenoylcarnitine", "very long-chain acyl-CoA dehydrogenase",
    ],
    "pompe_screen": [
        "Pompe disease", "GAA", "acid alpha-glucosidase activity",
        "dried blood spot GAA", "glycogen storage disease type II",
    ],
    "mps_i_screen": [
        "mucopolysaccharidosis type I", "MPS I", "IDUA",
        "alpha-L-iduronidase activity", "Hurler syndrome",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 12. CLINICAL FEATURES MAP
# ═══════════════════════════════════════════════════════════════════════

CLINICAL_FEATURES_MAP: Dict[str, List[str]] = {
    "storage_features": [
        "coarse facies", "hepatosplenomegaly", "dysostosis multiplex",
        "corneal clouding", "inguinal hernia", "gibbus deformity",
        "vacuolated lymphocytes", "cherry red spot",
    ],
    "neurodegenerative_features": [
        "regression", "seizures", "vision loss", "ataxia",
        "spasticity", "cognitive decline", "myoclonus",
        "white matter abnormalities", "brain atrophy",
    ],
    "connective_tissue_features": [
        "hypermobility", "skin hyperextensibility", "easy bruising",
        "scoliosis", "joint dislocations", "pes planus",
        "aortic root dilation", "mitral valve prolapse",
    ],
    "metabolic_crisis_features": [
        "vomiting", "lethargy", "hyperammonemia", "metabolic acidosis",
        "hypoglycemia", "encephalopathy", "tachypnea",
        "poor feeding", "dehydration",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# QUERY EXPANDER CLASS
# ═══════════════════════════════════════════════════════════════════════

class QueryExpander:
    """Expand rare-disease diagnostic queries with related domain terms.

    Orchestrates 8 domain-specific expansion maps, entity alias
    resolution, phenotype-to-HPO matching, and workflow-aware
    term boosting to maximise recall across rare disease databases.
    """

    def __init__(self) -> None:
        self.expansion_maps: List[Dict[str, List[str]]] = [
            DISEASE_SYNONYM_MAP,
            GENE_SYNONYM_MAP,
            PHENOTYPE_MAP,
            INHERITANCE_MAP,
            THERAPY_MAP,
            DIAGNOSTIC_MAP,
            METABOLIC_MAP,
            NEWBORN_SCREENING_MAP,
            CLINICAL_FEATURES_MAP,
        ]
        self.entity_aliases: Dict[str, str] = ENTITY_ALIASES
        self._comparative_re = re.compile(
            r"\b(?:vs\.?|versus|compared\s+to|compared\s+with|"
            r"differ(?:ent|ence)|distinguish)\b",
            re.IGNORECASE,
        )

    # ───────────────────────────────────────────────────────────────
    # Public API
    # ───────────────────────────────────────────────────────────────

    def expand(
        self,
        query: str,
        workflow: Optional[DiagnosticWorkflowType] = None,
    ) -> dict:
        """Expand a query with related rare-disease diagnostic terms.

        Parameters
        ----------
        query : str
            The raw user query.
        workflow : DiagnosticWorkflowType, optional
            If provided, additional workflow-specific terms are appended.

        Returns
        -------
        dict
            ``original``            - the raw query text
            ``expanded_terms``      - list of additional search terms
            ``detected_entities``   - dict of categorised entities found
            ``is_comparative``      - whether this is a comparison query
            ``workflow_hint``       - workflow type used for expansion
        """
        query_lower = query.lower().strip()

        # Resolve aliases first — inject canonical forms
        resolved_query = query_lower
        for abbr, canonical in self.entity_aliases.items():
            pattern = re.compile(r"\b" + re.escape(abbr) + r"\b")
            if pattern.search(resolved_query):
                resolved_query = resolved_query + " " + canonical.lower()

        # Collect expanded terms from all maps
        expanded_terms: List[str] = []
        for exp_map in self.expansion_maps:
            for trigger, terms in exp_map.items():
                if self._trigger_matches(trigger, resolved_query):
                    for term in terms:
                        if (
                            term.lower() not in query_lower
                            and term not in expanded_terms
                        ):
                            expanded_terms.append(term)

        # Add workflow-specific terms
        if workflow is not None:
            workflow_terms = self.get_workflow_terms(workflow)
            for term in workflow_terms:
                if (
                    term.lower() not in query_lower
                    and term not in expanded_terms
                ):
                    expanded_terms.append(term)

        # Detect entities
        detected_entities = self.detect_entities(query)

        # Detect comparative queries
        is_comparative = bool(self._comparative_re.search(query))

        return {
            "original": query,
            "expanded_terms": expanded_terms,
            "detected_entities": detected_entities,
            "is_comparative": is_comparative,
            "workflow_hint": workflow.value if workflow else None,
        }

    def detect_entities(self, text: str) -> dict:
        """Detect rare-disease entities in text by category.

        Parameters
        ----------
        text : str
            Input text to scan.

        Returns
        -------
        dict
            Keys are category names; values are lists of matched
            entity strings.
        """
        detected: Dict[str, List[str]] = {}

        for category, entities in _ENTITY_CATEGORIES.items():
            matches: List[str] = []
            for entity in entities:
                pattern = re.compile(
                    r"\b" + re.escape(entity) + r"\b", re.IGNORECASE
                )
                if pattern.search(text):
                    matches.append(entity)
            if matches:
                detected[category] = matches

        # Also detect aliases that appear directly
        alias_matches: List[str] = []
        for abbr in self.entity_aliases:
            pattern = re.compile(r"\b" + re.escape(abbr) + r"\b")
            if pattern.search(text):
                canonical = self.entity_aliases[abbr]
                entry = f"{abbr} ({canonical})"
                if entry not in alias_matches:
                    alias_matches.append(entry)

        if alias_matches:
            detected["resolved_aliases"] = alias_matches

        return detected

    def get_workflow_terms(self, workflow: DiagnosticWorkflowType) -> List[str]:
        """Return additional search terms for a specific workflow.

        Parameters
        ----------
        workflow : DiagnosticWorkflowType
            The diagnostic workflow.

        Returns
        -------
        list of str
            Terms relevant to the workflow.
        """
        return list(_WORKFLOW_TERMS.get(workflow, []))

    def expand_hpo_terms(self, term: str) -> dict:
        """Expand a term with broader, narrower, and related HPO-style terms.

        Parameters
        ----------
        term : str
            The term to expand.

        Returns
        -------
        dict
            ``broader``  - broader/parent terms
            ``narrower`` - narrower/child terms
            ``related``  - related terms at same level
        """
        term_lower = term.lower().strip()

        broader: List[str] = []
        narrower: List[str] = []
        related: List[str] = []

        # Check phenotype maps for hierarchy
        for phenotype, terms in PHENOTYPE_MAP.items():
            terms_lower = [t.lower() for t in terms]
            if term_lower == phenotype.replace("_", " "):
                narrower.extend(terms)
            elif term_lower in terms_lower:
                broader.append(phenotype.replace("_", " "))
                related.extend(
                    t for t in terms
                    if t.lower() != term_lower
                )

        # Check disease synonym maps
        for disease, synonyms in DISEASE_SYNONYM_MAP.items():
            synonyms_lower = [s.lower() for s in synonyms]
            if term_lower == disease.replace("_", " "):
                related.extend(synonyms)
            elif term_lower in synonyms_lower:
                broader.append(disease.replace("_", " "))
                related.extend(
                    s for s in synonyms
                    if s.lower() != term_lower
                )

        # Check gene synonym maps
        for gene, synonyms in GENE_SYNONYM_MAP.items():
            synonyms_lower = [s.lower() for s in synonyms]
            if term_lower == gene.lower():
                narrower.extend(synonyms)
                broader.append("gene")
            elif term_lower in synonyms_lower:
                broader.append(gene)
                related.extend(
                    s for s in synonyms
                    if s.lower() != term_lower
                )

        # Check metabolic maps
        for pathway, terms in METABOLIC_MAP.items():
            terms_lower = [t.lower() for t in terms]
            if term_lower == pathway.replace("_", " "):
                narrower.extend(terms)
                broader.append("metabolic disorder")
            elif term_lower in terms_lower:
                broader.append(pathway.replace("_", " "))
                related.extend(
                    t for t in terms
                    if t.lower() != term_lower
                )

        # Deduplicate while preserving order
        broader = list(dict.fromkeys(broader))
        narrower = list(dict.fromkeys(narrower))
        related = list(dict.fromkeys(related))

        return {
            "broader": broader,
            "narrower": narrower,
            "related": related,
        }

    # ───────────────────────────────────────────────────────────────
    # Internal helpers
    # ───────────────────────────────────────────────────────────────

    @staticmethod
    def _trigger_matches(trigger: str, query_lower: str) -> bool:
        """Check if a trigger term appears as a whole word in the query."""
        # Replace underscores with spaces for matching
        trigger_text = trigger.replace("_", " ").lower()
        pattern = re.compile(r"\b" + re.escape(trigger_text) + r"\b")
        return bool(pattern.search(query_lower))
