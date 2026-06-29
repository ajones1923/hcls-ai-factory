"""Tumor microenvironment (TME) parser for the Single-Cell Intelligence Agent.

Seeds TME profiles for major cancer types, capturing immune composition,
stromal components, spatial organization patterns, and therapeutic
implications from single-cell atlas studies.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA: TME PROFILES FOR MAJOR CANCER TYPES
# ===================================================================

TME_PROFILES: List[Dict[str, Any]] = [
    {
        "cancer_type": "Non-small cell lung cancer (NSCLC)",
        "cancer_code": "NSCLC",
        "immune_composition": {
            "T_cells": "CD8+ TILs often exhausted (PD-1+TIM-3+LAG-3+); Tregs enriched in tumor core",
            "macrophages": "M2-TAMs (CD163+CD206+) predominate; M1 macrophages in stroma",
            "dendritic_cells": "cDC1 scarcity limits cross-presentation; pDC produce type I IFN",
            "b_cells": "Tertiary lymphoid structures (TLS) with germinal centers in responders",
            "nk_cells": "CD56dim NK cells with reduced cytotoxicity in TME",
        },
        "stromal_components": ["Cancer-associated fibroblasts (iCAF, myCAF, apCAF)", "Tumor endothelium", "Pericytes"],
        "spatial_patterns": "Immune exclusion common; T cells restricted to invasive margin. TLS formation correlates with immunotherapy response.",
        "therapeutic_implications": "PD-1/PD-L1 checkpoint blockade standard of care. TIL exhaustion markers predict response. STK11/KRAS mutations associate with cold TME.",
        "key_references": ["Zilionis et al. Immunity 2019", "Wu et al. Nat Med 2021"],
    },
    {
        "cancer_type": "Breast cancer",
        "cancer_code": "BRCA",
        "immune_composition": {
            "T_cells": "TNBC: high TIL infiltration; HR+: lower immune infiltration. CD8+ T cell clonal expansion in responders",
            "macrophages": "TAMs promote angiogenesis and invasion. CD163+ M2 TAMs associate with poor prognosis",
            "dendritic_cells": "LAMP3+ mature DCs in TLS. DC-SIGN+ DCs in stroma",
            "b_cells": "Plasma cell enrichment in TNBC with IgG class switching",
            "nk_cells": "NK cell dysfunction via TGF-beta and PGE2 in TME",
        },
        "stromal_components": ["CAF subtypes (inflammatory, contractile, antigen-presenting)", "Adipocytes", "Endothelial cells"],
        "spatial_patterns": "Stromal TILs predict pCR in neoadjuvant chemotherapy. Immune hot/cold heterogeneity within same tumor.",
        "therapeutic_implications": "TNBC: pembrolizumab + chemo in PD-L1+ CPS>=10. HR+: lower immunogenicity. HER2+: trastuzumab-deruxtecan (T-DXd).",
        "key_references": ["Wu et al. Cell 2021", "Bassez et al. Nat Med 2021"],
    },
    {
        "cancer_type": "Colorectal cancer (CRC)",
        "cancer_code": "CRC",
        "immune_composition": {
            "T_cells": "MSI-H: dense CD8+ TIL infiltrate; MSS: immune desert. Th17 cells in CMS1",
            "macrophages": "CMS4 mesenchymal: high TAM density. SPP1+ macrophages promote fibrosis",
            "dendritic_cells": "DC maturation impaired by Wnt signaling",
            "b_cells": "Plasma cell infiltration in MSI-H tumors",
            "nk_cells": "NK dysfunction via NKG2D ligand shedding",
        },
        "stromal_components": ["Myofibroblastic CAFs", "Pericytes", "Cancer stem cells (Lgr5+)"],
        "spatial_patterns": "Immunoscore (CD3/CD8 at tumor center and invasive margin) predicts recurrence better than TNM staging.",
        "therapeutic_implications": "MSI-H/dMMR: pembrolizumab first-line. MSS: limited checkpoint benefit. Immunoscore guides adjuvant therapy.",
        "key_references": ["Zhang et al. Nature 2020", "Pelka et al. Cell 2021"],
    },
    {
        "cancer_type": "Melanoma",
        "cancer_code": "SKCM",
        "immune_composition": {
            "T_cells": "High mutational burden drives neoantigen-specific CD8+ T cells. TCF7+ progenitor exhausted T cells replenish effectors",
            "macrophages": "Melanoma-associated macrophages express PD-L1. TREM2+ macrophages immunosuppressive",
            "dendritic_cells": "cDC1 essential for anti-PD-1 response. Batf3-dependent DC priming",
            "b_cells": "TLS with B cell germinal centers in responders to immunotherapy",
            "nk_cells": "NK cells contribute to ADCC with anti-CTLA-4 antibodies",
        },
        "stromal_components": ["Tumor-associated fibroblasts", "Melanocyte stem cells", "Endothelial cells"],
        "spatial_patterns": "Brisk vs non-brisk TIL patterns. Spatial proximity of CD8 T cells to tumor cells predicts response.",
        "therapeutic_implications": "Ipilimumab + nivolumab combination. LAG-3 blockade (relatlimab). Adoptive TIL therapy showing 50%+ ORR.",
        "key_references": ["Sade-Feldman et al. Cell 2018", "Jerby-Arnon et al. Cell 2018"],
    },
    {
        "cancer_type": "Pancreatic ductal adenocarcinoma (PDAC)",
        "cancer_code": "PAAD",
        "immune_composition": {
            "T_cells": "Severely immune-excluded; CD8+ T cells trapped in stroma. Tregs predominate",
            "macrophages": "Dense M2-TAM infiltration. CSF1R+ macrophages maintain immunosuppression",
            "dendritic_cells": "DC paucity in tumor. Impaired antigen presentation",
            "b_cells": "IL-35-producing Breg cells contribute to immunosuppression",
            "nk_cells": "NK cells rare in PDAC TME",
        },
        "stromal_components": ["Dense desmoplastic stroma", "Stellate cells/CAFs", "Hyaluronan-rich ECM"],
        "spatial_patterns": "Immune-excluded phenotype dominant. Dense stroma acts as physical barrier to immune infiltration.",
        "therapeutic_implications": "Checkpoint monotherapy largely ineffective. Stromal targeting (anti-CD40, CXCR4 inhibitors) may reprogram TME. KRAS G12C inhibitors emerging.",
        "key_references": ["Peng et al. Cell Res 2019", "Steele et al. Nat Med 2020"],
    },
    {
        "cancer_type": "Glioblastoma (GBM)",
        "cancer_code": "GBM",
        "immune_composition": {
            "T_cells": "T cell exhaustion extreme. TIM-3 dominant exhaustion checkpoint. Few neoantigen-specific clones",
            "macrophages": "Microglia and bone marrow-derived macrophages (BMDM). MHCIIhi microglia vs MHCIIlo BMDM",
            "dendritic_cells": "Scarcity of conventional DCs in brain parenchyma",
            "b_cells": "Rare B cell infiltration",
            "nk_cells": "NK cells comprise <5% of immune infiltrate",
        },
        "stromal_components": ["Glioma stem cells (GSC)", "Astrocytes", "Pericytes", "Oligodendrocyte precursors"],
        "spatial_patterns": "Blood-brain barrier limits immune trafficking. Perivascular immune niches around blood vessels.",
        "therapeutic_implications": "Temozolomide standard. Anti-PD-1 monotherapy failed in CheckMate 143. CAR-T targeting EGFRvIII under investigation.",
        "key_references": ["Pombo Antunes et al. Nat Neurosci 2021", "Abdelfattah et al. Cell 2022"],
    },
    {
        "cancer_type": "Hepatocellular carcinoma (HCC)",
        "cancer_code": "LIHC",
        "immune_composition": {
            "T_cells": "LAYN+ exhausted CD8+ T cells. MAIT cells enriched in liver. Treg accumulation via CCL17/CCL22",
            "macrophages": "Kupffer cells and recruited monocyte-derived macrophages. TREM2+ macrophages",
            "dendritic_cells": "LAMP3+ DCs migrate to lymph nodes for T cell priming",
            "b_cells": "IgA+ plasma cells contribute to immunosuppression",
            "nk_cells": "Liver NK cells (CD56bright) abundant but functionally impaired",
        },
        "stromal_components": ["Hepatic stellate cells/CAFs", "Sinusoidal endothelium", "Cholangiocytes"],
        "spatial_patterns": "Immune-rich border vs immune-poor core. Portal tract immune aggregates.",
        "therapeutic_implications": "Atezolizumab + bevacizumab first-line. Tremelimumab + durvalumab alternative. HBV viral status influences TME.",
        "key_references": ["Zhang et al. Cell 2019", "Sun et al. Cell 2021"],
    },
    {
        "cancer_type": "Renal cell carcinoma (RCC)",
        "cancer_code": "KIRC",
        "immune_composition": {
            "T_cells": "High T cell infiltration but exhaustion. CD8+CXCL13+ T cells mark tumor-reactive clones",
            "macrophages": "M2 TAMs enriched. CD163+ macrophages correlate with grade",
            "dendritic_cells": "cDC1 presence correlates with response to nivolumab",
            "b_cells": "Tertiary lymphoid structures predict immunotherapy benefit",
            "nk_cells": "NK cell frequency reduced compared to normal kidney",
        },
        "stromal_components": ["CAFs with myofibroblast features", "Clear cell lipid-laden tumor cells", "Endothelial cells (VEGF-driven angiogenesis)"],
        "spatial_patterns": "Immune-inflamed subtype common. Angiogenic subtype has vessel-rich, immune-poor regions.",
        "therapeutic_implications": "Ipilimumab + nivolumab or pembrolizumab + axitinib first-line. VEGF-TKI + IO combinations. PBRM1 loss associates with IO benefit.",
        "key_references": ["Bi et al. Cancer Cell 2021", "Krishna et al. Science 2021"],
    },
    {
        "cancer_type": "Ovarian cancer",
        "cancer_code": "OV",
        "immune_composition": {
            "T_cells": "Intraepithelial CD8+ TILs predict survival. CD103+CD8+ tissue-resident memory T cells",
            "macrophages": "TAMs dominate immune infiltrate. Folate receptor beta (FR-beta)+ TAMs",
            "dendritic_cells": "cDC1 scarcity limits anti-tumor immunity",
            "b_cells": "B cells in TLS correlate with chemo-response",
            "nk_cells": "Ascites NK cells show impaired cytotoxicity",
        },
        "stromal_components": ["Mesothelial cells", "CAFs", "Ascites-associated cells"],
        "spatial_patterns": "Immune hot (intraepithelial TILs) vs immune cold (stromal restricted). Peritoneal metastases with unique TME.",
        "therapeutic_implications": "PARP inhibitors in BRCA1/2-mutated. Niraparib + pembrolizumab under investigation. VEGF blockade (bevacizumab) in maintenance.",
        "key_references": ["Vázquez-García et al. Nature 2022", "Hornburg et al. Cancer Cell 2021"],
    },
    {
        "cancer_type": "Head and neck squamous cell carcinoma (HNSCC)",
        "cancer_code": "HNSC",
        "immune_composition": {
            "T_cells": "HPV+: higher CD8+ TIL density. HPV-: more immunosuppressive. Exhausted T cells express PD-1, CTLA-4",
            "macrophages": "M2 TAMs enriched in tumor nests. CD68+ density varies by HPV status",
            "dendritic_cells": "Mature DCs in peritumoral regions of HPV+ tumors",
            "b_cells": "Plasma cells and TLS more common in HPV+ tumors",
            "nk_cells": "NKG2D downregulation in TME",
        },
        "stromal_components": ["Myofibroblastic CAFs", "Endothelial cells", "Salivary gland progenitors"],
        "spatial_patterns": "HPV status determines immune landscape. HPV+: inflamed. HPV-: variable exclusion patterns.",
        "therapeutic_implications": "Pembrolizumab first-line in PD-L1 CPS>=1. HPV+ tumors more immunogenic. De-escalation trials in HPV+ disease.",
        "key_references": ["Cillo et al. Immunity 2020", "Puram et al. Cell 2017"],
    },
    {
        "cancer_type": "Bladder cancer (urothelial carcinoma)",
        "cancer_code": "BLCA",
        "immune_composition": {
            "T_cells": "Luminal-infiltrated subtype has CD8+ TILs. Basal-squamous subtype more immune-excluded",
            "macrophages": "SPP1+ macrophages in muscle-invasive disease. M2-TAMs correlate with invasion",
            "dendritic_cells": "cDC1 density predicts atezolizumab response",
            "b_cells": "B cells and plasma cells in stroma",
            "nk_cells": "NK cell infiltration variable",
        },
        "stromal_components": ["Smooth muscle cells", "CAFs", "Urothelial progenitors"],
        "spatial_patterns": "Molecular subtypes (luminal, basal, neuronal-like) determine immune topology.",
        "therapeutic_implications": "Avelumab maintenance after platinum chemo. Enfortumab vedotin + pembrolizumab emerging. FGFR3 mutations in luminal subtypes.",
        "key_references": ["Chen et al. Cell 2020", "Luo et al. Cancer Cell 2022"],
    },
    {
        "cancer_type": "Prostate cancer",
        "cancer_code": "PRAD",
        "immune_composition": {
            "T_cells": "Low mutational burden limits neoantigen-driven immunity. Tregs prominent. CD8+ T cells sparse",
            "macrophages": "M2 TAMs predominate. TREM2+ macrophages immunosuppressive",
            "dendritic_cells": "DC infiltration low compared to other solid tumors",
            "b_cells": "B cell infiltration limited",
            "nk_cells": "NK cell frequency very low in prostate TME",
        },
        "stromal_components": ["Luminal and basal epithelium", "Neuroendocrine cells", "Smooth muscle stroma"],
        "spatial_patterns": "Immunologically cold tumor. Myeloid-dominated immune landscape.",
        "therapeutic_implications": "Ipilimumab + nivolumab in mCRPC with CDK12/MMR deficiency. Sipuleucel-T DC vaccine. PARP inhibitors in BRCA2.",
        "key_references": ["Chen et al. Cancer Cell 2021", "He et al. Cell 2021"],
    },
    {
        "cancer_type": "Acute myeloid leukemia (AML)",
        "cancer_code": "AML",
        "immune_composition": {
            "T_cells": "Exhausted CD8+ T cells with high PD-1/TIM-3. Tregs enriched in bone marrow. Reduced TCR diversity",
            "macrophages": "Leukemia-associated macrophages promote LSC survival. CD163+ macrophages in marrow niches",
            "dendritic_cells": "DC maturation impaired by AML blasts. Reduced antigen presentation capacity",
            "b_cells": "Normal B cell precursors suppressed by leukemic blasts",
            "nk_cells": "NK cell dysfunction via NKG2D ligand shedding by AML blasts",
        },
        "stromal_components": ["Bone marrow mesenchymal stromal cells", "Leukemic stem cell (LSC) niches", "Endosteal and perivascular niches", "Adipocytes"],
        "spatial_patterns": "LSCs reside in protective endosteal niches. Perivascular niches support blast survival. Normal hematopoiesis suppressed by TME remodeling.",
        "therapeutic_implications": "Venetoclax + azacitidine disrupts LSC metabolism. CD47 blockade (magrolimab) enhances phagocytosis. Flotetuzumab bispecific targets CD123. Midostaurin for FLT3-mutated AML.",
        "key_references": ["van Galen et al. Cell 2019", "Lasry et al. Nat Cancer 2023"],
    },
    {
        "cancer_type": "Chronic lymphocytic leukemia (CLL)",
        "cancer_code": "CLL",
        "immune_composition": {
            "T_cells": "CD8+ T cell exhaustion with PD-1 upregulation. Expanded Th2 and Treg populations. Impaired immune synapse formation",
            "macrophages": "Nurse-like cells (NLCs) are CLL-specific M2-like macrophages that protect CLL cells from apoptosis via BAFF and APRIL",
            "dendritic_cells": "DC dysfunction with impaired antigen presentation. CLL cells can differentiate into DC-like cells",
            "b_cells": "Malignant B cells dominate. Residual normal B cell immunity severely impaired",
            "nk_cells": "NK cell cytotoxicity reduced. Downregulation of activating receptors NKp30 and NKp44",
        },
        "stromal_components": ["Nurse-like cells", "Follicular dendritic cells", "Bone marrow stromal cells", "Lymph node proliferation centers"],
        "spatial_patterns": "Lymph node proliferation centers are sites of active CLL cell proliferation with T cell and stromal cell support. Bone marrow shows diffuse or nodular infiltration patterns.",
        "therapeutic_implications": "BTK inhibitors (ibrutinib, acalabrutinib) disrupt BCR signaling and TME interactions. Venetoclax targets BCL-2. Anti-CD20 (obinutuzumab) depletes CLL cells. CAR-T targeting CD19 in relapsed/refractory.",
        "key_references": ["Rendeiro et al. Nat Commun 2020", "Haselager et al. Blood 2023"],
    },
    {
        "cancer_type": "Multiple myeloma",
        "cancer_code": "MM",
        "immune_composition": {
            "T_cells": "T cell exhaustion and senescence. Reduced CD4:CD8 ratio. Treg expansion in bone marrow",
            "macrophages": "M2-polarized marrow macrophages support myeloma cell survival. CD163+ TAMs promote drug resistance",
            "dendritic_cells": "DC function impaired. Plasmacytoid DCs produce IL-6 that sustains myeloma growth",
            "b_cells": "Normal B cell and plasma cell populations suppressed. Immunoparesis with reduced uninvolved immunoglobulins",
            "nk_cells": "NK cell dysfunction but retained ADCC potential exploited by elotuzumab (anti-SLAMF7)",
        },
        "stromal_components": ["Bone marrow stromal cells", "Osteoclasts (activated by RANKL)", "Osteoblasts (suppressed)", "Bone marrow adipocytes", "Myeloma-associated fibroblasts"],
        "spatial_patterns": "Bone marrow plasma cell niche with osteoclast activation and osteolytic lesions. Extramedullary disease indicates TME-independent growth.",
        "therapeutic_implications": "Daratumumab (anti-CD38) + lenalidomide backbone. BCMA-targeting CAR-T (ide-cel, cilta-cel) and bispecifics (teclistamab). IMiDs modulate TME immunity. Elotuzumab enhances NK cell ADCC.",
        "key_references": ["Cohen et al. Nat Med 2021", "de Jong et al. Blood Cancer Discov 2022"],
    },
    {
        "cancer_type": "Thyroid cancer",
        "cancer_code": "THCA",
        "immune_composition": {
            "T_cells": "Immune-rich but functionally impaired. BRAF V600E-mutated tumors have higher CD8+ TIL density but exhausted phenotype",
            "macrophages": "Dense TAM infiltration correlates with dedifferentiation. M2 TAMs at invasive front",
            "dendritic_cells": "LAMP3+ mature DCs present but antigen presentation capacity reduced by tumor-derived TGF-beta",
            "b_cells": "B cell infiltrates associated with autoimmune thyroiditis background (Hashimoto's). TLS in PTC",
            "nk_cells": "NK cells present but with reduced NKG2D expression and cytotoxic function",
        },
        "stromal_components": ["CAFs with thyroid stroma features", "Endothelial cells", "Thyroid follicular epithelium"],
        "spatial_patterns": "Immune-rich but often functionally excluded. BRAF-mutated PTC has inflamed phenotype but immune escape via PD-L1 upregulation. Anaplastic thyroid cancer shifts to immune desert.",
        "therapeutic_implications": "RAI-refractory DTC: lenvatinib/sorafenib. BRAF V600E: dabrafenib + trametinib. Checkpoint inhibitors in ATC showing activity with pembrolizumab. RET inhibitors (selpercatinib) for RET-altered tumors.",
        "key_references": ["Lu et al. Cell Rep 2023", "Luo et al. Thyroid 2022"],
    },
    {
        "cancer_type": "Endometrial cancer",
        "cancer_code": "UCEC",
        "immune_composition": {
            "T_cells": "MSI-H/POLE ultramutated: dense CD8+ TIL infiltration ('hot'). Copy-number low: moderate immune infiltration. Copy-number high (serous-like): immune cold",
            "macrophages": "M2 TAMs enriched in serous-like subtype. CD68+ macrophage density varies by molecular subtype",
            "dendritic_cells": "cDC1 density higher in MSI-H subtype, correlating with neoantigen load",
            "b_cells": "Plasma cell infiltration in MSI-H tumors. TLS formation in immune-hot subtypes",
            "nk_cells": "NK cell activity variable across subtypes",
        },
        "stromal_components": ["Endometrial stromal cells", "CAFs", "Myometrial smooth muscle"],
        "spatial_patterns": "TCGA molecular subtypes define immune landscape. MSI-H and POLE subtypes are immunologically hot with high TMB. p53-abnormal serous-like subtype often immune-excluded.",
        "therapeutic_implications": "Dostarlimab-gxly for dMMR/MSI-H (first tissue-agnostic approval applies). Pembrolizumab + lenvatinib for MSS/pMMR advanced disease. Trastuzumab-deruxtecan for HER2+ subset.",
        "key_references": ["Sanchez-Vega et al. Cancer Cell 2023", "Makker et al. NEJM 2022"],
    },
    {
        "cancer_type": "Mesothelioma",
        "cancer_code": "MESO",
        "immune_composition": {
            "T_cells": "CD8+ T cells present but deeply exhausted (PD-1+LAG-3+TIM-3+). Tregs enriched. Low neoantigen burden despite asbestos exposure",
            "macrophages": "Highly immunosuppressive TAM population. M2-polarized macrophages dominate. Asbestos-driven chronic inflammation recruits monocytes",
            "dendritic_cells": "DC maturation severely impaired by tumor-derived IL-10 and TGF-beta",
            "b_cells": "B cell infiltrates present but limited germinal center activity",
            "nk_cells": "NK cell function suppressed by soluble NKG2D ligands and PGE2",
        },
        "stromal_components": ["Dense fibrous stroma", "Mesothelial-derived CAFs", "Asbestos fibers embedded in tissue", "Hyaluronan-rich ECM"],
        "spatial_patterns": "Highly immunosuppressive TME with asbestos-driven chronic inflammation. Dense fibrotic stroma limits immune cell access. Epithelioid subtype more immune-infiltrated than sarcomatoid.",
        "therapeutic_implications": "Ipilimumab + nivolumab first-line (CheckMate 743). BAP1 loss associates with distinct TME and immunotherapy response. Mesothelin-targeted therapies and anti-CTLA-4/PD-1 combinations under investigation.",
        "key_references": ["Bueno et al. Nat Genet 2016", "Alcala et al. Cancer Discov 2022"],
    },
    {
        "cancer_type": "Cholangiocarcinoma (bile duct cancer)",
        "cancer_code": "CHOL",
        "immune_composition": {
            "T_cells": "CD8+ T cells present but restricted to peritumoral stroma. Tregs accumulate at tumor-stroma interface",
            "macrophages": "Dense TAM infiltration with M2 polarization. SPP1+ macrophages promote desmoplastic reaction",
            "dendritic_cells": "DC paucity within tumor parenchyma. Impaired cross-presentation of tumor antigens",
            "b_cells": "B cell infiltrates in peribiliary regions. Tertiary lymphoid structures rare",
            "nk_cells": "NK cell infiltration very low",
        },
        "stromal_components": ["Dense desmoplastic stroma", "CAF-rich reactive fibrosis", "Hepatic stellate cell-derived CAFs", "Biliary epithelium"],
        "spatial_patterns": "Extremely dense desmoplastic stroma with CAF-rich microenvironment. Immune cells largely excluded from tumor core. Peribiliary immune aggregates at periphery.",
        "therapeutic_implications": "Gemcitabine + cisplatin + durvalumab first-line (TOPAZ-1). FGFR2 fusions: pemigatinib/futibatinib. IDH1 mutations: ivosidenib. HER2 amplification: trastuzumab-deruxtecan being evaluated.",
        "key_references": ["Job et al. Hepatology 2020", "Loeuillard et al. J Hepatol 2023"],
    },
    {
        "cancer_type": "Soft tissue sarcoma",
        "cancer_code": "SARC",
        "immune_composition": {
            "T_cells": "Highly heterogeneous by histological subtype. UPS: dense CD8+ TIL infiltration. Well-differentiated liposarcoma: immune sparse. GIST: variable",
            "macrophages": "TAMs often most abundant immune population. M2 macrophages correlate with poor outcome in leiomyosarcoma",
            "dendritic_cells": "DC infiltration generally low. Better DC presence in UPS compared to other subtypes",
            "b_cells": "B cell infiltration variable. TLS present in some UPS and dedifferentiated liposarcoma",
            "nk_cells": "NK cell infiltration low across most subtypes",
        },
        "stromal_components": ["Tumor-associated vasculature", "Heterogeneous mesenchymal stroma", "Tumor-specific ECM (myxoid, lipomatous, fibrous)"],
        "spatial_patterns": "TME varies dramatically by subtype. Undifferentiated pleomorphic sarcoma (UPS) is most immune-infiltrated. Synovial sarcoma and Ewing sarcoma are typically immune cold. Well-differentiated liposarcoma has adipocyte-rich TME.",
        "therapeutic_implications": "UPS and dedifferentiated liposarcoma: pembrolizumab showing responses. GIST: imatinib/sunitinib for KIT/PDGFRA-mutated. Synovial sarcoma: NY-ESO-1 TCR-T cell therapy. Trabectedin modulates TAMs in L-sarcoma.",
        "key_references": ["Petitprez et al. Nature 2020", "Cancer Genome Atlas Research Network, Cell 2017"],
    },
]


def get_tme_profile_count() -> int:
    """Return the number of seed TME profiles."""
    return len(TME_PROFILES)


def get_cancer_types() -> List[str]:
    """Return cancer type names from seed data."""
    return [p["cancer_type"] for p in TME_PROFILES]


# ===================================================================
# TME PARSER
# ===================================================================


class TMEParser(BaseIngestParser):
    """Ingest parser for tumor microenvironment profiles.

    Seeds the knowledge base with TME composition data for major cancer
    types, including immune cell populations, stromal components,
    spatial patterns, and therapeutic implications.

    Usage::

        parser = TMEParser()
        records, stats = parser.run()
    """

    def __init__(
        self,
        collection_manager=None,
        embedder=None,
    ) -> None:
        super().__init__(
            source_name="tme_atlas",
            collection_manager=collection_manager,
            embedder=embedder,
        )

    def fetch(self, **kwargs) -> List[Dict[str, Any]]:
        """Return seed TME profile records."""
        return list(TME_PROFILES)

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse TME profile dictionaries into IngestRecord objects."""
        records = []
        for entry in raw_data:
            immune = entry.get("immune_composition", {})
            immune_summary = "; ".join(
                f"{k}: {v}" for k, v in immune.items()
            )
            stromal_str = ", ".join(entry.get("stromal_components", []))
            refs_str = ", ".join(entry.get("key_references", []))

            text = (
                f"Tumor microenvironment: {entry['cancer_type']}. "
                f"Immune composition: {immune_summary}. "
                f"Stromal components: {stromal_str}. "
                f"Spatial patterns: {entry.get('spatial_patterns', '')}. "
                f"Therapeutic implications: {entry.get('therapeutic_implications', '')}. "
                f"Key references: {refs_str}."
            )

            record = IngestRecord(
                text=text,
                metadata={
                    "cancer_type": entry.get("cancer_type", ""),
                    "cancer_code": entry.get("cancer_code", ""),
                    "immune_composition": immune,
                    "stromal_components": entry.get("stromal_components", []),
                    "spatial_patterns": entry.get("spatial_patterns", ""),
                    "therapeutic_implications": entry.get("therapeutic_implications", ""),
                },
                collection_name="sc_tme",
                record_id=f"tme_{entry.get('cancer_code', 'unknown')}",
                source="tme_atlas",
            )
            records.append(record)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate that a TME record has minimum required data."""
        if len(record.text) < 50:
            return False
        if not record.metadata.get("cancer_type"):
            return False
        if not record.metadata.get("immune_composition"):
            return False
        return True
