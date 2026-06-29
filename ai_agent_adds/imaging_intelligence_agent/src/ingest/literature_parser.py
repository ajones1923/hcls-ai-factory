"""PubMed literature ingest pipeline for Imaging Intelligence Agent.

Fetches medical imaging AI research papers via NCBI E-utilities
(esearch + efetch), parses PubMed XML into ImagingLiterature models,
and stores embeddings in the imaging_literature Milvus collection.

Author: Adam Jones
Date: February 2026
"""

import re
from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import (
    AITaskType,
    BodyRegion,
    ImagingLiterature,
    ImagingModality,
    SourceType,
)
from src.utils.pubmed_client import PubMedClient

from .base import BaseIngestPipeline


# ── Keyword sets for modality classification ─────────────────────────

_MODALITY_KEYWORDS: Dict[ImagingModality, List[str]] = {
    ImagingModality.CT: [
        "computed tomography",
        "ct scan",
        "ct imaging",
        "hounsfield",
        "ctdi",
        "ct angiography",
        "cta",
        "hrct",
        "high-resolution ct",
        "low-dose ct",
        "ldct",
    ],
    ImagingModality.MRI: [
        "magnetic resonance",
        "mri",
        "mr imaging",
        "flair",
        "diffusion weighted",
        "dwi",
        "t1-weighted",
        "t2-weighted",
        "gadolinium",
        "functional mri",
        "fmri",
    ],
    ImagingModality.XRAY: [
        "x-ray",
        "xray",
        "radiograph",
        "plain film",
        "projection radiography",
    ],
    ImagingModality.CXR: [
        "chest x-ray",
        "chest radiograph",
        "cxr",
        "chest x ray",
        "posteroanterior chest",
        "pa chest",
    ],
    ImagingModality.ULTRASOUND: [
        "ultrasound",
        "ultrasonography",
        "sonography",
        "doppler",
        "echocardiography",
    ],
    ImagingModality.PET: [
        "positron emission",
        "pet scan",
        "pet imaging",
        "fdg-pet",
        "fdg pet",
        "suv",
    ],
    ImagingModality.PET_CT: [
        "pet/ct",
        "pet-ct",
        "hybrid pet",
        "combined pet",
    ],
    ImagingModality.MAMMOGRAPHY: [
        "mammography",
        "mammogram",
        "digital breast",
        "tomosynthesis",
        "dbt",
    ],
    ImagingModality.FLUOROSCOPY: [
        "fluoroscopy",
        "fluoroscopic",
        "real-time imaging",
    ],
}

# ── Keyword sets for body region classification ──────────────────────

_BODY_REGION_KEYWORDS: Dict[BodyRegion, List[str]] = {
    BodyRegion.HEAD: [
        "head ct",
        "intracranial",
        "cranial",
        "skull",
        "head trauma",
    ],
    BodyRegion.BRAIN: [
        "brain",
        "cerebral",
        "cortical",
        "hippocampal",
        "white matter",
        "gray matter",
        "glioma",
        "glioblastoma",
        "meningioma",
        "alzheimer",
        "neurodegenerative",
        "stroke",
        "ischemic",
    ],
    BodyRegion.NECK: [
        "neck",
        "cervical",
        "thyroid",
        "laryngeal",
        "pharyngeal",
    ],
    BodyRegion.CHEST: [
        "chest",
        "thorax",
        "thoracic",
        "pulmonary",
        "lung",
        "pleural",
        "mediastinal",
        "pneumonia",
        "pneumothorax",
        "covid",
        "tuberculosis",
        "copd",
    ],
    BodyRegion.CARDIAC: [
        "cardiac",
        "heart",
        "coronary",
        "myocardial",
        "aortic",
        "echocardiography",
        "ejection fraction",
        "cardiomegaly",
    ],
    BodyRegion.BREAST: [
        "breast",
        "mammary",
        "bi-rads",
        "birads",
        "breast cancer screening",
        "breast lesion",
    ],
    BodyRegion.ABDOMEN: [
        "abdominal",
        "abdomen",
        "liver",
        "hepatic",
        "renal",
        "kidney",
        "pancreatic",
        "splenic",
        "bowel",
        "gallbladder",
    ],
    BodyRegion.PELVIS: [
        "pelvic",
        "pelvis",
        "prostate",
        "uterine",
        "ovarian",
        "bladder",
    ],
    BodyRegion.SPINE: [
        "spine",
        "spinal",
        "vertebral",
        "lumbar",
        "thoracolumbar",
        "disc herniation",
        "scoliosis",
    ],
    BodyRegion.EXTREMITY: [
        "extremity",
        "musculoskeletal",
        "fracture",
        "orthopedic",
        "knee",
        "hip",
        "shoulder",
        "wrist",
        "ankle",
        "bone",
    ],
    BodyRegion.WHOLE_BODY: [
        "whole body",
        "whole-body",
        "total body",
        "pan-scan",
    ],
}

# ── Keyword sets for AI task classification ──────────────────────────

_AI_TASK_KEYWORDS: Dict[str, List[str]] = {
    AITaskType.SEGMENTATION.value: [
        "segmentation",
        "segment",
        "delineation",
        "contouring",
        "region of interest",
        "roi",
        "voxel-wise",
        "pixel-wise",
        "u-net",
        "unet",
        "nnunet",
        "monai",
    ],
    AITaskType.DETECTION.value: [
        "detection",
        "detect",
        "localization",
        "computer-aided detection",
        "cad",
        "bounding box",
        "object detection",
        "lesion detection",
        "nodule detection",
    ],
    AITaskType.CLASSIFICATION.value: [
        "classification",
        "classify",
        "diagnosis",
        "diagnostic",
        "pathology classification",
        "malignancy prediction",
        "benign vs malignant",
        "grading",
        "staging",
    ],
    AITaskType.REGISTRATION.value: [
        "registration",
        "alignment",
        "atlas-based",
        "deformable registration",
        "rigid registration",
        "image fusion",
    ],
    AITaskType.RECONSTRUCTION.value: [
        "reconstruction",
        "super-resolution",
        "denoising",
        "artifact reduction",
        "image enhancement",
        "sparse-view",
        "low-dose reconstruction",
    ],
    AITaskType.TRIAGE.value: [
        "triage",
        "prioritization",
        "worklist",
        "urgent finding",
        "critical finding",
        "alert",
    ],
    AITaskType.QUANTIFICATION.value: [
        "quantification",
        "volumetric",
        "measurement",
        "biomarker",
        "radiomics",
        "texture analysis",
        "shape analysis",
    ],
    AITaskType.REPORT_GENERATION.value: [
        "report generation",
        "automated report",
        "natural language generation",
        "radiology report",
        "impression generation",
        "structured reporting",
    ],
}

# Default PubMed query for medical imaging AI literature
DEFAULT_QUERY = (
    "medical imaging AI OR radiology artificial intelligence OR "
    "CT segmentation OR chest X-ray classification OR "
    "MRI brain segmentation OR MONAI OR DICOM AI"
)


def _truncate_utf8(text: str, max_bytes: int) -> str:
    """Truncate a string to fit within max_bytes when UTF-8 encoded.

    Milvus VARCHAR max_length counts bytes, not characters, so we need
    byte-aware truncation for text containing multi-byte characters
    (Greek letters, mathematical symbols, CJK characters, etc.).
    """
    encoded = text.encode("utf-8")
    if len(encoded) <= max_bytes:
        return text
    return encoded[:max_bytes].decode("utf-8", errors="ignore")


class PubMedImagingIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for PubMed medical imaging AI literature.

    Fetches abstracts from PubMed using NCBI E-utilities, classifies each
    paper by imaging modality, body region, and AI task, and stores
    the results in the imaging_literature Milvus collection.

    Usage:
        client = PubMedClient(api_key="...")
        pipeline = PubMedImagingIngestPipeline(collection_manager, embedder, client)
        count = pipeline.run(query="chest CT segmentation", max_results=500)
    """

    COLLECTION_NAME = "imaging_literature"

    def __init__(
        self,
        collection_manager,
        embedder,
        pubmed_client: Optional[PubMedClient] = None,
    ):
        """Initialize the PubMed imaging ingest pipeline.

        Args:
            collection_manager: ImagingCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method.
            pubmed_client: PubMedClient instance.  If None, a default client
                is created.
        """
        super().__init__(collection_manager, embedder)
        self.pubmed_client = pubmed_client or PubMedClient()

    def fetch(
        self,
        query: str = DEFAULT_QUERY,
        max_results: int = 5000,
    ) -> List[Dict[str, Any]]:
        """Fetch abstracts from PubMed via NCBI E-utilities.

        Performs a two-step retrieval:
          1. esearch -- get PMIDs matching the query
          2. efetch  -- retrieve full abstract records for those PMIDs

        Args:
            query: PubMed search query string.
            max_results: Maximum number of articles to retrieve.

        Returns:
            List of dicts with keys: pmid, title, abstract, authors,
            journal, year, mesh_terms.
        """
        pmids = self.pubmed_client.search(query, max_results)
        logger.info(f"Found {len(pmids)} PMIDs for imaging query")
        articles = self.pubmed_client.fetch_abstracts(pmids)
        return articles

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[ImagingLiterature]:
        """Parse PubMed article dicts into ImagingLiterature models.

        For each article, classifies the imaging modality, body region,
        and AI task based on keyword analysis of the title and abstract.

        Args:
            raw_data: List of dicts from fetch(), each containing:
                pmid, title, abstract, authors, journal, year, mesh_terms.

        Returns:
            List of validated ImagingLiterature model instances.
        """
        records = []
        for article in raw_data:
            try:
                pmid = article.get("pmid", "")
                title = article.get("title", "")
                abstract = article.get("abstract", "")
                journal = article.get("journal", "")
                year = article.get("year", None)
                mesh_terms = article.get("mesh_terms", [])

                # Combine title + abstract as the text chunk
                text_chunk = _truncate_utf8(f"{title} {abstract}".strip(), 2990)

                # Classify modality, body region, and AI task
                modality = self._classify_modality(text_chunk)
                body_region = self._classify_body_region(text_chunk)
                ai_task = self._classify_ai_task(text_chunk)

                # Extract disease mentions from text
                disease = self._extract_disease(text_chunk)

                # Convert year to int, default to 2020 if missing/invalid
                try:
                    year = int(year)
                except (TypeError, ValueError):
                    year = 2020

                # Join mesh terms list into a semicolon-separated string
                keywords = "; ".join(mesh_terms) if mesh_terms else ""

                record = ImagingLiterature(
                    id=_truncate_utf8(pmid, 95),
                    title=_truncate_utf8(title, 490),
                    text_chunk=text_chunk,
                    source_type=SourceType.PUBMED,
                    year=year,
                    modality=modality,
                    body_region=body_region,
                    ai_task=ai_task,
                    disease=_truncate_utf8(disease, 195),
                    keywords=_truncate_utf8(keywords, 950),
                    journal=_truncate_utf8(journal, 190),
                )
                records.append(record)
            except Exception as e:
                logger.warning(f"Failed to parse article {article.get('pmid', '?')}: {e}")
                continue

        return records

    @staticmethod
    def _classify_modality(text: str) -> ImagingModality:
        """Classify an article's imaging modality by keyword matching.

        Counts keyword hits for each ImagingModality and returns the
        modality with the highest count.  Falls back to CT if no
        keywords match.

        Args:
            text: The title + abstract text to classify.

        Returns:
            The ImagingModality with the most keyword matches.
        """
        if not text:
            return ImagingModality.CT

        text_lower = text.lower()
        scores: Dict[ImagingModality, int] = {}

        for modality, keywords in _MODALITY_KEYWORDS.items():
            count = sum(1 for kw in keywords if kw.lower() in text_lower)
            scores[modality] = count

        best = max(scores, key=scores.get)  # type: ignore[arg-type]
        if scores[best] == 0:
            return ImagingModality.CT

        return best

    @staticmethod
    def _classify_body_region(text: str) -> BodyRegion:
        """Classify an article's body region by keyword matching.

        Args:
            text: The title + abstract text to classify.

        Returns:
            The BodyRegion with the most keyword matches.
        """
        if not text:
            return BodyRegion.CHEST

        text_lower = text.lower()
        scores: Dict[BodyRegion, int] = {}

        for region, keywords in _BODY_REGION_KEYWORDS.items():
            count = sum(1 for kw in keywords if kw.lower() in text_lower)
            scores[region] = count

        best = max(scores, key=scores.get)  # type: ignore[arg-type]
        if scores[best] == 0:
            return BodyRegion.CHEST

        return best

    @staticmethod
    def _classify_ai_task(text: str) -> str:
        """Classify an article's AI task by keyword matching.

        Args:
            text: The title + abstract text to classify.

        Returns:
            The AI task string value with the most keyword matches.
        """
        if not text:
            return ""

        text_lower = text.lower()
        scores: Dict[str, int] = {}

        for task, keywords in _AI_TASK_KEYWORDS.items():
            count = sum(1 for kw in keywords if kw.lower() in text_lower)
            scores[task] = count

        best = max(scores, key=scores.get)  # type: ignore[arg-type]
        if scores[best] == 0:
            return ""

        return best

    @staticmethod
    def _extract_disease(text: str) -> str:
        """Extract disease/condition mentions from abstract text.

        Uses regex patterns to find common disease terms in medical
        imaging literature.

        Args:
            text: The abstract text to search.

        Returns:
            Semicolon-separated string of found diseases, or "".
        """
        if not text:
            return ""

        disease_patterns = [
            r"\bpneumonia\b",
            r"\bcovid[-\s]?19\b",
            r"\btuberculosis\b",
            r"\blung cancer\b",
            r"\bpulmonary nodule\b",
            r"\bpleural effusion\b",
            r"\bpneumothorax\b",
            r"\bglioma\b",
            r"\bglioblastoma\b",
            r"\bmeningioma\b",
            r"\bstroke\b",
            r"\bischemic\b",
            r"\bhemorrhage\b",
            r"\balzheimer\b",
            r"\bbreast cancer\b",
            r"\bcolorectal\b",
            r"\bhepato(?:cellular)?\s*carcinoma\b",
            r"\bpancreatic cancer\b",
            r"\bprostate cancer\b",
            r"\bfracture\b",
            r"\bosteoporosis\b",
            r"\bstenosis\b",
            r"\baneurysm\b",
            r"\bcardiomegaly\b",
            r"\batelectasis\b",
        ]

        found = set()
        for pattern in disease_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                found.add(match.group(0).strip().lower())

        return "; ".join(sorted(found)) if found else ""

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        **fetch_kwargs,
    ) -> int:
        """Execute the full PubMed imaging ingest pipeline.

        Args:
            collection_name: Target collection (defaults to 'imaging_literature').
            batch_size: Batch size for embedding and insertion.
            **fetch_kwargs: Passed to fetch() (query, max_results).

        Returns:
            Total number of records ingested.
        """
        target = collection_name or self.COLLECTION_NAME
        logger.info(f"Starting PubMed imaging ingest pipeline -> {target}")

        raw = self.fetch(**fetch_kwargs)
        logger.info(f"Fetched {len(raw)} raw articles from PubMed")

        records = self.parse(raw)
        logger.info(f"Parsed {len(records)} ImagingLiterature records")

        count = self.embed_and_store(records, target, batch_size)
        logger.info(f"Ingested {count} records into {target}")
        return count
