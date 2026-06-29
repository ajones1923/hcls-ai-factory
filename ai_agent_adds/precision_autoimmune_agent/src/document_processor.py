"""
Precision Autoimmune Agent — Document Processor

Ingests clinical PDFs (progress notes, lab reports, imaging, pathology,
genetic tests, referral letters, medication lists) into Milvus collections.

Pipeline: PDF → text extraction → NLP chunking → entity extraction → embedding → insert

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from loguru import logger

try:
    from PyPDF2 import PdfReader
except ImportError:
    PdfReader = None


# ── Document type classification patterns ─────────────────────────────────

DOC_TYPE_PATTERNS = {
    "lab_report": [
        r"lab(?:oratory)?\s+report", r"test\s+results?", r"cbc\b", r"cmp\b",
        r"reference\s+range", r"specimen", r"panel",
    ],
    "progress_note": [
        r"progress\s+note", r"office\s+visit", r"follow[\s-]?up",
        r"chief\s+complaint", r"assessment\s+(?:and|&)\s+plan", r"(?:hpi|history\s+of\s+present)",
    ],
    "imaging_report": [
        r"radiology", r"imaging\s+report", r"x[\s-]?ray", r"mri\b", r"ct\s+scan",
        r"ultrasound", r"impression:", r"findings:",
    ],
    "pathology_report": [
        r"pathology", r"biopsy", r"histolog", r"microscopic", r"gross\s+description",
        r"specimen\s+received",
    ],
    "genetic_report": [
        r"hla\s+typing", r"genetic\s+test", r"genotyp", r"allele",
        r"pharmacogenomic", r"dna\s+analysis",
    ],
    "referral_letter": [
        r"referral", r"dear\s+(?:dr|doctor)", r"requesting\s+(?:evaluation|consultation)",
        r"please\s+(?:evaluate|see)",
    ],
    "medication_list": [
        r"medication\s+(?:list|reconciliation)", r"current\s+medications",
        r"prescription", r"dosage",
    ],
}

SPECIALTY_PATTERNS = {
    "rheumatology": [r"rheumatol", r"arthritis", r"lupus", r"sle\b", r"autoimmun"],
    "neurology": [r"neurolog", r"ms\b", r"multiple\s+sclerosis", r"neuropath"],
    "dermatology": [r"dermatol", r"skin", r"rash", r"psoriasis", r"lesion"],
    "nephrology": [r"nephrol", r"kidney", r"renal", r"proteinuria", r"gfr\b"],
    "gastroenterology": [r"gastro", r"gi\b", r"bowel", r"crohn", r"colitis"],
    "ophthalmology": [r"ophthal", r"eye", r"uveitis", r"schirmer"],
    "endocrinology": [r"endocrin", r"thyroid", r"diabetes", r"graves", r"hashimoto"],
    "cardiology": [r"cardiol", r"heart", r"pots\b", r"tachycardia", r"echo(?:cardiogram)?"],
    "allergy_immunology": [r"allerg", r"immunol", r"mast\s+cell", r"anaphylax"],
    "pulmonology": [r"pulmon", r"lung", r"ild\b", r"pfts?\b", r"spirometr"],
    "primary_care": [r"primary\s+care", r"pcp\b", r"family\s+medicine", r"internal\s+medicine"],
}

# Autoantibody names for entity extraction
AUTOANTIBODY_NAMES = [
    "ANA", "anti-dsDNA", "anti-Smith", "RF", "anti-CCP", "anti-Scl-70",
    "anti-centromere", "anti-SSA", "anti-SSB", "anti-Ro", "anti-La",
    "anti-Jo-1", "AChR", "anti-tTG", "TSI", "anti-TPO", "ANCA",
    "anti-cardiolipin", "lupus anticoagulant", "anti-beta2-glycoprotein",
    "anti-RNP", "anti-histone", "anti-Pm-Scl", "anti-RNA Polymerase III",
    "anti-MuSK", "c-ANCA", "p-ANCA", "PR3", "MPO",
]

# Lab test names for extraction
LAB_TEST_PATTERNS = {
    "CRP": r"c[\s-]?reactive\s+protein|crp\b",
    "ESR": r"(?:erythrocyte\s+)?sed(?:imentation)?\s+rate|esr\b",
    "complement_C3": r"complement\s+c3|c3\s+level|c3\b",
    "complement_C4": r"complement\s+c4|c4\s+level|c4\b",
    "WBC": r"white\s+blood\s+cell|wbc\b",
    "hemoglobin": r"hemoglobin|hgb\b|hb\b",
    "platelet": r"platelet|plt\b",
    "creatinine": r"creatinine|cr\b",
    "albumin": r"albumin\b",
    "TSH": r"thyroid\s+stimulat|tsh\b",
    "IL-6": r"interleukin[\s-]?6|il[\s-]?6\b",
    "calprotectin": r"calprotectin|fecal\s+calprotectin",
    "ALT": r"alanine\s+(?:amino)?transaminase|alt\b|sgpt\b",
    "AST": r"aspartate\s+(?:amino)?transaminase|ast\b|sgot\b",
    "ALP": r"alkaline\s+phosphatase|alp\b|alk\s*phos\b",
    "bilirubin": r"bilirubin\b|bili\b",
    "BUN": r"blood\s+urea\s+nitrogen|bun\b",
    "ferritin": r"ferritin\b",
    "iron": r"(?:serum\s+)?iron\b|fe\b",
    "TIBC": r"total\s+iron\s+binding|tibc\b",
    "vitamin_D": r"vitamin\s+d|25[\s-]?oh[\s-]?d|cholecalciferol",
    "HbA1c": r"hba1c\b|hemoglobin\s+a1c|glycated\s+hemoglobin|a1c\b",
    "uric_acid": r"uric\s+acid|urate\b",
    "LDH": r"lactate\s+dehydrogenase|ldh\b",
    "CK": r"creatine\s+kinase|ck\b|cpk\b",
    "aldolase": r"aldolase\b",
    "procalcitonin": r"procalcitonin\b|pct\b",
    "IgG": r"immunoglobulin\s+g|igg\b",
    "IgA": r"immunoglobulin\s+a|iga\b",
    "IgM": r"immunoglobulin\s+m|igm\b",
    "complement_CH50": r"ch50\b|total\s+(?:hemolytic\s+)?complement",
    "anti_dsDNA_titer": r"anti[\s-]?ds[\s-]?dna\s+(?:titer|level|quantitative)",
    "RF_quantitative": r"rheumatoid\s+factor\s+(?:quantitative|level|titer)",
    "free_T4": r"free\s+t4|ft4\b|free\s+thyroxine",
    "free_T3": r"free\s+t3|ft3\b|free\s+triiodothyronine",
    "C_peptide": r"c[\s-]?peptide\b",
    "neurofilament_light": r"neurofilament\s+light|nfl\b|nf[\s-]?l\b",
    "beta2_microglobulin": r"beta[\s-]?2[\s-]?microglobulin|b2m\b",
    "D_dimer": r"d[\s-]?dimer\b",
    "fibrinogen": r"fibrinogen\b",
    "haptoglobin": r"haptoglobin\b",
    "reticulocyte": r"reticulocyte\b|retic\b",
    "troponin": r"troponin\b|tnni\b|tnt\b",
    "NT_proBNP": r"nt[\s-]?pro[\s-]?bnp\b|n[\s-]?terminal[\s-]?pro[\s-]?bnp",
    "proteinuria": r"protein(?:uria)?[\s/]creatinine|upcr\b|24[\s-]?(?:hr|hour)\s+(?:urine\s+)?protein",
}


class DocumentProcessor:
    """Process clinical PDFs into Milvus-ready chunks."""

    def __init__(
        self,
        collection_manager=None,
        embedder=None,
        max_chunk_size: int = 2500,
        chunk_overlap: int = 200,
    ):
        self.collection_manager = collection_manager
        self.embedder = embedder
        self.max_chunk_size = max_chunk_size
        self.chunk_overlap = chunk_overlap

    # ── PDF extraction ────────────────────────────────────────────────
    def extract_text_from_pdf(self, pdf_path: Path) -> str:
        """Extract text from a PDF file."""
        if PdfReader is None:
            raise ImportError("PyPDF2 is required for PDF processing")
        reader = PdfReader(str(pdf_path))
        pages = []
        for page in reader.pages:
            text = page.extract_text()
            if text:
                pages.append(text.strip())
        return "\n\n".join(pages)

    def extract_pages_from_pdf(self, pdf_path: Path) -> List[Tuple[int, str]]:
        """Extract text from each page individually."""
        if PdfReader is None:
            raise ImportError("PyPDF2 is required for PDF processing")
        reader = PdfReader(str(pdf_path))
        pages = []
        for i, page in enumerate(reader.pages):
            text = page.extract_text()
            if text and text.strip():
                pages.append((i + 1, text.strip()))
        return pages

    # ── Text chunking ─────────────────────────────────────────────────
    def chunk_text(self, text: str) -> List[str]:
        """Split text into overlapping chunks at sentence boundaries."""
        sentences = re.split(r'(?<=[.!?])\s+', text)
        chunks = []
        current = []
        current_len = 0

        for sentence in sentences:
            if current_len + len(sentence) > self.max_chunk_size and current:
                chunks.append(" ".join(current))
                # Keep overlap
                overlap_text = " ".join(current)
                overlap_start = max(0, len(overlap_text) - self.chunk_overlap)
                overlap = overlap_text[overlap_start:]
                current = [overlap] if overlap else []
                current_len = len(overlap)
            current.append(sentence)
            current_len += len(sentence) + 1

        if current:
            chunks.append(" ".join(current))

        return chunks

    # ── Document classification ───────────────────────────────────────
    def classify_document_type(self, text: str) -> str:
        """Classify document type from text content."""
        text_lower = text.lower()
        scores = {}
        for doc_type, patterns in DOC_TYPE_PATTERNS.items():
            score = sum(1 for p in patterns if re.search(p, text_lower))
            if score > 0:
                scores[doc_type] = score
        if scores:
            return max(scores, key=scores.get)
        return "clinical_note"

    def detect_specialty(self, text: str) -> str:
        """Detect medical specialty from document content."""
        text_lower = text.lower()
        scores = {}
        for specialty, patterns in SPECIALTY_PATTERNS.items():
            score = sum(1 for p in patterns if re.search(p, text_lower))
            if score > 0:
                scores[specialty] = score
        if scores:
            return max(scores, key=scores.get)
        return "general"

    def extract_date(self, text: str) -> str:
        """Try to extract a visit/collection date from text."""
        patterns = [
            r"(?:date|visit|collected|drawn)[\s:]*(\d{1,2}[/-]\d{1,2}[/-]\d{2,4})",
            r"(\d{4}-\d{2}-\d{2})",
            r"((?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\w*\s+\d{1,2},?\s+\d{4})",
        ]
        for pattern in patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                return match.group(1)
        return ""

    def extract_provider(self, text: str) -> str:
        """Extract provider name from document."""
        patterns = [
            r"(?:provider|physician|doctor|attending|seen\s+by)[\s:]*(?:dr\.?\s*)?([A-Z][a-z]+(?:\s+[A-Z][a-z]+)+)",
            r"(?:Dr\.?\s+)([A-Z][a-z]+(?:\s+[A-Z][a-z]+)+)",
        ]
        for pattern in patterns:
            match = re.search(pattern, text)
            if match:
                return match.group(1).strip()
        return ""

    def extract_patient_id_from_path(self, pdf_path: Path) -> str:
        """Derive patient ID from file path (e.g., demo_data/sarah_mitchell/...)."""
        parts = pdf_path.parts
        for i, part in enumerate(parts):
            if part == "demo_data" and i + 1 < len(parts):
                return parts[i + 1]  # e.g., "sarah_mitchell"
        return pdf_path.parent.name

    # ── Entity extraction ─────────────────────────────────────────────
    def extract_autoantibodies(self, text: str) -> List[Dict[str, Any]]:
        """Extract autoantibody mentions and their values."""
        found = []
        text_lower = text.lower()
        for ab in AUTOANTIBODY_NAMES:
            if ab.lower() in text_lower:
                # Try to find value
                value = None
                pattern = rf"{re.escape(ab)}[:\s]*([<>]?\s*\d+\.?\d*)"
                match = re.search(pattern, text, re.IGNORECASE)
                if match:
                    try:
                        cleaned = re.sub(r"[<>\s]", "", match.group(1))
                        if cleaned:
                            value = float(cleaned)
                    except (ValueError, TypeError):
                        pass

                # Check positivity
                positive = False
                pos_pattern = rf"{re.escape(ab)}[^.]*(?:positive|detected|reactive|elevated|abnormal)"
                if re.search(pos_pattern, text, re.IGNORECASE):
                    positive = True

                # Check for titer
                titer = None
                titer_pattern = rf"{re.escape(ab)}[^.]*?(1:\d+)"
                titer_match = re.search(titer_pattern, text, re.IGNORECASE)
                if titer_match:
                    titer = titer_match.group(1)
                    positive = True

                found.append({
                    "antibody": ab,
                    "value": value,
                    "positive": positive,
                    "titer": titer,
                })
        return found

    def extract_lab_values(self, text: str) -> List[Dict[str, Any]]:
        """Extract lab test values from text."""
        found = []
        for test_name, pattern in LAB_TEST_PATTERNS.items():
            matches = re.finditer(
                rf"({pattern})[:\s]*([<>]?\s*\d+\.?\d*)\s*([a-zA-Z/%]+)?",
                text,
                re.IGNORECASE,
            )
            for match in matches:
                try:
                    value = float(re.sub(r"[<>\s]", "", match.group(2)))
                    unit = match.group(3) or ""
                    found.append({
                        "test_name": test_name,
                        "value": value,
                        "unit": unit,
                    })
                except (ValueError, IndexError):
                    continue
        return found

    # ── Full ingestion pipeline ───────────────────────────────────────
    def process_pdf(self, pdf_path: Path, patient_id: Optional[str] = None) -> List[Dict[str, Any]]:
        """Process a single PDF into Milvus-ready records.

        Returns list of records (dicts) ready for insertion into
        autoimmune_clinical_documents collection.
        """
        pdf_path = Path(pdf_path)
        if not pdf_path.exists():
            logger.warning(f"PDF not found: {pdf_path}")
            return []

        if patient_id is None:
            patient_id = self.extract_patient_id_from_path(pdf_path)

        logger.debug(f"Processing: {pdf_path.name} (patient: {patient_id})")

        pages = self.extract_pages_from_pdf(pdf_path)
        if not pages:
            logger.warning(f"No text extracted from: {pdf_path.name}")
            return []

        full_text = "\n\n".join(text for _, text in pages)
        doc_type = self.classify_document_type(full_text)
        specialty = self.detect_specialty(full_text)
        visit_date = self.extract_date(full_text)
        provider = self.extract_provider(full_text)

        records = []
        for page_num, page_text in pages:
            chunks = self.chunk_text(page_text)
            for chunk_idx, chunk in enumerate(chunks):
                record_id = hashlib.md5(
                    f"{pdf_path.name}:{page_num}:{chunk_idx}".encode()
                ).hexdigest()[:16]

                records.append({
                    "id": f"doc_{patient_id}_{record_id}",
                    "text_chunk": chunk[:3000],
                    "patient_id": patient_id,
                    "doc_type": doc_type,
                    "specialty": specialty,
                    "provider": provider,
                    "visit_date": visit_date,
                    "source_file": pdf_path.name,
                    "page_number": page_num,
                    "chunk_index": chunk_idx,
                })

        logger.info(
            f"Processed {pdf_path.name}: {len(records)} chunks "
            f"(type={doc_type}, specialty={specialty})"
        )
        return records

    def process_directory(
        self,
        directory: Path,
        patient_id: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Process all PDFs in a directory."""
        directory = Path(directory)
        all_records = []
        pdf_files = sorted(directory.glob("*.pdf"))
        if not pdf_files:
            pdf_files = sorted(directory.glob("**/*.pdf"))

        for pdf_path in pdf_files:
            pid = patient_id or self.extract_patient_id_from_path(pdf_path)
            records = self.process_pdf(pdf_path, pid)
            all_records.extend(records)

        logger.info(f"Processed directory {directory.name}: {len(all_records)} total chunks from {len(pdf_files)} PDFs")
        return all_records

    def embed_records(self, records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Add embeddings to records using the configured embedder."""
        if self.embedder is None:
            raise RuntimeError("No embedder configured")

        texts = [r["text_chunk"] for r in records]
        embeddings = self.embedder.encode(texts, batch_size=32, show_progress_bar=False)

        for rec, emb in zip(records, embeddings):
            rec["embedding"] = emb.tolist()

        return records

    def ingest_patient(
        self,
        patient_dir: Path,
        patient_id: Optional[str] = None,
        collection_name: str = "autoimmune_clinical_documents",
    ) -> int:
        """Full pipeline: PDF dir → extract → chunk → embed → insert."""
        records = self.process_directory(patient_dir, patient_id)
        if not records:
            return 0

        records = self.embed_records(records)
        count = self.collection_manager.insert_batch(collection_name, records)
        logger.info(f"Ingested {count} records for patient {patient_id or patient_dir.name}")
        return count

    def ingest_demo_data(self, demo_dir: Path, collection_name: str = "autoimmune_clinical_documents") -> Dict[str, int]:
        """Ingest all demo patient directories."""
        demo_dir = Path(demo_dir)
        results = {}
        for patient_dir in sorted(demo_dir.iterdir()):
            if patient_dir.is_dir() and not patient_dir.name.startswith("."):
                count = self.ingest_patient(patient_dir, patient_dir.name, collection_name)
                results[patient_dir.name] = count
        return results
