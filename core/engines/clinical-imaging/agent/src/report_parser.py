"""Radiology Report NLP parser for Imaging Intelligence Agent.

Parses free-text radiology reports into structured data: sections,
entities (findings, anatomy, laterality, measurements), and critical
finding detection.  Output is suitable for embedding with BGE-small-en-v1.5
and storage in the imaging_reports Milvus collection.

Follows the same regex + Pydantic pattern as:
  - src/ingest/literature_parser.py (keyword-based classification)
  - src/ingest/finding_parser.py (finding extraction)

Author: Adam Jones
Date: April 2026
"""

import re
from typing import Dict, List, Optional, Tuple

from loguru import logger

from src.models import (
    Measurement,
    ParsedReport,
    ParsedSection,
    ReportEntity,
)


# ═══════════════════════════════════════════════════════════════════════
# SECTION HEADER PATTERNS
# ═══════════════════════════════════════════════════════════════════════

# Common radiology report section headers.  Each pattern maps to a
# canonical section name used as a dict key in ParsedReport.sections.
_SECTION_HEADERS: List[Tuple[str, re.Pattern]] = [
    ("examination", re.compile(
        r"^\s*(EXAMINATION|EXAM|PROCEDURE|TYPE OF EXAMINATION|STUDY)\s*[:.]?\s*",
        re.IGNORECASE | re.MULTILINE,
    )),
    ("clinical_history", re.compile(
        r"^\s*(CLINICAL\s+HISTORY|HISTORY|CLINICAL\s+INFORMATION|"
        r"INDICATION|CLINICAL\s+INDICATION|REASON\s+FOR\s+EXAM|"
        r"REASON\s+FOR\s+STUDY)\s*[:.]?\s*",
        re.IGNORECASE | re.MULTILINE,
    )),
    ("technique", re.compile(
        r"^\s*(TECHNIQUE|PROTOCOL|ACQUISITION)\s*[:.]?\s*",
        re.IGNORECASE | re.MULTILINE,
    )),
    ("comparison", re.compile(
        r"^\s*(COMPARISON|COMPARISONS|PRIOR\s+STUDIES?|PRIOR\s+EXAM)\s*[:.]?\s*",
        re.IGNORECASE | re.MULTILINE,
    )),
    ("findings", re.compile(
        r"^\s*(FINDINGS|BODY\s+OF\s+REPORT|DESCRIPTION)\s*[:.]?\s*",
        re.IGNORECASE | re.MULTILINE,
    )),
    ("impression", re.compile(
        r"^\s*(IMPRESSION|CONCLUSION|CONCLUSIONS|SUMMARY|ASSESSMENT)\s*[:.]?\s*",
        re.IGNORECASE | re.MULTILINE,
    )),
]

# ═══════════════════════════════════════════════════════════════════════
# MEASUREMENT PATTERNS
# ═══════════════════════════════════════════════════════════════════════

# Pattern: "8mm nodule", "3.2 cm mass", "15 mL effusion", "density 45 HU"
# Also handles adjectives between measurement and finding: "4.5 cm hypodense lesion"
_MEASUREMENT_FINDING_RE = re.compile(
    r"(\d+\.?\d*)\s*(mm|cm|mL|ml|HU|cc)\s+"
    r"(?:\w+\s+)?"  # Optional adjective (e.g., "hypodense", "hyperdense")
    r"(nodule|mass|lesion|effusion|thickness|diameter|volume|shift|measurement|"
    r"cyst|calcification|opacity|collection|hematoma|abscess|aneurysm|"
    r"lymph\s*node)",
    re.IGNORECASE,
)

# Pattern: "density 45 HU", "attenuation 35 HU"
_MEASUREMENT_DENSITY_RE = re.compile(
    r"(?:density|attenuation)\s+(\d+\.?\d*)\s*(HU)",
    re.IGNORECASE,
)

# Pattern: finding word BEFORE measurement: "shift 5.2 mm", "volume 35 mL"
_MEASUREMENT_REVERSED_RE = re.compile(
    r"(shift|volume|thickness|diameter)\s+(?:measures?\s+|of\s+|approximately\s+)?"
    r"(\d+\.?\d*)\s*(mm|cm|mL|ml|cc)",
    re.IGNORECASE,
)

# Pattern: "measures 3.2 x 2.1 x 1.8 cm" or "measuring 12 x 8 mm"
_MEASUREMENT_DIMENSIONS_RE = re.compile(
    r"(?:measur(?:es|ing)|dimensions?)\s+"
    r"(\d+\.?\d*)\s*x\s*(\d+\.?\d*)"
    r"(?:\s*x\s*(\d+\.?\d*))?"
    r"\s*(mm|cm)",
    re.IGNORECASE,
)

# Pattern: "8 mm" or "3.2cm" standalone measurement near a finding word
_MEASUREMENT_STANDALONE_RE = re.compile(
    r"(\d+\.?\d*)\s*-?\s*(mm|cm)\s+"
    r"(part-solid|ground-glass|solid|calcified|cystic|spiculated|"
    r"lobulated|irregular|well-circumscribed|subcentimeter)?\s*"
    r"(nodule|mass|lesion|opacity|density|calcification|cyst)",
    re.IGNORECASE,
)

# ═══════════════════════════════════════════════════════════════════════
# ENTITY EXTRACTION PATTERNS
# ═══════════════════════════════════════════════════════════════════════

# Finding keywords (nouns describing pathological observations)
_FINDING_PATTERNS = re.compile(
    r"\b(nodule|mass|lesion|effusion|consolidation|opacity|"
    r"pneumothorax|atelectasis|calcification|fracture|hemorrhage|"
    r"stenosis|aneurysm|thrombus|thrombosis|embolism|infarct|"
    r"edema|fibrosis|emphysema|bronchiectasis|lymphadenopathy|"
    r"cardiomegaly|hernia|herniation|obstruction|perforation|"
    r"dissection|abscess|hematoma|cyst|polyp|diverticulum|"
    r"metastasis|metastases|tumor|neoplasm|malignancy|"
    r"ground-glass|ground.glass|tree-in-bud|air.trapping|"
    r"honeycombing|traction bronchiectasis|mosaic attenuation|"
    r"interstitial thickening|septal thickening)\b",
    re.IGNORECASE,
)

# Anatomical structures
_ANATOMY_PATTERNS = re.compile(
    r"\b(right upper lobe|right middle lobe|right lower lobe|"
    r"left upper lobe|left lower lobe|lingula|"
    r"right lung|left lung|bilateral lungs|"
    r"right hilum|left hilum|"
    r"aorta|aortic arch|ascending aorta|descending aorta|"
    r"pulmonary artery|pulmonary arteries|"
    r"trachea|main.?stem bronchus|carina|"
    r"esophagus|thyroid|"
    r"right atrium|left atrium|right ventricle|left ventricle|"
    r"pericardium|pericardial|"
    r"liver|spleen|gallbladder|pancreas|"
    r"right kidney|left kidney|adrenal|"
    r"thoracic spine|lumbar spine|cervical spine|"
    r"sternum|ribs|clavicle|scapula|"
    r"axilla|axillary|mediastinum|mediastinal|"
    r"pleura|pleural space|"
    r"diaphragm|peritoneum|retroperitoneum|"
    r"basal ganglia|cerebellum|cerebrum|brainstem|"
    r"frontal lobe|temporal lobe|parietal lobe|occipital lobe|"
    r"lateral ventricle|third ventricle|fourth ventricle|"
    r"carotid|vertebral artery|circle of Willis|"
    r"right breast|left breast|"
    r"prostate|bladder|uterus|ovary|ovaries)\b",
    re.IGNORECASE,
)

# Laterality markers
_LATERALITY_PATTERNS = re.compile(
    r"\b(right|left|bilateral|midline|right-sided|left-sided)\b",
    re.IGNORECASE,
)

# Recommendation keywords
_RECOMMENDATION_PATTERNS = re.compile(
    r"\b(recommend|recommended|suggest|suggested|advise|advised|"
    r"consider|should be|follow.?up|correlate|clinical correlation|"
    r"further evaluation|additional imaging|biopsy|tissue sampling|"
    r"short.?interval|repeat|return|surveillance|screening|"
    r"PET/?CT|MRI|CT|ultrasound)\b",
    re.IGNORECASE,
)

# ═══════════════════════════════════════════════════════════════════════
# CRITICAL FINDING KEYWORDS
# ═══════════════════════════════════════════════════════════════════════

_CRITICAL_KEYWORDS = [
    # Acute life-threatening
    "acute intracranial hemorrhage",
    "acute hemorrhage",
    "intracranial hemorrhage",
    "subarachnoid hemorrhage",
    "subdural hematoma",
    "epidural hematoma",
    "intraparenchymal hemorrhage",
    "hemorrhagic stroke",
    "ischemic stroke",
    "acute stroke",
    "acute infarct",
    "large vessel occlusion",
    # Cardiothoracic emergencies
    "aortic dissection",
    "aortic rupture",
    "pulmonary embolism",
    "tension pneumothorax",
    "cardiac tamponade",
    "pericardial tamponade",
    "massive pleural effusion",
    "acute myocardial infarction",
    # Abdominal emergencies
    "bowel perforation",
    "free air",
    "pneumoperitoneum",
    "bowel obstruction",
    "mesenteric ischemia",
    "ruptured aneurysm",
    "active extravasation",
    "active hemorrhage",
    "active bleeding",
    # Spine emergencies
    "cord compression",
    "spinal cord compression",
    "cauda equina",
    "unstable fracture",
    # Oncology critical
    "highly suspicious for malignancy",
    "lung-rads 4b",
    "lung-rads 4x",
    "bi-rads 5",
    "bi-rads 6",
    "li-rads 5",
    "li-rads lr-5",
    # General critical
    "emergent",
    "critical finding",
    "urgent finding",
    "immediate attention",
    "life-threatening",
    "requires immediate",
    "stat notification",
]

# ═══════════════════════════════════════════════════════════════════════
# MODALITY DETECTION
# ═══════════════════════════════════════════════════════════════════════

_MODALITY_KEYWORDS: Dict[str, List[str]] = {
    "ct": ["ct ", "ct,", "computed tomography", "helical ct", "non-contrast ct",
           "contrast-enhanced ct", "ctdi", "ct scan"],
    "mri": ["mri", "magnetic resonance", "mr imaging", "flair", "t1-weighted",
            "t2-weighted", "dwi", "diffusion", "gadolinium"],
    "xray": ["x-ray", "xray", "radiograph", "plain film", "chest film"],
    "cxr": ["chest x-ray", "chest radiograph", "pa and lateral", "portable chest"],
    "ultrasound": ["ultrasound", "ultrasonography", "sonography", "doppler"],
    "pet_ct": ["pet/ct", "pet-ct", "fdg-pet/ct"],
    "pet": ["pet ", "pet,", "pet scan", "fdg-pet"],
    "mammography": ["mammography", "mammogram", "tomosynthesis"],
    "fluoroscopy": ["fluoroscopy", "fluoroscopic"],
}

# ═══════════════════════════════════════════════════════════════════════
# BODY REGION DETECTION
# ═══════════════════════════════════════════════════════════════════════

_BODY_REGION_KEYWORDS: Dict[str, List[str]] = {
    "head": ["head ct", "head ", "cranial", "skull", "intracranial"],
    "brain": ["brain", "cerebral", "cortical", "hippocampal", "glioma",
              "meningioma", "ventricle"],
    "neck": ["neck", "cervical", "thyroid", "laryngeal"],
    "chest": ["chest", "thorax", "thoracic", "pulmonary", "lung", "pleural",
              "mediastinal", "lobe"],
    "cardiac": ["cardiac", "heart", "coronary", "myocardial", "aortic",
                "pericardial", "ejection fraction"],
    "breast": ["breast", "mammary", "bi-rads", "birads"],
    "abdomen": ["abdominal", "abdomen", "liver", "hepatic", "renal", "kidney",
                "pancreatic", "splenic", "bowel", "gallbladder"],
    "pelvis": ["pelvic", "pelvis", "prostate", "uterine", "ovarian", "bladder"],
    "spine": ["spine", "spinal", "vertebral", "lumbar", "thoracolumbar",
              "disc herniation", "scoliosis"],
    "extremity": ["extremity", "musculoskeletal", "fracture", "orthopedic",
                  "knee", "hip", "shoulder", "wrist", "ankle"],
    "whole_body": ["whole body", "whole-body", "total body"],
}

# ═══════════════════════════════════════════════════════════════════════
# COMPARISON DATE EXTRACTION
# ═══════════════════════════════════════════════════════════════════════

_DATE_PATTERNS = [
    # ISO: 2025-06-15
    re.compile(r"(\d{4}-\d{2}-\d{2})"),
    # US: 06/15/2025, 6/15/2025
    re.compile(r"(\d{1,2}/\d{1,2}/\d{4})"),
    # Written: June 15, 2025
    re.compile(
        r"((?:January|February|March|April|May|June|July|August|September|"
        r"October|November|December)\s+\d{1,2},?\s+\d{4})",
        re.IGNORECASE,
    ),
    # Abbreviated: Jun 15, 2025
    re.compile(
        r"((?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z]*\.?\s+"
        r"\d{1,2},?\s+\d{4})",
        re.IGNORECASE,
    ),
]

# ═══════════════════════════════════════════════════════════════════════
# QUALIFIER EXTRACTION
# ═══════════════════════════════════════════════════════════════════════

_QUALIFIER_RE = re.compile(
    r"\b(part-solid|part solid|ground-glass|ground glass|solid|"
    r"calcified|cystic|spiculated|lobulated|irregular|"
    r"well-circumscribed|well circumscribed|subcentimeter|"
    r"sub-centimeter|fat-containing|fat containing|"
    r"heterogeneous|homogeneous|enhancing|non-enhancing|"
    r"cavitary|necrotic|hemorrhagic|rim-enhancing)\b",
    re.IGNORECASE,
)

# ═══════════════════════════════════════════════════════════════════════
# LOCATION EXTRACTION (for measurements)
# ═══════════════════════════════════════════════════════════════════════

_LOCATION_RE = re.compile(
    r"\b(?:in\s+the\s+|of\s+the\s+|within\s+the\s+|at\s+the\s+)?"
    r"(right upper lobe|right middle lobe|right lower lobe|"
    r"left upper lobe|left lower lobe|lingula|"
    r"right lung|left lung|"
    r"right hilum|left hilum|"
    r"right lobe of (?:the )?liver|left lobe of (?:the )?liver|"
    r"caudate lobe|"
    r"right kidney|left kidney|"
    r"right adrenal|left adrenal|"
    r"right breast|left breast|"
    r"right ventricle|left ventricle|"
    r"right atrium|left atrium|"
    r"right basal ganglia|left basal ganglia|"
    r"right frontal lobe|left frontal lobe|"
    r"right temporal lobe|left temporal lobe|"
    r"right parietal lobe|left parietal lobe|"
    r"right occipital lobe|left occipital lobe|"
    r"head of (?:the )?pancreas|body of (?:the )?pancreas|"
    r"tail of (?:the )?pancreas|"
    r"mediastinum|anterior mediastinum|posterior mediastinum|"
    r"upper abdomen|lower abdomen|pelvis|"
    r"right axilla|left axilla|"
    r"supraclavicular|infraclavicular)\b",
    re.IGNORECASE,
)


# ═══════════════════════════════════════════════════════════════════════
# RADIOLOGY REPORT PARSER
# ═══════════════════════════════════════════════════════════════════════


class RadiologyReportParser:
    """NLP parser for free-text radiology reports.

    Extracts structured data from unstructured radiology report text:
      - Section segmentation (EXAMINATION, FINDINGS, IMPRESSION, etc.)
      - Entity extraction (findings, anatomy, laterality, recommendations)
      - Measurement extraction with units and qualifiers
      - Critical finding detection
      - Modality and body region classification
      - Comparison date extraction

    All extraction uses regex patterns — no external NLP models required.
    Set REPORT_EXTRACT_WITH_LLM=True in settings to use an LLM for
    entity extraction (not yet implemented; placeholder for future use).

    Usage:
        parser = RadiologyReportParser()
        parsed = parser.parse(report_text)
        embedding_text = parser.generate_embedding_text(parsed)
    """

    def __init__(self) -> None:
        """Initialize the parser.

        Pre-compiles and stores references to all regex patterns used
        for section splitting, entity extraction, measurement parsing,
        and critical finding detection.
        """
        self._section_headers = _SECTION_HEADERS
        self._critical_keywords = [kw.lower() for kw in _CRITICAL_KEYWORDS]
        logger.debug("RadiologyReportParser initialized")

    # ── Public API ──────────────────────────────────────────────────

    def parse(self, report_text: str) -> ParsedReport:
        """Parse a radiology report into structured components.

        This is the main entry point.  It segments the report into
        named sections, extracts entities and measurements from the
        findings and impression sections, detects critical findings,
        and classifies the modality and body region.

        Args:
            report_text: The raw radiology report text (free-text).

        Returns:
            ParsedReport with sections, entities, measurements,
            critical_finding flag, modality, body_region, and
            comparison_date.
        """
        if not report_text or not report_text.strip():
            logger.warning("Empty report text provided")
            return ParsedReport(
                sections={},
                entities=[],
                measurements=[],
                critical_finding=False,
                modality=None,
                body_region=None,
                comparison_date=None,
            )

        # Normalize whitespace
        text = self._normalize_text(report_text)

        # Step 1: Segment into sections
        sections = self._segment_sections(text)
        logger.debug(f"Segmented {len(sections)} sections: {list(sections.keys())}")

        # Step 2: Extract entities from findings and impression
        entities = self._extract_entities(sections)
        logger.debug(f"Extracted {len(entities)} entities")

        # Step 3: Extract measurements from all text
        measurements = self._extract_measurements(text)
        logger.debug(f"Extracted {len(measurements)} measurements")

        # Step 4: Detect critical findings
        critical = self._detect_critical_findings(sections)

        # Step 5: Classify modality
        modality = self._detect_modality(sections)

        # Step 6: Classify body region
        body_region = self._detect_body_region(sections)

        # Step 7: Extract comparison date
        comparison_date = self._extract_comparison_date(sections)

        parsed = ParsedReport(
            sections=sections,
            entities=entities,
            measurements=measurements,
            critical_finding=critical,
            modality=modality,
            body_region=body_region,
            comparison_date=comparison_date,
        )

        logger.info(
            f"Parsed report: {len(sections)} sections, "
            f"{len(entities)} entities, {len(measurements)} measurements, "
            f"critical={critical}, modality={modality}, region={body_region}"
        )
        return parsed

    def generate_embedding_text(self, parsed: ParsedReport) -> str:
        """Generate a text summary suitable for BGE-small-en-v1.5 embedding.

        Combines the most informative sections (examination, findings,
        impression) with extracted measurements and critical finding
        status into a single string optimized for semantic search.

        Args:
            parsed: A ParsedReport from self.parse().

        Returns:
            A text string suitable for embedding (typically 100-500 tokens).
        """
        parts: List[str] = []

        # Modality and body region context
        if parsed.modality:
            parts.append(f"Modality: {parsed.modality}")
        if parsed.body_region:
            parts.append(f"Region: {parsed.body_region}")

        # Examination line
        if "examination" in parsed.sections:
            parts.append(f"Examination: {parsed.sections['examination']}")

        # Clinical history (brief)
        if "clinical_history" in parsed.sections:
            history = parsed.sections["clinical_history"]
            # Truncate long histories
            if len(history) > 200:
                history = history[:200] + "..."
            parts.append(f"History: {history}")

        # Findings (the most important section for embedding)
        if "findings" in parsed.sections:
            findings = parsed.sections["findings"]
            if len(findings) > 1000:
                findings = findings[:1000] + "..."
            parts.append(f"Findings: {findings}")

        # Impression (concise clinical summary)
        if "impression" in parsed.sections:
            parts.append(f"Impression: {parsed.sections['impression']}")

        # Key measurements
        if parsed.measurements:
            meas_strs = []
            for m in parsed.measurements[:5]:  # Limit to top 5
                s = f"{m.value}{m.unit} {m.finding}"
                if m.qualifier:
                    s = f"{m.qualifier} {s}"
                if m.location:
                    s += f" in {m.location}"
                meas_strs.append(s)
            parts.append(f"Measurements: {'; '.join(meas_strs)}")

        # Critical finding flag
        if parsed.critical_finding:
            parts.append("CRITICAL FINDING PRESENT")

        # Comparison date
        if parsed.comparison_date:
            parts.append(f"Comparison: {parsed.comparison_date}")

        return " | ".join(parts)

    # ── Section Segmentation ────────────────────────────────────────

    def _segment_sections(self, text: str) -> Dict[str, str]:
        """Split report text into named sections using regex header matching.

        Identifies section boundaries by matching common radiology
        report headers (EXAMINATION, FINDINGS, IMPRESSION, etc.) and
        extracts the text between consecutive headers.

        If no section headers are found, the entire text is placed
        under a "findings" key as a fallback.

        Args:
            text: Normalized report text.

        Returns:
            Dict mapping canonical section names to their text content.
            Keys are lowercase: "examination", "clinical_history",
            "technique", "comparison", "findings", "impression".
        """
        # Find all header positions
        header_positions: List[Tuple[int, int, str]] = []  # (start, end, section_name)

        for section_name, pattern in self._section_headers:
            for match in pattern.finditer(text):
                header_positions.append((match.start(), match.end(), section_name))

        # Sort by position in the text
        header_positions.sort(key=lambda x: x[0])

        if not header_positions:
            # No recognized headers found — treat entire text as findings
            logger.debug("No section headers found, using full text as findings")
            return {"findings": text.strip()}

        sections: Dict[str, str] = {}

        for i, (start, end, section_name) in enumerate(header_positions):
            # Text runs from end of this header to start of next header
            if i + 1 < len(header_positions):
                section_text = text[end:header_positions[i + 1][0]]
            else:
                section_text = text[end:]

            section_text = section_text.strip()

            # If the same section appears multiple times, append
            if section_name in sections and section_text:
                sections[section_name] += "\n" + section_text
            elif section_text:
                sections[section_name] = section_text

        return sections

    # ── Entity Extraction ───────────────────────────────────────────

    def _extract_entities(self, sections: Dict[str, str]) -> List[ReportEntity]:
        """Extract clinical entities from findings and impression sections.

        Looks for five entity types:
          - finding: pathological observations (nodule, mass, effusion, etc.)
          - anatomy: anatomical structures (right upper lobe, liver, etc.)
          - laterality: spatial qualifiers (right, left, bilateral)
          - recommendation: suggested follow-up actions
          - measurement: numeric measurements (delegated to _extract_measurements)

        Args:
            sections: Dict of section_name -> text from _segment_sections().

        Returns:
            List of ReportEntity instances, deduplicated by (type, value).
        """
        entities: List[ReportEntity] = []
        seen: set = set()  # (entity_type, value_lower) for dedup

        # Focus on findings and impression — the most entity-rich sections
        target_sections = ["findings", "impression", "clinical_history"]
        target_text = ""
        for sec_name in target_sections:
            if sec_name in sections:
                target_text += " " + sections[sec_name]

        if not target_text.strip():
            # Fall back to all available sections
            target_text = " ".join(sections.values())

        # Extract findings
        for match in _FINDING_PATTERNS.finditer(target_text):
            value = match.group(1).strip()
            key = ("finding", value.lower())
            if key not in seen:
                seen.add(key)
                context = self._get_context(target_text, match.start(), match.end())
                entities.append(ReportEntity(
                    entity_type="finding",
                    value=value,
                    context=context,
                ))

        # Extract anatomy
        for match in _ANATOMY_PATTERNS.finditer(target_text):
            value = match.group(1).strip()
            key = ("anatomy", value.lower())
            if key not in seen:
                seen.add(key)
                context = self._get_context(target_text, match.start(), match.end())
                entities.append(ReportEntity(
                    entity_type="anatomy",
                    value=value,
                    context=context,
                ))

        # Extract laterality
        for match in _LATERALITY_PATTERNS.finditer(target_text):
            value = match.group(1).strip()
            key = ("laterality", value.lower())
            if key not in seen:
                seen.add(key)
                context = self._get_context(target_text, match.start(), match.end())
                entities.append(ReportEntity(
                    entity_type="laterality",
                    value=value,
                    context=context,
                ))

        # Extract recommendations (from impression only)
        impression_text = sections.get("impression", "")
        for match in _RECOMMENDATION_PATTERNS.finditer(impression_text):
            value = match.group(1).strip()
            key = ("recommendation", value.lower())
            if key not in seen:
                seen.add(key)
                context = self._get_context(impression_text, match.start(), match.end())
                entities.append(ReportEntity(
                    entity_type="recommendation",
                    value=value,
                    context=context,
                ))

        return entities

    # ── Measurement Extraction ──────────────────────────────────────

    def _extract_measurements(self, text: str) -> List[Measurement]:
        """Extract numeric measurements with units from report text.

        Handles several patterns:
          - "8mm nodule" -> Measurement(8.0, "mm", "nodule")
          - "3.2 cm mass" -> Measurement(3.2, "cm", "mass")
          - "15 mL effusion" -> Measurement(15.0, "mL", "effusion")
          - "density 45 HU" -> Measurement(45.0, "HU", "density")
          - "measures 3.2 x 2.1 x 1.8 cm" -> Measurement(3.2, "cm", "dimension")

        Also extracts qualifiers (part-solid, ground-glass, etc.) and
        anatomical locations from surrounding context.

        Args:
            text: Full report text to search for measurements.

        Returns:
            List of Measurement instances, deduplicated.
        """
        measurements: List[Measurement] = []
        seen: set = set()

        # Pattern 1: "8mm nodule", "3.2 cm mass"
        for match in _MEASUREMENT_FINDING_RE.finditer(text):
            value = float(match.group(1))
            unit = match.group(2)
            finding = match.group(3).lower()
            key = (value, unit.lower(), finding)
            if key not in seen:
                seen.add(key)
                # Look for qualifier and location in surrounding context
                context = self._get_context(text, match.start(), match.end(), window=100)
                qualifier = self._extract_qualifier(context)
                location = self._extract_location(context)
                measurements.append(Measurement(
                    value=value,
                    unit=unit,
                    finding=finding,
                    qualifier=qualifier,
                    location=location,
                ))

        # Pattern 2: standalone "8 mm part-solid nodule"
        for match in _MEASUREMENT_STANDALONE_RE.finditer(text):
            value = float(match.group(1))
            unit = match.group(2)
            qualifier_str = match.group(3)
            finding = match.group(4).lower()
            key = (value, unit.lower(), finding)
            if key not in seen:
                seen.add(key)
                context = self._get_context(text, match.start(), match.end(), window=100)
                qualifier = qualifier_str.strip() if qualifier_str else self._extract_qualifier(context)
                location = self._extract_location(context)
                measurements.append(Measurement(
                    value=value,
                    unit=unit,
                    finding=finding,
                    qualifier=qualifier,
                    location=location,
                ))

        # Pattern 3: reversed — "shift 5.2 mm", "volume 35 mL"
        for match in _MEASUREMENT_REVERSED_RE.finditer(text):
            finding = match.group(1).lower()
            value = float(match.group(2))
            unit = match.group(3)
            key = (value, unit.lower(), finding)
            if key not in seen:
                seen.add(key)
                context = self._get_context(text, match.start(), match.end(), window=100)
                location = self._extract_location(context)
                measurements.append(Measurement(
                    value=value,
                    unit=unit,
                    finding=finding,
                    qualifier=None,
                    location=location,
                ))

        # Pattern 4: "density 45 HU"
        for match in _MEASUREMENT_DENSITY_RE.finditer(text):
            value = float(match.group(1))
            unit = match.group(2)
            key = (value, unit.lower(), "density")
            if key not in seen:
                seen.add(key)
                context = self._get_context(text, match.start(), match.end(), window=100)
                location = self._extract_location(context)
                measurements.append(Measurement(
                    value=value,
                    unit=unit,
                    finding="density",
                    qualifier=None,
                    location=location,
                ))

        # Pattern 5: "measures 3.2 x 2.1 x 1.8 cm"
        for match in _MEASUREMENT_DIMENSIONS_RE.finditer(text):
            dim1 = float(match.group(1))
            dim2 = float(match.group(2))
            dim3_str = match.group(3)
            unit = match.group(4)

            # Use the largest dimension as the primary measurement
            dims = [dim1, dim2]
            if dim3_str:
                dims.append(float(dim3_str))
            max_dim = max(dims)

            # Build a descriptive finding string
            if dim3_str:
                finding_desc = f"dimension ({dim1} x {dim2} x {dim3_str} {unit})"
            else:
                finding_desc = f"dimension ({dim1} x {dim2} {unit})"

            key = (max_dim, unit.lower(), "dimension")
            if key not in seen:
                seen.add(key)
                context = self._get_context(text, match.start(), match.end(), window=100)
                qualifier = self._extract_qualifier(context)
                location = self._extract_location(context)
                measurements.append(Measurement(
                    value=max_dim,
                    unit=unit,
                    finding=finding_desc,
                    qualifier=qualifier,
                    location=location,
                ))

        return measurements

    # ── Critical Finding Detection ──────────────────────────────────

    def _detect_critical_findings(self, sections: Dict[str, str]) -> bool:
        """Check for critical/urgent findings in impression and findings sections.

        Scans the impression section (primary) and findings section
        (secondary) for keywords indicating critical or urgent pathology
        that may require immediate clinical notification.

        Args:
            sections: Dict of section_name -> text.

        Returns:
            True if any critical finding keyword is detected.
        """
        # Check impression first (most important for critical findings)
        check_text = ""
        for sec_name in ["impression", "findings"]:
            if sec_name in sections:
                check_text += " " + sections[sec_name]

        if not check_text:
            # Fall back to full text
            check_text = " ".join(sections.values())

        check_lower = check_text.lower()

        for keyword in self._critical_keywords:
            if keyword in check_lower:
                logger.warning(f"Critical finding detected: '{keyword}'")
                return True

        return False

    # ── Modality Detection ──────────────────────────────────────────

    def _detect_modality(self, sections: Dict[str, str]) -> Optional[str]:
        """Classify imaging modality from examination and technique sections.

        Uses keyword matching against the examination line and technique
        description to determine the imaging modality (CT, MRI, etc.).

        Args:
            sections: Dict of section_name -> text.

        Returns:
            Modality string (e.g., "ct", "mri") or None if undetectable.
        """
        # Prioritize examination and technique sections
        check_text = ""
        for sec_name in ["examination", "technique", "clinical_history"]:
            if sec_name in sections:
                check_text += " " + sections[sec_name]

        if not check_text:
            check_text = " ".join(sections.values())

        check_lower = check_text.lower()

        best_modality: Optional[str] = None
        best_count = 0

        for modality, keywords in _MODALITY_KEYWORDS.items():
            count = sum(1 for kw in keywords if kw in check_lower)
            if count > best_count:
                best_count = count
                best_modality = modality

        return best_modality

    # ── Body Region Detection ───────────────────────────────────────

    def _detect_body_region(self, sections: Dict[str, str]) -> Optional[str]:
        """Classify anatomical body region from examination and findings sections.

        Args:
            sections: Dict of section_name -> text.

        Returns:
            Body region string (e.g., "chest", "brain") or None.
        """
        check_text = ""
        for sec_name in ["examination", "findings", "clinical_history"]:
            if sec_name in sections:
                check_text += " " + sections[sec_name]

        if not check_text:
            check_text = " ".join(sections.values())

        check_lower = check_text.lower()

        best_region: Optional[str] = None
        best_count = 0

        for region, keywords in _BODY_REGION_KEYWORDS.items():
            count = sum(1 for kw in keywords if kw in check_lower)
            if count > best_count:
                best_count = count
                best_region = region

        return best_region

    # ── Comparison Date Extraction ──────────────────────────────────

    def _extract_comparison_date(self, sections: Dict[str, str]) -> Optional[str]:
        """Extract comparison study date from the comparison section.

        Looks for date patterns in the COMPARISON section, falling back
        to searching the full report text.

        Args:
            sections: Dict of section_name -> text.

        Returns:
            Date string as found in the report, or None.
        """
        # Check comparison section first
        check_text = sections.get("comparison", "")
        if not check_text:
            # Fall back to searching for "dated" or "prior" in all text
            all_text = " ".join(sections.values())
            dated_match = re.search(r"dated\s+(.{10,30})", all_text, re.IGNORECASE)
            if dated_match:
                check_text = dated_match.group(1)
            else:
                return None

        for pattern in _DATE_PATTERNS:
            match = pattern.search(check_text)
            if match:
                return match.group(1)

        return None

    # ── Helper Methods ──────────────────────────────────────────────

    @staticmethod
    def _normalize_text(text: str) -> str:
        """Normalize whitespace and line endings in report text.

        Collapses multiple whitespace characters into single spaces
        while preserving newlines that separate sections.

        Args:
            text: Raw report text.

        Returns:
            Normalized text with consistent whitespace.
        """
        # Normalize line endings
        text = text.replace("\r\n", "\n").replace("\r", "\n")
        # Collapse multiple blank lines into single newline
        text = re.sub(r"\n{3,}", "\n\n", text)
        # Collapse multiple spaces (but not newlines) into single space
        text = re.sub(r"[^\S\n]+", " ", text)
        return text.strip()

    @staticmethod
    def _get_context(text: str, start: int, end: int, window: int = 60) -> str:
        """Extract surrounding context for an entity match.

        Args:
            text: The full text being searched.
            start: Start index of the match.
            end: End index of the match.
            window: Number of characters before and after the match.

        Returns:
            Substring with context around the match, stripped.
        """
        ctx_start = max(0, start - window)
        ctx_end = min(len(text), end + window)
        return text[ctx_start:ctx_end].strip()

    @staticmethod
    def _extract_qualifier(context: str) -> Optional[str]:
        """Extract a finding qualifier from surrounding context.

        Args:
            context: Text surrounding a measurement or finding.

        Returns:
            Qualifier string (e.g., "part-solid") or None.
        """
        match = _QUALIFIER_RE.search(context)
        return match.group(1).strip() if match else None

    @staticmethod
    def _extract_location(context: str) -> Optional[str]:
        """Extract anatomical location from surrounding context.

        Args:
            context: Text surrounding a measurement or finding.

        Returns:
            Location string (e.g., "right upper lobe") or None.
        """
        match = _LOCATION_RE.search(context)
        return match.group(1).strip() if match else None
