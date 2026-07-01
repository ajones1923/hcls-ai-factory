"""
Precision Autoimmune Agent — Timeline Builder

Constructs patient diagnostic timelines from ingested clinical documents.
Extracts temporal events, builds chronological narratives, and identifies
patterns in the diagnostic odyssey.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import re
from datetime import datetime
from typing import Any, Dict, List, Optional

# Event type classification patterns
EVENT_PATTERNS = {
    "symptom_onset": [
        r"(?:first|initial|new)\s+(?:complaint|symptom|presentation)",
        r"(?:onset|began|started|developed)\s+(?:of\s+)?(?:pain|swelling|fatigue|rash)",
        r"chief\s+complaint",
    ],
    "diagnosis": [
        r"(?:diagnosed|diagnosis)\s+(?:of|with|:)",
        r"confirmed\s+(?:diagnosis|finding)",
        r"meets?\s+(?:criteria|classification)\s+for",
    ],
    "misdiagnosis": [
        r"(?:previously|initially)\s+(?:diagnosed|labeled|treated)\s+(?:as|for|with)",
        r"(?:revised|changed|updated)\s+diagnosis",
        r"(?:incorrect|wrong|erroneous)\s+diagnosis",
    ],
    "lab_result": [
        r"lab(?:oratory)?\s+result",
        r"(?:positive|negative|elevated|low|abnormal)\s+(?:for\s+)?(?:ANA|anti-|RF|CRP|ESR|complement)",
    ],
    "imaging": [
        r"(?:x[\s-]?ray|mri|ct|ultrasound|echo)\s+(?:shows?|reveals?|demonstrates?)",
        r"imaging\s+(?:findings?|results?)",
        r"radiology\s+report",
    ],
    "biopsy": [
        r"biopsy\s+(?:shows?|reveals?|confirms?)",
        r"pathology\s+(?:report|findings?)",
        r"histolog(?:y|ical)\s+(?:findings?|examination)",
    ],
    "genetic_test": [
        r"hla\s+(?:typing|test|result)",
        r"genetic\s+(?:test|panel|analysis)",
        r"pharmacogenomic",
    ],
    "treatment_start": [
        r"(?:started|initiated|began|prescribed)\s+(?:on\s+)?(?:\w+(?:mab|cept|nib|tinib))",
        r"commenced\s+(?:treatment|therapy)",
        r"trial\s+of\s+\w+",
    ],
    "treatment_change": [
        r"(?:switched|changed|transitioned)\s+(?:to|from)",
        r"(?:added|augmented)\s+(?:with|by)",
        r"(?:discontinued|stopped|tapered)",
    ],
    "flare": [
        r"flare[\s-]?up",
        r"(?:disease|symptom)\s+(?:exacerbation|worsening|recurrence)",
        r"acute\s+(?:episode|relapse)",
    ],
    "referral": [
        r"referred?\s+to\s+(?:\w+\s+)?(?:rheumatol|neurolog|dermatol|nephrolog|gastro|ophthal|immunol)",
        r"consultation\s+(?:with|requested)",
    ],
    "er_visit": [
        r"emergency\s+(?:room|department|visit)",
        r"(?:er|ed)\s+(?:visit|presentation)",
        r"acute\s+(?:visit|presentation|care)",
    ],
}

DATE_PATTERNS = [
    (r"(\d{4}-\d{2}-\d{2})", "%Y-%m-%d"),
    (r"(\d{1,2}/\d{1,2}/\d{4})", "%m/%d/%Y"),
    (r"(\d{1,2}/\d{1,2}/\d{2})", "%m/%d/%y"),
    (r"((?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\w*\s+\d{1,2},?\s+\d{4})", "%B %d, %Y"),
]


class TimelineBuilder:
    """Builds patient diagnostic timelines from clinical document chunks."""

    def __init__(self, rag_engine=None):
        self.rag = rag_engine

    def extract_date(self, text: str) -> Optional[str]:
        """Extract the most likely date from text."""
        for pattern, fmt in DATE_PATTERNS:
            match = re.search(pattern, text)
            if match:
                try:
                    date_str = match.group(1).replace(",", "")
                    dt = datetime.strptime(date_str, fmt.replace(",", ""))
                    return dt.strftime("%Y-%m-%d")
                except ValueError:
                    continue
        return None

    def classify_event(self, text: str) -> str:
        """Classify the type of clinical event described in text."""
        text_lower = text.lower()
        scores = {}
        for event_type, patterns in EVENT_PATTERNS.items():
            score = sum(1 for p in patterns if re.search(p, text_lower))
            if score > 0:
                scores[event_type] = score
        if scores:
            return max(scores, key=scores.get)
        return "clinical_note"

    def extract_events_from_chunks(
        self,
        chunks: List[Dict[str, Any]],
        patient_id: str = "",
    ) -> List[Dict[str, Any]]:
        """Extract timeline events from document chunks."""
        events = []

        for chunk in chunks:
            text = chunk.get("text_chunk", chunk.get("text", ""))
            if not text:
                continue

            event_type = self.classify_event(text)
            event_date = self.extract_date(text) or chunk.get("visit_date", "")
            specialty = chunk.get("specialty", "")
            provider = chunk.get("provider", "")

            # Extract a brief description
            description = self._summarize_event(text, event_type)

            events.append({
                "patient_id": patient_id or chunk.get("patient_id", ""),
                "event_type": event_type,
                "event_date": event_date,
                "description": description,
                "provider": provider,
                "specialty": specialty,
                "source_file": chunk.get("source_file", ""),
                "text_chunk": text[:500],
            })

        # Sort by date
        events.sort(key=lambda e: e.get("event_date", "9999"))

        # Assign days from first symptom
        first_date = None
        for e in events:
            if e["event_type"] == "symptom_onset" and e["event_date"]:
                first_date = e["event_date"]
                break
        if first_date is None and events and events[0].get("event_date"):
            first_date = events[0]["event_date"]

        if first_date:
            try:
                d0 = datetime.strptime(first_date[:10], "%Y-%m-%d")
                for e in events:
                    if e.get("event_date"):
                        try:
                            dt = datetime.strptime(e["event_date"][:10], "%Y-%m-%d")
                            e["days_from_first_symptom"] = (dt - d0).days
                        except ValueError:
                            e["days_from_first_symptom"] = 0
                    else:
                        e["days_from_first_symptom"] = 0
            except ValueError:
                for e in events:
                    e["days_from_first_symptom"] = 0

        return events

    def _summarize_event(self, text: str, event_type: str) -> str:
        """Create a brief summary of the event."""
        # Take first 2 meaningful sentences
        sentences = re.split(r'(?<=[.!?])\s+', text.strip())
        meaningful = [s for s in sentences if len(s) > 20][:2]
        if meaningful:
            return " ".join(meaningful)[:300]
        return text[:300]

    def build_timeline(
        self,
        patient_id: str,
        document_chunks: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Build a complete patient timeline from document chunks."""
        events = self.extract_events_from_chunks(document_chunks, patient_id)

        # Aggregate statistics
        specialties = set()
        event_types = {}
        date_range = {"first": None, "last": None}

        for e in events:
            if e.get("specialty"):
                specialties.add(e["specialty"])
            et = e["event_type"]
            event_types[et] = event_types.get(et, 0) + 1
            date = e.get("event_date")
            if date:
                if date_range["first"] is None or date < date_range["first"]:
                    date_range["first"] = date
                if date_range["last"] is None or date > date_range["last"]:
                    date_range["last"] = date

        return {
            "patient_id": patient_id,
            "total_events": len(events),
            "events": events,
            "specialties_seen": sorted(specialties),
            "event_type_counts": event_types,
            "date_range": date_range,
        }

    def build_timeline_for_milvus(
        self,
        patient_id: str,
        events: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """Convert timeline events into Milvus-ready records."""
        records = []
        for i, event in enumerate(events):
            record = {
                "id": f"tl_{patient_id}_{i:04d}",
                "text_chunk": (
                    f"[{event.get('event_date', 'unknown')}] "
                    f"{event['event_type'].replace('_', ' ').title()}: "
                    f"{event.get('description', '')}"
                )[:3000],
                "patient_id": patient_id,
                "event_type": event.get("event_type", ""),
                "event_date": event.get("event_date", ""),
                "description": event.get("description", "")[:2000],
                "provider": event.get("provider", ""),
                "specialty": event.get("specialty", ""),
                "days_from_first_symptom": event.get("days_from_first_symptom", 0),
            }
            records.append(record)
        return records
