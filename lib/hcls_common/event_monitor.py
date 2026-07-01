"""Live Event Bus Monitor for cross-agent intelligence visualization.

Reads the event bus log files and provides a real-time view of cross-agent
communication for demo and debugging purposes.
"""

from __future__ import annotations

import json
import os
import glob
from datetime import datetime, date
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, field


@dataclass
class EventRecord:
    """A single event from the event bus log."""
    timestamp: str
    event_type: str
    source_stage: str
    patient_id: Optional[str]
    priority: str
    payload: Dict[str, Any]

    @property
    def short_summary(self) -> str:
        """One-line summary for display."""
        patient = f" [{self.patient_id}]" if self.patient_id else ""
        return f"{self.timestamp} | {self.source_stage} -> {self.event_type}{patient}"


def read_event_log(log_dir: str = "data/events", max_events: int = 100) -> List[EventRecord]:
    """Read recent events from the event bus log files.

    The event bus writes JSONL files named events_YYYYMMDD.jsonl.
    """
    events = []

    # Find log files
    patterns = [
        os.path.join(log_dir, "events_*.jsonl"),
        os.path.join(log_dir, "*.jsonl"),
    ]

    log_files = []
    for pattern in patterns:
        log_files.extend(glob.glob(pattern))

    if not log_files:
        # Try common locations
        for base in [
            "/home/adam/projects/hcls-ai-factory/core/agents/precision-biomarker/data/events",
            "/home/adam/projects/hcls-ai-factory/data/events",
        ]:
            for pattern in patterns:
                log_files.extend(glob.glob(os.path.join(base, "*.jsonl")))

    # Sort by filename (date-based) and read most recent
    log_files.sort(reverse=True)

    for log_file in log_files[:5]:  # Last 5 days
        try:
            with open(log_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        data = json.loads(line)
                        events.append(EventRecord(
                            timestamp=data.get("timestamp", ""),
                            event_type=data.get("event_type", data.get("type", "UNKNOWN")),
                            source_stage=data.get("source_stage", data.get("source", "UNKNOWN")),
                            patient_id=data.get("patient_id"),
                            priority=data.get("priority", "NORMAL"),
                            payload=data.get("payload", data.get("data", {})),
                        ))
                    except json.JSONDecodeError:
                        continue
        except FileNotFoundError:
            continue

    # Sort by timestamp descending, limit
    events.sort(key=lambda e: e.timestamp, reverse=True)
    return events[:max_events]


def get_event_stats(events: List[EventRecord]) -> Dict[str, Any]:
    """Compute aggregate statistics from event records."""
    if not events:
        return {"total": 0, "by_type": {}, "by_source": {}, "by_priority": {}}

    by_type: Dict[str, int] = {}
    by_source: Dict[str, int] = {}
    by_priority: Dict[str, int] = {}
    patients: set = set()

    for e in events:
        by_type[e.event_type] = by_type.get(e.event_type, 0) + 1
        by_source[e.source_stage] = by_source.get(e.source_stage, 0) + 1
        by_priority[e.priority] = by_priority.get(e.priority, 0) + 1
        if e.patient_id:
            patients.add(e.patient_id)

    return {
        "total": len(events),
        "unique_patients": len(patients),
        "by_type": dict(sorted(by_type.items(), key=lambda x: -x[1])),
        "by_source": dict(sorted(by_source.items(), key=lambda x: -x[1])),
        "by_priority": by_priority,
    }


def render_event_flow_ascii(events: List[EventRecord]) -> str:
    """Generate ASCII art showing event flow between agents."""
    agents = set()
    flows: Dict[str, Dict[str, int]] = {}

    # Map source stages to agent names
    stage_to_agent = {
        "ANNOTATION": "Biomarker",
        "BIOMARKER_ANALYSIS": "Biomarker",
        "ONCOLOGY_ANALYSIS": "Oncology",
        "CART_ANALYSIS": "CAR-T",
        "IMAGING_ANALYSIS": "Imaging",
        "AUTOIMMUNE_ANALYSIS": "Autoimmune",
        "DRUG_DISCOVERY": "Drug Disc.",
        "RAG_QUERY": "RAG/Chat",
        "VARIANT_CALLING": "Genomics",
    }

    # Map event types to target agents
    event_targets = {
        "BIOLOGICAL_AGE_COMPUTED": "Patient 360",
        "BIOMARKER_ALERT": "Oncology",
        "PGX_RESULT_READY": "Drug Disc.",
        "DISEASE_TRAJECTORY_ALERT": "Oncology",
        "ONCOLOGY_CASE_CREATED": "CAR-T",
        "THERAPY_RANKED": "Drug Disc.",
        "CART_MANUFACTURING_READY": "Patient 360",
        "PGX_DRUG_FILTER": "Drug Disc.",
    }

    for e in events:
        source = stage_to_agent.get(e.source_stage, e.source_stage)
        target = event_targets.get(e.event_type, "Patient 360")
        agents.add(source)
        agents.add(target)

        key = f"{source}->{target}"
        if key not in flows:
            flows[key] = {"count": 0, "types": set()}
        flows[key]["count"] += 1
        flows[key]["types"].add(e.event_type)

    # Build ASCII
    lines = ["Cross-Agent Event Flow", "=" * 50]
    for flow_key, info in sorted(flows.items(), key=lambda x: -x[1]["count"]):
        types_str = ", ".join(sorted(info["types"]))
        lines.append(f"  {flow_key} ({info['count']}x)")
        lines.append(f"    Events: {types_str}")

    if not flows:
        lines.append("  No cross-agent events recorded yet.")
        lines.append("  Run an analysis to generate events.")

    return "\n".join(lines)
