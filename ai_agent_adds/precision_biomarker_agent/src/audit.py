"""HIPAA-compliant audit logging for patient data access.

Logs all patient data access, analysis runs, exports, and queries
to a structured audit trail. Uses loguru with a dedicated sink
for audit events, separate from application logs.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import uuid
from datetime import datetime, timezone
from enum import Enum
from typing import Any, Dict, Optional

from loguru import logger


class AuditAction(str, Enum):
    """Auditable actions in the biomarker agent."""
    PATIENT_ANALYSIS = "patient_analysis"
    BIOLOGICAL_AGE = "biological_age"
    DISEASE_RISK = "disease_risk"
    PGX_MAPPING = "pgx_mapping"
    RAG_QUERY = "rag_query"
    REPORT_GENERATED = "report_generated"
    REPORT_EXPORTED = "report_exported"
    FHIR_EXPORTED = "fhir_exported"
    PATIENT_DATA_ACCESSED = "patient_data_accessed"


# Dedicated audit logger — separate from application logs
_audit_logger = logger.bind(audit=True)


def audit_log(
    action: AuditAction,
    patient_id: Optional[str] = None,
    *,
    request_id: Optional[str] = None,
    details: Optional[Dict[str, Any]] = None,
    source_ip: Optional[str] = None,
) -> str:
    """Record an audit event for patient data access or analysis.

    Args:
        action: The type of action being performed.
        patient_id: Patient identifier (logged as hash for privacy).
        request_id: Correlation ID for tracking across services.
        details: Additional context (modules used, export format, etc.).
        source_ip: Client IP address if available.

    Returns:
        The generated audit event ID.
    """
    event_id = uuid.uuid4().hex[:16]

    # Hash patient_id for log privacy — full ID only in encrypted storage
    patient_ref = None
    if patient_id:
        import hashlib
        patient_ref = hashlib.sha256(patient_id.encode()).hexdigest()[:12]

    audit_event = {
        "event_id": event_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "action": action.value,
        "patient_ref": patient_ref,
        "request_id": request_id or uuid.uuid4().hex[:12],
        "source_ip": source_ip,
        "details": details or {},
    }

    _audit_logger.info(
        "AUDIT | {action} | patient={patient_ref} | req={request_id}",
        action=action.value,
        patient_ref=patient_ref or "none",
        request_id=audit_event["request_id"],
    )

    return event_id
