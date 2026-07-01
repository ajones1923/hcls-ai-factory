"""
Bidirectional Pipeline Triggers for the HCLS AI Factory.

Implements cross-stage trigger rules that enable automatic, event-driven
handoffs between the three pipeline stages and the CAR-T Intelligence Agent.
Addresses Section 10.3 Gap 1 (bidirectional triggers), Gap 2 (dynamic target
identification), and Section 11.0.3 (event-driven coordination).

Trigger rules:
    1. VCF -> RAG Ingest:       VCF ready => auto-trigger annotation + Milvus embedding
    2. RAG -> Drug Discovery:   Pathogenic variant => auto-trigger target hypothesis
                                => MolMIM molecule generation
    3. RAG <-> CAR-T:           Cell-therapy query => CAR-T queries genomic_evidence
                                for patient variants; CAR-T identifies targets =>
                                queries RAG for variant context
    4. Drug Discovery -> Report: Top candidates ranked => auto-trigger report generation
    5. Dynamic Target ID:       Replace hardcoded 11-gene lookup with dynamic query:
                                find top pathogenic/likely-pathogenic variants by gene,
                                cross-reference with knowledge graph, return
                                variant-driven targets

Hardware target: NVIDIA DGX Spark (GB10 GPU, 128 GB unified memory, $4,699).
"""

from __future__ import annotations

import time
import threading
import uuid
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum, unique
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from loguru import logger

from .event_bus import (
    EventBus,
    EventHandler,
    EventPriority,
    EventStatus,
    EventType,
    PipelineEvent,
    PipelineStage,
    get_event_bus,
)

try:
    from prometheus_client import Counter, Histogram

    TRIGGERS_FIRED = Counter(
        "hcls_triggers_fired_total",
        "Total bidirectional triggers fired",
        ["trigger_name"],
    )
    TRIGGER_LATENCY = Histogram(
        "hcls_trigger_action_seconds",
        "Trigger action execution latency",
        ["trigger_name"],
        buckets=[0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 30.0, 60.0],
    )
except ImportError:
    TRIGGERS_FIRED = None
    TRIGGER_LATENCY = None


# ---------------------------------------------------------------------------
# Trigger rule data class
# ---------------------------------------------------------------------------

@dataclass
class TriggerRule:
    """
    A rule that fires an action when matching events arrive on the EventBus.

    Attributes
    ----------
    name : str
        Human-readable rule name for logging and audit.
    event_types : set[EventType]
        The event types that can trigger this rule.
    condition : callable
        A predicate ``(event: PipelineEvent) -> bool``.  The rule fires
        only when the condition returns True.
    action : callable
        The action ``(event: PipelineEvent) -> None`` executed when the
        rule fires.
    cooldown_seconds : float
        Minimum time between consecutive firings of the same rule.
        Prevents cascading / duplicate triggers.
    enabled : bool
        Whether the rule is currently active.
    description : str
        Documentation for the rule.
    rule_id : str
        Unique identifier.
    """

    name: str
    event_types: set[EventType]
    condition: Callable[[PipelineEvent], bool]
    action: Callable[[PipelineEvent], None]
    cooldown_seconds: float = 10.0
    enabled: bool = True
    description: str = ""
    rule_id: str = field(default_factory=lambda: str(uuid.uuid4())[:8])

    def __post_init__(self) -> None:
        if not self.description:
            self.description = f"Trigger rule: {self.name}"


@dataclass
class TriggerExecution:
    """Audit record of a single trigger firing."""

    rule_name: str
    rule_id: str
    event_id: str
    event_type: str
    fired_at: str
    duration_seconds: float
    success: bool
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "rule_name": self.rule_name,
            "rule_id": self.rule_id,
            "event_id": self.event_id,
            "event_type": self.event_type,
            "fired_at": self.fired_at,
            "duration_seconds": round(self.duration_seconds, 4),
            "success": self.success,
            "error": self.error,
        }


# ---------------------------------------------------------------------------
# Bidirectional Trigger Manager
# ---------------------------------------------------------------------------

class BidirectionalTriggerManager:
    """
    Manages trigger rules that enable bidirectional, event-driven handoffs
    across the HCLS AI Factory pipeline stages.

    The manager registers itself as an EventBus subscriber and evaluates
    all registered rules against incoming events.  When a rule's condition
    is met, its action fires (subject to cooldown enforcement).

    Usage::

        from hcls_common.bidirectional_triggers import (
            BidirectionalTriggerManager,
            register_builtin_triggers,
        )

        manager = BidirectionalTriggerManager()
        register_builtin_triggers(manager)
        manager.start()

        # Now the manager listens on the EventBus and fires triggers
        # automatically as events arrive.

    Parameters
    ----------
    event_bus : EventBus, optional
        Shared event bus instance.  Default: singleton.
    enable_audit : bool
        Whether to record trigger executions for audit trail.
    max_audit_size : int
        Maximum number of audit records to keep in memory.
    """

    def __init__(
        self,
        event_bus: EventBus | None = None,
        enable_audit: bool = True,
        max_audit_size: int = 5000,
    ):
        self._bus = event_bus or get_event_bus()
        self._rules: dict[str, TriggerRule] = {}
        self._lock = threading.RLock()
        self._running = False

        # Cooldown tracking: rule_id -> last_fired_timestamp
        self._last_fired: dict[str, float] = {}

        # Audit trail
        self._enable_audit = enable_audit
        self._audit: list[TriggerExecution] = []
        self._max_audit_size = max_audit_size

        # Statistics
        self._stats: dict[str, int] = defaultdict(int)

        logger.info("BidirectionalTriggerManager initialized")

    # ---- lifecycle ---------------------------------------------------------

    def start(self) -> None:
        """
        Register the manager as an EventBus subscriber and begin
        evaluating triggers against incoming events.
        """
        if self._running:
            logger.warning("TriggerManager already running")
            return

        self._running = True

        # Subscribe to ALL event types -- we evaluate rules internally
        self._bus.subscribe(
            handler=self._on_event,
            event_types=None,  # all types
            name="bidirectional_trigger_manager",
        )

        logger.info(
            f"TriggerManager started: {len(self._rules)} rules registered"
        )

    def stop(self) -> None:
        """Unsubscribe from the EventBus and stop evaluating triggers."""
        self._running = False
        self._bus.unsubscribe("bidirectional_trigger_manager")
        logger.info("TriggerManager stopped")

    # ---- rule management ---------------------------------------------------

    def register_rule(self, rule: TriggerRule) -> str:
        """
        Register a trigger rule.

        Returns the rule_id.
        """
        with self._lock:
            self._rules[rule.rule_id] = rule
        logger.info(
            f"Registered trigger rule: '{rule.name}' (id={rule.rule_id}), "
            f"events={[e.value for e in rule.event_types]}"
        )
        return rule.rule_id

    def unregister_rule(self, rule_id: str) -> bool:
        """
        Remove a trigger rule by ID.

        Returns True if a rule was removed.
        """
        with self._lock:
            if rule_id in self._rules:
                name = self._rules[rule_id].name
                del self._rules[rule_id]
                logger.info(f"Unregistered trigger rule: '{name}' (id={rule_id})")
                return True
        return False

    def enable_rule(self, rule_id: str) -> bool:
        """Enable a previously disabled rule."""
        with self._lock:
            if rule_id in self._rules:
                self._rules[rule_id].enabled = True
                return True
        return False

    def disable_rule(self, rule_id: str) -> bool:
        """Disable a rule without removing it."""
        with self._lock:
            if rule_id in self._rules:
                self._rules[rule_id].enabled = False
                return True
        return False

    def list_rules(self) -> list[dict[str, Any]]:
        """List all registered rules with their status."""
        with self._lock:
            return [
                {
                    "rule_id": r.rule_id,
                    "name": r.name,
                    "event_types": [e.value for e in r.event_types],
                    "enabled": r.enabled,
                    "cooldown_seconds": r.cooldown_seconds,
                    "description": r.description,
                    "last_fired": self._last_fired.get(r.rule_id),
                }
                for r in self._rules.values()
            ]

    # ---- event handler -----------------------------------------------------

    def _on_event(self, event: PipelineEvent) -> None:
        """
        EventBus handler: evaluate all registered rules against the event.
        """
        if not self._running:
            return

        with self._lock:
            candidate_rules = [
                r for r in self._rules.values()
                if r.enabled and event.event_type in r.event_types
            ]

        for rule in candidate_rules:
            self._evaluate_and_fire(rule, event)

    def _evaluate_and_fire(self, rule: TriggerRule, event: PipelineEvent) -> None:
        """
        Evaluate a rule's condition and fire its action if matched.

        Enforces cooldown to prevent cascading triggers.
        """
        # Cooldown check
        now = time.time()
        last = self._last_fired.get(rule.rule_id, 0.0)
        if (now - last) < rule.cooldown_seconds:
            logger.debug(
                f"Trigger '{rule.name}' skipped (cooldown: "
                f"{rule.cooldown_seconds - (now - last):.1f}s remaining)"
            )
            return

        # Condition check
        try:
            condition_met = rule.condition(event)
        except Exception as exc:
            logger.error(
                f"Trigger '{rule.name}' condition evaluation failed: {exc}"
            )
            self._stats["condition_errors"] += 1
            return

        if not condition_met:
            return

        # Fire the action
        logger.info(
            f"Trigger '{rule.name}' FIRING for event {event.event_type.value} "
            f"(event_id={event.event_id})"
        )

        t0 = time.time()
        success = True
        error_msg: str | None = None

        try:
            rule.action(event)
        except Exception as exc:
            success = False
            error_msg = str(exc)
            logger.error(f"Trigger '{rule.name}' action failed: {exc}")
            self._stats["action_errors"] += 1

        elapsed = time.time() - t0

        # Update cooldown timestamp
        self._last_fired[rule.rule_id] = time.time()

        # Statistics
        self._stats["total_fired"] += 1
        if success:
            self._stats["successful_fires"] += 1
        else:
            self._stats["failed_fires"] += 1

        # Prometheus metrics
        if TRIGGERS_FIRED is not None:
            TRIGGERS_FIRED.labels(trigger_name=rule.name).inc()
        if TRIGGER_LATENCY is not None:
            TRIGGER_LATENCY.labels(trigger_name=rule.name).observe(elapsed)

        # Audit record
        if self._enable_audit:
            execution = TriggerExecution(
                rule_name=rule.name,
                rule_id=rule.rule_id,
                event_id=event.event_id,
                event_type=event.event_type.value,
                fired_at=datetime.now().isoformat(),
                duration_seconds=elapsed,
                success=success,
                error=error_msg,
            )
            self._add_audit_record(execution)

    def _add_audit_record(self, execution: TriggerExecution) -> None:
        """Add an execution record to the bounded audit trail."""
        with self._lock:
            self._audit.append(execution)
            if len(self._audit) > self._max_audit_size:
                self._audit = self._audit[-self._max_audit_size:]

    # ---- query -------------------------------------------------------------

    def get_audit_trail(
        self,
        rule_name: str | None = None,
        limit: int = 100,
    ) -> list[dict[str, Any]]:
        """
        Return recent trigger execution records.

        Parameters
        ----------
        rule_name : str, optional
            Filter by rule name.
        limit : int
            Maximum records to return.
        """
        with self._lock:
            records = list(self._audit)

        if rule_name:
            records = [r for r in records if r.rule_name == rule_name]

        records.sort(key=lambda r: r.fired_at, reverse=True)
        return [r.to_dict() for r in records[:limit]]

    def get_stats(self) -> dict[str, Any]:
        """Return trigger manager statistics."""
        with self._lock:
            return {
                "total_rules": len(self._rules),
                "enabled_rules": sum(1 for r in self._rules.values() if r.enabled),
                "total_fired": self._stats.get("total_fired", 0),
                "successful_fires": self._stats.get("successful_fires", 0),
                "failed_fires": self._stats.get("failed_fires", 0),
                "condition_errors": self._stats.get("condition_errors", 0),
                "action_errors": self._stats.get("action_errors", 0),
                "audit_records": len(self._audit),
                "running": self._running,
            }


# ===========================================================================
# Dynamic Target Identifier (Section 10.3 Gap 2, Section 11.2.1)
# ===========================================================================

# Knowledge-graph gene list with druggability and therapeutic context.
# Replaces hardcoded 11-gene lookup with a configurable, extensible registry.
KNOWLEDGE_GRAPH_GENES: dict[str, dict[str, Any]] = {
    # Neurodegeneration / FTD
    "VCP": {
        "protein": "p97/VCP ATPase",
        "therapeutic_area": "neurodegeneration",
        "druggable": True,
        "reference_drug": "CB-5083",
        "pdb_ids": ["5FTK", "7K56", "8OOI", "9DIL"],
        "diseases": ["FTD", "ALS", "IBMPFD"],
    },
    "GRN": {
        "protein": "Progranulin",
        "therapeutic_area": "neurodegeneration",
        "druggable": True,
        "reference_drug": "Latozinemab",
        "pdb_ids": ["2JYE"],
        "diseases": ["FTD"],
    },
    "C9orf72": {
        "protein": "C9orf72",
        "therapeutic_area": "neurodegeneration",
        "druggable": False,
        "reference_drug": None,
        "pdb_ids": ["6LT0"],
        "diseases": ["FTD", "ALS"],
    },
    "MAPT": {
        "protein": "Tau",
        "therapeutic_area": "neurodegeneration",
        "druggable": True,
        "reference_drug": "Semorinemab",
        "pdb_ids": ["6QJH", "5O3L"],
        "diseases": ["FTD", "Alzheimer"],
    },

    # Oncology
    "BRCA1": {
        "protein": "BRCA1 DNA repair",
        "therapeutic_area": "oncology",
        "druggable": True,
        "reference_drug": "Olaparib",
        "pdb_ids": ["1JNX"],
        "diseases": ["Breast cancer", "Ovarian cancer"],
    },
    "BRCA2": {
        "protein": "BRCA2 DNA repair",
        "therapeutic_area": "oncology",
        "druggable": True,
        "reference_drug": "Olaparib",
        "pdb_ids": ["1MIU"],
        "diseases": ["Breast cancer", "Ovarian cancer"],
    },
    "TP53": {
        "protein": "p53 tumor suppressor",
        "therapeutic_area": "oncology",
        "druggable": True,
        "reference_drug": "APR-246",
        "pdb_ids": ["1TSR", "2XWR"],
        "diseases": ["Li-Fraumeni syndrome", "Multiple cancers"],
    },
    "EGFR": {
        "protein": "Epidermal Growth Factor Receptor",
        "therapeutic_area": "oncology",
        "druggable": True,
        "reference_drug": "Osimertinib",
        "pdb_ids": ["1M17", "4ZAU"],
        "diseases": ["NSCLC", "Glioblastoma"],
    },
    "KRAS": {
        "protein": "KRAS GTPase",
        "therapeutic_area": "oncology",
        "druggable": True,
        "reference_drug": "Sotorasib",
        "pdb_ids": ["4OBE", "6GJ8"],
        "diseases": ["NSCLC", "Pancreatic cancer", "Colorectal cancer"],
    },
    "BRAF": {
        "protein": "BRAF kinase",
        "therapeutic_area": "oncology",
        "druggable": True,
        "reference_drug": "Vemurafenib",
        "pdb_ids": ["1UWH"],
        "diseases": ["Melanoma", "Colorectal cancer"],
    },
    "ALK": {
        "protein": "ALK receptor tyrosine kinase",
        "therapeutic_area": "oncology",
        "druggable": True,
        "reference_drug": "Crizotinib",
        "pdb_ids": ["2XP2"],
        "diseases": ["NSCLC", "Neuroblastoma"],
    },
    "PIK3CA": {
        "protein": "PI3K catalytic subunit alpha",
        "therapeutic_area": "oncology",
        "druggable": True,
        "reference_drug": "Alpelisib",
        "pdb_ids": ["4JPS"],
        "diseases": ["Breast cancer"],
    },

    # Pharmacogenomics
    "CYP2D6": {
        "protein": "Cytochrome P450 2D6",
        "therapeutic_area": "pharmacogenomics",
        "druggable": False,
        "reference_drug": None,
        "pdb_ids": ["2F9Q"],
        "diseases": ["Drug metabolism"],
    },
    "CYP2C19": {
        "protein": "Cytochrome P450 2C19",
        "therapeutic_area": "pharmacogenomics",
        "druggable": False,
        "reference_drug": None,
        "pdb_ids": ["4GQS"],
        "diseases": ["Drug metabolism"],
    },

    # Metabolic
    "GLP1R": {
        "protein": "GLP-1 Receptor",
        "therapeutic_area": "metabolic",
        "druggable": True,
        "reference_drug": "Semaglutide",
        "pdb_ids": ["6B3J"],
        "diseases": ["Type 2 Diabetes", "Obesity"],
    },
    "PCSK9": {
        "protein": "PCSK9",
        "therapeutic_area": "metabolic",
        "druggable": True,
        "reference_drug": "Evolocumab",
        "pdb_ids": ["2P4E"],
        "diseases": ["Hypercholesterolemia"],
    },

    # CAR-T relevant
    "CD19": {
        "protein": "CD19 B-cell antigen",
        "therapeutic_area": "immunology",
        "druggable": True,
        "reference_drug": "Tisagenlecleucel",
        "pdb_ids": ["6AL5"],
        "diseases": ["B-ALL", "DLBCL"],
    },
    "BCMA": {
        "protein": "B-cell maturation antigen",
        "therapeutic_area": "immunology",
        "druggable": True,
        "reference_drug": "Idecabtagene vicleucel",
        "pdb_ids": ["1XU2"],
        "diseases": ["Multiple Myeloma"],
    },

    # Rare disease
    "CFTR": {
        "protein": "Cystic fibrosis transmembrane conductance regulator",
        "therapeutic_area": "rare_disease",
        "druggable": True,
        "reference_drug": "Elexacaftor/Tezacaftor/Ivacaftor",
        "pdb_ids": ["5UAK"],
        "diseases": ["Cystic Fibrosis"],
    },

    # Hematology
    "FLT3": {
        "protein": "FMS-like tyrosine kinase 3",
        "therapeutic_area": "oncology",
        "druggable": True,
        "reference_drug": "Midostaurin",
        "pdb_ids": ["1RJB"],
        "diseases": ["AML"],
    },
    "BCL2": {
        "protein": "B-cell lymphoma 2",
        "therapeutic_area": "oncology",
        "druggable": True,
        "reference_drug": "Venetoclax",
        "pdb_ids": ["1G5M"],
        "diseases": ["CLL", "AML"],
    },
}

# Clinical significance values considered actionable for target identification
ACTIONABLE_SIGNIFICANCE: set[str] = {
    "pathogenic",
    "likely_pathogenic",
    "likely pathogenic",
    "Pathogenic",
    "Likely_pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely_pathogenic",
    "pathogenic/likely_pathogenic",
}

# AlphaMissense class values considered actionable
ACTIONABLE_AM_CLASSES: set[str] = {
    "likely_pathogenic",
    "ambiguous",
}

# AlphaMissense score threshold for actionable variants
AM_PATHOGENICITY_THRESHOLD: float = 0.564


class DynamicTargetIdentifier:
    """
    Dynamic target identification from patient variant data.

    Replaces hardcoded gene lookups (Section 10.3 Gap 2) with a dynamic
    approach that queries the Milvus ``genomic_evidence`` collection for
    pathogenic/likely-pathogenic variants, groups them by gene, and
    cross-references against the knowledge graph to produce variant-driven
    target hypotheses.

    Parameters
    ----------
    milvus_client : Any, optional
        Milvus client for querying genomic_evidence.  Must have a
        ``query()`` or ``search_by_gene()`` method.
    knowledge_genes : dict, optional
        Override the knowledge graph gene registry.
    max_targets : int
        Maximum number of target hypotheses to return.
    min_variant_count : int
        Minimum pathogenic variant count in a gene to consider it.
    """

    def __init__(
        self,
        milvus_client: Any | None = None,
        knowledge_genes: dict[str, dict[str, Any]] | None = None,
        max_targets: int = 20,
        min_variant_count: int = 1,
    ):
        self._milvus = milvus_client
        self._knowledge = knowledge_genes or KNOWLEDGE_GRAPH_GENES
        self._max_targets = max_targets
        self._min_variant_count = min_variant_count

        logger.info(
            f"DynamicTargetIdentifier initialized: "
            f"{len(self._knowledge)} knowledge-graph genes, "
            f"max_targets={max_targets}"
        )

    def identify_targets(
        self,
        patient_id: str | None = None,
        filter_expr: str | None = None,
    ) -> list[dict[str, Any]]:
        """
        Identify drug targets dynamically from patient variant data.

        Workflow:
            1. Query genomic_evidence for pathogenic/likely-pathogenic variants
            2. Group variants by gene
            3. Cross-reference each gene with the knowledge graph
            4. Score and rank targets by variant evidence + druggability
            5. Return top-N target hypotheses

        Parameters
        ----------
        patient_id : str, optional
            Patient identifier for audit trail.
        filter_expr : str, optional
            Additional Milvus filter expression.

        Returns
        -------
        list[dict]
            Ranked list of target hypotheses with variant evidence.
        """
        logger.info(
            f"Dynamic target identification started"
            + (f" for patient {patient_id}" if patient_id else "")
        )

        # Step 1: Query for pathogenic variants
        pathogenic_variants = self._query_pathogenic_variants(filter_expr)

        if not pathogenic_variants:
            logger.warning("No pathogenic variants found in genomic_evidence")
            return []

        logger.info(f"Found {len(pathogenic_variants)} actionable variants")

        # Step 2: Group by gene
        gene_groups = self._group_by_gene(pathogenic_variants)
        logger.info(
            f"Variants span {len(gene_groups)} genes: "
            f"{list(gene_groups.keys())[:10]}..."
        )

        # Step 3: Cross-reference with knowledge graph and score
        targets = self._score_and_rank(gene_groups)

        # Step 4: Return top-N
        targets = targets[: self._max_targets]
        logger.info(
            f"Identified {len(targets)} target hypotheses"
            + (f" for patient {patient_id}" if patient_id else "")
        )
        return targets

    def identify_targets_for_gene(self, gene: str) -> dict[str, Any] | None:
        """
        Identify target potential for a specific gene.

        Combines variant evidence from Milvus with knowledge-graph context.
        """
        gene_upper = gene.upper()

        # Get knowledge graph entry
        kg_entry = self._knowledge.get(gene_upper)

        # Query variants for this gene
        variants = self._query_gene_variants(gene_upper)

        if not variants and not kg_entry:
            return None

        pathogenic_count = sum(
            1 for v in variants
            if self._is_actionable_variant(v)
        )

        target = {
            "gene": gene_upper,
            "variant_count": len(variants),
            "pathogenic_count": pathogenic_count,
            "variants": variants[:10],  # top 10 variants
            "in_knowledge_graph": kg_entry is not None,
        }

        if kg_entry:
            target.update({
                "protein": kg_entry.get("protein"),
                "therapeutic_area": kg_entry.get("therapeutic_area"),
                "druggable": kg_entry.get("druggable", False),
                "reference_drug": kg_entry.get("reference_drug"),
                "pdb_ids": kg_entry.get("pdb_ids", []),
                "diseases": kg_entry.get("diseases", []),
            })

        # Compute a priority score
        target["priority_score"] = self._compute_priority(
            pathogenic_count=pathogenic_count,
            druggable=kg_entry.get("druggable", False) if kg_entry else False,
            has_pdb=bool(kg_entry.get("pdb_ids")) if kg_entry else False,
        )

        return target

    # ---- private helpers ---------------------------------------------------

    def _query_pathogenic_variants(
        self,
        filter_expr: str | None = None,
    ) -> list[dict[str, Any]]:
        """
        Query Milvus genomic_evidence for pathogenic/likely-pathogenic variants.
        """
        if self._milvus is None:
            logger.warning("No Milvus client -- returning empty variant list")
            return []

        results: list[dict[str, Any]] = []

        # Strategy 1: Query by clinical significance
        significance_values = [
            "pathogenic", "likely_pathogenic", "Pathogenic",
            "Likely_pathogenic", "Pathogenic/Likely_pathogenic",
        ]

        for sig in significance_values:
            try:
                expr = f'clinical_significance == "{sig}"'
                if filter_expr:
                    expr = f"({expr}) and ({filter_expr})"

                if hasattr(self._milvus, "get_collection"):
                    collection = self._milvus.get_collection()
                    batch = collection.query(
                        expr=expr,
                        output_fields=[
                            "id", "gene", "chrom", "pos", "ref", "alt",
                            "consequence", "impact", "clinical_significance",
                            "rsid", "text_summary", "am_pathogenicity",
                            "am_class",
                        ],
                        limit=500,
                    )
                    results.extend(batch)
                elif hasattr(self._milvus, "query_collection"):
                    batch = self._milvus.query_collection(
                        collection_name="genomic_evidence",
                        expr=expr,
                        limit=500,
                    )
                    results.extend(batch)
            except Exception as exc:
                logger.debug(
                    f"Query for significance='{sig}' failed: {exc}"
                )

        # Strategy 2: Query by AlphaMissense pathogenicity score
        try:
            am_expr = f"am_pathogenicity >= {AM_PATHOGENICITY_THRESHOLD}"
            if filter_expr:
                am_expr = f"({am_expr}) and ({filter_expr})"

            if hasattr(self._milvus, "get_collection"):
                collection = self._milvus.get_collection()
                am_results = collection.query(
                    expr=am_expr,
                    output_fields=[
                        "id", "gene", "chrom", "pos", "ref", "alt",
                        "consequence", "impact", "clinical_significance",
                        "rsid", "text_summary", "am_pathogenicity",
                        "am_class",
                    ],
                    limit=500,
                )
                # Merge, avoiding duplicates
                existing_ids = {r.get("id") for r in results}
                for r in am_results:
                    if r.get("id") not in existing_ids:
                        results.append(r)
                        existing_ids.add(r.get("id"))
        except Exception as exc:
            logger.debug(f"AlphaMissense query failed: {exc}")

        return results

    def _query_gene_variants(self, gene: str) -> list[dict[str, Any]]:
        """Query all variants for a specific gene."""
        if self._milvus is None:
            return []

        try:
            if hasattr(self._milvus, "search_by_gene"):
                return self._milvus.search_by_gene(gene, top_k=100)
            elif hasattr(self._milvus, "get_collection"):
                collection = self._milvus.get_collection()
                return collection.query(
                    expr=f'gene == "{gene}"',
                    output_fields=[
                        "id", "gene", "chrom", "pos", "ref", "alt",
                        "consequence", "impact", "clinical_significance",
                        "rsid", "text_summary", "am_pathogenicity",
                        "am_class",
                    ],
                    limit=100,
                )
        except Exception as exc:
            logger.warning(f"Failed to query variants for gene {gene}: {exc}")

        return []

    @staticmethod
    def _is_actionable_variant(variant: dict[str, Any]) -> bool:
        """Check if a variant is clinically actionable."""
        # Check clinical significance
        clin_sig = variant.get("clinical_significance", "")
        if clin_sig and clin_sig.lower().replace(" ", "_") in {
            s.lower().replace(" ", "_") for s in ACTIONABLE_SIGNIFICANCE
        }:
            return True

        # Check AlphaMissense
        am_score = variant.get("am_pathogenicity")
        if am_score is not None and am_score >= AM_PATHOGENICITY_THRESHOLD:
            return True

        am_class = variant.get("am_class", "")
        if am_class in ACTIONABLE_AM_CLASSES:
            return True

        # Check impact
        impact = variant.get("impact", "")
        if impact and impact.upper() == "HIGH":
            return True

        return False

    @staticmethod
    def _group_by_gene(
        variants: list[dict[str, Any]],
    ) -> dict[str, list[dict[str, Any]]]:
        """Group variants by gene symbol."""
        groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for v in variants:
            gene = v.get("gene", "")
            if gene:
                groups[gene.upper()].append(v)
        return dict(groups)

    def _score_and_rank(
        self,
        gene_groups: dict[str, list[dict[str, Any]]],
    ) -> list[dict[str, Any]]:
        """
        Score each gene as a potential target and return ranked list.

        Scoring factors:
            - Number of pathogenic variants (higher = better)
            - Presence in knowledge graph
            - Druggability
            - Availability of PDB structures
            - Variant impact severity
        """
        targets: list[dict[str, Any]] = []

        for gene, variants in gene_groups.items():
            if len(variants) < self._min_variant_count:
                continue

            pathogenic_count = sum(
                1 for v in variants if self._is_actionable_variant(v)
            )

            kg_entry = self._knowledge.get(gene)

            target: dict[str, Any] = {
                "gene": gene,
                "variant_count": len(variants),
                "pathogenic_count": pathogenic_count,
                "variants": sorted(
                    variants,
                    key=lambda v: v.get("am_pathogenicity", 0) or 0,
                    reverse=True,
                )[:10],
                "in_knowledge_graph": kg_entry is not None,
            }

            if kg_entry:
                target.update({
                    "protein": kg_entry.get("protein"),
                    "therapeutic_area": kg_entry.get("therapeutic_area"),
                    "druggable": kg_entry.get("druggable", False),
                    "reference_drug": kg_entry.get("reference_drug"),
                    "pdb_ids": kg_entry.get("pdb_ids", []),
                    "diseases": kg_entry.get("diseases", []),
                })

            priority = self._compute_priority(
                pathogenic_count=pathogenic_count,
                druggable=kg_entry.get("druggable", False) if kg_entry else False,
                has_pdb=bool(kg_entry.get("pdb_ids")) if kg_entry else False,
            )
            target["priority_score"] = priority

            targets.append(target)

        # Sort by priority descending
        targets.sort(key=lambda t: t["priority_score"], reverse=True)
        return targets

    @staticmethod
    def _compute_priority(
        pathogenic_count: int,
        druggable: bool,
        has_pdb: bool,
    ) -> float:
        """
        Compute a priority score for a candidate target.

        Score range: 0.0 -- 10.0
        """
        score = 0.0

        # Variant evidence (0--4 points)
        score += min(pathogenic_count, 10) * 0.4

        # Druggability (0 or 3 points)
        if druggable:
            score += 3.0

        # Structural data available (0 or 2 points)
        if has_pdb:
            score += 2.0

        # Knowledge graph membership (1 point implicit via druggable/pdb)
        return round(score, 2)


# ===========================================================================
# Built-in trigger rules
# ===========================================================================

def _condition_vcf_ready(event: PipelineEvent) -> bool:
    """Condition: VCF file path present in payload."""
    return bool(event.payload.get("vcf_path"))


def _action_vcf_to_rag_ingest(event: PipelineEvent) -> None:
    """
    Action: VCF Ready -> trigger RAG ingest (annotation + Milvus embedding).

    Publishes an ANNOTATION_COMPLETE event that downstream handlers can
    use to begin the embedding pipeline.
    """
    vcf_path = event.payload.get("vcf_path", "")
    patient_id = event.patient_id

    logger.info(
        f"[Trigger] VCF -> RAG Ingest: starting annotation for "
        f"patient={patient_id}, vcf={vcf_path}"
    )

    # Publish downstream event for the RAG ingest pipeline
    bus = get_event_bus()
    bus.publish(PipelineEvent(
        event_type=EventType.ANNOTATION_COMPLETE,
        source_stage=PipelineStage.ANNOTATION,
        target_stage=PipelineStage.RAG_INGEST,
        patient_id=patient_id,
        payload={
            "vcf_path": vcf_path,
            "triggered_by": "vcf_to_rag_ingest",
            "auto_triggered": True,
        },
        priority=EventPriority.HIGH,
    ))


def _condition_pathogenic_variant(event: PipelineEvent) -> bool:
    """Condition: payload contains identified targets with pathogenic variants."""
    targets = event.payload.get("targets", [])
    if not targets:
        return False
    # At least one target must be druggable
    return any(t.get("druggable", False) for t in targets)


def _action_rag_to_drug_discovery(event: PipelineEvent) -> None:
    """
    Action: Pathogenic variant identified -> trigger drug discovery pipeline.

    Selects the highest-priority druggable target and publishes a
    DRUG_DESIGN_REQUEST event.
    """
    targets = event.payload.get("targets", [])
    patient_id = event.patient_id

    # Find best druggable target
    druggable = [t for t in targets if t.get("druggable", False)]
    if not druggable:
        return

    best_target = max(druggable, key=lambda t: t.get("priority_score", 0))

    logger.info(
        f"[Trigger] RAG -> Drug Discovery: targeting {best_target['gene']} "
        f"(score={best_target.get('priority_score', 0)}) "
        f"for patient={patient_id}"
    )

    bus = get_event_bus()
    bus.publish(PipelineEvent(
        event_type=EventType.DRUG_DESIGN_REQUEST,
        source_stage=PipelineStage.RAG_CHAT,
        target_stage=PipelineStage.DRUG_DISCOVERY,
        patient_id=patient_id,
        payload={
            "gene": best_target["gene"],
            "protein": best_target.get("protein"),
            "pdb_ids": best_target.get("pdb_ids", []),
            "reference_drug": best_target.get("reference_drug"),
            "therapeutic_area": best_target.get("therapeutic_area"),
            "variant_count": best_target.get("pathogenic_count", 0),
            "triggered_by": "rag_to_drug_discovery",
            "auto_triggered": True,
        },
        priority=EventPriority.HIGH,
    ))


def _condition_cart_query(event: PipelineEvent) -> bool:
    """Condition: deep analysis request targeting CAR-T domain."""
    payload = event.payload
    query_text = payload.get("query", "").lower()

    cart_keywords = [
        "car-t", "car t", "cart", "cell therapy", "cell-therapy",
        "chimeric antigen", "adoptive cell", "cd19", "cd22", "bcma",
    ]
    return any(kw in query_text for kw in cart_keywords)


def _action_cart_cross_query(event: PipelineEvent) -> None:
    """
    Action: CAR-T query -> publish VARIANT_CONTEXT_REQUEST so the CAR-T
    agent can access patient variant data from genomic_evidence.
    """
    patient_id = event.patient_id
    query = event.payload.get("query", "")

    logger.info(
        f"[Trigger] RAG <-> CAR-T: cross-query for patient={patient_id}, "
        f"query='{query[:80]}...'"
    )

    bus = get_event_bus()
    bus.publish(PipelineEvent(
        event_type=EventType.VARIANT_CONTEXT_REQUEST,
        source_stage=PipelineStage.CART_ANALYSIS,
        target_stage=PipelineStage.RAG_CHAT,
        patient_id=patient_id,
        payload={
            "query": query,
            "requested_collections": ["genomic_evidence"],
            "triggered_by": "cart_cross_query",
            "auto_triggered": True,
        },
        priority=EventPriority.NORMAL,
    ))


def _condition_candidates_ready(event: PipelineEvent) -> bool:
    """Condition: drug candidates present with docking scores."""
    candidates = event.payload.get("candidates", [])
    return len(candidates) >= 1


def _action_drug_to_report(event: PipelineEvent) -> None:
    """
    Action: Drug candidates ranked -> trigger report generation.
    """
    candidates = event.payload.get("candidates", [])
    patient_id = event.patient_id
    gene = event.payload.get("gene", "unknown")

    logger.info(
        f"[Trigger] Drug Discovery -> Report: {len(candidates)} candidates "
        f"for {gene}, patient={patient_id}"
    )

    bus = get_event_bus()
    bus.publish(PipelineEvent(
        event_type=EventType.REPORT_GENERATED,
        source_stage=PipelineStage.DRUG_DISCOVERY,
        target_stage=PipelineStage.REPORTING,
        patient_id=patient_id,
        payload={
            "gene": gene,
            "candidate_count": len(candidates),
            "top_candidate": candidates[0] if candidates else None,
            "triggered_by": "drug_to_report",
            "auto_triggered": True,
        },
        priority=EventPriority.NORMAL,
    ))


# ===========================================================================
# Built-in trigger definitions
# ===========================================================================

BUILTIN_TRIGGERS: list[TriggerRule] = [
    TriggerRule(
        name="vcf_to_rag_ingest",
        event_types={EventType.VCF_READY},
        condition=_condition_vcf_ready,
        action=_action_vcf_to_rag_ingest,
        cooldown_seconds=60.0,
        description=(
            "VCF Ready -> auto-trigger variant annotation and Milvus embedding. "
            "When Parabricks produces a VCF, this trigger initiates the RAG "
            "ingest pipeline to annotate variants with ClinVar/AlphaMissense "
            "and embed them into the genomic_evidence collection."
        ),
    ),
    TriggerRule(
        name="rag_to_drug_discovery",
        event_types={EventType.TARGETS_IDENTIFIED},
        condition=_condition_pathogenic_variant,
        action=_action_rag_to_drug_discovery,
        cooldown_seconds=30.0,
        description=(
            "Pathogenic variant identified -> auto-trigger target hypothesis "
            "generation and MolMIM molecule generation.  Selects the highest-"
            "priority druggable target and initiates the drug discovery pipeline."
        ),
    ),
    TriggerRule(
        name="cart_cross_query",
        event_types={EventType.DEEP_ANALYSIS_REQUEST},
        condition=_condition_cart_query,
        action=_action_cart_cross_query,
        cooldown_seconds=5.0,
        description=(
            "RAG <-> CAR-T bidirectional cross-query.  When a user asks about "
            "cell therapy, the CAR-T agent queries genomic_evidence for patient "
            "variants.  When CAR-T identifies targets, it queries RAG for "
            "variant context."
        ),
    ),
    TriggerRule(
        name="drug_to_report",
        event_types={EventType.DRUG_CANDIDATES_READY},
        condition=_condition_candidates_ready,
        action=_action_drug_to_report,
        cooldown_seconds=30.0,
        description=(
            "Drug Discovery -> Report.  When top candidates are ranked with "
            "docking scores, auto-trigger report generation."
        ),
    ),
]


def register_builtin_triggers(
    manager: BidirectionalTriggerManager,
) -> list[str]:
    """
    Register all built-in trigger rules with a TriggerManager.

    Returns a list of rule_ids.
    """
    rule_ids: list[str] = []
    for rule in BUILTIN_TRIGGERS:
        rid = manager.register_rule(rule)
        rule_ids.append(rid)
    logger.info(
        f"Registered {len(rule_ids)} built-in trigger rules: "
        f"{[r.name for r in BUILTIN_TRIGGERS]}"
    )
    return rule_ids


# ---------------------------------------------------------------------------
# Module-level convenience
# ---------------------------------------------------------------------------

def create_trigger_manager(
    event_bus: EventBus | None = None,
    register_builtins: bool = True,
    auto_start: bool = True,
    **kwargs: Any,
) -> BidirectionalTriggerManager:
    """
    Factory function to create and configure a BidirectionalTriggerManager.

    Parameters
    ----------
    event_bus : EventBus, optional
        Shared event bus.  Default: singleton.
    register_builtins : bool
        Whether to register the 4 built-in trigger rules.
    auto_start : bool
        Whether to call ``start()`` immediately.
    **kwargs
        Forwarded to BidirectionalTriggerManager constructor.
    """
    manager = BidirectionalTriggerManager(event_bus=event_bus, **kwargs)

    if register_builtins:
        register_builtin_triggers(manager)

    if auto_start:
        manager.start()

    return manager


def create_dynamic_target_identifier(
    milvus_client: Any | None = None,
    **kwargs: Any,
) -> DynamicTargetIdentifier:
    """
    Factory function to create a DynamicTargetIdentifier.
    """
    return DynamicTargetIdentifier(milvus_client=milvus_client, **kwargs)
