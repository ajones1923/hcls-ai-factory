"""
hcls_common -- Shared library for the HCLS AI Factory precision medicine platform.

Consolidates duplicate code across the RAG Chat pipeline, Drug Discovery pipeline,
and CAR-T Intelligence Agent. Part of the end-to-end Patient DNA -> Drug Candidates
workflow running on a single NVIDIA DGX Spark.

Modules:
    config          Centralized Pydantic settings (Milvus, embedding, LLM, NIM, observability)
    milvus_client   Unified Milvus vector-database wrapper with connection pooling
    embedder        Embedding providers (local BGE, TEI) with caching and Prometheus metrics
    llm_client      Multi-provider LLM client with retry, budget tracking, response cache
    circuit_breaker General-purpose circuit breaker (sync/async, decorator, context manager)
    tracing         OpenTelemetry integration with graceful degradation
    security        Input validation, sanitization, rate limiting for genomic data
    enums           Domain enumerations for the precision medicine pipeline
    query_router    Cross-collection query routing (genomic_evidence + 11 CAR-T collections)
    event_bus       Event-driven pipeline coordination with file-based audit trail
    bidirectional_triggers  Bidirectional triggers and dynamic target identification
    meta_agent      Meta-agent orchestration using Claude tool-use for clinical intelligence
    report_generator  Unified clinical report generation (Markdown, JSON, PDF)
    knowledge_version  Semantic versioning for all agent knowledge bases
    reproducibility    Reproducibility manifest for FDA 21 CFR Part 11 compliance
    event_monitor      Live event bus monitor for cross-agent intelligence visualization
"""

__version__ = "0.1.0"
__author__ = "Adam Jones"
__license__ = "Apache-2.0"

# ---------------------------------------------------------------------------
# Lazy re-exports -- import the most commonly used symbols so downstream code
# can do ``from hcls_common import get_settings, UnifiedMilvusClient, ...``
# ---------------------------------------------------------------------------

from hcls_common.config import HCLSSettings, get_settings
from hcls_common.milvus_client import UnifiedMilvusClient
from hcls_common.embedder import (
    BaseEmbedder,
    LocalEmbedder,
    TEIEmbedder,
    get_embedder,
)
from hcls_common.llm_client import (
    BaseLLMClient,
    AnthropicClient,
    OpenAIClient,
    OllamaClient,
    VLLMClient,
    LLMClientFactory,
    get_llm_client,
)
from hcls_common.circuit_breaker import (
    CircuitBreaker,
    CircuitBreakerOpen,
    CircuitState,
)
from hcls_common.tracing import init_tracing, traced
from hcls_common.security import (
    sanitize_search_query,
    validate_milvus_filter,
    sanitize_gene_name,
    sanitize_chromosome,
    validate_patient_id,
    validate_smiles,
    validate_pdb_id,
    rate_limit_check,
    add_security_headers,
    sanitize_filename,
)
from hcls_common.enums import (
    PipelineStage,
    ClinicalSignificance,
    VariantImpact,
    TherapeuticArea,
    DrugLikenessFilter,
    DockingMethod,
    LLMProvider,
    EmbeddingProvider,
    NIMService,
    AlertSeverity,
)
from hcls_common.query_router import (
    QueryRouter,
    QueryPlan,
    QueryIntent,
    IntentClassification,
    MergedResult,
    classify_intent,
    create_query_router,
)
from hcls_common.event_bus import (
    EventBus,
    EventType,
    EventStatus,
    EventPriority,
    PipelineEvent,
    PatientCase,
    get_event_bus,
    publish_event,
    create_pipeline_event,
)
from hcls_common.bidirectional_triggers import (
    BidirectionalTriggerManager,
    TriggerRule,
    DynamicTargetIdentifier,
    register_builtin_triggers,
    create_trigger_manager,
    create_dynamic_target_identifier,
)
from hcls_common.meta_agent import (
    MetaAgent,
    AgentResponse,
    ToolCallRecord,
    ToolName,
    create_meta_agent,
)
from hcls_common.report_generator import (
    ClinicalReportGenerator,
    ReportData,
    VariantRecord,
    TargetRecord,
    StructureRecord,
    DrugCandidateRecord,
    CARTEvaluationRecord,
    create_report_generator,
)
from hcls_common.knowledge_version import (
    KnowledgeManifest,
    KnowledgeSource,
    BIOMARKER_KNOWLEDGE,
    ONCOLOGY_KNOWLEDGE,
    CART_KNOWLEDGE,
    IMAGING_KNOWLEDGE,
    AUTOIMMUNE_KNOWLEDGE,
)
from hcls_common.reproducibility import (
    ReproducibilityManifest,
    HardwareSpec,
)
from hcls_common.event_monitor import (
    EventRecord,
    read_event_log,
    get_event_stats,
    render_event_flow_ascii,
)
from hcls_common.demo_data import (
    DEMO_PATIENT_ID,
    DEMO_PATIENT_AGE,
    DEMO_PATIENT_SEX,
    DEMO_BIOMARKERS,
    DEMO_GENOTYPES,
    DEMO_STAR_ALLELES,
    DEMO_ONCOLOGY,
    DEMO_CART,
    DEMO_IMAGING,
    DEMO_AUTOIMMUNE,
    get_demo_patient_summary,
)
from hcls_common.verify_gate import (
    honesty_check,
    verify_claims,
    verify_text,
    is_publishable,
)
from hcls_common.capability_registry import (
    Capability,
    CapabilityRegistry,
    CapabilityType,
    Port,
    Serving,
    Status,
    ValueShape,
    ArtifactShape,
    get_registry,
    validate_inputs,
    inputs_ok,
)
from hcls_common.mcp_server import FactoryTools, build_server
from hcls_common.workflow_composer import (
    WorkflowComposer, Pipeline, Node, NodeInput,
)
from hcls_common.mlops import MLOpsStore, get_mlops_store
from hcls_common.governance import second_opinion, dual_predict
from hcls_common.external_tools import ingest_external_tools, external_tool_to_capability
from hcls_common.license_gate import audit as license_audit, classify_license, check_requirements
from hcls_common.deploy_wizard import teardown_order, deploy_order
from hcls_common.multiomics import (
    MultiOmicsStore, MultiOmicsRecord, PatientContext, PatientContextStore,
)
from hcls_common.biokey import (
    BioKey, BioKeyResolver, EntityKind, resolve, DEFAULT_RESOLVER,
)
from hcls_common.artifact import (
    Artifact,
    Honesty,
    Maturity,
    Provenance,
    new_artifact,
    combine_honesty,
    weakest_maturity,
    non_inflation_issues,
    derive_artifact,
)

__all__ = [
    # capability_registry (A1)
    "Capability",
    "CapabilityRegistry",
    "CapabilityType",
    "Port",
    "Serving",
    "Status",
    "ValueShape",
    "ArtifactShape",
    "get_registry",
    "validate_inputs",
    "inputs_ok",
    # mcp tool-surface (A2)
    "FactoryTools",
    "build_server",
    # workflow composer (G1)
    "WorkflowComposer",
    "Pipeline",
    "Node",
    "NodeInput",
    # mlops (A3)
    "MLOpsStore",
    "get_mlops_store",
    # governance (A5)
    "second_opinion",
    "dual_predict",
    # external tools (A3)
    "ingest_external_tools",
    "external_tool_to_capability",
    # A7/A8
    "license_audit",
    "classify_license",
    "check_requirements",
    "teardown_order",
    "deploy_order",
    # F1 multi-omics
    "MultiOmicsStore",
    "MultiOmicsRecord",
    # patient context (PF-1)
    "PatientContext",
    "PatientContextStore",
    # biokey resolver (PF-5)
    "BioKey",
    "BioKeyResolver",
    "EntityKind",
    "resolve",
    "DEFAULT_RESOLVER",
    # artifact envelope (PF-3)
    "Artifact",
    "Honesty",
    "Maturity",
    "Provenance",
    "new_artifact",
    "combine_honesty",
    "weakest_maturity",
    "non_inflation_issues",
    "derive_artifact",
    # verify_gate (ORCH-11)
    "honesty_check",
    "verify_claims",
    "verify_text",
    "is_publishable",
    # config
    "HCLSSettings",
    "get_settings",
    # milvus_client
    "UnifiedMilvusClient",
    # embedder
    "BaseEmbedder",
    "LocalEmbedder",
    "TEIEmbedder",
    "get_embedder",
    # llm_client
    "BaseLLMClient",
    "AnthropicClient",
    "OpenAIClient",
    "OllamaClient",
    "VLLMClient",
    "LLMClientFactory",
    "get_llm_client",
    # circuit_breaker
    "CircuitBreaker",
    "CircuitBreakerOpen",
    "CircuitState",
    # tracing
    "init_tracing",
    "traced",
    # security
    "sanitize_search_query",
    "validate_milvus_filter",
    "sanitize_gene_name",
    "sanitize_chromosome",
    "validate_patient_id",
    "validate_smiles",
    "validate_pdb_id",
    "rate_limit_check",
    "add_security_headers",
    "sanitize_filename",
    # enums
    "PipelineStage",
    "ClinicalSignificance",
    "VariantImpact",
    "TherapeuticArea",
    "DrugLikenessFilter",
    "DockingMethod",
    "LLMProvider",
    "EmbeddingProvider",
    "NIMService",
    "AlertSeverity",
    # query_router
    "QueryRouter",
    "QueryPlan",
    "QueryIntent",
    "IntentClassification",
    "MergedResult",
    "classify_intent",
    "create_query_router",
    # event_bus
    "EventBus",
    "EventType",
    "EventStatus",
    "EventPriority",
    "PipelineEvent",
    "PatientCase",
    "get_event_bus",
    "publish_event",
    "create_pipeline_event",
    # bidirectional_triggers
    "BidirectionalTriggerManager",
    "TriggerRule",
    "DynamicTargetIdentifier",
    "register_builtin_triggers",
    "create_trigger_manager",
    "create_dynamic_target_identifier",
    # meta_agent
    "MetaAgent",
    "AgentResponse",
    "ToolCallRecord",
    "ToolName",
    "create_meta_agent",
    # report_generator
    "ClinicalReportGenerator",
    "ReportData",
    "VariantRecord",
    "TargetRecord",
    "StructureRecord",
    "DrugCandidateRecord",
    "CARTEvaluationRecord",
    "create_report_generator",
    # knowledge_version
    "KnowledgeManifest",
    "KnowledgeSource",
    "BIOMARKER_KNOWLEDGE",
    "ONCOLOGY_KNOWLEDGE",
    "CART_KNOWLEDGE",
    "IMAGING_KNOWLEDGE",
    "AUTOIMMUNE_KNOWLEDGE",
    # reproducibility
    "ReproducibilityManifest",
    "HardwareSpec",
    # event_monitor
    "EventRecord",
    "read_event_log",
    "get_event_stats",
    "render_event_flow_ascii",
    # demo_data
    "DEMO_PATIENT_ID",
    "DEMO_PATIENT_AGE",
    "DEMO_PATIENT_SEX",
    "DEMO_BIOMARKERS",
    "DEMO_GENOTYPES",
    "DEMO_STAR_ALLELES",
    "DEMO_ONCOLOGY",
    "DEMO_CART",
    "DEMO_IMAGING",
    "DEMO_AUTOIMMUNE",
    "get_demo_patient_summary",
]
