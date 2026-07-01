"""Tests for hcls_common.event_bus — PipelineEvent, EventBus pub/sub."""

import pytest

from hcls_common.event_bus import (
    EventBus,
    EventPriority,
    EventStatus,
    EventType,
    PatientCase,
    PipelineEvent,
    PipelineStage,
    Subscription,
)


class TestEventType:
    def test_all_values_are_strings(self):
        for et in EventType:
            assert isinstance(et.value, str)

    def test_display_name(self):
        assert EventType.VCF_READY.display_name == "Vcf Ready"
        assert EventType.FASTQ_ARRIVED.display_name == "Fastq Arrived"


class TestEventStatus:
    def test_terminal_states(self):
        assert EventStatus.COMPLETED.is_terminal is True
        assert EventStatus.FAILED.is_terminal is True
        assert EventStatus.EXPIRED.is_terminal is True

    def test_non_terminal_states(self):
        assert EventStatus.PENDING.is_terminal is False
        assert EventStatus.DISPATCHED.is_terminal is False
        assert EventStatus.PROCESSING.is_terminal is False


class TestEventPriority:
    def test_ordering(self):
        assert EventPriority.CRITICAL.value < EventPriority.HIGH.value
        assert EventPriority.HIGH.value < EventPriority.NORMAL.value
        assert EventPriority.NORMAL.value < EventPriority.LOW.value
        assert EventPriority.LOW.value < EventPriority.BACKGROUND.value


class TestPipelineEvent:
    def test_creation_with_defaults(self):
        event = PipelineEvent(
            event_type=EventType.VCF_READY,
            source_stage=PipelineStage.GENOMICS,
        )
        assert event.status == EventStatus.PENDING
        assert event.priority == EventPriority.NORMAL
        assert event.event_id  # should be a UUID string
        assert event.created_at  # should be an ISO timestamp

    def test_to_dict_roundtrip(self):
        event = PipelineEvent(
            event_type=EventType.VCF_READY,
            source_stage=PipelineStage.GENOMICS,
            patient_id="HG002",
            payload={"vcf_path": "/data/HG002.vcf.gz"},
        )
        d = event.to_dict()
        assert d["event_type"] == "vcf_ready"
        assert d["source_stage"] == "genomics"
        assert d["patient_id"] == "HG002"

        restored = PipelineEvent.from_dict(d)
        assert restored.event_type == EventType.VCF_READY
        assert restored.patient_id == "HG002"

    def test_mark_dispatched(self):
        event = PipelineEvent(
            event_type=EventType.FASTQ_ARRIVED,
            source_stage=PipelineStage.GENOMICS,
        )
        event.mark_dispatched()
        assert event.status == EventStatus.DISPATCHED
        assert event.dispatched_at is not None

    def test_mark_completed(self):
        event = PipelineEvent(
            event_type=EventType.FASTQ_ARRIVED,
            source_stage=PipelineStage.GENOMICS,
        )
        event.mark_completed()
        assert event.status == EventStatus.COMPLETED
        assert event.completed_at is not None

    def test_mark_failed(self):
        event = PipelineEvent(
            event_type=EventType.FASTQ_ARRIVED,
            source_stage=PipelineStage.GENOMICS,
        )
        event.mark_failed("Something went wrong")
        assert event.status == EventStatus.FAILED
        assert event.error == "Something went wrong"

    def test_priority_comparison(self):
        critical = PipelineEvent(
            event_type=EventType.VCF_READY,
            source_stage=PipelineStage.GENOMICS,
            priority=EventPriority.CRITICAL,
        )
        low = PipelineEvent(
            event_type=EventType.VCF_READY,
            source_stage=PipelineStage.GENOMICS,
            priority=EventPriority.LOW,
        )
        assert critical < low


class TestSubscription:
    def test_matches_correct_event_type(self):
        sub = Subscription(
            handler=lambda e: None,
            event_types={EventType.VCF_READY},
        )
        event = PipelineEvent(
            event_type=EventType.VCF_READY,
            source_stage=PipelineStage.GENOMICS,
        )
        assert sub.matches(event) is True

    def test_rejects_wrong_event_type(self):
        sub = Subscription(
            handler=lambda e: None,
            event_types={EventType.VCF_READY},
        )
        event = PipelineEvent(
            event_type=EventType.FASTQ_ARRIVED,
            source_stage=PipelineStage.GENOMICS,
        )
        assert sub.matches(event) is False

    def test_source_filter(self):
        sub = Subscription(
            handler=lambda e: None,
            event_types={EventType.VCF_READY},
            source_filter=PipelineStage.GENOMICS,
        )
        wrong_source = PipelineEvent(
            event_type=EventType.VCF_READY,
            source_stage=PipelineStage.RAG_CHAT,
        )
        assert sub.matches(wrong_source) is False


class TestPatientCase:
    def test_creation(self):
        case = PatientCase(patient_id="HG002")
        assert case.patient_id == "HG002"
        assert case.case_id  # UUID
        assert case.variants == []
        assert case.drug_candidates == []

    def test_add_processing_step(self):
        case = PatientCase(patient_id="HG002")
        case.add_processing_step("genomics", "vcf_generated", {"path": "/data/HG002.vcf"})
        assert len(case.processing_log) == 1
        assert case.processing_log[0]["agent"] == "genomics"
        assert case.processing_log[0]["action"] == "vcf_generated"

    def test_to_dict_roundtrip(self):
        case = PatientCase(patient_id="HG002", age=45, sex="male")
        d = case.to_dict()
        assert d["patient_id"] == "HG002"
        assert d["age"] == 45

        restored = PatientCase.from_dict(d)
        assert restored.patient_id == "HG002"
        assert restored.age == 45


class TestEventBus:
    def setup_method(self):
        EventBus.reset_instance()

    def teardown_method(self):
        EventBus.reset_instance()

    def test_singleton(self):
        bus1 = EventBus.get_instance(enable_audit=False)
        bus2 = EventBus.get_instance(enable_audit=False)
        assert bus1 is bus2

    def test_publish_and_subscribe(self):
        bus = EventBus.get_instance(enable_audit=False)
        received = []

        def handler(event):
            received.append(event)

        bus.subscribe(
            handler=handler,
            event_types={EventType.VCF_READY},
            name="test_handler",
        )

        event = PipelineEvent(
            event_type=EventType.VCF_READY,
            source_stage=PipelineStage.GENOMICS,
            patient_id="HG002",
        )
        bus.publish(event)

        assert len(received) == 1
        assert received[0].patient_id == "HG002"

    def test_unsubscribe(self):
        bus = EventBus.get_instance(enable_audit=False)
        received = []

        bus.subscribe(
            handler=lambda e: received.append(e),
            event_types={EventType.VCF_READY},
            name="removable",
        )
        assert bus.unsubscribe("removable") is True
        assert bus.unsubscribe("nonexistent") is False

    def test_get_stats(self):
        bus = EventBus.get_instance(enable_audit=False)
        stats = bus.get_stats()
        assert "published" in stats
        assert "queue_size" in stats
        assert stats["running"] is True

    def test_handler_error_does_not_crash_bus(self):
        bus = EventBus.get_instance(enable_audit=False)
        received_ok = []

        def bad_handler(event):
            raise RuntimeError("handler exploded")

        def good_handler(event):
            received_ok.append(event)

        bus.subscribe(handler=bad_handler, event_types={EventType.VCF_READY}, name="bad")
        bus.subscribe(handler=good_handler, event_types={EventType.VCF_READY}, name="good")

        event = PipelineEvent(
            event_type=EventType.VCF_READY,
            source_stage=PipelineStage.GENOMICS,
        )
        bus.publish(event)

        # Good handler still received the event despite bad handler crashing
        assert len(received_ok) == 1

    def test_get_history(self):
        bus = EventBus.get_instance(enable_audit=False)
        event = PipelineEvent(
            event_type=EventType.VCF_READY,
            source_stage=PipelineStage.GENOMICS,
            patient_id="HG003",
        )
        bus.publish(event)
        history = bus.get_history(patient_id="HG003")
        assert len(history) >= 1
        assert history[0].patient_id == "HG003"
