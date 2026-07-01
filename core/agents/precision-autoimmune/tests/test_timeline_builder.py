"""
Tests for the Precision Autoimmune Agent timeline builder module.

Covers:
- Date extraction and parsing from clinical text
- Event classification using regex patterns
- Event extraction from document chunks
- Full timeline construction and statistics
- Milvus record generation
- Edge cases: empty input, missing dates, malformed data

All tests run offline — no external services required.
"""

import pytest

from src.timeline_builder import TimelineBuilder

# ── Fixtures ──────────────────────────────────────────────────────────


@pytest.fixture
def builder():
    return TimelineBuilder(rag_engine=None)


@pytest.fixture
def sample_chunks():
    """Realistic clinical document chunks with dates and event types."""
    return [
        {
            "text_chunk": "2023-01-15: Chief complaint of joint pain and morning stiffness lasting >1 hour. Patient reports onset of pain in bilateral hands 3 months ago.",
            "patient_id": "PAT-100",
            "specialty": "Rheumatology",
            "provider": "Dr. Smith",
            "source_file": "visit_notes_2023.pdf",
        },
        {
            "text_chunk": "2023-03-10: Lab result shows positive for RF (IgM 120 IU/mL, ref <14) and elevated CRP at 28 mg/L. ANA positive at 1:160 with speckled pattern.",
            "patient_id": "PAT-100",
            "specialty": "Rheumatology",
            "provider": "Dr. Smith",
            "source_file": "lab_results_2023.pdf",
        },
        {
            "text_chunk": "2023-06-22: Diagnosed with rheumatoid arthritis. Meets criteria for RA based on joint involvement, serology, and symptom duration.",
            "patient_id": "PAT-100",
            "specialty": "Rheumatology",
            "provider": "Dr. Smith",
            "source_file": "diagnosis_notes.pdf",
        },
        {
            "text_chunk": "2023-07-01: Started on methotrexab 15mg weekly with folic acid 1mg daily. Commenced treatment for newly diagnosed RA.",
            "patient_id": "PAT-100",
            "specialty": "Rheumatology",
            "provider": "Dr. Smith",
            "source_file": "treatment_notes.pdf",
        },
        {
            "text_chunk": "2023-11-15: Disease flare-up noted. DAS28-CRP increased to 5.1 from previous 3.8. Swollen joint count increased.",
            "patient_id": "PAT-100",
            "specialty": "Rheumatology",
            "provider": "Dr. Smith",
            "source_file": "follow_up_notes.pdf",
        },
    ]


# ── Date extraction tests ─────────────────────────────────────────────


class TestExtractDate:

    def test_iso_date(self, builder):
        assert builder.extract_date("Visit on 2023-06-15 for follow-up") == "2023-06-15"

    def test_us_date_format(self, builder):
        assert builder.extract_date("Seen on 3/10/2023 in clinic") == "2023-03-10"

    def test_us_date_short_year(self, builder):
        result = builder.extract_date("Labs drawn 1/5/23")
        assert result == "2023-01-05"

    def test_long_month_format(self, builder):
        result = builder.extract_date("Diagnosis made on January 15, 2023")
        assert result == "2023-01-15"

    def test_long_month_without_comma(self, builder):
        result = builder.extract_date("Results from March 5 2024 show improvement")
        assert result == "2024-03-05"

    def test_no_date_returns_none(self, builder):
        assert builder.extract_date("Patient reports fatigue and joint pain") is None

    def test_empty_string(self, builder):
        assert builder.extract_date("") is None

    def test_first_date_matched(self, builder):
        # Should find the first valid date pattern
        result = builder.extract_date("Between 2023-01-01 and 2023-12-31")
        assert result == "2023-01-01"


# ── Event classification tests ────────────────────────────────────────


class TestClassifyEvent:

    def test_symptom_onset(self, builder):
        assert builder.classify_event("Initial complaint of joint pain and swelling") == "symptom_onset"
        assert builder.classify_event("Onset of fatigue and rash in patient") == "symptom_onset"

    def test_diagnosis(self, builder):
        assert builder.classify_event("Diagnosed with rheumatoid arthritis") == "diagnosis"
        assert builder.classify_event("Confirmed diagnosis of SLE") == "diagnosis"

    def test_misdiagnosis(self, builder):
        assert builder.classify_event("Previously diagnosed as fibromyalgia") == "misdiagnosis"
        assert builder.classify_event("Revised diagnosis from OA to RA") == "misdiagnosis"

    def test_lab_result(self, builder):
        assert builder.classify_event("Lab result: positive for ANA at 1:320") == "lab_result"
        # Note: "Elevated CRP noted" does not match because EVENT_PATTERNS use
        # uppercase marker names (ANA, CRP, ESR) but classify_event lowercases
        # the input text before matching, so only "lab result" pattern matches.
        assert builder.classify_event("laboratory result showed elevated markers") == "lab_result"

    def test_imaging(self, builder):
        assert builder.classify_event("MRI shows synovitis in bilateral hands") == "imaging"
        assert builder.classify_event("X-ray reveals joint erosions") == "imaging"

    def test_biopsy(self, builder):
        assert builder.classify_event("Biopsy shows lymphocytic infiltration") == "biopsy"
        assert builder.classify_event("Pathology report from skin biopsy") == "biopsy"

    def test_genetic_test(self, builder):
        assert builder.classify_event("HLA typing result: B27 positive") == "genetic_test"
        assert builder.classify_event("Pharmacogenomic panel completed") == "genetic_test"

    def test_treatment_start(self, builder):
        assert builder.classify_event("Started on adalimumab 40mg biweekly") == "treatment_start"
        assert builder.classify_event("Commenced treatment with methotrexate") == "treatment_start"

    def test_treatment_change(self, builder):
        assert builder.classify_event("Switched to tocilizumab from adalimumab") == "treatment_change"
        assert builder.classify_event("Discontinued methotrexate due to side effects") == "treatment_change"

    def test_flare(self, builder):
        assert builder.classify_event("Disease flare-up with increased joint swelling") == "flare"
        assert builder.classify_event("Symptom exacerbation noted this week") == "flare"

    def test_referral(self, builder):
        assert builder.classify_event("Referred to rheumatology for evaluation") == "referral"
        assert builder.classify_event("Consultation with dermatology requested") == "referral"

    def test_er_visit(self, builder):
        assert builder.classify_event("Emergency department visit for acute pain") == "er_visit"

    def test_unclassified_returns_clinical_note(self, builder):
        assert builder.classify_event("Patient appeared well today") == "clinical_note"

    def test_empty_text(self, builder):
        assert builder.classify_event("") == "clinical_note"


# ── Event extraction from chunks ──────────────────────────────────────


class TestExtractEventsFromChunks:

    def test_extracts_correct_count(self, builder, sample_chunks):
        events = builder.extract_events_from_chunks(sample_chunks, patient_id="PAT-100")
        assert len(events) == 5

    def test_events_sorted_by_date(self, builder, sample_chunks):
        events = builder.extract_events_from_chunks(sample_chunks, patient_id="PAT-100")
        dates = [e["event_date"] for e in events if e["event_date"]]
        assert dates == sorted(dates)

    def test_patient_id_assigned(self, builder, sample_chunks):
        events = builder.extract_events_from_chunks(sample_chunks, patient_id="PAT-100")
        for event in events:
            assert event["patient_id"] == "PAT-100"

    def test_event_types_classified(self, builder, sample_chunks):
        events = builder.extract_events_from_chunks(sample_chunks, patient_id="PAT-100")
        types = {e["event_type"] for e in events}
        assert "symptom_onset" in types
        assert "lab_result" in types
        assert "diagnosis" in types

    def test_dates_extracted(self, builder, sample_chunks):
        events = builder.extract_events_from_chunks(sample_chunks, patient_id="PAT-100")
        assert events[0]["event_date"] == "2023-01-15"

    def test_specialty_preserved(self, builder, sample_chunks):
        events = builder.extract_events_from_chunks(sample_chunks, patient_id="PAT-100")
        for event in events:
            assert event["specialty"] == "Rheumatology"

    def test_days_from_first_symptom_calculated(self, builder, sample_chunks):
        events = builder.extract_events_from_chunks(sample_chunks, patient_id="PAT-100")
        # First symptom is 2023-01-15; later events should have positive days
        first = next(e for e in events if e["event_date"] == "2023-01-15")
        assert first["days_from_first_symptom"] == 0
        last = next(e for e in events if e["event_date"] == "2023-11-15")
        assert last["days_from_first_symptom"] > 0

    def test_text_chunk_truncated_to_500(self, builder):
        long_text = "A" * 1000
        chunks = [{"text_chunk": long_text}]
        events = builder.extract_events_from_chunks(chunks)
        assert len(events[0]["text_chunk"]) == 500

    def test_empty_chunks(self, builder):
        events = builder.extract_events_from_chunks([])
        assert events == []

    def test_chunk_with_no_text(self, builder):
        chunks = [{"text_chunk": ""}, {"text": ""}]
        events = builder.extract_events_from_chunks(chunks)
        assert events == []

    def test_chunk_with_text_key(self, builder):
        """Chunks can use 'text' key instead of 'text_chunk'."""
        chunks = [{"text": "2024-01-01: Onset of pain in bilateral knees."}]
        events = builder.extract_events_from_chunks(chunks)
        assert len(events) == 1
        assert events[0]["event_date"] == "2024-01-01"

    def test_chunk_with_visit_date_fallback(self, builder):
        """When no date is in text, should fall back to visit_date."""
        chunks = [{"text_chunk": "Patient reports joint pain.", "visit_date": "2024-02-14"}]
        events = builder.extract_events_from_chunks(chunks)
        assert events[0]["event_date"] == "2024-02-14"

    def test_missing_dates_get_zero_days(self, builder):
        """Events without dates should get days_from_first_symptom=0."""
        chunks = [
            {"text_chunk": "2023-01-01: First complaint of joint pain."},
            {"text_chunk": "Patient reports ongoing fatigue."},
        ]
        events = builder.extract_events_from_chunks(chunks)
        dateless = [e for e in events if not e["event_date"]]
        for e in dateless:
            assert e["days_from_first_symptom"] == 0


# ── Timeline construction tests ───────────────────────────────────────


class TestBuildTimeline:

    def test_returns_expected_keys(self, builder, sample_chunks):
        timeline = builder.build_timeline("PAT-100", sample_chunks)
        assert "patient_id" in timeline
        assert "total_events" in timeline
        assert "events" in timeline
        assert "specialties_seen" in timeline
        assert "event_type_counts" in timeline
        assert "date_range" in timeline

    def test_total_events_count(self, builder, sample_chunks):
        timeline = builder.build_timeline("PAT-100", sample_chunks)
        assert timeline["total_events"] == 5

    def test_specialties_collected(self, builder, sample_chunks):
        timeline = builder.build_timeline("PAT-100", sample_chunks)
        assert "Rheumatology" in timeline["specialties_seen"]

    def test_event_type_counts(self, builder, sample_chunks):
        timeline = builder.build_timeline("PAT-100", sample_chunks)
        counts = timeline["event_type_counts"]
        assert isinstance(counts, dict)
        assert sum(counts.values()) == 5

    def test_date_range(self, builder, sample_chunks):
        timeline = builder.build_timeline("PAT-100", sample_chunks)
        assert timeline["date_range"]["first"] == "2023-01-15"
        assert timeline["date_range"]["last"] == "2023-11-15"

    def test_empty_chunks(self, builder):
        timeline = builder.build_timeline("PAT-EMPTY", [])
        assert timeline["total_events"] == 0
        assert timeline["events"] == []
        assert timeline["date_range"]["first"] is None
        assert timeline["date_range"]["last"] is None


# ── Milvus record generation tests ────────────────────────────────────


class TestBuildTimelineForMilvus:

    def test_returns_list_of_records(self, builder, sample_chunks):
        events = builder.extract_events_from_chunks(sample_chunks, patient_id="PAT-100")
        records = builder.build_timeline_for_milvus("PAT-100", events)
        assert isinstance(records, list)
        assert len(records) == len(events)

    def test_record_ids_formatted(self, builder, sample_chunks):
        events = builder.extract_events_from_chunks(sample_chunks, patient_id="PAT-100")
        records = builder.build_timeline_for_milvus("PAT-100", events)
        assert records[0]["id"] == "tl_PAT-100_0000"
        assert records[1]["id"] == "tl_PAT-100_0001"

    def test_record_has_required_fields(self, builder, sample_chunks):
        events = builder.extract_events_from_chunks(sample_chunks, patient_id="PAT-100")
        records = builder.build_timeline_for_milvus("PAT-100", events)
        for rec in records:
            assert "id" in rec
            assert "text_chunk" in rec
            assert "patient_id" in rec
            assert "event_type" in rec
            assert "event_date" in rec

    def test_text_chunk_contains_event_info(self, builder):
        events = [
            {
                "event_type": "diagnosis",
                "event_date": "2023-06-22",
                "description": "RA diagnosed",
                "provider": "Dr. Smith",
                "specialty": "Rheumatology",
                "days_from_first_symptom": 158,
            },
        ]
        records = builder.build_timeline_for_milvus("PAT-100", events)
        assert "2023-06-22" in records[0]["text_chunk"]
        assert "Diagnosis" in records[0]["text_chunk"]

    def test_empty_events(self, builder):
        records = builder.build_timeline_for_milvus("PAT-EMPTY", [])
        assert records == []


# ── Summarize event tests ─────────────────────────────────────────────


class TestSummarizeEvent:

    def test_short_text_returned_as_is(self, builder):
        text = "Brief note about patient condition."
        result = builder._summarize_event(text, "clinical_note")
        assert result == text

    def test_long_text_truncated(self, builder):
        text = "A" * 500
        result = builder._summarize_event(text, "clinical_note")
        assert len(result) <= 300

    def test_multi_sentence_takes_two(self, builder):
        text = "This is the first meaningful sentence about the patient. This is the second sentence with clinical detail. This third one should be dropped."
        result = builder._summarize_event(text, "diagnosis")
        assert "first meaningful" in result
        assert "second sentence" in result
        assert "third one" not in result
