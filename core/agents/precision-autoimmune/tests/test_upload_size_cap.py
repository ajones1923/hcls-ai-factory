"""Security regression tests for /ingest/upload's size cap (2026-11-03 audit, S15).

The route did `content = await file.read()` — the whole upload into memory, with no
ceiling, on an endpoint that is unauthenticated in the platform's default posture, on a
box whose 128 GB of memory is also the GPU's. The request-size middleware only inspects a
declared Content-Length, so a chunked upload never met it.
"""
import io

import pytest
from fastapi.testclient import TestClient

import api.main as main_mod
from config.settings import settings


class _FakeDocProcessor:
    """Stands in for the real processor: the cap must be hit before any of this runs."""

    def __init__(self):
        self.calls = 0
        self.collection_manager = self

    def process_pdf(self, path, patient_id):
        self.calls += 1
        return [{"doc_type": "note", "specialty": "rheumatology"}]

    def embed_records(self, records):
        return records

    def insert_batch(self, collection, records):
        return len(records)


@pytest.fixture(autouse=True)
def quiet_env(monkeypatch):
    monkeypatch.delenv("HCLS_API_KEY", raising=False)
    store = getattr(main_mod, "_rate_limit_store", None)
    if store is not None:
        store.clear()
    original = settings.MAX_REQUEST_SIZE_MB
    yield
    settings.MAX_REQUEST_SIZE_MB = original


@pytest.fixture
def client_and_proc():
    """The app's lifespan replaces _state on startup, so the fake goes in after that."""
    fake = _FakeDocProcessor()
    with TestClient(main_mod.app) as c:
        previous = main_mod._state.get("doc_processor")
        main_mod._state["doc_processor"] = fake
        try:
            yield c, fake
        finally:
            if previous is None:
                main_mod._state.pop("doc_processor", None)
            else:
                main_mod._state["doc_processor"] = previous


def _upload(client, payload, name="note.pdf"):
    return client.post(
        "/ingest/upload",
        files={"file": (name, io.BytesIO(payload), "application/pdf")},
        data={"patient_id": "SYN-001"},
    )


def test_oversize_upload_is_refused(client_and_proc):
    c, proc = client_and_proc
    settings.MAX_REQUEST_SIZE_MB = 1
    r = _upload(c, b"x" * (2 * 1024 * 1024))
    assert r.status_code == 413
    # Nothing was embedded or inserted: the refusal happens before any processing.
    assert proc.calls == 0


def test_in_route_counter_refuses_a_body_with_no_declared_length(client_and_proc):
    """The middleware only sees a declared Content-Length; this is the other path.

    Called directly rather than over HTTP, because the test client always declares a
    length — which is exactly why the cap cannot live in the middleware alone. The fake
    file hands back 1 MB per read and never ends, so a route without a running total
    would read until the box ran out of memory.
    """
    import asyncio

    from fastapi import HTTPException

    _, proc = client_and_proc
    settings.MAX_REQUEST_SIZE_MB = 2

    class _EndlessUpload:
        filename = "note.pdf"

        def __init__(self):
            self.reads = 0

        async def read(self, size=-1):
            assert size and size > 0, "the route must ask for a bounded chunk"
            self.reads += 1
            assert self.reads < 64, "the route did not stop reading"
            return b"x" * (1024 * 1024)

    endless = _EndlessUpload()
    with pytest.raises(HTTPException) as err:
        asyncio.run(main_mod.upload_document(file=endless, patient_id="SYN-001"))
    assert err.value.status_code == 413
    assert "limit" in str(err.value.detail).lower()
    assert proc.calls == 0


def test_empty_upload_is_refused(client_and_proc):
    c, proc = client_and_proc
    r = _upload(c, b"")
    assert r.status_code == 400
    assert proc.calls == 0


def test_upload_within_the_cap_still_works(client_and_proc):
    c, proc = client_and_proc
    settings.MAX_REQUEST_SIZE_MB = 5
    r = _upload(c, b"%PDF-1.4 tiny")
    assert r.status_code == 200, r.text
    assert proc.calls == 1


def test_non_pdf_is_still_refused(client_and_proc):
    c, proc = client_and_proc
    r = _upload(c, b"x" * 16, name="payload.exe")
    assert r.status_code == 400
    assert proc.calls == 0


def test_the_route_no_longer_reads_the_whole_body_at_once():
    """Guard against a regression back to an unbounded read."""
    import inspect
    source = inspect.getsource(main_mod.upload_document)
    code = "\n".join(l for l in source.splitlines() if not l.strip().startswith("#"))
    assert "await file.read()" not in code
    assert "await file.read(" in code          # a bounded, chunked read
    assert "MAX_REQUEST_SIZE_MB" in code
