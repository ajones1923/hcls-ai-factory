"""Security regression tests for the variant store service (2026-11-03 audit, S5).

`POST /load` passed a caller-supplied filesystem path straight into the VCF reader, on a
service that installed neither the shared API-key gate nor governance. On a box doing
genomics, "any file the service user can read, parsed and returned" means other patients'
VCFs. These tests assert the refusals.
"""
import os

import pytest
from fastapi.testclient import TestClient

from variant_store import VariantStore
import variant_store_service as svc

VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t50\tPASS\t.
chr1\t200\t.\tC\tT\t50\tPASS\t.
"""


@pytest.fixture
def rooted(tmp_path, monkeypatch):
    """A data root with one legitimate VCF in it, and a secret outside it."""
    root = tmp_path / "data"
    (root / "output").mkdir(parents=True)
    good = root / "output" / "cohort.vcf"
    good.write_text(VCF)
    outside = tmp_path / "other-patient.vcf"
    outside.write_text(VCF)
    monkeypatch.setenv("VARIANT_STORE_DATA_ROOT", str(root))
    return root, good, outside


@pytest.fixture
def client(rooted):
    return TestClient(svc.create_app(store=VariantStore()))


def test_absolute_path_outside_the_root_is_refused(client, rooted):
    _, _, outside = rooted
    r = client.post("/load", json={"path": str(outside)})
    assert r.status_code == 400
    assert "data root" in r.json()["detail"]


def test_traversal_out_of_the_root_is_refused(client, rooted):
    r = client.post("/load", json={"path": "output/../../other-patient.vcf"})
    assert r.status_code == 400


def test_system_file_is_refused(client):
    r = client.post("/load", json={"path": "/etc/passwd"})
    assert r.status_code == 400


def test_non_vcf_inside_the_root_is_refused(client, rooted):
    root, _, _ = rooted
    (root / "notes.txt").write_text("not a vcf")
    r = client.post("/load", json={"path": "notes.txt"})
    assert r.status_code == 400
    assert "must name" in r.json()["detail"]


def test_symlink_out_of_the_root_is_refused(client, rooted):
    root, _, outside = rooted
    link = root / "output" / "sneaky.vcf"
    link.symlink_to(outside)
    r = client.post("/load", json={"path": "output/sneaky.vcf"})
    assert r.status_code == 400


def test_missing_file_is_a_refusal_not_a_crash(client):
    r = client.post("/load", json={"path": "output/nope.vcf"})
    assert r.status_code == 400


def test_empty_path_is_refused(client):
    assert client.post("/load", json={"path": ""}).status_code == 400


def test_error_does_not_echo_a_host_path(client, rooted):
    _, _, outside = rooted
    detail = client.post("/load", json={"path": str(outside)}).json()["detail"]
    assert str(outside) not in detail


def test_a_legitimate_relative_path_still_loads(client, rooted):
    r = client.post("/load", json={"path": "output/cohort.vcf"})
    assert r.status_code == 200, r.text
    assert r.json()["loaded"] == 2


def test_a_legitimate_absolute_path_inside_the_root_still_loads(client, rooted):
    _, good, _ = rooted
    r = client.post("/load", json={"path": str(good)})
    assert r.status_code == 200, r.text
    assert r.json()["loaded"] == 2


def test_data_root_defaults_to_the_engine_data_dir(monkeypatch):
    monkeypatch.delenv("VARIANT_STORE_DATA_ROOT", raising=False)
    assert svc.data_root().name == "data"


def test_service_installs_the_shared_api_key_gate(rooted, monkeypatch):
    """The service used to sit outside the platform's gate entirely."""
    monkeypatch.setenv("HCLS_API_KEY", "vs-test-key")
    app = svc.create_app(store=VariantStore())
    with TestClient(app) as c:
        assert c.post("/load", json={"path": "output/cohort.vcf"}).status_code == 401
        assert c.get("/stats").status_code == 401
        ok = c.post("/load", json={"path": "output/cohort.vcf"},
                    headers={"X-API-Key": "vs-test-key"})
        assert ok.status_code == 200, ok.text
        # Probes must not fail closed with the service.
        assert c.get("/health").status_code == 200


def test_governance_middleware_is_installed(monkeypatch):
    monkeypatch.delenv("HCLS_API_KEY", raising=False)
    app = svc.create_app(store=VariantStore())
    names = [getattr(m.kwargs.get("dispatch", None), "__name__", "") for m in app.user_middleware]
    assert "_governance" in names
