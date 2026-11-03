"""The demo dataset is downloaded, not bundled, and the image could not be built.

Two defects in one place:

* The Dockerfile said `COPY data/ /app/data/`. There is no data/ directory in this engine --
  not gitignored, simply absent -- so `docker build` failed with "no source files were
  specified" and the image could never be built from a checkout. The directory is genuinely
  needed, because scanpy writes downloaded datasets to `settings.datasetdir` (default ./data)
  and the non-root runtime user cannot create it, so it is now created rather than copied.

* `demo=true` called `sc.datasets.pbmc3k()`, which fetches the dataset over the network on
  first use. The Dockerfile called it "the bundled PBMC 3k demo dataset"; it is not bundled.
  On a host with no egress the request surfaced as a bare 500 from a scanpy traceback.
"""
import pathlib
import sys

import pytest

SRC = pathlib.Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

import single_cell_service as svc  # noqa: E402

DOCKERFILE = (pathlib.Path(__file__).resolve().parents[1] / "Dockerfile").read_text()


class TestTheImageCanBeBuilt:
    def test_no_copy_of_a_directory_that_does_not_exist(self):
        # Instructions only: the comment explaining the fix names the old line.
        instructions = [ln.strip() for ln in DOCKERFILE.splitlines()
                        if ln.strip() and not ln.strip().startswith("#")]
        assert not any(ln.startswith("COPY data/") for ln in instructions), instructions
        assert not (pathlib.Path(__file__).resolve().parents[1] / "data").exists()

    def test_the_dataset_directory_is_created_instead(self):
        """scanpy needs it writable, and the runtime user cannot create it."""
        assert "mkdir -p /app/data" in DOCKERFILE

    def test_the_dockerfile_does_not_call_the_dataset_bundled(self):
        assert "bundled PBMC 3k demo dataset" not in DOCKERFILE
        assert "not shipped" in DOCKERFILE or "fetched over the network" in DOCKERFILE


class TestAnUnavailableDatasetSaysWhy:
    def test_a_failed_download_raises_a_runtime_error_that_explains_itself(self, monkeypatch):
        class FakeSettings:
            datasetdir = "/nonexistent"

        class FakeDatasets:
            @staticmethod
            def pbmc3k():
                raise OSError("Network is unreachable")

        fake_sc = type("sc", (), {"settings": FakeSettings(), "datasets": FakeDatasets})
        monkeypatch.setitem(sys.modules, "scanpy", fake_sc)

        with pytest.raises(RuntimeError) as exc:
            svc.load_demo_dataset()
        msg = str(exc.value)
        assert "downloaded on first use rather than shipped" in msg
        assert "SC_DATASET_DIR" in msg
        assert "h5ad_path" in msg
        assert "Network is unreachable" in msg      # the real cause is not swallowed

    def test_the_route_answers_503_not_500(self, monkeypatch):
        """The service is working; the dataset is not available here."""
        from fastapi import HTTPException
        from fastapi.testclient import TestClient

        monkeypatch.setattr(svc, "load_demo_dataset",
                            lambda: (_ for _ in ()).throw(RuntimeError("no dataset here")))
        monkeypatch.setitem(sys.modules, "scanpy", type("sc", (), {}))

        app = svc.create_app(analysis=object())
        client = TestClient(app, raise_server_exceptions=False)
        r = client.post("/analyze", json={"demo": True})
        assert r.status_code == 503
        assert "no dataset here" in r.json()["detail"]


class TestTheDatasetDirIsConfigurable:
    def test_sc_dataset_dir_is_honoured(self, monkeypatch, tmp_path):
        seen = {}

        class FakeSettings:
            datasetdir = "./data"

        class FakeDatasets:
            @staticmethod
            def pbmc3k():
                seen["datasetdir"] = fake_sc.settings.datasetdir
                return "adata"

        fake_sc = type("sc", (), {"settings": FakeSettings(), "datasets": FakeDatasets})
        monkeypatch.setitem(sys.modules, "scanpy", fake_sc)
        monkeypatch.setattr(svc, "DATASET_DIR", str(tmp_path))

        assert svc.load_demo_dataset() == "adata"
        assert seen["datasetdir"] == str(tmp_path)

    def test_the_default_is_left_to_scanpy(self, monkeypatch):
        """Unset means scanpy's own default, which the Dockerfile provisions."""
        monkeypatch.setattr(svc, "DATASET_DIR", "")

        class FakeSettings:
            datasetdir = "./data"

        fake_sc = type("sc", (), {"settings": FakeSettings(),
                                  "datasets": type("d", (), {"pbmc3k": staticmethod(lambda: "a")})})
        monkeypatch.setitem(sys.modules, "scanpy", fake_sc)
        svc.load_demo_dataset()
        assert fake_sc.settings.datasetdir == "./data"
