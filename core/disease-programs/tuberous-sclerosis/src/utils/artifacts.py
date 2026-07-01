"""
Artifact store (PRD §2.5.5). Content-addressed buckets for cohort artifacts, generated
reports, and large provenance blobs. A MinIO backend (S3-compatible) and a local-filesystem
backend share one interface; `make_artifact_store()` picks MinIO when TSC_USE_MINIO=1, else
the filesystem store (the default, runs with no server).
"""
from __future__ import annotations

import os
from abc import ABC, abstractmethod
from pathlib import Path

from config.settings import settings

BUCKETS = ("tsc-cohort", "tsc-reports", "tsc-provenance")


class ArtifactStore(ABC):
    @abstractmethod
    def put(self, bucket: str, key: str, data: bytes) -> str: ...
    @abstractmethod
    def get(self, bucket: str, key: str) -> bytes: ...
    @abstractmethod
    def uri(self, bucket: str, key: str) -> str: ...


class FsArtifactStore(ArtifactStore):
    def __init__(self, root: Path | None = None) -> None:
        self.root = Path(root or settings.DATA_DIR / "artifacts")
        for b in BUCKETS:
            (self.root / b).mkdir(parents=True, exist_ok=True)

    def _path(self, bucket: str, key: str) -> Path:
        p = self.root / bucket / key
        p.parent.mkdir(parents=True, exist_ok=True)
        return p

    def put(self, bucket, key, data):
        self._path(bucket, key).write_bytes(data)
        return self.uri(bucket, key)

    def get(self, bucket, key):
        return self._path(bucket, key).read_bytes()

    def uri(self, bucket, key):
        return f"file://{self._path(bucket, key)}"


class MinioArtifactStore(ArtifactStore):  # pragma: no cover - needs a MinIO/S3 endpoint
    def __init__(self, endpoint: str = "localhost:9000",
                 access_key: str | None = None, secret_key: str | None = None) -> None:
        from minio import Minio  # noqa: WPS433

        # No hardcoded default credentials — require real secrets from the environment
        # (TSC_MINIO_USER / TSC_MINIO_PASSWORD), fail-closed if absent.
        access_key = access_key or os.environ.get("TSC_MINIO_USER")
        secret_key = secret_key or os.environ.get("TSC_MINIO_PASSWORD")
        if not access_key or not secret_key:
            raise RuntimeError(
                "MinIO credentials missing: set TSC_MINIO_USER and TSC_MINIO_PASSWORD "
                "(no default 'minioadmin' fallback)."
            )
        self.client = Minio(endpoint, access_key=access_key, secret_key=secret_key, secure=False)
        for b in BUCKETS:
            if not self.client.bucket_exists(b):
                self.client.make_bucket(b)

    def put(self, bucket, key, data):
        import io

        self.client.put_object(bucket, key, io.BytesIO(data), length=len(data))
        return self.uri(bucket, key)

    def get(self, bucket, key):
        return self.client.get_object(bucket, key).read()

    def uri(self, bucket, key):
        return f"s3://{bucket}/{key}"


def make_artifact_store() -> ArtifactStore:
    import os

    if os.environ.get("TSC_USE_MINIO") == "1":
        return MinioArtifactStore()
    return FsArtifactStore()
