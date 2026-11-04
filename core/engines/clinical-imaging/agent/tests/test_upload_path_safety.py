"""An uploaded filename must not decide where on disk the service writes.

`/live/upload` and `/live/upload-series` built their temp path with
`os.path.join(tmp_dir, file.filename)`. The name comes from whoever is uploading, so a
filename containing "../" segments escaped the temp directory and the DICOM bytes landed
wherever the service user could write.
"""
import os
import tempfile

import pytest
from fastapi import HTTPException

from api.routes.live_analysis import _safe_upload_path


@pytest.fixture
def tmp_dir():
    with tempfile.TemporaryDirectory() as d:
        yield d


class TestEscapesAreRefused:
    @pytest.mark.parametrize("name", [
        "../../../../etc/cron.d/pwn",
        "../outside.dcm",
        "a/../../b.dcm",
        "/etc/passwd",
        "..",
    ])
    def test_traversal_is_contained_or_refused(self, tmp_dir, name):
        try:
            path = _safe_upload_path(tmp_dir, name, "upload.dcm")
        except HTTPException as exc:
            assert exc.status_code == 400
            return
        root = os.path.realpath(tmp_dir)
        assert path.startswith(root + os.sep), f"{name!r} escaped to {path}"

    def test_absolute_path_does_not_write_outside(self, tmp_dir):
        path = _safe_upload_path(tmp_dir, "/etc/passwd", "upload.dcm")
        assert os.path.realpath(tmp_dir) == os.path.dirname(path)
        assert os.path.basename(path) == "passwd"


class TestOrdinaryNamesStillWork:
    def test_plain_name_is_kept(self, tmp_dir):
        path = _safe_upload_path(tmp_dir, "series1.dcm", "upload.dcm")
        assert os.path.basename(path) == "series1.dcm"
        assert path.startswith(os.path.realpath(tmp_dir) + os.sep)

    def test_missing_name_uses_the_fallback(self, tmp_dir):
        assert os.path.basename(_safe_upload_path(tmp_dir, None, "slice_0001.dcm")) == "slice_0001.dcm"

    def test_empty_name_uses_the_fallback(self, tmp_dir):
        assert os.path.basename(_safe_upload_path(tmp_dir, "", "slice_0002.dcm")) == "slice_0002.dcm"

    def test_nul_byte_is_stripped(self, tmp_dir):
        path = _safe_upload_path(tmp_dir, "scan\x00.dcm", "upload.dcm")
        assert "\x00" not in path
