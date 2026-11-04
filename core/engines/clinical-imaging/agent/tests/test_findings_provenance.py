"""Canned imaging findings must never be served as if they were measured.

Until 2026-11-03 `WorkflowResult.is_mock` carried `mock_mode` — what the CALLER asked for —
while every workflow but cxr_rapid_findings returned `_mock_inference()` from `infer()` even in
real mode. Three of them (CT head, CT chest, MRI brain) ran a genuine forward pass over a
`torch.randn` volume and then returned canned findings stamped `real_model_used=True`, so a
fixed 12.5 mL hemorrhage reached the API with `is_mock=False`. These tests pin the labelling.
"""
import pytest

from src.models import WorkflowResult, WorkflowStatus
from src.workflows.base import BaseImagingWorkflow


class _Canned(BaseImagingWorkflow):
    """Minimal workflow whose infer() returns canned findings, as the real ones do."""

    WORKFLOW_NAME = "canned_probe"

    def preprocess(self, input_path: str):
        return input_path

    def infer(self, preprocessed):
        return self._mock_inference()

    def postprocess(self, inference_result):
        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.COMPLETED,
            measurements={"hemorrhage_volume_ml": 12.5},
        )

    def _mock_inference(self):
        return {"hemorrhage_volume_ml": 12.5}


class _Measured(_Canned):
    """A workflow that genuinely derives its findings from the input."""

    WORKFLOW_NAME = "measured_probe"

    def infer(self, preprocessed):
        return {"hemorrhage_volume_ml": 3.0}


class _SyntheticPass(_Canned):
    """A real model ran — over synthetic input — and the findings stayed canned."""

    WORKFLOW_NAME = "synthetic_pass_probe"

    def infer(self, preprocessed):
        out = self._mock_inference()
        out[self.SYNTHETIC_PASS_KEY] = True
        return out


def test_canned_findings_are_labelled_mock_even_in_real_mode():
    result = _Canned(mock_mode=False).run("scan.dcm")
    assert result.is_mock is True
    assert result.findings_provenance == "canned"


def test_mock_mode_is_still_labelled_mock():
    result = _Canned(mock_mode=True).run("scan.dcm")
    assert result.is_mock is True
    assert result.findings_provenance == "canned"


def test_measured_findings_are_not_labelled_mock():
    result = _Measured(mock_mode=False).run("scan.dcm")
    assert result.is_mock is False
    assert result.findings_provenance == "measured"


def test_forward_pass_on_synthetic_input_is_not_called_measured():
    """A real model on a torch.randn volume is not a read of the patient's image."""
    result = _SyntheticPass(mock_mode=False).run("scan.dcm")
    assert result.is_mock is True
    assert result.findings_provenance == "canned_after_synthetic_forward_pass"


def test_mock_payload_is_marked_at_the_source():
    """_mock_inference() marks its payload however it is reached, so infer() cannot forget."""
    assert _Canned(mock_mode=False)._mock_inference()[BaseImagingWorkflow.CANNED_KEY] is True


@pytest.mark.parametrize(
    "module_name,workflow_name",
    [
        ("src.workflows.ct_head_hemorrhage", "CTHeadHemorrhageWorkflow"),
        ("src.workflows.ct_chest_lung_nodule", "CTChestLungNoduleWorkflow"),
        ("src.workflows.mri_brain_ms_lesion", "MRIBrainMSLesionWorkflow"),
    ],
)
def test_real_model_used_stamp_is_gone(module_name, workflow_name):
    """`real_model_used` read as "measured here"; it meant "ran on synthetic input"."""
    import importlib
    import inspect

    source = inspect.getsource(importlib.import_module(module_name))
    assert "real_model_used" not in source, (
        f"{workflow_name} still stamps real_model_used; use SYNTHETIC_PASS_KEY"
    )
