"""Warfarin dosing: a refusal has to reach the caller, and an unreadable genotype is not wild-type.

Found 2026-11-04, hand-reviewing pharmacogenomics for check 4. Three findings, all in the one
endpoint, and the first one is the reason the other two were invisible.

1. THERE ARE TWO IWPC IMPLEMENTATIONS, and only one of them was disabled.

   `src/dosing.py` carries the IWPC algorithm with a mis-transcribed coefficient block, and
   refuses to return a dose because the error is race-specific (CLINICAL_REVIEW §21).
   `api/routes/pgx_clinical.py::_warfarin_dosing_impl` carries a SECOND, independent copy and
   computes the dose inline -- it never calls the library. So "IWPC warfarin dosing is
   DISABLED" was true of the library and false of the HTTP endpoint, which is the surface a
   caller actually touches.

   The route's copy does NOT have the §21 defect: it assigns -0.2760 to Black/African American
   race, which is exactly the term the broken copy misapplies. That is why this test does not
   disable the endpoint. It pins the divergence instead, so the next person to read §21 learns
   that a second copy exists before concluding the capability is off.

2. THE LIBRARY'S REFUSAL WAS SWALLOWED. It raised a bare `RuntimeError`; the route caught
   `Exception` and returned HTTP 500 "Internal processing error". A careful explanation of a
   race-specific dosing error reached the server log and never the caller -- who would read a
   500 as a transient fault and retry. It now raises `CapabilityWithdrawn` and the route
   answers 501 with the reason.

3. AN UNREADABLE GENOTYPE WAS SCORED AS REFERENCE. An unrecognised VKORC1 became GG and an
   unrecognised CYP2C9 became *1/*1, with a warning appended beside a confidently-typed
   `predicted_weekly_dose_mg`. A warning next to a number is not a refusal; the number is what
   gets used. Both defaults bias the dose UPWARD, because every scoreable alternative reduces
   it (VKORC1 AA -1.6974, CYP2C9 *3/*3 -2.3312) -- the most dose-inflating default available,
   on a narrow-therapeutic-index anticoagulant.

   `src/dosing.py` already states the principle this violated, in `IndeterminateGenotype`: such
   a genotype "is not a normal metabolizer, not wild-type, and not a dose." The module had the
   rule; the HTTP path did not follow it.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

ROUTE_SRC = (ROOT / "api" / "routes" / "pgx_clinical.py").read_text()
LIB_SRC = (ROOT / "src" / "dosing.py").read_text()


def client():
    from fastapi.testclient import TestClient

    from api.main import app
    return TestClient(app)


def body(**over):
    payload = {"age": 65, "height_cm": 170, "weight_kg": 80,
               "vkorc1_genotype": "GG", "cyp2c9_genotype": "*1/*1",
               "race": "white", "amiodarone": False, "enzyme_inducer": False}
    payload.update(over)
    return payload


class TestTheWithdrawalIsNotA500:
    def test_the_library_raises_a_named_exception(self):
        from src.dosing import CapabilityWithdrawn
        assert issubclass(CapabilityWithdrawn, RuntimeError)
        assert "CapabilityWithdrawn(" in LIB_SRC, "the IWPC refusal no longer raises it"

    def test_the_route_maps_it_to_501_and_passes_the_reason(self):
        assert "except CapabilityWithdrawn" in ROUTE_SRC
        assert "status_code=501" in ROUTE_SRC
        assert 'detail=str(exc)' in ROUTE_SRC, (
            "the reason must reach the caller; a generic detail string is what this fixed")

    def test_the_generic_handler_no_longer_sees_it_first(self):
        """Ordering is the whole fix: a broad `except Exception` placed first eats it again."""
        # Scoped to the warfarin handler: the file has several generic 500 handlers for other
        # endpoints, and an unscoped index() finds the first one in the file instead of this one.
        block = ROUTE_SRC.split("async def warfarin_dosing", 1)[1].split("def _warfarin_dosing_impl", 1)[0]
        assert "except CapabilityWithdrawn" in block
        assert block.index("except CapabilityWithdrawn") < block.index("except Exception")


class TestAnUnreadableGenotypeIsDeclined:
    @pytest.mark.parametrize("field,value", [
        ("vkorc1_genotype", "XX"),
        ("vkorc1_genotype", ""),
        ("cyp2c9_genotype", "*1/*9"),
        ("cyp2c9_genotype", "unknown"),
    ])
    def test_no_dose_is_returned(self, field, value):
        r = client().post("/v1/pgx/dosing/warfarin", json=body(**{field: value}))
        assert r.status_code == 422, (
            f"{field}={value!r} returned {r.status_code}; an unscoreable genotype must not "
            f"produce a dose")
        assert "predicted_weekly_dose_mg" not in r.text

    @pytest.mark.parametrize("field,value", [
        ("vkorc1_genotype", "XX"), ("cyp2c9_genotype", "*1/*9")])
    def test_the_reason_says_why_and_where_to_read_more(self, field, value):
        r = client().post("/v1/pgx/dosing/warfarin", json=body(**{field: value}))
        detail = r.json().get("detail", "")
        assert "cannot be scored" in detail
        assert "§21" in detail or "CLINICAL_REVIEW" in detail

    def test_the_reference_substitution_is_gone_from_the_source(self):
        """Guards the specific wording that used to make this a warning instead of a refusal."""
        for gone in ("Using GG (reference)", "Using *1/*1 (reference)"):
            assert gone not in ROUTE_SRC, f"{gone!r} is back -- the dose is being guessed again"

    def test_a_scoreable_request_still_works(self):
        """The refusal must not have turned into a blanket outage."""
        r = client().post("/v1/pgx/dosing/warfarin", json=body())
        assert r.status_code == 200, r.text
        assert r.json()["predicted_weekly_dose_mg"] > 0


class TestTheTwoCopiesAreBothStillHere:
    """Pinned so §21's "disabled" is not read as covering the HTTP endpoint.

    If either of these fails, the duplication was resolved -- which is a good outcome, and the
    right response is to update CLINICAL_REVIEW §21 and delete this class, not to restore a
    second copy of a warfarin dosing algorithm.
    """

    def test_the_library_copy_refuses(self):
        from src.dosing import CapabilityWithdrawn, DosingCalculator
        calc = DosingCalculator()
        with pytest.raises(CapabilityWithdrawn):
            calc.warfarin_dose(cyp2c9_diplotype="*1/*1", vkorc1_genotype="GG")

    def test_the_route_has_its_own_inline_copy(self):
        assert "_warfarin_dosing_impl" in ROUTE_SRC
        assert "5.6044" in ROUTE_SRC, "the route's own IWPC constant is gone"

    def test_the_route_copy_assigns_the_race_coefficient_to_race(self):
        """The one difference that matters: this is the term §21 says the library misapplies."""
        race_block = ROUTE_SRC.split("# Race contribution", 1)[1][:600]
        assert "0.2760" in race_block, (
            "the -0.2760 Black/African American coefficient is no longer in the race block; "
            "if it moved, re-read CLINICAL_REVIEW §21 before trusting this endpoint")
