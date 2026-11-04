"""The trajectory prior must not claim a real registry it was not fitted to.

Found by the claim-integrity review of this subject, 2026-11-04, check 2.

The Trajectory Modeler shrinks a patient's observed growth slope toward a population prior.
It reported that prior as `method: "population-prior shrinkage (NHD mixed-effects)"` with
`population_informed: true`. NHD is the TSC Alliance Natural History Database -- a real
registry, named in the program's own documentation as the data source.

The prior is not fitted to it. It is fitted to `_synthetic_nhd()`: 24 simulated patients
whose slopes are drawn from `rng.normal(0.026, 0.006)`. The mixed-effects model then recovers
0.0249 -- the constant it was generated from, plus sampling noise -- so the fit carries no
information the hardcoded number did not already have, and the statsmodels call does not even
converge cleanly ("MLE may be on the boundary", "Hessian not positive definite").

Same class as the fabricated GPU telemetry removed in v2.2.0: a number that reads as derived
from a real source and is invented. It is worse here only in that it feeds a clinical
trajectory a clinician could act on.

The shrinkage arithmetic is deliberately unchanged -- moving a clinical value needs a primary
source, and there is none. What changed is the claim.
"""
import numpy as np
import pytest

from src.agents.trajectory_modeler.population import population_prior, shrink_slope
from src.agents.trajectory_modeler.runner import _bayesian_slope


@pytest.fixture
def series():
    """A growing lesion: four observations over 18 months."""
    t = np.array([0.0, 6.0, 12.0, 18.0])
    y = np.array([1.10, 1.28, 1.44, 1.63])
    return t, y, float(np.polyfit(t, y, 1)[0])


class TestItDoesNotClaimTheRealRegistry:
    def test_the_method_string_does_not_present_itself_as_the_tsc_alliance_nhd(self, series):
        t, y, ols = series
        _, info = _bayesian_slope(t, y, ols, True, {"slope_mean": 0.0, "slope_variance": 1.0})
        method = info["method"]
        assert "SYNTHETIC" in method, method
        assert "not the TSC Alliance NHD" in method, method

    def test_it_does_not_report_being_population_informed(self, series):
        """No population informed it. 24 simulated patients are not a population."""
        t, y, ols = series
        _, info = _bayesian_slope(t, y, ols, True, {"slope_mean": 0.0, "slope_variance": 1.0})
        assert info["population_informed"] is False

    def test_the_synthetic_basis_is_flagged_in_the_payload(self, series):
        """A caller must be able to see this without reading the module docstring."""
        t, y, ols = series
        _, info = _bayesian_slope(t, y, ols, True, {"slope_mean": 0.0, "slope_variance": 1.0})
        assert info["prior_is_synthetic"] is True


class TestWhyTheClaimWasWrong:
    def test_the_fit_recovers_the_constant_it_was_generated_from(self):
        """Pinning the circularity, so nobody re-reads this as a real fit later.

        `_synthetic_nhd` draws slopes from N(0.026, 0.006). If the fitted prior lands on
        0.026 the mixed-effects model has told us only what we told it.
        """
        slope, _ = population_prior()
        assert abs(slope - 0.026) < 0.005, (
            f"fitted prior {slope} vs the 0.026 its own generator uses")

    def test_the_shrinkage_arithmetic_is_unchanged(self, series):
        """The fix corrects a claim, not a clinical value. The posterior must not move."""
        t, y, ols = series
        shrunk, _ = shrink_slope(t, y, ols)
        _, info = _bayesian_slope(t, y, ols, True, {"slope_mean": 0.0, "slope_variance": 1.0})
        assert info["posterior_slope"] == pytest.approx(round(shrunk, 4))


class TestTheHonestPathIsStillHonest:
    def test_the_weak_prior_path_still_reports_itself_correctly(self, series):
        """Quantities with population_prior: false were always labelled truthfully."""
        t, y, ols = series
        _, info = _bayesian_slope(t, y, ols, False, {"slope_mean": 0.0, "slope_variance": 1.0})
        assert info["method"] == "conjugate-normal update (weak prior)"
        assert info["population_informed"] is False
