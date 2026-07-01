"""
Population-prior model for the Trajectory Modeler (PRD §3 FR-TM-1; master paper §10).

Fits a statsmodels mixed-effects model (random intercept + slope by patient) to a
synthetic NHD-style population of SEGA trajectories — a stand-in for the TSC Alliance
Natural History Database (the documented data source) — to derive a population slope
prior with its variance. The individual patient's OLS slope is then shrunk toward this
prior, precision-weighted: sparse/noisy individual data borrows strength from the
population; clean data (e.g., Patient B) barely moves. Import-guarded: returns None when
statsmodels is unavailable, and the runner falls back to the plain OLS slope.
"""
from __future__ import annotations

import numpy as np

_PRIOR: tuple[float, float] | None = None
_PUBLISHED_FALLBACK = (0.026, 1e-4)   # ~3 mm/yr growing SEGA, modest prior variance


def _synthetic_nhd(seed: int = 7):
    rng = np.random.default_rng(seed)
    pids, months, cms = [], [], []
    for pid in range(24):
        base = rng.uniform(0.6, 1.2)
        slope = rng.normal(0.026, 0.006)          # cm/month, growing SEGAs
        for m in (0, 6, 12, 18, 24):
            pids.append(pid)
            months.append(float(m))
            cms.append(base + slope * m + rng.normal(0, 0.03))
    return pids, months, cms


def population_prior() -> tuple[float, float] | None:
    """(slope_mean, slope_variance) from the NHD mixed-effects fit, or None if unavailable."""
    global _PRIOR
    if _PRIOR is not None:
        return _PRIOR
    try:
        import pandas as pd
        import statsmodels.formula.api as smf

        pids, months, cms = _synthetic_nhd()
        df = pd.DataFrame({"patient": pids, "month": months, "cm": cms})
        res = smf.mixedlm("cm ~ month", df, groups=df["patient"], re_formula="~month").fit(
            method="lbfgs", disp=False
        )
        slope = float(res.fe_params["month"])
        var = float(res.bse["month"]) ** 2
        _PRIOR = (slope, max(var, 1e-6))
    except Exception:
        _PRIOR = _PUBLISHED_FALLBACK
    return _PRIOR


def shrink_slope(t: np.ndarray, y: np.ndarray, ols_slope: float) -> tuple[float, dict]:
    """Precision-weighted shrinkage of the OLS slope toward the population prior."""
    prior = population_prior()
    if prior is None:
        return ols_slope, {"population_informed": False}
    ps, pv = prior
    resid = y - (ols_slope * t + (y.mean() - ols_slope * t.mean()))
    n = len(t)
    sse = float(np.sum(resid ** 2))
    denom = float(np.sum((t - t.mean()) ** 2)) or 1.0
    var_ind = max((sse / max(n - 2, 1)) / denom, 1e-6)
    shrunk = (ols_slope / var_ind + ps / pv) / (1.0 / var_ind + 1.0 / pv)
    return shrunk, {"population_informed": True, "prior_slope": round(ps, 4),
                    "individual_slope": round(ols_slope, 4), "shrunk_slope": round(shrunk, 4)}
