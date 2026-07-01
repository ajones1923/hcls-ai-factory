"""
Agent 3 — TSC-Trajectory Modeler (PRD §3 FR-TM-1..6; master paper §10).

The deliberately CLASSICAL agent (NOT an LLM). Forecasts every monitored quantity —
SEGA and AML size (cm), renal function (eGFR), and seizure frequency — at 6/12/18 months
with 50%/90% prediction intervals, a Bayesian growth slope (population-prior shrinkage for
SEGA, weak-prior conjugate update otherwise), a survival-style probability of crossing the
clinical threshold by each horizon, a graded crossing (likely/possible/unlikely), and a
surveillance-cadence recommendation that tightens follow-up vs. the ITSC floor when risk is
elevated. Config-driven via config/trajectory_config.yaml. Only the prose summary uses an
LLM tier (Haiku); the summary introduces no new number.
"""
from __future__ import annotations

import math

import numpy as np
import yaml

from config.settings import settings
from src.agents.base import Agent, AgentContext, AgentOutput
from src.agents.trajectory_modeler.population import shrink_slope
from src.utils.model_router import get_router
from src.utils.provenance import Provenance, stable_hash

SEGA_DISCUSSION_THRESHOLD_CM = 1.8   # kept for back-compat imports (eval harness)
AML_DISCUSSION_THRESHOLD_CM = 4.0

# cohort keys per quantity: (series_key, months_key, location_key, default_location)
_SERIES_KEYS = {
    "SEGA": ("sega_series_cm", "sega_months", "sega_location", "foramen of Monro"),
    "AML": ("aml_series_cm", "aml_months", "aml_location", "kidney"),
    "eGFR": ("egfr_series", "egfr_months", None, "kidney"),
    "seizure_frequency": ("seizure_freq_series", "seizure_months", None, "clinical"),
}

_CFG: dict | None = None


def _config() -> dict:
    global _CFG
    if _CFG is None:
        with open(settings.TRAJECTORY_CONFIG_PATH) as f:
            _CFG = yaml.safe_load(f)
    return _CFG


def _norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _gp_posterior(t: np.ndarray, y: np.ndarray, horizons: list[int]):
    """Gaussian-process mean+std at each horizon; falls back to OLS residual std."""
    try:
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

        kernel = ConstantKernel(1.0) * RBF(length_scale=12.0) + WhiteKernel(noise_level=0.01)
        gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True, alpha=1e-6)
        gp.fit(t.reshape(-1, 1), y)
        mu, std = gp.predict(np.array(horizons).reshape(-1, 1), return_std=True)
        return mu, np.maximum(std, 0.03)
    except Exception:
        slope, intercept = np.polyfit(t, y, 1)
        resid = y - (slope * t + intercept)
        sigma = max(0.05, float(np.std(resid)))
        mu = slope * np.array(horizons) + intercept
        return mu, np.array([sigma * np.sqrt(1 + h / 12.0) for h in horizons])


def _slope_var(t: np.ndarray, y: np.ndarray, ols_slope: float) -> float:
    resid = y - (ols_slope * t + (y.mean() - ols_slope * t.mean()))
    n = len(t)
    sse = float(np.sum(resid ** 2))
    denom = float(np.sum((t - t.mean()) ** 2)) or 1.0
    return max((sse / max(n - 2, 1)) / denom, 1e-6)


def _bayesian_slope(t, y, ols_slope, use_prior, weak_prior) -> tuple[float, dict]:
    """Bayesian growth-slope estimate. Population-prior shrinkage when available (SEGA);
    otherwise a conjugate-normal update against a weakly-informative prior."""
    if use_prior:
        shrunk, info = shrink_slope(t, y, float(ols_slope))
        return shrunk, {"method": "population-prior shrinkage (NHD mixed-effects)",
                        "population_informed": info.get("population_informed", True),
                        "prior_slope": info.get("prior_slope"),
                        "individual_slope": info.get("individual_slope", round(float(ols_slope), 4)),
                        "posterior_slope": info.get("shrunk_slope", round(shrunk, 4))}
    var_ind = _slope_var(t, y, ols_slope)
    pm, pv = weak_prior["slope_mean"], weak_prior["slope_variance"]
    post_var = 1.0 / (1.0 / var_ind + 1.0 / pv)
    post_mean = (ols_slope / var_ind + pm / pv) * post_var
    return post_mean, {"method": "conjugate-normal update (weak prior)",
                       "population_informed": False, "prior_slope": pm,
                       "individual_slope": round(float(ols_slope), 4),
                       "posterior_slope": round(post_mean, 4),
                       "posterior_sd": round(math.sqrt(post_var), 4)}


def _model_quantity(name: str, qcfg: dict, series, months, location, cfg: dict) -> dict:
    t = np.array(months[: len(series)], dtype=float)
    y = np.array(series, dtype=float)
    direction = qcfg["direction"]
    threshold = float(qcfg["threshold"])
    horizons = cfg["horizons_months"]
    z = cfg["prediction_intervals"]
    weak = cfg.get("weak_prior", {"slope_mean": 0.0, "slope_variance": 1.0})

    ols_slope = float(np.polyfit(t, y, 1)[0])
    slope, bayes = _bayesian_slope(t, y, ols_slope, qcfg.get("population_prior", False), weak)
    intercept = float(y.mean() - slope * t.mean())
    _, std = _gp_posterior(t, y, horizons)

    forecast, crossing_prob = {}, {}
    for i, h in enumerate(horizons):
        mu = slope * h + intercept
        sd = float(std[i])
        forecast[f"m{h}"] = {
            "mean_cm": round(mu, 3),
            "pi50": [round(mu - z["pi50"] * sd, 3), round(mu + z["pi50"] * sd, 3)],
            "pi90": [round(mu - z["pi90"] * sd, 3), round(mu + z["pi90"] * sd, 3)],
        }
        zc = (mu - threshold) / sd if direction == "increasing" else (threshold - mu) / sd
        crossing_prob[f"m{h}"] = round(_norm_cdf(zc), 3)

    heading = (slope > 0) if direction == "increasing" else (slope < 0)
    t_cross = round(float((threshold - intercept) / slope), 1) if heading and slope != 0 else None
    if t_cross is not None and t_cross < 0:
        t_cross = None
    past_now = (y[-1] >= threshold) if direction == "increasing" else (y[-1] <= threshold)

    p18 = crossing_prob[f"m{horizons[-1]}"]
    g = cfg["crossing_grades"]
    grade = "likely" if (p18 >= g["likely"] or past_now) else ("possible" if p18 >= g["possible"] else "unlikely")

    cad = cfg["surveillance_cadence"]
    floor = qcfg["surveillance"]["itsc_floor_months"]
    rec = max(round(floor * cad[grade]), cad["min_interval_months"])
    surveillance = {
        "modality": qcfg["surveillance"]["modality"],
        "itsc_floor_months": floor, "recommended_interval_months": rec,
        "rationale": (f"{grade} threshold crossing by 18mo (p={p18}) — "
                      + ("tighten vs ITSC floor" if rec < floor else "hold ITSC floor")),
    }

    return {
        "lesion": qcfg.get("label", name), "quantity": name, "unit": qcfg["unit"],
        "direction": direction, "is_lesion": bool(qcfg.get("is_lesion", False)),
        "location": location, "observed_cm": list(series),
        "observed_months": [float(m) for m in months[: len(series)]],
        "model": ("mixed-effects population-prior shrinkage + Gaussian-process intervals"
                  if qcfg.get("population_prior") else "Bayesian (weak-prior) trend + GP intervals"),
        "slope_cm_per_month": round(slope, 4), "population_prior": bayes, "bayesian": bayes,
        "forecast": forecast, "threshold_cm": threshold, "months_to_threshold": t_cross,
        "median_time_to_threshold_months": t_cross,
        "crossing_probability": crossing_prob, "crossing_grade": grade,
        "crosses_in_12_18mo_window": t_cross is not None and 12 <= t_cross <= 18,
        "at_or_above_threshold": bool(past_now) or (t_cross is not None and t_cross <= 0),
        "surveillance_recommendation": surveillance,
    }


def _model_lesion(name, series, months, location, threshold, use_prior) -> dict:
    """Back-compat wrapper (used by older callers/tests): model a single cm-lesion."""
    qcfg = {"label": name, "unit": "cm", "direction": "increasing", "threshold": threshold,
            "is_lesion": True, "population_prior": use_prior,
            "surveillance": {"modality": "MRI", "itsc_floor_months": 12}}
    return _model_quantity(name, qcfg, series, months, location, _config())


class TrajectoryModeler(Agent):
    name = "trajectory_modeler"

    def run(self, ctx: AgentContext) -> AgentOutput:
        cfg = _config()
        quantities, modeled = [], {}
        for qname, qcfg in cfg["quantities"].items():
            skey, mkey, lkey, ldef = _SERIES_KEYS[qname]
            series = ctx.cohort.get(skey)
            if series and len(series) >= 2:
                months = ctx.cohort.get(mkey, [-24, -12, -6])
                location = ctx.cohort.get(lkey, ldef) if lkey else ldef
                quantities.append(_model_quantity(qname, qcfg, series, months, location, cfg))
                modeled[qname] = series

        lesions = [q for q in quantities if q["is_lesion"]]   # back-compat view (SEGA/AML)

        prov_model = Provenance(
            agent=self.name, step="forecasting", tier="classical",
            model_id="statsmodels MixedLM population prior + Bayesian slope + sklearn GaussianProcessRegressor",
            input_hash=stable_hash(modeled),
        )
        router = get_router()
        _, prov_summary = router.call(
            self.name, "prose_summary",
            system="Summarize the structured trajectory forecast in one sentence. Introduce no new number.",
            prompt=str(quantities), prompt_template_version="prose_summary/v0",
        )
        data = {"quantities": quantities, "lesions": lesions,
                "_build_status": "config-driven multi-quantity (SEGA/AML/eGFR/seizure); "
                                 "Bayesian slope + survival-style crossing probability + cadence"}
        return self.ok(self.name, data, [prov_model, prov_summary])
