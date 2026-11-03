"""The scene must show every horizon the modeller produced.

`_lesion_scene` looped over the literal `(6, 12, 18)` while `horizons_months` is configurable
in `config/trajectory_config.yaml`. Change the config and a horizon disappears from the
volumetric scene without a word — the viz would quietly disagree with the forecast it is
drawing.
"""
from src.viz.build_spec import _lesion_scene


def _quantity(horizon_months):
    return {
        "lesion": "SEGA", "location": "foramen of Monro", "threshold_cm": 1.8,
        "direction": "increasing",
        "observed_months": [-24.0, -12.0, -6.0], "observed_cm": [1.1, 1.3, 1.5],
        "forecast": {f"m{h}": {"mean_cm": 1.6 + 0.1 * i,
                               "pi50": [1.5, 1.7 + 0.1 * i],
                               "pi90": [1.4, 1.9 + 0.1 * i]}
                     for i, h in enumerate(horizon_months)},
        "crossing_probability": {f"m{h}": 0.5 for h in horizon_months},
    }


def _forecast_months(q):
    scene = _lesion_scene(q, provenance={})
    return [kf.month for kf in scene.keyframes if kf.forecast]


def test_the_default_three_horizons_are_drawn():
    assert _forecast_months(_quantity([6, 12, 18])) == [6.0, 12.0, 18.0]


def test_a_fourth_horizon_is_not_dropped():
    """This is what the hardcoded tuple lost."""
    assert _forecast_months(_quantity([6, 12, 18, 24])) == [6.0, 12.0, 18.0, 24.0]


def test_a_reconfigured_set_is_followed():
    assert _forecast_months(_quantity([3, 9])) == [3.0, 9.0]


def test_horizons_are_ordered_even_if_the_dict_is_not():
    q = _quantity([18, 6, 12])
    assert _forecast_months(q) == [6.0, 12.0, 18.0]


def test_observed_history_is_still_drawn_and_not_marked_forecast():
    scene = _lesion_scene(_quantity([6, 12, 18]), provenance={})
    observed = [kf.month for kf in scene.keyframes if not kf.forecast]
    assert observed == [-24.0, -12.0, -6.0]


def test_a_quantity_with_no_forecast_yields_no_forecast_keyframes():
    q = _quantity([6])
    q["forecast"] = {}
    q["crossing_probability"] = {}
    assert _forecast_months(q) == []
