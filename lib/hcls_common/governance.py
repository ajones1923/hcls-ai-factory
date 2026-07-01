"""
Dual-model "second opinion" governance (A5).

For high-stakes predictions (toxicity, cell-type, pathogenicity), run two *independent*
models and surface both plus their agreement — turning model disagreement into a visible,
actionable confidence signal. Pairs with the output honesty gate (verify_gate) and the
input-validation gate (A2): inputs validated, two models consulted, claims refuted.
"""
from __future__ import annotations

from typing import Any, Callable


def second_opinion(
    value_a: Any,
    value_b: Any,
    *,
    label: str = "prediction",
    model_a: str = "model_a",
    model_b: str = "model_b",
    numeric_tol: float | None = None,
) -> dict[str, Any]:
    """Compare two independent model outputs and report agreement as a confidence signal."""
    if (numeric_tol is not None
            and isinstance(value_a, (int, float)) and not isinstance(value_a, bool)
            and isinstance(value_b, (int, float)) and not isinstance(value_b, bool)):
        agree = abs(value_a - value_b) <= numeric_tol
    else:
        agree = value_a == value_b
    return {
        "label": label,
        model_a: value_a,
        model_b: value_b,
        "agree": agree,
        "confidence": "high" if agree else "low",
        "note": None if agree
        else f"independent models disagree on {label} — low confidence, recommend review",
    }


def dual_predict(
    predictor_a: Callable[[Any], Any],
    predictor_b: Callable[[Any], Any],
    x: Any,
    *,
    extract: Callable[[Any], Any] = lambda r: r,
    label: str = "prediction",
    model_a: str = "model_a",
    model_b: str = "model_b",
    numeric_tol: float | None = None,
) -> dict[str, Any]:
    """Run two predictors on the same input, compare the extracted value, keep both raw results."""
    ra, rb = predictor_a(x), predictor_b(x)
    out = second_opinion(extract(ra), extract(rb), label=label,
                         model_a=model_a, model_b=model_b, numeric_tol=numeric_tol)
    out["raw_a"], out["raw_b"] = ra, rb
    return out
