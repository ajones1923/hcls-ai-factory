"""Threshold-state palette (FR-VZ-5). Colorblind-safe; state is also carried by a text
label and opacity so colour is never the sole channel."""
from __future__ import annotations

# Okabe-Ito-derived, deuteranopia/protanopia-safe: teal (clear) -> amber (approaching) -> vermillion (at/above)
STATE_RGB = {
    "below": (0.0, 0.62, 0.45),
    "approaching": (0.90, 0.62, 0.0),
    "at_or_above": (0.84, 0.0, 0.10),
}
ENVELOPE_RGB = (0.50, 0.69, 1.0)        # the uncertainty cloud
THRESHOLD_RGB = (0.84, 0.20, 0.20)      # the intervention-threshold membrane
ANATOMY_RGB = (0.20, 0.25, 0.40)        # stylized context anatomy
RECOVERY_RGB = (1.0, 0.78, 0.20)        # mosaic-recovery highlight (gold)

# ACMG classification -> body colour for the population view
CLASS_RGB = {
    "Pathogenic": (0.84, 0.19, 0.17),
    "Likely Pathogenic": (0.91, 0.51, 0.13),
    "Variant of Uncertain Significance": (0.55, 0.58, 0.62),
    "Likely Benign": (0.23, 0.65, 0.45),
    "Benign": (0.18, 0.55, 0.36),
    "No variant identified": (0.30, 0.33, 0.40),
}


def class_rgb(classification: str | None) -> tuple:
    return CLASS_RGB.get(classification or "", (0.40, 0.43, 0.50))


def state_for_value(value: float, threshold: float, direction: str) -> str:
    """Threshold state for an OBSERVED value."""
    if direction == "increasing":
        if value >= threshold:
            return "at_or_above"
        return "approaching" if value >= 0.85 * threshold else "below"
    if value <= threshold:
        return "at_or_above"
    return "approaching" if value <= 1.15 * threshold else "below"


def state_for_probability(p: float) -> str:
    """Threshold state for a FORECAST horizon, from its crossing probability."""
    if p >= 0.5:
        return "at_or_above"
    return "approaching" if p >= 0.15 else "below"
