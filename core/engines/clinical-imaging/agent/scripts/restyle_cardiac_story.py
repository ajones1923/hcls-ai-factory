#!/usr/bin/env python3
"""Bring the cardiac story illustration into agreement with the measured coronary analysis.

`cardiac_story.png` is a generated line illustration (kept pristine as `cardiac_story_source.png`).
It shipped with two things that contradicted the rest of the demo:

  1. it read "72% stenosis", a placeholder, while every other panel reports the MEASURED
     77.9% diameter stenosis from `coronary_analysis.json`; and
  2. its blockage glow sat on the mid vessel, while the measured lesion -- and therefore the
     marker in the 3D tree, the rotating model and the close-up -- is on the proximal segment.

A viewer comparing the four panels would see the illustration disagree with the measurement about
both the number and the place. This script rewrites the number from the analysis JSON and moves
the glow to the proximal segment so the illustration follows the measurement instead of leading it.

The illustration remains an ILLUSTRATION: it is not this patient's anatomy and is not drawn to
scale. It is corrected here only so it does not assert a different finding than the data.

Run:  python3 scripts/restyle_cardiac_story.py
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

BASE = pathlib.Path(__file__).resolve().parents[1]
OUT_DIR = BASE / "data" / "demo" / "coronary"
SRC = OUT_DIR / "cardiac_story_source.png"
DST = OUT_DIR / "cardiac_story.png"
DST_DARK = OUT_DIR / "cardiac_story_dark.png"
ANALYSIS = OUT_DIR / "coronary_analysis.json"

# Measured off the source image rather than guessed -- see the module docstring.
PAPER = np.array([249.0, 247.0, 243.0])     # cream stock
INK = np.array([25.0, 46.0, 61.0])          # navy line work

OLD_GLOW_BOX = (340, 285, 610, 585)         # x0, y0, x1, y1 of the mid-vessel bloom
OLD_LEADER_SEED = (533, 435)                # a point on the callout leader (it is a bent polyline,
                                            # so most straight-line samples miss it)
OLD_LEADER_EXPECT = (200, 700)              # sanity bounds on that component's pixel count
OLD_LEADER_WINDOW = (500, 350, 625, 480)    # search window for that stroke's component
NEW_LESION = (452, 348)                     # proximal segment; matches the measured marker
DIAGRAM_EDGE = (609, 424)                   # left edge of the cross-section inset

PCT_BOX = (645, 437, 722, 472)              # bbox of the "72%" numerals, plus bleed
PCT_CENTRE_X = 680                          # the cross-section inset's centre line
PCT_BASELINE = 466

GLOW_RGB = (247, 148, 42)                   # orange, per the request; warmer than the old amber
GLOW_SRC_RGB = (255, 165, 40)               # the ORIGINAL bloom colour, fitted by residual
                                            # search so it can be inverted rather than covered
GLOW_RADIUS = 78

# Panel D sits beside three dark 3D panels, so it also ships in a dark theme.
DARK_BG = np.array([13.0, 16.0, 21.0])      # matches the portal's panel background
DARK_FG = np.array([214.0, 226.0, 238.0])   # line work, light enough to read


def _font(px: int) -> ImageFont.FreeTypeFont:
    """Quicksand matches the illustration's geometric sans more closely than DejaVu."""
    for p in ("/usr/share/fonts/truetype/quicksand/Quicksand-Medium.ttf",
              "/usr/share/fonts/truetype/quicksand/Quicksand-Regular.ttf",
              "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        if pathlib.Path(p).exists():
            return ImageFont.truetype(p, px)
    return ImageFont.load_default()


def unglow(a: np.ndarray, box, colour=GLOW_SRC_RGB, sigma: float = 25.0) -> None:
    """Undo the old bloom by inverting the composite that produced it, in place.

    The bloom is a normal alpha composite of an orange over the artwork -- verified, not assumed:
    clean line work (25,46,61) appears under the bloom as (92,89,75), i.e. LIGHTER, which rules
    out a multiply, while paper (249,247,243) appears as (248,225,190), which rules out a plain
    add. Both are reproduced by `out = base*(1-a) + colour*a` at a ≈ 0.3.

    So rather than approximate the artwork underneath, recover it: the blue channel gives the
    local alpha wherever the underlying pixel is known to be paper, that sparse estimate is
    smoothed into a continuous field (the bloom is smooth by construction, and the line work it
    crosses is thin enough for the surrounding paper to carry the field across), and the composite
    is then inverted everywhere in the box.

    Two earlier approaches failed here and are worth not repeating: filling with paper erased the
    vessels the bloom covered, and remapping luminance onto a paper/ink ramp left a pale grey haze
    exactly where the bloom had been strongest.
    """
    from scipy import ndimage

    x0, y0, x1, y1 = box
    sub = a[y0:y1, x0:x1].astype(float)
    C = np.asarray(colour, float)

    # Alpha is only identifiable where the underlying pixel is paper; the orange keeps red high,
    # so a bright red channel is a reliable "this was paper" test even deep inside the bloom.
    known = (sub[:, :, 0] > 205).astype(float)
    alpha = np.clip((PAPER[2] - sub[:, :, 2]) / (PAPER[2] - C[2]), 0.0, 0.75)
    num = ndimage.gaussian_filter(alpha * known, sigma)
    den = ndimage.gaussian_filter(known, sigma)
    field = np.clip(num / np.maximum(den, 1e-6), 0.0, 0.75)[:, :, None]

    rebuilt = np.clip((sub - field * C[None, None, :]) / (1.0 - field), 0.0, 255.0)
    a[y0:y1, x0:x1] = rebuilt.round().astype(a.dtype)

    resid = np.abs(rebuilt[sub[:, :, 0] > 205] - PAPER[None, :]).mean()
    print(f"  bloom removed: mean paper residual {resid:.1f}/255 "
          f"(peak alpha {float(field.max()):.2f})")


def denature_stroke(a: np.ndarray, box) -> None:
    """Recolour the illustration's amber "diseased segment" stroke back to ordinary navy, in place.

    Removing the bloom is not enough. Underneath it the artist also drew that stretch of vessel in
    warm brown rather than the navy used everywhere else -- it is line work, not glow, so unglow()
    correctly leaves it alone. But it is the illustration's second way of saying "the blockage is
    HERE", and it points at the mid vessel, so it has to move with the marker.

    The stroke is re-expressed as INK COVERAGE rather than remapped by luminance. The amber is
    lighter than the navy it must become, so a luminance-preserving remap turns the vessel into a
    grey ghost; normalising against the stroke's own darkest tone instead makes its core fully navy
    and its antialiased edge fade to paper, matching the line work either side of it.

    Two gates keep this off everything that is not the stroke: warmth, and actual ink coverage.
    Warmth alone is not enough -- reconstructed paper carries a few units of residual warmth, so a
    warmth-only rule dulled the whole treated area to grey.
    """
    x0, y0, x1, y1 = box
    sub = a[y0:y1, x0:x1].astype(float)
    lum = sub @ np.array([0.299, 0.587, 0.114])
    l_pap = float(PAPER @ np.array([0.299, 0.587, 0.114]))

    warmth = sub[:, :, 0] - sub[:, :, 2]
    sel = (warmth > 20) & (lum < 215)
    if not sel.any():
        print("  recoloured warm stroke: nothing matched")
        return
    l_core = float(np.percentile(lum[sel], 3))

    cov = np.clip((l_pap - lum) / max(l_pap - l_core, 1.0), 0.0, 1.0)
    # The illustrator drew this segment with a soft warm halo around the stroke as well as the
    # stroke itself. Mapping that halo to ink coverage turned it into a grey smudge, so low
    # coverage is pushed to zero: the halo clears to paper, the stroke keeps its full weight.
    cov = np.clip((cov - 0.22) / 0.78, 0.0, 1.0)[:, :, None]
    inked = PAPER[None, None, :] + cov * (INK - PAPER)[None, None, :]

    # No darkness gate. An earlier version faded the effect out above luminance 215 to protect the
    # paper, which instead left the stroke's own light antialiased edge as a tan halo around a navy
    # core. The coverage formula already protects paper: at paper luminance cov is 0, so a paper
    # pixel that slips past the warmth test is rewritten to paper and nothing changes.
    w = np.clip((warmth - 10.0) / 12.0, 0.0, 1.0)[:, :, None]
    a[y0:y1, x0:x1] = (sub * (1 - w) + inked * w).round().astype(a.dtype)
    print(f"  recoloured warm stroke: {int((w > 0.5).sum())} px moved to navy "
          f"(stroke core luminance {l_core:.0f})")


def erase_component(a: np.ndarray, seed, window, expect, thresh: int = 430, grow: int = 2) -> None:
    """Erase exactly the connected stroke that `seed` lands on, and nothing else.

    Painting a fat paper-coloured line along the old leader's endpoints was the obvious approach
    and the wrong one: the leader crosses the heart's right border, so a 9 px brush punched a gap
    in the outline. Selecting the leader's own connected component removes the stroke and leaves
    every line it merely passes over intact.
    """
    from scipy import ndimage

    x0, y0, x1, y1 = window
    sub = a[y0:y1, x0:x1]
    dark = sub.sum(2) < thresh
    lab, _ = ndimage.label(dark, structure=np.ones((3, 3)))
    tag = lab[seed[1] - y0, seed[0] - x0]
    if tag == 0:
        raise SystemExit(f"erase_component: seed {seed} is not on a dark stroke")
    sel = lab == tag
    n = int(sel.sum())
    # Guard, because the failure is silent and destructive: if the seed drifted onto the heart
    # outline or a vessel, that component is far larger and would be wiped out with no error.
    lo, hi = expect
    if not lo <= n <= hi:
        raise SystemExit(f"erase_component: component at {seed} has {n} px, expected {lo}-{hi} "
                         f"-- refusing to erase (the seed is probably on the wrong stroke)")
    ys, xs = np.nonzero(sel)
    print(f"  erasing leader stroke: {n} px, bbox x {xs.min()+x0}-{xs.max()+x0} "
          f"y {ys.min()+y0}-{ys.max()+y0}")
    # Fill from the surrounding paper rather than with a flat constant: the paper carries a faint
    # texture and vignette, and a uniform patch reads as a bright band even when its mean tone is
    # correct. A masked blur inherits the local tone; a little noise restores the grain.
    mask = ndimage.binary_dilation(sel, iterations=grow)
    keep = (~mask).astype(float)[:, :, None]
    src = sub.astype(float) * keep
    num = np.stack([ndimage.gaussian_filter(src[:, :, c], 6.0) for c in range(3)], -1)
    den = ndimage.gaussian_filter(keep[:, :, 0], 6.0)[:, :, None]
    local = num / np.maximum(den, 1e-6)
    rng = np.random.default_rng(7)
    local = local + rng.normal(0.0, 1.8, local.shape)
    sub[mask] = np.clip(local, 0, 255).round().astype(a.dtype)[mask]


def add_glow(img: Image.Image, centre, radius: int, rgb) -> Image.Image:
    """Soft radial bloom, composited the way the source illustration's glow behaves.

    Same compositing model unglow() inverts -- a normal alpha blend of an orange, peaking near the
    0.30 alpha the original bloom used -- so the new marker looks native to the drawing rather than
    like a sticker pasted on top.
    """
    w, h = img.size
    cx, cy = centre
    yy, xx = np.mgrid[0:h, 0:w]
    d = np.hypot(xx - cx, yy - cy) / float(radius)
    field = (np.clip(1.0 - d, 0.0, 1.0) ** 1.9) * 0.34         # smooth falloff to nothing

    base = np.asarray(img).astype(float)
    C = np.array(rgb, float)[None, None, :]
    out = base * (1.0 - field[:, :, None]) + C * field[:, :, None]
    return Image.fromarray(np.clip(out, 0, 255).round().astype(np.uint8))


def to_dark(img: Image.Image) -> Image.Image:
    """Re-render the illustration for a dark UI: dark stock, light line work, accents kept.

    A plain lightness inversion is the obvious approach and is wrong here. It would fix the
    achromatic parts -- cream to near-black, navy to pale blue -- but it also flips the mid-tone
    accents the wrong way: the amber marker sits near lightness 0.63, so inverting sends it to
    0.37, i.e. a dark brown that disappears against a dark background.

    So the two kinds of pixel are treated differently. Near-neutral pixels (paper and navy line
    work) are re-expressed as ink coverage and painted from DARK_BG to DARK_FG, which puts light
    lines on a dark ground with the drawing's weights intact. Saturated pixels (the teal ribbon,
    the amber glow, the gauge) keep their hue and saturation and are only lifted in lightness
    enough to read against the new background.
    """
    a = np.asarray(img).astype(float) / 255.0
    mx, mn = a.max(2), a.min(2)
    chroma = mx - mn
    lum = a @ np.array([0.299, 0.587, 0.114])

    l_pap = float(PAPER @ np.array([0.299, 0.587, 0.114])) / 255.0
    l_ink = float(INK @ np.array([0.299, 0.587, 0.114])) / 255.0
    cov = np.clip((l_pap - lum) / (l_pap - l_ink), 0.0, 1.0)[:, :, None]
    neutral = (DARK_BG[None, None, :] + cov * (DARK_FG - DARK_BG)[None, None, :]) / 255.0

    # Accents: rebuild in HSL with hue and saturation intact and only the lightness retargeted.
    # Done as arithmetic on the hue-normalised colour rather than per-pixel colorsys, both for
    # speed over 1.6 M pixels and so the mapping is one readable expression.
    l = (mx + mn) / 2.0
    denom = np.maximum(1.0 - np.abs(2.0 * l - 1.0), 1e-6)
    sat = np.clip(chroma / denom, 0.0, 1.0)
    shape = (a - mn[:, :, None]) / np.maximum(chroma, 1e-6)[:, :, None]   # min 0, max 1

    # Compress rather than floor. Flooring dark accents to a fixed lightness also dragged the soft
    # amber HALOS up to it, which turned every glow into a pale blob; a linear compression lifts
    # the dark cores into view while leaving the already-light halos roughly where they were, so
    # they can fall away into the background instead of spreading.
    l2 = np.clip(0.32 + 0.55 * l, 0.32, 0.80)
    c2 = (1.0 - np.abs(2.0 * l2 - 1.0)) * np.minimum(1.0, sat * 1.10)
    accent = (l2 - c2 / 2.0)[:, :, None] + c2[:, :, None] * shape

    # Feathered handover, so an accent's antialiased edge does not get a hard outline against the
    # neutral field. Below 0.20 chroma a pixel is treated as line work, which is what the soft
    # halos are closest to.
    blend = np.clip((chroma - 0.14) / 0.06, 0.0, 1.0)[:, :, None]
    out = neutral * (1.0 - blend) + accent * blend
    print(f"  dark theme: {int((chroma > 0.20).sum())} accent px kept in colour, "
          f"{int((chroma <= 0.20).sum())} px remapped to light-on-dark")
    return Image.fromarray(np.clip(out * 255.0, 0, 255).round().astype(np.uint8))


def main() -> int:
    if not SRC.exists():
        sys.exit(f"missing pristine source: {SRC}\n"
                 f"copy the original illustration there before running this script")
    pct = 77.9
    if ANALYSIS.exists():
        pct = round(float(json.loads(ANALYSIS.read_text())["max_diameter_stenosis_pct"]), 1)

    img = Image.open(SRC).convert("RGB")
    a = np.asarray(img).copy()

    # Order matters. Undo the composite first, on pixels exactly as the illustrator left them;
    # then neutralise the warm stroke the bloom was sitting on; erase the leader last, because
    # filling it with final paper BEFORE the inversion made unglow() divide that fill by (1 - a)
    # and blow it out to a white streak.
    unglow(a, OLD_GLOW_BOX)
    erase_component(a, OLD_LEADER_SEED, OLD_LEADER_WINDOW, OLD_LEADER_EXPECT)
    denature_stroke(a, OLD_GLOW_BOX)
    img = Image.fromarray(a)

    img = add_glow(img, NEW_LESION, GLOW_RADIUS, GLOW_RGB)

    # Leader and ring are drawn 4x oversampled then downsampled: the illustration is smooth
    # vector-style artwork, and hard-aliased overlays on top of it read as pasted on.
    ss = 4
    ang = math.atan2(DIAGRAM_EDGE[1] - NEW_LESION[1], DIAGRAM_EDGE[0] - NEW_LESION[0])
    start = (NEW_LESION[0] + math.cos(ang) * 46, NEW_LESION[1] + math.sin(ang) * 46)
    for kind, colour in (("leader", (38, 52, 66)), ("ring", (201, 98, 18))):
        mask = Image.new("L", (img.size[0] * ss, img.size[1] * ss), 0)
        md = ImageDraw.Draw(mask)
        if kind == "leader":
            md.line([(start[0] * ss, start[1] * ss),
                     (DIAGRAM_EDGE[0] * ss, DIAGRAM_EDGE[1] * ss)], fill=255, width=3 * ss)
        else:
            r = 26 * ss
            md.ellipse([NEW_LESION[0] * ss - r, NEW_LESION[1] * ss - r,
                        NEW_LESION[0] * ss + r, NEW_LESION[1] * ss + r],
                       outline=255, width=4 * ss)
        m = np.asarray(mask.resize(img.size, Image.LANCZOS)).astype(float)[:, :, None] / 255.0
        base = np.asarray(img).astype(float)
        img = Image.fromarray(np.clip(
            base * (1 - m) + np.array(colour, float)[None, None, :] * m, 0, 255
        ).round().astype(np.uint8))

    # The number, re-typeset from the measurement.
    d = ImageDraw.Draw(img)
    d.rectangle(PCT_BOX, fill=tuple(PAPER.astype(int)))
    label = f"{pct:g}%"
    f = _font(44)   # matches the cap height of the "72%" it replaces (22 px in Quicksand)
    bb = d.textbbox((0, 0), label, font=f)
    d.text((PCT_CENTRE_X - (bb[2] - bb[0]) / 2 - bb[0], PCT_BASELINE - bb[3]),
           label, font=f, fill=tuple(INK.astype(int)))

    img.save(DST)
    # No dark variant is derived here any more. Inverting the light illustration got panel D
    # onto a dark ground but still read as an inversion; it was replaced by a purpose-drawn
    # dark illustration (cardiac_pathway_dark.png). to_dark() is kept because it is the
    # fallback if that asset is ever lost, but it is no longer part of the pipeline.
    print(f"  wrote {DST.relative_to(BASE)}  ({label} stenosis, lesion moved to {NEW_LESION})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
