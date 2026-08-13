#!/usr/bin/env python3
"""
Generate AI-annotated versions of medical images for demo purposes.

Creates realistic radiology AI detection overlays including bounding boxes,
heatmap overlays, measurement lines, severity badges, and diagnostic labels.
"""

import os
import math
from PIL import Image, ImageDraw, ImageFont, ImageFilter

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(BASE_DIR, "data", "sample_images")
FULLRES_DIR = os.path.join(SRC_DIR, "fullres")
OUT_DIR = os.path.join(SRC_DIR, "annotated")

# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------
NVIDIA_GREEN = (118, 185, 0)
RED = (255, 51, 51)
ORANGE = (255, 140, 0)
CYAN = (0, 212, 255)
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
DARK_SHADOW = (0, 0, 0, 180)

# ---------------------------------------------------------------------------
# Font helpers
# ---------------------------------------------------------------------------
FONT_PATHS = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
]

_font_cache: dict = {}


def _find_system_font() -> str | None:
    for p in FONT_PATHS:
        if os.path.isfile(p):
            return p
    return None


def get_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    if size in _font_cache:
        return _font_cache[size]
    path = _find_system_font()
    if path:
        font = ImageFont.truetype(path, size)
    else:
        font = ImageFont.load_default()
    _font_cache[size] = font
    return font


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def draw_text_with_shadow(draw: ImageDraw.ImageDraw, xy, text: str, font,
                          fill=WHITE, shadow_fill=BLACK, offset=2):
    """Draw text with a dark shadow for readability on grayscale backgrounds."""
    x, y = xy
    # Draw shadow in 4 directions for thickness
    for dx, dy in [(-offset, -offset), (offset, -offset),
                   (-offset, offset), (offset, offset),
                   (0, offset), (0, -offset), (offset, 0), (-offset, 0)]:
        draw.text((x + dx, y + dy), text, font=font, fill=shadow_fill)
    draw.text((x, y), text, font=font, fill=fill)


def draw_badge(draw: ImageDraw.ImageDraw, xy, text: str, bg_color, font,
               text_color=WHITE, padding=8):
    """Draw a rounded-rectangle badge with text."""
    x, y = xy
    bbox = font.getbbox(text)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    rx0, ry0 = x, y
    rx1 = x + tw + padding * 2
    ry1 = y + th + padding * 2
    draw.rounded_rectangle([rx0, ry0, rx1, ry1], radius=6, fill=bg_color + (220,))
    draw.text((x + padding, y + padding - bbox[1]), text, font=font, fill=text_color)
    return rx1, ry1


def draw_ai_badge(draw: ImageDraw.ImageDraw, img_w: int):
    """Draw small 'AI' badge in NVIDIA green at top-right."""
    font = get_font(20)
    text = "AI"
    bbox = font.getbbox(text)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    pad = 6
    x = img_w - tw - pad * 2 - 10
    y = 10
    draw.rounded_rectangle([x, y, x + tw + pad * 2, y + th + pad * 2],
                           radius=4, fill=NVIDIA_GREEN + (230,))
    draw.text((x + pad, y + pad - bbox[1]), text, font=font, fill=WHITE)


def add_green_border(img: Image.Image, width: int = 2) -> Image.Image:
    """Add a thin NVIDIA green border around the image."""
    draw = ImageDraw.Draw(img)
    w, h = img.size
    for i in range(width):
        draw.rectangle([i, i, w - 1 - i, h - 1 - i], outline=NVIDIA_GREEN + (255,))
    return img


def load_source(path: str, target_size: int | None = None) -> Image.Image:
    """Load a grayscale source image and convert to RGBA."""
    img = Image.open(path).convert("RGBA")
    if target_size and img.size[0] != target_size:
        img = img.resize((target_size, target_size), Image.LANCZOS)
    return img


def save_annotated(img: Image.Image, name: str):
    """Save the annotated image as PNG."""
    # Flatten to RGB (drop alpha) then save
    out = Image.new("RGB", img.size, (0, 0, 0))
    out.paste(img, mask=img.split()[3])
    add_green_border(out)
    out.save(os.path.join(OUT_DIR, name), "PNG")
    print(f"  -> Saved {name} ({out.size[0]}x{out.size[1]})")


def draw_dashed_line(draw, start, end, fill, width=2, dash_len=12, gap_len=8):
    """Draw a dashed line between two points."""
    x0, y0 = start
    x1, y1 = end
    length = math.hypot(x1 - x0, y1 - y0)
    if length == 0:
        return
    dx = (x1 - x0) / length
    dy = (y1 - y0) / length
    pos = 0
    while pos < length:
        seg_end = min(pos + dash_len, length)
        draw.line([(x0 + dx * pos, y0 + dy * pos),
                    (x0 + dx * seg_end, y0 + dy * seg_end)],
                   fill=fill, width=width)
        pos = seg_end + gap_len


# ===========================================================================
# Image 1: Normal CXR
# ===========================================================================
def annotate_cxr_normal():
    print("1. fullres_cxr_synth_000 — Normal CXR")
    img = load_source(os.path.join(FULLRES_DIR, "fullres_cxr_synth_000.png"))
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    w, h = img.size

    # --- Measurement grid: CTR lines ---
    # Thorax line (full width at mid-height)
    thorax_y = int(h * 0.48)
    thorax_x0 = int(w * 0.18)
    thorax_x1 = int(w * 0.82)
    draw_dashed_line(draw, (thorax_x0, thorax_y), (thorax_x1, thorax_y),
                     fill=CYAN + (180,), width=2, dash_len=14, gap_len=8)
    # Small end-caps
    for x in [thorax_x0, thorax_x1]:
        draw.line([(x, thorax_y - 10), (x, thorax_y + 10)], fill=CYAN + (180,), width=2)

    # Heart line (~40% of thorax width, centred slightly left)
    heart_cx = int(w * 0.47)
    heart_half = int((thorax_x1 - thorax_x0) * 0.24)
    heart_y = int(h * 0.52)
    draw_dashed_line(draw, (heart_cx - heart_half, heart_y),
                     (heart_cx + heart_half, heart_y),
                     fill=NVIDIA_GREEN + (200,), width=2, dash_len=10, gap_len=6)
    for x in [heart_cx - heart_half, heart_cx + heart_half]:
        draw.line([(x, heart_y - 10), (x, heart_y + 10)],
                  fill=NVIDIA_GREEN + (200,), width=2)

    # CTR label near measurement
    font_sm = get_font(18)
    draw_text_with_shadow(draw, (heart_cx + heart_half + 8, heart_y - 10),
                          "CTR 0.48", font_sm, fill=NVIDIA_GREEN)

    # --- NORMAL badge top-left ---
    font_badge = get_font(28)
    draw_badge(draw, (20, 20), "NORMAL", NVIDIA_GREEN, font_badge)

    # --- Bottom label ---
    font_label = get_font(22)
    draw_text_with_shadow(draw, (20, h - 50),
                          "AI Analysis: No acute findings. CTR 0.48 (normal)",
                          font_label, fill=NVIDIA_GREEN)

    # --- AI badge ---
    draw_ai_badge(draw, w)

    img = Image.alpha_composite(img, overlay)
    save_annotated(img, "fullres_cxr_synth_000_annotated.png")


# ===========================================================================
# Image 2: Consolidation
# ===========================================================================
def annotate_cxr_consolidation():
    print("2. fullres_cxr_synth_001 — Consolidation")
    img = load_source(os.path.join(FULLRES_DIR, "fullres_cxr_synth_001.png"))
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    w, h = img.size

    # --- Semi-transparent red/orange region on left lower lobe ---
    # Approximate ellipse in lower-left lung field
    region_cx, region_cy = int(w * 0.62), int(h * 0.68)
    region_rx, region_ry = int(w * 0.14), int(h * 0.13)
    region_overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    rdraw = ImageDraw.Draw(region_overlay)
    rdraw.ellipse([region_cx - region_rx, region_cy - region_ry,
                   region_cx + region_rx, region_cy + region_ry],
                  fill=(255, 80, 20, 100))
    # Blur for heat-map feel
    region_overlay = region_overlay.filter(ImageFilter.GaussianBlur(radius=18))
    img = Image.alpha_composite(img, region_overlay)
    # Re-create overlay for crisp annotations
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)

    # Bounding box
    bx0 = region_cx - region_rx - 20
    by0 = region_cy - region_ry - 20
    bx1 = region_cx + region_rx + 20
    by1 = region_cy + region_ry + 20
    draw.rectangle([bx0, by0, bx1, by1], outline=RED + (220,), width=3)

    # Label above box
    font_label = get_font(22)
    draw_text_with_shadow(draw, (bx0, by0 - 30),
                          "Consolidation \u2014 85% confidence", font_label, fill=RED)

    # Measurement
    font_sm = get_font(18)
    draw_text_with_shadow(draw, (bx0, by1 + 6), "Area: 42 cm\u00b2", font_sm, fill=CYAN)

    # CRITICAL badge
    font_badge = get_font(28)
    draw_badge(draw, (20, 20), "CRITICAL", RED, font_badge)

    draw_ai_badge(draw, w)
    img = Image.alpha_composite(img, overlay)
    save_annotated(img, "fullres_cxr_synth_001_annotated.png")


# ===========================================================================
# Image 3: Pleural Effusion
# ===========================================================================
def annotate_cxr_effusion():
    print("3. fullres_cxr_synth_002 — Pleural Effusion")
    img = load_source(os.path.join(FULLRES_DIR, "fullres_cxr_synth_002.png"))
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    w, h = img.size

    # --- Semi-transparent blue/cyan meniscus along bottom-right lung ---
    # Build a polygon representing a meniscus shape
    meniscus_pts = []
    cx_right = int(w * 0.35)  # right lung is on left side of image (PA view)
    base_y = int(h * 0.88)
    top_y = int(h * 0.72)
    left_x = int(w * 0.15)
    right_x = int(w * 0.50)
    # Bottom edge
    meniscus_pts.append((left_x, base_y))
    meniscus_pts.append((right_x, base_y))
    # Curved top (meniscus)
    steps = 20
    for i in range(steps + 1):
        t = i / steps
        x = right_x + (left_x - right_x) * t
        # Parabolic curve: higher at edges, lower at centre
        offset = 4 * t * (1 - t)  # peaks at 0.5
        y = base_y - (base_y - top_y) * (0.3 + 0.7 * offset)
        meniscus_pts.append((x, y))

    region_overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    rdraw = ImageDraw.Draw(region_overlay)
    rdraw.polygon(meniscus_pts, fill=(0, 180, 255, 90))
    region_overlay = region_overlay.filter(ImageFilter.GaussianBlur(radius=10))
    img = Image.alpha_composite(img, region_overlay)

    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)

    # Label
    font_label = get_font(22)
    draw_text_with_shadow(draw, (left_x, top_y - 35),
                          "Pleural Effusion \u2014 91% confidence", font_label, fill=CYAN)

    # Depth measurement line
    meas_x = int(w * 0.33)
    draw.line([(meas_x, base_y), (meas_x, top_y + 20)], fill=CYAN + (200,), width=2)
    draw.line([(meas_x - 8, base_y), (meas_x + 8, base_y)], fill=CYAN + (200,), width=2)
    draw.line([(meas_x - 8, top_y + 20), (meas_x + 8, top_y + 20)],
              fill=CYAN + (200,), width=2)
    font_sm = get_font(18)
    draw_text_with_shadow(draw, (meas_x + 12, int((base_y + top_y) / 2)),
                          "Depth: 4.2 cm", font_sm, fill=CYAN)

    # Moderate severity badge
    font_badge = get_font(28)
    draw_badge(draw, (20, 20), "MODERATE", ORANGE, font_badge)

    draw_ai_badge(draw, w)
    img = Image.alpha_composite(img, overlay)
    save_annotated(img, "fullres_cxr_synth_002_annotated.png")


# ===========================================================================
# Image 4: Cardiomegaly
# ===========================================================================
def annotate_cxr_cardiomegaly():
    print("4. fullres_cxr_synth_003 — Cardiomegaly")
    img = load_source(os.path.join(FULLRES_DIR, "fullres_cxr_synth_003.png"))
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    w, h = img.size

    # --- Cardiac silhouette outline (ellipse) ---
    card_cx, card_cy = int(w * 0.48), int(h * 0.55)
    card_rx, card_ry = int(w * 0.22), int(h * 0.18)
    # Draw the outline as a thick ellipse border
    for offset in range(-2, 3):
        draw.ellipse([card_cx - card_rx + offset, card_cy - card_ry + offset,
                       card_cx + card_rx + offset, card_cy + card_ry + offset],
                      outline=(255, 100, 30, 180), width=2)

    # --- CTR measurement lines ---
    # Heart width line
    heart_y = card_cy
    heart_x0 = card_cx - card_rx
    heart_x1 = card_cx + card_rx
    draw.line([(heart_x0, heart_y), (heart_x1, heart_y)],
              fill=RED + (220,), width=3)
    for x in [heart_x0, heart_x1]:
        draw.line([(x, heart_y - 12), (x, heart_y + 12)], fill=RED + (220,), width=2)

    # Thorax width line
    thorax_y = int(h * 0.48)
    thorax_x0 = int(w * 0.16)
    thorax_x1 = int(w * 0.84)
    draw_dashed_line(draw, (thorax_x0, thorax_y), (thorax_x1, thorax_y),
                     fill=CYAN + (180,), width=2)
    for x in [thorax_x0, thorax_x1]:
        draw.line([(x, thorax_y - 10), (x, thorax_y + 10)], fill=CYAN + (180,), width=2)

    # Labels near lines
    font_sm = get_font(18)
    draw_text_with_shadow(draw, (heart_x1 + 10, heart_y - 10),
                          "Heart", font_sm, fill=RED)
    draw_text_with_shadow(draw, (thorax_x1 + 6, thorax_y - 10),
                          "Thorax", font_sm, fill=CYAN)

    # CTR label
    font_label = get_font(22)
    draw_text_with_shadow(draw, (20, h - 80),
                          "Cardiomegaly \u2014 CTR 0.62 (>0.50)", font_label, fill=RED)
    draw_text_with_shadow(draw, (20, h - 50),
                          "Cardiac-to-Thoracic Ratio exceeds normal threshold",
                          get_font(18), fill=ORANGE)

    # URGENT badge
    font_badge = get_font(28)
    draw_badge(draw, (20, 20), "URGENT", ORANGE, font_badge)

    draw_ai_badge(draw, w)
    img = Image.alpha_composite(img, overlay)
    save_annotated(img, "fullres_cxr_synth_003_annotated.png")


# ===========================================================================
# Image 5: Pneumothorax
# ===========================================================================
def annotate_cxr_pneumothorax():
    print("5. fullres_cxr_synth_004 — Pneumothorax")
    img = load_source(os.path.join(FULLRES_DIR, "fullres_cxr_synth_004.png"))
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    w, h = img.size

    # --- Visceral pleural line on left side (right side of image in PA view) ---
    # Bright thin line
    pleural_pts = []
    line_x_base = int(w * 0.72)
    for i in range(20):
        t = i / 19
        y = int(h * 0.15 + t * h * 0.50)
        # Slight curvature
        x = line_x_base + int(15 * math.sin(t * math.pi))
        pleural_pts.append((x, y))

    for i in range(len(pleural_pts) - 1):
        draw.line([pleural_pts[i], pleural_pts[i + 1]],
                  fill=WHITE + (240,), width=3)

    # --- Semi-transparent cyan shading between pleural line and chest wall ---
    # Build polygon: pleural line + chest wall edge
    wall_pts = []
    wall_x = int(w * 0.88)
    for pt in reversed(pleural_pts):
        wall_pts.append((wall_x, pt[1]))
    shade_polygon = pleural_pts + wall_pts

    shade_overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    sdraw = ImageDraw.Draw(shade_overlay)
    sdraw.polygon(shade_polygon, fill=(0, 212, 255, 80))
    shade_overlay = shade_overlay.filter(ImageFilter.GaussianBlur(radius=6))
    img = Image.alpha_composite(img, shade_overlay)

    overlay2 = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw2 = ImageDraw.Draw(overlay2)
    # Re-draw the pleural line on top
    for i in range(len(pleural_pts) - 1):
        draw2.line([pleural_pts[i], pleural_pts[i + 1]],
                   fill=WHITE + (240,), width=3)

    # Depth measurement at apex
    apex_pt = pleural_pts[0]
    meas_end = (wall_x, apex_pt[1])
    draw2.line([apex_pt, meas_end], fill=CYAN + (200,), width=2)
    draw2.line([(apex_pt[0], apex_pt[1] - 6), (apex_pt[0], apex_pt[1] + 6)],
               fill=CYAN + (200,), width=2)
    draw2.line([(meas_end[0], meas_end[1] - 6), (meas_end[0], meas_end[1] + 6)],
               fill=CYAN + (200,), width=2)
    font_sm = get_font(18)
    draw_text_with_shadow(draw2, (apex_pt[0], apex_pt[1] - 25),
                          "Depth: 2.8 cm at apex", font_sm, fill=CYAN)

    # Label
    font_label = get_font(22)
    draw_text_with_shadow(draw2, (20, h - 50),
                          "Pneumothorax \u2014 94% confidence", font_label, fill=RED)

    # CRITICAL badge
    font_badge = get_font(28)
    draw_badge(draw2, (20, 20), "CRITICAL", RED, font_badge)

    draw_ai_badge(draw2, w)

    # Composite the original overlay too
    img = Image.alpha_composite(img, overlay)
    img = Image.alpha_composite(img, overlay2)
    save_annotated(img, "fullres_cxr_synth_004_annotated.png")


# ===========================================================================
# Image 6: CT Head Hemorrhage
# ===========================================================================
def annotate_ct_head_hemorrhage():
    print("6. organ_ct_000 — CT Head Hemorrhage")
    img = load_source(os.path.join(SRC_DIR, "organ_ct_000.png"), target_size=512)
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    w, h = img.size

    # --- Red/orange elliptical overlay in right basal ganglia region ---
    hem_cx, hem_cy = int(w * 0.40), int(h * 0.45)
    hem_rx, hem_ry = 40, 35

    region_overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    rdraw = ImageDraw.Draw(region_overlay)
    rdraw.ellipse([hem_cx - hem_rx, hem_cy - hem_ry,
                   hem_cx + hem_rx, hem_cy + hem_ry],
                  fill=(255, 70, 20, 100))
    region_overlay = region_overlay.filter(ImageFilter.GaussianBlur(radius=12))
    img = Image.alpha_composite(img, region_overlay)

    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)

    # Crosshair measurement lines
    draw.line([(hem_cx - 55, hem_cy), (hem_cx + 55, hem_cy)],
              fill=CYAN + (200,), width=2)
    draw.line([(hem_cx, hem_cy - 50), (hem_cx, hem_cy + 50)],
              fill=CYAN + (200,), width=2)
    # Small circle at center
    draw.ellipse([hem_cx - 3, hem_cy - 3, hem_cx + 3, hem_cy + 3],
                 fill=CYAN + (220,))

    # Ellipse outline
    draw.ellipse([hem_cx - hem_rx, hem_cy - hem_ry,
                  hem_cx + hem_rx, hem_cy + hem_ry],
                 outline=RED + (200,), width=2)

    # Label
    font_label = get_font(18)
    draw_text_with_shadow(draw, (20, h - 70),
                          "ICH \u2014 Volume: 28.5 mL", font_label, fill=RED)
    draw_text_with_shadow(draw, (20, h - 45),
                          "Midline Shift: 4.8 mm", font_label, fill=ORANGE)

    # CRITICAL badge
    font_badge = get_font(22)
    draw_badge(draw, (15, 15), "CRITICAL", RED, font_badge)

    draw_ai_badge(draw, w)
    img = Image.alpha_composite(img, overlay)
    save_annotated(img, "ct_head_hemorrhage_annotated.png")


# ===========================================================================
# Image 7: CT Chest Nodule
# ===========================================================================
def annotate_ct_chest_nodule():
    print("7. organ_ct_002 — CT Chest Nodule")
    img = load_source(os.path.join(SRC_DIR, "organ_ct_002.png"), target_size=512)
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    w, h = img.size

    # --- Small circle with crosshairs in upper right area ---
    nod_cx, nod_cy = int(w * 0.35), int(h * 0.35)
    nod_r = 22

    draw.ellipse([nod_cx - nod_r, nod_cy - nod_r,
                  nod_cx + nod_r, nod_cy + nod_r],
                 outline=RED + (220,), width=3)
    # Crosshairs extending beyond circle
    ext = 20
    draw.line([(nod_cx - nod_r - ext, nod_cy), (nod_cx + nod_r + ext, nod_cy)],
              fill=CYAN + (180,), width=1)
    draw.line([(nod_cx, nod_cy - nod_r - ext), (nod_cx, nod_cy + nod_r + ext)],
              fill=CYAN + (180,), width=1)

    # Labels
    font_label = get_font(18)
    draw_text_with_shadow(draw, (nod_cx + nod_r + ext + 5, nod_cy - 10),
                          "18 mm", font_label, fill=CYAN)

    draw_text_with_shadow(draw, (20, h - 70),
                          "Part-Solid Nodule \u2014 18mm, Lung-RADS 4B",
                          font_label, fill=RED)
    font_sm = get_font(16)
    draw_text_with_shadow(draw, (20, h - 45),
                          "VDT: 245 days", font_sm, fill=ORANGE)

    # CRITICAL badge
    font_badge = get_font(22)
    draw_badge(draw, (15, 15), "CRITICAL", RED, font_badge)

    draw_ai_badge(draw, w)
    img = Image.alpha_composite(img, overlay)
    save_annotated(img, "ct_chest_nodule_annotated.png")


# ===========================================================================
# Image 8: CT Coronary
# ===========================================================================
def annotate_ct_coronary():
    """DISABLED — this drew coronary anatomy onto a slice of SPINE.

    The source, `dicom_ct_sample_preview.png`, is a pydicom public-domain axial CT slice of a
    VERTEBRA: vertebral body, pedicles, lamina, spinous process, spinal canal. This function drew
    stylised LAD/LCx/RCA paths, stenosis markers and an Agatston overlay on top of it and saved the
    result as `ct_coronary_annotated.png`, which the demo captioned "CTA Coronary — AI Stenosis
    Analysis". Any radiologist or cardiologist identifies that in about a second.

    The coronary panel is now produced by `scripts/render_coronary_mesh.py` from real coronary
    artery meshes (`data/cardiac_ct/CoronariesNC6/`). Do not re-enable this without a genuine
    cardiac source image.
    """
    print("8. CT Coronary — SKIPPED (source is a vertebra; see render_coronary_mesh.py)")
    return


def _disabled_annotate_ct_coronary_vertebra_source():
    print("8. dicom_ct_sample_preview — CT Coronary")
    img = load_source(os.path.join(SRC_DIR, "dicom_ct_sample_preview.png"),
                      target_size=512)
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    w, h = img.size

    # --- Stylised coronary artery paths ---
    # LAD (Left Anterior Descending) — runs from upper-center downward-left
    lad_pts = [(int(w * 0.48), int(h * 0.32)),
               (int(w * 0.46), int(h * 0.40)),
               (int(w * 0.44), int(h * 0.48)),
               (int(w * 0.43), int(h * 0.56)),
               (int(w * 0.42), int(h * 0.64)),
               (int(w * 0.41), int(h * 0.70))]

    # LCx (Left Circumflex) — curves leftward
    lcx_pts = [(int(w * 0.48), int(h * 0.32)),
               (int(w * 0.52), int(h * 0.38)),
               (int(w * 0.56), int(h * 0.44)),
               (int(w * 0.58), int(h * 0.50)),
               (int(w * 0.58), int(h * 0.56))]

    # RCA (Right Coronary Artery) — curves rightward then down
    rca_pts = [(int(w * 0.45), int(h * 0.30)),
               (int(w * 0.40), int(h * 0.34)),
               (int(w * 0.36), int(h * 0.40)),
               (int(w * 0.34), int(h * 0.48)),
               (int(w * 0.35), int(h * 0.56)),
               (int(w * 0.37), int(h * 0.62))]

    # Draw artery paths
    for pts, color, label in [(lad_pts, NVIDIA_GREEN, "LAD"),
                               (lcx_pts, NVIDIA_GREEN, "LCx"),
                               (rca_pts, NVIDIA_GREEN, "RCA")]:
        for i in range(len(pts) - 1):
            draw.line([pts[i], pts[i + 1]], fill=color + (180,), width=3)
        # Label at end of path
        font_sm = get_font(14)
        ex, ey = pts[-1]
        draw_text_with_shadow(draw, (ex + 6, ey - 8), label, font_sm,
                              fill=NVIDIA_GREEN)

    # --- Stenosis marker on LAD (between points 1 and 2) ---
    stenosis_pt = lad_pts[1]
    sx, sy = stenosis_pt
    # Red marker
    draw.ellipse([sx - 8, sy - 8, sx + 8, sy + 8], fill=RED + (220,))
    draw.ellipse([sx - 10, sy - 10, sx + 10, sy + 10],
                 outline=RED + (255,), width=2)
    # Arrow line pointing to label
    font_label = get_font(16)
    draw.line([(sx + 10, sy), (sx + 35, sy - 20)], fill=RED + (220,), width=2)
    draw_text_with_shadow(draw, (sx + 38, sy - 30),
                          "72% Stenosis", font_label, fill=RED)

    # Main label at bottom
    font_main = get_font(17)
    # 77.9% is the measured value from coronary_analysis.json; this label read 72% while every
    # other surface said 77.9%. Keep it in step with the analysis if that number ever changes.
    draw_text_with_shadow(draw, (15, h - 70),
                          "LAD Proximal: 77.9% Stenosis \u2014 CAD-RADS 4A/P3/HRP",
                          font_main, fill=RED)

    # Calcium score badge. Agatston needs Hounsfield units a surface mesh does not carry, so the
    # score is representative; an exact MESA percentile would also need race/ethnicity.
    font_badge = get_font(18)
    draw_badge(draw, (15, h - 40), "Ca Score: 385 (repr.)",
               ORANGE, font_badge)

    # CRITICAL badge top-left
    font_crit = get_font(22)
    draw_badge(draw, (15, 15), "CRITICAL", RED, font_crit)

    draw_ai_badge(draw, w)
    img = Image.alpha_composite(img, overlay)
    save_annotated(img, "ct_coronary_annotated.png")


# ===========================================================================
# Main
# ===========================================================================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"Source directory : {SRC_DIR}")
    print(f"Full-res source : {FULLRES_DIR}")
    print(f"Output directory : {OUT_DIR}")
    print("=" * 60)

    annotate_cxr_normal()
    annotate_cxr_consolidation()
    annotate_cxr_effusion()
    annotate_cxr_cardiomegaly()
    annotate_cxr_pneumothorax()
    annotate_ct_head_hemorrhage()
    annotate_ct_chest_nodule()
    annotate_ct_coronary()

    print("=" * 60)
    print(f"Done. {8} annotated images written to:\n  {OUT_DIR}")


if __name__ == "__main__":
    main()
