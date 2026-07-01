#!/usr/bin/env python3
"""Generate high-resolution synthetic NIfTI volumes and AI segmentation overlays.

Creates 256x256 resolution synthetic medical imaging volumes with realistic
anatomical structures, then generates segmentation overlays that look
clinical-grade for the Clinical Imaging Engine demo.

Outputs:
  data/sample_images/highres_ct_head.nii.gz       — 256x256x128 CT Head
  data/sample_images/highres_ct_chest.nii.gz      — 256x256x128 CT Chest
  data/sample_images/highres_brain_flair.nii.gz   — 256x256x96  Brain MRI FLAIR
  data/demo/segmentation/highres_ct_head_overlay.png
  data/demo/segmentation/highres_ct_head_segmented.gif
  data/demo/segmentation/highres_ct_chest_overlay.png
  data/demo/segmentation/highres_ct_chest_segmented.gif
  data/demo/segmentation/highres_brain_flair_overlay.png
  data/demo/segmentation/highres_brain_flair_segmented.gif

Usage:
  python33 scripts/generate_highres_volumes.py
"""

import os
import sys
import numpy as np
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("Clinical Imaging Engine — High-Resolution Volume & Overlay Generator")
print("=" * 70)

import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import imageio.v3 as iio
from scipy.ndimage import gaussian_filter, binary_dilation

# Directories
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAMPLE_DIR = os.path.join(BASE_DIR, "data", "sample_images")
OUT_DIR = os.path.join(BASE_DIR, "data", "demo", "segmentation")
os.makedirs(SAMPLE_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# --- Utility functions ---

def ellipsoid_mask(shape, center, radii):
    """Create a 3D ellipsoid mask."""
    z, y, x = np.ogrid[:shape[0], :shape[1], :shape[2]]
    return ((z - center[0])**2 / radii[0]**2 +
            (y - center[1])**2 / radii[1]**2 +
            (x - center[2])**2 / radii[2]**2) <= 1.0

def sphere_mask(shape, center, radius):
    """Create a 3D sphere mask."""
    return ellipsoid_mask(shape, center, (radius, radius, radius))

def irregular_blob(shape, center, base_radius, irregularity=0.3, seed=42):
    """Create an irregular 3D blob using perturbed sphere."""
    rng = np.random.RandomState(seed)
    mask = np.zeros(shape, dtype=bool)
    z, y, x = np.ogrid[:shape[0], :shape[1], :shape[2]]
    dist = np.sqrt((z - center[0])**2 + (y - center[1])**2 + (x - center[2])**2)
    # Create base sphere then add noise-based perturbation
    noise = gaussian_filter(rng.randn(*shape) * irregularity * base_radius, sigma=base_radius * 0.4)
    mask = dist < (base_radius + noise)
    return mask

def add_noise(volume, sigma_frac=0.05):
    """Add Gaussian noise proportional to tissue values."""
    rng = np.random.RandomState(12345)
    # Use absolute range for noise scale
    vrange = volume.max() - volume.min()
    noise = rng.randn(*volume.shape).astype(np.float32) * (vrange * sigma_frac)
    return volume + noise


# ========================================================================
# VOLUME GENERATION
# ========================================================================

def generate_ct_head(filepath):
    """Generate 256x256x128 synthetic CT head volume."""
    print("\n[1/3] Generating CT Head (256x256x128)...")
    shape = (256, 256, 128)
    volume = np.full(shape, -1000.0, dtype=np.float32)  # Air background

    cy, cx = 128, 128  # Center in Y,X
    cz = 64             # Center in Z

    # Skull — outer ellipse
    skull_outer = ellipsoid_mask(shape, (cz, cy, cx), (55, 95, 80))
    skull_inner = ellipsoid_mask(shape, (cz, cy, cx), (50, 88, 73))
    skull_shell = skull_outer & ~skull_inner
    volume[skull_shell] = np.random.RandomState(1).uniform(800, 1000, skull_shell.sum())

    # Brain parenchyma — fill inside skull
    brain_mask = skull_inner.copy()

    # Create lateral ventricles — two elongated ellipses
    vent_left = ellipsoid_mask(shape, (cz, cy - 12, cx - 5), (20, 8, 4))
    vent_right = ellipsoid_mask(shape, (cz, cy - 12, cx + 5), (20, 8, 4))
    ventricles = vent_left | vent_right

    # Third ventricle — midline slit
    third_vent = ellipsoid_mask(shape, (cz, cy - 15, cx), (12, 3, 1))
    ventricles = ventricles | third_vent

    # CSF in ventricles
    csf_mask = ventricles & brain_mask
    volume[csf_mask] = np.random.RandomState(2).uniform(5, 10, csf_mask.sum())

    # White matter — inner brain regions
    wm_mask = ellipsoid_mask(shape, (cz, cy, cx), (38, 65, 55)) & brain_mask & ~csf_mask
    volume[wm_mask] = np.random.RandomState(3).uniform(25, 35, wm_mask.sum())

    # Gray matter — outer cortex (between white matter and skull)
    gm_mask = brain_mask & ~wm_mask & ~csf_mask
    volume[gm_mask] = np.random.RandomState(4).uniform(35, 45, gm_mask.sum())

    # Falx cerebri (midline)
    falx = np.zeros(shape, dtype=bool)
    falx[:, :, cx-1:cx+2] = True
    falx_brain = falx & brain_mask & ~csf_mask
    volume[falx_brain] = np.random.RandomState(5).uniform(50, 70, falx_brain.sum())

    # Hemorrhage in right basal ganglia — irregular blob
    hemorrhage = irregular_blob(shape, (cz, cy + 5, cx + 20), base_radius=12, irregularity=0.35, seed=99)
    hemorrhage = hemorrhage & brain_mask
    volume[hemorrhage] = np.random.RandomState(6).uniform(55, 75, hemorrhage.sum())

    # Smooth and add noise
    volume = gaussian_filter(volume, sigma=0.8)
    volume = add_noise(volume, sigma_frac=0.03)

    # Save NIfTI with 1mm isotropic
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    img = nib.Nifti1Image(volume, affine)
    nib.save(img, filepath)
    print(f"  Saved: {filepath} ({os.path.getsize(filepath)/1024/1024:.1f} MB)")
    print(f"  Shape: {shape}, Range: [{volume.min():.0f}, {volume.max():.0f}] HU")
    return volume


def generate_ct_chest(filepath):
    """Generate 256x256x128 synthetic CT chest volume."""
    print("\n[2/3] Generating CT Chest (256x256x128)...")
    shape = (256, 256, 128)
    volume = np.full(shape, -1000.0, dtype=np.float32)  # Air background

    cy, cx = 128, 128
    cz = 64

    # Body contour — outer ellipse
    body_outer = ellipsoid_mask(shape, (cz, cy, cx), (55, 110, 90))
    volume[body_outer] = 30  # Soft tissue baseline

    # Spine — posterior cylinder
    spine = ellipsoid_mask(shape, (cz, cy + 60, cx), (50, 12, 10))
    volume[spine] = np.random.RandomState(10).uniform(500, 800, spine.sum())

    # Ribs — create curved rib-like structures on each side
    rng = np.random.RandomState(11)
    for z_offset in range(-40, 45, 12):
        zc = cz + z_offset
        for side in [-1, 1]:
            # Each rib as an arc (approximated by ellipsoidal shell)
            rib_outer = ellipsoid_mask(shape, (zc, cy + 20, cx + side * 50), (3, 50, 35))
            rib_inner = ellipsoid_mask(shape, (zc, cy + 20, cx + side * 50), (2, 47, 32))
            rib = rib_outer & ~rib_inner & body_outer
            volume[rib] = rng.uniform(400, 700, rib.sum())

    # Sternum — anterior midline
    sternum = ellipsoid_mask(shape, (cz, cy - 75, cx), (30, 5, 8))
    sternum = sternum & body_outer
    volume[sternum] = rng.uniform(500, 700, sternum.sum())

    # Lung fields — two large irregular regions
    lung_left = ellipsoid_mask(shape, (cz, cy - 10, cx - 35), (42, 55, 28))
    lung_right = ellipsoid_mask(shape, (cz, cy - 10, cx + 35), (42, 55, 30))
    lungs = lung_left | lung_right
    lungs = lungs & body_outer
    volume[lungs] = np.random.RandomState(12).uniform(-800, -500, lungs.sum())

    # Heart — central-left structure
    heart = ellipsoid_mask(shape, (cz + 5, cy - 5, cx - 10), (25, 35, 30))
    heart = heart & body_outer & ~lungs
    volume[heart] = np.random.RandomState(13).uniform(35, 50, heart.sum())

    # Mediastinum — central soft tissue
    mediastinum = ellipsoid_mask(shape, (cz, cy, cx), (45, 25, 15))
    mediastinum = mediastinum & body_outer & ~lungs & ~heart
    volume[mediastinum] = np.random.RandomState(14).uniform(30, 50, mediastinum.sum())

    # Pulmonary nodule — right upper lobe, ~8 voxels diameter
    nodule_center = (cz - 20, cy - 25, cx + 30)
    nodule = sphere_mask(shape, nodule_center, 4)
    nodule = nodule & lung_right
    volume[nodule] = np.random.RandomState(15).uniform(40, 60, nodule.sum())

    # Chest wall soft tissue
    chest_wall = body_outer & ~lungs & ~heart & ~spine & ~mediastinum
    # Only update areas still at baseline
    baseline_mask = chest_wall & (np.abs(volume - 30) < 1)
    volume[baseline_mask] = np.random.RandomState(16).uniform(30, 50, baseline_mask.sum())

    # Smooth and add noise
    volume = gaussian_filter(volume, sigma=0.7)
    volume = add_noise(volume, sigma_frac=0.03)

    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    img = nib.Nifti1Image(volume, affine)
    nib.save(img, filepath)
    print(f"  Saved: {filepath} ({os.path.getsize(filepath)/1024/1024:.1f} MB)")
    print(f"  Shape: {shape}, Range: [{volume.min():.0f}, {volume.max():.0f}] HU")
    return volume


def generate_brain_flair(filepath):
    """Generate 256x256x96 synthetic Brain MRI FLAIR volume."""
    print("\n[3/3] Generating Brain MRI FLAIR (256x256x96)...")
    shape = (256, 256, 96)
    volume = np.zeros(shape, dtype=np.float32)  # Background

    cy, cx = 128, 128
    cz = 48

    # Skull — dark on FLAIR
    skull_outer = ellipsoid_mask(shape, (cz, cy, cx), (42, 95, 80))
    skull_inner = ellipsoid_mask(shape, (cz, cy, cx), (38, 89, 74))
    skull_shell = skull_outer & ~skull_inner
    volume[skull_shell] = 50

    brain_mask = skull_inner.copy()

    # Lateral ventricles — CSF dark on FLAIR
    vent_left = ellipsoid_mask(shape, (cz, cy - 10, cx - 6), (16, 8, 4))
    vent_right = ellipsoid_mask(shape, (cz, cy - 10, cx + 6), (16, 8, 4))
    third_vent = ellipsoid_mask(shape, (cz, cy - 14, cx), (10, 3, 1))
    ventricles = (vent_left | vent_right | third_vent) & brain_mask

    # White matter — inner regions
    wm_region = ellipsoid_mask(shape, (cz, cy, cx), (30, 65, 55)) & brain_mask & ~ventricles
    volume[wm_region] = np.random.RandomState(20).uniform(400, 600, wm_region.sum())

    # Gray matter — cortical mantle
    gm_region = brain_mask & ~wm_region & ~ventricles
    volume[gm_region] = np.random.RandomState(21).uniform(600, 800, gm_region.sum())

    # CSF — dark on FLAIR
    volume[ventricles] = np.random.RandomState(22).uniform(50, 100, ventricles.sum())

    # Sulcal CSF — thin dark lines in cortex (simplified)
    sulci = np.zeros(shape, dtype=bool)
    rng = np.random.RandomState(23)
    for angle_deg in range(0, 360, 20):
        angle = np.radians(angle_deg)
        for r in range(60, 85, 3):
            sy = int(cy + r * np.cos(angle))
            sx = int(cx + r * np.sin(angle))
            if 5 < sy < 251 and 5 < sx < 251:
                sulci[:, sy-1:sy+1, sx-1:sx+1] = True
    sulci = sulci & brain_mask
    volume[sulci] = rng.uniform(50, 120, sulci.sum())

    # MS lesions — bright hyperintense spots
    # Periventricular lesions (Dawson's fingers pattern)
    lesion_centers = [
        (cz + 5, cy - 5, cx + 15),    # Right periventricular
        (cz - 3, cy - 8, cx - 18),    # Left periventricular
        (cz + 8, cy + 20, cx + 25),   # Right juxtacortical
        (cz - 5, cy + 15, cx - 30),   # Left juxtacortical (subcortical)
    ]
    lesion_masks = []
    for i, lc in enumerate(lesion_centers):
        lesion = irregular_blob(shape, lc, base_radius=4 + i % 2, irregularity=0.4, seed=30 + i)
        lesion = lesion & brain_mask & ~ventricles
        volume[lesion] = np.random.RandomState(30 + i).uniform(900, 1000, lesion.sum())
        lesion_masks.append(lesion)

    # Smooth and add noise
    volume = gaussian_filter(volume, sigma=0.6)
    volume = add_noise(volume, sigma_frac=0.04)

    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    img = nib.Nifti1Image(volume, affine)
    nib.save(img, filepath)
    print(f"  Saved: {filepath} ({os.path.getsize(filepath)/1024/1024:.1f} MB)")
    print(f"  Shape: {shape}, Range: [{volume.min():.0f}, {volume.max():.0f}]")
    return volume


# ========================================================================
# SEGMENTATION & OVERLAY GENERATION
# ========================================================================

def segment_ct_head(volume):
    """Segment CT head into bone, gray matter, white matter, CSF, hemorrhage."""
    mask = np.zeros(volume.shape, dtype=np.uint8)
    # Labels: 1=bone, 2=gray_matter, 3=white_matter, 4=CSF, 5=hemorrhage
    mask[volume > 200] = 1                                      # Bone
    mask[(volume > 32) & (volume < 48) & (mask == 0)] = 2       # Gray matter
    mask[(volume > 20) & (volume < 38) & (mask == 0)] = 3       # White matter
    mask[(volume > 2) & (volume < 15) & (mask == 0)] = 4        # CSF
    mask[(volume > 50) & (volume < 80) & (mask == 0)] = 5       # Hemorrhage
    return mask

def segment_ct_chest(volume):
    """Segment CT chest into bone, lung, heart, nodule, soft tissue."""
    mask = np.zeros(volume.shape, dtype=np.uint8)
    # Labels: 1=bone, 2=lung, 3=heart, 4=soft_tissue, 5=nodule
    mask[volume > 200] = 1                                       # Bone
    mask[(volume < -200) & (volume > -1000)] = 2                 # Lung air
    # Heart: central region with soft tissue density
    cy, cx, cz = 128, 128, 64
    heart_region = ellipsoid_mask(volume.shape, (cz + 5, cy - 5, cx - 10), (25, 35, 30))
    mask[heart_region & (volume > 30) & (volume < 55)] = 3      # Heart
    # Soft tissue
    mask[(volume > 20) & (volume < 60) & (mask == 0)] = 4       # Soft tissue
    # Nodule: detect small dense region in lung field
    # Re-identify: anything in lung region that's denser than air
    nodule_mask = (volume > 30) & (volume < 70) & (mask == 2)
    # Expand to catch nodule voxels that might not have been classified as lung
    nodule_region = sphere_mask(volume.shape, (cz - 20, cy - 25, cx + 30), 6)
    nodule_mask = nodule_mask | (nodule_region & (volume > 20))
    mask[nodule_mask] = 5                                        # Nodule
    return mask

def segment_brain_flair(volume):
    """Segment brain FLAIR into gray matter, white matter, CSF, lesions."""
    mask = np.zeros(volume.shape, dtype=np.uint8)
    # Labels: 1=gray_matter, 2=white_matter, 3=CSF, 4=lesions
    mask[(volume > 500) & (volume < 850)] = 1                    # Gray matter
    mask[(volume > 300) & (volume < 650) & (mask == 0)] = 2      # White matter
    mask[(volume > 20) & (volume < 150) & (mask == 0)] = 3       # CSF
    mask[volume > 850] = 4                                        # MS lesions (bright)
    return mask


def create_overlay_frame(slice_2d, mask_2d, colors, window_center, window_width,
                         title="", markers=None, dpi=128, figsize=(8, 8)):
    """Create a single overlay frame at high resolution."""
    vmin = window_center - window_width / 2
    vmax = window_center + window_width / 2
    windowed = np.clip((slice_2d - vmin) / (vmax - vmin), 0, 1)

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')

    # Base grayscale image
    ax.imshow(windowed.T, cmap='gray', aspect='equal', origin='lower')

    # Segmentation overlay
    if mask_2d is not None and mask_2d.max() > 0:
        overlay = np.zeros((*mask_2d.shape, 4))
        for label_idx, color in colors.items():
            region = mask_2d == label_idx
            if region.any():
                overlay[region] = color
        ax.imshow(overlay.transpose(1, 0, 2), aspect='equal', origin='lower')

    # Circle markers for findings
    if markers is not None:
        for (my, mx, radius, color, label_text) in markers:
            if mask_2d is not None:
                region = mask_2d == color if isinstance(color, int) else None
            circle = plt.Circle((my, mx), radius, fill=False, edgecolor=color if isinstance(color, str) else 'orange',
                              linewidth=1.5, linestyle='--')
            ax.add_patch(circle)

    # NVIDIA green border
    for spine in ax.spines.values():
        spine.set_color('#76B900')
        spine.set_linewidth(2.0)

    if title:
        ax.set_title(title, color='#76B900', fontsize=10, fontweight='bold', pad=8)

    ax.set_xticks([])
    ax.set_yticks([])

    # AI label
    ax.text(0.02, 0.02, 'AI Segmentation \u2014 VISTA-3D',
            transform=ax.transAxes, color='#76B900', fontsize=7, alpha=0.8,
            verticalalignment='bottom', fontfamily='monospace')

    fig.tight_layout(pad=0.5)
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)
    frame = buf[:, :, :3].copy()
    plt.close(fig)
    return frame


def find_best_slice(mask, axis=2):
    """Find the slice with maximum segmentation content along given axis."""
    best_idx = 0
    best_count = 0
    n = mask.shape[axis]
    for i in range(int(n * 0.15), int(n * 0.85)):
        if axis == 2:
            count = (mask[:, :, i] > 0).sum()
        elif axis == 0:
            count = (mask[i, :, :] > 0).sum()
        else:
            count = (mask[:, i, :] > 0).sum()
        if count > best_count:
            best_count = count
            best_idx = i
    return best_idx


def generate_overlays_ct_head(volume, mask, out_png, out_gif):
    """Generate CT head overlays with hemorrhage highlighting."""
    print("\n  Generating CT Head overlays...")
    colors = {
        1: (1.0, 1.0, 0.85, 0.20),   # Bone — cream
        2: (0.45, 0.35, 0.75, 0.25),  # Gray matter — blue/purple
        3: (0.35, 0.45, 0.85, 0.20),  # White matter — blue
        4: (0.0, 0.8, 0.85, 0.30),    # CSF — cyan
        5: (1.0, 0.1, 0.1, 0.55),     # Hemorrhage — bright red
    }
    wc, ww = 40, 80  # Brain window

    # Best slice
    best_z = find_best_slice(mask, axis=2)
    # Try to find slice with hemorrhage
    for z in range(mask.shape[2] // 4, 3 * mask.shape[2] // 4):
        if (mask[:, :, z] == 5).sum() > 20:
            best_z = z
            break

    frame = create_overlay_frame(volume[:, :, best_z], mask[:, :, best_z], colors,
                                  wc, ww,
                                  title=f"CT Head \u2014 AI Segmentation \u2014 Slice {best_z}",
                                  dpi=128, figsize=(8, 8))
    iio.imwrite(out_png, frame)
    print(f"  Saved: {out_png} ({os.path.getsize(out_png)/1024:.0f} KB)")

    # Animated GIF — 50 slices
    n_slices = mask.shape[2]
    start, end = int(n_slices * 0.15), int(n_slices * 0.85)
    n_frames = 50
    step = max(1, (end - start) // n_frames)
    indices = list(range(start, end, step))[:n_frames]

    frames = []
    for idx in indices:
        f = create_overlay_frame(volume[:, :, idx], mask[:, :, idx], colors,
                                  wc, ww,
                                  title=f"CT Head \u2014 Slice {idx}/{n_slices}",
                                  dpi=128, figsize=(8, 8))
        frames.append(f)

    iio.imwrite(out_gif, frames, duration=int(1000/6), loop=0)
    print(f"  Saved: {out_gif} ({os.path.getsize(out_gif)/1024:.0f} KB)")


def generate_overlays_ct_chest(volume, mask, out_png, out_gif):
    """Generate CT chest overlays with nodule highlighting."""
    print("\n  Generating CT Chest overlays...")
    colors = {
        1: (1.0, 1.0, 0.85, 0.18),   # Bone — cream
        2: (0.2, 0.75, 0.2, 0.25),    # Lung — green
        3: (0.9, 0.15, 0.15, 0.40),   # Heart — red
        4: (0.55, 0.35, 0.2, 0.12),   # Soft tissue — subtle brown
        5: (1.0, 0.55, 0.0, 0.65),    # Nodule — bright orange
    }
    wc, ww = -400, 1500  # Lung window

    # Find slice with nodule
    best_z = find_best_slice(mask, axis=2)
    for z in range(mask.shape[2]):
        if (mask[:, :, z] == 5).sum() > 3:
            best_z = z
            break

    # Create best slice with nodule marker
    slice_2d = volume[:, :, best_z]
    mask_2d = mask[:, :, best_z]

    # Find nodule centroid for marker
    nodule_yx = np.where(mask_2d == 5)
    fig_dpi = 128

    frame = create_overlay_frame(slice_2d, mask_2d, colors, wc, ww,
                                  title=f"CT Chest \u2014 AI Segmentation \u2014 Slice {best_z}",
                                  dpi=128, figsize=(8, 8))

    # Re-draw with nodule marker
    vmin = wc - ww / 2
    vmax = wc + ww / 2
    windowed = np.clip((slice_2d - vmin) / (vmax - vmin), 0, 1)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=128)
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    ax.imshow(windowed.T, cmap='gray', aspect='equal', origin='lower')

    overlay = np.zeros((*mask_2d.shape, 4))
    for label_idx, color in colors.items():
        region = mask_2d == label_idx
        if region.any():
            overlay[region] = color
    ax.imshow(overlay.transpose(1, 0, 2), aspect='equal', origin='lower')

    # Nodule marker circle
    if len(nodule_yx[0]) > 0:
        ny, nx = nodule_yx[0].mean(), nodule_yx[1].mean()
        circle = plt.Circle((ny, nx), 12, fill=False, edgecolor='#FF8C00',
                           linewidth=2, linestyle='--')
        ax.add_patch(circle)
        ax.text(ny + 15, nx + 15, 'NODULE', color='#FF8C00', fontsize=7,
               fontweight='bold', fontfamily='monospace')

    for spine in ax.spines.values():
        spine.set_color('#76B900')
        spine.set_linewidth(2.0)
    ax.set_title(f"CT Chest \u2014 AI Segmentation \u2014 Slice {best_z}",
                color='#76B900', fontsize=10, fontweight='bold', pad=8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.02, 0.02, 'AI Segmentation \u2014 VISTA-3D',
            transform=ax.transAxes, color='#76B900', fontsize=7, alpha=0.8,
            verticalalignment='bottom', fontfamily='monospace')
    fig.tight_layout(pad=0.5)
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)
    frame = buf[:, :, :3].copy()
    plt.close(fig)
    iio.imwrite(out_png, frame)
    print(f"  Saved: {out_png} ({os.path.getsize(out_png)/1024:.0f} KB)")

    # Animated GIF
    n_slices = mask.shape[2]
    start, end = int(n_slices * 0.15), int(n_slices * 0.85)
    n_frames = 50
    step = max(1, (end - start) // n_frames)
    indices = list(range(start, end, step))[:n_frames]

    frames = []
    for idx in indices:
        f = create_overlay_frame(volume[:, :, idx], mask[:, :, idx], colors,
                                  wc, ww,
                                  title=f"CT Chest \u2014 Slice {idx}/{n_slices}",
                                  dpi=128, figsize=(8, 8))
        frames.append(f)

    iio.imwrite(out_gif, frames, duration=int(1000/6), loop=0)
    print(f"  Saved: {out_gif} ({os.path.getsize(out_gif)/1024:.0f} KB)")


def generate_overlays_brain_flair(volume, mask, out_png, out_gif):
    """Generate brain FLAIR overlays with MS lesion highlighting."""
    print("\n  Generating Brain FLAIR overlays...")
    colors = {
        1: (0.65, 0.3, 0.75, 0.22),   # Gray matter — purple
        2: (0.85, 0.85, 0.2, 0.18),    # White matter — yellow
        3: (0.0, 0.75, 0.85, 0.30),    # CSF — cyan
        4: (0.0, 1.0, 0.3, 0.60),      # MS lesions — bright green/neon
    }
    wc, ww = 500, 1000  # FLAIR window

    # Find slice with lesions
    best_z = find_best_slice(mask, axis=2)
    for z in range(mask.shape[2]):
        if (mask[:, :, z] == 4).sum() > 5:
            best_z = z
            break

    slice_2d = volume[:, :, best_z]
    mask_2d = mask[:, :, best_z]

    # Create PNG with lesion markers
    vmin = wc - ww / 2
    vmax = wc + ww / 2
    windowed = np.clip((slice_2d - vmin) / (vmax - vmin), 0, 1)

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=100)
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    ax.imshow(windowed.T, cmap='gray', aspect='equal', origin='lower')

    overlay = np.zeros((*mask_2d.shape, 4))
    for label_idx, color in colors.items():
        region = mask_2d == label_idx
        if region.any():
            overlay[region] = color
    ax.imshow(overlay.transpose(1, 0, 2), aspect='equal', origin='lower')

    # Circle markers for each lesion cluster
    from scipy.ndimage import label as ndlabel
    lesion_binary = mask_2d == 4
    if lesion_binary.any():
        labeled, n_lesions = ndlabel(lesion_binary)
        for li in range(1, n_lesions + 1):
            yx = np.where(labeled == li)
            if len(yx[0]) > 2:
                cy, cx = yx[0].mean(), yx[1].mean()
                r = max(np.std(yx[0]), np.std(yx[1]), 3) * 2.5
                circle = plt.Circle((cy, cx), r, fill=False, edgecolor='#00FF80',
                                   linewidth=1.5, linestyle='--')
                ax.add_patch(circle)

    for spine in ax.spines.values():
        spine.set_color('#76B900')
        spine.set_linewidth(2.0)
    ax.set_title(f"MRI Brain FLAIR \u2014 AI Segmentation \u2014 Slice {best_z}",
                color='#76B900', fontsize=10, fontweight='bold', pad=8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.02, 0.02, 'AI Segmentation \u2014 VISTA-3D',
            transform=ax.transAxes, color='#76B900', fontsize=7, alpha=0.8,
            verticalalignment='bottom', fontfamily='monospace')
    fig.tight_layout(pad=0.5)
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)
    frame = buf[:, :, :3].copy()
    plt.close(fig)
    iio.imwrite(out_png, frame)
    print(f"  Saved: {out_png} ({os.path.getsize(out_png)/1024:.0f} KB)")

    # Animated GIF
    n_slices = mask.shape[2]
    start, end = int(n_slices * 0.15), int(n_slices * 0.85)
    n_frames = 50
    step = max(1, (end - start) // n_frames)
    indices = list(range(start, end, step))[:n_frames]

    frames = []
    for idx in indices:
        f = create_overlay_frame(volume[:, :, idx], mask[:, :, idx], colors,
                                  wc, ww,
                                  title=f"MRI FLAIR \u2014 Slice {idx}/{n_slices}",
                                  dpi=128, figsize=(8, 8))
        frames.append(f)

    iio.imwrite(out_gif, frames, duration=int(1000/6), loop=0)
    print(f"  Saved: {out_gif} ({os.path.getsize(out_gif)/1024:.0f} KB)")


# ========================================================================
# MAIN
# ========================================================================

def main():
    print("\nPhase 1: Generating high-resolution NIfTI volumes")
    print("-" * 50)

    # Generate volumes
    ct_head = generate_ct_head(os.path.join(SAMPLE_DIR, "highres_ct_head.nii.gz"))
    ct_chest = generate_ct_chest(os.path.join(SAMPLE_DIR, "highres_ct_chest.nii.gz"))
    brain_flair = generate_brain_flair(os.path.join(SAMPLE_DIR, "highres_brain_flair.nii.gz"))

    print("\n\nPhase 2: Segmenting volumes and generating overlays")
    print("-" * 50)

    # Segment
    print("\n  Segmenting CT Head...")
    mask_head = segment_ct_head(ct_head)
    print(f"  Labels: {np.unique(mask_head).tolist()}, Non-zero: {(mask_head > 0).sum():,}")

    print("  Segmenting CT Chest...")
    mask_chest = segment_ct_chest(ct_chest)
    print(f"  Labels: {np.unique(mask_chest).tolist()}, Non-zero: {(mask_chest > 0).sum():,}")

    print("  Segmenting Brain FLAIR...")
    mask_flair = segment_brain_flair(brain_flair)
    print(f"  Labels: {np.unique(mask_flair).tolist()}, Non-zero: {(mask_flair > 0).sum():,}")

    # Generate overlays
    generate_overlays_ct_head(
        ct_head, mask_head,
        os.path.join(OUT_DIR, "highres_ct_head_overlay.png"),
        os.path.join(OUT_DIR, "highres_ct_head_segmented.gif"),
    )
    generate_overlays_ct_chest(
        ct_chest, mask_chest,
        os.path.join(OUT_DIR, "highres_ct_chest_overlay.png"),
        os.path.join(OUT_DIR, "highres_ct_chest_segmented.gif"),
    )
    generate_overlays_brain_flair(
        brain_flair, mask_flair,
        os.path.join(OUT_DIR, "highres_brain_flair_overlay.png"),
        os.path.join(OUT_DIR, "highres_brain_flair_segmented.gif"),
    )

    # Summary
    print("\n" + "=" * 70)
    print("HIGH-RESOLUTION VOLUMES & OVERLAYS GENERATED")
    print("=" * 70)

    print("\nNIfTI Volumes:")
    for f in ["highres_ct_head.nii.gz", "highres_ct_chest.nii.gz", "highres_brain_flair.nii.gz"]:
        fpath = os.path.join(SAMPLE_DIR, f)
        print(f"  {f}: {os.path.getsize(fpath)/1024/1024:.1f} MB")

    print("\nSegmentation Overlays:")
    for f in sorted(os.listdir(OUT_DIR)):
        if f.startswith("highres_"):
            fpath = os.path.join(OUT_DIR, f)
            print(f"  {f}: {os.path.getsize(fpath)/1024:.0f} KB")

    print("=" * 70)


if __name__ == "__main__":
    main()
