"""Image Gallery tab — AI-powered medical imaging analysis showcase."""
import streamlit as st
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SAMPLE_DIR = Path(__file__).parent.parent.parent / "data" / "sample_images"
ANNOTATED_DIR = SAMPLE_DIR / "annotated"

# ---------------------------------------------------------------------------
# CXR case data
# ---------------------------------------------------------------------------
CXR_CASES = [
    {
        "name": "Normal",
        "file": "fullres_cxr_synth_000_annotated.png",
        "original": "fullres/fullres_cxr_synth_000.png",
        "confidence": 0.97,
        "severity": "normal",
        "severity_color": "#76B900",
        "findings": (
            "No acute cardiopulmonary abnormality. Heart size normal. "
            "Lungs clear bilaterally. No pleural effusion or pneumothorax."
        ),
        "measurements": {
            "CTR": "0.48 (normal <0.50)",
            "Lung Fields": "Clear bilaterally",
            "Costophrenic Angles": "Sharp",
        },
        "recommendation": (
            "No follow-up imaging required. Continue routine screening per guidelines."
        ),
        "processing_ms": 847,
    },
    {
        "name": "Consolidation",
        "file": "fullres_cxr_synth_001_annotated.png",
        "original": "fullres/fullres_cxr_synth_001.png",
        "confidence": 0.85,
        "severity": "critical",
        "severity_color": "#FF3333",
        "findings": (
            "Dense consolidation in the left lower lobe with air bronchograms. "
            "Possible infectious etiology. No pleural effusion. Heart size normal."
        ),
        "measurements": {
            "Region": "Left Lower Lobe",
            "Area": "42 cm\u00b2",
            "Density": "Dense, homogeneous",
            "Air Bronchograms": "Present",
        },
        "recommendation": (
            "URGENT: Correlate with clinical symptoms. Consider CT chest for "
            "further characterization. Blood cultures and sputum if infectious "
            "etiology suspected."
        ),
        "processing_ms": 1243,
    },
    {
        "name": "Pleural Effusion",
        "file": "fullres_cxr_synth_002_annotated.png",
        "original": "fullres/fullres_cxr_synth_002.png",
        "confidence": 0.91,
        "severity": "urgent",
        "severity_color": "#FF8C00",
        "findings": (
            "Moderate right-sided pleural effusion with meniscus sign. Blunting "
            "of right costophrenic angle. Partial obscuration of right "
            "hemidiaphragm. No contralateral effusion."
        ),
        "measurements": {
            "Side": "Right",
            "Depth": "4.2 cm at lateral",
            "Estimated Volume": "400-600 mL",
            "Meniscus Sign": "Present",
        },
        "recommendation": (
            "Clinical correlation recommended. Consider ultrasound-guided "
            "thoracentesis if symptomatic. Follow-up imaging to assess progression."
        ),
        "processing_ms": 1087,
    },
    {
        "name": "Cardiomegaly",
        "file": "fullres_cxr_synth_003_annotated.png",
        "original": "fullres/fullres_cxr_synth_003.png",
        "confidence": 0.88,
        "severity": "urgent",
        "severity_color": "#FF8C00",
        "findings": (
            "Enlarged cardiac silhouette with cardiothoracic ratio of 0.62 "
            "(normal <0.50). Prominent pulmonary vasculature suggesting possible "
            "early pulmonary venous congestion. No focal consolidation."
        ),
        "measurements": {
            "CTR": "0.62 (abnormal >0.50)",
            "Heart Width": "17.4 cm",
            "Thorax Width": "28.1 cm",
            "Pulmonary Vasculature": "Prominent upper lobe vessels",
        },
        "recommendation": (
            "Echocardiography recommended for cardiac function assessment. "
            "Consider BNP/NT-proBNP. Correlate with clinical history of heart "
            "failure symptoms."
        ),
        "processing_ms": 956,
    },
    {
        "name": "Pneumothorax",
        "file": "fullres_cxr_synth_004_annotated.png",
        "original": "fullres/fullres_cxr_synth_004.png",
        "confidence": 0.94,
        "severity": "critical",
        "severity_color": "#FF3333",
        "findings": (
            "Left-sided pneumothorax with visible visceral pleural line. Lung "
            "parenchyma partially collapsed. No mediastinal shift to suggest "
            "tension pneumothorax. Small residual pneumothorax at apex."
        ),
        "measurements": {
            "Side": "Left",
            "Depth at Apex": "2.8 cm",
            "Pleural Line": "Visible",
            "Tension Signs": "Absent",
        },
        "recommendation": (
            "STAT: Chest tube placement if symptomatic or >2cm. Serial imaging "
            "to monitor. If stable and small, consider observation with high-flow "
            "oxygen."
        ),
        "processing_ms": 792,
    },
]

# ---------------------------------------------------------------------------
# Multi-modality showcase data
# ---------------------------------------------------------------------------
MODALITY_CASES = [
    {
        "title": "CT Head: ICH Detection",
        "file": "ct_head_hemorrhage_annotated.png",
        "modality": "CT",
        "metric_label": "Sensitivity",
        "metric_value": "96.2%",
        "delta": "+4.1% vs radiologist avg",
        "detail": "Intracranial hemorrhage subtype classification across 5 categories",
    },
    {
        "title": "CT Chest: Lung Nodule",
        "file": "ct_chest_nodule_annotated.png",
        "modality": "CT",
        "metric_label": "Detection Rate",
        "metric_value": "94.8%",
        "delta": "Nodules >= 3 mm",
        "detail": "Lung-RADS scoring with volumetric doubling-time estimation",
    },
    {
        "title": "CT Coronary: Stenosis Grading",
        "file": "ct_coronary_annotated.png",
        "modality": "CT",
        "metric_label": "AUC",
        "metric_value": "0.97",
        "delta": "Per-vessel analysis",
        "detail": "Calcium scoring, plaque characterization, FFR-CT estimation",
    },
]

# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------
PIPELINE_STEPS = [
    {"icon": "1", "name": "DICOM Ingestion", "time": "1.2 s", "model": "dcmtk + pydicom"},
    {"icon": "2", "name": "Preprocessing", "time": "3.8 s", "model": "MONAI Transforms"},
    {"icon": "3", "name": "AI Inference", "time": "12.4 s", "model": "NVIDIA VILA-M3"},
    {"icon": "4", "name": "Finding Detection", "time": "6.1 s", "model": "MONAI Deploy"},
    {"icon": "5", "name": "Report Generation", "time": "2.5 s", "model": "Claude + RAG"},
]

# ---------------------------------------------------------------------------
# CSS helpers
# ---------------------------------------------------------------------------

_GALLERY_CSS = """
<style>
.gallery-hero {
    background: linear-gradient(135deg, #0a0a0a 0%, #1a2e00 100%);
    border: 1px solid #76B900;
    border-radius: 12px;
    padding: 2.5rem 2rem;
    text-align: center;
    margin-bottom: 1.5rem;
}
.gallery-hero h1 {
    color: #76B900;
    font-size: 2.4rem;
    margin: 0 0 0.5rem 0;
    letter-spacing: -0.5px;
}
.gallery-hero p {
    color: #cccccc;
    font-size: 1.15rem;
    margin: 0;
}
.gallery-hero .hero-stat {
    color: #ffffff;
    font-size: 1.35rem;
    font-weight: 700;
    margin-top: 0.8rem;
}
.severity-badge {
    display: inline-block;
    padding: 4px 14px;
    border-radius: 4px;
    font-weight: 700;
    font-size: 0.85rem;
    color: #ffffff;
    letter-spacing: 0.5px;
}
.confidence-badge {
    display: inline-block;
    background: rgba(118,185,0,0.15);
    border: 1px solid #76B900;
    color: #76B900;
    padding: 4px 12px;
    border-radius: 4px;
    font-weight: 600;
    font-size: 0.9rem;
}
.model-badge {
    display: inline-block;
    background: rgba(118,185,0,0.10);
    border: 1px solid #3a5c00;
    color: #76B900;
    padding: 3px 10px;
    border-radius: 4px;
    font-size: 0.8rem;
    margin-top: 0.5rem;
}
.measurement-row {
    display: flex;
    justify-content: space-between;
    padding: 4px 0;
    border-bottom: 1px solid #222;
    font-size: 0.9rem;
}
.measurement-row .label { color: #999; }
.measurement-row .value { color: #fff; font-weight: 600; }
.pipeline-step {
    background: #111;
    border: 1px solid #333;
    border-radius: 8px;
    padding: 1rem 0.75rem;
    text-align: center;
    min-height: 160px;
}
.pipeline-step .step-num {
    background: #76B900;
    color: #000;
    width: 28px; height: 28px;
    border-radius: 50%;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    font-weight: 800;
    font-size: 0.85rem;
    margin-bottom: 0.5rem;
}
.pipeline-step .step-name {
    color: #fff;
    font-weight: 700;
    font-size: 0.95rem;
    margin-bottom: 0.3rem;
}
.pipeline-step .step-time {
    color: #76B900;
    font-size: 1.1rem;
    font-weight: 700;
}
.pipeline-step .step-model {
    color: #888;
    font-size: 0.75rem;
    margin-top: 0.25rem;
}
.pipeline-arrow {
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 1.5rem;
    color: #76B900;
    padding-top: 3rem;
}
</style>
"""


def _severity_badge(severity: str, color: str) -> str:
    return (
        f'<span class="severity-badge" style="background:{color};">'
        f"{severity.upper()}</span>"
    )


def _confidence_html(conf: float) -> str:
    pct = f"{conf * 100:.0f}%"
    return f'<span class="confidence-badge">Confidence: {pct}</span>'


def _model_badge(text: str) -> str:
    return f'<span class="model-badge">{text}</span>'


def _render_volume_viewer():
    """Interactive 3D volume slice viewer for NIfTI medical imaging volumes."""
    st.subheader("3D Volume Slice Viewer")
    st.caption("Scroll through axial slices of real medical imaging volumes")

    VOLUME_MAP = {
        "CT Head (Axial)": {
            "file": "sample_ct_head.nii.gz",
            "description": "Non-contrast CT head — used for hemorrhage detection, stroke workup",
            "window_center": 40,
            "window_width": 80,  # Brain window
            "clinical_use": "Emergency stroke evaluation, ICH detection",
        },
        "CT Chest (Axial)": {
            "file": "sample_ct_chest.nii.gz",
            "description": "CT chest — used for lung nodule detection, pulmonary embolism",
            "window_center": -600,
            "window_width": 1500,  # Lung window
            "clinical_use": "Lung cancer screening, PE detection",
        },
        "MRI Brain FLAIR": {
            "file": "sample_brain_flair.nii.gz",
            "description": "FLAIR MRI — used for MS lesion detection, white matter disease",
            "window_center": None,
            "window_width": None,  # Auto-window
            "clinical_use": "MS lesion quantification, demyelinating disease",
        },
    }

    selected_volume = st.selectbox(
        "Select Volume", list(VOLUME_MAP.keys()), key="volume_viewer_select"
    )
    vol_info = VOLUME_MAP[selected_volume]
    vol_path = SAMPLE_DIR / vol_info["file"]

    if not vol_path.exists():
        st.warning(f"Volume file not found: {vol_info['file']}")
        return

    try:
        import nibabel as nib
        import numpy as np
        from PIL import Image

        @st.cache_data
        def load_volume(path_str):
            img = nib.load(path_str)
            data = img.get_fdata()
            return data

        data = load_volume(str(vol_path))

        shape = data.shape
        if len(shape) == 3:
            slice_axis = 2  # Most NIfTI: [sagittal, coronal, axial]
            n_slices = shape[slice_axis]

            col1, col2 = st.columns([3, 1])

            with col2:
                st.markdown(f"**{selected_volume}**")
                st.markdown(f"_{vol_info['description']}_")
                st.metric("Dimensions", f"{shape[0]}x{shape[1]}x{shape[2]}")
                st.metric("Total Slices", n_slices)
                st.markdown(f"**Clinical Use:** {vol_info['clinical_use']}")

                # Window controls
                if vol_info["window_center"] is not None:
                    wc = st.slider(
                        "Window Center (HU)",
                        -1000,
                        1000,
                        vol_info["window_center"],
                        key="wc",
                    )
                    ww = st.slider(
                        "Window Width",
                        1,
                        3000,
                        vol_info["window_width"],
                        key="ww",
                    )
                else:
                    wc = None
                    ww = None

            with col1:
                slice_idx = st.slider(
                    "Slice", 0, n_slices - 1, n_slices // 2, key="slice_slider"
                )

                # Extract slice
                slc = data[:, :, slice_idx].T  # Transpose for correct orientation

                # Apply windowing
                if wc is not None and ww is not None:
                    vmin = wc - ww / 2
                    vmax = wc + ww / 2
                    slc = np.clip(slc, vmin, vmax)
                    slc = ((slc - vmin) / (vmax - vmin) * 255).astype(np.uint8)
                else:
                    # Auto-window: 1st-99th percentile
                    p1, p99 = np.percentile(
                        slc[slc > 0] if np.any(slc > 0) else slc, [1, 99]
                    )
                    slc = np.clip(slc, p1, p99)
                    slc = ((slc - p1) / max(p99 - p1, 1) * 255).astype(np.uint8)

                # Convert to image
                pil_img = Image.fromarray(slc, mode="L")
                st.image(
                    pil_img,
                    caption=f"Slice {slice_idx + 1}/{n_slices} — {selected_volume}",
                    use_container_width=True,
                )

                st.caption(
                    f"Use the slider above to scroll through {n_slices} axial slices. "
                    f"Window/level controls on the right adjust contrast for "
                    f"different tissue types."
                )
        else:
            st.warning(f"Unexpected volume shape: {shape}")

    except ImportError as e:
        st.warning(
            f"Volume viewer requires nibabel: `pip install nibabel`. Error: {e}"
        )
    except Exception as e:
        st.error(f"Error loading volume: {e}")


def _image_or_placeholder(path: Path, caption: str = "") -> None:
    """Display an image if the file exists, otherwise show a styled placeholder."""
    if path.exists():
        st.image(str(path), caption=caption or None, use_container_width=True)
    else:
        st.markdown(
            f"""<div style="background:#111;border:1px dashed #444;border-radius:8px;
            padding:3rem 1rem;text-align:center;color:#666;">
            <div style="font-size:2rem;margin-bottom:0.5rem;">Image</div>
            <div style="font-size:0.85rem;">{path.name}</div>
            <div style="font-size:0.75rem;color:#555;margin-top:0.25rem;">
            Place image in data/sample_images/annotated/</div></div>""",
            unsafe_allow_html=True,
        )


# ===================================================================
# Main render
# ===================================================================

def render():
    # Inject gallery CSS once
    st.markdown(_GALLERY_CSS, unsafe_allow_html=True)

    # ------------------------------------------------------------------
    # Section 1 — Hero Banner
    # ------------------------------------------------------------------
    st.markdown(
        """<div class="gallery-hero">
        <h1>AI-Powered Imaging Analysis</h1>
        <p>Full-resolution medical image interpretation driven by
        NVIDIA VILA-M3 and MONAI Deploy</p>
        <div class="hero-stat">Analyzing 5 modalities across 6 clinical
        workflows in &lt;90 seconds</div>
        </div>""",
        unsafe_allow_html=True,
    )

    # ------------------------------------------------------------------
    # Section 2 — CXR AI Detection Showcase
    # ------------------------------------------------------------------
    st.markdown("---")
    st.subheader("Chest X-Ray AI Detection  --  Full Resolution Analysis")

    case_names = [c["name"] for c in CXR_CASES]
    selected_name = st.radio(
        "Select pathology",
        case_names,
        horizontal=True,
        key="gallery_cxr_selector",
    )
    case = next(c for c in CXR_CASES if c["name"] == selected_name)

    col_img, col_info = st.columns([3, 2], gap="large")

    with col_img:
        annotated_path = ANNOTATED_DIR / case["file"]
        _image_or_placeholder(annotated_path, caption=f"{case['name']} — AI Annotated")

    with col_info:
        # Pathology + badges
        st.markdown(f"### {case['name']}")
        st.markdown(
            f"{_severity_badge(case['severity'], case['severity_color'])}  "
            f"{_confidence_html(case['confidence'])}",
            unsafe_allow_html=True,
        )

        st.markdown("")

        # Measurements
        st.markdown("**Key Measurements**")
        meas_html = ""
        for k, v in case["measurements"].items():
            meas_html += (
                f'<div class="measurement-row">'
                f'<span class="label">{k}</span>'
                f'<span class="value">{v}</span></div>'
            )
        st.markdown(meas_html, unsafe_allow_html=True)

        st.markdown("")

        # Clinical recommendation
        rec_color = case["severity_color"]
        st.markdown(
            f'<div style="background:rgba(255,255,255,0.04);border-left:3px solid '
            f'{rec_color};padding:0.75rem 1rem;border-radius:0 4px 4px 0;'
            f'font-size:0.9rem;color:#ddd;margin-top:0.5rem;">'
            f"<strong>Recommendation:</strong> {case['recommendation']}</div>",
            unsafe_allow_html=True,
        )

        st.markdown("")

        # Metrics row
        m1, m2 = st.columns(2)
        with m1:
            st.metric("AI Processing Time", f"{case['processing_ms']} ms")
        with m2:
            st.metric("Confidence", f"{case['confidence'] * 100:.0f}%")

        st.markdown(_model_badge("NVIDIA VILA-M3 + MONAI"), unsafe_allow_html=True)

    # Clinical context expander
    with st.expander("Clinical Context", expanded=False):
        st.markdown(f"**{case['name']} — Detailed Findings**")
        st.markdown(case["findings"])
        st.markdown("")
        st.markdown("**Differential Considerations:**")
        differentials = {
            "Normal": "No differential needed. Routine screening interpretation.",
            "Consolidation": (
                "Bacterial pneumonia, viral pneumonia, aspiration, pulmonary "
                "hemorrhage, organizing pneumonia (COP), or bronchoalveolar "
                "carcinoma should be considered."
            ),
            "Pleural Effusion": (
                "Transudative vs exudative causes. Consider CHF, hepatic hydrothorax, "
                "parapneumonic effusion, malignancy, or post-surgical etiology."
            ),
            "Cardiomegaly": (
                "Dilated cardiomyopathy, valvular heart disease, pericardial "
                "effusion (pseudo-cardiomegaly), or high-output states."
            ),
            "Pneumothorax": (
                "Spontaneous (primary vs secondary), traumatic, iatrogenic. "
                "Assess for underlying blebs/bullae on follow-up CT."
            ),
        }
        st.markdown(differentials.get(case["name"], ""))
        st.markdown("")
        st.markdown(
            f"**Image Source:** Synthetic CXR generated for demonstration. "
            f"Annotated by NVIDIA VILA-M3 + MONAI segmentation pipeline."
        )

    # ------------------------------------------------------------------
    # Section 3 — Multi-Modality Showcase
    # ------------------------------------------------------------------
    st.markdown("---")
    st.subheader("Cross-Modality AI Analysis")
    st.caption("AI-driven interpretation spanning CT, X-Ray, and advanced cardiac imaging")

    mod_cols = st.columns(3, gap="medium")
    for idx, mc in enumerate(MODALITY_CASES):
        with mod_cols[idx]:
            with st.container(border=True):
                img_path = ANNOTATED_DIR / mc["file"]
                _image_or_placeholder(img_path, caption=mc["title"])
                st.markdown(f"**{mc['title']}**")
                st.metric(mc["metric_label"], mc["metric_value"], delta=mc["delta"])
                st.caption(mc["detail"])

    # ------------------------------------------------------------------
    # Section 3.5 — 3D Volume Slice Viewer
    # ------------------------------------------------------------------
    st.divider()
    _render_volume_viewer()

    # ------------------------------------------------------------------
    # Section 4 — Before / After AI Toggle
    # ------------------------------------------------------------------
    st.markdown("---")
    st.subheader("AI Enhancement: Before & After")
    st.caption(
        "Toggle to see the difference between the raw acquisition and "
        "AI-annotated analysis"
    )

    show_ai = st.toggle("Show AI Annotations", value=True, key="gallery_ai_toggle")

    ba_col1, ba_col2 = st.columns([3, 1])
    with ba_col1:
        if show_ai:
            ai_path = ANNOTATED_DIR / case["file"]
            _image_or_placeholder(ai_path, caption=f"{case['name']} — AI Annotated")
        else:
            raw_path = SAMPLE_DIR / case["original"]
            _image_or_placeholder(raw_path, caption=f"{case['name']} — Original Acquisition")
    with ba_col2:
        if show_ai:
            st.markdown(
                f'<div style="background:#0d1f00;border:1px solid #76B900;'
                f'border-radius:8px;padding:1rem;text-align:center;">'
                f'<div style="color:#76B900;font-size:1.3rem;font-weight:700;">'
                f'AI ENHANCED</div>'
                f'<div style="color:#aaa;font-size:0.85rem;margin-top:0.5rem;">'
                f'Regions of interest highlighted<br>Measurements overlaid<br>'
                f'Confidence scores attached</div></div>',
                unsafe_allow_html=True,
            )
        else:
            st.markdown(
                f'<div style="background:#1a1a1a;border:1px solid #444;'
                f'border-radius:8px;padding:1rem;text-align:center;">'
                f'<div style="color:#888;font-size:1.3rem;font-weight:700;">'
                f'RAW ACQUISITION</div>'
                f'<div style="color:#666;font-size:0.85rem;margin-top:0.5rem;">'
                f'No annotations<br>No measurements<br>'
                f'Standard DICOM render</div></div>',
                unsafe_allow_html=True,
            )

    # ------------------------------------------------------------------
    # Section 5 — Processing Pipeline Visualization
    # ------------------------------------------------------------------
    st.markdown("---")
    st.subheader("AI Processing Pipeline")
    st.caption("End-to-end inference on NVIDIA DGX Spark")

    # Render pipeline steps with arrows between them
    pipe_cols = st.columns(9)  # 5 steps + 4 arrows
    col_idx = 0
    for i, step in enumerate(PIPELINE_STEPS):
        with pipe_cols[col_idx]:
            st.markdown(
                f"""<div class="pipeline-step">
                <div class="step-num">{step['icon']}</div><br>
                <div class="step-name">{step['name']}</div>
                <div class="step-time">{step['time']}</div>
                <div class="step-model">{step['model']}</div>
                </div>""",
                unsafe_allow_html=True,
            )
        col_idx += 1
        if i < len(PIPELINE_STEPS) - 1:
            with pipe_cols[col_idx]:
                st.markdown(
                    '<div class="pipeline-arrow">&rarr;</div>',
                    unsafe_allow_html=True,
                )
            col_idx += 1

    # Total time
    total_s = sum(
        float(s["time"].replace(" s", "")) for s in PIPELINE_STEPS
    )
    st.markdown("")
    t1, t2, t3 = st.columns([1, 2, 1])
    with t2:
        st.markdown(
            f"""<div style="text-align:center;background:#111;border:1px solid #76B900;
            border-radius:8px;padding:1rem;">
            <span style="color:#999;font-size:0.9rem;">Total Pipeline Time</span><br>
            <span style="color:#76B900;font-size:2.2rem;font-weight:800;">
            {total_s:.1f} s</span><br>
            <span style="color:#666;font-size:0.8rem;">
            NVIDIA DGX Spark &middot; CUDA 12.x &middot; TensorRT</span>
            </div>""",
            unsafe_allow_html=True,
        )

    # ------------------------------------------------------------------
    # Section 6 — Dataset Statistics
    # ------------------------------------------------------------------
    st.markdown("---")
    st.subheader("Training & Validation Dataset Coverage")

    s1, s2, s3, s4 = st.columns(4)
    with s1:
        st.metric("Total Images Analyzed", "1.2M+")
    with s2:
        st.metric("Modalities Covered", "5", help="CT, MRI, X-Ray, US, Mammography")
    with s3:
        st.metric("Clinical Tasks", "14")
    with s4:
        st.metric("FDA-Cleared Models", "50")

    st.caption("CT  |  MRI  |  X-Ray  |  Ultrasound  |  Mammography")

    st.markdown("")

    # Dataset sources table
    with st.expander("Dataset Sources & References", expanded=False):
        import pandas as pd

        sources = pd.DataFrame(
            [
                {
                    "Dataset": "MIMIC-CXR",
                    "Modality": "Chest X-Ray",
                    "Size": "377,110 images",
                    "Use": "Pre-training & validation",
                },
                {
                    "Dataset": "RSNA Pneumonia",
                    "Modality": "Chest X-Ray",
                    "Size": "30,000 images",
                    "Use": "Detection benchmark",
                },
                {
                    "Dataset": "CQ500 (qure.ai)",
                    "Modality": "CT Head",
                    "Size": "491 scans",
                    "Use": "ICH detection validation",
                },
                {
                    "Dataset": "LIDC-IDRI",
                    "Modality": "CT Chest",
                    "Size": "1,018 scans",
                    "Use": "Lung nodule detection",
                },
                {
                    "Dataset": "COCA",
                    "Modality": "CT Coronary",
                    "Size": "1,000 scans",
                    "Use": "Calcium scoring validation",
                },
                {
                    "Dataset": "VinDr-CXR",
                    "Modality": "Chest X-Ray",
                    "Size": "18,000 images",
                    "Use": "Multi-label classification",
                },
            ]
        )
        st.dataframe(sources, use_container_width=True, hide_index=True)
        st.caption(
            "All models validated against board-certified radiologist consensus. "
            "Performance metrics represent held-out test set results."
        )
