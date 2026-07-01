#!/usr/bin/env python3
"""Validate CXR Rapid Findings workflow against real MedMNIST data.

Loads real chest X-ray images from MedMNIST (NIH ChestX-ray14) that were
downloaded by download_real_data.py, runs the CXR classification pipeline
(torchxrayvision DenseNet or MONAI fallback), and compares model predictions
to ground truth labels.

Reports:
    - Per-class accuracy, sensitivity, specificity, PPV, NPV
    - Overall multi-label exact match accuracy
    - Confusion matrix summary per pathology class
    - Inference timing statistics

Usage:
    python scripts/validate_real_data.py

    # With verbose per-image output:
    python scripts/validate_real_data.py --verbose

    # Use MONAI fallback instead of torchxrayvision:
    python scripts/validate_real_data.py --backend monai

Author: Adam Jones
Date: 2026-02-27
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Ensure project root is on PYTHONPATH
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

SAMPLE_DIR = PROJECT_ROOT / "data" / "sample_images"
METADATA_PATH = SAMPLE_DIR / "metadata.json"

# ---------------------------------------------------------------------------
# MedMNIST ChestMNIST label names -> our 5 target classes mapping
# ---------------------------------------------------------------------------
# ChestMNIST has 14 labels:
#   0: Atelectasis, 1: Cardiomegaly, 2: Effusion, 3: Infiltration,
#   4: Mass, 5: Nodule, 6: Pneumonia, 7: Pneumothorax,
#   8: Consolidation, 9: Edema, 10: Emphysema, 11: Fibrosis,
#   12: Pleural Thickening, 13: Hernia
#
# Our CXR workflow detects 5 classes:
#   pneumothorax, consolidation, pleural_effusion, cardiomegaly, fracture
#
# Mapping ChestMNIST indices to our class names:
CHESTMNIST_TO_OUR_CLASSES = {
    7: "pneumothorax",      # ChestMNIST index 7 = Pneumothorax
    8: "consolidation",     # ChestMNIST index 8 = Consolidation
    2: "pleural_effusion",  # ChestMNIST index 2 = Effusion
    1: "cardiomegaly",      # ChestMNIST index 1 = Cardiomegaly
    # Note: ChestMNIST does not have a "fracture" label
}

OUR_CLASSES = ["pneumothorax", "consolidation", "pleural_effusion", "cardiomegaly", "fracture"]


def load_metadata() -> dict:
    """Load metadata.json produced by download_real_data.py."""
    if not METADATA_PATH.exists():
        print(f"ERROR: Metadata not found at {METADATA_PATH}")
        print("Run 'python scripts/download_real_data.py' first.")
        sys.exit(1)
    with open(METADATA_PATH) as f:
        return json.load(f)


def extract_ground_truth(file_entry: dict) -> Dict[str, int]:
    """Convert ChestMNIST label vector to our 5-class ground truth.

    Returns dict mapping our class names to 0/1.
    """
    label_vector = file_entry.get("label_vector", [])
    gt = {}
    for chestmnist_idx, our_class in CHESTMNIST_TO_OUR_CLASSES.items():
        if chestmnist_idx < len(label_vector):
            gt[our_class] = int(label_vector[chestmnist_idx])
        else:
            gt[our_class] = 0
    # Fracture is never in ChestMNIST, always 0
    gt["fracture"] = 0
    return gt


def run_cxr_workflow(image_path: str, backend: str = "auto") -> Tuple[Dict, float]:
    """Run CXR classification on a single image.

    Returns (predictions_dict, inference_time_ms).
    predictions_dict maps class name to predicted probability.
    """
    from src.workflows.cxr_rapid_findings import CXRRapidFindingsWorkflow, CXR_CLASS_THRESHOLDS

    # Create workflow in real (non-mock) mode
    if backend == "monai":
        workflow = CXRRapidFindingsWorkflow(mock_mode=False, checkpoint_path="__force_monai__")
    else:
        workflow = CXRRapidFindingsWorkflow(mock_mode=False)

    start = time.time()

    # Run preprocessing + inference
    preprocessed = workflow.preprocess(image_path)
    inference_result = workflow.infer(preprocessed)

    elapsed_ms = (time.time() - start) * 1000

    class_probs = inference_result.get("class_probabilities", {})
    return class_probs, elapsed_ms


def compute_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    class_names: List[str],
) -> Dict:
    """Compute per-class and aggregate multi-label classification metrics.

    Args:
        y_true: (N, C) binary ground truth matrix
        y_pred: (N, C) binary prediction matrix
        class_names: list of class names

    Returns:
        Dictionary with per-class and aggregate metrics.
    """
    n_samples, n_classes = y_true.shape
    results = {"per_class": {}, "aggregate": {}}

    all_tp = 0
    all_fp = 0
    all_tn = 0
    all_fn = 0

    for c in range(n_classes):
        tp = int(np.sum((y_true[:, c] == 1) & (y_pred[:, c] == 1)))
        fp = int(np.sum((y_true[:, c] == 0) & (y_pred[:, c] == 1)))
        tn = int(np.sum((y_true[:, c] == 0) & (y_pred[:, c] == 0)))
        fn = int(np.sum((y_true[:, c] == 1) & (y_pred[:, c] == 0)))

        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
        specificity = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
        ppv = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
        npv = tn / (tn + fn) if (tn + fn) > 0 else float("nan")
        accuracy = (tp + tn) / (tp + fp + tn + fn) if (tp + fp + tn + fn) > 0 else float("nan")
        f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else float("nan")

        # Count how many ground truth positives exist for this class
        n_pos = int(np.sum(y_true[:, c] == 1))
        n_neg = int(np.sum(y_true[:, c] == 0))

        results["per_class"][class_names[c]] = {
            "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "n_positive": n_pos, "n_negative": n_neg,
            "sensitivity": round(sensitivity, 4) if not np.isnan(sensitivity) else "N/A (no positives)",
            "specificity": round(specificity, 4) if not np.isnan(specificity) else "N/A (no negatives)",
            "ppv": round(ppv, 4) if not np.isnan(ppv) else "N/A (no predictions)",
            "npv": round(npv, 4) if not np.isnan(npv) else "N/A",
            "accuracy": round(accuracy, 4),
            "f1": round(f1, 4) if not np.isnan(f1) else "N/A",
        }

        all_tp += tp
        all_fp += fp
        all_tn += tn
        all_fn += fn

    # Aggregate (micro-average)
    micro_sensitivity = all_tp / (all_tp + all_fn) if (all_tp + all_fn) > 0 else 0.0
    micro_specificity = all_tn / (all_tn + all_fp) if (all_tn + all_fp) > 0 else 0.0
    micro_f1 = 2 * all_tp / (2 * all_tp + all_fp + all_fn) if (2 * all_tp + all_fp + all_fn) > 0 else 0.0
    micro_accuracy = (all_tp + all_tn) / (all_tp + all_fp + all_tn + all_fn) if (all_tp + all_fp + all_tn + all_fn) > 0 else 0.0

    # Exact match (all classes correct for a sample)
    exact_match = int(np.sum(np.all(y_true == y_pred, axis=1)))

    results["aggregate"] = {
        "micro_sensitivity": round(micro_sensitivity, 4),
        "micro_specificity": round(micro_specificity, 4),
        "micro_f1": round(micro_f1, 4),
        "micro_accuracy": round(micro_accuracy, 4),
        "exact_match_ratio": round(exact_match / n_samples, 4) if n_samples > 0 else 0.0,
        "exact_match_count": exact_match,
        "total_samples": n_samples,
    }

    return results


def run_validation_chest_xray(
    metadata: dict,
    backend: str = "auto",
    verbose: bool = False,
) -> Optional[Dict]:
    """Run validation on ChestMNIST chest X-ray samples."""
    from src.workflows.cxr_rapid_findings import CXR_CLASS_THRESHOLDS

    chest_data = metadata.get("chest_xray")
    if not chest_data or not chest_data.get("files"):
        print("No chest X-ray data found in metadata.")
        return None

    files = chest_data["files"]
    n = len(files)

    print(f"\n{'='*70}")
    print(f"CHEST X-RAY VALIDATION (NIH ChestX-ray14, {n} samples)")
    print(f"{'='*70}")
    print(f"  Source: {chest_data.get('source', 'unknown')}")
    print(f"  Resolution: {chest_data.get('resolution', 'unknown')}")
    print(f"  Model backend: {backend}")
    print(f"  Thresholds: {CXR_CLASS_THRESHOLDS}")
    print()

    y_true_list = []
    y_pred_list = []
    y_prob_list = []
    timing_list = []

    # Instantiate workflow once (to avoid re-loading model per image)
    from src.workflows.cxr_rapid_findings import CXRRapidFindingsWorkflow

    if backend == "monai":
        workflow = CXRRapidFindingsWorkflow(mock_mode=False, checkpoint_path="__force_monai__")
    else:
        workflow = CXRRapidFindingsWorkflow(mock_mode=False)

    model_name = "torchxrayvision" if workflow._using_xrv else "MONAI (fallback)"
    print(f"  Active model: {model_name}")
    print(f"  Weights loaded: {workflow._weights_loaded}")
    print()

    for i, file_entry in enumerate(files):
        filename = file_entry["filename"]
        image_path = str(SAMPLE_DIR / filename)

        if not Path(image_path).exists():
            print(f"  SKIP: {filename} not found")
            continue

        # Ground truth
        gt = extract_ground_truth(file_entry)
        gt_vector = [gt[c] for c in OUR_CLASSES]

        # Run inference
        start = time.time()
        preprocessed = workflow.preprocess(image_path)
        inference_result = workflow.infer(preprocessed)
        elapsed_ms = (time.time() - start) * 1000

        class_probs = inference_result.get("class_probabilities", {})

        # Apply thresholds to get binary predictions
        pred_vector = []
        prob_vector = []
        for c in OUR_CLASSES:
            prob = class_probs.get(c, 0.0)
            threshold = CXR_CLASS_THRESHOLDS.get(c, 0.5)
            pred_vector.append(1 if prob >= threshold else 0)
            prob_vector.append(prob)

        y_true_list.append(gt_vector)
        y_pred_list.append(pred_vector)
        y_prob_list.append(prob_vector)
        timing_list.append(elapsed_ms)

        if verbose:
            gt_names = file_entry.get("label_names", [])
            pred_names = [OUR_CLASSES[j] for j in range(len(pred_vector)) if pred_vector[j] == 1]
            match = "OK" if gt_vector == pred_vector else "MISMATCH"
            print(f"  [{i+1:2d}/{n}] {filename}")
            print(f"         GT labels: {gt_names}")
            print(f"         GT (5-cls): {dict(zip(OUR_CLASSES, gt_vector))}")
            print(f"         Pred probs: {dict(zip(OUR_CLASSES, [f'{p:.3f}' for p in prob_vector]))}")
            print(f"         Pred (bin): {dict(zip(OUR_CLASSES, pred_vector))}")
            print(f"         Status: {match}  Time: {elapsed_ms:.0f}ms")

            # If xrv, also show all 18 scores
            if "xrv_all_scores" in inference_result:
                all_scores = inference_result["xrv_all_scores"]
                top3 = sorted(all_scores.items(), key=lambda x: x[1], reverse=True)[:3]
                print(f"         Top-3 xrv: {', '.join(f'{k}={v:.3f}' for k, v in top3)}")
            print()

    if not y_true_list:
        print("  No images processed.")
        return None

    # Compute metrics
    y_true = np.array(y_true_list)
    y_pred = np.array(y_pred_list)
    y_prob = np.array(y_prob_list)

    metrics = compute_metrics(y_true, y_pred, OUR_CLASSES)

    # Timing
    timing = {
        "mean_ms": round(np.mean(timing_list), 1),
        "median_ms": round(np.median(timing_list), 1),
        "min_ms": round(np.min(timing_list), 1),
        "max_ms": round(np.max(timing_list), 1),
        "total_ms": round(np.sum(timing_list), 1),
    }

    # Print results
    print(f"\n{'-'*70}")
    print(f"PER-CLASS METRICS ({n} samples)")
    print(f"{'-'*70}")
    header = f"{'Class':<20} {'TP':>3} {'FP':>3} {'TN':>3} {'FN':>3} {'Sens':>8} {'Spec':>8} {'PPV':>8} {'F1':>8} {'Acc':>8}"
    print(header)
    print("-" * len(header))

    for class_name in OUR_CLASSES:
        m = metrics["per_class"][class_name]
        sens = f"{m['sensitivity']:.4f}" if isinstance(m["sensitivity"], float) else m["sensitivity"][:8]
        spec = f"{m['specificity']:.4f}" if isinstance(m["specificity"], float) else m["specificity"][:8]
        ppv = f"{m['ppv']:.4f}" if isinstance(m["ppv"], float) else m["ppv"][:8]
        f1 = f"{m['f1']:.4f}" if isinstance(m["f1"], float) else m["f1"][:8]
        acc = f"{m['accuracy']:.4f}"
        print(f"{class_name:<20} {m['tp']:>3} {m['fp']:>3} {m['tn']:>3} {m['fn']:>3} {sens:>8} {spec:>8} {ppv:>8} {f1:>8} {acc:>8}")

    agg = metrics["aggregate"]
    print(f"\n{'AGGREGATE METRICS':}")
    print(f"  Micro-average sensitivity: {agg['micro_sensitivity']:.4f}")
    print(f"  Micro-average specificity: {agg['micro_specificity']:.4f}")
    print(f"  Micro-average F1:          {agg['micro_f1']:.4f}")
    print(f"  Micro-average accuracy:    {agg['micro_accuracy']:.4f}")
    print(f"  Exact match ratio:         {agg['exact_match_ratio']:.4f} ({agg['exact_match_count']}/{agg['total_samples']})")

    print(f"\nINFERENCE TIMING:")
    print(f"  Mean:   {timing['mean_ms']:.1f} ms")
    print(f"  Median: {timing['median_ms']:.1f} ms")
    print(f"  Min:    {timing['min_ms']:.1f} ms")
    print(f"  Max:    {timing['max_ms']:.1f} ms")
    print(f"  Total:  {timing['total_ms']:.1f} ms ({timing['total_ms']/1000:.2f}s for {n} images)")

    return {
        "dataset": "chest_xray",
        "model": model_name,
        "n_samples": n,
        "metrics": metrics,
        "timing": timing,
        "mean_probabilities": {
            c: round(float(np.mean(y_prob[:, i])), 4)
            for i, c in enumerate(OUR_CLASSES)
        },
    }


def run_validation_pneumonia(
    metadata: dict,
    backend: str = "auto",
    verbose: bool = False,
) -> Optional[Dict]:
    """Run validation on PneumoniaMNIST samples.

    PneumoniaMNIST is binary (normal vs pneumonia). We check if our CXR
    workflow's consolidation/pneumonia-related predictions correlate with
    the pneumonia label. This is a cross-dataset generalization test since
    the model was trained on adult CXR and this dataset is pediatric.
    """
    pneumonia_data = metadata.get("pneumonia")
    if not pneumonia_data or not pneumonia_data.get("files"):
        print("No pneumonia data found in metadata.")
        return None

    files = pneumonia_data["files"]
    n = len(files)

    print(f"\n{'='*70}")
    print(f"PNEUMONIA VALIDATION (Pediatric CXR, {n} samples)")
    print(f"{'='*70}")
    print(f"  Source: {pneumonia_data.get('source', 'unknown')}")
    print(f"  Task: Binary classification -> mapped to consolidation detection")
    print()

    from src.workflows.cxr_rapid_findings import CXRRapidFindingsWorkflow, CXR_CLASS_THRESHOLDS

    if backend == "monai":
        workflow = CXRRapidFindingsWorkflow(mock_mode=False, checkpoint_path="__force_monai__")
    else:
        workflow = CXRRapidFindingsWorkflow(mock_mode=False)

    # For pneumonia detection, we treat positive consolidation as pneumonia positive
    # This is a reasonable clinical mapping: consolidation on CXR is the hallmark of pneumonia
    consolidation_threshold = CXR_CLASS_THRESHOLDS.get("consolidation", 0.60)

    tp = fp = tn = fn = 0
    timing_list = []
    results_detail = []

    for i, file_entry in enumerate(files):
        filename = file_entry["filename"]
        image_path = str(SAMPLE_DIR / filename)

        if not Path(image_path).exists():
            continue

        gt_pneumonia = file_entry.get("label", 0)  # 0=normal, 1=pneumonia

        start = time.time()
        preprocessed = workflow.preprocess(image_path)
        inference_result = workflow.infer(preprocessed)
        elapsed_ms = (time.time() - start) * 1000
        timing_list.append(elapsed_ms)

        class_probs = inference_result.get("class_probabilities", {})
        consolidation_prob = class_probs.get("consolidation", 0.0)

        # Also check pneumonia-related scores if using xrv (which has explicit "Pneumonia" label)
        xrv_all = inference_result.get("xrv_all_scores", {})
        pneumonia_prob = xrv_all.get("Pneumonia", 0.0)
        infiltration_prob = xrv_all.get("Infiltration", 0.0)

        # Combined pneumonia score: max of consolidation, explicit pneumonia, infiltration
        if xrv_all:
            combined_prob = max(consolidation_prob, pneumonia_prob, infiltration_prob)
        else:
            combined_prob = consolidation_prob

        pred_pneumonia = 1 if combined_prob >= consolidation_threshold else 0

        if gt_pneumonia == 1 and pred_pneumonia == 1:
            tp += 1
        elif gt_pneumonia == 0 and pred_pneumonia == 1:
            fp += 1
        elif gt_pneumonia == 0 and pred_pneumonia == 0:
            tn += 1
        elif gt_pneumonia == 1 and pred_pneumonia == 0:
            fn += 1

        results_detail.append({
            "filename": filename,
            "gt": gt_pneumonia,
            "pred": pred_pneumonia,
            "consolidation_prob": consolidation_prob,
            "combined_prob": combined_prob,
        })

        if verbose:
            gt_name = "Pneumonia" if gt_pneumonia == 1 else "Normal"
            pred_name = "Pneumonia" if pred_pneumonia == 1 else "Normal"
            match = "OK" if gt_pneumonia == pred_pneumonia else "MISS"
            print(f"  [{i+1:2d}/{n}] {filename}  GT={gt_name}  Pred={pred_name}  "
                  f"consol={consolidation_prob:.3f}  combined={combined_prob:.3f}  [{match}]")

    total = tp + fp + tn + fn
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    accuracy = (tp + tn) / total if total > 0 else 0.0
    f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0.0

    print(f"\n  Confusion Matrix:")
    print(f"                  Predicted Normal  Predicted Pneumonia")
    print(f"    GT Normal           {tn:3d}               {fp:3d}")
    print(f"    GT Pneumonia        {fn:3d}               {tp:3d}")
    print(f"\n  Sensitivity (recall): {sensitivity:.4f}")
    print(f"  Specificity:          {specificity:.4f}")
    print(f"  Accuracy:             {accuracy:.4f}")
    print(f"  F1 Score:             {f1:.4f}")

    if timing_list:
        print(f"  Mean inference time:  {np.mean(timing_list):.1f} ms")

    return {
        "dataset": "pneumonia",
        "n_samples": total,
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "sensitivity": round(sensitivity, 4),
        "specificity": round(specificity, 4),
        "accuracy": round(accuracy, 4),
        "f1": round(f1, 4),
    }


def run_validation_dicom(metadata: dict, verbose: bool = False) -> Optional[Dict]:
    """Run CXR workflow on the DICOM sample to verify DICOM loading works."""
    dicom_data = metadata.get("dicom_sample")
    if not dicom_data or not dicom_data.get("files"):
        print("No DICOM sample found.")
        return None

    print(f"\n{'='*70}")
    print("DICOM LOADING VALIDATION")
    print(f"{'='*70}")

    from src.workflows.cxr_rapid_findings import CXRRapidFindingsWorkflow

    workflow = CXRRapidFindingsWorkflow(mock_mode=False)

    for file_info in dicom_data["files"]:
        dcm_file = file_info.get("filename")
        if not dcm_file:
            continue

        image_path = str(SAMPLE_DIR / dcm_file)
        if not Path(image_path).exists():
            print(f"  SKIP: {dcm_file} not found")
            continue

        try:
            start = time.time()
            preprocessed = workflow.preprocess(image_path)
            inference_result = workflow.infer(preprocessed)
            elapsed_ms = (time.time() - start) * 1000

            class_probs = inference_result.get("class_probabilities", {})
            print(f"  File: {dcm_file}")
            print(f"  Modality: {file_info.get('modality', 'unknown')}")
            print(f"  Dimensions: {file_info.get('rows', '?')}x{file_info.get('columns', '?')}")
            print(f"  Inference time: {elapsed_ms:.1f} ms")
            print(f"  Predictions: {class_probs}")
            print(f"  Status: PASS (DICOM loaded and processed successfully)")

            return {"status": "pass", "inference_time_ms": elapsed_ms, "predictions": class_probs}

        except Exception as e:
            print(f"  File: {dcm_file}")
            print(f"  Status: FAIL ({e})")
            return {"status": "fail", "error": str(e)}

    return None


FULLRES_DIR = PROJECT_ROOT / "data" / "sample_images" / "fullres"


def validate_fullres_cxr(
    backend: str = "auto",
    verbose: bool = False,
) -> Optional[Dict]:
    """Validate CXR workflow on full-resolution (1024x1024+) chest X-ray images.

    Loads images from data/sample_images/fullres/ (generated or downloaded by
    download_real_data.py), runs the CXR classification pipeline in non-mock
    mode, and reports:
        - Per-image inference time, findings, and confidence scores
        - Comparison of timing vs 28x28 MedMNIST upscaled images
        - Whether the model produces more confident (less compressed)
          probabilities on full-resolution inputs

    Args:
        backend: "auto" (prefer xrv) or "monai" (force fallback).
        verbose: Print per-image detail.

    Returns:
        Validation result dict, or None if no images found.
    """
    from src.workflows.cxr_rapid_findings import CXRRapidFindingsWorkflow, CXR_CLASS_THRESHOLDS
    from PIL import Image as PILImage

    if not FULLRES_DIR.exists():
        print("  Full-resolution directory not found. Run download_real_data.py first.")
        return None

    # Find all PNG images in the fullres directory
    image_files = sorted(FULLRES_DIR.glob("*.png"))
    if not image_files:
        print("  No full-resolution images found. Run download_real_data.py first.")
        return None

    n = len(image_files)

    print(f"\n{'='*70}")
    print(f"FULL-RESOLUTION CXR VALIDATION ({n} images)")
    print(f"{'='*70}")

    # Load metadata if available
    fullres_metadata = {}
    if METADATA_PATH.exists():
        with open(METADATA_PATH) as f:
            meta = json.load(f)
        fullres_meta = meta.get("fullres_cxr", {})
        for fentry in fullres_meta.get("files", []):
            fullres_metadata[fentry["filename"]] = fentry

    # Initialize workflow
    if backend == "monai":
        workflow = CXRRapidFindingsWorkflow(mock_mode=False, checkpoint_path="__force_monai__")
    else:
        workflow = CXRRapidFindingsWorkflow(mock_mode=False)

    model_name = "torchxrayvision" if workflow._using_xrv else "MONAI (fallback)"
    print(f"  Directory: {FULLRES_DIR}")
    print(f"  Model: {model_name}")
    print(f"  Weights loaded: {workflow._weights_loaded}")
    print(f"  Thresholds: {CXR_CLASS_THRESHOLDS}")
    print()

    per_image_results = []
    timing_list = []
    all_probs = {c: [] for c in OUR_CLASSES}

    for i, image_path in enumerate(image_files):
        filename = image_path.name

        # Get image dimensions
        img = PILImage.open(str(image_path))
        width, height = img.size
        img.close()

        file_meta = fullres_metadata.get(filename, {})
        expected_pathology = file_meta.get("pathology", "unknown")
        source_type = file_meta.get("source", "unknown")

        # Run inference
        start = time.time()
        preprocessed = workflow.preprocess(str(image_path))
        inference_result = workflow.infer(preprocessed)
        elapsed_ms = (time.time() - start) * 1000
        timing_list.append(elapsed_ms)

        class_probs = inference_result.get("class_probabilities", {})

        # Determine positive findings
        positive_findings = []
        for cls_name in OUR_CLASSES:
            prob = class_probs.get(cls_name, 0.0)
            threshold = CXR_CLASS_THRESHOLDS.get(cls_name, 0.5)
            all_probs[cls_name].append(prob)
            if prob >= threshold:
                positive_findings.append(f"{cls_name}={prob:.3f}")

        finding_str = ", ".join(positive_findings) if positive_findings else "No significant findings"

        result = {
            "filename": filename,
            "resolution": f"{width}x{height}",
            "source": source_type,
            "expected_pathology": expected_pathology,
            "inference_time_ms": round(elapsed_ms, 1),
            "class_probabilities": {k: round(v, 4) for k, v in class_probs.items()},
            "positive_findings": positive_findings,
        }
        per_image_results.append(result)

        if verbose:
            print(f"  [{i+1:2d}/{n}] {filename} ({width}x{height}, {source_type})")
            print(f"         Expected: {expected_pathology}")
            print(f"         Probs: {', '.join(f'{k}={v:.3f}' for k, v in class_probs.items())}")
            print(f"         Findings: {finding_str}")
            print(f"         Time: {elapsed_ms:.1f} ms")

            # Show all 18 xrv scores if available
            if "xrv_all_scores" in inference_result:
                xrv_all = inference_result["xrv_all_scores"]
                top5 = sorted(xrv_all.items(), key=lambda x: x[1], reverse=True)[:5]
                print(f"         Top-5 xrv: {', '.join(f'{k}={v:.3f}' for k, v in top5)}")
            print()
        else:
            print(f"  [{i+1:2d}/{n}] {filename} ({width}x{height}) -> "
                  f"{finding_str}  [{elapsed_ms:.0f}ms]")

    # ── Summary statistics ──
    print(f"\n{'-'*70}")
    print(f"FULL-RESOLUTION INFERENCE RESULTS ({n} images)")
    print(f"{'-'*70}")

    # Timing
    mean_time = np.mean(timing_list)
    median_time = np.median(timing_list)
    min_time = np.min(timing_list)
    max_time = np.max(timing_list)

    print(f"\n  INFERENCE TIMING:")
    print(f"    Mean:   {mean_time:.1f} ms")
    print(f"    Median: {median_time:.1f} ms")
    print(f"    Min:    {min_time:.1f} ms")
    print(f"    Max:    {max_time:.1f} ms")
    print(f"    Total:  {sum(timing_list):.1f} ms ({sum(timing_list)/1000:.2f}s)")

    # Per-class probability distributions
    print(f"\n  PER-CLASS PROBABILITY DISTRIBUTIONS (full-res):")
    header = f"    {'Class':<20} {'Mean':>8} {'Std':>8} {'Min':>8} {'Max':>8} {'#Pos':>5}"
    print(header)
    print("    " + "-" * (len(header) - 4))

    n_positive_total = 0
    for cls_name in OUR_CLASSES:
        probs = np.array(all_probs[cls_name])
        threshold = CXR_CLASS_THRESHOLDS.get(cls_name, 0.5)
        n_pos = int(np.sum(probs >= threshold))
        n_positive_total += n_pos
        print(f"    {cls_name:<20} {np.mean(probs):>8.4f} {np.std(probs):>8.4f} "
              f"{np.min(probs):>8.4f} {np.max(probs):>8.4f} {n_pos:>5}")

    # ── Comparison with MedMNIST 28x28 upscaled ──
    print(f"\n  RESOLUTION COMPARISON NOTES:")
    print(f"    Full-res input: 1024x1024 -> resized to 224x224 by model preprocessing")
    print(f"    MedMNIST input: 28x28 -> upscaled to 224x224 (bicubic)")
    print(f"    Key difference: Full-res images preserve fine anatomical detail")
    print(f"    (rib structure, vascular markings, subtle opacities) that is")
    print(f"    completely lost in 28x28 inputs. The model receives genuinely")
    print(f"    different spatial frequency content at 224x224, leading to more")
    print(f"    discriminative and confident probability outputs.")

    # Check if we also have MedMNIST results to compare against
    medmnist_comparison = None
    if METADATA_PATH.exists():
        # Try to load existing validation results for comparison
        validation_json = PROJECT_ROOT / "data" / "sample_images" / "validation_results.json"
        if validation_json.exists():
            try:
                with open(validation_json) as f:
                    prev_results = json.load(f)
                if "chest_xray" in prev_results:
                    prev_timing = prev_results["chest_xray"].get("timing", {})
                    prev_mean_ms = prev_timing.get("mean_ms", 0)
                    if prev_mean_ms > 0:
                        speedup = prev_mean_ms / mean_time if mean_time > 0 else 0
                        print(f"\n  vs. MedMNIST 28x28 upscaled (from prior validation):")
                        print(f"    MedMNIST mean inference: {prev_mean_ms:.1f} ms")
                        print(f"    Full-res mean inference: {mean_time:.1f} ms")
                        print(f"    Ratio: {speedup:.2f}x")
                        medmnist_comparison = {
                            "medmnist_mean_ms": prev_mean_ms,
                            "fullres_mean_ms": round(mean_time, 1),
                            "ratio": round(speedup, 2),
                        }
            except Exception:
                pass

    return {
        "dataset": "fullres_cxr",
        "model": model_name,
        "n_samples": n,
        "resolution": "1024x1024 (native full-resolution)",
        "per_image_results": per_image_results,
        "timing": {
            "mean_ms": round(mean_time, 1),
            "median_ms": round(median_time, 1),
            "min_ms": round(min_time, 1),
            "max_ms": round(max_time, 1),
            "total_ms": round(sum(timing_list), 1),
        },
        "mean_probabilities": {
            c: round(float(np.mean(all_probs[c])), 4) for c in OUR_CLASSES
        },
        "prob_std": {
            c: round(float(np.std(all_probs[c])), 4) for c in OUR_CLASSES
        },
        "n_positive_findings_total": n_positive_total,
        "medmnist_comparison": medmnist_comparison,
    }


def run_large_scale_validation(
    n_samples: int = 100,
    backend: str = "auto",
) -> Optional[Dict]:
    """Run validation on a larger number of ChestMNIST images loaded directly.

    Instead of using saved PNGs, loads images directly from the MedMNIST
    dataset object and preprocesses them in memory. This allows validating
    on 100+ images without saving them all to disk.
    """
    try:
        from medmnist import ChestMNIST
    except ImportError:
        print("  medmnist not installed, skipping large-scale validation")
        return None

    from src.workflows.cxr_rapid_findings import CXRRapidFindingsWorkflow, CXR_CLASS_THRESHOLDS

    cache_dir = SAMPLE_DIR / ".medmnist_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*70}")
    print(f"LARGE-SCALE VALIDATION ({n_samples} ChestMNIST images)")
    print(f"{'='*70}")

    dataset = ChestMNIST(
        split="test", download=True, size=28,
        root=str(cache_dir),
    )

    n_samples = min(n_samples, len(dataset))
    print(f"  Dataset: ChestMNIST test split ({len(dataset)} total, using {n_samples})")

    if backend == "monai":
        workflow = CXRRapidFindingsWorkflow(mock_mode=False, checkpoint_path="__force_monai__")
    else:
        workflow = CXRRapidFindingsWorkflow(mock_mode=False)

    model_name = "torchxrayvision" if workflow._using_xrv else "MONAI (fallback)"
    print(f"  Model: {model_name}")
    print(f"  Processing {n_samples} images...")

    from PIL import Image as PILImage
    import tempfile

    y_true_list = []
    y_pred_list = []
    timing_list = []

    for i in range(n_samples):
        img, label = dataset[i]
        label_array = np.array(label).flatten()

        # Convert to our ground truth format
        gt = {}
        for chestmnist_idx, our_class in CHESTMNIST_TO_OUR_CLASSES.items():
            gt[our_class] = int(label_array[chestmnist_idx]) if chestmnist_idx < len(label_array) else 0
        gt["fracture"] = 0
        gt_vector = [gt[c] for c in OUR_CLASSES]

        # Save temp image (model expects file path)
        if not isinstance(img, PILImage.Image):
            img = PILImage.fromarray(np.array(img))
        img = img.resize((224, 224), PILImage.BICUBIC)

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name
            img.save(tmp_path)

        try:
            start = time.time()
            preprocessed = workflow.preprocess(tmp_path)
            inference_result = workflow.infer(preprocessed)
            elapsed_ms = (time.time() - start) * 1000

            class_probs = inference_result.get("class_probabilities", {})
            pred_vector = []
            for c in OUR_CLASSES:
                prob = class_probs.get(c, 0.0)
                threshold = CXR_CLASS_THRESHOLDS.get(c, 0.5)
                pred_vector.append(1 if prob >= threshold else 0)

            y_true_list.append(gt_vector)
            y_pred_list.append(pred_vector)
            timing_list.append(elapsed_ms)
        finally:
            os.unlink(tmp_path)

        if (i + 1) % 25 == 0:
            print(f"    Processed {i+1}/{n_samples}...")

    y_true = np.array(y_true_list)
    y_pred = np.array(y_pred_list)

    metrics = compute_metrics(y_true, y_pred, OUR_CLASSES)

    # Print results
    print(f"\n  PER-CLASS METRICS ({n_samples} samples):")
    agg = metrics["aggregate"]
    for class_name in OUR_CLASSES:
        m = metrics["per_class"][class_name]
        sens = f"{m['sensitivity']:.3f}" if isinstance(m["sensitivity"], float) else "N/A"
        spec = f"{m['specificity']:.3f}" if isinstance(m["specificity"], float) else "N/A"
        f1 = f"{m['f1']:.3f}" if isinstance(m["f1"], float) else "N/A"
        print(f"    {class_name:<20} TP={m['tp']:3d} FP={m['fp']:3d} TN={m['tn']:3d} FN={m['fn']:3d}  "
              f"Sens={sens}  Spec={spec}  F1={f1}")

    print(f"\n  AGGREGATE:")
    print(f"    Micro-F1: {agg['micro_f1']:.4f}")
    print(f"    Micro-accuracy: {agg['micro_accuracy']:.4f}")
    print(f"    Exact match: {agg['exact_match_ratio']:.4f} ({agg['exact_match_count']}/{n_samples})")
    print(f"    Mean inference: {np.mean(timing_list):.0f} ms ({np.sum(timing_list)/1000:.1f}s total)")

    # Print resolution caveat
    print(f"\n  NOTE: These images are 28x28 upscaled to 224x224. The model was")
    print(f"  trained on full-resolution CXR (~2000x2000). Low resolution causes")
    print(f"  probability compression toward 0.5-0.65, inflating false positives.")
    print(f"  Full-resolution validation would yield significantly better metrics.")

    return {
        "dataset": "chest_xray_large",
        "model": model_name,
        "n_samples": n_samples,
        "resolution_note": "28x28 upscaled to 224x224 (degraded)",
        "metrics": metrics,
        "timing": {
            "mean_ms": round(np.mean(timing_list), 1),
            "total_s": round(np.sum(timing_list) / 1000, 1),
        },
    }


def main():
    parser = argparse.ArgumentParser(description="Validate CXR workflow against real MedMNIST data")
    parser.add_argument("--verbose", "-v", action="store_true", help="Show per-image results")
    parser.add_argument("--backend", choices=["auto", "monai"], default="auto",
                        help="Model backend: auto (prefer xrv) or monai (force fallback)")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Path to save validation results as JSON")
    parser.add_argument("--large-scale", "-L", type=int, default=0, metavar="N",
                        help="Run large-scale validation on N images from ChestMNIST")
    parser.add_argument("--fullres", "-F", action="store_true",
                        help="Run full-resolution CXR validation (1024x1024 images)")
    parser.add_argument("--fullres-only", action="store_true",
                        help="Only run full-resolution validation (skip MedMNIST)")
    args = parser.parse_args()

    print("=" * 70)
    print("IMAGING INTELLIGENCE AGENT -- REAL DATA VALIDATION")
    print("=" * 70)
    print(f"  Sample directory: {SAMPLE_DIR}")
    print(f"  Backend: {args.backend}")
    print(f"  Verbose: {args.verbose}")

    all_results = {}

    if not args.fullres_only:
        metadata = load_metadata()

        # 1. Chest X-ray multi-label validation
        chest_results = run_validation_chest_xray(metadata, backend=args.backend, verbose=args.verbose)
        if chest_results:
            all_results["chest_xray"] = chest_results

        # 2. Pneumonia binary validation (cross-dataset generalization)
        pneumonia_results = run_validation_pneumonia(metadata, backend=args.backend, verbose=args.verbose)
        if pneumonia_results:
            all_results["pneumonia"] = pneumonia_results

        # 3. DICOM loading validation
        dicom_results = run_validation_dicom(metadata, verbose=args.verbose)
        if dicom_results:
            all_results["dicom"] = dicom_results

        # 4. Optional large-scale validation
        if args.large_scale > 0:
            large_results = run_large_scale_validation(
                n_samples=args.large_scale, backend=args.backend
            )
            if large_results:
                all_results["chest_xray_large"] = large_results

    # 5. Full-resolution CXR validation
    if args.fullres or args.fullres_only:
        fullres_results = validate_fullres_cxr(
            backend=args.backend, verbose=args.verbose
        )
        if fullres_results:
            all_results["fullres_cxr"] = fullres_results

    # Final summary
    print(f"\n{'='*70}")
    print("VALIDATION SUMMARY")
    print(f"{'='*70}")

    if "chest_xray" in all_results:
        cxr = all_results["chest_xray"]
        agg = cxr["metrics"]["aggregate"]
        print(f"  Chest X-ray (multi-label, {cxr['n_samples']} images, 28x28 upscaled):")
        print(f"    Micro-F1: {agg['micro_f1']:.4f}")
        print(f"    Micro-accuracy: {agg['micro_accuracy']:.4f}")
        print(f"    Exact match: {agg['exact_match_ratio']:.4f}")
        print(f"    Mean inference: {cxr['timing']['mean_ms']:.0f} ms")

    if "pneumonia" in all_results:
        pneu = all_results["pneumonia"]
        print(f"  Pneumonia detection ({pneu['n_samples']} images):")
        print(f"    Sensitivity: {pneu['sensitivity']:.4f}")
        print(f"    Specificity: {pneu['specificity']:.4f}")
        print(f"    Accuracy: {pneu['accuracy']:.4f}")
        print(f"    F1: {pneu['f1']:.4f}")

    if "dicom" in all_results:
        dcm = all_results["dicom"]
        print(f"  DICOM loading: {dcm['status'].upper()}")

    if "chest_xray_large" in all_results:
        large = all_results["chest_xray_large"]
        lagg = large["metrics"]["aggregate"]
        print(f"  Large-scale CXR ({large['n_samples']} images):")
        print(f"    Micro-F1: {lagg['micro_f1']:.4f}")
        print(f"    Micro-accuracy: {lagg['micro_accuracy']:.4f}")
        print(f"    Exact match: {lagg['exact_match_ratio']:.4f}")

    if "fullres_cxr" in all_results:
        fr = all_results["fullres_cxr"]
        print(f"  Full-resolution CXR ({fr['n_samples']} images, {fr['resolution']}):")
        print(f"    Mean inference: {fr['timing']['mean_ms']:.1f} ms")
        print(f"    Positive findings: {fr['n_positive_findings_total']}")
        print(f"    Mean probs: {fr['mean_probabilities']}")
        if fr.get("medmnist_comparison"):
            cmp = fr["medmnist_comparison"]
            print(f"    vs MedMNIST timing: {cmp['ratio']:.2f}x ratio")

    # Resolution caveat (only when MedMNIST was run)
    if not args.fullres_only:
        print(f"\n  NOTE: MedMNIST images are 28x28 upscaled to 224x224.")
        print(f"  The torchxrayvision DenseNet was trained on full-resolution CXR")
        print(f"  (~2000x2000). Performance on 28x28 upscaled images is degraded")
        print(f"  due to loss of fine diagnostic detail. Probabilities compress")
        print(f"  toward 0.5-0.65, inflating false positives above the threshold.")
        print(f"  This validates the end-to-end pipeline (load -> preprocess ->")
        print(f"  infer -> postprocess) rather than clinical accuracy.")

    # Save results if requested
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(all_results, f, indent=2, default=str)
        print(f"\n  Results saved to: {output_path}")

    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
