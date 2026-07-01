"""TensorRT model optimization for MONAI clinical workflows.

Converts PyTorch models to TensorRT for 2-5x inference speedup on
NVIDIA GPUs. Supports FP16 and INT8 precision with dynamic input
shapes for variable CT/MRI volume dimensions.

Free SDK, available on DGX Spark ARM64.
"""

import hashlib
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from loguru import logger

from .triton_config import MONAI_MODEL_CONFIGS

# ── Dynamic input shape profiles for medical imaging ──
# Each model has min/opt/max shapes to handle variable volume sizes.
DYNAMIC_SHAPE_PROFILES: Dict[str, Dict[str, List[int]]] = {
    "segresnet_hemorrhage": {
        "min": [1, 1, 128, 128, 128],
        "opt": [1, 1, 256, 256, 256],
        "max": [1, 1, 512, 512, 512],
    },
    "unest_brain_ms": {
        "min": [1, 1, 64, 64, 64],
        "opt": [1, 1, 96, 96, 96],
        "max": [1, 1, 128, 128, 128],
    },
    "retinanet_lung_nodule": {
        "min": [1, 1, 64, 64, 64],
        "opt": [1, 1, 128, 128, 128],
        "max": [1, 1, 256, 256, 256],
    },
    "densenet121_cxr": {
        "min": [1, 1, 128, 128],
        "opt": [1, 1, 224, 224],
        "max": [1, 1, 512, 512],
    },
}

# ── Estimated speedups by model architecture and precision ──
SPEEDUP_ESTIMATES: Dict[str, Dict[str, float]] = {
    "segresnet_hemorrhage": {
        "fp16": 3.0,
        "int8": 4.5,
        "fp32": 1.0,
    },
    "unest_brain_ms": {
        "fp16": 2.5,
        "int8": 3.5,
        "fp32": 1.0,
    },
    "retinanet_lung_nodule": {
        "fp16": 3.0,
        "int8": 4.0,
        "fp32": 1.0,
    },
    "densenet121_cxr": {
        "fp16": 4.0,
        "int8": 5.5,
        "fp32": 1.0,
    },
}


class TensorRTOptimizer:
    """TensorRT model optimization for MONAI clinical workflows.

    Converts PyTorch models to TensorRT for 2-5x inference speedup.
    Supports FP16 and INT8 precision with dynamic input shapes for
    variable CT/MRI volume dimensions.

    Free SDK, available on DGX Spark ARM64.
    """

    def __init__(
        self,
        cache_dir: str = "data/cache/tensorrt",
        precision: str = "fp16",
    ):
        self.cache_dir = Path(cache_dir)
        self.precision = precision  # "fp32", "fp16", "int8"
        self._available = self._check_tensorrt()

        if self._available:
            logger.info(
                f"TensorRT optimizer ready (precision={precision}, "
                f"cache={self.cache_dir})"
            )
        else:
            logger.info(
                "TensorRT not available. Models will use standard "
                "PyTorch inference. Install with: pip install torch-tensorrt"
            )

    def _check_tensorrt(self) -> bool:
        """Check if torch_tensorrt is available."""
        try:
            import torch_tensorrt  # noqa: F401

            return True
        except ImportError:
            return False

    @property
    def available(self) -> bool:
        """Whether TensorRT optimization is available."""
        return self._available

    # ── Optimization ──────────────────────────────────────────────────

    def optimize_model(
        self,
        model_path: str,
        model_name: str,
        input_shapes: Optional[Dict[str, List[int]]] = None,
    ) -> str:
        """Convert a PyTorch model to TensorRT.

        Args:
            model_path: Path to PyTorch .pt model file.
            model_name: Key in MONAI_MODEL_CONFIGS (e.g., "segresnet_hemorrhage").
            input_shapes: Dict with 'min', 'opt', 'max' shape lists.
                          If None, uses DYNAMIC_SHAPE_PROFILES for the model.

        Returns:
            Path to the optimized TensorRT engine file.

        Raises:
            RuntimeError: If TensorRT is not available or optimization fails.
        """
        if not self._available:
            raise RuntimeError(
                "TensorRT not available. Install torch-tensorrt to enable "
                "model optimization."
            )

        import torch
        import torch_tensorrt

        model_path_obj = Path(model_path)
        if not model_path_obj.exists():
            raise FileNotFoundError(f"Model file not found: {model_path}")

        # Resolve input shapes
        shapes = input_shapes or DYNAMIC_SHAPE_PROFILES.get(model_name)
        if not shapes:
            raise ValueError(
                f"No input shape profile for '{model_name}'. Provide "
                f"input_shapes manually."
            )

        # Generate cache key from model path + shapes + precision
        cache_key = self._cache_key(model_path, shapes)
        output_path = self.cache_dir / f"{model_name}_{cache_key}.ts"

        # Return cached version if available
        if output_path.exists():
            logger.info(
                f"Using cached TensorRT model: {output_path.name}"
            )
            return str(output_path)

        logger.info(
            f"Optimizing '{model_name}' with TensorRT "
            f"(precision={self.precision})"
        )
        start = time.time()

        # Load PyTorch model
        model = torch.jit.load(str(model_path_obj))
        model.eval()
        model.cuda()

        # Build TensorRT compile spec
        compile_inputs = [
            torch_tensorrt.Input(
                min_shape=shapes["min"],
                opt_shape=shapes["opt"],
                max_shape=shapes["max"],
                dtype=torch.float32,
            )
        ]

        enabled_precisions = {torch.float32}
        if self.precision in ("fp16", "int8"):
            enabled_precisions.add(torch.float16)
        if self.precision == "int8":
            enabled_precisions.add(torch.int8)

        # Compile
        trt_model = torch_tensorrt.compile(
            model,
            inputs=compile_inputs,
            enabled_precisions=enabled_precisions,
            workspace_size=1 << 30,  # 1 GB workspace
            truncate_long_and_double=True,
        )

        # Save optimized model
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        torch.jit.save(trt_model, str(output_path))

        # Write metadata
        meta_path = output_path.with_suffix(".json")
        meta = {
            "model_name": model_name,
            "source_model": str(model_path),
            "precision": self.precision,
            "input_shapes": shapes,
            "optimization_time_sec": round(time.time() - start, 1),
            "output_path": str(output_path),
        }
        meta_path.write_text(json.dumps(meta, indent=2))

        elapsed = time.time() - start
        logger.info(
            f"TensorRT optimization complete for '{model_name}' "
            f"in {elapsed:.1f}s -> {output_path.name}"
        )
        return str(output_path)

    def load_optimized(self, model_name: str) -> Any:
        """Load a previously optimized TensorRT model.

        Args:
            model_name: Key in MONAI_MODEL_CONFIGS.

        Returns:
            The loaded TorchScript TensorRT model, or None if not found.

        Raises:
            RuntimeError: If TensorRT/PyTorch is not available.
        """
        if not self._available:
            raise RuntimeError("TensorRT not available")

        import torch

        # Find cached engine
        pattern = f"{model_name}_*.ts"
        matches = sorted(self.cache_dir.glob(pattern))
        if not matches:
            logger.warning(
                f"No cached TensorRT model found for '{model_name}'"
            )
            return None

        latest = matches[-1]
        logger.info(f"Loading TensorRT model: {latest.name}")
        model = torch.jit.load(str(latest))
        model.eval()
        model.cuda()
        return model

    # ── Query ─────────────────────────────────────────────────────────

    def is_optimized(self, model_name: str) -> bool:
        """Check if an optimized version exists in cache.

        Args:
            model_name: Key in MONAI_MODEL_CONFIGS.

        Returns:
            True if a cached TensorRT engine exists for this model.
        """
        if not self.cache_dir.exists():
            return False

        pattern = f"{model_name}_*.ts"
        return any(self.cache_dir.glob(pattern))

    def get_speedup_estimate(self, model_name: str) -> Dict[str, Any]:
        """Return estimated speedup for a model based on architecture.

        Args:
            model_name: Key in MONAI_MODEL_CONFIGS.

        Returns:
            Dict with speedup estimates per precision level, input shape
            profile, and the currently configured precision.
        """
        estimates = SPEEDUP_ESTIMATES.get(model_name, {})
        shapes = DYNAMIC_SHAPE_PROFILES.get(model_name, {})
        config = MONAI_MODEL_CONFIGS.get(model_name, {})

        return {
            "model_name": model_name,
            "display_name": config.get("display_name", model_name),
            "tensorrt_available": self._available,
            "current_precision": self.precision,
            "estimated_speedup": estimates.get(self.precision, 1.0),
            "speedup_by_precision": estimates,
            "input_shape_profile": shapes,
            "is_cached": self.is_optimized(model_name),
        }

    def get_all_speedup_estimates(self) -> List[Dict[str, Any]]:
        """Return speedup estimates for all known MONAI models."""
        return [
            self.get_speedup_estimate(name)
            for name in MONAI_MODEL_CONFIGS
        ]

    def list_cached_models(self) -> List[Dict[str, Any]]:
        """List all cached TensorRT models with metadata.

        Returns:
            List of dicts with model name, path, and metadata.
        """
        if not self.cache_dir.exists():
            return []

        cached: List[Dict[str, Any]] = []
        for meta_file in sorted(self.cache_dir.glob("*.json")):
            try:
                meta = json.loads(meta_file.read_text())
                engine_path = Path(meta.get("output_path", ""))
                meta["engine_exists"] = engine_path.exists()
                cached.append(meta)
            except (json.JSONDecodeError, KeyError) as e:
                logger.debug(f"Skipping invalid metadata {meta_file}: {e}")

        return cached

    def clear_cache(self, model_name: Optional[str] = None) -> int:
        """Remove cached TensorRT engines.

        Args:
            model_name: If provided, only clear cache for this model.
                        If None, clear all cached engines.

        Returns:
            Number of files removed.
        """
        if not self.cache_dir.exists():
            return 0

        if model_name:
            pattern_ts = f"{model_name}_*.ts"
            pattern_json = f"{model_name}_*.json"
            files = list(self.cache_dir.glob(pattern_ts)) + list(
                self.cache_dir.glob(pattern_json)
            )
        else:
            files = list(self.cache_dir.glob("*.ts")) + list(
                self.cache_dir.glob("*.json")
            )

        count = 0
        for f in files:
            try:
                f.unlink()
                count += 1
            except OSError as e:
                logger.warning(f"Failed to remove {f}: {e}")

        logger.info(f"Cleared {count} cached TensorRT file(s)")
        return count

    # ── Status ────────────────────────────────────────────────────────

    def get_status(self) -> Dict[str, Any]:
        """Return comprehensive status of TensorRT optimization."""
        return {
            "available": self._available,
            "precision": self.precision,
            "cache_dir": str(self.cache_dir),
            "cached_models": len(self.list_cached_models()),
            "known_models": list(MONAI_MODEL_CONFIGS.keys()),
            "speedup_estimates": {
                name: SPEEDUP_ESTIMATES.get(name, {}).get(self.precision, 1.0)
                for name in MONAI_MODEL_CONFIGS
            },
        }

    # ── Internals ─────────────────────────────────────────────────────

    def _cache_key(
        self,
        model_path: str,
        shapes: Dict[str, List[int]],
    ) -> str:
        """Generate a deterministic cache key for a model + config."""
        h = hashlib.sha256()
        h.update(model_path.encode())
        h.update(self.precision.encode())
        h.update(json.dumps(shapes, sort_keys=True).encode())
        return h.hexdigest()[:12]
