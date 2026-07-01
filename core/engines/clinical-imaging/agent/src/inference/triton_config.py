"""Triton Inference Server integration for MONAI clinical workflows.

Manages model deployment on NVIDIA Triton Inference Server with dynamic
model loading/unloading for GPU memory management, health monitoring,
and graceful fallback to direct PyTorch inference when Triton is
unavailable.

BSD 3-Clause (Triton) + Free SDK (TensorRT). Both available on DGX
Spark ARM64.
"""

import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from loguru import logger

# ── Model repository configurations for the 4 MONAI clinical models ──
# Each config maps to a Triton model_repository/ layout with config.pbtxt
MONAI_MODEL_CONFIGS: Dict[str, Dict[str, Any]] = {
    "segresnet_hemorrhage": {
        "display_name": "SegResNet - CT Head Hemorrhage",
        "framework": "pytorch",
        "max_batch_size": 1,
        "input_name": "INPUT__0",
        "input_dtype": "FP32",
        "input_shape": [1, 1, 256, 256, 256],
        "output_name": "OUTPUT__0",
        "output_dtype": "FP32",
        "output_shape": [1, 5, 256, 256, 256],  # 5-class hemorrhage
        "instance_group_count": 1,
        "dynamic_batching_delay_ms": 0,
        "preferred_batch_sizes": [1],
        "model_file": "model.pt",
        "optimization": {
            "tensorrt": True,
            "precision": "fp16",
            "estimated_speedup": 3.0,
        },
    },
    "unest_brain_ms": {
        "display_name": "UNEST - MRI Brain MS Lesion",
        "framework": "pytorch",
        "max_batch_size": 1,
        "input_name": "INPUT__0",
        "input_dtype": "FP32",
        "input_shape": [1, 1, 96, 96, 96],
        "output_name": "OUTPUT__0",
        "output_dtype": "FP32",
        "output_shape": [1, 2, 96, 96, 96],  # Binary MS lesion
        "instance_group_count": 1,
        "dynamic_batching_delay_ms": 0,
        "preferred_batch_sizes": [1],
        "model_file": "model.pt",
        "optimization": {
            "tensorrt": True,
            "precision": "fp16",
            "estimated_speedup": 2.5,
        },
    },
    "retinanet_lung_nodule": {
        "display_name": "RetinaNet - CT Lung Nodule Detection",
        "framework": "pytorch",
        "max_batch_size": 1,
        "input_name": "INPUT__0",
        "input_dtype": "FP32",
        "input_shape": [1, 1, 128, 128, 128],
        "output_name": "OUTPUT__0",
        "output_dtype": "FP32",
        "output_shape": [1, -1, 7],  # Variable detections: [x,y,z,w,h,d,score]
        "instance_group_count": 1,
        "dynamic_batching_delay_ms": 0,
        "preferred_batch_sizes": [1],
        "model_file": "model.pt",
        "optimization": {
            "tensorrt": True,
            "precision": "fp16",
            "estimated_speedup": 3.0,
        },
    },
    "densenet121_cxr": {
        "display_name": "DenseNet-121 - Chest X-Ray Classification",
        "framework": "pytorch",
        "max_batch_size": 8,
        "input_name": "INPUT__0",
        "input_dtype": "FP32",
        "input_shape": [1, 1, 224, 224],  # 2D, fastest to optimize
        "output_name": "OUTPUT__0",
        "output_dtype": "FP32",
        "output_shape": [1, 14],  # 14 CheXpert pathology classes
        "instance_group_count": 1,
        "dynamic_batching_delay_ms": 50,
        "preferred_batch_sizes": [1, 4, 8],
        "model_file": "model.pt",
        "optimization": {
            "tensorrt": True,
            "precision": "fp16",
            "estimated_speedup": 4.0,
        },
    },
}


class TritonModelManager:
    """Manages model deployment on NVIDIA Triton Inference Server.

    Provides:
    - Model repository configuration generation
    - Dynamic model loading/unloading for memory management
    - Health checking and status monitoring
    - TensorRT optimization wrapper
    - Fallback to direct PyTorch inference when Triton unavailable

    BSD 3-Clause (Triton) + Free SDK (TensorRT). Both available on
    DGX Spark ARM64.
    """

    def __init__(self, triton_url: str = "localhost:8000", enabled: bool = True):
        self.triton_url = triton_url
        self.enabled = enabled
        self._client: Any = None
        self._last_health_check: float = 0
        self._health_cache: Optional[bool] = None
        self._health_ttl: float = 10.0  # seconds
        self._init_client()

    def _init_client(self) -> None:
        """Try to import tritonclient; graceful fallback if unavailable."""
        if not self.enabled:
            logger.info("Triton integration disabled by configuration")
            return

        try:
            import tritonclient.http as httpclient

            self._client = httpclient.InferenceServerClient(
                url=self.triton_url,
                verbose=False,
            )
            logger.info(f"Triton client initialized at {self.triton_url}")
        except ImportError:
            logger.warning(
                "tritonclient not installed. Triton features disabled. "
                "Install with: pip install tritonclient[http]"
            )
            self.enabled = False
        except Exception as e:
            logger.warning(f"Failed to initialize Triton client: {e}")
            self.enabled = False

    # ── Health & Status ──────────────────────────────────────────────

    def is_healthy(self) -> bool:
        """Check Triton server health with TTL-based caching."""
        if not self.enabled or self._client is None:
            return False

        now = time.time()
        if (
            self._health_cache is not None
            and (now - self._last_health_check) < self._health_ttl
        ):
            return self._health_cache

        try:
            healthy = self._client.is_server_live() and self._client.is_server_ready()
            self._health_cache = healthy
            self._last_health_check = now
            if healthy:
                logger.debug("Triton server is healthy")
            else:
                logger.warning("Triton server is not ready")
            return healthy
        except Exception as e:
            logger.warning(f"Triton health check failed: {e}")
            self._health_cache = False
            self._last_health_check = now
            return False

    def get_server_metadata(self) -> Dict[str, Any]:
        """Return Triton server metadata (name, version, extensions)."""
        if not self.is_healthy():
            return {"status": "unavailable", "enabled": self.enabled}

        try:
            meta = self._client.get_server_metadata()
            return {
                "status": "healthy",
                "name": meta.get("name", "triton"),
                "version": meta.get("version", "unknown"),
                "extensions": meta.get("extensions", []),
            }
        except Exception as e:
            logger.error(f"Failed to get Triton metadata: {e}")
            return {"status": "error", "error": str(e)}

    # ── Model Management ─────────────────────────────────────────────

    def list_models(self) -> List[Dict[str, Any]]:
        """List loaded models with status."""
        if not self.is_healthy():
            return []

        try:
            repo_index = self._client.get_model_repository_index()
            models = []
            for entry in repo_index:
                models.append(
                    {
                        "name": entry.get("name", ""),
                        "version": entry.get("version", ""),
                        "state": entry.get("state", "UNKNOWN"),
                        "reason": entry.get("reason", ""),
                    }
                )
            return models
        except Exception as e:
            logger.error(f"Failed to list Triton models: {e}")
            return []

    def load_model(self, model_name: str) -> bool:
        """Dynamically load a model into Triton.

        Args:
            model_name: Name matching a directory in the model repository.

        Returns:
            True if model was loaded successfully.
        """
        if not self.is_healthy():
            logger.warning(
                f"Cannot load model '{model_name}': Triton not available"
            )
            return False

        try:
            self._client.load_model(model_name)
            logger.info(f"Loaded model '{model_name}' into Triton")
            return True
        except Exception as e:
            logger.error(f"Failed to load model '{model_name}': {e}")
            return False

    def unload_model(self, model_name: str) -> bool:
        """Unload a model to free GPU memory.

        Args:
            model_name: Name of the loaded model to unload.

        Returns:
            True if model was unloaded successfully.
        """
        if not self.is_healthy():
            logger.warning(
                f"Cannot unload model '{model_name}': Triton not available"
            )
            return False

        try:
            self._client.unload_model(model_name)
            logger.info(f"Unloaded model '{model_name}' from Triton")
            return True
        except Exception as e:
            logger.error(f"Failed to unload model '{model_name}': {e}")
            return False

    def is_model_ready(self, model_name: str) -> bool:
        """Check if a specific model is ready for inference."""
        if not self.is_healthy():
            return False

        try:
            return self._client.is_model_ready(model_name)
        except Exception:
            return False

    def get_model_config(self, model_name: str) -> Dict[str, Any]:
        """Get model configuration from Triton.

        Args:
            model_name: Name of the model.

        Returns:
            Model configuration dict, or static config from
            MONAI_MODEL_CONFIGS if Triton is unavailable.
        """
        # Try live Triton first
        if self.is_healthy():
            try:
                config = self._client.get_model_config(model_name)
                return dict(config)
            except Exception as e:
                logger.debug(
                    f"Triton config unavailable for '{model_name}': {e}"
                )

        # Fall back to static config
        if model_name in MONAI_MODEL_CONFIGS:
            return dict(MONAI_MODEL_CONFIGS[model_name])

        return {}

    # ── Inference ─────────────────────────────────────────────────────

    def infer(
        self,
        model_name: str,
        inputs: Dict[str, np.ndarray],
        outputs: Optional[List[str]] = None,
    ) -> Dict[str, np.ndarray]:
        """Run inference through Triton.

        Args:
            model_name: Name of the loaded model.
            inputs: Dict mapping input tensor names to numpy arrays.
            outputs: Optional list of output tensor names to request.
                     If None, all outputs are returned.

        Returns:
            Dict mapping output tensor names to numpy arrays.

        Raises:
            RuntimeError: If Triton is unavailable or inference fails.
        """
        if not self.is_healthy():
            raise RuntimeError(
                f"Triton not available for inference on '{model_name}'. "
                "Use direct PyTorch inference as fallback."
            )

        try:
            import tritonclient.http as httpclient

            # Build Triton input objects
            triton_inputs = []
            for name, data in inputs.items():
                inp = httpclient.InferInput(
                    name, list(data.shape), _numpy_to_triton_dtype(data.dtype)
                )
                inp.set_data_from_numpy(data)
                triton_inputs.append(inp)

            # Build Triton output requests
            triton_outputs = None
            if outputs:
                triton_outputs = [
                    httpclient.InferRequestedOutput(name) for name in outputs
                ]

            # Run inference
            result = self._client.infer(
                model_name=model_name,
                inputs=triton_inputs,
                outputs=triton_outputs,
            )

            # Extract output tensors
            output_dict: Dict[str, np.ndarray] = {}
            config = MONAI_MODEL_CONFIGS.get(model_name, {})
            output_names = outputs or [config.get("output_name", "OUTPUT__0")]
            for name in output_names:
                try:
                    output_dict[name] = result.as_numpy(name)
                except Exception:
                    pass

            logger.debug(
                f"Triton inference on '{model_name}' completed: "
                f"{len(output_dict)} output(s)"
            )
            return output_dict

        except ImportError:
            raise RuntimeError("tritonclient not installed")
        except Exception as e:
            logger.error(f"Triton inference failed for '{model_name}': {e}")
            raise RuntimeError(f"Triton inference failed: {e}") from e

    # ── Config Generation ─────────────────────────────────────────────

    def generate_config_pbtxt(self, model_name: str) -> str:
        """Generate a Triton config.pbtxt for a MONAI model.

        Args:
            model_name: Key in MONAI_MODEL_CONFIGS.

        Returns:
            String contents for config.pbtxt, or empty string if
            model_name is not recognized.
        """
        config = MONAI_MODEL_CONFIGS.get(model_name)
        if not config:
            logger.warning(f"No config template for model '{model_name}'")
            return ""

        input_dims = ", ".join(str(d) for d in config["input_shape"][1:])
        output_dims = ", ".join(str(d) for d in config["output_shape"][1:])

        batch_sizes_str = ", ".join(
            str(b) for b in config.get("preferred_batch_sizes", [1])
        )

        pbtxt = f"""\
name: "{model_name}"
platform: "pytorch_libtorch"
max_batch_size: {config["max_batch_size"]}

input [
  {{
    name: "{config["input_name"]}"
    data_type: TYPE_{config["input_dtype"]}
    dims: [ {input_dims} ]
  }}
]

output [
  {{
    name: "{config["output_name"]}"
    data_type: TYPE_{config["output_dtype"]}
    dims: [ {output_dims} ]
  }}
]

instance_group [
  {{
    count: {config["instance_group_count"]}
    kind: KIND_GPU
  }}
]

dynamic_batching {{
  max_queue_delay_microseconds: {config["dynamic_batching_delay_ms"] * 1000}
  preferred_batch_size: [ {batch_sizes_str} ]
}}
"""
        return pbtxt

    def setup_model_repository(
        self, repo_dir: str, model_names: Optional[List[str]] = None
    ) -> List[str]:
        """Create Triton model repository directory structure.

        Generates config.pbtxt files for each model. Model weights
        (model.pt) must be placed manually or via TensorRTOptimizer.

        Args:
            repo_dir: Root directory for the model repository.
            model_names: List of model names to set up. If None, all
                         MONAI models are configured.

        Returns:
            List of model directories created.
        """
        repo_path = Path(repo_dir)
        names = model_names or list(MONAI_MODEL_CONFIGS.keys())
        created: List[str] = []

        for name in names:
            model_dir = repo_path / name / "1"
            model_dir.mkdir(parents=True, exist_ok=True)

            config_pbtxt = self.generate_config_pbtxt(name)
            if config_pbtxt:
                config_path = repo_path / name / "config.pbtxt"
                config_path.write_text(config_pbtxt)
                created.append(str(repo_path / name))
                logger.info(f"Created model config: {config_path}")

        return created

    # ── Status ────────────────────────────────────────────────────────

    def get_status(self) -> Dict[str, Any]:
        """Return comprehensive status of the Triton integration."""
        status: Dict[str, Any] = {
            "enabled": self.enabled,
            "triton_url": self.triton_url,
            "healthy": self.is_healthy(),
            "models_loaded": [],
            "available_configs": list(MONAI_MODEL_CONFIGS.keys()),
        }
        if status["healthy"]:
            status["models_loaded"] = self.list_models()
            status["server"] = self.get_server_metadata()
        return status


# ── Utilities ─────────────────────────────────────────────────────────


def _numpy_to_triton_dtype(np_dtype: np.dtype) -> str:
    """Convert numpy dtype to Triton dtype string."""
    mapping = {
        np.float32: "FP32",
        np.float16: "FP16",
        np.int32: "INT32",
        np.int64: "INT64",
        np.int8: "INT8",
        np.uint8: "UINT8",
        np.bool_: "BOOL",
    }
    return mapping.get(np_dtype.type, "FP32")


def get_model_configs() -> Dict[str, Dict[str, Any]]:
    """Return all MONAI model configurations.

    Convenience function for external modules that need model metadata
    without instantiating a TritonModelManager.
    """
    return dict(MONAI_MODEL_CONFIGS)
