"""
Evidence Embedder - Generate vector embeddings for variant evidence.
"""
from typing import List, Optional
import numpy as np
from pathlib import Path
from loguru import logger

try:
    from sentence_transformers import SentenceTransformer
    HAS_SENTENCE_TRANSFORMERS = True
except ImportError:
    HAS_SENTENCE_TRANSFORMERS = False
    logger.warning("sentence-transformers not installed.")

from .vcf_parser import VariantEvidence


class EvidenceEmbedder:
    """
    Generate embeddings for variant evidence text summaries.

    Uses sentence-transformers models optimized for semantic search.
    Recommended models:
    - BAAI/bge-small-en-v1.5 (384 dim, fast, good quality)
    - BAAI/bge-base-en-v1.5 (768 dim, better quality)
    - all-MiniLM-L6-v2 (384 dim, very fast)
    """

    def __init__(
        self,
        model_name: str = "BAAI/bge-small-en-v1.5",
        cache_dir: Optional[Path] = None,
        device: Optional[str] = None,
        batch_size: int = 32,
    ):
        self.model_name = model_name
        self.batch_size = batch_size

        if not HAS_SENTENCE_TRANSFORMERS:
            raise ImportError(
                "sentence-transformers is required. "
                "Install with: pip install sentence-transformers"
            )

        logger.info(f"Loading embedding model: {model_name}")

        # Set cache directory if specified
        cache_folder = str(cache_dir) if cache_dir else None

        # Initialize model
        self.model = SentenceTransformer(
            model_name,
            cache_folder=cache_folder,
            device=device,
        )

        self.dimension = self.model.get_sentence_embedding_dimension()
        logger.info(f"Embedding dimension: {self.dimension}")

    def embed_text(self, text: str) -> np.ndarray:
        """Embed a single text string."""
        embedding = self.model.encode(
            text,
            normalize_embeddings=True,
            show_progress_bar=False,
        )
        return embedding

    def embed_texts(self, texts: List[str]) -> np.ndarray:
        """Embed multiple text strings."""
        embeddings = self.model.encode(
            texts,
            normalize_embeddings=True,
            batch_size=self.batch_size,
            show_progress_bar=len(texts) > 100,
        )
        return embeddings

    def embed_evidence(self, evidence: VariantEvidence) -> np.ndarray:
        """Embed a single variant evidence object."""
        return self.embed_text(evidence.text_summary)

    def embed_evidence_batch(self, evidence_list: List[VariantEvidence]) -> np.ndarray:
        """Embed a batch of variant evidence objects."""
        texts = [e.text_summary for e in evidence_list]
        return self.embed_texts(texts)

    def embed_query(self, query: str) -> np.ndarray:
        """
        Embed a search query.

        For BGE models, queries should be prefixed for better retrieval.
        """
        # BGE models benefit from instruction prefix for queries
        if "bge" in self.model_name.lower():
            query = f"Represent this sentence for searching relevant passages: {query}"

        return self.embed_text(query)


class CachedEmbedder(EvidenceEmbedder):
    """
    Embedder with disk caching for large datasets.
    Avoids re-computing embeddings for variants already processed.
    """

    def __init__(
        self,
        cache_path: Path,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.cache_path = Path(cache_path)
        self.cache_path.mkdir(parents=True, exist_ok=True)

        # In-memory cache for current session
        self._memory_cache: dict = {}

        # Load existing cache
        self._load_cache()

    def _cache_file(self) -> Path:
        """Path to cache file."""
        model_slug = self.model_name.replace("/", "_")
        return self.cache_path / f"embeddings_{model_slug}.npz"

    def _load_cache(self) -> None:
        """Load cached embeddings from disk."""
        cache_file = self._cache_file()
        if cache_file.exists():
            try:
                data = np.load(cache_file, allow_pickle=True)
                self._memory_cache = dict(data['cache'].item())
                logger.info(f"Loaded {len(self._memory_cache)} cached embeddings")
            except Exception as e:
                logger.warning(f"Could not load cache: {e}")
                self._memory_cache = {}

    def _save_cache(self) -> None:
        """Save cached embeddings to disk."""
        try:
            np.savez_compressed(
                self._cache_file(),
                cache=self._memory_cache
            )
            logger.debug(f"Saved {len(self._memory_cache)} embeddings to cache")
        except Exception as e:
            logger.warning(f"Could not save cache: {e}")

    def embed_evidence(self, evidence: VariantEvidence) -> np.ndarray:
        """Embed with caching."""
        cache_key = evidence.variant_id

        if cache_key in self._memory_cache:
            return self._memory_cache[cache_key]

        embedding = super().embed_evidence(evidence)
        self._memory_cache[cache_key] = embedding
        return embedding

    def embed_evidence_batch(self, evidence_list: List[VariantEvidence]) -> np.ndarray:
        """Embed batch with caching."""
        embeddings = []
        to_compute = []
        to_compute_idx = []

        # Check cache
        for i, evidence in enumerate(evidence_list):
            cache_key = evidence.variant_id
            if cache_key in self._memory_cache:
                embeddings.append((i, self._memory_cache[cache_key]))
            else:
                to_compute.append(evidence)
                to_compute_idx.append(i)

        # Compute missing embeddings
        if to_compute:
            texts = [e.text_summary for e in to_compute]
            new_embeddings = self.embed_texts(texts)

            for evidence, embedding, idx in zip(to_compute, new_embeddings, to_compute_idx):
                self._memory_cache[evidence.variant_id] = embedding
                embeddings.append((idx, embedding))

        # Sort by original index and return
        embeddings.sort(key=lambda x: x[0])
        return np.array([e[1] for e in embeddings])

    def flush_cache(self) -> None:
        """Save cache to disk."""
        self._save_cache()

    def clear_cache(self) -> None:
        """Clear all cached embeddings."""
        self._memory_cache = {}
        cache_file = self._cache_file()
        if cache_file.exists():
            cache_file.unlink()
        logger.info("Embedding cache cleared")
