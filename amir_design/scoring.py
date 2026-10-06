"""Verified inference adapter for the published shRNAI models.

Scores are model outputs for shRNAI's supplied hairpin, not calibrated
knockdown percentages or predictions validated for another backbone.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
from typing import Any

import numpy as np


SHRNAI_COMMIT = "69f8576e027b8036eaeb2d039c40a2ccdbabf373"
SHRNAI_URL = "https://github.com/ParkSJ-91/shRNAI"
EXPECTED_SHA256 = {
    "shRNAI/module_simple.py": "b07a46728948099888915c1f38eb515030ab01bf3868a682374eafc231db4b96",
    "models/22nt.h5": "5650a6b618afc0c6993c6ec3c84c1d50716902f2fc42462948a3e9f897556471",
    "models/pri.h5": "6eb73b65a93fd2bee613d2b0ae06dc131e761705ef3e71aec9bf5c1f7a62dcdd",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_legacy_config(config: dict[str, Any]) -> dict[str, Any]:
    """Convert legacy integral float pooling dimensions in memory only.

    The saved Keras 2.3.0-tf models use 9.0/12.0 for pool sizes and strides.
    Their numeric values and weights remain unchanged. Nonintegral dimensions
    are rejected rather than rounded or otherwise altering the architecture.
    """
    normalized = deepcopy(config)

    def visit(value: Any) -> None:
        if isinstance(value, dict):
            for name, child in value.items():
                if name in {"pool_size", "strides"} and isinstance(child, list):
                    for index, dimension in enumerate(child):
                        if isinstance(dimension, float):
                            if not dimension.is_integer():
                                raise ValueError(f"Nonintegral {name} dimension: {dimension}")
                            child[index] = int(dimension)
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)

    visit(normalized)
    return normalized


class ShRNAIScorer:
    """Score 22-nt target-sense DNA windows with unchanged upstream weights."""

    def __init__(self, repo_path: Path, *, batch_size: int = 64) -> None:
        self.repo_path = Path(repo_path).resolve()
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        self.batch_size = batch_size
        hashes: dict[str, str] = {}
        for relative, expected in EXPECTED_SHA256.items():
            path = self.repo_path / relative
            if not path.is_file():
                raise FileNotFoundError(
                    f"Missing shRNAI artifact: {path}. Run the repository setup script."
                )
            actual = _sha256(path)
            if actual != expected:
                raise ValueError(f"shRNAI artifact checksum mismatch: {relative}")
            hashes[relative] = actual

        if (self.repo_path / ".git").exists():
            head = subprocess.run(
                ["git", "-C", str(self.repo_path), "rev-parse", "HEAD"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
            if head != SHRNAI_COMMIT:
                raise ValueError(f"Expected shRNAI commit {SHRNAI_COMMIT}; found {head}")

        self.provenance = {
            "repository": SHRNAI_URL,
            "commit": SHRNAI_COMMIT,
            "sha256": hashes,
            "scoring_scaffold": "shRNAI supplied pri-miR-30 template (113 nt)",
            "score_interpretation": "uncalibrated model activity score",
            "compatibility_adjustment": "integral float pooling dimensions normalized in memory",
        }
        source_path = self.repo_path / "shRNAI/module_simple.py"
        spec = importlib.util.spec_from_file_location("_verified_shrnai", source_path)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Cannot import verified shRNAI helper: {source_path}")
        self._helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self._helper)
        self._models: tuple[Any, Any] | None = None

    def encode_targets(
        self, target_sequences: list[str]
    ) -> tuple[list[str], list[str], np.ndarray, np.ndarray]:
        """Use the original encoder, maintaining candidate order and identity."""
        guides: list[str] = []
        hairpins: list[str] = []
        mature_inputs: list[np.ndarray] = []
        pri_inputs: list[np.ndarray] = []
        for sequence in target_sequences:
            if not isinstance(sequence, str) or len(sequence) != 22:
                raise ValueError("shRNAI inputs must be 22-nt target-sense DNA strings")
            if set(sequence) - set("ACGT"):
                raise ValueError("shRNAI inputs require uppercase A/C/G/T without ambiguous bases")
            upstream_guides, upstream_hairpins, mature, pri = self._helper.convert(sequence)
            if len(upstream_guides) != 1 or len(upstream_hairpins) != 1:
                raise RuntimeError("Unexpected candidate count from verified shRNAI encoder")
            guides.extend(upstream_guides)
            hairpins.extend(upstream_hairpins)
            mature_inputs.append(mature)
            pri_inputs.append(pri)
        if not guides:
            return (
                [],
                [],
                np.empty((0, 22, 4, 1), dtype=np.float32),
                np.empty((0, 56, 9, 1), dtype=np.float32),
            )
        return (
            guides,
            hairpins,
            np.concatenate(mature_inputs).astype(np.float32),
            np.concatenate(pri_inputs).astype(np.float32),
        )

    def _load_models(self) -> tuple[Any, Any]:
        if self._models is None:
            import h5py
            import tensorflow as tf

            def load(relative: str) -> Any:
                path = self.repo_path / relative
                with h5py.File(path, "r") as artifact:
                    stored_config = artifact.attrs["model_config"]
                if isinstance(stored_config, bytes):
                    stored_config = stored_config.decode("utf-8")
                config = normalize_legacy_config(json.loads(stored_config))
                with tf.device("/CPU:0"):
                    model = tf.keras.models.model_from_json(json.dumps(config))
                    model.load_weights(str(path))
                return model

            self._models = (load("models/pri.h5"), load("models/22nt.h5"))
        return self._models

    def score(self, target_sequences: list[str]) -> list[float]:
        """Return one raw final shRNAI score for every supplied candidate."""
        _, _, mature, pri = self.encode_targets(target_sequences)
        if not target_sequences:
            return []
        import tensorflow as tf

        pri_model, mature_model = self._load_models()
        scores: list[float] = []
        with tf.device("/CPU:0"):
            for start in range(0, len(target_sequences), self.batch_size):
                stop = start + self.batch_size
                processing = pri_model(pri[start:stop], training=False)
                output = mature_model([mature[start:stop], processing], training=False)
                values = np.asarray(output).reshape(-1)
                if not np.isfinite(values).all():
                    raise RuntimeError("shRNAI returned a nonfinite activity score")
                scores.extend(float(value) for value in values)
        if len(scores) != len(target_sequences):
            raise RuntimeError("shRNAI returned an unexpected number of scores")
        return scores
