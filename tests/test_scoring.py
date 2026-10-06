"""Regression checks against upstream encoding and published notebook results."""

from pathlib import Path

import numpy as np
import pytest

from amir_design.scoring import ShRNAIScorer, normalize_legacy_config


REPO = Path(__file__).resolve().parents[1] / "tools" / "shrnai"
TARGETS = ["ACAAGATGATGTTTGAAACTAT", "ACCAGATGTTAGTGACAATGAA"]
GUIDES = ["ATAGTTTCAAACATCATCTTGT", "TTCATTGTCACTAACATCTGGT"]


@pytest.fixture(scope="module")
def scorer():
    return ShRNAIScorer(REPO)


def test_upstream_encoding_parity_and_strand(scorer):
    guides, hairpins, mature, pri = scorer.encode_targets(TARGETS)
    assert guides == GUIDES
    assert all(len(hairpin) == 113 for hairpin in hairpins)
    assert mature.shape == (2, 22, 4, 1)
    assert pri.shape == (2, 56, 9, 1)
    for index, target in enumerate(TARGETS):
        upstream = scorer._helper.convert(target)
        assert hairpins[index] == upstream[1][0]
        np.testing.assert_array_equal(mature[index], upstream[2][0])
        np.testing.assert_array_equal(pri[index], upstream[3][0])


def test_real_cpu_scores_reproduce_notebook_and_preserve_order(scorer):
    first = scorer.score(TARGETS)
    second = scorer.score(TARGETS)
    assert len(first) == 2
    assert np.isfinite(first).all()
    np.testing.assert_allclose(first, [0.7972, 0.7761], atol=0.00015, rtol=0)
    np.testing.assert_allclose(first, second, atol=1e-7, rtol=0)
    np.testing.assert_allclose(scorer.score(TARGETS[::-1]), first[::-1], atol=1e-6, rtol=0)


@pytest.mark.parametrize("target", ["A" * 21, "A" * 23, "N" * 22, "a" * 22, "U" * 22])
def test_reject_invalid_model_input(scorer, target):
    with pytest.raises(ValueError):
        scorer.score([target])


def test_legacy_normalization_is_exact_and_does_not_mutate_source():
    original = {"config": {"pool_size": [9.0, 1], "strides": [12.0, 1]}}
    normalized = normalize_legacy_config(original)
    assert normalized == {"config": {"pool_size": [9, 1], "strides": [12, 1]}}
    assert isinstance(original["config"]["pool_size"][0], float)
    assert isinstance(normalized["config"]["pool_size"][0], int)
    with pytest.raises(ValueError):
        normalize_legacy_config({"pool_size": [9.5, 1]})


def test_empty_candidate_list_is_empty(scorer):
    assert scorer.score([]) == []
