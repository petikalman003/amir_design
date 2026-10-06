import pytest

from amir_design.candidates import (
    DESIGN_RULES, FIGURE_S1_RULES, normalize_cds, passes_design_rules,
    reverse_complement, select_candidates,
)


TARGET = "CGGCTGCGCATATACATATATA"


def test_guide_orientation_seed_and_every_isoform_coverage():
    selection = select_candidates(TARGET, {"tx1": "ATG" + TARGET, "tx2": TARGET + "TAA"})
    assert len(selection.candidates) == 1
    candidate = selection.candidates[0]
    assert candidate.target_sequence == TARGET
    assert candidate.guide_dna == reverse_complement(TARGET)
    assert reverse_complement(candidate.guide_dna) == TARGET
    assert candidate.guide_dna[1:8] == reverse_complement(TARGET[14:21])
    assert candidate.isoform_sites == {"tx1": (4,), "tx2": (1,)}


def test_missing_site_in_one_isoform_rejects_single_guide():
    selection = select_candidates(TARGET, {"tx1": TARGET, "tx2": "ATG" * 30})
    assert selection.sequence_filtered_count == 1
    assert selection.rejected_isoform_count == 1
    assert selection.candidates == []


def test_duplicate_windows_keep_one_candidate_and_all_origins():
    selection = select_candidates(TARGET + TARGET, {"tx": TARGET + TARGET})
    matches = [candidate for candidate in selection.candidates if candidate.target_sequence == TARGET]
    assert len(matches) == 1
    assert matches[0].cds_starts == (1, 23)
    assert matches[0].isoform_sites["tx"] == (1, 23)


def test_n3_profiles_are_explicit_and_have_different_strand_semantics():
    target_with_a = TARGET[:2] + "A" + TARGET[3:]
    target_with_t = TARGET[:2] + "T" + TARGET[3:]
    assert not passes_design_rules(target_with_a, DESIGN_RULES)
    assert passes_design_rules(target_with_a, FIGURE_S1_RULES)
    assert passes_design_rules(target_with_t, DESIGN_RULES)
    assert not passes_design_rules(target_with_t, FIGURE_S1_RULES)
    assert reverse_complement(target_with_t)[19] == "A"


@pytest.mark.parametrize("text", ["", "N" * 22, "U" * 22, "ACGT", ">one\nACGT\n>two\nACGT"])
def test_unusable_inputs_fail_before_external_operations(text):
    with pytest.raises(ValueError):
        normalize_cds(text)


def test_plain_sequence_and_single_fasta_normalize_identically():
    assert normalize_cds("\n" + TARGET.lower() + "\t") == TARGET
    assert normalize_cds(">my_cds\n" + TARGET[:10] + "\n" + TARGET[10:]) == TARGET
