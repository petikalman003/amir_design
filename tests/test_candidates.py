import pytest

from amir_design.candidates import (
    DESIGN_RULES, FIGURE_S1_RULES, normalize_cds, passes_design_rules,
    reverse_complement, select_candidates, select_region_candidates,
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


def test_partial_coverage_keeps_at_least_one_isoform_and_reports_missing_ids():
    transcripts = {"tx1": TARGET, "tx2": "ATG" * 30}
    regions = {"tx1": [(1, 22)], "tx2": []}
    partial = select_region_candidates(transcripts, regions)
    assert partial.sequence_filtered_count == 1
    assert partial.rejected_isoform_count == 0
    assert partial.missing_region_transcript_ids == ("tx2",)
    candidate = partial.candidates[0]
    assert candidate.isoform_sites == {"tx1": (1,), "tx2": ()}
    assert candidate.covered_isoforms == ("tx1",)
    assert candidate.missed_isoforms == ("tx2",)
    assert candidate.coverage_fraction == 0.5
    complete = select_region_candidates(transcripts, regions, require_all_isoforms=True)
    assert complete.sequence_filtered_count == 1
    assert complete.rejected_isoform_count == 1
    assert complete.candidates == []


def test_gene_generation_deduplicates_all_transcript_origins():
    transcripts = {"tx1": TARGET + "NNN" + TARGET, "tx2": "NN" + TARGET}
    regions = {"tx1": [(1, 22), (26, 47)], "tx2": [(3, 24)]}
    selection = select_region_candidates(transcripts, regions, region="utr3")
    assert len(selection.candidates) == 1
    candidate = selection.candidates[0]
    assert candidate.source_sites == {"tx1": (1, 26), "tx2": (3,)}
    assert candidate.source_region == "utr3"
    assert candidate.cds_starts == ()
    assert candidate.coverage_fraction == 1


def test_transcript_generation_keeps_full_gene_denominator_and_cross_region_hits():
    other = "CCGCTGCGCATATACATATATA"
    transcripts = {"tx1": TARGET, "tx2": other + "NN" + TARGET, "tx3": "ATG" * 10}
    # tx2's selected CDS contains OTHER. TARGET there lies outside that region.
    regions = {"tx1": [(1, 22)], "tx2": [(1, 22)], "tx3": []}
    selected = select_region_candidates(
        transcripts, regions, generation_transcript_ids=["tx1"], region="cds",
    )
    assert len(selected.candidates) == 1
    candidate = selected.candidates[0]
    assert candidate.source_sites == {"tx1": (1,)}
    assert candidate.covered_isoforms == ("tx1", "tx2")
    assert candidate.missed_isoforms == ("tx3",)
    assert candidate.isoform_sites["tx2"] == (25,)
    assert candidate.coverage_fraction == pytest.approx(2 / 3)
    assert selected.target_isoforms == ("tx1", "tx2", "tx3")
    all_sources = select_region_candidates(transcripts, regions)
    assert {candidate.target_sequence for candidate in all_sources.candidates} == {TARGET, other}


def test_region_boundaries_exclude_partial_windows_and_include_last_full_window():
    assert select_region_candidates({"tx": TARGET}, {"tx": [(2, 22)]}).candidates == []
    assert select_region_candidates({"tx": TARGET}, {"tx": [(1, 21)]}).candidates == []
    assert select_region_candidates(
        {"tx": TARGET}, {"tx": [(1, 11), (12, 22)]},
    ).candidates == []
    selection = select_region_candidates({"tx": "N" + TARGET}, {"tx": [(2, 23)]})
    assert selection.candidates[0].source_sites == {"tx": (2,)}


def test_ambiguous_windows_only_are_skipped_and_overlapping_annotations_deduplicate():
    sequence = "N" + TARGET
    selection = select_region_candidates(
        {"tx": sequence}, {"tx": [(1, 23), (1, 23)]}, region="whole",
    )
    assert selection.ambiguous_windows_skipped == 1
    assert len(selection.candidates) == 1
    assert selection.candidates[0].source_sites == {"tx": (2,)}


@pytest.mark.parametrize("regions", [{"tx": [(0, 22)]}, {"tx": [(1, 23)]},
                                   {"tx": [(20, 1)]}])
def test_invalid_annotation_coordinates_fail_instead_of_silently_clipping(regions):
    with pytest.raises(ValueError, match="Invalid cds coordinates"):
        select_region_candidates({"tx": TARGET}, regions)


def test_unknown_transcript_or_region_cannot_silently_change_generation_scope():
    with pytest.raises(ValueError, match="Unknown generation transcript"):
        select_region_candidates({"tx": TARGET}, {"tx": [(1, 22)]},
                                 generation_transcript_ids=["missing"])
    with pytest.raises(ValueError, match="Unknown search region"):
        select_region_candidates({"tx": TARGET}, {"tx": [(1, 22)]}, region="intron")
