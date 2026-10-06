from pathlib import Path
import os
import shutil

import pytest

from amir_design.offtargets import SeedMatchRunner, summarize_sites


GUIDE = "TTATAGAGCAAGAACACTGTTT"


def test_class_overlap_and_target_isoforms_do_not_inflate_counts():
    def hit(candidate, transcript, seed, start, end, anchor):
        return {
            "candidate_id": candidate, "transcript_id": transcript,
            "seed_class": seed, "start": start, "end": end,
            "width": end - start + 1, "g2_position": anchor,
            "target_sequence": "CTCTATAA",
        }
    mapping = {"T1": "TARGET", "T2": "TARGET", "O1": "OTHER", "O2": "OTHER"}
    rows = [
        hit("a", "T1", "mer8", 4, 11, 10),
        hit("a", "T2", "mer8", 4, 11, 10),
        hit("a", "O1", "mer7m8", 4, 10, 10),
        hit("a", "O1", "mer7A1", 5, 11, 10),
        hit("a", "O1", "mer8", 4, 11, 10),
        hit("a", "O2", "mer8", 4, 11, 10),
    ]
    summary, sites = summarize_sites(["a", "zero"], rows, mapping, "TARGET")
    assert summary["a"]["off_target_genes"] == 1
    assert summary["a"]["off_target_transcripts"] == 2
    assert summary["a"]["off_target_sites"] == 2
    assert summary["a"]["mer8_sites"] == 2
    assert summary["a"]["mer7m8_sites"] == 0
    assert summary["a"]["mer7A1_sites"] == 0
    assert "mer6_sites" not in summary["a"]
    assert all(value == 0 for value in summary["zero"].values())
    assert {row["transcript_id"] for row in sites} == {"O1", "O2"}


@pytest.fixture
def reference(tmp_path):
    records = {
        "T1": ("TARGET", "GGGCTCTATAAGGG"),
        "T2": ("TARGET", "GGGCTCTATAAGGG"),
        "O1": ("GENE_A", "GGGCTCTATAAGGG"),
        "O2": ("GENE_A", "AAACTCTATAAGGG"),
        "O3": ("GENE_B", "GGGTCTATAAGGG"),
        "O4": ("GENE_C", "GGGCTCTATAGGG"),
        "N1": ("GENE_D", "GGGTATAGAGGGG"),  # seed guide sequence, not its target
    }
    fasta = tmp_path / "reference.fa"
    fasta.write_text("".join(f">{tx} annotation description\n{seq}\n" for tx, (_, seq) in records.items()))
    return fasta, {tx: gene for tx, (gene, _) in records.items()}


@pytest.fixture
def runner():
    repo = Path(__file__).resolve().parents[1] / "tools/seedmatchr"
    local = Path("/workspace/.runtime/bin/Rscript")
    rscript = os.environ.get("AMIR_RSCRIPT") or (str(local) if local.is_file() else shutil.which("Rscript"))
    if not rscript or not repo.is_dir():
        pytest.skip("Run environment setup to provide pinned SeedMatchR and R/Biostrings")
    return SeedMatchRunner(repo, rscript)


def test_actual_seedmatchr_orientation_counts_shared_seeds_and_zero_hits(runner, reference):
    fasta, mapping = reference
    alternative_same_seed = GUIDE[:8] + "ACGTACGTACGTAC"
    no_hits = "CGGGGGGGACGTACGTACGTAC"
    result = runner.screen(
        {"dna": GUIDE, "rna": GUIDE.replace("T", "U"), "shared": alternative_same_seed, "zero": no_hits},
        fasta, "TARGET", mapping,
    )
    for candidate in ("dna", "rna", "shared"):
        assert result.summary[candidate] == {
            "off_target_genes": 3, "off_target_transcripts": 4, "off_target_sites": 4,
            "mer8_sites": 2, "mer7m8_sites": 1, "mer7A1_sites": 1,
        }
    assert all(value == 0 for value in result.summary["zero"].values())
    assert len(result.sites) == 12
    assert not {"T1", "T2", "N1"}.intersection(row["transcript_id"] for row in result.sites)
    assert result.provenance["engine"] == "SeedMatchR"
    assert result.provenance["r_version"]
    assert result.provenance["biostrings_version"]
    assert result.provenance["coordinates"].startswith("1-based inclusive")
    one = next(row for row in result.sites if row["candidate_id"] == "dna" and row["transcript_id"] == "O1")
    assert one["start"] == 4
    assert one["end"] == 11
    assert one["g2_position"] == 10
    assert one["target_sequence"] == "CTCTATAA"


def test_actual_seedmatchr_complete_zero_match_output(runner, tmp_path):
    fasta = tmp_path / "no-matches.fa"
    fasta.write_text(">target\nACGTACGTACGTACGTACGTACGT\n>other\nACGTACGTACGTACGTACGTACGT\n")
    result = runner.screen({"a": GUIDE}, fasta, "TARGET", {"target": "TARGET", "other": "OTHER"})
    assert result.sites == []
    assert all(value == 0 for value in result.summary["a"].values())
    assert "mer6_sites" not in result.summary["a"]


def test_mer6_counts_are_only_reported_when_explicitly_searched(runner, tmp_path):
    fasta = tmp_path / "six-mer-only.fa"
    fasta.write_text(">target\nGGGCTCTATAAGGG\n>other\nGGGTCTATAGGG\n")
    mapping = {"target": "TARGET", "other": "OTHER"}
    default = runner.screen({"a": GUIDE}, fasta, "TARGET", mapping)
    assert default.summary["a"]["off_target_genes"] == 0
    assert "mer6_sites" not in default.summary["a"]
    explicit = runner.screen(
        {"a": GUIDE}, fasta, "TARGET", mapping,
        seed_classes=("mer8", "mer7m8", "mer7A1", "mer6"),
    )
    assert explicit.summary["a"]["mer6_sites"] == 1
    assert explicit.summary["a"]["off_target_genes"] == 1
    assert explicit.summary["a"]["off_target_transcripts"] == 1
    assert explicit.summary["a"]["off_target_sites"] == 1
    assert explicit.sites[0]["seed_class"] == "mer6"
    assert explicit.sites[0]["g2_position"] == 9


def test_incomplete_mapping_is_rejected_before_starting_r(reference):
    fasta, mapping = reference
    mapping.pop("O2")
    with pytest.raises(ValueError, match="missing gene mappings.*O2"):
        SeedMatchRunner(Path("absent"), "absent-rscript").screen({"a": GUIDE}, fasta, "TARGET", mapping)


def test_guide_must_be_mature_22_nt(reference):
    fasta, mapping = reference
    with pytest.raises(ValueError, match="22-nt guide"):
        SeedMatchRunner(Path("absent")).screen({"a": "A" * 113}, fasta, "TARGET", mapping)
