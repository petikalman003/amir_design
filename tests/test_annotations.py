"""Strand-aware transcript regions and safe, exact target lookup."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from amir_design.annotations import (
    map_transcript_regions,
    merge_intervals,
    prepare_reference_annotations,
)
from amir_design.references import (
    AmbiguousGeneError,
    ReferenceError,
    prepare_reference_bundle,
)


def _fake_blast(command):
    if "-version" in command:
        return "makeblastdb: test"
    prefix = Path(command[command.index("-out") + 1])
    for suffix in (".nhr", ".nin", ".nsq"):
        Path(str(prefix) + suffix).write_bytes(b"index")
    return ""


def _gtf_line(transcript_id, gene_id, feature, start, end, strand="+"):
    return f'chr1\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene_id}"; transcript_id "{transcript_id}";\n'


@pytest.fixture
def reference(tmp_path):
    fasta = tmp_path / "input.fa"
    fasta.write_text(
        ">ENSG100.2|ENST100.4|Example\n" + "ACGT" * 7 + "AC\n"
        ">ENSG100.2|ENST101.5|Example\n" + "TGCA" * 7 + "TG\n"
        ">ENSG200.1|ENST200.1|Other\n" + "A" * 30 + "\n"
    )
    with patch("amir_design.references._run", side_effect=_fake_blast):
        bundle = prepare_reference_bundle(tmp_path / "bundle", species="human", source_fasta=fasta)
    gtf = tmp_path / "input.gtf"
    plus = "".join(_gtf_line("ENST100.4", "ENSG100.2", feature, start, end) for feature, start, end in (
        ("exon", 100, 109), ("exon", 200, 219),
        ("CDS", 105, 109), ("CDS", 200, 211),
        ("start_codon", 105, 107), ("stop_codon", 212, 214),
        ("UTR", 100, 104), ("UTR", 215, 219),
    ))
    minus = "".join(_gtf_line("ENST101.5", "ENSG100.2", feature, start, end, "-") for feature, start, end in (
        ("exon", 100, 109), ("exon", 200, 219),
        ("CDS", 200, 214), ("CDS", 105, 109),
        ("stop_codon", 102, 104), ("UTR", 215, 219), ("UTR", 100, 101),
    ))
    gtf.write_text(plus + minus + _gtf_line("ENST200.1", "ENSG200.1", "exon", 300, 329))
    return bundle, gtf


def test_plus_minus_splicing_and_stop_codon_are_in_mature_coordinates(reference):
    bundle, gtf = reference
    base_manifest = (bundle.root / "manifest.json").read_bytes()
    prepare_reference_annotations(bundle, source_gtf=gtf)
    regions = bundle.regions_for_gene("ENSG100.2")
    assert regions["ENST100.4"] == {"cds": [(6, 25)], "utr5": [(1, 5)], "utr3": [(26, 30)]}
    assert regions["ENST101.5"] == {"cds": [(6, 28)], "utr5": [(1, 5)], "utr3": [(29, 30)]}
    assert bundle.regions_for_gene("ENSG200.1")["ENST200.1"] == {"cds": [], "utr5": [], "utr3": []}
    assert (bundle.root / "manifest.json").read_bytes() == base_manifest
    bundle.verify_integrity()
    bundle.verify_annotations()


def test_true_region_gaps_are_not_joined():
    regions = map_transcript_regions(
        transcript_id="test", length=30, strand="+", exons=[(100, 109), (200, 219)],
        features={"cds": [(105, 107), (200, 211)], "stop": [(212, 214)]},
    )
    assert regions["cds"] == [(6, 8), (11, 25)]
    assert merge_intervals([(1, 3), (4, 7), (9, 10)]) == [(1, 7), (9, 10)]


@pytest.mark.parametrize("strand,features,expected", [
    ("+", {"cds": [(105, 109), (200, 211)], "stop": [(212, 214)],
            "utr": [(100, 104), (212, 219)]},
     {"cds": [(6, 25)], "utr5": [(1, 5)], "utr3": [(26, 30)]}),
    ("-", {"cds": [(200, 214), (105, 109)], "stop": [(102, 104)],
            "utr": [(215, 219), (100, 104)]},
     {"cds": [(6, 28)], "utr5": [(1, 5)], "utr3": [(29, 30)]}),
])
def test_gencode_utr_includes_stop_bases_which_belong_only_to_cds(strand, features, expected):
    # The minus case follows the actual Xkr4 GENCODE M36 convention: generic
    # 3prime UTR reaches through all three separately annotated stop bases.
    regions = map_transcript_regions(
        transcript_id="test", length=30, strand=strand,
        exons=[(100, 109), (200, 219)], features=features,
    )
    assert regions == expected


def test_gencode_stop_only_utr_is_empty_and_split_stop_is_spliced():
    # A stop codon can straddle an exon junction. Its mature bases still form
    # one codon, and a UTR consisting only of those bases contains no UTR site.
    regions = map_transcript_regions(
        transcript_id="test", length=30, strand="+", exons=[(100, 109), (200, 219)],
        features={"cds": [(100, 108)], "stop": [(109, 109), (200, 201)],
                  "utr": [(109, 109), (200, 201)]},
    )
    assert regions == {"cds": [(1, 12)], "utr5": [], "utr3": []}


@pytest.mark.parametrize("exons,features,message", [
    ([(100, 108), (200, 219)], {}, "length"),
    ([(100, 110), (110, 128)], {}, "Overlapping"),
    ([(100, 109), (200, 219)], {"cds": [(105, 204)]}, "Non-exonic"),
    ([(100, 109), (200, 219)], {"cds": [(105, 109)], "utr": [(105, 106)]}, "UTR overlaps"),
    ([(100, 109), (200, 219)], {"stop": [(105, 107)]}, "without CDS"),
])
def test_invalid_mappings_fail_closed(exons, features, message):
    with pytest.raises(ReferenceError, match=message):
        map_transcript_regions(transcript_id="test", length=30, strand="+", exons=exons, features=features)


def test_identifier_lookup_gene_transcript_versions_and_symbol(reference):
    bundle, _ = reference
    for text in ("ENSG100.2", "ENSG100", "Example", "example"):
        resolved = bundle.resolve_identifier(text)
        assert resolved["input_type"] == "gene"
        assert resolved["gene_id"] == "ENSG100.2"
        assert resolved["matching_transcript_ids"] == ["ENST100.4", "ENST101.5"]
    for text in ("ENST100.4", "ENST100"):
        resolved = bundle.resolve_identifier(text)
        assert resolved["input_type"] == "transcript"
        assert resolved["matching_transcript_ids"] == ["ENST100.4"]
        assert resolved["transcript_ids"] == ["ENST100.4", "ENST101.5"]
    for text in ("ENST100.3", "ENSG100.1", "missing"):
        with pytest.raises(ReferenceError, match="not found"):
            bundle.resolve_identifier(text)


def test_ambiguous_symbols_and_stable_ids_are_rejected(tmp_path):
    fasta = tmp_path / "ambiguous.fa"
    fasta.write_text(
        ">ENSG100.1|ENST100.1|DUP\nACGT\n"
        ">ENSG100.2|ENST100.2|dup\nACGT\n"
    )
    with patch("amir_design.references._run", side_effect=_fake_blast):
        bundle = prepare_reference_bundle(tmp_path / "bundle", species="human", source_fasta=fasta)
    for text in ("ENSG100", "ENST100", "DUP", "Dup"):
        with pytest.raises(AmbiguousGeneError, match="multiple"):
            bundle.resolve_identifier(text)
    assert bundle.resolve_identifier("ENSG100.1")["gene_id"] == "ENSG100.1"


def test_missing_annotations_are_actionable_and_identifier_lookup_still_works(reference):
    bundle, _ = reference
    assert bundle.annotation_provenance == {}
    assert bundle.resolve_identifier("ENST100")["gene_id"] == "ENSG100.2"
    with pytest.raises(ReferenceError, match="setup.sh --references"):
        bundle.regions_for_gene("ENSG100.2")


def test_sidecar_reuse_and_sha_validation_do_not_download(reference):
    bundle, gtf = reference
    sha = hashlib.sha256(gtf.read_bytes()).hexdigest()
    first = prepare_reference_annotations(bundle, source_gtf=gtf, expected_sha256=sha)
    with patch("amir_design.annotations._download") as download:
        assert prepare_reference_annotations(bundle, source_gtf=gtf, expected_sha256=sha) == first
    download.assert_not_called()
    assert first["source_sha256"] == sha
    assert first["reference"]["source_sha256"] == bundle.manifest["source_sha256"]
    with pytest.raises(ReferenceError, match="SHA-256"):
        prepare_reference_annotations(bundle, source_gtf=gtf, expected_sha256="0" * 64)
    gtf.write_text(gtf.read_text() + "#changed\n")
    with pytest.raises(ReferenceError, match="differs"):
        prepare_reference_annotations(bundle, source_gtf=gtf)


@pytest.mark.parametrize("change,message", [
    (lambda text: text.replace('gene_id "ENSG100.2"', 'gene_id "WRONG"', 1), "gene assignment"),
    (lambda text: text.replace("\t100\t109\t", "\t100\t108\t", 1), "length"),
    (lambda text: "".join(line for line in text.splitlines(True) if "ENST101.5" not in line), "lacks matching"),
])
def test_bad_gtf_does_not_publish_partial_sidecar(reference, change, message):
    bundle, gtf = reference
    gtf.write_text(change(gtf.read_text()))
    with pytest.raises(ReferenceError, match=message):
        prepare_reference_annotations(bundle, source_gtf=gtf)
    assert not (bundle.root / "annotations").exists()
    assert not list(bundle.root.glob(".annotations-*"))


def test_compatible_gtf_extra_transcripts_are_ignored_and_recorded(reference):
    bundle, gtf = reference
    gtf.write_text(gtf.read_text() + _gtf_line("not_in_fasta", "other_gene", "exon", 900, 950))
    manifest = prepare_reference_annotations(bundle, source_gtf=gtf)
    assert manifest["ignored_gtf_transcripts_absent_from_fasta"] == 1


def test_annotation_tampering_or_other_reference_binding_is_rejected(reference):
    bundle, gtf = reference
    prepare_reference_annotations(bundle, source_gtf=gtf)
    database = bundle.root / "annotations" / "regions.sqlite3"
    with database.open("r+b") as handle:
        handle.seek(100)
        handle.write(b"changed")
    with pytest.raises(ReferenceError, match="integrity"):
        bundle.verify_annotations()
    manifest_path = bundle.root / "annotations" / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["reference"]["source_sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ReferenceError, match="different transcript"):
        _ = bundle.annotation_provenance


def test_gzipped_gtf_is_supported(reference):
    bundle, gtf = reference
    compressed = gtf.with_suffix(".gtf.gz")
    with gzip.open(compressed, "wb") as source:
        source.write(gtf.read_bytes())
    prepare_reference_annotations(bundle, source_gtf=compressed)
    assert bundle.regions_for_gene("ENSG100.2")["ENST100.4"]["cds"] == [(6, 25)]


def test_sidecar_manifest_cannot_reference_outside_files(reference):
    bundle, gtf = reference
    prepare_reference_annotations(bundle, source_gtf=gtf)
    manifest_path = bundle.root / "annotations" / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["sqlite"] = "../transcripts.sqlite3"
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ReferenceError, match="Invalid region"):
        bundle.regions_for_gene("ENSG100.2")


def test_official_gtf_checksum_failure_never_publishes(reference):
    bundle, gtf = reference
    # Bind this synthetic reference to the official source solely to exercise
    # download/checksum control flow, not to claim it is a real GENCODE bundle.
    from amir_design.references import GENCODE_SOURCES
    bundle.manifest["source"] = GENCODE_SOURCES["human"]
    archive = gzip.compress(gtf.read_bytes())
    md5 = hashlib.md5(archive).hexdigest()

    def download(url, destination):
        if url.endswith("MD5SUMS"):
            destination.write_text(f"{md5}  gencode.v47.annotation.gtf.gz\n")
        else:
            destination.write_bytes(archive + b"corrupt")

    with patch("amir_design.annotations._download", side_effect=download):
        with pytest.raises(ReferenceError, match="authoritative MD5"):
            prepare_reference_annotations(bundle)
    assert not (bundle.root / "annotations").exists()
