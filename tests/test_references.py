"""Reference contracts: preserve annotation, fail closed, and assign one gene."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from amir_design.references import (
    AmbiguousGeneError,
    ReferenceBundle,
    ReferenceError,
    _HTTPSOnlyRedirect,
    _download,
    parse_transcript_header,
    prepare_reference_bundle,
)


CDS = "ACGTCAGTTCGATGACCTAGTCAA"


def _legacy_fasta() -> str:
    return (
        f">GENE1.2|TX1.3|Example\nAAA{CDS}TTT\n"
        f">GENE1.2|TX2.8|Example\nCCC{CDS}AAAA\n"
        ">GENE1.2|TX3.1|Example\nACGTACGTACGTACGTACGTACGT\n"
        f">GENE2.1|TX4.2|Other\nTTT{CDS}CCC\n"
    )


def _mock_makeblastdb(command: list[str]) -> str:
    if "-version" in command:
        return "makeblastdb: 2.17.0+"
    prefix = Path(command[command.index("-out") + 1])
    for suffix in (".nhr", ".nin", ".nsq"):
        Path(str(prefix) + suffix).write_bytes(b"test index")
    return ""


class HeaderTests(unittest.TestCase):
    def test_http_sources_and_redirects_are_rejected_before_fetch(self):
        with self.assertRaisesRegex(ReferenceError, "HTTPS"):
            _download("http://example.org/reference.fa", Path("unused.fa"))
        with self.assertRaisesRegex(ReferenceError, "HTTPS"):
            _HTTPSOnlyRedirect().redirect_request(
                None, None, 302, "Found", {}, "http://example.org/reference.fa"
            )

    def test_gencode_header_retains_versions_and_biotype(self):
        header = (
            ">ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|"
            "OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|lncRNA|"
        )
        record = parse_transcript_header(header)
        self.assertEqual(record["transcript_id"], "ENST00000456328.2")
        self.assertEqual(record["gene_id"], "ENSG00000223972.5")
        self.assertEqual(record["gene_name"], "DDX11L1")
        self.assertEqual(record["biotype"], "lncRNA")

    def test_gencode_protein_coding_region_headers_have_correct_biotype(self):
        headers = (
            ">ENST00000335674.1|ENSG00000185962.1|OTTHUMG00000012397.3|OTTHUMT00000034517.3|LCE3A-201|LCE3A|270|CDS:1-270|",
            ">ENST00000641515.2|ENSG00000186092.7|OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-201|OR4F5|2618|UTR5:1-60|CDS:61-1041|UTR3:1042-2618|",
            ">ENSMUST00000070533.5|ENSMUSG00000051951.6|OTTMUSG00000026353|OTTMUST00000065166.1|Xkr4-201|Xkr4|3634|UTR5:1-150|CDS:151-2094|UTR3:2095-3634|",
        )
        for header in headers:
            with self.subTest(header=header):
                record = parse_transcript_header(header)
                self.assertEqual(record["biotype"], "protein_coding")
                self.assertEqual(record["transcript_id"], header[1:].split("|")[0])

    def test_gencode_region_tags_must_include_valid_cds_coordinates(self):
        prefix = ">ENST1.1|ENSG1.1|OTTG|OTTT|EX-201|EX|24|"
        for tags in ("UTR5:1-24|", "CDS:1-25|", "CDS:2-1|", "CDS:1-24|CDS:1-24|", "CDS:one-many|"):
            with self.subTest(tags=tags):
                with self.assertRaises(ReferenceError):
                    parse_transcript_header(prefix + tags)

    def test_mouse_and_par_y_headers(self):
        mouse = parse_transcript_header(
            ">ENSMUST00000193812.1|ENSMUSG00000102693.2|OTT1|OTT2|4933401J01Rik-201|4933401J01Rik|24|protein_coding|"
        )
        par_y = parse_transcript_header(
            ">ENST000001.4_PAR_Y|ENSG000001.3_PAR_Y|OTT1|OTT2|TEST-201|TEST|24|protein_coding|"
        )
        self.assertEqual(mouse["gene_id"], "ENSMUSG00000102693.2")
        self.assertEqual(par_y["transcript_id"], "ENST000001.4_PAR_Y")

    def test_legacy_header_and_unannotated_rejection(self):
        record = parse_transcript_header(">gene.1|transcript.2|Symbol")
        self.assertEqual(record["biotype"], "unspecified")
        with self.assertRaisesRegex(ReferenceError, "lack gene annotation"):
            parse_transcript_header(">transcript.2")


class BundleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)
        self.source = self.directory / "custom.fa"
        self.source.write_text(_legacy_fasta(), encoding="ascii")
        self.root = self.directory / "bundle"

    def tearDown(self):
        self.temporary.cleanup()

    def build(self):
        with patch("amir_design.references._run", side_effect=_mock_makeblastdb):
            return prepare_reference_bundle(self.root, species="human", source_fasta=self.source)

    def test_preserves_all_isoforms_full_sequences_and_normalized_ids(self):
        bundle = self.build()
        transcripts = bundle.transcripts_for_gene("GENE1.2")
        self.assertEqual(set(transcripts), {"TX1.3", "TX2.8", "TX3.1"})
        self.assertEqual(transcripts["TX1.3"], f"AAA{CDS}TTT")
        self.assertEqual(bundle.transcript_gene_map()["TX4.2"], "GENE2.1")
        self.assertIn(">TX1.3\n", bundle.fasta_path.read_text())
        self.assertNotIn("GENE1.2|", bundle.fasta_path.read_text())
        self.assertEqual(bundle.manifest["transcript_count"], 4)
        bundle.verify_integrity()

    def test_gzipped_gencode_and_biotypes_are_not_filtered(self):
        self.source = self.directory / "custom.fa.gz"
        with gzip.open(self.source, "wt") as handle:
            for tx_id, biotype in (
                ("ENST1.3", "protein_coding"),
                ("ENST2.5", "nonsense_mediated_decay"),
                ("ENST3.1", "retained_intron"),
            ):
                handle.write(
                    f">{tx_id}|ENSG1.2|OTTG|OTTT|EXAMPLE-201|EXAMPLE|{len(CDS)}|{biotype}|\n{CDS}\n"
                )
        bundle = self.build()
        self.assertEqual(len(bundle.transcripts_for_gene("ENSG1.2")), 3)
        self.assertEqual(
            set(bundle.manifest["transcript_biotype_counts"]),
            {"protein_coding", "nonsense_mediated_decay", "retained_intron"},
        )

    def test_pc_reference_metadata_does_not_label_coordinates_as_biotype(self):
        self.source.write_text(
            f">ENST1.1|ENSG1.1|OTTG|OTTT|EX-201|EX|{len(CDS)}|UTR5:1-3|CDS:4-21|UTR3:22-24|\n{CDS}\n"
        )
        bundle = self.build()
        self.assertEqual(bundle.manifest["transcript_biotype_counts"], {"protein_coding": 1})
        self.assertEqual(bundle.transcripts_for_gene("ENSG1.1"), {"ENST1.1": CDS})

    def test_bad_checksum_and_duplicate_ids_do_not_publish(self):
        with self.assertRaisesRegex(ReferenceError, "SHA-256"):
            prepare_reference_bundle(
                self.root, species="human", source_fasta=self.source, expected_sha256="0" * 64
            )
        self.assertFalse(self.root.exists())
        self.source.write_text(_legacy_fasta() + f">GENE1.2|TX1.3|Example\n{CDS}\n")
        with self.assertRaisesRegex(ReferenceError, "Duplicate transcript"):
            self.build()
        self.assertFalse(self.root.exists())
        self.assertFalse(list(self.directory.glob(".bundle-*")))

    def test_declared_length_must_match(self):
        self.source.write_text(f">ENST1.1|ENSG1.1|OTTG|OTTT|EX-201|EX|100|protein_coding|\n{CDS}\n")
        with self.assertRaisesRegex(ReferenceError, "length disagrees"):
            self.build()
        self.assertFalse(self.root.exists())

    def test_reuses_only_identical_reference_and_checks_integrity(self):
        bundle = self.build()
        with patch("amir_design.references._run") as run:
            reused = prepare_reference_bundle(self.root, species="human", source_fasta=self.source)
        run.assert_not_called()
        self.assertEqual(reused.manifest, bundle.manifest)
        self.source.write_text(_legacy_fasta().replace("AAA", "CCC", 1))
        with self.assertRaisesRegex(ReferenceError, "differs"):
            prepare_reference_bundle(self.root, species="human", source_fasta=self.source)
        with bundle.fasta_path.open("r+") as handle:
            handle.write("!")
        with self.assertRaisesRegex(ReferenceError, "integrity"):
            bundle.verify_integrity()

    def test_missing_blast_index_is_actionable(self):
        bundle = self.build()
        missing = Path(str(bundle.blast_db_path) + ".nin")
        missing.unlink()
        with self.assertRaisesRegex(ReferenceError, "missing|Missing"):
            ReferenceBundle.open(self.root)

    def test_manifest_cannot_escape_bundle(self):
        self.build()
        manifest_path = self.root / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["fasta"] = "../custom.fa"
        manifest_path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ReferenceError, "outside"):
            ReferenceBundle.open(self.root)

    def test_unique_full_cds_gene_assignment_includes_nonmatching_isoforms(self):
        bundle = self.build()
        output = f"TX1.3\t1\t{len(CDS)}\t{len(CDS)}\t0\t0\t50\nTX2.8\t1\t{len(CDS)}\t{len(CDS)}\t0\t0\t50\nTX4.2\t1\t20\t20\t0\t0\t30\n"
        result = subprocess.CompletedProcess([], 0, stdout=output, stderr="")
        with patch("amir_design.references.subprocess.run", return_value=result) as run:
            gene = bundle.gene_for_cds(CDS)
        command = run.call_args.args[0]
        self.assertEqual(command[command.index("-strand") + 1], "plus")
        self.assertEqual(gene["gene_id"], "GENE1.2")
        self.assertEqual(gene["transcript_ids"], ["TX1.3", "TX2.8", "TX3.1"])
        self.assertEqual(gene["matching_transcript_ids"], ["TX1.3", "TX2.8"])

    def test_exact_matches_to_two_genes_are_rejected(self):
        bundle = self.build()
        output = "".join(
            f"{tx}\t1\t{len(CDS)}\t{len(CDS)}\t0\t0\t50\n" for tx in ("TX1.3", "TX4.2")
        )
        result = subprocess.CompletedProcess([], 0, stdout=output, stderr="")
        with patch("amir_design.references.subprocess.run", return_value=result):
            with self.assertRaisesRegex(AmbiguousGeneError, "GENE1.2.*GENE2.1"):
                bundle.gene_for_cds(CDS)

    def test_partial_matches_do_not_assign_a_gene(self):
        bundle = self.build()
        result = subprocess.CompletedProcess(
            [], 0, stdout="TX1.3\t2\t20\t19\t0\t0\t30\n", stderr=""
        )
        with patch("amir_design.references.subprocess.run", return_value=result):
            with self.assertRaisesRegex(ReferenceError, "No exact full-length"):
                bundle.gene_for_cds(CDS)

    def test_official_download_requires_authoritative_checksum(self):
        compressed = gzip.compress(self.source.read_bytes())
        good_md5 = hashlib.md5(compressed).hexdigest()

        def download(url, destination):
            if url.endswith("MD5SUMS"):
                destination.write_text(f"{good_md5}  gencode.v47.transcripts.fa.gz\n")
            else:
                destination.write_bytes(compressed + b"corruption")

        with patch("amir_design.references._download", side_effect=download):
            with self.assertRaisesRegex(ReferenceError, "authoritative MD5"):
                prepare_reference_bundle(self.root, species="human")
        self.assertFalse(self.root.exists())

    @unittest.skipUnless(
        shutil.which("makeblastdb") and shutil.which("blastn"), "NCBI BLAST+ is not on PATH"
    )
    def test_real_blast_maps_unique_cds_on_plus_strand(self):
        self.source.write_text(
            f">GENE1.2|TX1.3|Example\nAAA{CDS}TTT\n"
            ">GENE1.2|TX2.8|Example\nACGTACGTACGTACGTACGTACGT\n"
        )
        bundle = prepare_reference_bundle(self.root, species="human", source_fasta=self.source)
        self.assertEqual(bundle.gene_for_cds(CDS)["gene_id"], "GENE1.2")
        self.assertEqual(len(bundle.transcripts_for_gene("GENE1.2")), 2)

    @unittest.skipUnless(
        shutil.which("makeblastdb") and shutil.which("blastn"), "NCBI BLAST+ is not on PATH"
    )
    def test_real_blast_rejects_cross_gene_ambiguity_and_partial_cds(self):
        bundle = prepare_reference_bundle(self.root, species="human", source_fasta=self.source)
        with self.assertRaises(AmbiguousGeneError):
            bundle.gene_for_cds(CDS)
        with self.assertRaisesRegex(ReferenceError, "No exact full-length"):
            bundle.gene_for_cds(CDS + "GGGGGGGG")


if __name__ == "__main__":
    unittest.main()
