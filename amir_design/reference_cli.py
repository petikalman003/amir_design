"""Prepare versioned human/mouse references for the Streamlit pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from .references import ReferenceError, prepare_reference_bundle
from .annotations import prepare_reference_annotations


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build checked, reusable full-transcript FASTA/SQLite/BLAST reference bundles."
    )
    parser.add_argument("--species", choices=("human", "mouse", "both"), required=True)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("data/references"),
        help="Parent directory; bundles are stored beneath human/ and mouse/.",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        help="Local custom FASTA or FASTA.gz with GENCODE or gene|transcript|symbol headers.",
    )
    parser.add_argument(
        "--sha256", help="Optional expected SHA-256 of the original custom FASTA/archive."
    )
    parser.add_argument("--release", help="Custom reference release label (requires --fasta).")
    parser.add_argument("--assembly", help="Custom assembly label (requires --fasta).")
    parser.add_argument("--gtf", type=Path, help="Matching local GTF/GTF.gz for custom transcript regions.")
    parser.add_argument("--gtf-sha256", help="Optional expected SHA-256 of the original GTF/archive.")
    args = parser.parse_args(argv)
    if args.fasta and args.species == "both":
        parser.error(
            "--fasta requires a single species; prepare custom human and mouse inputs separately."
        )
    if not args.fasta and (args.release or args.assembly):
        parser.error("--release and --assembly apply to custom --fasta inputs.")
    if args.gtf and args.species == "both":
        parser.error("--gtf requires a single species.")
    if args.gtf_sha256 and not args.gtf:
        parser.error("--gtf-sha256 requires --gtf.")
    species_list = ("human", "mouse") if args.species == "both" else (args.species,)
    failed = False
    for species in species_list:
        print(f"Preparing {species} reference at {args.root / species}", flush=True)
        try:
            bundle = prepare_reference_bundle(
                args.root / species,
                species=species,
                source_fasta=args.fasta,
                expected_sha256=args.sha256,
                release=args.release,
                assembly=args.assembly,
            )
            if args.gtf or not args.fasta:
                print(f"Preparing matching {species} CDS/UTR annotations", flush=True)
                prepare_reference_annotations(
                    bundle, source_gtf=args.gtf, expected_sha256=args.gtf_sha256,
                )
        except ReferenceError as exc:
            print(f"{species}: {exc}", file=sys.stderr, flush=True)
            failed = True
            continue
        source = bundle.manifest["source"]
        print(
            f"Ready: {source['release']} / {source['assembly']}; {bundle.manifest['transcript_count']:,} transcripts, {bundle.manifest['gene_count']:,} genes; checksums verified.",
            flush=True,
        )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
