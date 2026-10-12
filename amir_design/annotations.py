"""Checked GTF annotations mapped to mature-transcript sequence coordinates.

Region data are an atomic sidecar. Installing them does not rewrite an existing
full-transcript FASTA, SQLite database, BLAST index, or version-1 manifest.
"""

from __future__ import annotations

from contextlib import closing
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import tempfile

from .references import (
    GENCODE_SOURCES,
    ReferenceBundle,
    ReferenceError,
    _download,
    _hash_file,
    _official_md5,
)


ANNOTATION_SOURCES = {
    species: {
        "species": species,
        "release": source["release"],
        "assembly": source["assembly"],
        "gtf_url": source["fasta_url"].replace(".transcripts.fa.gz", ".annotation.gtf.gz"),
        "checksum_url": source["checksum_url"],
    }
    for species, source in GENCODE_SOURCES.items()
}
_FEATURES = {
    "exon": "exon", "CDS": "cds", "start_codon": "start", "stop_codon": "stop",
    "UTR": "utr", "five_prime_utr": "utr5", "three_prime_utr": "utr3",
    "5UTR": "utr5", "3UTR": "utr3",
}
_ATTRIBUTE = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)\s+"([^"\n]*)"\s*(?:;|$)')
_MISSING = (
    "Matching transcript-region annotations are missing. Run bash scripts/setup.sh --references "
    "to add them; existing FASTA and database downloads are reused. Custom references need a matching --gtf. "
    "Whole-transcript searches can run without region annotations."
)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Union overlapping/adjacent bases without joining a real sequence gap."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _subtract_intervals(
    intervals: list[tuple[int, int]], excluded: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Remove only explicitly annotated bases, preserving all remaining gaps."""
    remaining = list(intervals)
    for excluded_start, excluded_end in excluded:
        next_remaining: list[tuple[int, int]] = []
        for start, end in remaining:
            if excluded_end < start or excluded_start > end:
                next_remaining.append((start, end))
                continue
            if start < excluded_start:
                next_remaining.append((start, excluded_start - 1))
            if excluded_end < end:
                next_remaining.append((excluded_end + 1, end))
        remaining = next_remaining
    return remaining


def map_transcript_regions(
    *, transcript_id: str, length: int, strand: str,
    exons: list[tuple[int, int]], features: dict[str, list[tuple[int, int]]],
) -> dict[str, list[tuple[int, int]]]:
    """Map exonic genomic features onto a 1-based, spliced 5′→3′ transcript."""
    exons = sorted(exons)
    if not exons or strand not in ("+", "-"):
        raise ReferenceError(f"Missing exons or invalid strand for {transcript_id}.")
    previous_end = 0
    for start, end in exons:
        if start < 1 or end < start or start <= previous_end:
            raise ReferenceError(f"Overlapping or invalid exons for {transcript_id}.")
        previous_end = end
    if sum(end - start + 1 for start, end in exons) != length:
        raise ReferenceError(f"GTF exon length disagrees with stored transcript length for {transcript_id}.")
    oriented = exons if strand == "+" else list(reversed(exons))
    offsets: list[tuple[int, int, int]] = []
    offset = 1
    for start, end in oriented:
        offsets.append((start, end, offset))
        offset += end - start + 1

    def mapped(feature: str) -> list[tuple[int, int]]:
        result: list[tuple[int, int]] = []
        for start, end in features.get(feature, []):
            if start < 1 or end < start:
                raise ReferenceError(f"Invalid {feature} feature for {transcript_id}.")
            represented = 0
            for exon_start, exon_end, transcript_start in offsets:
                left, right = max(start, exon_start), min(end, exon_end)
                if left > right:
                    continue
                represented += right - left + 1
                if strand == "+":
                    result.append((transcript_start + left - exon_start,
                                   transcript_start + right - exon_start))
                else:
                    result.append((transcript_start + exon_end - right,
                                   transcript_start + exon_end - left))
            if represented != end - start + 1:
                raise ReferenceError(f"Non-exonic {feature} feature for {transcript_id}.")
        return merge_intervals(result)

    cds = mapped("cds")
    start_codons, stop_codons = mapped("start"), mapped("stop")
    if (start_codons or stop_codons) and not cds:
        raise ReferenceError(f"Codon annotation without CDS for {transcript_id}.")
    cds = merge_intervals(cds + start_codons + stop_codons)
    # GENCODE GTF generic UTR features include the stop codon, while its
    # transcript FASTA CDS tags include that same codon in the coding region.
    # Assign those explicitly annotated bases to CDS alone. Do not clip an
    # unannotated or genuine CDS overlap: it must still fail the checks below.
    utr5, utr3, utr = (
        _subtract_intervals(mapped(feature), stop_codons)
        for feature in ("utr5", "utr3", "utr")
    )
    if (utr5 or utr3 or utr) and not cds:
        raise ReferenceError(f"UTR annotation without CDS for {transcript_id}.")
    if cds:
        first_coding, last_coding = cds[0][0], cds[-1][1]
        for start, end in utr:
            if end < first_coding:
                utr5.append((start, end))
            elif start > last_coding:
                utr3.append((start, end))
            else:
                raise ReferenceError(f"UTR overlaps or falls inside the coding span for {transcript_id}.")
        if any(end >= first_coding for _, end in utr5) or any(start <= last_coding for start, _ in utr3):
            raise ReferenceError(f"UTR strand/coordinates disagree with CDS for {transcript_id}.")
    return {"cds": cds, "utr5": merge_intervals(utr5), "utr3": merge_intervals(utr3)}


def _parse_gtf(path: Path, records: dict[str, tuple[str, int]]) -> tuple[dict, int]:
    annotations: dict[str, dict] = {}
    ignored_transcripts: set[str] = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as source:
        for line_number, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ReferenceError(f"GTF line {line_number} does not contain nine columns.")
            chromosome, _, feature, start, end, _, strand, _, raw_attributes = fields
            if feature not in _FEATURES:
                continue
            attribute_pairs = _ATTRIBUTE.findall(raw_attributes)
            for identifier_key in ("transcript_id", "gene_id"):
                values = {value for key, value in attribute_pairs if key == identifier_key}
                if len(values) > 1:
                    raise ReferenceError(f"Conflicting {identifier_key} attributes at GTF line {line_number}.")
            attributes = dict(attribute_pairs)
            transcript_id, gene_id = attributes.get("transcript_id"), attributes.get("gene_id")
            if not transcript_id or not gene_id:
                raise ReferenceError(f"GTF line {line_number} lacks transcript_id or gene_id.")
            if transcript_id not in records:
                ignored_transcripts.add(transcript_id)
                continue
            if gene_id != records[transcript_id][0]:
                raise ReferenceError(f"GTF gene assignment disagrees with reference for {transcript_id}.")
            try:
                interval = (int(start), int(end))
            except ValueError as exc:
                raise ReferenceError(f"Invalid GTF coordinates at line {line_number}.") from exc
            annotation = annotations.setdefault(transcript_id, {
                "chromosome": chromosome, "strand": strand, "features": {},
            })
            if annotation["chromosome"] != chromosome or annotation["strand"] != strand:
                raise ReferenceError(f"Conflicting chromosome or strand for {transcript_id}.")
            annotation["features"].setdefault(_FEATURES[feature], []).append(interval)
    missing = sorted(records.keys() - annotations.keys())
    if missing:
        raise ReferenceError(
            f"GTF lacks matching transcript annotations for {len(missing):,} reference records: "
            + ", ".join(missing[:5]) + ". Use the same annotation release as the FASTA."
        )
    return annotations, len(ignored_transcripts)


def _binding(bundle: ReferenceBundle) -> dict:
    return {
        "source_sha256": bundle.manifest["source_sha256"],
        "transcript_sqlite_sha256": bundle.manifest["files"][bundle.manifest["sqlite"]]["sha256"],
        "transcript_count": bundle.manifest["transcript_count"],
    }


def annotation_manifest(bundle: ReferenceBundle, *, _root: Path | None = None) -> dict:
    root = _root or bundle.root / "annotations"
    if not root.exists():
        return {}
    try:
        manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
        if manifest.get("format_version") != 1 or manifest.get("reference") != _binding(bundle):
            raise ReferenceError("Region annotations belong to a different transcript reference.")
        files = manifest.get("files")
        if manifest.get("sqlite") != "regions.sqlite3" or not isinstance(files, dict) or set(files) != {"regions.sqlite3"}:
            raise ReferenceError("Invalid region-annotation manifest.")
        metadata = files["regions.sqlite3"]
        if (
            not isinstance(metadata, dict)
            or not isinstance(metadata.get("size"), int)
            or metadata["size"] < 1
            or not isinstance(metadata.get("sha256"), str)
            or not re.fullmatch(r"[a-f0-9]{64}", metadata["sha256"])
            or not isinstance(manifest.get("source"), dict)
            or not isinstance(manifest.get("source_sha256"), str)
            or not re.fullmatch(r"[a-f0-9]{64}", manifest["source_sha256"])
        ):
            raise ReferenceError("Invalid region-annotation file metadata.")
        return manifest
    except ReferenceError:
        raise
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        raise ReferenceError(f"Cannot open transcript-region annotations: {exc}") from exc


def verify_annotations(bundle: ReferenceBundle, *, _root: Path | None = None) -> None:
    root = _root or bundle.root / "annotations"
    manifest = annotation_manifest(bundle, _root=root)
    if not manifest:
        raise ReferenceError(_MISSING)
    for name, metadata in manifest["files"].items():
        if name != "regions.sqlite3":
            raise ReferenceError("Region manifest contains an unsupported file path.")
        path = root / name
        if not path.is_file() or path.stat().st_size != metadata["size"] or _hash_file(path) != metadata["sha256"]:
            raise ReferenceError(f"Region-annotation integrity verification failed for {name}.")
    try:
        path = root / "regions.sqlite3"
        with closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)) as connection:
            count = connection.execute("SELECT COUNT(*) FROM transcript_regions").fetchone()[0]
        if count != manifest["reference"]["transcript_count"]:
            raise ReferenceError("Region-annotation transcript count disagrees with reference.")
    except sqlite3.Error as exc:
        raise ReferenceError(f"Cannot verify region-annotation database: {exc}") from exc


def regions_for_gene(bundle: ReferenceBundle, gene_id: str) -> dict[str, dict[str, list[tuple[int, int]]]]:
    manifest = annotation_manifest(bundle)
    if not manifest:
        raise ReferenceError(_MISSING)
    path = bundle.root / "annotations" / manifest["sqlite"]
    try:
        with closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)) as connection:
            rows = connection.execute(
                "SELECT transcript_id, cds, utr5, utr3 FROM transcript_regions WHERE gene_id = ? ORDER BY transcript_id",
                (gene_id,),
            )
            result = {
                transcript_id: {
                    name: [tuple(interval) for interval in json.loads(encoded)]
                    for name, encoded in zip(("cds", "utr5", "utr3"), (cds, utr5, utr3), strict=True)
                }
                for transcript_id, cds, utr5, utr3 in rows
            }
        expected = bundle.transcripts_for_gene(gene_id)
        if result.keys() != expected.keys():
            raise ReferenceError("Region annotations omit target-gene transcripts.")
        for transcript_id, regions in result.items():
            for intervals in regions.values():
                if any(not 1 <= start <= end <= len(expected[transcript_id]) for start, end in intervals):
                    raise ReferenceError(f"Invalid stored region coordinates for {transcript_id}.")
        return result
    except ReferenceError:
        raise
    except (sqlite3.Error, ValueError, TypeError) as exc:
        raise ReferenceError(f"Cannot read transcript-region annotations: {exc}") from exc


def prepare_reference_annotations(
    bundle: ReferenceBundle, *, source_gtf: str | Path | None = None,
    expected_sha256: str | None = None,
) -> dict:
    """Install/reuse checked GTF coordinates bound to an existing FASTA bundle."""
    if expected_sha256 and not re.fullmatch(r"[a-fA-F0-9]{64}", expected_sha256):
        raise ReferenceError("Expected GTF SHA-256 must contain exactly 64 hexadecimal characters.")
    local_path = Path(source_gtf).resolve() if source_gtf else None
    species = bundle.manifest["source"]["species"]
    if local_path:
        if not local_path.is_file():
            raise ReferenceError(f"Local annotation GTF does not exist: {local_path}")
        source = {
            "species": species, "release": bundle.manifest["source"]["release"],
            "assembly": bundle.manifest["source"]["assembly"], "local_filename": local_path.name,
        }
    else:
        if bundle.manifest["source"] != GENCODE_SOURCES[species]:
            raise ReferenceError("Custom reference regions require a matching local GTF via --gtf.")
        source = dict(ANNOTATION_SOURCES[species])
    # Region coordinates must bind to the checked sequence database we actually
    # read, rather than merely trusting an unchanged original manifest.
    if _hash_file(bundle.sqlite_path) != bundle.manifest["files"][bundle.manifest["sqlite"]]["sha256"]:
        raise ReferenceError("Transcript database integrity failed before preparing annotations.")
    if (bundle.root / "annotations").exists():
        existing = annotation_manifest(bundle)
        if existing["source"] != source:
            raise ReferenceError("Existing region annotations use a different source; preserved without replacement.")
        verify_annotations(bundle)
        if local_path and _hash_file(local_path) != existing["source_sha256"]:
            raise ReferenceError("Local GTF differs from existing region annotations; preserved without replacement.")
        if expected_sha256 and expected_sha256.lower() != existing["source_sha256"]:
            raise ReferenceError("Expected GTF SHA-256 does not match existing annotations.")
        return existing
    temporary = Path(tempfile.mkdtemp(prefix=".annotations-", dir=bundle.root))
    try:
        authoritative_md5 = None
        if local_path:
            input_path = local_path
        else:
            checksum_path = temporary / "MD5SUMS"
            _download(source["checksum_url"], checksum_path)
            filename = source["gtf_url"].rsplit("/", 1)[-1]
            authoritative_md5 = _official_md5(checksum_path, filename)
            input_path = temporary / filename
            _download(source["gtf_url"], input_path)
            if _hash_file(input_path, "md5") != authoritative_md5:
                raise ReferenceError("Downloaded GTF failed the authoritative MD5 checksum; annotations were not published.")
        source_sha256 = _hash_file(input_path)
        if expected_sha256 and source_sha256 != expected_sha256.lower():
            raise ReferenceError("GTF failed the supplied SHA-256 checksum; annotations were not published.")
        with bundle._connection() as reference:
            records = {transcript_id: (gene_id, length) for transcript_id, gene_id, length in reference.execute(
                "SELECT transcript_id, gene_id, LENGTH(sequence) FROM transcripts"
            )}
        annotations, ignored = _parse_gtf(input_path, records)
        database_path = temporary / "regions.sqlite3"
        counts = {"cds": 0, "utr5": 0, "utr3": 0}
        with closing(sqlite3.connect(database_path)) as connection, connection:
            connection.execute(
                "CREATE TABLE transcript_regions (transcript_id TEXT PRIMARY KEY, gene_id TEXT NOT NULL, chromosome TEXT NOT NULL, strand TEXT NOT NULL, transcript_length INTEGER NOT NULL, cds TEXT NOT NULL, utr5 TEXT NOT NULL, utr3 TEXT NOT NULL)"
            )
            connection.execute("CREATE INDEX region_gene ON transcript_regions(gene_id)")
            for transcript_id in sorted(records):
                gene_id, length = records[transcript_id]
                annotation = annotations.pop(transcript_id)
                features = annotation["features"]
                regions = map_transcript_regions(
                    transcript_id=transcript_id, length=length, strand=annotation["strand"],
                    exons=features.pop("exon", []), features=features,
                )
                for name in counts:
                    counts[name] += bool(regions[name])
                connection.execute("INSERT INTO transcript_regions VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (
                    transcript_id, gene_id, annotation["chromosome"], annotation["strand"], length,
                    *(json.dumps(regions[name], separators=(",", ":")) for name in ("cds", "utr5", "utr3")),
                ))
        manifest = {
            "format_version": 1, "source": source, "source_sha256": source_sha256,
            "authoritative_md5": authoritative_md5, "reference": _binding(bundle),
            "coordinate_system": "1-based inclusive, mature spliced transcript 5prime-to-3prime",
            "mapping_policy": "GTF exons define transcript coordinates; CDS includes annotated start/stop codons; explicitly annotated stop-codon bases are assigned to CDS and removed from GENCODE UTR features; only remaining annotated UTR bases are classified against the coding span; separated intervals remain separated",
            "transcripts_with_regions": counts, "matched_gtf_transcripts": len(records),
            "ignored_gtf_transcripts_absent_from_fasta": ignored,
            "sqlite": "regions.sqlite3", "files": {"regions.sqlite3": {
                "sha256": _hash_file(database_path), "size": database_path.stat().st_size,
            }},
        }
        (temporary / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        # Keep just the compact reproducible sidecar, not the large source archive.
        if not local_path:
            input_path.unlink()
            checksum_path.unlink()
        verify_annotations(bundle, _root=temporary)
        os.rename(temporary, bundle.root / "annotations")
        return manifest
    except ReferenceError:
        raise
    except (OSError, UnicodeError, sqlite3.Error) as exc:
        raise ReferenceError(f"Unable to prepare transcript-region annotations: {exc}") from exc
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
