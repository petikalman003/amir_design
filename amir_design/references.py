"""Versioned transcript references shared by coverage, BLAST and seed screening.

The sequence database deliberately keeps every transcript in the supplied FASTA.
It never selects a canonical transcript or strips transcript/gene version suffixes.
"""

from __future__ import annotations

from contextlib import closing, contextmanager
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import tempfile
from typing import Iterator
from urllib.parse import urlparse
from urllib.request import HTTPRedirectHandler, build_opener


class ReferenceError(RuntimeError):
    """A reference could not be prepared, opened, or used safely."""


class AmbiguousGeneError(ReferenceError):
    """The input CDS matches more than one annotated gene."""


GENCODE_SOURCES = {
    "human": {
        "species": "human",
        "release": "v47",
        "assembly": "GRCh38",
        "fasta_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.transcripts.fa.gz",
        "checksum_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/MD5SUMS",
    },
    "mouse": {
        "species": "mouse",
        "release": "M36",
        "assembly": "GRCm39",
        "fasta_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz",
        "checksum_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/MD5SUMS",
    },
}

_DNA = re.compile(r"^[ACGTRYSWKMBDHVN]+$")
_ID = re.compile(r"^[A-Za-z0-9_.:-]+$")
_FORMAT_VERSION = 1


def _hash_file(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_transcript_header(header: str) -> dict[str, str]:
    """Read official GENCODE or legacy ``gene|transcript|symbol`` headers.

    GENCODE headers contain transcript and gene IDs first, gene symbol at field
    six, and declared sequence length at field seven. Full-transcript headers
    give a biotype next; protein-coding FASTA headers instead give CDS/UTR tags.
    Version suffixes, including GENCODE's PAR_Y suffix, remain intact.
    """
    fields = header.lstrip(">").strip().split("|")
    if (
        len(fields) >= 8
        and fields[0].startswith(("ENST", "ENSMUST"))
        and fields[1].startswith(("ENSG", "ENSMUSG"))
    ):
        transcript_id, gene_id = fields[:2]
        gene_name, length, biotype = fields[5:8]
        try:
            declared_length = str(int(length))
        except ValueError as exc:
            raise ReferenceError(
                f"Invalid declared length in FASTA header for {transcript_id}."
            ) from exc
        if biotype.startswith(("CDS:", "UTR5:", "UTR3:")):
            regions: set[str] = set()
            for tag in filter(None, fields[7:]):
                match = re.fullmatch(r"(CDS|UTR5|UTR3):(\d+)-(\d+)", tag)
                if match is None:
                    raise ReferenceError(
                        f"Invalid GENCODE protein-coding region tag for {transcript_id}: {tag}."
                    )
                region, start, end = match.groups()
                if region in regions or not 1 <= int(start) <= int(end) <= int(length):
                    raise ReferenceError(
                        f"Invalid or duplicate GENCODE sequence region for {transcript_id}: {tag}."
                    )
                regions.add(region)
            if "CDS" not in regions:
                raise ReferenceError(
                    f"Protein-coding header lacks a CDS region for {transcript_id}."
                )
            biotype = "protein_coding"
        elif ":" in biotype:
            raise ReferenceError(f"Unrecognized GENCODE biotype or region tag for {transcript_id}.")
    elif len(fields) == 3 and all(fields):
        gene_id, transcript_id, gene_name = fields
        biotype, declared_length = "unspecified", ""
    else:
        raise ReferenceError(
            "Transcript FASTA headers must be GENCODE headers or "
            "gene_id|transcript_id|gene_symbol; plain transcript IDs lack gene annotation."
        )
    if not _ID.fullmatch(transcript_id) or not _ID.fullmatch(gene_id):
        raise ReferenceError(f"Unsupported transcript or gene identifier in header: {header[:160]}")
    if not gene_name or not biotype:
        raise ReferenceError(f"Missing gene symbol or transcript biotype for {transcript_id}.")
    return {
        "transcript_id": transcript_id,
        "gene_id": gene_id,
        "gene_name": gene_name,
        "biotype": biotype,
        "declared_length": declared_length,
    }


def _fasta_records(path: Path) -> Iterator[tuple[dict[str, str], str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as source:
        header: str | None = None
        pieces: list[str] = []
        for line_number, line in enumerate(source, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield _validated_record(header, pieces)
                header, pieces = line, []
            elif header is None:
                raise ReferenceError(f"Sequence before first FASTA header at line {line_number}.")
            else:
                pieces.append("".join(line.split()))
        if header is not None:
            yield _validated_record(header, pieces)


def _validated_record(header: str, pieces: list[str]) -> tuple[dict[str, str], str]:
    record = parse_transcript_header(header)
    sequence = "".join(pieces).upper().replace("U", "T")
    if not sequence or not _DNA.fullmatch(sequence):
        raise ReferenceError(f"Empty or invalid nucleotide sequence for {record['transcript_id']}.")
    if record["declared_length"] and int(record["declared_length"]) != len(sequence):
        raise ReferenceError(
            f"FASTA sequence length disagrees with header for {record['transcript_id']}."
        )
    return record, sequence


class _HTTPSOnlyRedirect(HTTPRedirectHandler):
    def redirect_request(self, request, response, code, message, headers, new_url):
        if urlparse(new_url).scheme != "https":
            raise ReferenceError("Reference download redirected away from HTTPS.")
        return super().redirect_request(request, response, code, message, headers, new_url)


def _download(url: str, destination: Path) -> None:
    if urlparse(url).scheme != "https":
        raise ReferenceError("Reference downloads require HTTPS with certificate verification.")
    try:
        opener = build_opener(_HTTPSOnlyRedirect())
        with opener.open(url, timeout=60) as response, destination.open("wb") as target:
            if urlparse(response.geturl()).scheme != "https":
                raise ReferenceError("Reference download redirected away from HTTPS.")
            shutil.copyfileobj(response, target, length=1024 * 1024)
    except ReferenceError:
        raise
    except Exception as exc:
        raise ReferenceError(f"Unable to download {url}: {exc}") from exc


def _official_md5(checksum_path: Path, filename: str) -> str:
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) == 2 and parts[1].lstrip("*").removeprefix("./") == filename:
            if re.fullmatch(r"[a-fA-F0-9]{32}", parts[0]):
                return parts[0].lower()
    raise ReferenceError(f"Authoritative MD5SUMS does not list {filename}.")


def _run(command: list[str]) -> str:
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
    except FileNotFoundError as exc:
        raise ReferenceError(
            f"Required program is missing: {command[0]}. Install NCBI BLAST+."
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise ReferenceError(
            f"{command[0]} failed: {(exc.stderr or exc.stdout or '').strip()[-3000:]}"
        ) from exc
    return result.stdout.strip()


def prepare_reference_bundle(
    root: str | Path,
    *,
    species: str,
    source_fasta: str | Path | None = None,
    expected_sha256: str | None = None,
    release: str | None = None,
    assembly: str | None = None,
) -> "ReferenceBundle":
    """Build an atomic bundle from a checked official download or local FASTA.

    Existing bundles are verified and reused only when they match the requested
    source. A different reference belongs in a different root; existing files are
    never silently replaced. Local inputs receive a recorded SHA-256; callers can
    additionally supply the checksum they obtained from their source.
    """
    if species not in GENCODE_SOURCES:
        raise ReferenceError("Species must be human or mouse.")
    if expected_sha256 and not re.fullmatch(r"[a-fA-F0-9]{64}", expected_sha256):
        raise ReferenceError("Expected SHA-256 must contain exactly 64 hexadecimal characters.")
    root = Path(root).resolve()
    local_path = Path(source_fasta).resolve() if source_fasta else None
    if local_path and not local_path.is_file():
        raise ReferenceError(f"Local reference FASTA does not exist: {local_path}")
    source = dict(GENCODE_SOURCES[species])
    if local_path:
        source = {
            "species": species,
            "release": release or "custom",
            "assembly": assembly or "custom",
            "local_filename": local_path.name,
        }
    if root.exists():
        bundle = ReferenceBundle.open(root)
        if bundle.manifest["source"] != source:
            raise ReferenceError(
                "Existing bundle uses a different reference. Choose a different bundle root."
            )
        bundle.verify_integrity()
        if local_path and _hash_file(local_path) != bundle.manifest["source_sha256"]:
            raise ReferenceError(
                "Local FASTA differs from the existing bundle. Choose a different bundle root."
            )
        if expected_sha256 and expected_sha256.lower() != bundle.manifest["source_sha256"]:
            raise ReferenceError("Expected SHA-256 does not match the existing bundle source.")
        return bundle
    root.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{root.name}-", dir=root.parent))
    try:
        authoritative_md5 = None
        if local_path:
            input_path = local_path
        else:
            checksum_path = temporary / "MD5SUMS"
            _download(source["checksum_url"], checksum_path)
            filename = Path(urlparse(source["fasta_url"]).path).name
            authoritative_md5 = _official_md5(checksum_path, filename)
            input_path = temporary / filename
            _download(source["fasta_url"], input_path)
            if _hash_file(input_path, "md5") != authoritative_md5:
                raise ReferenceError(
                    "Downloaded reference failed the authoritative MD5 checksum. Bundle was not published."
                )
        source_sha256 = _hash_file(input_path)
        if expected_sha256 and source_sha256 != expected_sha256.lower():
            raise ReferenceError(
                "Reference FASTA failed the supplied SHA-256 checksum. Bundle was not published."
            )

        database_path = temporary / "transcripts.sqlite3"
        fasta_path = temporary / "transcripts.fa"
        biotype_counts: dict[str, int] = {}
        count = 0
        with (
            closing(sqlite3.connect(database_path)) as connection,
            connection,
            fasta_path.open("w", encoding="ascii") as fasta,
        ):
            connection.execute(
                "CREATE TABLE transcripts (transcript_id TEXT PRIMARY KEY, gene_id TEXT NOT NULL, gene_name TEXT NOT NULL, biotype TEXT NOT NULL, sequence TEXT NOT NULL)"
            )
            connection.execute("CREATE INDEX transcript_gene ON transcripts(gene_id)")
            connection.execute(
                "CREATE TABLE genes (gene_id TEXT PRIMARY KEY, gene_name TEXT NOT NULL, biotype TEXT NOT NULL)"
            )
            for record, sequence in _fasta_records(input_path):
                existing_gene = connection.execute(
                    "SELECT gene_name, biotype FROM genes WHERE gene_id = ?", (record["gene_id"],)
                ).fetchone()
                if existing_gene and existing_gene[0] != record["gene_name"]:
                    raise ReferenceError(f"Conflicting gene symbols for {record['gene_id']}.")
                if existing_gene and existing_gene[1] != record["biotype"]:
                    connection.execute(
                        "UPDATE genes SET biotype = 'mixed' WHERE gene_id = ?", (record["gene_id"],)
                    )
                elif not existing_gene:
                    connection.execute(
                        "INSERT INTO genes VALUES (?, ?, ?)",
                        (record["gene_id"], record["gene_name"], record["biotype"]),
                    )
                try:
                    connection.execute(
                        "INSERT INTO transcripts VALUES (?, ?, ?, ?, ?)",
                        (
                            record["transcript_id"],
                            record["gene_id"],
                            record["gene_name"],
                            record["biotype"],
                            sequence,
                        ),
                    )
                except sqlite3.IntegrityError as exc:
                    raise ReferenceError(
                        f"Duplicate transcript ID: {record['transcript_id']}."
                    ) from exc
                fasta.write(f">{record['transcript_id']}\n")
                for offset in range(0, len(sequence), 80):
                    fasta.write(sequence[offset : offset + 80] + "\n")
                count += 1
                biotype_counts[record["biotype"]] = biotype_counts.get(record["biotype"], 0) + 1
            if not count:
                raise ReferenceError("Reference FASTA contains no transcripts.")
            gene_count = connection.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
        blast_dir = temporary / "blast"
        blast_dir.mkdir()
        blast_prefix = blast_dir / "transcripts"
        _run(
            [
                "makeblastdb",
                "-in",
                str(fasta_path),
                "-dbtype",
                "nucl",
                "-parse_seqids",
                "-out",
                str(blast_prefix),
            ]
        )
        blast_version = _run(["makeblastdb", "-version"]).splitlines()[0]
        stored_paths = [fasta_path, database_path, *sorted(blast_dir.iterdir())]
        manifest = {
            "format_version": _FORMAT_VERSION,
            "source": source,
            "source_sha256": source_sha256,
            "authoritative_md5": authoritative_md5,
            "transcript_policy": "all transcript records in source FASTA; no canonical filtering; versioned identifiers preserved",
            "transcript_count": count,
            "gene_count": gene_count,
            "transcript_biotype_counts": biotype_counts,
            "fasta": "transcripts.fa",
            "sqlite": "transcripts.sqlite3",
            "blast_db": "blast/transcripts",
            "blast_version": blast_version,
            "files": {
                str(path.relative_to(temporary)): {
                    "sha256": _hash_file(path),
                    "size": path.stat().st_size,
                }
                for path in stored_paths
            },
        }
        (temporary / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        # Check the complete staged result before making it visible to the UI.
        ReferenceBundle.open(temporary).verify_integrity()
        os.rename(temporary, root)
        return ReferenceBundle.open(root)
    except ReferenceError:
        raise
    except (OSError, UnicodeError, sqlite3.Error) as exc:
        raise ReferenceError(f"Unable to prepare reference bundle: {exc}") from exc
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


class ReferenceBundle:
    def __init__(self, root: Path, manifest: dict):
        self.root = root
        self.manifest = manifest
        self.fasta_path = self._path(manifest["fasta"])
        self.sqlite_path = self._path(manifest["sqlite"])
        self.blast_db_path = self._path(manifest["blast_db"])

    def _path(self, relative: str) -> Path:
        path = (self.root / relative).resolve()
        if not path.is_relative_to(self.root):
            raise ReferenceError("Reference manifest contains a path outside its bundle.")
        return path

    @classmethod
    def open(cls, root: str | Path) -> "ReferenceBundle":
        root = Path(root).resolve()
        try:
            manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
            if not isinstance(manifest, dict):
                raise ReferenceError("Reference manifest must contain an object.")
            if manifest.get("format_version") != _FORMAT_VERSION:
                raise ReferenceError("Unsupported reference manifest format.")
            bundle = cls(root, manifest)
            for key in ("fasta", "sqlite"):
                if manifest[key] not in manifest["files"]:
                    raise ReferenceError(
                        f"Reference manifest omits required file: {manifest[key]}."
                    )
            for filename, metadata in manifest["files"].items():
                path = bundle._path(filename)
                if not path.is_file() or path.stat().st_size != metadata["size"]:
                    raise ReferenceError(
                        f"Reference file is missing or has changed size: {filename}."
                    )
            for suffix in (".nhr", ".nin", ".nsq"):
                if not Path(str(bundle.blast_db_path) + suffix).is_file():
                    raise ReferenceError(
                        f"Missing BLAST nucleotide index: {bundle.blast_db_path}{suffix}"
                    )
            with bundle._connection() as connection:
                count = connection.execute("SELECT COUNT(*) FROM transcripts").fetchone()[0]
                if count != manifest["transcript_count"]:
                    raise ReferenceError("Transcript count disagrees with reference manifest.")
            return bundle
        except ReferenceError:
            raise
        except (OSError, ValueError, KeyError, TypeError, sqlite3.Error) as exc:
            raise ReferenceError(f"Cannot open reference bundle at {root}: {exc}") from exc

    @contextmanager
    def _connection(self) -> Iterator[sqlite3.Connection]:
        connection = sqlite3.connect(self.sqlite_path.as_uri() + "?mode=ro", uri=True)
        try:
            yield connection
        finally:
            connection.close()

    def verify_integrity(self) -> None:
        for filename, metadata in self.manifest["files"].items():
            path = self._path(filename)
            if not path.is_file() or _hash_file(path) != metadata["sha256"]:
                raise ReferenceError(f"Reference integrity verification failed for {filename}.")

    def transcripts_for_gene(self, gene_id: str) -> dict[str, str]:
        with self._connection() as connection:
            return dict(
                connection.execute(
                    "SELECT transcript_id, sequence FROM transcripts WHERE gene_id = ? ORDER BY transcript_id",
                    (gene_id,),
                )
            )

    def transcript_gene_map(self) -> dict[str, str]:
        with self._connection() as connection:
            return dict(
                connection.execute(
                    "SELECT transcript_id, gene_id FROM transcripts ORDER BY transcript_id"
                )
            )

    def gene_for_cds(self, cds: str) -> dict:
        """Require an exact full-length CDS match to a single annotated gene.

        Every full-CDS hit is considered. A highest-scoring homolog is never
        silently selected, and minus-strand matches cannot assign a target gene.
        """
        cds = "".join(cds.split()).upper()
        if len(cds) < 22 or not re.fullmatch(r"[ACGT]+", cds):
            raise ReferenceError(
                "Target CDS must contain at least 22 unambiguous DNA bases (A/C/G/T)."
            )
        command = [
            "blastn",
            "-query",
            "-",
            "-db",
            str(self.blast_db_path),
            "-task",
            "blastn-short" if len(cds) < 50 else "blastn",
            "-strand",
            "plus",
            "-dust",
            "no",
            "-perc_identity",
            "100",
            "-max_target_seqs",
            "1000000",
            "-outfmt",
            "6 sseqid qstart qend length mismatch gapopen bitscore",
        ]
        try:
            result = subprocess.run(
                command, input=f">target_cds\n{cds}\n", capture_output=True, text=True, check=True
            )
        except FileNotFoundError as exc:
            raise ReferenceError("blastn is missing. Install NCBI BLAST+.") from exc
        except subprocess.CalledProcessError as exc:
            raise ReferenceError(
                f"CDS mapping failed: {(exc.stderr or '').strip()[-3000:]}"
            ) from exc
        transcript_ids: set[str] = set()
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) != 7:
                raise ReferenceError("Unexpected BLAST output while identifying the target gene.")
            try:
                start, end, length, mismatches, gaps = map(int, fields[1:6])
            except ValueError as exc:
                raise ReferenceError(
                    "Invalid BLAST coordinates while identifying the target gene."
                ) from exc
            if (
                start == 1
                and end == len(cds)
                and length == len(cds)
                and mismatches == 0
                and gaps == 0
            ):
                transcript_ids.add(fields[0].removeprefix("lcl|"))
        if not transcript_ids:
            raise ReferenceError(
                "No exact full-length CDS match in this reference. Check species, reference release, sequence orientation, and whether the input is a CDS rather than a whole transcript."
            )
        matches: dict[str, tuple[str, str]] = {}
        with self._connection() as connection:
            for transcript_id in transcript_ids:
                row = connection.execute(
                    "SELECT gene_id, gene_name FROM transcripts WHERE transcript_id = ?",
                    (transcript_id,),
                ).fetchone()
                if row is None:
                    raise ReferenceError(
                        f"BLAST returned an unknown transcript ID: {transcript_id}."
                    )
                matches[row[0]] = (row[0], row[1])
            if len(matches) != 1:
                names = ", ".join(
                    f"{name} ({gene_id})" for gene_id, name in sorted(matches.values())
                )
                raise AmbiguousGeneError(
                    f"Input CDS exactly matches multiple genes: {names}. Supply a gene-specific CDS."
                )
            gene_id = next(iter(matches))
            gene_name, biotype = connection.execute(
                "SELECT gene_name, biotype FROM genes WHERE gene_id = ?", (gene_id,)
            ).fetchone()
        return {
            "gene_id": gene_id,
            "gene_name": gene_name,
            "biotype": biotype,
            "transcript_ids": sorted(self.transcripts_for_gene(gene_id)),
            "matching_transcript_ids": sorted(transcript_ids),
        }
