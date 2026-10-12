"""Run the pinned SeedMatchR sequence API against a species transcript bundle.

The upstream functions are sourced, rather than replacing their search with a
Python implementation or installing the unrelated RNA-seq reporting stack.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
from typing import Callable, Iterator, Mapping, Sequence


SEEDMATCHR_COMMIT = "2fc0c4fba31e7306377488ac8b4c81d6445df5e2"
DEFAULT_SEED_CLASSES = ("mer8", "mer7m8", "mer7A1")
SEED_PRIORITY = {"mer8": 0, "mer7m8": 1, "mer7A1": 2, "mer6": 3}


@dataclass(frozen=True)
class OffTargetResult:
    summary: dict[str, dict[str, int]]
    sites: list[dict[str, str | int]] | OffTargetSiteStore
    provenance: dict[str, str]


SITE_COLUMNS = (
    "candidate_id", "transcript_id", "start", "end", "width", "g2_position",
    "seed_class", "target_sequence", "gene_id",
)
_INTEGER_SITE_COLUMNS = frozenset(("start", "end", "width", "g2_position"))


class OffTargetSiteStore:
    """Complete site rows held on disk, with bounded-memory preview and export.

    The temporary directory remains alive for as long as this result exists.
    Iteration opens the TSV afresh, so sorting summaries or exporting again does
    not consume the store. No search result is truncated.
    """

    columns = SITE_COLUMNS

    def __init__(
        self, path: Path, row_count: int, *, owner: tempfile.TemporaryDirectory | None = None,
    ):
        self.path = Path(path)
        self.row_count = row_count
        self._owner = owner

    def __len__(self) -> int:
        return self.row_count

    @property
    def empty(self) -> bool:
        return self.row_count == 0

    def iter_dicts(self) -> Iterator[dict[str, str | int]]:
        with self.path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != self.columns:
                raise RuntimeError("Off-target site file has unexpected columns")
            for raw in reader:
                yield {field: int(value) if field in _INTEGER_SITE_COLUMNS else value
                       for field, value in raw.items()}

    def iter_rows(self) -> Iterator[tuple[str | int, ...]]:
        for row in self.iter_dicts():
            yield tuple(row[column] for column in self.columns)

    def __iter__(self) -> Iterator[dict[str, str | int]]:
        return self.iter_dicts()

    def __getitem__(self, index: int) -> dict[str, str | int]:
        if index < 0:
            index += self.row_count
        if not 0 <= index < self.row_count:
            raise IndexError(index)
        for offset, row in enumerate(self.iter_dicts()):
            if offset == index:
                return row
        raise RuntimeError("Off-target site file is shorter than its recorded count")

    def __eq__(self, other: object) -> bool:
        if isinstance(other, list):
            return len(other) == self.row_count and all(
                row == expected for row, expected in zip(self, other, strict=True)
            )
        return self is other

    def to_frame(self, limit: int | None = None):
        import itertools
        import pandas as pd

        if limit is not None and limit < 0:
            raise ValueError("Preview limit must be nonnegative")
        rows = self.iter_rows()
        if limit is not None:
            rows = itertools.islice(rows, limit)
        return pd.DataFrame.from_records(rows, columns=self.columns)

    def head(self, n: int = 5):
        return self.to_frame(limit=n)


def _identifier(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or re.search(r"[\t\r\n\x00]", value):
        raise ValueError(f"{name} must be a nonempty identifier without control characters")
    return value


def summarize_sites(
    candidate_ids: Sequence[str],
    sites: Sequence[Mapping[str, str | int]],
    transcript_gene_map: Mapping[str, str],
    target_gene_id: str,
    *,
    seed_classes: Sequence[str] = DEFAULT_SEED_CLASSES,
) -> tuple[dict[str, dict[str, int]], list[dict[str, str | int]]]:
    """Exclude the intended gene and retain one strongest class per seed site.

    The position opposite guide nucleotide 2 identifies one binding location
    across the overlapping canonical seed definitions. Sites in two transcript
    isoforms remain two transcript sites; the gene count counts that gene once.
    """
    selected_classes = tuple(dict.fromkeys(seed_classes))
    if not selected_classes or any(seed not in SEED_PRIORITY for seed in selected_classes):
        raise ValueError("Seed classes must be selected from mer8, mer7m8, mer7A1, mer6")
    candidates = set(candidate_ids)
    strongest: dict[tuple[str, str, int], dict[str, str | int]] = {}
    for raw in sites:
        candidate = str(raw["candidate_id"])
        transcript = str(raw["transcript_id"])
        seed_class = str(raw["seed_class"])
        if candidate not in candidates:
            raise ValueError(f"SeedMatchR returned an unknown candidate: {candidate}")
        if transcript not in transcript_gene_map:
            raise ValueError(f"Transcript is missing its gene mapping: {transcript}")
        if seed_class not in SEED_PRIORITY:
            raise ValueError(f"Unsupported seed class: {seed_class}")
        if seed_class not in selected_classes:
            raise ValueError(f"SeedMatchR returned an unsearched seed class: {seed_class}")
        gene = transcript_gene_map[transcript]
        if gene == target_gene_id:
            continue
        row = dict(raw)
        row["gene_id"] = gene
        for field in ("start", "end", "width", "g2_position"):
            row[field] = int(row[field])
        key = candidate, transcript, int(row["g2_position"])
        previous = strongest.get(key)
        if previous is None or SEED_PRIORITY[seed_class] < SEED_PRIORITY[str(previous["seed_class"])]:
            strongest[key] = row
    result = sorted(
        strongest.values(),
        key=lambda row: (str(row["candidate_id"]), str(row["transcript_id"]), int(row["g2_position"])),
    )
    grouped: dict[str, list[dict[str, str | int]]] = {candidate: [] for candidate in candidate_ids}
    for row in result:
        grouped[str(row["candidate_id"])].append(row)
    summary = {}
    for candidate, rows in grouped.items():
        summary[candidate] = {
            "off_target_genes": len({str(row["gene_id"]) for row in rows}),
            "off_target_transcripts": len({str(row["transcript_id"]) for row in rows}),
            "off_target_sites": len(rows),
            **{
                f"{seed}_sites": sum(row["seed_class"] == seed for row in rows)
                for seed in selected_classes
            },
        }
    return summary, result


class SeedMatchRunner:
    def __init__(self, repo_path: Path, rscript: str = "Rscript"):
        self.repo_path = Path(repo_path).resolve()
        self.rscript = str(rscript)
        self.bridge = Path(__file__).resolve().parents[1] / "scripts" / "seedmatch.R"

    def _check_runtime(self) -> None:
        if shutil.which(self.rscript) is None:
            raise RuntimeError(f"Rscript is unavailable: {self.rscript}. Run scripts/install_r.sh first.")
        if not self.bridge.is_file():
            raise RuntimeError(f"SeedMatchR bridge is missing: {self.bridge}")
        for name in ("get_seed.R", "SeedMatchR.R"):
            if not (self.repo_path / "R" / name).is_file():
                raise RuntimeError(f"SeedMatchR source is missing: {self.repo_path / 'R' / name}")
        try:
            commit = subprocess.run(
                ["git", "-C", str(self.repo_path), "rev-parse", "HEAD"],
                check=True, capture_output=True, text=True, timeout=30,
            ).stdout.strip()
            subprocess.run(
                ["git", "-C", str(self.repo_path), "diff", "--quiet", "HEAD", "--", "R/get_seed.R", "R/SeedMatchR.R"],
                check=True, capture_output=True, text=True, timeout=30,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise RuntimeError("Cannot verify the pinned, unmodified SeedMatchR source") from exc
        if commit != SEEDMATCHR_COMMIT:
            raise RuntimeError(f"SeedMatchR must be pinned to {SEEDMATCHR_COMMIT}; found {commit}")

    def screen(
        self,
        guides: Mapping[str, str],
        reference_fasta: Path,
        target_gene_id: str,
        transcript_gene_map: Mapping[str, str],
        *,
        seed_classes: Sequence[str] = DEFAULT_SEED_CLASSES,
        include_sites: bool = True,
        timeout: int = 3600,
        progress: Callable[[str], None] | None = None,
    ) -> OffTargetResult:
        """Search exact canonical guide seeds over full spliced transcripts.

        Guides are mature 22-nt strands in 5' to 3' orientation; DNA or RNA text
        is accepted. All transcripts of the intended gene are excluded. The
        returned site count is a transcript-site count, not a genomic-site or
        experimentally established off-target count. Counts-only mode still searches
        every requested seed against the complete reference; it only avoids writing
        the potentially millions of detailed site rows. Progress is emitted as
        each guide completes.
        """
        _identifier(target_gene_id, "Target gene")
        seed_classes = tuple(dict.fromkeys(seed_classes))
        if not seed_classes or any(seed not in SEED_PRIORITY for seed in seed_classes):
            raise ValueError("Seed classes must be selected from mer8, mer7m8, mer7A1, mer6")
        normalized = {}
        for candidate, guide in guides.items():
            _identifier(candidate, "Candidate")
            dna = str(guide).upper().replace("U", "T")
            if len(dna) != 22 or re.fullmatch(r"[ACGT]+", dna) is None:
                raise ValueError(f"Candidate {candidate} requires a mature, unambiguous 22-nt guide")
            normalized[candidate] = dna
        provenance = {
            "engine": "SeedMatchR",
            "commit": SEEDMATCHR_COMMIT,
            "interface": "upstream get_seed and SeedMatchR iranges functions",
            "seed_classes": ",".join(seed_classes),
            "matching": "exact; no indels; no wobbles",
            "background": "full spliced transcript sequences",
            "primary_metric": "unique non-target genes",
            "coordinates": "1-based inclusive positions within each full spliced transcript",
        }
        if not normalized:
            return OffTargetResult({}, [], provenance)
        reference_fasta = Path(reference_fasta).resolve()
        if not reference_fasta.is_file():
            raise ValueError(f"Reference FASTA does not exist: {reference_fasta}")
        fasta_ids = []
        with reference_fasta.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    words = line[1:].strip().split()
                    if not words:
                        raise ValueError("Reference FASTA contains an empty transcript identifier")
                    fasta_ids.append(words[0])
        if not fasta_ids or len(set(fasta_ids)) != len(fasta_ids):
            raise ValueError("Reference FASTA must contain unique, nonempty transcript identifiers")
        missing = set(fasta_ids).difference(transcript_gene_map)
        if missing:
            raise ValueError(f"Reference transcripts missing gene mappings: {', '.join(sorted(missing)[:5])}")
        for transcript in fasta_ids:
            _identifier(transcript, "Transcript")
            _identifier(transcript_gene_map[transcript], "Mapped gene")
        if target_gene_id not in {transcript_gene_map[tx] for tx in fasta_ids}:
            raise ValueError(f"Target gene is absent from the reference: {target_gene_id}")
        self._check_runtime()
        bridge_source = self.bridge.read_bytes()
        provenance.update({
            "get_seed_sha256": hashlib.sha256((self.repo_path / "R/get_seed.R").read_bytes()).hexdigest(),
            "seedmatchr_api_sha256": hashlib.sha256((self.repo_path / "R/SeedMatchR.R").read_bytes()).hexdigest(),
            "bridge_sha256": hashlib.sha256(bridge_source).hexdigest(),
        })
        owner = tempfile.TemporaryDirectory(prefix="amir-seedmatch-")
        temporary = Path(owner.name)
        candidates_file = temporary / "candidates.tsv"
        sites_file = temporary / "sites.tsv"
        info_file = temporary / "provenance.tsv"
        mapping_file = temporary / "mapping.tsv"
        summary_file = temporary / "summary.tsv"
        progress_file = temporary / "progress.txt"
        output_file = temporary / "runtime.log"
        bridge_file = temporary / "seedmatch.R"
        try:
            bridge_file.write_bytes(bridge_source)
            with candidates_file.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(("candidate_id", "guide"))
                writer.writerows(normalized.items())
            with mapping_file.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(("transcript_id", "gene_id"))
                writer.writerows((tx, transcript_gene_map[tx]) for tx in fasta_ids)
            command = [
                self.rscript, "--vanilla", str(bridge_file), str(self.repo_path),
                str(candidates_file), str(reference_fasta), str(sites_file),
                str(info_file), ",".join(seed_classes), str(mapping_file), target_gene_id,
                str(summary_file), str(progress_file), "sites" if include_sites else "summary",
            ]
            self._run_bridge(command, output_file, progress_file, timeout, progress)
            if not all(path.is_file() for path in (sites_file, info_file, summary_file)):
                raise RuntimeError("SeedMatchR did not produce its expected output files")
            summary = {}
            with summary_file.open(newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    candidate = row.pop("candidate_id")
                    if candidate not in normalized or candidate in summary:
                        raise RuntimeError("SeedMatchR returned an inconsistent candidate summary")
                    counts = {field: int(value) for field, value in row.items()}
                    expected = {"off_target_genes", "off_target_transcripts", "off_target_sites",
                                *(f"{seed}_sites" for seed in seed_classes)}
                    if set(counts) != expected or any(value < 0 for value in counts.values()):
                        raise RuntimeError("SeedMatchR returned invalid off-target counts")
                    summary[candidate] = counts
            if set(summary) != set(normalized):
                raise RuntimeError("SeedMatchR omitted a candidate summary")
            with info_file.open(newline="") as handle:
                provenance.update({row["key"]: row["value"]
                                   for row in csv.DictReader(handle, delimiter="\t")})
            provenance["site_details"] = "complete" if include_sites else "not materialized; exact counts retained"
            row_count = sum(counts["off_target_sites"] for counts in summary.values())
            sites = OffTargetSiteStore(sites_file, row_count, owner=owner)
            if not include_sites:
                owner.cleanup()
                return OffTargetResult(summary, [], provenance)
            return OffTargetResult(summary, sites, provenance)
        except BaseException:
            owner.cleanup()
            raise

    @staticmethod
    def _run_bridge(
        command: list[str], output_file: Path, progress_file: Path, timeout: int,
        progress: Callable[[str], None] | None,
    ) -> None:
        """Drain progress without buffering the process's output or blocking the UI."""
        progress_file.touch()
        started = time.monotonic()
        try:
            with output_file.open("w") as output, progress_file.open() as updates:
                with subprocess.Popen(command, stdout=output, stderr=output, text=True) as process:
                    try:
                        while process.poll() is None:
                            for message in updates:
                                if progress is not None:
                                    progress(message.rstrip("\n"))
                            if time.monotonic() - started > timeout:
                                raise RuntimeError(f"SeedMatchR exceeded its {timeout}-second time limit")
                            time.sleep(0.15)
                        for message in updates:
                            if progress is not None:
                                progress(message.rstrip("\n"))
                    except BaseException:
                        process.kill()
                        process.wait()
                        raise
            if process.returncode:
                detail = output_file.read_text(errors="replace").strip()[-4000:]
                raise RuntimeError(f"SeedMatchR failed: {detail}")
        except OSError as exc:
            raise RuntimeError(f"Could not start SeedMatchR: {exc}") from exc
