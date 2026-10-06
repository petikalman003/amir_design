"""Guide selection and strict coverage of every annotated target isoform."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Mapping


DESIGN_RULES = "uploaded_design"
FIGURE_S1_RULES = "figure_s1_guide"


def normalize_cds(text: str) -> str:
    """Accept plain DNA or a single FASTA entry; reject ambiguous model inputs."""
    text = text.strip()
    if text.startswith(">"):
        lines = text.splitlines()
        if sum(line.startswith(">") for line in lines) != 1:
            raise ValueError("Submit one CDS sequence at a time.")
        text = "".join(lines[1:])
    sequence = re.sub(r"\s+", "", text).upper()
    if not sequence:
        raise ValueError("Paste a target CDS sequence.")
    if set(sequence) - set("ACGT"):
        raise ValueError("Use DNA bases A, C, G and T; resolve ambiguous bases before design.")
    if len(sequence) < 22:
        raise ValueError("The target CDS must contain at least 22 bases.")
    return sequence


def reverse_complement(sequence: str) -> str:
    if set(sequence) - set("ACGT"):
        raise ValueError("Reverse complementation requires uppercase A/C/G/T.")
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def passes_design_rules(target: str, profile: str = DESIGN_RULES) -> bool:
    """Keep uploaded Design.py's biological filters with explicit strand semantics.

    Most Figure S1 constraints are equivalent under reverse complementation.
    Its literal red-guide 'no A at N3' is not equivalent to the uploaded
    target-sense exclusion. Never change that one rule without recording it.
    """
    if profile not in {DESIGN_RULES, FIGURE_S1_RULES}:
        raise ValueError(f"Unknown guide-selection profile: {profile}")
    if len(target) != 22 or set(target) - set("ACGT"):
        return False
    excluded_base = "A" if profile == DESIGN_RULES else "T"
    if target[0] not in "CG" or target[2] == excluded_base:
        return False
    if target[9] not in "AT" or target[21] not in "AT":
        return False
    if any(run in target for run in ("AAAAA", "TTTTT", "CCCC", "GGGG")):
        return False
    first_gc = sum(base in "GC" for base in target[:9])
    last_gc = sum(base in "GC" for base in target[9:])
    return first_gc / 9 > 0.5 and last_gc / 13 < 0.5


def positions(sequence: str, target: str) -> tuple[int, ...]:
    """All overlapping occurrences as 1-based starts in a transcript or CDS."""
    starts: list[int] = []
    offset = sequence.find(target)
    while offset != -1:
        starts.append(offset + 1)
        offset = sequence.find(target, offset + 1)
    return tuple(starts)


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    target_sequence: str
    guide_dna: str
    cds_starts: tuple[int, ...]
    isoform_sites: dict[str, tuple[int, ...]]

    @property
    def guide_rna(self) -> str:
        return self.guide_dna.replace("T", "U")

    @property
    def seed_rna(self) -> str:
        return self.guide_rna[1:8]


@dataclass(frozen=True)
class CandidateSelection:
    candidates: list[Candidate]
    sequence_filtered_count: int
    rejected_isoform_count: int
    target_isoforms: tuple[str, ...]
    profile: str


def select_candidates(
    cds: str, target_transcripts: Mapping[str, str], *, profile: str = DESIGN_RULES
) -> CandidateSelection:
    """Retain unique target windows present exactly in every target transcript."""
    cds = normalize_cds(cds)
    if not target_transcripts:
        raise ValueError("The reference contains no isoforms for the target gene.")
    if any(not sequence for sequence in target_transcripts.values()):
        raise ValueError("The reference contains an empty target transcript.")
    windows: dict[str, list[int]] = {}
    for start in range(len(cds) - 21):
        window = cds[start : start + 22]
        if passes_design_rules(window, profile):
            windows.setdefault(window, []).append(start + 1)
    candidates: list[Candidate] = []
    for window, starts in windows.items():
        isoform_sites = {
            transcript_id: positions(sequence, window)
            for transcript_id, sequence in target_transcripts.items()
        }
        if not all(isoform_sites.values()):
            continue
        candidates.append(
            Candidate(
                candidate_id=f"amiR_{len(candidates) + 1:04d}",
                target_sequence=window,
                guide_dna=reverse_complement(window),
                cds_starts=tuple(starts),
                isoform_sites=isoform_sites,
            )
        )
    return CandidateSelection(
        candidates=candidates,
        sequence_filtered_count=len(windows),
        rejected_isoform_count=len(windows) - len(candidates),
        target_isoforms=tuple(target_transcripts),
        profile=profile,
    )
