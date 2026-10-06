"""Orchestrate Design -> shRNAI -> SeedMatchR without losing candidate identity."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from amir_design import __version__
from amir_design.candidates import DESIGN_RULES, normalize_cds, select_candidates


SUMMARY_COLUMNS = [
    "candidate_id", "predicted_efficiency_score", "off_target_genes",
    "off_target_transcripts", "off_target_sites", "target_isoforms_covered",
    "target_isoforms_total", "guide_rna_5to3", "guide_dna_5to3",
    "seed_rna_positions_2to8", "target_site_dna_5to3", "cds_start_positions_1based",
    "mire_insert_dna_5to3", "synthesis_status",
]


@dataclass
class PipelineResult:
    status: str
    summary: pd.DataFrame
    sites: pd.DataFrame
    on_target_sites: pd.DataFrame
    metadata: dict[str, Any]
    notices: list[str] = field(default_factory=list)


def run_pipeline(
    cds_text: str,
    reference: Any,
    scorer: Any,
    seed_runner: Any,
    *,
    profile: str = DESIGN_RULES,
    progress: Callable[[str], None] | None = None,
    insert_builder: Callable[[str], str] | None = None,
    scaffold_provenance: dict[str, Any] | None = None,
    template_path: Path | None = None,
) -> PipelineResult:
    """Run all stages; unavailable or failed stages raise rather than invent data.

    The scorer and R runner are injectable to enable isolated stage tests;
    production always supplies the actual upstream model and SeedMatchR adapter.
    """
    from amir_design.scaffold import (
        DEFAULT_MIRE_TEMPLATE, ScaffoldUnavailableError, build_mire_insert, load_template,
    )

    notify = progress or (lambda message: None)
    cds = normalize_cds(cds_text)
    notify("Identifying the target gene in the selected reference")
    gene = reference.gene_for_cds(cds)
    transcripts = reference.transcripts_for_gene(gene["gene_id"])
    notify(f"Checking shared target sites across {len(transcripts)} isoforms")
    selection = select_candidates(cds, transcripts, profile=profile)
    metadata: dict[str, Any] = {
        "application_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_cds_sha256": hashlib.sha256(cds.encode("ascii")).hexdigest(),
        "input_cds_length": len(cds),
        "species_reference": reference.manifest,
        "target_gene_id": gene["gene_id"],
        "target_gene_name": gene["gene_name"],
        "target_isoforms": list(selection.target_isoforms),
        "design_rule_profile": selection.profile,
        "sequence_filtered_candidates": selection.sequence_filtered_count,
        "rejected_incomplete_isoform_coverage": selection.rejected_isoform_count,
        "on_target_matching": "exact 22-nt target site in every annotated gene transcript",
        "off_target_scope": "full spliced reference transcripts; all target-gene isoforms excluded",
        "score_interpretation": "uncalibrated shRNAI activity score; higher is predicted stronger",
        "sort_order": "efficiency score descending, off-target genes ascending; no combined score",
        "eligible_candidates": len(selection.candidates),
    }
    if not selection.candidates:
        metadata["stages"] = {
            "candidate_selection": "completed",
            "shrnai": "not run: no eligible candidate",
            "seedmatchr": "not run: no eligible candidate",
            "synthesis": "not run: no eligible candidate",
        }
        return PipelineResult(
            status="no_eligible_candidate",
            summary=pd.DataFrame(columns=SUMMARY_COLUMNS),
            sites=pd.DataFrame(),
            on_target_sites=pd.DataFrame(),
            metadata=metadata,
            notices=["No eligible candidate covers every target isoform with the selected rules."],
        )

    candidates = selection.candidates
    if template_path is not None:
        template = load_template(template_path)

        def insert_builder(guide: str) -> str:
            return build_mire_insert(guide, template=template)

        scaffold_provenance = template.provenance
    elif scaffold_provenance is None and insert_builder is None and DEFAULT_MIRE_TEMPLATE is not None:
        scaffold_provenance = DEFAULT_MIRE_TEMPLATE.provenance
    notify(f"Predicting efficiency for {len(candidates)} eligible guides")
    scores = scorer.score([candidate.target_sequence for candidate in candidates])
    if len(scores) != len(candidates):
        raise RuntimeError("Efficiency scoring returned an inconsistent candidate count.")
    metadata["shrnai"] = scorer.provenance
    notify("Screening guide seeds across full transcripts")
    off_targets = seed_runner.screen(
        {candidate.candidate_id: candidate.guide_dna for candidate in candidates},
        reference.fasta_path,
        gene["gene_id"],
        reference.transcript_gene_map(),
    )
    metadata["seedmatchr"] = off_targets.provenance
    metadata["site_coordinates"] = "1-based inclusive positions within full spliced transcripts"
    if scaffold_provenance is not None:
        metadata["synthesis_scaffold"] = scaffold_provenance
    builder = insert_builder or build_mire_insert
    rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    synthesis_missing = False
    for candidate, score in zip(candidates, scores, strict=True):
        counts = off_targets.summary.get(candidate.candidate_id)
        if counts is None:
            raise RuntimeError(f"No off-target summary returned for {candidate.candidate_id}.")
        try:
            insert = builder(candidate.guide_dna)
            synthesis_status = "miR-E insert; cloning junctions are vector-specific"
        except ScaffoldUnavailableError:
            insert = ""
            synthesis_status = "Requires a verified complete miR-E template"
            synthesis_missing = True
        row = {
            "candidate_id": candidate.candidate_id,
            "predicted_efficiency_score": float(score),
            **counts,
            "target_isoforms_covered": len(candidate.isoform_sites),
            "target_isoforms_total": len(transcripts),
            "guide_rna_5to3": candidate.guide_rna,
            "guide_dna_5to3": candidate.guide_dna,
            "seed_rna_positions_2to8": candidate.seed_rna,
            "target_site_dna_5to3": candidate.target_sequence,
            "cds_start_positions_1based": ";".join(map(str, candidate.cds_starts)),
            "mire_insert_dna_5to3": insert,
            "synthesis_status": synthesis_status,
        }
        rows.append(row)
        for transcript_id, starts in candidate.isoform_sites.items():
            for start in starts:
                target_rows.append({
                    "candidate_id": candidate.candidate_id,
                    "transcript_id": transcript_id,
                    "gene_id": gene["gene_id"],
                    "start_1based": start,
                    "end_1based": start + 21,
                })
    summary = pd.DataFrame(rows).sort_values(
        ["predicted_efficiency_score", "off_target_genes", "candidate_id"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    metadata["eligible_candidates"] = len(candidates)
    metadata["synthesis_status"] = "template_required" if synthesis_missing else "insert_available"
    metadata["stages"] = {
        "candidate_selection": "completed", "shrnai": "completed", "seedmatchr": "completed",
        "synthesis": "template required" if synthesis_missing else "declared template assembled",
    }
    return PipelineResult(
        status="analysis_complete" if synthesis_missing else "complete",
        summary=summary,
        sites=pd.DataFrame(off_targets.sites),
        on_target_sites=pd.DataFrame(target_rows),
        metadata=metadata,
        notices=[
            "Efficiency scores are model predictions, not knockdown percentages.",
            "Seed sites indicate potential binding; counts do not establish off-target knockdown.",
            *(["Complete miR-E inserts require a verified template; mature guides are provided."]
              if synthesis_missing else []),
        ],
    )
