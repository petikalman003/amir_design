"""Orchestrate Design -> shRNAI -> SeedMatchR without losing candidate identity."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
import hashlib
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from amir_design import __version__
from amir_design.candidates import (
    DESIGN_RULES, CandidateSelection, normalize_cds, select_candidates, select_region_candidates,
)


SUMMARY_COLUMNS = [
    "candidate_id", "predicted_efficiency_score", "off_target_genes",
    "off_target_transcripts", "off_target_sites", "target_isoforms_covered",
    "target_isoforms_total", "target_isoform_coverage_fraction", "covered_transcript_ids",
    "missed_transcript_ids", "source_region", "source_transcript_ids",
    "source_transcript_start_positions_1based", "guide_rna_5to3", "guide_dna_5to3",
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
    candidate_origins: pd.DataFrame = field(default_factory=pd.DataFrame)


_SORT_ORDERS = {
    "activity": (
        ["predicted_efficiency_score", "off_target_genes", "target_isoform_coverage_fraction",
         "candidate_id"], [False, True, False, True],
        "activity descending, off-target genes ascending, coverage descending",
    ),
    "coverage": (
        ["target_isoform_coverage_fraction", "predicted_efficiency_score", "off_target_genes",
         "candidate_id"], [False, False, True, True],
        "isoform coverage descending, activity descending, off-target genes ascending",
    ),
    "off_targets": (
        ["off_target_genes", "predicted_efficiency_score", "target_isoform_coverage_fraction",
         "candidate_id"], [True, False, False, True],
        "off-target genes ascending, activity descending, coverage descending",
    ),
}


def sort_candidates(summary: pd.DataFrame, sort_by: str = "activity") -> pd.DataFrame:
    """Sort independent metrics without constructing a combined activity score."""
    if sort_by not in _SORT_ORDERS:
        raise ValueError(f"Unknown candidate sort: {sort_by}")
    if summary.empty:
        return summary.copy()
    columns, ascending, _ = _SORT_ORDERS[sort_by]
    return summary.sort_values(columns, ascending=ascending, kind="stable").reset_index(drop=True)


def sort_result(result: PipelineResult, order: str = "activity") -> PipelineResult:
    """Return a sorted report without repeating model prediction or seed search."""
    summary = sort_candidates(result.summary, order)
    metadata = {**result.metadata, "sort_by": order, "sort_order": _SORT_ORDERS[order][2]}
    return replace(result, summary=summary, metadata=metadata)


def _selection_metadata(reference: Any, gene: dict, selection: CandidateSelection) -> dict[str, Any]:
    return {
        "application_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "species_reference": reference.manifest,
        "target_gene_id": gene["gene_id"],
        "target_gene_name": gene["gene_name"],
        "target_isoforms": list(selection.target_isoforms),
        "design_rule_profile": selection.profile,
        "sequence_filtered_candidates": selection.sequence_filtered_count,
        "rejected_incomplete_isoform_coverage": selection.rejected_isoform_count,
        "require_all_isoforms": selection.require_all_isoforms,
        "on_target_matching": "exact 22-nt target site anywhere in each full gene transcript",
        "isoform_coverage_denominator": "all annotated transcripts of the target gene",
        "off_target_scope": "full spliced reference transcripts; all target-gene isoforms excluded",
        "score_interpretation": "uncalibrated shRNAI activity score; higher is predicted stronger",
        "eligible_candidates": len(selection.candidates),
    }


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
    require_all_isoforms: bool = True,
    sort_by: str = "activity",
) -> PipelineResult:
    """Keep the original pasted-CDS entry point and complete-coverage default."""
    notify = progress or (lambda message: None)
    cds = normalize_cds(cds_text)
    notify("Identifying the target gene in the selected reference")
    gene = reference.gene_for_cds(cds)
    transcripts = reference.transcripts_for_gene(gene["gene_id"])
    notify(f"Checking target sites across {len(transcripts)} isoforms")
    selection = select_candidates(
        cds, transcripts, profile=profile, require_all_isoforms=require_all_isoforms,
    )
    metadata = _selection_metadata(reference, gene, selection)
    metadata.update({
        "input_cds_sha256": hashlib.sha256(cds.encode("ascii")).hexdigest(),
        "input_cds_length": len(cds),
        "input_identifier_kind": "pasted_cds",
        "search_region": "cds",
    })
    return _run_selected_pipeline(
        selection, transcripts, gene, reference, scorer, seed_runner, metadata,
        notify=notify, insert_builder=insert_builder, scaffold_provenance=scaffold_provenance,
        template_path=template_path, sort_by=sort_by,
    )


def run_identifier_pipeline(
    identifier: str,
    reference: Any,
    scorer: Any,
    seed_runner: Any,
    *,
    region: str = "cds",
    require_all_isoforms: bool = False,
    sort_by: str = "activity",
    profile: str = DESIGN_RULES,
    progress: Callable[[str], None] | None = None,
    insert_builder: Callable[[str], str] | None = None,
    scaffold_provenance: dict[str, Any] | None = None,
    template_path: Path | None = None,
) -> PipelineResult:
    """Generate from a gene or transcript ID, retaining the full-gene denominator."""
    if region not in {"cds", "utr3", "utr5", "whole"}:
        raise ValueError(f"Unknown search region: {region}")
    notify = progress or (lambda message: None)
    notify("Finding the gene or transcript in the selected reference")
    gene = reference.resolve_identifier(identifier)
    transcripts = reference.transcripts_for_gene(gene["gene_id"])
    generation_ids = (
        gene["matching_transcript_ids"] if gene["input_type"] == "transcript"
        else list(transcripts)
    )
    if region == "whole":
        regions = {transcript: [(1, len(sequence))] for transcript, sequence in transcripts.items()}
    else:
        annotation = reference.regions_for_gene(gene["gene_id"])
        regions = {transcript: annotation.get(transcript, {}).get(region, [])
                   for transcript in transcripts}
    notify(f"Selecting {region} target sites and checking {len(transcripts)} gene isoforms")
    selection = select_region_candidates(
        transcripts, regions, generation_transcript_ids=generation_ids,
        require_all_isoforms=require_all_isoforms, profile=profile, region=region,
    )
    metadata = _selection_metadata(reference, gene, selection)
    metadata.update({
        "input_identifier": identifier.strip(),
        "input_identifier_kind": gene["input_type"],
        "search_region": region,
        "source_transcript_ids": list(selection.generation_transcript_ids),
        "missing_region_transcript_ids": list(selection.missing_region_transcript_ids),
        "gene_transcripts_without_selected_region": [transcript for transcript in transcripts
                                                    if not regions.get(transcript)],
        "ambiguous_source_windows_skipped": selection.ambiguous_windows_skipped,
        "candidate_region_policy": "each complete 22-nt source window lies within selected region",
        "annotation_provenance": getattr(reference, "annotation_provenance", {}),
    })
    return _run_selected_pipeline(
        selection, transcripts, gene, reference, scorer, seed_runner, metadata,
        notify=notify, insert_builder=insert_builder, scaffold_provenance=scaffold_provenance,
        template_path=template_path, sort_by=sort_by,
    )


def _run_selected_pipeline(
    selection: CandidateSelection,
    transcripts: dict[str, str],
    gene: dict,
    reference: Any,
    scorer: Any,
    seed_runner: Any,
    metadata: dict[str, Any],
    *,
    notify: Callable[[str], None],
    insert_builder: Callable[[str], str] | None,
    scaffold_provenance: dict[str, Any] | None,
    template_path: Path | None,
    sort_by: str,
) -> PipelineResult:
    """Run all stages; unavailable or failed stages raise rather than invent data.

    The scorer and R runner are injectable to enable isolated stage tests;
    production always supplies the actual upstream model and SeedMatchR adapter.
    """
    from amir_design.scaffold import (
        DEFAULT_MIRE_TEMPLATE, ScaffoldUnavailableError, build_mire_insert, load_template,
    )

    if sort_by not in _SORT_ORDERS:
        raise ValueError(f"Unknown candidate sort: {sort_by}")
    metadata.update({"sort_by": sort_by, "sort_order": _SORT_ORDERS[sort_by][2]})
    if not selection.candidates:
        metadata["stages"] = {
            "candidate_selection": "completed",
            "shrnai": "not run: no eligible candidate",
            "seedmatchr": "not run: no eligible candidate",
            "synthesis": "not run: no eligible candidate",
        }
        if (selection.generation_transcript_ids and len(selection.missing_region_transcript_ids)
                == len(selection.generation_transcript_ids)):
            notice = "The selected source transcripts have no annotated sequence in this region."
        elif selection.sequence_filtered_count == 0:
            notice = "No 22-nt candidate in the selected sequence region passes the design rules."
        elif selection.require_all_isoforms:
            notice = "No eligible candidate covers every target isoform with the selected rules."
        else:
            notice = "No eligible candidate has an exact site in any target-gene isoform."
        return PipelineResult(
            status="no_eligible_candidate",
            summary=pd.DataFrame(columns=SUMMARY_COLUMNS),
            sites=pd.DataFrame(),
            on_target_sites=pd.DataFrame(),
            metadata=metadata,
            notices=[notice],
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
    origin_rows: list[dict[str, Any]] = []
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
            "target_isoforms_covered": len(candidate.covered_isoforms),
            "target_isoforms_total": len(transcripts),
            "target_isoform_coverage_fraction": candidate.coverage_fraction,
            "covered_transcript_ids": ";".join(candidate.covered_isoforms),
            "missed_transcript_ids": ";".join(candidate.missed_isoforms),
            "source_region": candidate.source_region,
            "source_transcript_ids": ";".join(candidate.source_sites),
            "source_transcript_start_positions_1based": ";".join(
                f"{transcript}:{','.join(map(str, starts))}"
                for transcript, starts in candidate.source_sites.items()
            ),
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
        for transcript_id, starts in candidate.source_sites.items():
            for start in starts:
                origin_rows.append({
                    "candidate_id": candidate.candidate_id,
                    "transcript_id": transcript_id,
                    "gene_id": gene["gene_id"],
                    "region": candidate.source_region,
                    "source_region": candidate.source_region,
                    "target_site_dna_5to3": candidate.target_sequence,
                    "start_1based": start,
                    "end_1based": start + 21,
                })
    # Keep SeedMatchR's per-seed-class counters as well as the base report columns.
    summary = sort_candidates(pd.DataFrame(rows), sort_by)
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
        candidate_origins=pd.DataFrame(origin_rows),
        notices=[
            "Efficiency scores are model predictions, not knockdown percentages.",
            "Isoform coverage counts transcripts, not their expression or expected knockdown.",
            "Seed sites indicate potential binding; counts do not establish off-target knockdown.",
            *(["Some source transcripts lack the selected annotated region; their IDs are "
               "recorded in Run information and they remain in the coverage denominator."]
              if selection.missing_region_transcript_ids else []),
            *(["Complete miR-E inserts require a verified template; mature guides are provided."]
              if synthesis_missing else []),
        ],
    )
