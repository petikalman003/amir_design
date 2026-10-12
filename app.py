"""Run with ``python -m streamlit run app.py`` from the repository root."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import tempfile
from time import monotonic

import streamlit as st
from xlsxwriter.exceptions import XlsxWriterException

from amir_design.candidates import DESIGN_RULES, FIGURE_S1_RULES
from amir_design.offtargets import SeedMatchRunner
from amir_design.pipeline import run_identifier_pipeline, sort_result, with_site_details
from amir_design.references import ReferenceBundle, ReferenceError
from amir_design.reports import workbook_bytes, write_workbook
from amir_design.scoring import ShRNAIScorer


PROJECT_ROOT = Path(__file__).resolve().parent
REGION_LABELS = {
    "CDS only": "cds",
    "3′ UTR only": "utr3",
    "5′ UTR only": "utr5",
    "Whole transcript": "whole",
}
SORT_LABELS = {
    "Predicted activity": "activity",
    "Isoform coverage": "coverage",
    "Lowest off-target gene count": "off_targets",
}


def elapsed_text(seconds: float) -> str:
    seconds = int(seconds)
    return f"{seconds // 60}m {seconds % 60}s" if seconds >= 60 else f"{seconds}s"


def discard_detailed_workbook() -> None:
    path = st.session_state.pop("detailed_workbook", None)
    if path:
        Path(path).unlink(missing_ok=True)
    st.session_state.pop("detailed_workbook_sort", None)


@st.cache_resource
def get_scorer(path: str) -> ShRNAIScorer:
    return ShRNAIScorer(Path(path))


@st.cache_resource
def get_reference(path: str) -> ReferenceBundle:
    reference = ReferenceBundle.open(Path(path))
    reference.verify_integrity()
    if reference.annotation_provenance:
        reference.verify_annotations()
    return reference


class CachedScorer:
    """Defer model prerequisites until at least one eligible candidate exists."""

    def __init__(self, path: str):
        self.path = path

    def score(self, targets: list[str]) -> list[float]:
        return get_scorer(self.path).score(targets)

    @property
    def provenance(self):
        return get_scorer(self.path).provenance


def rscript_path() -> str:
    configured = os.environ.get("AMIR_RSCRIPT")
    if configured:
        return configured
    installed = shutil.which("Rscript")
    if installed:
        return installed
    local = PROJECT_ROOT.parent / ".runtime/bin/Rscript"
    return str(local) if local.is_file() else "Rscript"


def main() -> None:
    st.set_page_config(page_title="amiR Design", page_icon="🧬", layout="wide")
    st.title("amiR Design")
    st.write("Design 22-base miR-E guides, predict activity, and screen full transcripts for seed sites.")
    reference_root = Path(os.environ.get("AMIR_REFERENCE_ROOT", PROJECT_ROOT / "data/references"))
    tools_root = Path(os.environ.get("AMIR_TOOLS_ROOT", PROJECT_ROOT / "tools"))

    with st.form("design_form"):
        species_label = st.radio("Reference organism", ["Human", "Mouse"], horizontal=True)
        identifier = st.text_input(
            "Target gene or transcript",
            placeholder="Enter a gene ID, transcript ID, or gene symbol",
            help=("A gene selects candidate sites from all its annotated isoforms. A transcript "
                  "selects sites from that transcript. Both report coverage across every "
                  "annotated transcript of the gene. Gene symbols must identify one gene "
                  "in the reference."),
        )
        region_label = st.selectbox("Candidate search region", list(REGION_LABELS))
        require_all = st.checkbox(
            "Require all isoforms",
            value=False,
            help=("If checked, keep only guides with an exact 22-base target site in every "
                  "annotated transcript of the target gene. Otherwise keep guides covering "
                  "at least one transcript."),
        )
        with st.expander("Guide-selection rules"):
            rule_label = st.selectbox(
                "Position 3 rule",
                ["Uploaded Design.py rules", "Figure S1 literal guide-strand rule"],
                help=("The upload excludes A at target-site position 3. A literal reading of "
                      "Figure S1 excludes A at guide position 20, equivalent to excluding T "
                      "at target-site position 3. Other sequence filters are the same."),
            )
        submitted = st.form_submit_button("Design and screen", type="primary")

    st.caption("The search region controls where candidate sites originate. Isoform coverage "
               "checks their exact 22-base sites anywhere in each mature transcript. "
               "Off-target screening searches exact 7- and 8-base seed sites across full transcripts.")

    if submitted:
        # A failed or changed submission must not leave a previous result visible.
        discard_detailed_workbook()
        st.session_state.pop("result", None)
        st.session_state.pop("workbook", None)
        st.session_state.pop("result_species", None)
        st.session_state.pop("result_context", None)
        st.session_state.pop("workbook_sort", None)
        species = species_label.lower()
        profile = DESIGN_RULES if rule_label.startswith("Uploaded") else FIGURE_S1_RULES
        try:
            if not identifier.strip():
                raise ValueError("Enter a target gene or transcript before designing guides.")
            scorer = CachedScorer(str(tools_root / "shrnai"))
            runner = SeedMatchRunner(tools_root / "seedmatchr", rscript=rscript_path())
            template_path = os.environ.get("AMIR_MIRE_TEMPLATE")
            with st.status("Preparing your candidates", expanded=True) as status:
                started = monotonic()
                progress_line = st.empty()

                def show_progress(message: str) -> None:
                    progress_line.write(f"{message} · {elapsed_text(monotonic() - started)} elapsed")

                show_progress("Loading the selected reference")
                reference = get_reference(str(reference_root / species))
                result = run_identifier_pipeline(
                    identifier, reference, scorer, runner,
                    region=REGION_LABELS[region_label], require_all_isoforms=require_all,
                    profile=profile,
                    progress=show_progress,
                    template_path=Path(template_path) if template_path else None,
                    include_site_details=False,
                )
                label = ("No eligible candidate" if result.status == "no_eligible_candidate"
                         else f"Analysis complete in {elapsed_text(monotonic() - started)}")
                status.update(label=label, state="complete", expanded=False)
            st.session_state.result = result
            st.session_state.result_species = species
            st.session_state.result_context = {
                "region_label": region_label,
                "require_all_isoforms": require_all,
            }
        except ReferenceError as exc:
            st.error(str(exc))
        except (ValueError, RuntimeError, FileNotFoundError) as exc:
            st.error(str(exc))

    result = st.session_state.get("result")
    if result is None:
        return
    gene_name = result.metadata["target_gene_name"]
    gene_id = result.metadata["target_gene_id"]
    species = st.session_state.get("result_species", "reference")
    st.subheader(f"{species.title()} · {gene_name} · {gene_id}")
    context = st.session_state.result_context
    isoform_count = len(result.metadata["target_isoforms"])
    policy = "every isoform required" if context["require_all_isoforms"] else "at least one isoform required"
    st.caption(f"Search region: {context['region_label']} · {isoform_count} annotated gene "
               f"transcripts in the coverage denominator · {policy}.")
    if result.metadata.get("input_identifier_kind") == "transcript":
        source_transcripts = ", ".join(result.metadata.get("source_transcript_ids", []))
        st.caption(f"Candidates generated from the selected transcript: {source_transcripts}.")
    else:
        st.caption("Candidates generated from the selected region across the gene's isoforms.")
    missing_regions = result.metadata.get("missing_region_transcript_ids", [])
    if missing_regions:
        st.info(f"{len(missing_regions)} source transcripts have no annotated selected region: "
                f"{', '.join(missing_regions)}. They remain in the isoform-coverage denominator.")
    st.caption("Isoform coverage is a transcript count, not a prediction of expression knockdown. "
               "The denominator includes noncoding and other annotated transcripts.")
    if result.status == "no_eligible_candidate":
        st.warning(result.notices[0])
        st.caption(f"{result.metadata['sequence_filtered_candidates']} unique candidates passed "
                   f"sequence filters; {result.metadata['rejected_incomplete_isoform_coverage']} "
                   "were rejected by the isoform-coverage requirement.")
    else:
        cols = st.columns(3)
        cols[0].metric("Eligible candidates", len(result.summary))
        cols[1].metric("Target isoforms", isoform_count)
        cols[2].metric("Lowest off-target gene count", int(result.summary.off_target_genes.min()))
    sort_label = st.selectbox("Sort results by", list(SORT_LABELS), key="result_sort")
    order = SORT_LABELS[sort_label]
    displayed_result = sort_result(result, order=order)
    if result.status != "no_eligible_candidate":
        st.dataframe(
            displayed_result.summary, hide_index=True, use_container_width=True,
            column_config={
                "predicted_efficiency_score": st.column_config.NumberColumn(
                    "Predicted activity", format="%.4f",
                ),
                "target_isoforms_covered": "Isoforms covered",
                "target_isoforms_total": "Total gene isoforms",
                "target_isoform_coverage_fraction": st.column_config.NumberColumn(
                    "Coverage fraction", format="%.3f",
                    help="Fraction of annotated transcripts covered, not expression knockdown.",
                ),
                "off_target_genes": "Off-target genes",
                "off_target_transcripts": "Off-target transcripts",
                "off_target_sites": "Off-target sites",
                "covered_transcript_ids": "Covered transcript IDs",
                "missed_transcript_ids": "Missed transcript IDs",
            },
        )
        st.caption("Scores are not knockdown percentages; seed counts indicate potential binding sites. "
                   "The summary download follows the selected sorting order.")
        if result.metadata.get("synthesis_status") == "template_required":
            st.warning("Mature guide sequences are available. Complete synthesis inserts require "
                       "a verified miR-E scaffold template.")
        with st.expander("Target isoform sites and candidate origins"):
            st.write("Exact target-site occurrences across the gene's transcripts")
            st.dataframe(result.on_target_sites, hide_index=True, use_container_width=True)
            st.write("Sites used to generate candidates in the selected search region")
            st.dataframe(result.candidate_origins, hide_index=True, use_container_width=True)
    if st.session_state.get("workbook_sort") != order:
        st.session_state.workbook = workbook_bytes(
            displayed_result, include_off_target_sites=False,
        )
        st.session_state.workbook_sort = order
    st.caption("The summary contains every candidate, exact off-target counts, isoform coverage, "
               "and synthesis sequences. Individual off-target positions are an optional export.")
    st.download_button(
        "Download results and sequences", data=st.session_state.workbook,
        file_name=f"amir_{species}_{gene_name}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
    if result.status != "no_eligible_candidate":
        site_count = int(result.summary.off_target_sites.sum())
        st.caption(f"The detailed workbook adds all {site_count:,} individual off-target "
                   "sites. Preparing it can take longer and create a large file.")
        if st.button("Prepare detailed workbook"):
            discard_detailed_workbook()
            descriptor, name = tempfile.mkstemp(prefix="amir-detailed-", suffix=".xlsx")
            os.close(descriptor)
            export_ready = False
            try:
                with st.status("Preparing detailed workbook", expanded=True) as report_status:
                    report_progress = st.empty()
                    report_started = monotonic()

                    def show_report_progress(message: str) -> None:
                        report_progress.write(
                            f"{message} · {elapsed_text(monotonic() - report_started)} elapsed"
                        )

                    reference = get_reference(str(reference_root / species))
                    runner = SeedMatchRunner(tools_root / "seedmatchr", rscript=rscript_path())
                    result = with_site_details(
                        result, reference, runner, progress=show_report_progress,
                    )
                    st.session_state.result = result
                    detailed_result = sort_result(result, order=order)
                    write_workbook(detailed_result, Path(name), progress=show_report_progress)
                    report_status.update(label="Detailed workbook ready", state="complete",
                                         expanded=False)
                st.session_state.detailed_workbook = name
                st.session_state.detailed_workbook_sort = sort_label
                export_ready = True
            except (ValueError, RuntimeError, OSError, XlsxWriterException) as exc:
                st.error(f"Could not prepare the detailed workbook: {exc}")
            finally:
                # Streamlit cancellation and rerun exceptions must also remove partial files.
                if not export_ready:
                    Path(name).unlink(missing_ok=True)
        detailed_path = st.session_state.get("detailed_workbook")
        if detailed_path:
            st.caption(f"Detailed workbook order: {st.session_state.detailed_workbook_sort}. "
                       "Changing the table sorting keeps this prepared file; prepare it again "
                       "to change its order.")
            with Path(detailed_path).open("rb") as report_file:
                st.download_button(
                    "Download detailed workbook", data=report_file,
                    file_name=f"amir_{species}_{gene_name}_detailed.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                )
        if result.metadata.get("off_target_site_details_retained"):
            with st.expander("Off-target site details"):
                st.dataframe(result.sites.head(10000), hide_index=True, use_container_width=True)
                if len(result.sites) > 10000:
                    st.caption("Showing the first 10,000 sites. The detailed workbook contains "
                               "every site.")


if __name__ == "__main__":
    main()
