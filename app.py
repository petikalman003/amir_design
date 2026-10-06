"""Run with ``streamlit run app.py`` from the repository root."""

from __future__ import annotations

import os
from pathlib import Path
import shutil

import streamlit as st

from amir_design.candidates import DESIGN_RULES, FIGURE_S1_RULES, normalize_cds
from amir_design.offtargets import SeedMatchRunner
from amir_design.pipeline import run_pipeline
from amir_design.references import ReferenceBundle, ReferenceError
from amir_design.reports import workbook_bytes
from amir_design.scoring import ShRNAIScorer


PROJECT_ROOT = Path(__file__).resolve().parent


@st.cache_resource
def get_scorer(path: str) -> ShRNAIScorer:
    return ShRNAIScorer(Path(path))


@st.cache_resource
def get_reference(path: str) -> ReferenceBundle:
    reference = ReferenceBundle.open(Path(path))
    reference.verify_integrity()
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
        cds = st.text_area("Target CDS", height=200, placeholder="Paste your DNA CDS sequence here")
        with st.expander("Guide-selection rules"):
            rule_label = st.selectbox(
                "Position 3 rule",
                ["Uploaded Design.py rules", "Figure S1 literal guide-strand rule"],
                help=("The upload excludes A at target-site position 3. A literal reading of "
                      "Figure S1 excludes A at guide position 20, equivalent to excluding T "
                      "at target-site position 3. Other sequence filters are the same."),
            )
        submitted = st.form_submit_button("Design and screen", type="primary")

    st.caption("Candidates must cover every annotated isoform of the target gene. "
               "Off-target screening uses exact 7- and 8-base seed matches across complete transcripts.")

    if submitted:
        # A failed or changed submission must not leave a previous result visible.
        st.session_state.pop("result", None)
        st.session_state.pop("workbook", None)
        st.session_state.pop("result_species", None)
        species = species_label.lower()
        profile = DESIGN_RULES if rule_label.startswith("Uploaded") else FIGURE_S1_RULES
        try:
            cds = normalize_cds(cds)
            reference = get_reference(str(reference_root / species))
            scorer = CachedScorer(str(tools_root / "shrnai"))
            runner = SeedMatchRunner(tools_root / "seedmatchr", rscript=rscript_path())
            template_path = os.environ.get("AMIR_MIRE_TEMPLATE")
            with st.status("Preparing your candidates", expanded=True) as status:
                result = run_pipeline(
                    cds, reference, scorer, runner, profile=profile,
                    progress=st.write,
                    template_path=Path(template_path) if template_path else None,
                )
                report = workbook_bytes(result)
                label = "No eligible candidate" if result.status == "no_eligible_candidate" else "Analysis complete"
                status.update(label=label, state="complete", expanded=False)
            st.session_state.result = result
            st.session_state.workbook = report
            st.session_state.result_species = species
        except ReferenceError as exc:
            st.error(str(exc))
            st.info("Prepare the selected species reference using the setup instructions in README.md.")
        except (ValueError, RuntimeError, FileNotFoundError) as exc:
            st.error(str(exc))

    result = st.session_state.get("result")
    if result is None:
        return
    gene_name = result.metadata["target_gene_name"]
    gene_id = result.metadata["target_gene_id"]
    species = st.session_state.get("result_species", "reference")
    st.subheader(f"{species.title()} · {gene_name} · {gene_id}")
    if result.status == "no_eligible_candidate":
        st.warning(result.notices[0])
        st.caption(f"{result.metadata['sequence_filtered_candidates']} unique candidates passed "
                   f"sequence filters; {result.metadata['rejected_incomplete_isoform_coverage']} "
                   "were rejected for incomplete isoform coverage.")
    else:
        cols = st.columns(3)
        cols[0].metric("Eligible candidates", len(result.summary))
        cols[1].metric("Target isoforms", len(result.metadata["target_isoforms"]))
        cols[2].metric("Lowest off-target gene count", int(result.summary.off_target_genes.min()))
        st.dataframe(result.summary, hide_index=True, use_container_width=True)
        st.caption("Sorted by predicted activity, then off-target gene count. "
                   "Scores are not knockdown percentages; seed counts indicate potential binding sites.")
        if result.metadata.get("synthesis_status") == "template_required":
            st.warning("Mature guide sequences are available. Complete synthesis inserts require "
                       "a verified miR-E scaffold template.")
        with st.expander("Off-target site details"):
            st.dataframe(result.sites.head(10000), hide_index=True, use_container_width=True)
            if len(result.sites) > 10000:
                st.caption("Showing the first 10,000 sites. The download contains all sites.")
    st.download_button(
        "Download results and sequences", data=st.session_state.workbook,
        file_name=f"amir_{species}_{gene_name}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )


if __name__ == "__main__":
    main()
