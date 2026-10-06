from io import BytesIO
from pathlib import Path
from types import SimpleNamespace

import openpyxl
import pandas as pd
import pytest

from amir_design.candidates import reverse_complement
from amir_design.pipeline import run_pipeline
from amir_design.reports import workbook_bytes


TARGET = "CGGCTGCGCATATACATATATA"


class ReferenceFixture:
    manifest = {"source": {"species": "test", "release": "fixture"}}
    fasta_path = Path("unused-fixture.fa")

    def __init__(self, missing_isoform=False):
        self.missing_isoform = missing_isoform

    def gene_for_cds(self, cds):
        assert cds == TARGET
        return {"gene_id": "target", "gene_name": "FixtureTarget"}

    def transcripts_for_gene(self, gene):
        return {"tx1": TARGET, "tx2": "ATG" * 30 if self.missing_isoform else "ATG" + TARGET}

    def transcript_gene_map(self):
        return {"tx1": "target", "tx2": "target", "off1": "other"}


class ScorerFixture:
    provenance = {"kind": "test stub"}

    def score(self, targets):
        assert targets == [TARGET]
        return [0.8]


class RunnerFixture:
    def screen(self, guides, fasta, gene, mapping):
        assert guides == {"amiR_0001": reverse_complement(TARGET)}
        assert gene == "target"
        return SimpleNamespace(
            summary={"amiR_0001": {"off_target_genes": 1, "off_target_transcripts": 1,
                                    "off_target_sites": 1}},
            sites=[{"candidate_id": "amiR_0001", "gene_id": "other", "transcript_id": "off1"}],
            provenance={"kind": "test stub"},
        )


def test_no_eligible_candidate_does_not_score_or_screen(tmp_path):
    result = run_pipeline(
        TARGET, ReferenceFixture(missing_isoform=True), object(), object(),
        template_path=tmp_path / "missing-template.json",
    )
    assert result.status == "no_eligible_candidate"
    assert result.summary.empty
    assert len(result.metadata["target_isoforms"]) == 2
    assert result.metadata["rejected_incomplete_isoform_coverage"] == 1
    workbook = openpyxl.load_workbook(BytesIO(workbook_bytes(result)))
    assert "Candidates" in workbook.sheetnames
    assert workbook["Candidates"].max_row == 1


def test_missing_template_preserves_real_stage_results_without_claiming_insert(monkeypatch):
    import amir_design.scaffold as scaffold

    monkeypatch.setattr(scaffold, "DEFAULT_MIRE_TEMPLATE", None)
    result = run_pipeline(TARGET, ReferenceFixture(), ScorerFixture(), RunnerFixture())
    assert result.summary.iloc[0].predicted_efficiency_score == 0.8
    assert result.summary.iloc[0].off_target_genes == 1
    assert result.summary.iloc[0].mire_insert_dna_5to3 == ""
    assert result.metadata["synthesis_status"] == "template_required"
    assert len(result.on_target_sites) == 2


def test_published_default_insert_contains_exact_scored_and_screened_guide():
    result = run_pipeline(TARGET, ReferenceFixture(), ScorerFixture(), RunnerFixture())
    row = result.summary.iloc[0]
    assert result.status == "complete"
    assert len(row.mire_insert_dna_5to3) == 119
    assert row.mire_insert_dna_5to3[68:90] == row.guide_dna_5to3
    assert result.metadata["synthesis_scaffold"]["template_origin"] == "verified published miR-E fragment"


def test_report_records_insert_provenance_and_literal_external_strings(monkeypatch):
    import amir_design.reports as reports

    result = run_pipeline(
        TARGET, ReferenceFixture(), ScorerFixture(), RunnerFixture(),
        insert_builder=lambda guide: "ACGT" + guide,
        scaffold_provenance={"source": "test fixture only"},
    )
    result.metadata["target_gene_name"] = "=SOME_EXTERNAL_IDENTIFIER"
    result.sites = pd.DataFrame({"gene_id": ["=literal", "+literal", "@literal", "-literal"]})
    monkeypatch.setattr(reports, "MAX_DATA_ROWS", 2)
    workbook = openpyxl.load_workbook(BytesIO(workbook_bytes(result)))
    assert "Off-target sites 2" in workbook.sheetnames
    cell = workbook["Off-target sites 1"]["A2"]
    assert cell.value == "=literal" and cell.data_type == "s"
    assert result.metadata["synthesis_scaffold"]["source"] == "test fixture only"


def test_stage_failure_is_not_reported_as_zero_offtargets():
    class BrokenRunner:
        def screen(self, *args):
            raise RuntimeError("R search failed")

    with pytest.raises(RuntimeError, match="R search failed"):
        run_pipeline(TARGET, ReferenceFixture(), ScorerFixture(), BrokenRunner())


def test_malformed_insert_is_not_treated_as_an_unavailable_template():
    def builder(guide):
        raise ValueError("invalid template")

    with pytest.raises(ValueError, match="invalid template"):
        run_pipeline(TARGET, ReferenceFixture(), ScorerFixture(), RunnerFixture(), insert_builder=builder)
