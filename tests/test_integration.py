"""Real BLAST -> model -> R tests on deliberately synthetic transcript references."""

from io import BytesIO
from pathlib import Path

import numpy as np
import openpyxl
import pytest
from streamlit.testing.v1 import AppTest

from amir_design.offtargets import SeedMatchRunner
from amir_design.pipeline import run_pipeline
from amir_design.references import prepare_reference_bundle
from amir_design.reports import workbook_bytes
from amir_design.scoring import ShRNAIScorer


ROOT = Path(__file__).resolve().parents[1]
TARGET = "CGGCTGCGCATATACATATATA"


@pytest.fixture(scope="module")
def real_reference(tmp_path_factory):
    root = tmp_path_factory.mktemp("integration") / "human"
    return prepare_reference_bundle(
        root, species="human", source_fasta=ROOT / "tests/fixtures/integration_transcripts.fa",
        release="synthetic-fixture", assembly="synthetic-fixture",
    )


@pytest.fixture(scope="module")
def real_dependencies():
    scorer = ShRNAIScorer(ROOT / "tools/shrnai")
    runner = SeedMatchRunner(ROOT / "tools/seedmatchr", rscript="/workspace/.runtime/bin/Rscript")
    return scorer, runner


def test_real_three_stage_pipeline_and_download(real_reference, real_dependencies):
    result = run_pipeline(TARGET, real_reference, *real_dependencies)
    assert result.status == "complete"
    assert len(result.summary) == 1
    row = result.summary.iloc[0]
    assert np.isfinite(row.predicted_efficiency_score)
    assert row.target_isoforms_covered == row.target_isoforms_total == 2
    assert row.off_target_genes == 2
    assert row.off_target_transcripts == row.off_target_sites == 3
    assert row.mer8_sites == 3
    assert set(result.sites.gene_id) == {"TEST_OFF_GENE1", "TEST_OFF_GENE2"}
    assert result.metadata["synthesis_status"] == "insert_available"
    assert len(row.mire_insert_dna_5to3) == 119
    assert row.mire_insert_dna_5to3[68:90] == row.guide_dna_5to3
    assert result.metadata["species_reference"]["source"]["release"] == "synthetic-fixture"
    report = openpyxl.load_workbook(BytesIO(workbook_bytes(result)))
    assert report["Candidates"].max_row == 2
    assert report["Off-target sites"].max_row == 4
    assert report["Target isoform sites"].max_row == 3


def test_real_reference_without_shared_isoform_returns_no_candidate(tmp_path):
    reference = prepare_reference_bundle(
        tmp_path / "no-shared", species="mouse",
        source_fasta=ROOT / "tests/fixtures/no_shared_isoform_transcripts.fa",
        release="synthetic-fixture", assembly="synthetic-fixture",
    )
    result = run_pipeline(TARGET, reference, object(), object())
    assert result.status == "no_eligible_candidate"
    assert result.summary.empty


def test_streamlit_real_submission_and_failed_submission_clear_old_result(
    real_reference, monkeypatch,
):
    monkeypatch.setenv("AMIR_REFERENCE_ROOT", str(real_reference.root.parent))
    monkeypatch.setenv("AMIR_TOOLS_ROOT", str(ROOT / "tools"))
    monkeypatch.setenv("AMIR_RSCRIPT", "/workspace/.runtime/bin/Rscript")
    monkeypatch.delenv("AMIR_MIRE_TEMPLATE", raising=False)
    app = AppTest.from_file(str(ROOT / "app.py"), default_timeout=30).run()
    assert not app.exception
    app.text_area[0].set_value(TARGET)
    app.button[0].click().run()
    assert not app.exception and not app.error
    assert app.session_state["result"].summary.iloc[0].off_target_genes == 2
    assert len(app.session_state["workbook"]) > 1000
    app.text_area[0].set_value("invalid sequence")
    app.button[0].click().run()
    assert not app.exception
    assert len(app.error) == 1
    assert "result" not in app.session_state
