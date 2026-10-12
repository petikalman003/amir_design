"""The single workbook retains large gene-wide fields without Excel truncation."""

from io import BytesIO
import json

import openpyxl
import pandas as pd
import pytest

from amir_design.pipeline import PipelineResult
from amir_design.reports import workbook_bytes, write_workbook


def test_large_candidate_fields_are_preserved_in_ordered_literal_parts():
    positions = "=" + ";".join(f"TEST_TX{i}:4,5,6" for i in range(5000))
    result = PipelineResult(
        status="complete",
        summary=pd.DataFrame({
            "candidate_id": ["amiR_0001"],
            "source_transcript_start_positions_1based": [positions],
        }),
        sites=pd.DataFrame(),
        on_target_sites=pd.DataFrame(),
        metadata={"target_isoforms": [f"TEST_TX{i}" for i in range(5000)]},
    )
    workbook = openpyxl.load_workbook(BytesIO(workbook_bytes(result)))
    assert "Full text recorded in Long values" in workbook["Candidates"]["B2"].value
    values = workbook["Long values"]
    rows = list(values.iter_rows(min_row=2, values_only=True))
    assert "".join(row[-1] for row in rows) == positions
    assert {row[1] for row in rows} == {"Candidates"}
    assert {row[2] for row in rows} == {"B2"}
    assert [row[4] for row in rows] == list(range(1, len(rows) + 1))
    assert values["G2"].data_type == "s"
    metadata = list(workbook["Run information"].iter_rows(min_row=2, values_only=True))
    isoform_parts = [value for key, value in metadata if key.startswith("target_isoforms")]
    assert json.loads("".join(isoform_parts)) == result.metadata["target_isoforms"]


class StreamingSites:
    columns = ("candidate_id", "transcript_id", "start", "target_sequence")

    def __init__(self, records):
        self.records = records
        self.iterations = 0

    def __len__(self):
        return len(self.records)

    def iter_rows(self):
        self.iterations += 1
        yield from self.records


def test_streamed_sites_split_sheets_without_losing_rows_or_literal_identifiers(monkeypatch):
    monkeypatch.setattr("amir_design.reports.MAX_DATA_ROWS", 2)
    sites = StreamingSites([
        ("a", "=1+1", 1, "ACGT"),
        ("a", "https://example.org/gene", 2, "ACGT"),
        ("a", "00123", 3, "ACGT"),
        ("a", "+SUM(A1:A2)", 4, "ACGT"),
        ("a", "last", 5, "ACGT"),
    ])
    result = PipelineResult(
        "complete", pd.DataFrame({"candidate_id": ["a"], "off_target_sites": [5]}),
        sites, pd.DataFrame(), {},
    )
    workbook = openpyxl.load_workbook(BytesIO(workbook_bytes(result)))
    sheets = [workbook[f"Off-target sites {index}"] for index in range(1, 4)]
    rows = [row for sheet in sheets for row in sheet.iter_rows(min_row=2, values_only=True)]
    assert rows == sites.records
    assert sites.iterations == 1
    assert [sheet.max_row for sheet in sheets] == [3, 3, 2]
    for sheet in sheets:
        assert sheet.freeze_panes == "A2"
        assert sheet["B2"].data_type == "s"
        assert sheet["B2"].hyperlink is None
    assert sheets[1]["B2"].value == "00123"
    metadata = dict(workbook["Run information"].iter_rows(min_row=2, values_only=True))
    assert metadata["report_detail_level"] == "all_sites"
    assert metadata["report_off_target_site_rows_included"] == "5"


def test_summary_export_does_not_read_positions_and_records_exact_counts():
    class UnreadSites(StreamingSites):
        def iter_rows(self):
            raise AssertionError("A summary must not read off-target positions")

    result = PipelineResult(
        "complete", pd.DataFrame({
            "candidate_id": ["a", "b"], "off_target_sites": [100000, 23456],
        }), UnreadSites([]), pd.DataFrame(), {"off_target_sites_retained": False},
    )
    workbook = openpyxl.load_workbook(
        BytesIO(workbook_bytes(result, include_off_target_sites=False)),
    )
    assert "omitted" in workbook["Off-target sites"]["A1"].value
    assert workbook["Candidates"].max_row == 3
    metadata = dict(workbook["Run information"].iter_rows(min_row=2, values_only=True))
    assert metadata["report_detail_level"] == "summary"
    assert metadata["report_off_target_sites_total"] == "123456"
    assert metadata["report_off_target_site_rows_included"] == "0"
    assert metadata["off_target_sites_retained"] == "false"


def test_detailed_export_rejects_missing_positions_instead_of_reporting_zero_sites():
    result = PipelineResult(
        "complete", pd.DataFrame({"off_target_sites": [100]}),
        pd.DataFrame(), pd.DataFrame(), {},
    )
    with pytest.raises(ValueError, match="positions were not retained"):
        workbook_bytes(result)


def test_file_export_has_complete_rows_and_reports_progress(tmp_path):
    sites = StreamingSites([("a", "TX1", 1, "ACGT"), ("a", "TX2", 2, "ACGT")])
    result = PipelineResult(
        "complete", pd.DataFrame({"candidate_id": ["a"], "off_target_sites": [2]}),
        sites, pd.DataFrame(), {},
    )
    messages = []
    destination = tmp_path / "complete.xlsx"
    assert write_workbook(result, destination, progress=messages.append) == destination
    workbook = openpyxl.load_workbook(destination)
    assert list(workbook["Off-target sites"].iter_rows(min_row=2, values_only=True)) == sites.records
    assert "Wrote Off-target sites: 2 rows" in messages
    assert messages[-1] == "Compressing the workbook"


def test_streamed_source_declared_row_count_is_verified():
    class ShortSites(StreamingSites):
        def __len__(self):
            return len(self.records) + 1

    result = PipelineResult(
        "complete", pd.DataFrame(), ShortSites([("a", "TX1", 1, "ACGT")]),
        pd.DataFrame(), {},
    )
    with pytest.raises(ValueError, match="ended before its declared row count"):
        workbook_bytes(result)


def test_long_value_links_point_to_the_actual_split_sheet(monkeypatch):
    monkeypatch.setattr("amir_design.reports.MAX_DATA_ROWS", 2)
    long_text = "=" + "ACGT" * 10000
    result = PipelineResult(
        "complete", pd.DataFrame({"candidate_id": ["a", "b", "c"],
                                  "text": ["short", "short", long_text]}),
        pd.DataFrame(), pd.DataFrame(), {},
    )
    workbook = openpyxl.load_workbook(BytesIO(workbook_bytes(result)))
    assert "Full text recorded in Long values" in workbook["Candidates 2"]["B2"].value
    parts = list(workbook["Long values"].iter_rows(min_row=2, values_only=True))
    assert "".join(row[-1] for row in parts) == long_text
    assert {row[1] for row in parts} == {"Candidates 2"}
    assert {row[2] for row in parts} == {"B2"}


def test_long_interpretation_notice_is_preserved_as_literal_parts():
    notice = "=" + "interpretation " * 5000
    result = PipelineResult("complete", pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {},
                            notices=[notice])
    workbook = openpyxl.load_workbook(BytesIO(workbook_bytes(result)))
    rows = list(workbook["Run information"].iter_rows(min_row=2, values_only=True))
    parts = [value for key, value in rows if key.startswith("Interpretation")]
    assert "".join(parts) == notice
    for row in workbook["Run information"].iter_rows(min_row=2):
        assert row[1].data_type == "s"
