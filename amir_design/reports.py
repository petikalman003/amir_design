"""One workbook containing candidate sequences, scores, sites, and provenance."""

from __future__ import annotations

from io import BytesIO
import json

from openpyxl import Workbook
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter
import pandas as pd

from amir_design.pipeline import PipelineResult


MAX_DATA_ROWS = 1048575


def _literal(sheet, value):
    cell = WriteOnlyCell(sheet, value=value)
    if isinstance(value, str):
        if len(value) > 32767:
            raise ValueError(f"A report cell in {sheet.title} exceeds Excel's text limit.")
        cell.data_type = "s"
    return cell


def _write_sheet(workbook: Workbook, title: str, frame, long_values: list[dict] | None = None) -> None:
    # Keep one downloadable artifact even when full-transcript sites exceed a sheet.
    if len(frame) > MAX_DATA_ROWS:
        for index, start in enumerate(range(0, len(frame), MAX_DATA_ROWS), 1):
            _write_sheet(
                workbook, f"{title} {index}", frame.iloc[start:start + MAX_DATA_ROWS], long_values,
            )
        return
    sheet = workbook.create_sheet(title)
    if not len(frame.columns):
        sheet.append(["No records"])
        return
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = f"A1:{get_column_letter(len(frame.columns))}{len(frame) + 1}"
    headers = []
    for index, name in enumerate(frame.columns, 1):
        sheet.column_dimensions[get_column_letter(index)].width = min(55, max(16, len(name) + 2))
        cell = _literal(sheet, name)
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="245B55")
        headers.append(cell)
    sheet.append(headers)
    for row_index, record in enumerate(frame.itertuples(index=False, name=None), 2):
        cells = []
        for column_index, value in enumerate(record, 1):
            if hasattr(value, "item"):
                value = value.item()
            if isinstance(value, str) and len(value) > 32767 and long_values is not None:
                # Gene-wide searches can yield long lists of transcript IDs or positions.
                # Preserve every character in ordered parts and make the cell reference explicit.
                value_id = f"value_{len(long_values) + 1:06d}"
                parts = (len(value) + 29999) // 30000
                for part, start in enumerate(range(0, len(value), 30000), 1):
                    long_values.append({
                        "value_id": value_id,
                        "source_sheet": title,
                        "source_cell": f"{get_column_letter(column_index)}{row_index}",
                        "column_name": frame.columns[column_index - 1],
                        "part": part,
                        "parts": parts,
                        "value_text": value[start:start + 30000],
                    })
                value = f"Full text recorded in Long values ({value_id})."
            # Write all text as literal data, including external gene identifiers.
            cells.append(_literal(sheet, value))
        sheet.append(cells)


def workbook_bytes(result: PipelineResult) -> bytes:
    workbook = Workbook(write_only=True)
    long_values: list[dict] = []
    _write_sheet(workbook, "Candidates", result.summary, long_values)
    _write_sheet(workbook, "Off-target sites", result.sites, long_values)
    _write_sheet(workbook, "Target isoform sites", result.on_target_sites, long_values)
    if len(result.candidate_origins.columns):
        _write_sheet(workbook, "Candidate origins", result.candidate_origins, long_values)
    if long_values:
        _write_sheet(workbook, "Long values", pd.DataFrame(long_values))
    sheet = workbook.create_sheet("Run information")
    sheet.freeze_panes = "A2"
    sheet.column_dimensions["A"].width = 40
    sheet.column_dimensions["B"].width = 100
    sheet.append(["Field", "Value"])
    sheet.append(["status", _literal(sheet, result.status)])
    for key, value in result.metadata.items():
        text = value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)
        if len(text) <= 32767:
            sheet.append([key, _literal(sheet, text)])
        else:
            # Large provenance lists remain recoverable without silent truncation.
            for index, start in enumerate(range(0, len(text), 30000), 1):
                sheet.append([f"{key} (part {index})", _literal(sheet, text[start:start + 30000])])
    for notice in result.notices:
        sheet.append(["Interpretation", _literal(sheet, notice)])
    output = BytesIO()
    workbook.save(output)
    return output.getvalue()
