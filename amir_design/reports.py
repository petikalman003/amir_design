"""Workbooks with candidate summaries and optional complete, streamed site detail."""

from __future__ import annotations

from io import BytesIO
import json
from pathlib import Path
from typing import Callable

import pandas as pd
import xlsxwriter
from xlsxwriter.utility import xl_col_to_name

from amir_design.pipeline import PipelineResult


MAX_DATA_ROWS = 1048575
_TEXT_LIMIT = 32767
_TEXT_PART_LENGTH = 30000


def _rows(source):
    """Read either a small dataframe or the disk-backed off-target site store."""
    if hasattr(source, "iter_rows"):
        return iter(source.iter_rows())
    return iter(source.itertuples(index=False, name=None))


def _value(value):
    if hasattr(value, "item"):
        value = value.item()
    if value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    return value


def _write_sheet(
    workbook,
    title: str,
    frame,
    long_values: list[dict] | None = None,
    progress: Callable[[str], None] | None = None,
) -> None:
    """Write every row once, splitting sheets at Excel's actual row limit."""
    columns = list(frame.columns)
    count = len(frame)
    rows = _rows(frame)
    sheet_count = max(1, (count + MAX_DATA_ROWS - 1) // MAX_DATA_ROWS)
    header_format = workbook.add_format({
        "bold": True, "font_color": "FFFFFF", "bg_color": "245B55",
    })
    written = 0
    for sheet_index in range(sheet_count):
        sheet_title = f"{title} {sheet_index + 1}" if sheet_count > 1 else title
        sheet = workbook.add_worksheet(sheet_title)
        if not columns:
            sheet.write_string(0, 0, "No records")
            continue
        sheet.freeze_panes(1, 0)
        for column_index, name in enumerate(columns):
            sheet.set_column(column_index, column_index, min(55, max(16, len(name) + 2)))
        sheet.write_row(0, 0, columns, header_format)
        sheet_rows = min(MAX_DATA_ROWS, count - written)
        sheet.autofilter(0, 0, sheet_rows, len(columns) - 1)
        for row_index in range(1, sheet_rows + 1):
            try:
                record = next(rows)
            except StopIteration as exc:
                raise ValueError(f"The {title} row source ended before its declared row count.") from exc
            values = []
            for column_index, raw in enumerate(record):
                value = _value(raw)
                if isinstance(value, str) and len(value) > _TEXT_LIMIT:
                    if long_values is None:
                        raise ValueError(f"A report cell in {sheet_title} exceeds Excel's text limit.")
                    value_id = f"value_{len(long_values) + 1:06d}"
                    parts = (len(value) + _TEXT_PART_LENGTH - 1) // _TEXT_PART_LENGTH
                    for part, start in enumerate(range(0, len(value), _TEXT_PART_LENGTH), 1):
                        long_values.append({
                            "value_id": value_id,
                            "source_sheet": sheet_title,
                            "source_cell": f"{xl_col_to_name(column_index)}{row_index + 1}",
                            "column_name": columns[column_index],
                            "part": part,
                            "parts": parts,
                            "value_text": value[start:start + _TEXT_PART_LENGTH],
                        })
                    value = f"Full text recorded in Long values ({value_id})."
                values.append(value)
            # Automatic formula, URL, and number conversion are disabled globally.
            # write_row is faster than allocating a Cell for each value.
            sheet.write_row(row_index, 0, values)
            written += 1
            if progress and written % 25000 == 0:
                progress(f"Writing {title}: {written:,}/{count:,} rows")
    if next(rows, None) is not None:
        raise ValueError(f"The {title} row source exceeds its declared row count.")
    if progress and count:
        progress(f"Wrote {title}: {count:,} rows")


def _write_workbook(
    result: PipelineResult,
    destination,
    *,
    include_off_target_sites: bool,
    progress: Callable[[str], None] | None,
) -> None:
    total_sites = (
        int(result.summary["off_target_sites"].sum())
        if "off_target_sites" in result.summary.columns else len(result.sites)
    )
    if include_off_target_sites and len(result.sites) != total_sites:
        raise ValueError(
            "The complete off-target positions were not retained. Prepare site details "
            "before generating the detailed workbook."
        )
    options = {
        "constant_memory": True,
        "strings_to_formulas": False,
        "strings_to_urls": False,
        "strings_to_numbers": False,
    }
    with xlsxwriter.Workbook(destination, options) as workbook:
        # Complete whole-transcript reports can exceed ZIP's legacy 4 GB limit.
        workbook.use_zip64()
        long_values: list[dict] = []
        _write_sheet(workbook, "Candidates", result.summary, long_values, progress)
        if include_off_target_sites:
            _write_sheet(workbook, "Off-target sites", result.sites, long_values, progress)
        else:
            sheet = workbook.add_worksheet("Off-target sites")
            sheet.set_column(0, 0, 100)
            sheet.write_string(0, 0, "Individual off-target positions are omitted from this summary.")
            sheet.write_string(1, 0, "Candidate counts are exact. Prepare the detailed workbook "
                               "to download every off-target site.")
        _write_sheet(workbook, "Target isoform sites", result.on_target_sites, long_values, progress)
        if len(result.candidate_origins.columns):
            _write_sheet(workbook, "Candidate origins", result.candidate_origins,
                         long_values, progress)
        if long_values:
            _write_sheet(workbook, "Long values", pd.DataFrame(long_values), progress=progress)
        sheet = workbook.add_worksheet("Run information")
        sheet.freeze_panes(1, 0)
        sheet.set_column(0, 0, 40)
        sheet.set_column(1, 1, 100)
        sheet.write_row(0, 0, ["Field", "Value"])
        metadata = {
            **result.metadata,
            "report_detail_level": "all_sites" if include_off_target_sites else "summary",
            "report_off_target_sites_total": total_sites,
            "report_off_target_site_rows_included": len(result.sites)
            if include_off_target_sites else 0,
        }
        information = [("status", result.status), *metadata.items()]
        row_index = 1
        for key, value in information:
            text = value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)
            if len(text) <= _TEXT_LIMIT:
                sheet.write_row(row_index, 0, [key, text])
                row_index += 1
            else:
                for index, start in enumerate(range(0, len(text), _TEXT_PART_LENGTH), 1):
                    sheet.write_row(row_index, 0, [f"{key} (part {index})",
                                                   text[start:start + _TEXT_PART_LENGTH]])
                    row_index += 1
        for notice in result.notices:
            if len(notice) <= _TEXT_LIMIT:
                sheet.write_row(row_index, 0, ["Interpretation", notice])
                row_index += 1
            else:
                for index, start in enumerate(range(0, len(notice), _TEXT_PART_LENGTH), 1):
                    sheet.write_row(row_index, 0, [f"Interpretation (part {index})",
                                                   notice[start:start + _TEXT_PART_LENGTH]])
                    row_index += 1
        if progress:
            progress("Compressing the workbook")


def write_workbook(
    result: PipelineResult,
    path: Path,
    *,
    include_off_target_sites: bool = True,
    progress: Callable[[str], None] | None = None,
) -> Path:
    """Stream a report to disk without keeping its complete ZIP bytes in memory."""
    path = Path(path)
    _write_workbook(result, str(path), include_off_target_sites=include_off_target_sites,
                    progress=progress)
    return path


def workbook_bytes(
    result: PipelineResult,
    *,
    include_off_target_sites: bool = True,
    progress: Callable[[str], None] | None = None,
) -> bytes:
    """Return a downloadable workbook; use write_workbook for large detail exports."""
    output = BytesIO()
    _write_workbook(result, output, include_off_target_sites=include_off_target_sites,
                    progress=progress)
    return output.getvalue()
