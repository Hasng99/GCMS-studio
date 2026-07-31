from __future__ import annotations

import io

import pandas as pd

from .pipeline import PipelineResult


def dataframe_to_csv_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False).encode("utf-8-sig")


def results_to_xlsx_bytes(
    result: PipelineResult,
    standards: pd.DataFrame,
    profile: pd.DataFrame,
) -> bytes:
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        result.peak_summary.to_excel(writer, sheet_name="peak_summary", index=False)
        result.selected_hits.to_excel(writer, sheet_name="selected_hits", index=False)
        result.rejected_hits.to_excel(writer, sheet_name="rejected_hits", index=False)
        result.all_hits.to_excel(writer, sheet_name="all_hits", index=False)
        standards.to_excel(writer, sheet_name="standards", index=False)
        profile.to_excel(writer, sheet_name="profile", index=False)
    return output.getvalue()
