from __future__ import annotations

import io
from typing import Mapping

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


def multi_results_to_xlsx_bytes(
    results: Mapping[str, PipelineResult],
    standards: pd.DataFrame,
    profile: pd.DataFrame,
    *,
    comparison: pd.DataFrame | None = None,
    comparison_sheet: str = "comparison",
) -> bytes:
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        if comparison is not None:
            comparison.to_excel(writer, sheet_name=comparison_sheet[:31], index=False)
        pd.DataFrame({
            "sample_index": range(1, len(results) + 1),
            "sample_name": list(results),
        }).to_excel(writer, sheet_name="samples", index=False)
        for index, result in enumerate(results.values(), start=1):
            prefix = f"sample_{index}"
            result.peak_summary.to_excel(
                writer, sheet_name=f"{prefix}_summary", index=False
            )
            result.selected_hits.to_excel(
                writer, sheet_name=f"{prefix}_selected", index=False
            )
            result.rejected_hits.to_excel(
                writer, sheet_name=f"{prefix}_rejected", index=False
            )
            result.all_hits.to_excel(
                writer, sheet_name=f"{prefix}_all", index=False
            )
        standards.to_excel(writer, sheet_name="standards", index=False)
        profile.to_excel(writer, sheet_name="profile", index=False)
    return output.getvalue()
