import pandas as pd
import pytest

from src.multi_sample import (
    area_column,
    build_replicate_area_comparison,
    build_sample_presence_comparison,
    detected_column,
    replicate_area_view,
    unique_sample_labels,
)


def _frame(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_unique_sample_labels_suffix_duplicates() -> None:
    assert unique_sample_labels(["sample", "sample", "other", "sample"]) == [
        "sample", "sample (2)", "other", "sample (3)",
    ]


def test_replicate_comparison_chooses_lowest_area_cv_valid_combination() -> None:
    frames = {
        "rep-1": _frame([
            {
                "canonical_name": "Hexanal", "rt_min": 3.00, "ri": 690,
                "area": 100, "quality": 95,
            },
            {
                "canonical_name": "Hexanal", "rt_min": 3.04, "ri": 700,
                "area": 500, "quality": 90,
            },
        ]),
        "rep-2": _frame([
            {
                "canonical_name": "Hexanal", "rt_min": 3.08, "ri": 710,
                "area": 110, "quality": 92,
            },
            {
                "canonical_name": "Hexanal", "rt_min": 3.05, "ri": 705,
                "area": 900, "quality": 99,
            },
        ]),
    }

    comparison = build_replicate_area_comparison(frames)

    assert comparison["canonical_name"].tolist() == ["Hexanal"]
    row = comparison.iloc[0]
    assert row[area_column("rep-1")] == 100
    assert row[area_column("rep-2")] == 110
    assert row["area_mean"] == 105
    assert row["area_std"] == pytest.approx(7.0710678)
    assert row["sample_count"] == 2


def test_replicate_area_view_keeps_requested_columns_in_order() -> None:
    comparison = pd.DataFrame([{
        "canonical_name": "Hexanal",
        "sample_count": 2,
        "detected_samples": "rep-1, rep-2",
        "mean_rt": 3.04,
        "rt_range": 0.08,
        "mean_ri": 700,
        "ri_range": 20,
        area_column("rep-1"): 100,
        area_column("rep-2"): 110,
        "area_mean": 105,
        "area_std": 7.07,
        "area_cv_percent": 6.73,
    }])

    view = replicate_area_view(comparison)

    assert view.columns.tolist() == [
        "mean_rt",
        "canonical_name",
        "mean_ri",
        "sample_count",
        "area_mean",
        area_column("rep-1"),
        area_column("rep-2"),
        "detected_samples",
    ]


def test_replicate_comparison_requires_rt_and_ri_tolerances() -> None:
    frames = {
        "rep-1": _frame([{
            "canonical_name": "Hexanal", "rt_min": 3.00, "ri": 690,
            "area": 100, "quality": 95,
        }]),
        "rep-2": _frame([{
            "canonical_name": "Hexanal", "rt_min": 3.11, "ri": 710,
            "area": 105, "quality": 95,
        }]),
    }

    comparison = build_replicate_area_comparison(frames)

    assert comparison.empty


def test_different_sample_comparison_marks_common_and_unique_compounds() -> None:
    frames = {
        "sample-a": _frame([
            {
                "canonical_name": "Hexanal", "rt_min": 3.0, "ri": 690,
                "area": 100, "quality": 95,
            },
            {
                "canonical_name": "Octanal", "rt_min": 5.0, "ri": 890,
                "area": 200, "quality": 90,
            },
        ]),
        "sample-b": _frame([
            {
                "canonical_name": "Hexanal", "rt_min": 3.1, "ri": 700,
                "area": 120, "quality": 94,
            },
            {
                "canonical_name": "Nonanal", "rt_min": 6.0, "ri": 990,
                "area": 300, "quality": 91,
            },
        ]),
    }

    comparison = build_sample_presence_comparison(frames)

    status = dict(zip(comparison["canonical_name"], comparison["detection_status"]))
    assert status == {
        "Hexanal": "공통 검출",
        "Nonanal": "개별 검출",
        "Octanal": "개별 검출",
    }
    hexanal = comparison[comparison["canonical_name"] == "Hexanal"].iloc[0]
    assert hexanal[detected_column("sample-a")] == "●"
    assert hexanal[detected_column("sample-b")] == "●"
