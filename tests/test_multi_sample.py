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
                "canonical_name": "Hexanal", "rt_min": 3.05, "ri": 710,
                "area": 110, "quality": 92,
            },
            {
                "canonical_name": "Hexanal", "rt_min": 3.04, "ri": 705,
                "area": 900, "quality": 99,
            },
        ]),
    }

    comparison = build_replicate_area_comparison(frames, ri_tolerance=30)

    # 두 쌍(A-C, B-D) 모두 RT/RI 허용범위를 만족하므로 생략되지 않고 각각 별도 행으로 남는다.
    assert comparison["canonical_name"].tolist() == ["Hexanal", "Hexanal"]
    row = comparison.iloc[0]
    assert row[area_column("rep-1")] == 100
    assert row[area_column("rep-2")] == 110
    assert row["area_mean"] == 105
    assert row["area_std"] == pytest.approx(7.0710678)
    assert row["sample_count"] == 2
    second = comparison.iloc[1]
    assert second[area_column("rep-1")] == 500
    assert second[area_column("rep-2")] == 900


def test_replicate_comparison_default_ri_tolerance_is_15() -> None:
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
                "canonical_name": "Hexanal", "rt_min": 3.05, "ri": 710,
                "area": 110, "quality": 92,
            },
            {
                "canonical_name": "Hexanal", "rt_min": 3.04, "ri": 705,
                "area": 900, "quality": 99,
            },
        ]),
    }

    comparison = build_replicate_area_comparison(frames)

    # A-C 쌍은 RI 차이가 20으로 기본 허용범위(±15)를 벗어나 유효하지 않다.
    assert comparison["canonical_name"].tolist() == ["Hexanal"]
    row = comparison.iloc[0]
    assert row[area_column("rep-1")] == 500
    assert row[area_column("rep-2")] == 900


def test_replicate_comparison_reports_multiple_distinct_clusters_with_same_name() -> None:
    frames = {
        "rep-1": _frame([
            {"canonical_name": "Aldehyde X", "rt_min": 5.00, "ri": 800, "area": 100, "quality": 90},
            {"canonical_name": "Aldehyde X", "rt_min": 7.00, "ri": 1000, "area": 200, "quality": 91},
        ]),
        "rep-2": _frame([
            {"canonical_name": "Aldehyde X", "rt_min": 5.01, "ri": 805, "area": 110, "quality": 92},
            {"canonical_name": "Aldehyde X", "rt_min": 7.01, "ri": 1005, "area": 210, "quality": 93},
        ]),
        "rep-3": _frame([
            {"canonical_name": "Aldehyde X", "rt_min": 5.02, "ri": 810, "area": 120, "quality": 94},
            {"canonical_name": "Aldehyde X", "rt_min": 7.02, "ri": 1010, "area": 220, "quality": 95},
        ]),
    }

    comparison = build_replicate_area_comparison(frames)

    assert comparison["canonical_name"].tolist() == ["Aldehyde X", "Aldehyde X"]
    assert sorted(comparison["mean_rt"].round(2).tolist()) == [5.01, 7.01]
    assert comparison["sample_count"].tolist() == [3, 3]


def test_replicate_comparison_rescues_quality_excluded_sample_above_ratio() -> None:
    frames = {
        "rep-1": _frame([{"canonical_name": "Ketone Y", "rt_min": 4.00, "ri": 600, "area": 100, "quality": 90}]),
        "rep-2": _frame([{"canonical_name": "Ketone Y", "rt_min": 4.01, "ri": 603, "area": 105, "quality": 91}]),
        "rep-3": _frame([{"canonical_name": "Ketone Y", "rt_min": 4.02, "ri": 606, "area": 110, "quality": 92}]),
        "rep-4": _frame([]),
    }
    rejected = {
        "rep-1": _frame([]),
        "rep-2": _frame([]),
        "rep-3": _frame([]),
        "rep-4": _frame([{
            "canonical_name": "Ketone Y", "rt_min": 4.015, "ri": 604, "area": 95,
            "quality": 50, "inclusion_reason": "",
        }]),
    }

    comparison = build_replicate_area_comparison(frames, rejected)

    assert len(comparison) == 1
    row = comparison.iloc[0]
    assert row["sample_count"] == 4
    assert row["quality_flag_samples"] == "rep-4"
    assert row[area_column("rep-4")] == 95


def test_replicate_comparison_does_not_rescue_below_min_ratio() -> None:
    frames = {
        "rep-1": _frame([{"canonical_name": "Ketone Y", "rt_min": 4.00, "ri": 600, "area": 100, "quality": 90}]),
        "rep-2": _frame([{"canonical_name": "Ketone Y", "rt_min": 4.01, "ri": 603, "area": 105, "quality": 91}]),
        "rep-3": _frame([]),
        "rep-4": _frame([]),
    }
    rejected = {
        "rep-1": _frame([]),
        "rep-2": _frame([]),
        "rep-3": _frame([{
            "canonical_name": "Ketone Y", "rt_min": 4.015, "ri": 604, "area": 95,
            "quality": 50, "inclusion_reason": "",
        }]),
        "rep-4": _frame([]),
    }

    comparison = build_replicate_area_comparison(frames, rejected)

    assert len(comparison) == 1
    row = comparison.iloc[0]
    assert row["sample_count"] == 2
    assert row["quality_flag_samples"] == ""


def test_replicate_comparison_ignores_siloxane_excluded_rejects_for_rescue() -> None:
    frames = {
        "rep-1": _frame([{"canonical_name": "Ketone Y", "rt_min": 4.00, "ri": 600, "area": 100, "quality": 90}]),
        "rep-2": _frame([{"canonical_name": "Ketone Y", "rt_min": 4.01, "ri": 603, "area": 105, "quality": 91}]),
        "rep-3": _frame([{"canonical_name": "Ketone Y", "rt_min": 4.02, "ri": 606, "area": 110, "quality": 92}]),
        "rep-4": _frame([]),
    }
    rejected = {
        "rep-1": _frame([]),
        "rep-2": _frame([]),
        "rep-3": _frame([]),
        "rep-4": _frame([{
            "canonical_name": "Ketone Y", "rt_min": 4.015, "ri": 604, "area": 95,
            "quality": 99, "inclusion_reason": "SILOXANE FAMILY EXCLUDED",
        }]),
    }

    comparison = build_replicate_area_comparison(frames, rejected)

    row = comparison.iloc[0]
    assert row["sample_count"] == 3
    assert row["quality_flag_samples"] == ""


def test_replicate_area_view_orders_replicates_descending_and_keeps_columns() -> None:
    comparison = pd.DataFrame([
        {
            "canonical_name": "Hexanal",
            "sample_count": 2,
            "detected_samples": "rep-1, rep-2",
            "mean_rt": 3.04,
            "rt_range": 0.05,
            "mean_ri": 700,
            "ri_range": 20,
            area_column("rep-1"): 100,
            area_column("rep-2"): 110,
            "area_mean": 105,
            "area_std": 7.07,
            "area_cv_percent": 6.73,
        },
        {
            "canonical_name": "Octanal",
            "sample_count": 3,
            "detected_samples": "rep-1, rep-2, rep-3",
            "mean_rt": 5.02,
            "rt_range": 0.04,
            "mean_ri": 900,
            "ri_range": 18,
            area_column("rep-1"): 200,
            area_column("rep-2"): 220,
            "area_mean": 210,
            "area_std": 10,
            "area_cv_percent": 4.76,
        },
    ])

    view = replicate_area_view(comparison)

    assert view["sample_count"].tolist() == [3, 2]
    assert view["canonical_name"].tolist() == ["Octanal", "Hexanal"]
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
            "canonical_name": "Hexanal", "rt_min": 3.051, "ri": 710,
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
