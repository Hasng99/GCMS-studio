import pandas as pd

from app import summary_view

from src.pipeline import run_pipeline


def _inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    hits = pd.DataFrame({
        "compound_number": [1, 1, 2], "rt_min": [3.0, 3.0, 5.0],
        "hit_number": [1, 2, 1], "hit_name": ["Hexanal", "Unknown A", "Unknown B"],
        "quality": [40, 90, 30], "cas_number": ["000066-25-1", "", ""],
        "area": [100, 100, 200],
    })
    profile = pd.DataFrame({
        "parent_fatty_acid": ["Linoleic acid"], "source_name": ["Hexanal"],
        "canonical_name": ["Hexanal"], "aliases_semicolon": [""],
    })
    standards = pd.DataFrame({
        "carbon_number": [6, 7, 8], "alkane_name": ["Hexane", "Heptane", "Octane"],
        "ri": [600, 700, 800], "rt_min": [2.0, 4.0, 8.0],
    })
    return hits, profile, standards


def test_profile_below_threshold_and_nonprofile_above_threshold() -> None:
    result = run_pipeline(*_inputs(), quality_threshold=80)
    reasons = set(result.selected_hits["inclusion_reason"])
    assert reasons == {"PROFILE", "QUALITY"}
    assert len(result.selected_hits) == 2
    assert len(result.rejected_hits) == 1


def test_both_is_not_duplicated() -> None:
    hits, profile, standards = _inputs()
    hits.loc[0, "quality"] = 95
    result = run_pipeline(hits, profile, standards, quality_threshold=80)
    both = result.selected_hits[result.selected_hits["inclusion_reason"] == "BOTH"]
    assert len(both) == 1


def test_summary_column_priority_and_hidden_identifiers() -> None:
    result = run_pipeline(*_inputs(), quality_threshold=80)
    summary = summary_view(result.peak_summary)
    assert summary.columns[:5].tolist() == [
        "rt_min", "canonical_name", "quality", "ri", "area"
    ]
    assert "sample_name" not in summary.columns
    assert "compound_number" not in summary.columns
