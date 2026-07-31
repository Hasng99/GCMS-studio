from pathlib import Path

import pandas as pd

from src.parsers import parse_masshunter


def test_csv_parser_accepts_normalized_columns(tmp_path: Path) -> None:
    path = tmp_path / "sample.csv"
    pd.DataFrame({
        "compound_number": [1, None], "rt_min": [2.0, None], "area": [10, None],
        "hit_number": [1, 2], "hit_name": ["Hexanal", "Other"],
        "quality": [40, 90], "cas_number": ["66-25-1", "1-00-0"],
        "scan_number": [100, None], "baseline_height": [1, None],
        "absolute_height": [2, None], "peak_width_50_min": [0.1, None],
    }).to_csv(path, index=False)
    frame, _ = parse_masshunter(path)
    assert frame["compound_number"].tolist() == [1, 1]
    assert frame["rt_min"].tolist() == [2.0, 2.0]


def test_masshunter_xls_matches_reference() -> None:
    source = Path(__file__).parents[3] / "work" / "source" / "GCMS_RI_Codex_All_In_One"
    if not source.exists():
        return
    parsed, metadata = parse_masshunter(source / "Mix5_4.xls")
    expected = pd.read_csv(source / "masshunter_sample_normalized.csv", dtype={"cas_number": "string"})
    assert len(parsed) == 422
    assert metadata["header_row"] == "9"
    pd.testing.assert_series_equal(parsed["hit_name"].reset_index(drop=True), expected["hit_name"], check_names=False, check_dtype=False)
    pd.testing.assert_series_equal(parsed["rt_min"].reset_index(drop=True), expected["rt_min"], check_names=False, check_dtype=False)
