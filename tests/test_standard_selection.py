import pandas as pd

from src.standard_selection import apply_selected_candidate_rts


STANDARDS = pd.DataFrame({
    "carbon_number": [6, 7, 8],
    "alkane_name": ["Hexane", "Heptane", "Octane"],
    "ri": [600, 700, 800],
    "rt_min": [2.0, 4.0, 8.0],
    "source": ["default", "default", "default"],
    "confirmed": [True, True, True],
})


def test_selected_candidate_replaces_only_matching_standard() -> None:
    candidates = pd.DataFrame({
        "selected": [False, True, True],
        "carbon_number": [6, 6, 8],
        "rt_min": [2.1, 2.2, 8.4],
    })
    updated = apply_selected_candidate_rts(STANDARDS, candidates)

    assert updated["rt_min"].tolist() == [2.2, 4.0, 8.4]
    assert updated.loc[updated["carbon_number"] == 6, "source"].iloc[0] == "사용자 선택 MassHunter RT"


def test_multiple_rts_for_same_compound_are_rejected() -> None:
    candidates = pd.DataFrame({
        "selected": [True, True],
        "carbon_number": [6, 6],
        "rt_min": [2.1, 2.2],
    })
    try:
        apply_selected_candidate_rts(STANDARDS, candidates)
    except ValueError as exc:
        assert "C6" in str(exc)
    else:
        raise AssertionError("duplicate compound selection should fail")


def test_empty_selection_is_rejected() -> None:
    candidates = pd.DataFrame({
        "selected": [False],
        "carbon_number": [6],
        "rt_min": [2.1],
    })
    try:
        apply_selected_candidate_rts(STANDARDS, candidates)
    except ValueError as exc:
        assert "선택" in str(exc)
    else:
        raise AssertionError("empty selection should fail")
