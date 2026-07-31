import pandas as pd

from src.ri import calculate_ri


STANDARDS = pd.DataFrame({
    "carbon_number": [6, 7, 8],
    "alkane_name": ["Hexane", "Heptane", "Octane"],
    "ri": [600, 700, 800],
    "rt_min": [2.0, 4.0, 8.0],
})


def test_midpoint_ri() -> None:
    assert calculate_ri(3.0, STANDARDS).ri == 650.0


def test_exact_standard() -> None:
    result = calculate_ri(4.0, STANDARDS)
    assert result.ri == 700.0
    assert result.lower_alkane == "Heptane"


def test_out_of_range() -> None:
    assert calculate_ri(1.0, STANDARDS).ri_status == "OUT_OF_RANGE"


def test_invalid_standard() -> None:
    bad = STANDARDS.assign(rt_min=[2.0, 2.0, 8.0])
    assert calculate_ri(3.0, bad).ri_status == "INVALID_STANDARD"
