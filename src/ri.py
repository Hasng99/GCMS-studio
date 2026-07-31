from __future__ import annotations

import bisect
from dataclasses import dataclass, asdict

import pandas as pd


@dataclass(frozen=True)
class RIResult:
    lower_alkane: str = ""
    upper_alkane: str = ""
    lower_rt: float | None = None
    upper_rt: float | None = None
    ri: float | None = None
    ri_status: str = "INVALID_STANDARD"

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def validate_standards(standards: pd.DataFrame) -> pd.DataFrame:
    required = {"carbon_number", "alkane_name", "rt_min"}
    missing = required - set(standards.columns)
    if missing:
        raise ValueError("Standard 필수 열이 없습니다: " + ", ".join(sorted(missing)))
    clean = standards.copy()
    clean["carbon_number"] = pd.to_numeric(clean["carbon_number"], errors="coerce")
    clean["rt_min"] = pd.to_numeric(clean["rt_min"], errors="coerce")
    if clean[["carbon_number", "rt_min"]].isna().any().any() or len(clean) < 2:
        raise ValueError("Standard에는 유효한 carbon_number와 rt_min이 2개 이상 필요합니다.")
    clean = clean.sort_values("rt_min").reset_index(drop=True)
    if not clean["rt_min"].is_monotonic_increasing or clean["rt_min"].duplicated().any():
        raise ValueError("Standard RT는 중복 없이 엄격히 증가해야 합니다.")
    if (clean["carbon_number"].diff().dropna() <= 0).any():
        raise ValueError("Standard 탄소 수는 RT와 함께 엄격히 증가해야 합니다.")
    return clean


def calculate_ri(
    rt_min: float,
    standards: pd.DataFrame,
    *,
    allow_extrapolation: bool = False,
    round_digits: int = 1,
    exact_tolerance: float = 0.003,
) -> RIResult:
    try:
        clean = validate_standards(standards)
        rt = float(rt_min)
    except (TypeError, ValueError):
        return RIResult()
    rts = clean["rt_min"].astype(float).tolist()
    exact = min(range(len(rts)), key=lambda i: abs(rts[i] - rt))
    if abs(rts[exact] - rt) <= exact_tolerance:
        row = clean.iloc[exact]
        ri = float(row.get("ri", 100 * row["carbon_number"]))
        return RIResult(str(row["alkane_name"]), str(row["alkane_name"]), rts[exact], rts[exact], round(ri, round_digits), "OK")
    position = bisect.bisect_left(rts, rt)
    if position == 0 or position == len(rts):
        if not allow_extrapolation:
            return RIResult(ri_status="OUT_OF_RANGE")
        lower_i, upper_i = (0, 1) if position == 0 else (len(rts) - 2, len(rts) - 1)
    else:
        lower_i, upper_i = position - 1, position
    lower, upper = clean.iloc[lower_i], clean.iloc[upper_i]
    lower_rt, upper_rt = float(lower["rt_min"]), float(upper["rt_min"])
    lower_ri = float(lower.get("ri", 100 * lower["carbon_number"]))
    upper_ri = float(upper.get("ri", 100 * upper["carbon_number"]))
    value = lower_ri + (upper_ri - lower_ri) * (rt - lower_rt) / (upper_rt - lower_rt)
    return RIResult(
        str(lower["alkane_name"]),
        str(upper["alkane_name"]),
        lower_rt,
        upper_rt,
        round(value, round_digits),
        "OK",
    )
