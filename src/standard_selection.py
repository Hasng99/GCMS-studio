from __future__ import annotations

import pandas as pd

from src.ri import validate_standards


def apply_selected_candidate_rts(
    standards: pd.DataFrame,
    candidate_rows: pd.DataFrame,
) -> pd.DataFrame:
    if "selected" not in candidate_rows.columns:
        raise ValueError("RT 후보의 선택 정보가 없습니다.")
    selected = candidate_rows[candidate_rows["selected"].fillna(False).astype(bool)].copy()
    if selected.empty:
        raise ValueError("교체할 RT 후보를 하나 이상 선택하세요.")
    required = {"carbon_number", "rt_min"}
    missing = required - set(selected.columns)
    if missing:
        raise ValueError("RT 후보 필수 열이 없습니다: " + ", ".join(sorted(missing)))
    duplicated = selected["carbon_number"].duplicated(keep=False)
    if duplicated.any():
        carbons = sorted(selected.loc[duplicated, "carbon_number"].astype(int).unique().tolist())
        raise ValueError(
            "같은 물질에서는 RT를 하나만 선택하세요: "
            + ", ".join(f"C{carbon}" for carbon in carbons)
        )

    updated = standards.copy()
    known_carbons = set(pd.to_numeric(updated["carbon_number"], errors="coerce").dropna().astype(int))
    selected_carbons = set(pd.to_numeric(selected["carbon_number"], errors="coerce").dropna().astype(int))
    unknown = sorted(selected_carbons - known_carbons)
    if unknown:
        raise ValueError("현재 Standard에 없는 탄소 수입니다: " + ", ".join(f"C{value}" for value in unknown))

    for row in selected.itertuples(index=False):
        carbon = int(row.carbon_number)
        updated.loc[updated["carbon_number"].astype(int) == carbon, "rt_min"] = float(row.rt_min)
        if "source" in updated.columns:
            updated.loc[updated["carbon_number"].astype(int) == carbon, "source"] = "사용자 선택 MassHunter RT"
        if "confirmed" in updated.columns:
            updated.loc[updated["carbon_number"].astype(int) == carbon, "confirmed"] = True
    return validate_standards(updated)
