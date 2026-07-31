from __future__ import annotations

from itertools import combinations, product
from typing import Mapping

import pandas as pd

from .matching import normalize_name


AREA_PREFIX = "area::"
DETECTED_PREFIX = "detected::"
RT_PREFIX = "rt::"
RI_PREFIX = "ri::"


def area_column(sample_name: str) -> str:
    return f"{AREA_PREFIX}{sample_name}"


def detected_column(sample_name: str) -> str:
    return f"{DETECTED_PREFIX}{sample_name}"


def rt_column(sample_name: str) -> str:
    return f"{RT_PREFIX}{sample_name}"


def ri_column(sample_name: str) -> str:
    return f"{RI_PREFIX}{sample_name}"


def unique_sample_labels(sample_names: list[str]) -> list[str]:
    """중복 시료명을 화면과 비교표에서 구분 가능한 이름으로 만든다."""
    counts: dict[str, int] = {}
    used: set[str] = set()
    labels: list[str] = []
    for raw_name in sample_names:
        name = str(raw_name).strip() or "sample"
        next_count = counts.get(name, 0) + 1
        label = name if next_count == 1 else f"{name} ({next_count})"
        while label in used:
            next_count += 1
            label = f"{name} ({next_count})"
        counts[name] = next_count
        used.add(label)
        labels.append(label)
    return labels


def _comparison_candidates(
    frames: Mapping[str, pd.DataFrame],
) -> tuple[list[str], dict[str, list[dict[str, object]]]]:
    sample_names = list(frames)
    candidates: dict[str, list[dict[str, object]]] = {}
    row_order = 0
    for sample_name, frame in frames.items():
        if frame.empty:
            continue
        for record in frame.to_dict(orient="records"):
            name_key = normalize_name(record.get("canonical_name"))
            rt = pd.to_numeric(record.get("rt_min"), errors="coerce")
            ri = pd.to_numeric(record.get("ri"), errors="coerce")
            area = pd.to_numeric(record.get("area"), errors="coerce")
            if not name_key or pd.isna(rt) or pd.isna(ri) or pd.isna(area):
                continue
            record.update({
                "_sample_name": sample_name,
                "_name_key": name_key,
                "_rt": float(rt),
                "_ri": float(ri),
                "_area": float(area),
                "_row_order": row_order,
            })
            candidates.setdefault(name_key, []).append(record)
            row_order += 1
    return sample_names, candidates


def build_replicate_area_comparison(
    frames: Mapping[str, pd.DataFrame],
    *,
    rt_tolerance: float = 0.05,
    ri_tolerance: float = 30.0,
) -> pd.DataFrame:
    """반복시료에서 RT·RI 허용범위를 만족하며 Area CV가 가장 작은 조합을 찾는다."""
    sample_names, grouped_candidates = _comparison_candidates(frames)
    area_columns = [area_column(name) for name in sample_names]
    base_columns = [
        "canonical_name", "sample_count", "detected_samples",
        "mean_rt", "rt_range", "mean_ri", "ri_range",
        *area_columns, "area_mean", "area_std", "area_cv_percent",
    ]
    output_rows: list[dict[str, object]] = []

    for name_key, candidates in grouped_candidates.items():
        by_sample = {
            sample_name: [
                candidate
                for candidate in candidates
                if candidate["_sample_name"] == sample_name
            ]
            for sample_name in sample_names
        }
        available_samples = [
            sample_name for sample_name in sample_names if by_sample[sample_name]
        ]
        best_combo: tuple[dict[str, object], ...] | None = None
        best_score: tuple[object, ...] | None = None

        for sample_count in range(len(available_samples), 1, -1):
            for sample_subset in combinations(available_samples, sample_count):
                candidate_sets = [by_sample[sample_name] for sample_name in sample_subset]
                for combo in product(*candidate_sets):
                    rt_values = [float(item["_rt"]) for item in combo]
                    ri_values = [float(item["_ri"]) for item in combo]
                    rt_range = max(rt_values) - min(rt_values)
                    ri_range = max(ri_values) - min(ri_values)
                    if (
                        rt_range > rt_tolerance + 1e-12
                        or ri_range > ri_tolerance + 1e-12
                    ):
                        continue
                    area_values = pd.Series(
                        [float(item["_area"]) for item in combo],
                        dtype="float64",
                    )
                    area_mean = float(area_values.mean())
                    area_std = float(area_values.std(ddof=1))
                    if area_mean == 0:
                        area_cv = 0.0 if area_std == 0 else float("inf")
                    else:
                        area_cv = abs(area_std / area_mean) * 100
                    score: tuple[object, ...] = (
                        area_cv,
                        area_std,
                        rt_range,
                        ri_range,
                        tuple(int(item["_row_order"]) for item in combo),
                    )
                    if best_score is None or score < best_score:
                        best_combo = combo
                        best_score = score
            if best_combo is not None:
                break

        if best_combo is None:
            continue

        rt_values = pd.Series(
            [float(item["_rt"]) for item in best_combo],
            dtype="float64",
        )
        ri_values = pd.Series(
            [float(item["_ri"]) for item in best_combo],
            dtype="float64",
        )
        area_values = pd.Series(
            [float(item["_area"]) for item in best_combo],
            dtype="float64",
        )
        area_mean = float(area_values.mean())
        area_std = float(area_values.std(ddof=1))
        if area_mean == 0:
            area_cv = 0.0 if area_std == 0 else float("inf")
        else:
            area_cv = abs(area_std / area_mean) * 100
        selected_samples = [str(item["_sample_name"]) for item in best_combo]
        representative = max(
            best_combo,
            key=lambda item: (
                pd.to_numeric(item.get("quality"), errors="coerce")
                if pd.notna(pd.to_numeric(item.get("quality"), errors="coerce"))
                else float("-inf")
            ),
        )
        output_row: dict[str, object] = {
            "canonical_name": representative.get("canonical_name", name_key),
            "sample_count": len(best_combo),
            "detected_samples": ", ".join(selected_samples),
            "mean_rt": float(rt_values.mean()),
            "rt_range": float(rt_values.max() - rt_values.min()),
            "mean_ri": float(ri_values.mean()),
            "ri_range": float(ri_values.max() - ri_values.min()),
            **{column: pd.NA for column in area_columns},
            "area_mean": area_mean,
            "area_std": area_std,
            "area_cv_percent": area_cv,
        }
        for item in best_combo:
            output_row[area_column(str(item["_sample_name"]))] = float(item["_area"])
        output_rows.append(output_row)

    if not output_rows:
        return pd.DataFrame(columns=base_columns)
    return (
        pd.DataFrame(output_rows, columns=base_columns)
        .sort_values(["area_cv_percent", "canonical_name"], kind="stable")
        .reset_index(drop=True)
    )


def replicate_area_view(frame: pd.DataFrame) -> pd.DataFrame:
    """반복시료 Area 비교 화면에 필요한 열만 요청 순서대로 반환한다."""
    area_columns = [
        column for column in frame.columns if column.startswith(AREA_PREFIX)
    ]
    columns = [
        "mean_rt",
        "canonical_name",
        "mean_ri",
        "sample_count",
        "area_mean",
        *area_columns,
        "detected_samples",
    ]
    return frame[[column for column in columns if column in frame.columns]].copy()


def build_sample_presence_comparison(
    frames: Mapping[str, pd.DataFrame],
) -> pd.DataFrame:
    """서로 다른 샘플의 공통·부분 공통·개별 검출 물질을 비교한다."""
    sample_names = list(frames)
    representatives: dict[str, dict[str, dict[str, object]]] = {}

    for sample_name, frame in frames.items():
        if frame.empty:
            continue
        working = frame.copy()
        working["_name_key"] = working["canonical_name"].map(normalize_name)
        working["_quality_sort"] = pd.to_numeric(
            working.get("quality"), errors="coerce"
        )
        working["_area_sort"] = pd.to_numeric(
            working.get("area"), errors="coerce"
        )
        working = working[working["_name_key"] != ""].sort_values(
            ["_quality_sort", "_area_sort"],
            ascending=[False, False],
            na_position="last",
            kind="stable",
        )
        for record in working.drop_duplicates("_name_key").to_dict(orient="records"):
            representatives.setdefault(str(record["_name_key"]), {})[sample_name] = record

    columns = [
        "canonical_name", "detection_status", "detected_count", "detected_samples",
    ]
    for sample_name in sample_names:
        columns.extend([
            detected_column(sample_name),
            rt_column(sample_name),
            ri_column(sample_name),
            area_column(sample_name),
        ])

    rows: list[dict[str, object]] = []
    for records_by_sample in representatives.values():
        detected_samples = [
            sample_name for sample_name in sample_names if sample_name in records_by_sample
        ]
        detected_count = len(detected_samples)
        if detected_count == len(sample_names):
            status = "공통 검출"
        elif detected_count == 1:
            status = "개별 검출"
        else:
            status = "부분 공통"
        first_record = records_by_sample[detected_samples[0]]
        row: dict[str, object] = {
            "canonical_name": first_record.get("canonical_name", ""),
            "detection_status": status,
            "detected_count": detected_count,
            "detected_samples": ", ".join(detected_samples),
        }
        for sample_name in sample_names:
            record = records_by_sample.get(sample_name)
            row[detected_column(sample_name)] = "●" if record is not None else "—"
            row[rt_column(sample_name)] = (
                record.get("rt_min") if record is not None else pd.NA
            )
            row[ri_column(sample_name)] = (
                record.get("ri") if record is not None else pd.NA
            )
            row[area_column(sample_name)] = (
                record.get("area") if record is not None else pd.NA
            )
        rows.append(row)

    if not rows:
        return pd.DataFrame(columns=columns)
    status_order = {"공통 검출": 0, "부분 공통": 1, "개별 검출": 2}
    output = pd.DataFrame(rows, columns=columns)
    output["_status_order"] = output["detection_status"].map(status_order)
    return (
        output.sort_values(
            ["_status_order", "canonical_name"],
            kind="stable",
        )
        .drop(columns="_status_order")
        .reset_index(drop=True)
    )
