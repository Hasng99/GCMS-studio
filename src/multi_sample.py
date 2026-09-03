from __future__ import annotations

from itertools import combinations, product
from typing import Mapping

import pandas as pd

from .matching import normalize_name


AREA_PREFIX = "area::"
DETECTED_PREFIX = "detected::"
RT_PREFIX = "rt::"
RI_PREFIX = "ri::"

DEFAULT_RT_TOLERANCE = 0.05
DEFAULT_RI_TOLERANCE = 15.0
MIN_REPLICATE_DETECTION_RATIO = 0.7


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


def _rejected_pool(
    rejected_frames: Mapping[str, pd.DataFrame],
) -> dict[str, dict[str, list[dict[str, object]]]]:
    """구제 후보 풀: sample_name -> name_key -> quality 미달 등으로 제외된 hit 목록.

    siloxane 계열 제외처럼 quality와 무관한 사유로 빠진 hit은 구제 대상에서 뺀다.
    """
    pool: dict[str, dict[str, list[dict[str, object]]]] = {}
    for sample_name, frame in rejected_frames.items():
        if frame.empty:
            continue
        bucket: dict[str, list[dict[str, object]]] = {}
        for record in frame.to_dict(orient="records"):
            if record.get("inclusion_reason") == "SILOXANE FAMILY EXCLUDED":
                continue
            name_key = normalize_name(record.get("canonical_name"))
            rt = pd.to_numeric(record.get("rt_min"), errors="coerce")
            ri = pd.to_numeric(record.get("ri"), errors="coerce")
            area = pd.to_numeric(record.get("area"), errors="coerce")
            if not name_key or pd.isna(rt) or pd.isna(ri) or pd.isna(area):
                continue
            enriched = dict(record)
            enriched.update({
                "_sample_name": sample_name,
                "_rt": float(rt),
                "_ri": float(ri),
                "_area": float(area),
            })
            bucket.setdefault(name_key, []).append(enriched)
        pool[sample_name] = bucket
    return pool


def _search_best_combo(
    by_sample: dict[str, list[dict[str, object]]],
    available_samples: list[str],
    rt_tolerance: float,
    ri_tolerance: float,
) -> tuple[dict[str, object], ...] | None:
    """RT·RI 허용범위 안에서 Area CV가 가장 작은, 가능한 많은 샘플의 조합을 찾는다."""
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
                if rt_range > rt_tolerance + 1e-12 or ri_range > ri_tolerance + 1e-12:
                    continue
                area_values = pd.Series(
                    [float(item["_area"]) for item in combo], dtype="float64"
                )
                area_mean = float(area_values.mean())
                area_std = float(area_values.std(ddof=1))
                area_cv = (
                    (0.0 if area_std == 0 else float("inf"))
                    if area_mean == 0
                    else abs(area_std / area_mean) * 100
                )
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
    return best_combo


def build_replicate_area_comparison(
    frames: Mapping[str, pd.DataFrame],
    rejected_frames: Mapping[str, pd.DataFrame] | None = None,
    *,
    rt_tolerance: float = DEFAULT_RT_TOLERANCE,
    ri_tolerance: float = DEFAULT_RI_TOLERANCE,
    min_detection_ratio: float = MIN_REPLICATE_DETECTION_RATIO,
) -> pd.DataFrame:
    """반복시료에서 RT·RI 허용범위를 만족하는 동일 물질 조합을 모두 찾는다.

    같은 이름이라도 RT/RI가 서로 다른 군집이 여러 개면 하나만 남기지 않고
    각 군집을 별도 행으로 제시한다. 조합이 전체 반복 수의
    `min_detection_ratio` 이상을 만족하면, 나머지 샘플 중 quality 미달로
    제외됐던 동일 물질(rejected_frames)도 함께 표시하고
    `quality_flag_samples`에 어느 샘플이 구제됐는지 기록한다.
    """
    sample_names, grouped_candidates = _comparison_candidates(frames)
    rescue_pool = _rejected_pool(rejected_frames or {})
    sample_total = len(sample_names)
    area_columns = [area_column(name) for name in sample_names]
    base_columns = [
        "canonical_name", "sample_count", "detected_samples",
        "mean_rt", "rt_range", "mean_ri", "ri_range",
        *area_columns, "area_mean", "area_std", "area_cv_percent",
        "quality_flag_samples",
    ]
    output_rows: list[dict[str, object]] = []

    for name_key, candidates in grouped_candidates.items():
        by_sample = {
            sample_name: [
                candidate for candidate in candidates if candidate["_sample_name"] == sample_name
            ]
            for sample_name in sample_names
        }
        rejected_by_sample = {
            sample_name: list(rescue_pool.get(sample_name, {}).get(name_key, []))
            for sample_name in sample_names
        }

        while True:
            available_samples = [
                sample_name for sample_name in sample_names if by_sample[sample_name]
            ]
            if len(available_samples) < 2:
                break
            best_combo = _search_best_combo(by_sample, available_samples, rt_tolerance, ri_tolerance)
            if best_combo is None:
                break

            used_ids = {id(item) for item in best_combo}
            for sample_name in sample_names:
                by_sample[sample_name] = [
                    item for item in by_sample[sample_name] if id(item) not in used_ids
                ]

            selected_samples = {str(item["_sample_name"]) for item in best_combo}
            rescued: dict[str, dict[str, object]] = {}
            if sample_total and len(best_combo) >= min_detection_ratio * sample_total - 1e-9:
                rt_low = min(float(item["_rt"]) for item in best_combo) - rt_tolerance
                rt_high = max(float(item["_rt"]) for item in best_combo) + rt_tolerance
                ri_low = min(float(item["_ri"]) for item in best_combo) - ri_tolerance
                ri_high = max(float(item["_ri"]) for item in best_combo) + ri_tolerance
                for sample_name in sample_names:
                    if sample_name in selected_samples:
                        continue
                    pool = rejected_by_sample.get(sample_name, [])
                    matches = [
                        item for item in pool
                        if rt_low <= item["_rt"] <= rt_high and ri_low <= item["_ri"] <= ri_high
                    ]
                    if not matches:
                        continue
                    best_match = max(
                        matches,
                        key=lambda item: (
                            pd.to_numeric(item.get("quality"), errors="coerce")
                            if pd.notna(pd.to_numeric(item.get("quality"), errors="coerce"))
                            else float("-inf")
                        ),
                    )
                    rescued[sample_name] = best_match
                    rejected_by_sample[sample_name] = [
                        item for item in pool if item is not best_match
                    ]

            all_members = list(best_combo) + list(rescued.values())
            rt_values = pd.Series([float(item["_rt"]) for item in all_members], dtype="float64")
            ri_values = pd.Series([float(item["_ri"]) for item in all_members], dtype="float64")
            area_values = pd.Series([float(item["_area"]) for item in all_members], dtype="float64")
            area_mean = float(area_values.mean())
            area_std = float(area_values.std(ddof=1))
            area_cv = (
                (0.0 if area_std == 0 else float("inf"))
                if area_mean == 0
                else abs(area_std / area_mean) * 100
            )
            selected_sample_names = sorted({str(item["_sample_name"]) for item in all_members})
            representative = max(
                all_members,
                key=lambda item: (
                    pd.to_numeric(item.get("quality"), errors="coerce")
                    if pd.notna(pd.to_numeric(item.get("quality"), errors="coerce"))
                    else float("-inf")
                ),
            )
            output_row: dict[str, object] = {
                "canonical_name": representative.get("canonical_name", name_key),
                "sample_count": len(all_members),
                "detected_samples": ", ".join(selected_sample_names),
                "mean_rt": float(rt_values.mean()),
                "rt_range": float(rt_values.max() - rt_values.min()),
                "mean_ri": float(ri_values.mean()),
                "ri_range": float(ri_values.max() - ri_values.min()),
                **{column: pd.NA for column in area_columns},
                "area_mean": area_mean,
                "area_std": area_std,
                "area_cv_percent": area_cv,
                "quality_flag_samples": ";".join(sorted(rescued)),
            }
            for item in all_members:
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
    ordered = frame.sort_values(
        "sample_count",
        ascending=False,
        kind="stable",
        na_position="last",
    ).reset_index(drop=True)
    area_columns = [
        column for column in ordered.columns if column.startswith(AREA_PREFIX)
    ]
    columns = [
        "mean_rt",
        "canonical_name",
        "mean_ri",
        "sample_count",
        "area_mean",
        *area_columns,
        "detected_samples",
        "quality_flag_samples",
    ]
    return ordered[
        [column for column in columns if column in ordered.columns]
    ].copy()


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
