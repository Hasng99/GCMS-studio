from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from .matching import match_profile
from .nist_links import nist_gc_url
from .ri import calculate_ri, validate_standards


OUTPUT_COLUMNS = [
    "sample_name", "compound_number", "rt_min", "hit_number", "hit_name_original",
    "canonical_name", "cas_number", "quality", "profile_match", "parent_fatty_acid",
    "inclusion_reason", "lower_alkane", "upper_alkane", "lower_rt", "upper_rt",
    "ri", "ri_status", "nist_gc_url", "area", "selected_for_peak_summary",
]


@dataclass
class PipelineResult:
    all_hits: pd.DataFrame
    selected_hits: pd.DataFrame
    peak_summary: pd.DataFrame
    rejected_hits: pd.DataFrame
    metrics: dict[str, int]


def run_pipeline(
    hits: pd.DataFrame,
    profile: pd.DataFrame,
    standards: pd.DataFrame,
    *,
    sample_name: str = "",
    quality_threshold: float = 80,
    fuzzy: bool = False,
    allow_extrapolation: bool = False,
    round_digits: int = 1,
    exact_tolerance: float = 0.003,
) -> PipelineResult:
    standards = validate_standards(standards)
    rows: list[dict[str, object]] = []
    for _, hit in hits.iterrows():
        match = match_profile(hit.get("hit_name"), hit.get("cas_number"), profile, fuzzy=fuzzy)
        quality = pd.to_numeric(hit.get("quality"), errors="coerce")
        quality_pass = bool(pd.notna(quality) and quality >= quality_threshold)
        profile_match = bool(match["profile_match"])
        included = profile_match or quality_pass
        reason = "BOTH" if profile_match and quality_pass else ("PROFILE" if profile_match else ("QUALITY" if quality_pass else ""))
        ri_result = calculate_ri(
            hit.get("rt_min"), standards, allow_extrapolation=allow_extrapolation,
            round_digits=round_digits, exact_tolerance=exact_tolerance,
        ).to_dict()
        row = {
            "sample_name": sample_name,
            "compound_number": hit.get("compound_number"),
            "rt_min": hit.get("rt_min"),
            "hit_number": hit.get("hit_number"),
            "hit_name_original": hit.get("hit_name"),
            "canonical_name": match["canonical_name"],
            "cas_number": hit.get("cas_number"),
            "quality": quality,
            "profile_match": profile_match,
            "parent_fatty_acid": match["parent_fatty_acid"],
            "inclusion_reason": reason,
            **ri_result,
            "nist_gc_url": nist_gc_url(hit.get("cas_number"), hit.get("hit_name")),
            "area": hit.get("area"),
            "selected_for_peak_summary": False,
            "_included": included,
            "_quality_pass": quality_pass,
        }
        rows.append(row)
    all_hits = pd.DataFrame(rows)
    selected = all_hits[all_hits["_included"]].copy()
    if not selected.empty:
        winners = (
            selected.sort_values(["quality", "hit_number"], ascending=[False, True], na_position="last")
            .groupby(["compound_number", "inclusion_reason"], dropna=False, sort=False)
            .head(1)
            .index
        )
        selected.loc[winners, "selected_for_peak_summary"] = True
    summary = selected[selected["selected_for_peak_summary"]].copy()
    rejected = all_hits[~all_hits["_included"]].copy()
    metrics = {
        "total_peaks": int(hits["compound_number"].nunique()),
        "total_hits": int(len(hits)),
        "profile_match": int(all_hits["profile_match"].sum()),
        "quality_pass": int(all_hits["_quality_pass"].sum()),
        "both": int((all_hits["inclusion_reason"] == "BOTH").sum()),
        "ri_ok": int((selected["ri_status"] == "OK").sum()),
        "out_of_range": int((selected["ri_status"] == "OUT_OF_RANGE").sum()),
    }
    visible = OUTPUT_COLUMNS
    return PipelineResult(
        all_hits=all_hits[visible].reset_index(drop=True),
        selected_hits=selected[visible].reset_index(drop=True),
        peak_summary=summary[visible].reset_index(drop=True),
        rejected_hits=rejected[visible].reset_index(drop=True),
        metrics=metrics,
    )
