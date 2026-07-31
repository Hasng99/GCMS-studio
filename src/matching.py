from __future__ import annotations

import re
import unicodedata

import pandas as pd
from rapidfuzz import fuzz, process


def normalize_name(value: object) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).casefold().strip()
    text = text.replace("–", "-").replace("—", "-")
    text = re.sub(r"[()\[\]{}]", "", text)
    return re.sub(r"[\W_]+", "", text, flags=re.UNICODE)


def normalize_cas(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    digits = re.sub(r"\D", "", str(value))
    return digits.lstrip("0")


def build_profile_index(profile: pd.DataFrame) -> tuple[dict[str, set[int]], dict[str, set[int]]]:
    names: dict[str, set[int]] = {}
    cases: dict[str, set[int]] = {}
    for index, row in profile.reset_index(drop=True).iterrows():
        candidates = [row.get("canonical_name"), row.get("source_name")]
        aliases = row.get("aliases_semicolon")
        if pd.notna(aliases):
            candidates.extend(str(aliases).split(";"))
        for value in candidates:
            key = normalize_name(value)
            if key:
                names.setdefault(key, set()).add(index)
        cas = normalize_cas(row.get("cas_number"))
        if cas:
            cases.setdefault(cas, set()).add(index)
    return names, cases


def match_profile(
    hit_name: object,
    cas_number: object,
    profile: pd.DataFrame,
    fuzzy: bool = False,
    fuzzy_score: int = 95,
) -> dict[str, object]:
    names, cases = build_profile_index(profile)
    matched: set[int] = set()
    cas = normalize_cas(cas_number)
    if cas and cas in cases:
        matched |= cases[cas]
    name_key = normalize_name(hit_name)
    if not matched and name_key in names:
        matched |= names[name_key]
    if not matched and fuzzy and name_key:
        result = process.extractOne(name_key, names.keys(), scorer=fuzz.ratio, score_cutoff=fuzzy_score)
        if result:
            matched |= names[result[0]]
    if not matched:
        return {"profile_match": False, "canonical_name": str(hit_name), "parent_fatty_acid": ""}
    rows = profile.iloc[sorted(matched)]
    canonical = rows["canonical_name"].dropna().astype(str).iloc[0]
    parents = "; ".join(sorted(set(rows["parent_fatty_acid"].dropna().astype(str))))
    return {"profile_match": True, "canonical_name": canonical, "parent_fatty_acid": parents}
