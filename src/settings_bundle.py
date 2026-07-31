from __future__ import annotations

from datetime import datetime, timezone
import json
from typing import Any

import pandas as pd

from src.ri import validate_standards


SETTINGS_SCHEMA_VERSION = 1
REQUIRED_PROFILE_COLUMNS = {"canonical_name", "parent_fatty_acid"}


def validate_profile(profile: pd.DataFrame) -> pd.DataFrame:
    missing = REQUIRED_PROFILE_COLUMNS - set(profile.columns)
    if missing:
        raise ValueError("프로필 필수 열이 없습니다: " + ", ".join(sorted(missing)))
    clean = profile.copy()
    clean["canonical_name"] = clean["canonical_name"].fillna("").astype(str).str.strip()
    clean["parent_fatty_acid"] = clean["parent_fatty_acid"].fillna("").astype(str).str.strip()
    clean = clean[clean["canonical_name"] != ""].reset_index(drop=True)
    if clean.empty:
        raise ValueError("프로필에는 물질 이름이 1개 이상 필요합니다.")
    return clean


def _records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    serializable = frame.astype(object).where(pd.notna(frame), None)
    return serializable.to_dict(orient="records")


def settings_to_json_bytes(
    standards: pd.DataFrame,
    profile: pd.DataFrame,
    *,
    quality_threshold: float,
    fuzzy_matching: bool,
) -> bytes:
    clean_standards = validate_standards(standards)
    clean_profile = validate_profile(profile)
    payload = {
        "schema_version": SETTINGS_SCHEMA_VERSION,
        "app": "GC-MS Studio",
        "saved_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": {
            "quality_threshold": float(quality_threshold),
            "fuzzy_matching": bool(fuzzy_matching),
        },
        "standard_rt": _records(clean_standards),
        "volatile_profile": _records(clean_profile),
    }
    return json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False).encode("utf-8")


def settings_from_json(source: bytes | str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    raw = source.decode("utf-8-sig") if isinstance(source, bytes) else source
    payload = json.loads(raw)
    if payload.get("schema_version") != SETTINGS_SCHEMA_VERSION:
        raise ValueError("지원하지 않는 설정 파일 버전입니다.")
    standards = validate_standards(pd.DataFrame(payload.get("standard_rt", [])))
    profile = validate_profile(pd.DataFrame(payload.get("volatile_profile", [])))
    analysis = payload.get("analysis", {})
    threshold = float(analysis.get("quality_threshold", 80.0))
    if not 0 <= threshold <= 100:
        raise ValueError("Quality 기준은 0~100 사이여야 합니다.")
    return standards, profile, {
        "quality_threshold": threshold,
        "fuzzy_matching": bool(analysis.get("fuzzy_matching", False)),
    }
