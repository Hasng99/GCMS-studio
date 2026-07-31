import json

import pandas as pd

from src.settings_bundle import settings_from_json, settings_to_json_bytes, validate_profile


STANDARDS = pd.DataFrame({
    "carbon_number": [6, 7, 8],
    "alkane_name": ["Hexane", "Heptane", "Octane"],
    "ri": [600, 700, 800],
    "rt_min": [2.0, 4.0, 8.0],
})

PROFILE = pd.DataFrame({
    "canonical_name": ["Hexanal", "Octanal"],
    "parent_fatty_acid": ["Linoleic acid", "Oleic acid"],
    "note": [None, "confirmed"],
})


def test_settings_bundle_round_trip() -> None:
    payload = settings_to_json_bytes(
        STANDARDS,
        PROFILE,
        quality_threshold=85,
        fuzzy_matching=True,
    )
    standards, profile, analysis = settings_from_json(payload)

    assert standards["alkane_name"].tolist() == ["Hexane", "Heptane", "Octane"]
    assert profile["canonical_name"].tolist() == ["Hexanal", "Octanal"]
    assert analysis == {"quality_threshold": 85.0, "fuzzy_matching": True}
    assert json.loads(payload)["volatile_profile"][0]["note"] is None


def test_profile_requires_expected_columns() -> None:
    try:
        validate_profile(pd.DataFrame({"canonical_name": ["Hexanal"]}))
    except ValueError as exc:
        assert "parent_fatty_acid" in str(exc)
    else:
        raise AssertionError("missing profile column should fail")


def test_settings_bundle_rejects_unknown_schema() -> None:
    payload = settings_to_json_bytes(
        STANDARDS,
        PROFILE,
        quality_threshold=80,
        fuzzy_matching=False,
    )
    decoded = json.loads(payload)
    decoded["schema_version"] = 99

    try:
        settings_from_json(json.dumps(decoded))
    except ValueError as exc:
        assert "버전" in str(exc)
    else:
        raise AssertionError("unknown schema version should fail")
