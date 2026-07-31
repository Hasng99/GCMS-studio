import pandas as pd

from src.matching import match_profile
from src.nist_links import nist_gc_url


PROFILE = pd.DataFrame({
    "parent_fatty_acid": ["Linoleic acid"],
    "source_name": ["Hexanal"],
    "canonical_name": ["Hexanal"],
    "aliases_semicolon": [""],
})


def test_exact_profile_name() -> None:
    result = match_profile("(Hexanal)", "", PROFILE)
    assert result["profile_match"] is True
    assert result["canonical_name"] == "Hexanal"


def test_nist_cas_url_removes_padding() -> None:
    assert nist_gc_url("000124-38-9", "Carbon dioxide").endswith("ID=C124389&Mask=2000")


def test_nist_name_fallback() -> None:
    assert "Name=Hexanal" in nist_gc_url("", "Hexanal")
