from __future__ import annotations

import io
import re
from pathlib import Path
from typing import BinaryIO

import pandas as pd


COLUMN_ALIASES = {
    "compoundnumber": "compound_number",
    "compoundnumber#": "compound_number",
    "rtmin": "rt_min",
    "scannumber": "scan_number",
    "scannumber#": "scan_number",
    "areaabs": "area",
    "baselineheigthab": "baseline_height",
    "baselineheightab": "baseline_height",
    "absoluteheigthab": "absolute_height",
    "absoluteheightab": "absolute_height",
    "peakwidth50min": "peak_width_50_min",
    "hitnumber": "hit_number",
    "hitname": "hit_name",
    "quality": "quality",
    "molweightamu": "mol_weight",
    "casnumber": "cas_number",
    "library": "library",
    "entrynumberlibrary": "entry_number_library",
}

PEAK_COLUMNS = [
    "compound_number",
    "rt_min",
    "scan_number",
    "area",
    "baseline_height",
    "absolute_height",
    "peak_width_50_min",
]
REQUIRED_COLUMNS = {
    "compound_number",
    "rt_min",
    "hit_number",
    "hit_name",
    "quality",
    "cas_number",
    "area",
}
NUMERIC_COLUMNS = [
    "compound_number",
    "rt_min",
    "scan_number",
    "area",
    "baseline_height",
    "absolute_height",
    "peak_width_50_min",
    "hit_number",
    "quality",
    "mol_weight",
    "entry_number_library",
]


def _column_key(value: object) -> str:
    text = str(value).strip().lower()
    return re.sub(r"[\s()_%*]+", "", text)


def _read_bytes(source: str | Path | bytes | BinaryIO) -> tuple[bytes | None, str]:
    if isinstance(source, (str, Path)):
        path = Path(source)
        return None, path.name
    if isinstance(source, bytes):
        return source, "upload.xls"
    name = getattr(source, "name", "upload.xls")
    source.seek(0)
    return source.read(), str(name)


def _excel_source(source: str | Path | bytes | BinaryIO) -> tuple[object, str]:
    payload, name = _read_bytes(source)
    return (io.BytesIO(payload) if payload is not None else source), name


def _find_header(raw: pd.DataFrame) -> int:
    for index, row in raw.head(30).iterrows():
        keys = {_column_key(value) for value in row if pd.notna(value)}
        if "hitnumber" in keys and "hitname" in keys and "quality" in keys:
            return int(index)
    raise ValueError("MassHunter 헤더 행(Hit Number, Hit Name, Quality)을 찾지 못했습니다.")


def _metadata_from_rows(raw: pd.DataFrame) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for value in raw.iloc[:8, 0].dropna().astype(str):
        if ":" in value:
            key, item = value.split(":", 1)
            metadata[key.strip().lower().replace(" ", "_")] = item.strip()
    return metadata


def _normalize_frame(frame: pd.DataFrame) -> pd.DataFrame:
    rename: dict[object, str] = {}
    for column in frame.columns:
        key = _column_key(column)
        rename[column] = COLUMN_ALIASES.get(key, str(column).strip())
    frame = frame.rename(columns=rename)
    missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError("필수 열이 없습니다: " + ", ".join(missing))
    frame = frame.dropna(how="all").copy()
    frame[PEAK_COLUMNS] = frame[PEAK_COLUMNS].ffill()
    for column in NUMERIC_COLUMNS:
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["hit_name"] = frame["hit_name"].astype("string").str.strip()
    frame["cas_number"] = frame["cas_number"].astype("string").str.strip()
    frame = frame[frame["hit_number"].notna() & frame["hit_name"].notna()].reset_index(drop=True)
    for column in ("compound_number", "hit_number"):
        frame[column] = frame[column].astype("Int64")
    return frame


def parse_masshunter(
    source: str | Path | bytes | BinaryIO,
    sheet_name: str = "LibRes",
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Parse MassHunter LibRes from xls/xlsx/csv and forward-fill peak fields."""
    payload, name = _read_bytes(source)
    suffix = Path(name).suffix.lower()
    if suffix == ".csv":
        csv_source: object = io.BytesIO(payload) if payload is not None else source
        frame = pd.read_csv(csv_source, dtype={"cas_number": "string", "CAS Number": "string"})
        return _normalize_frame(frame), {"sample_name": Path(name).stem}
    if suffix not in {".xls", ".xlsx"}:
        raise ValueError("지원 형식은 .xls, .xlsx, .csv입니다.")
    engine = "xlrd" if suffix == ".xls" else "openpyxl"
    excel_source: object = io.BytesIO(payload) if payload is not None else source
    raw = pd.read_excel(excel_source, sheet_name=sheet_name, header=None, engine=engine, dtype=object)
    header_index = _find_header(raw)
    metadata = _metadata_from_rows(raw)
    columns = raw.iloc[header_index].tolist()
    frame = raw.iloc[header_index + 1 :].copy()
    frame.columns = columns
    normalized = _normalize_frame(frame)
    metadata["sample_name"] = metadata.get("sample_name", Path(name).stem)
    metadata["header_row"] = str(header_index + 1)
    return normalized, metadata
