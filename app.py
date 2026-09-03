from __future__ import annotations

from html import escape
from pathlib import Path
import hashlib

import pandas as pd
import streamlit as st
import yaml

from src.exporters import (
    dataframe_to_csv_bytes,
    multi_results_to_xlsx_bytes,
    results_to_xlsx_bytes,
    selected_hits_to_xlsx_bytes,
)
from src.nist_links import nist_gc_url
from src.matching import normalize_name
from src.multi_sample import (
    AREA_PREFIX,
    DEFAULT_RI_TOLERANCE,
    DETECTED_PREFIX,
    RI_PREFIX,
    RT_PREFIX,
    build_replicate_area_comparison,
    build_sample_presence_comparison,
    replicate_area_view,
    unique_sample_labels,
)
from src.parsers import parse_masshunter
from src.pipeline import PipelineResult, run_pipeline
from src.ri import validate_standards
from src.settings_bundle import settings_from_json, settings_to_json_bytes, validate_profile
from src.standard_selection import apply_selected_candidate_rts


BASE_DIR = Path(__file__).parent
APP_NAME = "GC-MS Studio"
EXCLUDE_SILOXANE_DEFAULT = True
RESULT_COLUMN_LABELS = {
    "sample_name": "Sample info.",
    "compound_number": "No.",
    "rt_min": "RT(min)",
    "hit_number": "Hit No.",
    "hit_name_original": "Hit name",
    "canonical_name": "Compound name",
    "cas_number": "CAS No.",
    "quality": "Quality",
    "profile_match": "Profile match",
    "parent_fatty_acid": "Parent FAs",
    "inclusion_reason": "Supporting",
    "lower_alkane": "Lower alkane",
    "upper_alkane": "Upper alkane",
    "lower_rt": "Lower RT",
    "upper_rt": "Upper RT",
    "ri": "RI",
    "ri_status": "RI status",
    "nist_gc_url": "NIST url",
    "area": "Area",
    "selected_for_peak_summary": "Selected for peak summary",
}


def apply_theme() -> None:
    st.markdown("""
    <style>
    .stApp{background:#f5f7fb;color:#172033}
    [data-testid="stHeader"]{background:rgba(245,247,251,.9)}
    [data-testid="stStatusWidget"]{display:none!important}
    [data-testid="stSpinner"]{
        position:fixed!important;inset:0!important;z-index:999999!important;
        display:flex!important;align-items:center!important;justify-content:center!important;
        padding:0!important;background:rgba(245,247,251,.78);backdrop-filter:blur(2px)
    }
    [data-testid="stSpinner"]>div{
        width:auto!important;padding:1rem 1.35rem!important;border:1px solid #d7dfec;
        border-radius:14px;background:#fff;box-shadow:0 12px 34px rgba(15,23,42,.16)
    }
    [data-testid="stSidebar"]{background:linear-gradient(180deg,#10233f,#22365f)}
    [data-testid="stSidebar"] h1,[data-testid="stSidebar"] h2,
    [data-testid="stSidebar"] h3,[data-testid="stSidebar"] label,
    [data-testid="stSidebar"] p,[data-testid="stSidebar"] span{color:#f8fafc}
    [data-testid="stSidebar"] input{color:#172033!important;background:#fff!important}
    [data-testid="stSidebar"] [data-testid="stFileUploaderDropzone"]{background:#fff;border:1px solid #d7dfec}
    [data-testid="stSidebar"] [data-testid="stFileUploaderDropzone"] *{color:#24324a!important}
    [data-testid="stSidebar"] [data-testid="stFileUploaderFile"]{background:#fff}
    [data-testid="stSidebar"] [data-testid="stFileUploaderFile"] *{color:#24324a!important}
    [data-testid="stSidebar"] button{border-color:#cbd5e1}
    .st-key-analysis_start button{
        color:#fff!important;border:0!important;font-weight:700!important;
        background:linear-gradient(120deg,#15b8a6,#2563eb 55%,#7c3aed)!important;
        box-shadow:0 8px 20px rgba(37,99,235,.28)!important;
        transition:transform .16s ease,box-shadow .16s ease,filter .16s ease!important
    }
    .st-key-analysis_start button:hover{
        transform:translateY(-1px);filter:brightness(1.06);
        box-shadow:0 11px 24px rgba(37,99,235,.36)!important
    }
    .st-key-analysis_start button:disabled{
        transform:none;filter:saturate(.25);box-shadow:none!important;opacity:.55
    }
    .block-container{max-width:1500px;padding-top:1.5rem}
    .hero{padding:2rem 2.2rem;border-radius:22px;color:white;background:linear-gradient(120deg,#10233f,#126e82);box-shadow:0 14px 38px rgba(16,35,63,.18);margin-bottom:1.25rem}
    .hero h1{margin:0 0 .5rem;font-size:2rem;color:white}.hero p{margin:0;color:#d9f3f5;font-size:1.02rem}
    [data-testid="stMetric"]{background:#fff;border:1px solid #e5eaf1;padding:1rem;border-radius:16px;box-shadow:0 5px 18px rgba(15,23,42,.05)}
    .stTabs [data-baseweb="tab-list"]{gap:.3rem}.stTabs [data-baseweb="tab"]{background:#fff;border-radius:10px 10px 0 0;padding:.5rem 1rem}
    .sample-info{display:flex;align-items:center;justify-content:flex-end;gap:.55rem;margin:.85rem .25rem -.15rem;color:#475569}
    .sample-info span{font-size:.78rem;font-weight:600;text-transform:uppercase;letter-spacing:.04em}
    .sample-info strong{color:#172033;font-size:.95rem}
    .status-card{background:#fff;border:1px solid #e5eaf1;border-radius:14px;padding:1rem 1.2rem;margin:.35rem 0}
    .nist-card{background:linear-gradient(135deg,#eef8f8,#fff);border:1px solid #c8e7e8;border-radius:16px;padding:1.2rem 1.4rem;margin:.5rem 0 1rem}
    </style>""", unsafe_allow_html=True)


@st.cache_data
def load_defaults() -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    with (BASE_DIR / "app_config.yaml").open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    return (
        config,
        pd.read_csv(BASE_DIR / "data" / "alkane_standard_rt.csv"),
        pd.read_csv(BASE_DIR / "data" / "fatty_acid_volatile_profile.csv"),
    )


def initialize_session(config: dict, default_standards: pd.DataFrame, default_profile: pd.DataFrame) -> None:
    if "active_standards" not in st.session_state:
        st.session_state["active_standards"] = validate_standards(default_standards)
    if "active_profile" not in st.session_state:
        st.session_state["active_profile"] = validate_profile(default_profile)
    st.session_state.setdefault("quality_threshold", float(config["filter"]["quality_threshold"]))
    st.session_state.setdefault("fuzzy_matching", False)
    st.session_state.setdefault("exclude_siloxane", EXCLUDE_SILOXANE_DEFAULT)


def active_reference_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    return st.session_state["active_standards"].copy(), st.session_state["active_profile"].copy()


def read_confirmed_standard(upload: object) -> pd.DataFrame:
    upload.seek(0)
    suffix = Path(upload.name).suffix.lower()
    if suffix == ".csv":
        return validate_standards(pd.read_csv(upload))
    if suffix == ".xlsx":
        return validate_standards(pd.read_excel(upload, engine="openpyxl"))
    raise ValueError("확정 Standard RT는 CSV 또는 XLSX 파일을 사용하세요.")


def read_profile(upload: object, default_profile: pd.DataFrame) -> pd.DataFrame:
    upload.seek(0)
    suffix = Path(upload.name).suffix.lower()
    if suffix == ".csv":
        return validate_profile(pd.read_csv(upload))
    if suffix == ".pdf":
        payload = upload.read()
        expected_sha256 = "a6251642a3f934be45a9d95963a2975bf1ea7e03891923f090f17148af41a516"
        if hashlib.sha256(payload).hexdigest() != expected_sha256:
            raise ValueError("이 PDF는 검증된 지방산 휘발성분 원자료와 다릅니다. 정규화 CSV를 업로드하세요.")
        return validate_profile(default_profile)
    raise ValueError("휘발성분 프로필은 CSV 또는 검증된 PDF를 사용하세요.")


def replace_active_settings(
    standards: pd.DataFrame,
    profile: pd.DataFrame,
    *,
    quality_threshold: float | None = None,
    fuzzy_matching: bool | None = None,
) -> None:
    st.session_state["active_standards"] = validate_standards(standards)
    st.session_state["active_profile"] = validate_profile(profile)
    if quality_threshold is not None:
        st.session_state["quality_threshold"] = float(quality_threshold)
        st.session_state["analysis_quality_widget"] = float(quality_threshold)
    if fuzzy_matching is not None:
        st.session_state["fuzzy_matching"] = bool(fuzzy_matching)
        st.session_state["analysis_fuzzy_widget"] = bool(fuzzy_matching)
    for key in ("standard_editor", "profile_editor", "raw_candidate_editor"):
        st.session_state.pop(key, None)


def raw_standard_candidates(upload: object, standards: pd.DataFrame) -> pd.DataFrame:
    upload.seek(0)
    hits, _ = parse_masshunter(upload)
    target_names = {
        normalize_name(name): int(carbon)
        for name, carbon in zip(standards["alkane_name"], standards["carbon_number"])
    }
    candidates = hits[hits["hit_name"].map(normalize_name).isin(target_names)].copy()
    if candidates.empty:
        return candidates
    candidates["carbon_number"] = candidates["hit_name"].map(
        lambda value: target_names[normalize_name(value)]
    )
    confirmed = standards[["carbon_number", "rt_min"]].rename(columns={"rt_min": "confirmed_rt"})
    candidates = candidates.merge(confirmed, on="carbon_number", how="left")
    candidates["rt_difference"] = (candidates["rt_min"] - candidates["confirmed_rt"]).abs()
    columns = [
        "carbon_number", "hit_name", "rt_min", "quality", "hit_number",
        "compound_number", "confirmed_rt", "rt_difference",
    ]
    return candidates[columns].sort_values(
        ["carbon_number", "rt_difference", "quality"],
        ascending=[True, True, False],
    )


def settings_page(config: dict, default_standards: pd.DataFrame, default_profile: pd.DataFrame) -> None:
    st.subheader("기준 데이터 설정")
    st.caption("표를 직접 편집하거나 파일을 불러올 수 있습니다. 적용한 설정은 분석 화면에 즉시 반영됩니다.")

    st.markdown("### 전체 설정 저장·복원")
    bundle_left, bundle_middle, bundle_right = st.columns([1.4, .8, .8])
    with bundle_left:
        bundle_upload = st.file_uploader(
            "저장된 설정 파일 불러오기",
            type=["json"],
            key="settings_bundle_upload",
            help="이 앱에서 내려받은 gcms_ri_settings.json 파일을 선택하세요.",
        )
    with bundle_middle:
        st.write("")
        st.write("")
        if st.button("설정 파일 적용", use_container_width=True, disabled=bundle_upload is None):
            try:
                bundle_upload.seek(0)
                standards, profile, analysis = settings_from_json(bundle_upload.read())
                replace_active_settings(
                    standards,
                    profile,
                    quality_threshold=analysis["quality_threshold"],
                    fuzzy_matching=analysis["fuzzy_matching"],
                )
                st.rerun()
            except Exception as exc:
                st.error(f"설정 파일을 적용하지 못했습니다: {exc}")
    standards, profile = active_reference_data()
    with bundle_right:
        st.write("")
        st.write("")
        st.download_button(
            "현재 설정 파일 저장",
            settings_to_json_bytes(
                standards,
                profile,
                quality_threshold=st.session_state["quality_threshold"],
                fuzzy_matching=st.session_state["fuzzy_matching"],
            ),
            "gcms_ri_settings.json",
            "application/json",
            use_container_width=True,
        )

    reset_col, status_col = st.columns([.8, 2.2])
    with reset_col:
        if st.button("내장 기본값으로 초기화", use_container_width=True):
            replace_active_settings(
                default_standards,
                default_profile,
                quality_threshold=float(config["filter"]["quality_threshold"]),
                fuzzy_matching=False,
            )
            st.rerun()
    with status_col:
        st.info(
            f"현재 적용값: Standard RT {len(standards)}행 · 프로필 {len(profile)}행 · "
            f"Quality {st.session_state['quality_threshold']:g}"
        )

    with st.expander("분석 기본 설정 수기 변경"):
        settings_quality = st.number_input(
            "Quality 기준",
            min_value=0.0,
            max_value=100.0,
            value=float(st.session_state["quality_threshold"]),
            step=1.0,
            key="settings_quality_widget",
        )
        settings_fuzzy = st.checkbox(
            "유사 이름 매칭 사용",
            value=bool(st.session_state["fuzzy_matching"]),
            key="settings_fuzzy_widget",
        )
        if st.button("분석 설정 적용", use_container_width=True):
            replace_active_settings(
                standards,
                profile,
                quality_threshold=settings_quality,
                fuzzy_matching=settings_fuzzy,
            )
            st.rerun()

    standard_tab, profile_tab = st.tabs([
        f"Standard RT ({len(standards)}행)",
        f"휘발성분 프로필 ({len(profile)}행)",
    ])
    try:
        with standard_tab:
            st.markdown("#### 파일에서 가져오기")
            standard_upload = st.file_uploader(
                "Standard RT 파일",
                type=["csv", "xlsx", "xls"],
                key="standard_rt_upload",
                help="확정표에는 carbon_number, alkane_name, ri, rt_min 열이 필요합니다.",
            )
            if standard_upload is not None and standard_upload.name.lower().endswith((".csv", ".xlsx")):
                if st.button("업로드한 Standard RT 적용", use_container_width=True):
                    new_standards = read_confirmed_standard(standard_upload)
                    replace_active_settings(new_standards, profile)
                    st.rerun()
            if standard_upload is not None and standard_upload.name.lower().endswith(".xls"):
                st.warning(
                    "같은 물질의 RT 후보가 여러 개일 수 있습니다. 사용할 후보의 '사용' 칸을 체크한 뒤 "
                    "'선택한 RT로 Standard 교체'를 누르세요. 물질별로 하나만 선택할 수 있습니다."
                )
                candidates = raw_standard_candidates(standard_upload, standards)
                if candidates.empty:
                    st.info("현재 Standard 물질과 일치하는 RT 후보를 찾지 못했습니다.")
                else:
                    selectable = candidates.copy()
                    selectable.insert(0, "selected", False)
                    selected_candidates = st.data_editor(
                        selectable,
                        use_container_width=True,
                        hide_index=True,
                        key="raw_candidate_editor",
                        disabled=[column for column in selectable.columns if column != "selected"],
                        column_config={
                            "selected": st.column_config.CheckboxColumn(
                                "사용", help="물질별로 RT 하나만 선택하세요."
                            ),
                            "carbon_number": st.column_config.NumberColumn("탄소 수", format="C%d"),
                            "rt_min": st.column_config.NumberColumn("후보 RT", format="%.4f"),
                            "confirmed_rt": st.column_config.NumberColumn("기존 RT", format="%.4f"),
                            "rt_difference": st.column_config.NumberColumn("차이", format="%.4f"),
                        },
                    )
                    if st.button("선택한 RT로 Standard 교체", type="primary", use_container_width=True):
                        updated_standards = apply_selected_candidate_rts(standards, selected_candidates)
                        replace_active_settings(updated_standards, profile)
                        st.rerun()

            st.markdown("#### 표에서 직접 수정")
            edited_standards = st.data_editor(
                standards,
                num_rows="dynamic",
                use_container_width=True,
                hide_index=True,
                key="standard_editor",
                column_config={
                    "carbon_number": st.column_config.NumberColumn("탄소 수", step=1),
                    "alkane_name": st.column_config.TextColumn("알케인 이름"),
                    "ri": st.column_config.NumberColumn("RI"),
                    "rt_min": st.column_config.NumberColumn("RT (min)", format="%.4f"),
                },
            )
            if st.button("수기 변경 Standard RT 적용", type="primary", use_container_width=True):
                replace_active_settings(edited_standards, profile)
                st.rerun()

        with profile_tab:
            st.markdown("#### 파일에서 가져오기")
            profile_upload = st.file_uploader(
                "휘발성분 프로필 파일",
                type=["csv", "pdf"],
                key="profile_upload",
                help="CSV에는 canonical_name, parent_fatty_acid 열이 필요합니다.",
            )
            if st.button("업로드한 프로필 적용", use_container_width=True, disabled=profile_upload is None):
                new_profile = read_profile(profile_upload, default_profile)
                replace_active_settings(standards, new_profile)
                st.rerun()

            st.markdown("#### 표에서 직접 수정")
            edited_profile = st.data_editor(
                profile,
                num_rows="dynamic",
                use_container_width=True,
                hide_index=True,
                key="profile_editor",
            )
            if st.button("수기 변경 프로필 적용", type="primary", use_container_width=True):
                replace_active_settings(standards, edited_profile)
                st.rerun()
    except Exception as exc:
        st.error(f"설정을 적용하지 못했습니다: {exc}")


def nist_search_tab(result: PipelineResult) -> None:
    st.markdown("""<div class="nist-card"><b>NIST Chemistry WebBook GC/RI 비교</b><br>
    계산된 RI와 NIST의 컬럼 종류, 고정상(active phase), 온도 프로그램이 비슷한 문헌값을 비교하세요.
    실험 조건이 다르면 RI가 달라질 수 있으므로 숫자만 단독 비교하지 않는 것이 중요합니다.</div>""", unsafe_allow_html=True)
    selected = result.selected_hits.drop_duplicates(subset=["canonical_name", "cas_number"]).copy()
    names = selected["canonical_name"].dropna().astype(str).tolist()
    if names:
        compound = st.selectbox("분석 결과에서 물질 선택", names)
        row = selected[selected["canonical_name"].astype(str) == compound].iloc[0]
        c1, c2, c3 = st.columns([1, 1, 1.2])
        c1.metric("내 계산 RI", "-" if pd.isna(row["ri"]) else row["ri"])
        c2.metric("RT (min)", row["rt_min"])
        with c3:
            st.write("NIST GC 데이터")
            st.link_button("선택 물질의 NIST RI 열기 ↗", row["nist_gc_url"], use_container_width=True)
    st.divider()
    manual_name = st.text_input("다른 물질 이름으로 검색", placeholder="예: Hexanal, 2-Octenal")
    left, right = st.columns(2)
    with left:
        if manual_name.strip():
            st.link_button("이름으로 NIST GC 검색 ↗", nist_gc_url("", manual_name.strip()), use_container_width=True)
        else:
            st.button("이름으로 NIST GC 검색 ↗", disabled=True, use_container_width=True)
    with right:
        st.link_button("NIST 상세 이름 검색 페이지 ↗", "https://webbook.nist.gov/chemistry/name-ser/", use_container_width=True)
    st.caption("NIST 페이지에서는 Gas Chromatography 표에서 temperature ramp/isothermal, active phase, 길이·내경·막 두께를 현재 실험 조건과 함께 비교하세요.")
    nist_result = selected[
        ["canonical_name", "cas_number", "rt_min", "ri", "ri_status", "nist_gc_url"]
    ]
    st.dataframe(
        nist_result,
        use_container_width=True,
        hide_index=True,
        column_config=result_column_config(nist_result),
    )


def result_column_config(frame: pd.DataFrame) -> dict[str, object]:
    """분석 결과의 내부 열 이름을 사용자용 표 머리글로 표시한다."""
    column_config: dict[str, object] = {
        column: label
        for column, label in RESULT_COLUMN_LABELS.items()
        if column in frame.columns
    }
    if "nist_gc_url" in frame.columns:
        column_config["nist_gc_url"] = st.column_config.LinkColumn(
            RESULT_COLUMN_LABELS["nist_gc_url"],
            display_text="열기 ↗",
        )
    return column_config


def summary_view(frame: pd.DataFrame) -> pd.DataFrame:
    """요약 화면에서 분석 우선순위에 맞춘 열만 반환한다."""
    priority = [
        "rt_min", "canonical_name", "quality", "ri", "nist_gc_url",
        "area", "profile_match", "inclusion_reason", "parent_fatty_acid",
    ]
    return frame[[column for column in priority if column in frame.columns]].copy()


def result_table_view(frame: pd.DataFrame) -> pd.DataFrame:
    """시료명은 표마다 반복하지 않고 별도의 단일 정보로 표시한다."""
    return frame.drop(columns=["sample_name"], errors="ignore").copy()


def analysis_input_signature(
    uploads: list[object],
    *,
    threshold: float,
    fuzzy: bool,
    exclude_siloxane: bool,
    relationship: str,
    standards: pd.DataFrame,
    profile: pd.DataFrame,
) -> str:
    """파일과 분석 설정이 마지막 실행 이후 바뀌었는지 확인한다."""
    digest = hashlib.sha256()
    for upload in uploads:
        payload = upload.getvalue()
        digest.update(str(upload.name).encode("utf-8"))
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(hashlib.sha256(payload).digest())
    digest.update(
        repr((
            float(threshold),
            bool(fuzzy),
            bool(exclude_siloxane),
            relationship,
        )).encode("utf-8")
    )
    for frame in (standards, profile):
        digest.update(repr(tuple(frame.columns)).encode("utf-8"))
        digest.update(pd.util.hash_pandas_object(frame, index=True).values.tobytes())
    return digest.hexdigest()


def comparison_column_config(frame: pd.DataFrame) -> dict[str, object]:
    config: dict[str, object] = {
        "canonical_name": "Compound name",
        "sample_count": "Replicates",
        "detected_samples": "Detected samples",
        "mean_rt": "Mean RT(min)",
        "rt_range": "RT range",
        "mean_ri": "Mean RI",
        "ri_range": "RI range",
        "area_mean": "Area mean",
        "area_std": "Area SD",
        "area_cv_percent": "Area CV(%)",
        "detection_status": "검출 구분",
        "detected_count": "검출 샘플 수",
        "quality_flag_samples": "Quality 유의 시료",
    }
    for column in frame.columns:
        if column.startswith(AREA_PREFIX):
            config[column] = f"Area · {column.removeprefix(AREA_PREFIX)}"
        elif column.startswith(DETECTED_PREFIX):
            config[column] = column.removeprefix(DETECTED_PREFIX)
        elif column.startswith(RT_PREFIX):
            config[column] = f"RT(min) · {column.removeprefix(RT_PREFIX)}"
        elif column.startswith(RI_PREFIX):
            config[column] = f"RI · {column.removeprefix(RI_PREFIX)}"
    return {column: label for column, label in config.items() if column in frame.columns}


QUALITY_FLAG_STYLE = "background-color:#fee2e2;color:#b91c1c;font-weight:600"


def style_replicate_comparison(frame: pd.DataFrame):
    """quality_flag_samples에 기록된 샘플의 Area 칸을 빨간색으로 강조한다.

    st.dataframe의 캔버스 그리드는 Styler의 배경색을 반영하지 않아 실제 DOM으로
    렌더링되는 st.table과 함께 사용한다. 열 이름은 표시용으로 먼저 바꾸고,
    강조 대상은 원래 열 순서(위치)로 판단해 이름 변경과 무관하게 동작하도록 한다.
    """
    columns = list(frame.columns)
    area_positions = {
        position: column.removeprefix(AREA_PREFIX)
        for position, column in enumerate(columns)
        if column.startswith(AREA_PREFIX)
    }
    flag_position = columns.index("quality_flag_samples") if "quality_flag_samples" in columns else None

    def highlight(row: pd.Series) -> list[str]:
        raw_flag = row.iloc[flag_position] if flag_position is not None else ""
        flagged = {name for name in str(raw_flag or "").split(";") if name}
        styles = [""] * len(row)
        for position, sample in area_positions.items():
            if sample in flagged:
                styles[position] = QUALITY_FLAG_STYLE
        return styles

    labels = comparison_column_config(frame)
    numeric_formats = {
        "mean_rt": "{:.3f}",
        "mean_ri": "{:.1f}",
        "area_mean": "{:.1f}",
        **{column: "{:.1f}" for column in columns if column.startswith(AREA_PREFIX)},
    }
    display = frame.rename(columns=labels)
    return (
        display.style
        .apply(highlight, axis=1)
        .format({labels.get(column, column): fmt for column, fmt in numeric_formats.items()}, na_rep="—")
        .hide(axis="index")
    )


def render_selected_hits_download(
    label: str,
    results_by_sample: dict[str, pd.DataFrame],
    file_stem: str,
    widget_key: str,
) -> None:
    """Selected hits 결과를 CSV 또는 XLSX 형식을 선택해 다운로드한다.

    XLSX는 샘플별로 시트를 나누고 시트 안에서는 sample_name 열을 생략한다.
    CSV는 하나의 파일로 합치므로 sample_name 열을 유지한다.
    """
    file_format = st.radio(
        f"{label} 파일 형식",
        ["CSV", "XLSX"],
        horizontal=True,
        key=f"{widget_key}_format",
    )
    if file_format == "CSV":
        combined = pd.concat(list(results_by_sample.values()), ignore_index=True)
        st.download_button(
            f"{label} 다운로드",
            dataframe_to_csv_bytes(combined),
            f"{file_stem}.csv",
            "text/csv",
            use_container_width=True,
            key=f"{widget_key}_download",
        )
    else:
        st.download_button(
            f"{label} 다운로드",
            selected_hits_to_xlsx_bytes(results_by_sample),
            f"{file_stem}.xlsx",
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            use_container_width=True,
            key=f"{widget_key}_download",
        )


def render_metrics(result: PipelineResult) -> None:
    labels = {
        "total_peaks": "전체 peak", "total_hits": "전체 hits", "profile_match": "Profile 일치",
        "quality_pass": "Quality 통과", "both": "BOTH", "ri_ok": "RI 성공", "out_of_range": "범위 밖",
    }
    columns = st.columns(len(labels))
    for column, (key, label) in zip(columns, labels.items()):
        column.metric(label, result.metrics[key])


def render_sample_info(sample_name: str) -> None:
    st.markdown(
        f'<div class="sample-info"><span>Sample name</span><strong>{escape(sample_name)}</strong></div>',
        unsafe_allow_html=True,
    )


def render_single_result(
    result: PipelineResult,
    sample_name: str,
    standards: pd.DataFrame,
    profile: pd.DataFrame,
) -> None:
    render_metrics(result)
    render_sample_info(sample_name)
    tabs = st.tabs([
        "요약", "Peak summary", "Selected hits", "Rejected hits",
        "NIST RI 검색", "Standards", "Profile",
    ])
    frames = [
        summary_view(result.peak_summary),
        result_table_view(result.peak_summary),
        result_table_view(result.selected_hits),
        result_table_view(result.rejected_hits),
    ]
    for tab, frame in zip(tabs[:4], frames):
        with tab:
            st.dataframe(
                frame,
                use_container_width=True,
                hide_index=True,
                column_config=result_column_config(frame),
            )
    with tabs[4]:
        nist_search_tab(result)
    with tabs[5]:
        st.dataframe(standards, use_container_width=True, hide_index=True)
    with tabs[6]:
        st.dataframe(profile, use_container_width=True, hide_index=True)
    left, right = st.columns(2)
    with left:
        render_selected_hits_download(
            "Selected hits",
            {sample_name: result.selected_hits},
            "selected_hits",
            "single_selected_hits",
        )
    right.download_button(
        "전체 결과 XLSX 다운로드",
        results_to_xlsx_bytes(result, standards, profile),
        "gcms_ri_results.xlsx",
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        use_container_width=True,
    )


def render_multi_results(
    results: dict[str, PipelineResult],
    relationship: str,
    standards: pd.DataFrame,
    profile: pd.DataFrame,
    *,
    ri_tolerance: float = DEFAULT_RI_TOLERANCE,
) -> None:
    summary_frames = {
        sample_name: result.peak_summary
        for sample_name, result in results.items()
    }
    rejected_frames = {
        sample_name: result.rejected_hits
        for sample_name, result in results.items()
    }
    same_sample = relationship == "동일 샘플 내 반복시료"
    if same_sample:
        comparison = replicate_area_view(
            build_replicate_area_comparison(
                summary_frames, rejected_frames, ri_tolerance=ri_tolerance,
            )
        )
        comparison_label = "Area 비교"
        comparison_sheet = "area_comparison"
    else:
        comparison = build_sample_presence_comparison(summary_frames)
        comparison_label = "샘플 비교"
        comparison_sheet = "sample_comparison"

    tabs = st.tabs([
        comparison_label,
        *[f"{sample_name} 요약" for sample_name in results],
        "Standards",
        "Profile",
    ])
    with tabs[0]:
        if same_sample:
            st.caption(
                f"같은 Compound name 중 RT 차이가 ±0.05분 이내이고 RI 차이가 ±{ri_tolerance:g} 이내인 "
                "조합을 동일 물질로 봅니다. 이 범위를 만족하는 군집이 여러 개면 생략하지 않고 모두 "
                "별도 행으로 표시합니다. 한 군집이 전체 반복 수의 70% 이상에서 검출되면, 나머지 "
                "샘플에서 Quality가 낮아 제외됐더라도 같은 물질로 함께 제시하며 해당 Area 칸을 "
                "빨간색으로 표시해 유의가 필요함을 나타냅니다."
            )
        else:
            st.caption(
                "Compound name을 기준으로 전체 샘플의 공통 검출, 일부 샘플의 부분 공통, "
                "한 샘플에서만 나온 개별 검출 물질을 구분합니다."
            )
            status_counts = comparison["detection_status"].value_counts()
            count_columns = st.columns(3)
            for column, status in zip(
                count_columns,
                ["공통 검출", "부분 공통", "개별 검출"],
            ):
                column.metric(status, int(status_counts.get(status, 0)))
        if comparison.empty:
            st.info("비교 조건을 만족하는 물질이 없습니다.")
        elif same_sample:
            st.table(style_replicate_comparison(comparison))
        else:
            st.dataframe(
                comparison,
                use_container_width=True,
                hide_index=True,
                column_config=comparison_column_config(comparison),
            )

    sample_tabs = tabs[1:1 + len(results)]
    for tab, (sample_name, result) in zip(sample_tabs, results.items()):
        with tab:
            render_metrics(result)
            render_sample_info(sample_name)
            frame = summary_view(result.peak_summary)
            st.dataframe(
                frame,
                use_container_width=True,
                hide_index=True,
                column_config=result_column_config(frame),
            )
    with tabs[-2]:
        st.dataframe(standards, use_container_width=True, hide_index=True)
    with tabs[-1]:
        st.dataframe(profile, use_container_width=True, hide_index=True)

    left, right = st.columns(2)
    with left:
        render_selected_hits_download(
            "전체 샘플 Selected hits",
            {sample_name: result.selected_hits for sample_name, result in results.items()},
            "multi_selected_hits",
            "multi_selected_hits",
        )
    right.download_button(
        "전체 샘플 결과 XLSX 다운로드",
        multi_results_to_xlsx_bytes(
            results,
            standards,
            profile,
            comparison=comparison,
            comparison_sheet=comparison_sheet,
        ),
        "gcms_multi_results.xlsx",
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        use_container_width=True,
    )


def analysis_page(config: dict, standards: pd.DataFrame, profile: pd.DataFrame) -> None:
    with st.sidebar:
        st.markdown("### 분석 입력")
        threshold = st.number_input(
            "Quality 기준",
            min_value=0.0,
            max_value=100.0,
            value=float(st.session_state["quality_threshold"]),
            step=1.0,
            key="analysis_quality_widget",
        )
        fuzzy = st.checkbox(
            "유사 이름 매칭 사용",
            value=bool(st.session_state["fuzzy_matching"]),
            key="analysis_fuzzy_widget",
            help="오탐 가능성이 있어 기본값은 꺼짐입니다.",
        )
        exclude_siloxane = st.checkbox(
            "siloxane 계열 제외",
            value=bool(st.session_state["exclude_siloxane"]),
            key="analysis_siloxane_widget",
            help="Compound name에 'siloxane', 'siloxyl' 또는 'siloxy'가 포함된 물질을 추천 결과에서 제외합니다.",
        )
        st.session_state["quality_threshold"] = float(threshold)
        st.session_state["fuzzy_matching"] = bool(fuzzy)
        st.session_state["exclude_siloxane"] = bool(exclude_siloxane)
        sample_uploads = st.file_uploader(
            "시료 결과 (.xls/.xlsx/.csv, 복수 선택 가능)",
            type=["xls", "xlsx", "csv"],
            key="sample_upload",
            accept_multiple_files=True,
        )
        relationship = "단일 시료"
        replicate_ri_tolerance = float(
            st.session_state.get("replicate_ri_tolerance", DEFAULT_RI_TOLERANCE)
        )
        if len(sample_uploads) > 1:
            relationship = st.radio(
                "업로드 파일 관계",
                ["동일 샘플 내 반복시료", "서로 다른 샘플"],
                key="multi_file_relationship",
            )
            if relationship == "동일 샘플 내 반복시료":
                replicate_ri_tolerance = st.number_input(
                    "반복시료 동일 물질 판정 RI 허용범위 (±)",
                    min_value=1.0,
                    max_value=200.0,
                    value=replicate_ri_tolerance,
                    step=1.0,
                    key="replicate_ri_tolerance_widget",
                    help="같은 이름의 물질을 동일 물질로 볼 때 허용하는 RI 차이입니다.",
                )
                st.session_state["replicate_ri_tolerance"] = float(replicate_ri_tolerance)
        start_clicked = st.button(
            "Start",
            key="analysis_start",
            type="primary",
            use_container_width=True,
            disabled=not sample_uploads,
        )
        st.caption(f"Standard RT: {len(standards)}행 · 프로필: {len(profile)}행")
    if not sample_uploads:
        st.session_state.pop("analysis_results", None)
        st.session_state.pop("analysis_signature", None)
        st.session_state.pop("analysis_relationship", None)
        st.info("왼쪽에서 MassHunter 시료 결과 파일을 하나 이상 업로드하세요.")
        st.markdown("**처리 흐름:** 업로드 → Start → 후보 선별 → RI 계산 → NIST 비교 → CSV/XLSX 다운로드")
        st.dataframe(standards[["alkane_name", "ri", "rt_min"]], use_container_width=True, hide_index=True)
        return

    input_signature = analysis_input_signature(
        sample_uploads,
        threshold=threshold,
        fuzzy=fuzzy,
        exclude_siloxane=exclude_siloxane,
        relationship=relationship,
        standards=standards,
        profile=profile,
    )
    if start_clicked:
        try:
            with st.spinner("loading …"):
                parsed_uploads: list[tuple[pd.DataFrame, str]] = []
                for upload in sample_uploads:
                    try:
                        hits, metadata = parse_masshunter(upload)
                    except Exception as exc:
                        raise ValueError(f"{upload.name}: {exc}") from exc
                    parsed_uploads.append((
                        hits,
                        str(metadata.get("sample_name") or Path(upload.name).stem),
                    ))
                sample_labels = unique_sample_labels([
                    sample_name for _, sample_name in parsed_uploads
                ])
                results: dict[str, PipelineResult] = {}
                for (hits, _), sample_name in zip(parsed_uploads, sample_labels):
                    results[sample_name] = run_pipeline(
                        hits,
                        profile,
                        standards,
                        sample_name=sample_name,
                        quality_threshold=threshold,
                        fuzzy=fuzzy,
                        exclude_siloxane=exclude_siloxane,
                        allow_extrapolation=bool(config["ri"]["allow_extrapolation"]),
                        round_digits=int(config["ri"]["round_digits"]),
                        exact_tolerance=float(config["ri"]["exact_standard_tolerance_min"]),
                    )
            st.session_state["analysis_results"] = results
            st.session_state["analysis_signature"] = input_signature
            st.session_state["analysis_relationship"] = relationship
        except Exception as exc:
            st.session_state.pop("analysis_results", None)
            st.session_state.pop("analysis_signature", None)
            st.session_state.pop("analysis_relationship", None)
            st.error(f"분석을 완료하지 못했습니다: {exc}")
            return

    if st.session_state.get("analysis_signature") != input_signature:
        st.info("파일과 분석 옵션을 확인한 뒤 왼쪽의 Start 버튼을 눌러주세요.")
        return
    results = st.session_state.get("analysis_results")
    if not results:
        st.info("왼쪽의 Start 버튼을 눌러 분석을 시작하세요.")
        return
    relationship = str(
        st.session_state.get("analysis_relationship", relationship)
    )
    if len(results) == 1:
        sample_name, result = next(iter(results.items()))
        render_single_result(result, sample_name, standards, profile)
    else:
        render_multi_results(
            results, relationship, standards, profile, ri_tolerance=replicate_ri_tolerance,
        )


def main() -> None:
    st.set_page_config(page_title=APP_NAME, page_icon="🧪", layout="wide")
    apply_theme()
    config, default_standards, default_profile = load_defaults()
    initialize_session(config, default_standards, default_profile)
    st.markdown(f"""<div class="hero"><h1>{APP_NAME}</h1><p>MassHunter hit 선별 · 지방산 산화 프로필 매칭 · RI 계산 · NIST 조건 비교</p></div>""", unsafe_allow_html=True)
    with st.sidebar:
        st.markdown(f"## {APP_NAME}")
        page = st.radio("화면 선택", ["분석", "기준 설정"], horizontal=True)
        st.divider()
    if page == "기준 설정":
        settings_page(config, default_standards, default_profile)
        return
    standards, profile = active_reference_data()
    analysis_page(config, standards, profile)


if __name__ == "__main__":
    main()
