# GC-MS Studio

HS-SPME/GC-MS로 얻어진 휘발성 지방산화생성물 프로필 또는 MassHunter library 기준으로 선별하고, Standard n-alkane RT로 RI를 판별하는 Streamlit 웹 앱입니다.

## 웹에서 사용

[GC-MS Studio 열기](https://lllstudio.streamlit.app/)

Windows, macOS, 모바일에서 별도 설치 없이 브라우저로 사용할 수 있습니다. 첫 접속은 Streamlit 서버가 깨어나는 데 수십 초가 걸릴 수 있습니다.

## 주요 기능

- MassHunter 추출 `.xls`, `.xlsx`, `.csv` 파일 업로드
- Profile match or Quality threshold로 후보 물질 선별
- Kovats RI 계산
- NIST Chemistry WebBook GC/RI 바로가기
- Standard RT와 휘발성분 프로필의 표 직접 편집
- MassHunter Standard `.xls`의 여러 RT 후보 중 사용할 값 선택
- 현재 Standard RT·프로필·Quality 설정을 단일 JSON 파일로 저장 및 복원
- CSV 및 다중 시트 XLSX 결과 다운로드

## Standard RT 변경

왼쪽 메뉴의 **기준 설정**에서 변경합니다.

### 표에서 직접 변경

1. **Standard RT** 탭의 표에서 값을 수정하거나 행을 추가·삭제합니다.
2. **수기 변경 Standard RT 적용**을 누릅니다.
3. RT와 탄소 수가 순서대로 증가하는지 자동 검증한 뒤 분석에 적용됩니다.

### MassHunter `.xls`에서 RT 선택

1. **Standard RT 파일**에 Standard 측정 결과 `.xls`를 업로드합니다.
2. 같은 물질의 여러 RT 후보 중 사용할 행의 **사용** 체크박스를 선택합니다.
3. **선택한 RT로 Standard 교체**를 누릅니다.

여러 물질을 한 번에 변경할 수 있지만 같은 탄소 수의 후보는 하나만 선택해야 합니다. 선택하지 않은 물질은 기존 RT를 유지합니다.

### 설정 저장 및 다음 접속에서 복원

1. 원하는 RT, 프로필, Quality 기준을 적용합니다.
2. **현재 설정 파일 저장**을 눌러 `gcms_ri_settings.json`을 받습니다.
3. 다음 접속에서는 **저장된 설정 파일 불러오기**에 이 JSON 하나를 올립니다.
4. **설정 파일 적용**을 누르면 모든 설정이 한 번에 복원됩니다.

설정은 사용자의 브라우저 세션 메모리에만 유지됩니다. 다른 사용자나 다음 세션에서 사용하려면 JSON 설정 파일을 저장하십시오.

## 분석 결과 사용

1. 왼쪽 메뉴에서 **분석**을 선택합니다.
2. Quality 기준과 유사 이름 매칭 여부를 확인합니다.
3. MassHunter 시료 결과를 업로드합니다.
4. 요약, 후보, NIST RI 검색 결과를 확인합니다.
5. CSV 또는 전체 XLSX를 다운로드합니다.

요약 탭은 `RT → Compound → Quality → RI → Area` 순서로 표시됩니다.

## 로컬 실행

```bash
python -m pip install -r requirements.txt
streamlit run app.py
```

업로드된 파일은 서버 메모리에서 처리됩니다. 연구 데이터가 민감하다면 공개 앱 대신 접근 제한이 있는 비공개 앱 또는 사내 서버를 사용하십시오.
