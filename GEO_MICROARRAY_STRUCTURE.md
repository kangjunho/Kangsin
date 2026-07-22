# GEO Microarray 교육과정 구조

GitHub 저장소의 루트를 기준으로 `public/geo/` 폴더를 생성합니다.

```text
public/
└── geo/
    ├── index.html
    ├── 00_installation.html
    ├── 01_r_basics.html
    ├── 02_geo_search.html
    ├── 03_geo_download.html
    ├── 04_preprocessing.html
    ├── 05_deg.html
    ├── 06_visualization.html
    ├── 07_enrichment.html
    ├── 08_interpretation.html
    ├── assets/
    │   ├── css/
    │   │   └── geo-course.css
    │   ├── images/
    │   └── data/
    └── scripts/
        ├── 00_setup.R
        ├── 03_download_geo.R
        ├── 04_preprocessing.R
        ├── 05_deg_limma.R
        ├── 06_visualization.R
        └── 07_enrichment.R
```

현재 `e-learning.html`의 GEO 카드는 이미 `public/geo/index.html`을 가리키고 있으므로 수정하지 않아도 됩니다.

## 작성 순서

1. `00_installation.html`: Windows용 R/RStudio 설치와 패키지 설정
2. `01_r_basics.html`: 실습에 필요한 최소 R 문법
3. `02_geo_search.html`: GEO 검색과 데이터셋 선정 기준
4. 예제 GSE 확정
5. `03`–`08` 페이지를 동일한 예제와 객체명으로 연결

## 운영 원칙

- 한 과정에서는 하나의 GEO 예제 데이터셋을 끝까지 사용합니다.
- 각 페이지는 `학습목표 → 개념 → 따라하기 → 코드 → 결과 확인 → 문제 해결 → 다음 단계` 순서를 유지합니다.
- 실행 코드는 페이지 안에 제시하고, 전체 코드는 `scripts/`에서 별도로 다운로드할 수 있게 합니다.
- Windows 경로와 한글 사용자명에서 발생하는 오류를 별도 안내합니다.
- 아직 작성하지 않은 강의는 목차에서 `준비 중`으로 표시합니다.
