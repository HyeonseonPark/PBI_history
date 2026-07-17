# Bulk mRNA-seq 차등발현 분석 (Panax_ginseng_123w_development_mrRNA-seq, Panax_ginseng_Heatstress_mrRNA-seq, Panax_ginseng_Light_mrRNA-seq, Zoysia_sinica_japonica_mrRNAseq)

## 1. 공통 파이프라인

```
[Raw FASTQ (paired-end)]
        │  QC/Trim  (프로젝트별 상이 — 2번 섹션 참고: Trimmomatic+PRINSEQ++ / fastp / Prinseq-lite)
        ▼
[Clean FASTQ]
        │  HISAT2 정렬 (2.2.1, 전 프로젝트 공통 확인)
        ▼
[Sorted BAM(.bai)]
        │  StringTie 정량 (3.0.0, -eB ballgown 모드, 전 프로젝트 공통 확인)
        ▼
[샘플별 transcripts.gtf, gene_abund.tab, *.ctab]
        │  prepDE 계열 스크립트(count matrix) + TPM 추출 스크립트(도구는 프로젝트별 상이)
        ▼
[gene/transcript count matrix, TPM matrix]
        │  edgeR TMM 정규화 (calcNormFactors) — 확인된 3/4 프로젝트에서 수행 (Panax_ginseng_Light_mrRNA-seq 제외, 2번 섹션 참고)
        ▼
[TMM 정규화 matrix]
        │  Trinity 스타일 edgeR 차등발현분석 (pairwise exactTest) — 확인된 3/4 프로젝트에서 수행 (Panax_ginseng_Light_mrRNA-seq 제외)
        ▼
[edgeR DE_results, 클러스터링/히트맵 등 최종 산출물]
```

| 스텝 | 소프트웨어(버전) | 핵심 파라미터 | Input | Output |
|---|---|---|---|---|
| QC/Trim | 프로젝트별 상이(아래 2번 섹션 참고) | 프로젝트별 상이 | Raw FASTQ R1/R2 | Clean FASTQ R1/R2 |
| HISAT2 정렬 | HISAT2 2.2.1 (전 프로젝트 공통 확인) | `--rna-strandness RF` (공통 확인), 스레드 수(`-p`)는 프로젝트별 상이 | Clean FASTQ + HISAT2 index | Sorted BAM(.bai) |
| StringTie 정량 | StringTie 3.0.0 (전 프로젝트 공통 확인) | `-eB --rf -G <참조 GTF/GFF3>` (ballgown 모드) | Sorted BAM + 참조 GTF/GFF3 | 샘플별 `transcripts.gtf`, `gene_abund.tab`, `*.ctab` |
| Count/TPM matrix 생성 | prepDE.py3 계열(공통, 세부 변형은 프로젝트별 상이) + TPM 추출 스크립트(TPM_extraction.py 또는 stringtie_expression_matrix.pl, 프로젝트별 상이) | 기본값 사용(확인 필요 항목 다수) | 샘플별 ctab/gene_abund.tab | `gene_count_matrix.csv`, `transcript_count_matrix.csv`, TPM matrix |
| TMM 정규화 | R edgeR `calcNormFactors` — 3/4 프로젝트 확인(Panax_ginseng_Light_mrRNA-seq은 관련 스크립트 자체가 프로젝트 내에 없음) | 기본 TMM 옵션 | count/TPM matrix | `TMM.EXPR.matrix`, `*.TMM_info.txt` |
| edgeR 차등발현분석 | Trinity util(run_DE_analysis.pl 계열) + R edgeR `exactTest` — 3/4 프로젝트 확인(Panax_ginseng_Light_mrRNA-seq은 DEG 분석 디렉토리가 비어 있어 미수행으로 확인됨) | `cpm(matrix)>1 in ≥2 samples` 필터, pairwise exactTest; P-value/log2FC cutoff는 프로젝트별 상이(아래 2번 섹션 참고) | count matrix + samples.file | `*.edgeR.DE_results`, `*.subset`, MA/Volcano plot, 클러스터링/히트맵 |

## 2. 프로젝트별 특이사항

### Panax_ginseng_123w_development_mrRNA-seq
- 연구목적/샘플 설계: 인삼(Panax ginseng)의 1년생 1·2·3주차 뿌리 이차 생장 발달 분자 메커니즘 연구, 9샘플(1wr/2wr/3wr × 3반복)
- 공통 파이프라인과의 차이점:
  - QC/Trim이 2단계(Trimmomatic 0.40으로 어댑터 제거 후 PRINSEQ++ 1.2로 품질필터/dedup)로 구성됨 — 다른 세 프로젝트와 상이한 QC 툴 조합
  - edgeR DE 분석은 1w/2w/3w_vs_REST + Default 조합의 pairwise 비교, P0.001/log2FC cutoff=2(C2) 기준
  - Trinity 고정높이 트리컷 클러스터링(P=55, 6개 클러스터) 수행
  - TMM 발현행렬에 대해 Z-score 행렬 생성(Z-scoring_TPM_matrix.py, log2(TMM+1) 후 유전자별 z-score) 후, 최상위 후처리 스크립트(filter_matrix.py, filter_deg_matrix_by_go.py)로 DEG 목록 및 GO term 7개 기준 세분화 matrix까지 추가 생성
- (해당 시) 프로젝트 고유 확장 분석: `functional_annotation/`의 Diamond 기반 인삼-애기장대(Arabidopsis) reciprocal best hit 오솔로그 매핑 → 이를 이용한 GO/KEGG enrichment(원본은 외부 DAVID 웹툴 분석 결과를 수동 입력한 것으로 확인됨) → 버블플롯 생성 → GSEA 분석(clusterProfiler, 오솔로그 매핑 기반 Custom TERM2GENE 구성) — 논문 투고 후 리비전 단계(1w/2w/3w 개별)에서 재작업되어 최종 보충표까지 작성됨. 단, `Gene_Regulatory_Network/`(HOMER 모티프·TF 조절망)와 `reference/split_subgenome.py`(서브게놈 비교) 분석은 최종 논문에 사용되지 않은 것으로 확인되어 제외됨.

### Panax_ginseng_Heatstress_mrRNA-seq
- 연구목적/샘플 설계: 인삼 CP/JW 품종을 24도(대조군)와 35도(1일/3일 열스트레스) 처리 후 잎 전사체 연구, 18샘플(CP_Con/CP_1d/CP_3d/JW_Con/JW_1d/JW_3d × 3반복)
- 공통 파이프라인과의 차이점:
  - QC/Trim은 fastp 1.0.1 단일 스텝(`--cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --dedup`) + MultiQC 통합 리포트로 구성 — 별도 PRINSEQ 단계 없음
  - 데이터 정합성 문제 확인됨: 실제 15개 pairwise edgeR DE 분석에 사용된 count matrix(`gene_count_matrix.tsv`)는 prepDE.py3 버그가 있는 원본에서 파생되어 84개 유전자/isoform이 누락된 상태로 산출된 것으로 확인됨. 버그 수정본(`gene_count_matrix_fixed.csv`)은 생성되었으나 어떤 하류 분석에서도 참조되지 않음(연구자 확인, gene_count_matrix.csv는 prepDE.py로 분석된 것으로 최종 확정 — 감사 기록에는 유지)
  - `matrix_preprocess.py` → `genes.counts.matrix` 경로는 산출물 mtime이 입력파일보다 16시간 이상 앞서는 타임스탬프 역행이 확인되어 폐기된/미완주 분기로 판단됨
  - edgeR DE 분석은 6그룹 전체 15개 pairwise 조합, P0.001/C2 cutoff, 클러스터링 P=30 기준
  - 다른 세 프로젝트와 달리 Z-score 행렬 생성 단계는 확인되지 않음
- (해당 시) 프로젝트 고유 확장 분석: `reference/diamond.sh` 기반 인삼-애기장대 오솔로그 매핑(reciprocal best hit)이 존재하나, DEG 분석(2026-04-07 완료) 완료 8일 후(2026-04-14) 실행되었고 파일 상호참조가 전혀 없어 DEG 파이프라인과 완전히 분리된 별개 산출물로 확인됨(연구자 확인: 최종 인수인계 리포트 범위에서 제외하기로 결정, 감사 기록에는 유지).

### Panax_ginseng_Light_mrRNA-seq
- 연구목적/샘플 설계: 인삼 CP/JW 품종의 광질 처리(Control/White/Blue/Red/Far-Red) 후 잎 전사체 연구로 추정, Pg1~Pg15가 Control/White/Blue/Red/Far-Red 순서로 3반복씩 배정됨(작성자 확인으로 확정)
- 공통 파이프라인과의 차이점:
  - QC/Trim은 fastp 1.0.1 단일 스텝(`cleaning/` 디렉토리) + MultiQC 통합이나, 이를 감싼 wrapper 스크립트 파일 자체가 저장소에 없어 `fastp.log`로만 커맨드가 재구성됨
  - **TMM 정규화 및 edgeR 차등발현분석 단계가 확인되지 않음**: `DEG_anaylsis/` 디렉토리가 비어 있고, 프로젝트 전체에서 `runTMM.R`이나 edgeR 관련 스크립트가 전혀 발견되지 않음 — count/TPM matrix 생성 이후 정규화·DE 분석이 수행되지 않았거나 본 저장소 범위 밖에서 수행된 것으로 추정됨(확인 필요)
  - `matrix_preprocess.py`가 의도한 출력(`genes.counts.matrix`)도 실제로는 존재하지 않고, 대신 `gene_count_matrix.tsv`가 생성되어 있으나 이 파일의 mtime이 스크립트 최종 수정시각보다 앞서는 시간 모순이 있어 생성 경위가 불명확함(확인 필요)
- (해당 시) 프로젝트 고유 확장 분석: 없음 — DE 분석 자체가 미완료로 확인되어 하류 확장 분석도 존재하지 않음

### Zoysia_sinica_japonica_mrRNAseq
- 연구목적/샘플 설계: 갯벌잔디(Z. sinica)와 육지잔디(Z. japonica)의 NaCl 처리 1일 후 뿌리 mRNA-seq을 통한 종간 전사체 비교, 조건(Control/Salt)당 3반복 × 2종
- 공통 파이프라인과의 차이점:
  - QC/Trim은 Prinseq-lite 0.20.4 단일 스텝(별도 Trimmomatic/fastp 없음)
  - 공통 백본(Track A)이 참조 게놈 3개 버전(Zj pseudomolecule / Zj_contig / Zs)에 대해 동일 스크립트 세트로 병렬 반복 실행됨 — 다른 세 프로젝트는 단일 참조 게놈만 사용
  - TPM matrix 생성에 Perl 스크립트(`stringtie_expression_matrix.pl`, griffithlab 계열)를 사용 — 다른 세 프로젝트는 Python `TPM_extraction.py` 사용
  - edgeR DE 분석 cutoff는 P0.01/C1(+ 일부 P0.001/C2도 병행), Trinity 버전이 v2.15.2로 로그에서 직접 확인됨(다른 프로젝트는 Trinity util로 추정되나 123w에서만 동일 버전 명시 확인됨)
  - Zs(Z. sinica)에 한해서만 Trinity 고정높이 서브클러스터링(P=30) 및 클러스터별 fpkm matrix 생성 수행(Zj/Zj_contig는 서브클러스터링 미실시로 확인)
  - Zj/Zs 각각에 대해 Z-score 행렬 생성(Z-scoring_TPM_matrix.py, log2(TPM+1) 후 유전자별 z-score) — Zj_contig는 스크립트만 존재하고 실제 실행/output 확인 안 됨(확인 필요)
- (해당 시) 프로젝트 고유 확장 분석: 이 프로젝트는 Track A(공통 RNA-seq 1차 처리) 외에 4개의 독립 확장 트랙(Track B~E)을 보유함
  - Track B (TF_Prediction): Track A의 Zs 서브클러스터1 결과를 기반으로 TSS 좌표 추출(±2000/100bp) 후 bed→fasta 추출(사용 툴 미확인) 및 TF 후보 유전자군 트리(OG0006158) 구성 — 다수 단계의 생성 스크립트/선정 기준이 확인 필요로 남아 있음
  - Track C (Zoysia_Interspecies_Transcriptomic_analysis): Track A의 Zj·Zs 단백질 서열을 reciprocal BLASTP로 매칭 → 매칭 서브셋 count/TPM matrix 병합(Matrix_merge.py) → TMM 정규화 → 그룹별(Control/Salt_Treatment/Whole) Trinity/edgeR 종간 비교 DE 분석(Whole은 4그룹 전체 6개 쌍별 비교 + 서브클러스터링) → Z-score 행렬 생성
  - Track D (Salt_transporters, 12sp/12sp_Zj_COM): 외부 12종 orthogroup 단백질 세트에 대해 InterProScan(v5.74-105.0) 도메인 스캔 → 종별/전체종 도메인 presence·count 통계 산출 — 본 프로젝트 rawdata와 직접 연결되지 않는 2차 비교분석으로 확인됨
  - Track E (Zoysia_ENTAP): Track A 유래 단백질 서열(ZJN pseudochromosome, ZJN contig)과 외부 3종(Oryza_sativa, ZMW=Z. matrella, ZPZ=Z. pacifica) 단백질에 대해 EnTAP v2.2.0(DIAMOND+EggNOG-mapper) 기능 어노테이션 수행 → 일부 종(ZJN_pseudochromosome)만 TBtools 포맷팅 추가 수행 → DE 결과와 ENTAP 주석 결합(결합 스크립트는 저장소에 없어 수동 join으로 추정)

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Panax_ginseng_123w_development_mrRNA-seq | Trimmomatic+PRINSEQ++ 2단계 QC; TMM z-score 행렬을 DEG 목록·GO term 기준으로 추가 세분화하는 후처리 단계 보유 | Trimmomatic 0.40, PRINSEQ++ 1.2 (타 프로젝트와 다른 QC 툴 조합) | Diamond 오솔로그 매핑 + DAVID 기반 GO/KEGG enrichment + GSEA(clusterProfiler), 논문 리비전 재작업 포함 |
| Panax_ginseng_Heatstress_mrRNA-seq | fastp 단일 QC; 실제 DE 분석에 쓰인 count matrix가 버그 있는 prepDE 원본(84개 유전자/isoform 누락)으로 확인됨 | fastp 1.0.1(단일 스텝, PRINSEQ 미사용) | 인삼-애기장대 오솔로그 매핑(diamond, RBH) — DEG 파이프라인과 파일시스템상 완전히 분리된 별개 산출물 |
| Panax_ginseng_Light_mrRNA-seq | fastp 단일 QC; **DE 분석(TMM 정규화·edgeR) 자체가 미수행**으로 확인됨(DEG_anaylsis 디렉토리 비어 있음) | fastp 1.0.1 | 없음(DE 분석 미완료로 하류 확장분석 부재) |
| Zoysia_sinica_japonica_mrRNAseq | Prinseq-lite 단일 QC; 참조 게놈 3버전(Zj/Zj_contig/Zs) 병렬 처리; TPM 추출에 Perl 스크립트 사용 | prinseq-lite 0.20.4, Trinity v2.15.2(로그로 직접 확인), InterProScan 5.74-105.0, EnTAP v2.2.0 | TF_Prediction(Track B), 종간비교분석(Track C), Salt_transporters 도메인분석(Track D), Zoysia_ENTAP 기능어노테이션(Track E) |
