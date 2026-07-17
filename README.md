# 생물정보학 분석 파이프라인 종합 정리

15개 프로젝트(acc1 10개 + acc2 4개 + NAS_Archive 1개)의 실측 기반 파이프라인 문서를 분석
목표/방법 기준 5개 카테고리로 분류하고, 카테고리 내 공통 파이프라인과 프로젝트별 특이사항을
비교 정리한 문서 세트입니다. 저장 위치가 아니라 분석 방법 기준으로 분류했습니다.

## 분류 체계

| 카테고리 | 프로젝트 수 | 설명 | 문서 |
|---|---|---|---|
| De novo genome assembly | 5 | 장리드(PacBio HiFi/ONT) 시퀀싱 기반 de novo 전장유전체 조립, scaffolding, QC, 유전자/반복서열 주석 | [de_novo_genome_assembly.md](de_novo_genome_assembly.md) |
| Resequencing 변이탐지 DNA-seq | 2 | 참조 기반 리시퀀싱(bwa-mem2/bwa+GATK)을 통한 QTL 후보 유전자 탐색 또는 분자마커 개발 | [resequencing_dnaseq.md](resequencing_dnaseq.md) |
| Bulk mRNA-seq 차등발현 | 4 | HISAT2+StringTie 정량 → TMM 정규화 → Trinity/edgeR 차등발현분석 표준 파이프라인 | [bulk_mRNAseq_DE.md](bulk_mRNAseq_DE.md) |
| Multi-omics/holobiome | 1 | 뿌리 전사체(mRNA-seq) + 근권 메타지놈(16S) 통합 분석(MOFA+) | [multi_omics_holobiome.md](multi_omics_holobiome.md) |
| Single-cell RNA-seq | 2 | 10X Genomics cellranger 기반 단일세포 전사체 분석(Seurat/R 또는 Scanpy/Python) | [single_cell_rnaseq.md](single_cell_rnaseq.md) |

## 전체 프로젝트 한눈에 보기

| 프로젝트 | 저장위치 | 카테고리 | 연구목표(한줄) | 핵심 소프트웨어 스택 | 프로젝트 고유 특징 |
|---|---|---|---|---|---|
| Cannabis_sativa_cv.ChungSam_dnDNAseq | acc1 | De novo genome assembly | 한국 재래종 대마(청삼) 자웅이주 개체의 haplotype-resolved de novo 전장유전체 조립 및 카나비노이드 합성 유전자 분석 | hifiasm, RagTag, BUSCO, QUAST, RepeatModeler/RepeatMasker, BRAKER3, Sourmash, chromeister/syri/plotsr | Oatk로 organelle 분리 후 hifiasm 2회 조립(1차 폐기 추정); Female/Male 각 hap1/hap2, 총 4 haplotype; 판게놈(98샘플) 비교 및 THCAS/CBDAS/CBGAS 유전자 분석 |
| Pichia_Pastoris_dnDNAseq_HanHwa | acc1 | De novo genome assembly | 한화솔루션 Pichia 균주(BG-10, X-33) 2종의 ONT+Illumina 하이브리드 de novo genome assembly | NextDenovo, NextPolish, Ratatosk, RagTag, BUSCO, QUAST, Merqury, mummer, LiftOff | trycycler multi-assembler 합의조립 시도 후 reconcile 실패로 폐기, NextDenovo+NextPolish 단일 경로 채택; 유전자 주석은 de novo 예측이 아닌 LiftOff 기반 liftover(GS115/CBS7435 참조) |
| Abies_koreana_dnDNAseq | acc2 | De novo genome assembly | 구상나무 표준목(제주도)의 전장 유전체 해독 | hifiasm, Oatk, minimap2, seqkit, QUAST | 2026-07 기준 진행 중 프로젝트 — assembly + 사후 organelle 필터링 + QUAST QC까지만 완료, scaffolding/비교유전체/유전자주석은 아직 미수행; BUSCO 미사용 |
| Zoysia_macrostachya_dnDNAseq | acc2 | De novo genome assembly | Zoysia macrostachya de novo 조립 → Omni-C 스캐폴딩 → 갭클로징 → 주석 → 비교유전체(분화시기/유전자가족진화/양성선택) | Ratatosk, NextDenovo/NextPolish/ntLink, HapHiC, Juicebox, YagCloser, BRAKER3+TSEBRA, OrthoFinder, MCMCTree, CAFE5 | Hi-C 기반 de novo 스캐폴딩 사용(HapHiC); pan-genome 비교 + 염스트레스 SRA 전사체 재분석까지 포함한 복합적인 파이프라인 |
| Zoysia_sinica_Genome_Assembly | NAS_Archive | De novo genome assembly | *Zoysia sinica*(갯벌잔디) de novo 전장유전체 조립 및 비교유전체(연구목적은 description.txt 부재로 폴더/파일명 기반 추정) | yahs/SALSA2/3D-DNA(Hi-C 스캐폴딩), RagTag, BUSCO, QUAST, Merqury, RepeatModeler/RepeatMasker, BRAKER3(v3.0.7), AGAT(v1.4.0), OrthoFinder, CAFE, MCMCTree | 원시 de novo assembly 단계가 접근 불가 디렉토리(`Original/`)에서 수행되어 검증 자체가 불가능(5개 프로젝트 중 불확실 항목이 가장 많음); Hi-C 스캐폴딩 툴 3종(yahs/SALSA2/3D-DNA) 병행 시도 후 yahs 채택; 별도 수동 큐레이션(Manual_curation) 단계 보유 |
| Capsella_rubella_rsDNAseq | acc1 | Resequencing 변이탐지 DNA-seq | *Capsella rubella* 1GR1×22.5 F2 100개체의 flowering time/cauline leaf number 형질 QTL 연구 | fastp, bwa-mem2, GATK4, vcftools, R/qtl2 | 두 참조 계통(Cr1GR1/Cr22.5)으로 alignment 단계부터 완전 병렬 수행, 서로 다른 지점(step09/step10)에서 종료; 최종 후보 유전자 리스트는 Cr22.5 계통에서만 기능주석까지 진행되어 산출 |
| Actinia_arguat_rsDNAseq | acc2 | Resequencing 변이탐지 DNA-seq | Actinidia arguta(다래) 12품종 구분용 분자마커(rhAmp/IDT genotyping) 개발 | prinseq-lite, bwa, GATK4, vcftools, bcftools, bedtools, PLINK | 단일 참조 게놈, MarkDuplicates 없이 RG 태깅만 수행; QTL이 아닌 SNP/InDel 조합 기반 최소 마커셋 탐색(greedy set-cover)으로 종결 |
| Panax_ginseng_123w_development_mrRNA-seq | acc1 | Bulk mRNA-seq 차등발현 | 인삼 1년생 1·2·3주차 뿌리 이차 생장 발달 분자 메커니즘 연구 | Trimmomatic, PRINSEQ++, HISAT2, StringTie, edgeR/Trinity, clusterProfiler | 2단계 QC(Trimmomatic+PRINSEQ++, 타 프로젝트와 상이); Diamond 오솔로그 매핑+GO/KEGG+GSEA 후속 분석까지 논문 리비전 단계에서 재작업 완료 |
| Panax_ginseng_Heatstress_mrRNA-seq | acc1 | Bulk mRNA-seq 차등발현 | 인삼 CP/JW 품종의 24도/35도(1일·3일) 열스트레스 잎 전사체 연구 | fastp, HISAT2, StringTie, edgeR/Trinity | 실제 DE 분석에 쓰인 count matrix가 prepDE.py3 버그로 84개 유전자/isoform 누락 상태(연구자 확인, 최종본으로 유지); 오솔로그 매핑은 DEG 파이프라인과 파일상 완전 분리된 별개 산출물 |
| Panax_ginseng_Light_mrRNA-seq | acc1 | Bulk mRNA-seq 차등발현 | 인삼 CP/JW 품종의 광질(Control/White/Blue/Red/Far-Red) 처리 잎 전사체 연구 | fastp, HISAT2, StringTie | **TMM 정규화/edgeR DE 분석 자체가 미수행**으로 확인(DEG_anaylsis 디렉토리 비어 있음) — 4개 프로젝트 중 유일하게 DE 결과 없음 |
| Zoysia_sinica_japonica_mrRNAseq | acc1 | Bulk mRNA-seq 차등발현 | 갯벌잔디(Z. sinica)-육지잔디(Z. japonica) NaCl 처리 뿌리 전사체 종간 비교 | Prinseq-lite, HISAT2, StringTie, Trinity v2.15.2/edgeR, reciprocal BLASTP, InterProScan, EnTAP | 참조 게놈 3버전(Zj/Zj_contig/Zs) 병렬 처리(Track A) + 4개 독립 확장 트랙(TF_Prediction/종간비교/Salt_transporters/ENTAP 주석) 보유 — 카테고리 내 가장 복합적인 DAG |
| Panax_ginseng_holobiome | acc1 | Multi-omics/holobiome | 인삼 2·3·4년생 뿌리 전사체-근권 미생물(16S) 상호관계 규명을 통한 연작장해·진세노사이드 합성 조절 메커니즘 연구 | QIIME2, phyloseq/ANCOM-BC2, HISAT2, StringTie2, edgeR, WGCNA, eggNOG-mapper, clusterProfiler, SpiecEasi/SparCC, MOFA+ | 16S+mRNA-seq 두 축을 MOFA+로 통합(유일한 멀티오믹스 프로젝트); mRNA-seq 정량 경로에서 featureCounts/DESeq2/구버전WGCNA 실패 이력 존재; "양자 영감 ML" 계획은 미착수 상태로 확인 |
| Homo_sapiens_BEST1_scRNAseq | acc2 | Single-cell RNA-seq | Human brain organoid 세포주(#2, #7)의 BEST1 knockdown 단일세포 전사체 프로파일링 | Cell Ranger 10.0(OCM 4-plex), Seurat 5, DoubletFinder, Harmony, clusterProfiler류(fgsea/msigdbr) | 10X OCM 4-plex 라이브러리 1회 cellranger multi 실행으로 4샘플 디멀티플렉싱; Harmony 기반 donor 통합 경로와 Harmony 미사용 donor 분리 경로가 병행 수행되는 유일한 프로젝트 |
| Homo_sapiens_Tcell_scRNAseq | acc1 (`2nd_generation/` 하위) | Single-cell RNA-seq | 인간 제대혈 T세포 anti-CD3/CD28 activation 단일세포 전사체 프로파일링(Control vs Treatment) | Cell Ranger 10.0, Scanpy 1.12.1, AnnData, Scrublet | 동일 라이브러리 3회 재시퀀싱 배치를 심볼릭 링크 lane 위장 트릭으로 병합; Seurat 대신 Scanpy/Python 생태계 사용(카테고리 내 유일); NgR1(RTN4R)/NogoA(RTN4) 발현 기반 T세포 서브셋 DEG 분석; DoubletFinder 교차검증은 사이드 브랜치로만 존재하고 최종 파이프라인엔 미반영 |

## 범위 및 원칙

- 대상은 acc1 10개 + acc2 4개(각자 폴더의 기존 `handover_docs`: 01_directory_map.md +
  02_pipeline_steps.md, 스크립트·로그 실측 교차검증 기반) + NAS_Archive 1개(핸드오버 문서는
  없었으나 이번 작업에서 직접 스크립트·로그를 실측 검증) 총 15개 프로젝트로 한정합니다.
  acc1/acc2 중 handover_docs가 없는 프로젝트(acc1: Arabidopsis_scRNAseq_Chungbuk_Prof.Jo,
  Camelina_sativa_genome, Capsella_rubella_dnDNAseq, Capsicum_annum_snRNAseq,
  Larix_kaempferi_mRNAseq, Panax_ginseng_snRNAseq / acc2: GBM_mRNAseq,
  korean_pine_genome_assembly, Lemna_Kinnex_Isoseq, scRNA-seq_data, Zoysia_Phenomics)는
  이번 범위에서 제외되었습니다.
- Homo_sapiens_Tcell_scRNAseq는 프로젝트 최상위가 아니라 `2nd_generation/` 하위 폴더에
  handover_docs가 있어 최초 스캔에서 누락되었다가 이후 추가되었습니다. 프로젝트 최상위의
  1차 분석 산출물(`rawdata`, `cellranger_output`, `reference`,
  `run_cellranger_analysis.sh` 등)은 해당 handover_docs 자체가 스코프 외로 명시한 부분이라
  이번 종합정리에도 포함하지 않았습니다 — 향후 필요 시 별도 확인 필요.
- Zoysia_sinica_Genome_Assembly(NAS_Archive)는 description.txt/README 자체가 없고 별도
  handover_docs도 작성되지 않은 상태였습니다. 이번 작업에서 실제 스크립트/로그/파일
  타임스탬프를 직접 읽어 검증한 뒤 바로 이 종합정리 문서에 반영했으며, 원본 위치에는
  별도 handover_docs를 남기지 않았습니다. 원시 de novo assembly 단계가 이루어진
  `Original/ONT_longread_DNAseq/` 디렉토리는 조사 시점에 접근 불가(간헐적 NFS 이슈)였어서
  검증하지 못했습니다 — 재확인 필요.
- 원본 문서에서 "확인 필요"로 표시된 불확실 항목은 새로 추정하지 않고 각 카테고리
  문서에 그대로 보존했습니다.
- 각 카테고리 문서의 파이프라인 스텝은 표가 아닌 스텝별 서술형(소프트웨어/파라미터/
  Input/Output + 실제 스크립트에서 발췌한 실행 커맨드)으로 정리했습니다. 실행 커맨드는
  각 프로젝트의 실제 `.sh`/`.py`/`.R`/`.pl` 스크립트를 읽고 발췌한 것이며, 반복문/로깅
  등 보일러플레이트는 제거하고 대표 1건으로 축약했습니다.
- 같은 카테고리에 여러 프로젝트가 속할 경우, 실제로 공유되는 부분만 "공통 파이프라인"
  섹션에 두고, 진행 단계 차이나 방법론 차이로 공유되지 않는 부분은 프로젝트별 섹션
  또는 별도 스텝으로 명확히 구분했습니다(예: de_novo_genome_assembly.md의 Abies_koreana는
  아직 조립 단계까지만 진행된 프로젝트로 별도 표기).
- 이 문서 세트는 자체 완결적으로 작성되어, 원본 `handover_docs` 경로에 대한 참조
  링크를 포함하지 않습니다.
