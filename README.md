# 생물정보학 분석 파이프라인 종합 정리

9개 프로젝트의 실측 기반 파이프라인 문서(`handover_docs`)를 분석 목표/방법 기준 4개
카테고리로 분류하고, 카테고리 내 공통 파이프라인과 프로젝트별 특이사항을 비교 정리한
문서 세트입니다.

## 분류 체계

| 카테고리 | 프로젝트 수 | 설명 | 문서 |
|---|---|---|---|
| De novo genome assembly | 2 | 장리드(PacBio HiFi/ONT) 시퀀싱 기반 de novo 전장유전체 조립, scaffolding, QC, 유전자/반복서열 주석 | [de_novo_genome_assembly.md](de_novo_genome_assembly.md) |
| Resequencing/QTL DNA-seq | 1 | 참조 기반 리시퀀싱(bwa-mem2+GATK) 및 형질 연관 QTL 후보 유전자 탐색 | [resequencing_dnaseq.md](resequencing_dnaseq.md) |
| Bulk mRNA-seq 차등발현 | 4 | HISAT2+StringTie 정량 → TMM 정규화 → Trinity/edgeR 차등발현분석 표준 파이프라인 | [bulk_mRNAseq_DE.md](bulk_mRNAseq_DE.md) |
| Multi-omics/holobiome | 1 | 뿌리 전사체(mRNA-seq) + 근권 메타지놈(16S) 통합 분석(MOFA+) | [multi_omics_holobiome.md](multi_omics_holobiome.md) |

## 전체 프로젝트 한눈에 보기

| 프로젝트 | 카테고리 | 연구목표(한줄) | 핵심 소프트웨어 스택 | 프로젝트 고유 특징 |
|---|---|---|---|---|
| Cannabis_sativa_cv.ChungSam_dnDNAseq | De novo genome assembly | 한국 재래종 대마(청삼) 자웅이주 개체의 haplotype-resolved de novo 전장유전체 조립 및 카나비노이드 합성 유전자 분석 | hifiasm, RagTag, BUSCO, QUAST, RepeatModeler/RepeatMasker, BRAKER3, Sourmash, chromeister/syri/plotsr | Oatk로 organelle 분리 후 hifiasm 2회 조립(1차 폐기 추정); Female/Male 각 hap1/hap2, 총 4 haplotype; 판게놈(98샘플) 비교 및 THCAS/CBDAS/CBGAS 유전자 분석 |
| Pichia_Pastoris_dnDNAseq_HanHwa | De novo genome assembly | 한화솔루션 Pichia 균주(BG-10, X-33) 2종의 ONT+Illumina 하이브리드 de novo genome assembly | NextDenovo, NextPolish, Ratatosk, RagTag, BUSCO, QUAST, Merqury, mummer, LiftOff | trycycler multi-assembler 합의조립 시도 후 reconcile 실패로 폐기, NextDenovo+NextPolish 단일 경로 채택; 유전자 주석은 de novo 예측이 아닌 LiftOff 기반 liftover(GS115/CBS7435 참조) |
| Capsella_rubella_rsDNAseq | Resequencing/QTL DNA-seq | *Capsella rubella* 1GR1×22.5 F2 100개체의 flowering time/cauline leaf number 형질 QTL 연구 | fastp, bwa-mem2, GATK4, vcftools, R/qtl2 | 두 참조 계통(Cr1GR1/Cr22.5)으로 alignment 단계부터 완전 병렬 수행, 서로 다른 지점(step09/step10)에서 종료; 최종 후보 유전자 리스트는 Cr22.5 계통에서만 기능주석까지 진행되어 산출 |
| Panax_ginseng_123w_development_mrRNA-seq | Bulk mRNA-seq 차등발현 | 인삼 1년생 1·2·3주차 뿌리 이차 생장 발달 분자 메커니즘 연구 | Trimmomatic, PRINSEQ++, HISAT2, StringTie, edgeR/Trinity, clusterProfiler | 2단계 QC(Trimmomatic+PRINSEQ++, 타 프로젝트와 상이); Diamond 오솔로그 매핑+GO/KEGG+GSEA 후속 분석까지 논문 리비전 단계에서 재작업 완료 |
| Panax_ginseng_Heatstress_mrRNA-seq | Bulk mRNA-seq 차등발현 | 인삼 CP/JW 품종의 24도/35도(1일·3일) 열스트레스 잎 전사체 연구 | fastp, HISAT2, StringTie, edgeR/Trinity | 실제 DE 분석에 쓰인 count matrix가 prepDE.py3 버그로 84개 유전자/isoform 누락 상태(연구자 확인, 최종본으로 유지); 오솔로그 매핑은 DEG 파이프라인과 파일상 완전 분리된 별개 산출물 |
| Panax_ginseng_Light_mrRNA-seq | Bulk mRNA-seq 차등발현 | 인삼 CP/JW 품종의 광질(Control/White/Blue/Red/Far-Red) 처리 잎 전사체 연구 | fastp, HISAT2, StringTie | **TMM 정규화/edgeR DE 분석 자체가 미수행**으로 확인(DEG_anaylsis 디렉토리 비어 있음) — 4개 프로젝트 중 유일하게 DE 결과 없음 |
| Zoysia_sinica_japonica_mrRNAseq | Bulk mRNA-seq 차등발현 | 갯벌잔디(Z. sinica)-육지잔디(Z. japonica) NaCl 처리 뿌리 전사체 종간 비교 | Prinseq-lite, HISAT2, StringTie, Trinity v2.15.2/edgeR, reciprocal BLASTP, InterProScan, EnTAP | 참조 게놈 3버전(Zj/Zj_contig/Zs) 병렬 처리(Track A) + 4개 독립 확장 트랙(TF_Prediction/종간비교/Salt_transporters/ENTAP 주석) 보유 — 카테고리 내 가장 복합적인 DAG |
| Panax_ginseng_holobiome | Multi-omics/holobiome | 인삼 2·3·4년생 뿌리 전사체-근권 미생물(16S) 상호관계 규명을 통한 연작장해·진세노사이드 합성 조절 메커니즘 연구 | QIIME2, phyloseq/ANCOM-BC2, HISAT2, StringTie2, edgeR, WGCNA, eggNOG-mapper, clusterProfiler, SpiecEasi/SparCC, MOFA+ | 16S+mRNA-seq 두 축을 MOFA+로 통합(유일한 멀티오믹스 프로젝트); mRNA-seq 정량 경로에서 featureCounts/DESeq2/구버전WGCNA 실패 이력 존재; "양자 영감 ML" 계획은 미착수 상태로 확인 |

## 범위 및 원칙

- 대상은 `handover_docs`(01_directory_map.md + 02_pipeline_steps.md, 스크립트·로그 실측
  교차검증 기반)가 이미 작성된 9개 프로젝트로 한정합니다. handover_docs가 없는 7개
  프로젝트(Arabidopsis_scRNAseq_Chungbuk_Prof.Jo, Camelina_sativa_genome,
  Capsella_rubella_dnDNAseq, Capsicum_annum_snRNAseq, Homo_sapiens_Tcell_scRNAseq,
  Larix_kaempferi_mRNAseq, Panax_ginseng_snRNAseq)는 이번 범위에서 제외되었습니다.
- 원본 문서에서 "확인 필요"로 표시된 불확실 항목은 새로 추정하지 않고 각 카테고리
  문서에 그대로 보존했습니다.
- 이 문서 세트는 자체 완결적으로 작성되어, 원본 `handover_docs` 경로에 대한 참조
  링크를 포함하지 않습니다.
