# Resequencing/QTL DNA-seq (Capsella_rubella_rsDNAseq)

## 1. 공통 파이프라인

두 개의 참조 계통(`Cr1GR1`, `Cr22.5` ragtag assembly)에 대해 alignment 단계부터 파이프라인이 병렬로 수행된다. 아래 DAG는 공통 상류(rawdata~trimming)와 이후 두 계통 병렬 분기(Cr1GR1 / Cr22.5)를 함께 표시한다.

```
[rawdata] Illumina PE fastq.gz (100 샘플 x R1/R2)
   │ (verify_md5.sh, check_reads.py, parse_summary.py = 사전 QC/검증, 산출물 없음, 비블로킹)
   ▼
[fastq.sh] fastp v1.0.1 트리밍 + dedup
   ├──▶ multiqc ──▶ multiqc_result/
   │
   ├───────────────────────────────┬────────────────────────────────┐
   ▼ (Cr1GR1 경로)                  ▼ (Cr22.5 경로)
[bwamem2-GATK-Cr1GR1.sh]            [bwamem2-GATK.sh]
 REF=Cr1GR1.ragtag.fasta             REF=Cr22.5.ragtag.fasta
 bwa-mem2 mem → samtools sort        bwa-mem2 mem → samtools sort
   ▼                                   ▼
 gatk MarkDuplicates                 gatk MarkDuplicates
   ▼                                   ▼
 gatk HaplotypeCaller -ERC GVCF      gatk HaplotypeCaller -ERC GVCF
   │                                   │
   ├──▶ [get_alignment_stats.sh] (samtools flagstat, 사후 QC, DAG 곁가지) ──▶ Cr1GR1/Cr22.5 flagstat/*.tsv
   ▼                                   ▼
[Merging_gvcf.sh] (Cr1GR1)           [Merging_gvcf.sh] (Cr22.5)
 gatk CombineGVCFs → GenotypeGVCFs   gatk CombineGVCFs → GenotypeGVCFs
   ▼                                   ▼
 SelectVariants(SNP only)            SelectVariants(SNP only)
 (Cr1GR1용 별도 .sh 없음 — Cr22.5     [GATK_filtering.sh] 1단계
  GATK_filtering.sh 코드를 동일
  파라미터로 재사용해 대화형 실행)
   ▼                                   ▼
 VariantFiltration(Basic_Filter)     VariantFiltration(Basic_Filter)
 (동일 재사용 원칙 적용)              [GATK_filtering.sh] 2단계
   ▼                                   ▼
 vcftools 추가 필터링                 vcftools 추가 필터링
 (동일 재사용 원칙 적용)              [GATK_filtering.sh] 3단계
   ▼                                   ▼
[vcf_to_qtl2_1GR1.py]                [vcf_to_qtl2.py]
 (qtl_analysis conda env)            (qtl_analysis conda env, 동일 로직)
   ▼                                   ▼
[reanalysis/run_all.sh]              [reanalysis/run_all.sh]
 step01_data_qc.R                    step01_data_qc.R
 step02_linkage_map.R                step02_linkage_map.R
 step03_genoprob.R                   step03_genoprob.R
 step04_interval_mapping.R           step04_interval_mapping.R
 step05_permutation.R (n=1000)       step05_permutation.R (n=1000)
 step06_cim_mqm.R                    step06_cim_mqm.R
 step07_qtl_interval.R               step07_qtl_interval.R
 step08_effects.R                    step08_effects.R
 step09_candidates.R                 step09_candidates.R
   ▼                                   ▼
[FINAL] results/09_candidate_genes_HIGH.csv   step10_candidate_annotation.R (BLAST+TAIR10+Araport11)
 (Cr1GR1은 step10 없이 여기서 종료)             ▼
                                     replot_figures.R (보고서용 최종 Figure)
                                       ▼
                                     [FINAL] results/10_candidate_annotation.csv
                                      ← 송부 완료된 최종 후보 유전자 리스트
```

| 스텝 | 소프트웨어(버전) | 핵심 파라미터 | Input | Output |
|---|---|---|---|---|
| 0. Raw data 검증 (비블로킹) | md5sum(coreutils), python3, zcat/wc | `xargs -n 2 -P 50` 병렬 50, `ThreadPoolExecutor(max_workers=10)`, batch_size 50 | 100 샘플 R1/R2 fastq.gz, md5/summary 텍스트 | stdout 리포트만 (파일 산출물 없음) |
| 1. Trimming/QC | fastp **v1.0.1**, multiqc | `--detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --dedup --compression 4 --thread 100`(실제 64 threads로 제한됨) | `rawdata/*_R1/R2_*.fastq.gz` | `clean_data/*_clean_R1/R2.fastq.gz`(200개), `qc_reports/*`, `multiqc_result/` |
| 2. Alignment + 전처리 (Cr1GR1/Cr22.5 병렬) | bwa-mem2 **2.2.1**, samtools **1.21**, GATK **4.6.2.0** | `BWA_THREADS=60(Cr1GR1)/80(Cr22.5) SORT_THREADS=20`; RG 태그; MarkDuplicates `--VALIDATION_STRINGENCY SILENT`; HaplotypeCaller `-ERC GVCF --native-pair-hmm-threads 2` | `clean_data/*_clean_R1/R2.fastq.gz`, 계통별 ragtag reference | `{Cr1GR1,Cr22.5}/bam/*.dedup.bam(.bai)`, `*.dedup_metrics.txt`, `{Cr1GR1,Cr22.5}/gvcf/*.g.vcf.gz(.tbi)` |
| 3. Alignment 통계 (곁가지, 비블로킹) | samtools flagstat (1.21) | `THREADS=8`, `samtools flagstat -@ 2` | `{Cr1GR1,Cr22.5}/bam/*.dedup.bam` | `{Cr1GR1,Cr22.5}/flagstat/*.flagstat`, `*_alignment_stats.tsv` |
| 4. gVCF 병합 + Genotyping | GATK **4.6.2.0** (CombineGVCFs, GenotypeGVCFs) | 기본 옵션 | `gvcf/*.g.vcf.gz` | `merged.gvcf.gz(.tbi)`, `raw_variants.vcf.gz(.tbi)` |
| 5. SNP 추출 + Hard filtering | GATK **4.6.2.0** (SelectVariants, VariantFiltration) | `-select-type SNP`; `--filter-expression "QD<2.0\|\|FS>60.0\|\|MQ<40.0\|\|SOR>3.0\|\|MQRankSum<-12.5\|\|ReadPosRankSum<-8.0" --filter-name Basic_Filter` (Cr1GR1은 Cr22.5 코드를 동일 파라미터로 재사용) | `raw_variants.vcf.gz` | `raw_snps*.vcf.gz(.tbi)` → `filtered_snps*.vcf.gz(.tbi)` |
| 6. VCF 추가 필터링 | vcftools **0.1.16** | `--remove-filtered-all --min-alleles 2 --max-alleles 2 --max-missing 0.7 --maf 0.05 --recode --recode-INFO-all` (Cr1GR1은 Cr22.5 코드 재사용) | `filtered_snps*.vcf.gz` | `clean_snps*.recode.vcf` |
| 7. VCF → R/qtl2 포맷 변환 | python **3.13.9**, `scikit-allel` **1.3.13** (qtl_analysis conda env) | Chi-square 임계값(분리비 왜곡) 20.0, 결측률 임계값 50%, 인코딩 1=A(Alt)/2=H/3=B(Ref) | `clean_snps*.recode.vcf`, phenotype csv(외부 제공 원시 표현형) | `qtl2*_geno/gmap/pheno.csv`, `qtl2*_control.yaml` |
| 8. QTL 재분석 파이프라인 (최종 채택) | R **4.5.2**, R 패키지 `qtl2` **v0.38** (qtl_analysis conda env); Cr22.5 step10 한정 ggplot2 | step01: missing rate 20% 개체 제거, 200kb bin thinning; step02: Kosambi cM; step04: HK 회귀(scan1); step05: permutation **n=1000**; step07: Bayes 95% CI; step08: Beavis 보정 PVE; step09: 피크 200kb 이내=HIGH/500kb 이내=MEDIUM; (Cr22.5 step10: BLAST bitscore 기준 best-hit, `dist_to_peak_kb` tier) | `qtl2*_control.yaml`(geno/gmap/pheno), GFF3 주석(공유), (Cr22.5 step10: 프로젝트 외부 `Functional_annotation/` 리소스) | Cr1GR1: `results/09_candidate_genes(_HIGH).csv`(최종); Cr22.5: `results/09_candidate_genes(_HIGH).csv` → step10 `results/10_candidate_annotation.csv`(최종), `10_top100kb_annotated.csv`, `10_gene_map.pdf` |

## 2. 프로젝트별 특이사항

### Capsella_rubella_rsDNAseq
- (단일 프로젝트 카테고리 — 위 공통 파이프라인이 곧 이 프로젝트의 파이프라인)
- 연구목적/샘플 설계: *Capsella rubella* 1GR1 x 22.5 ecotype F2 100개체를 이용한 flowering time, cauline leaf number 형질의 QTL 연구
- 프로젝트 내부 분기: 파이프라인이 alignment 단계부터 두 참조 계통(`Cr1GR1`, `Cr22.5`)으로 병렬 수행되며, 두 분기는 최종 QTL 재분석 파이프라인(step08) 이후 서로 다른 지점에서 종료된다 — Cr1GR1은 step09(`09_candidate_genes_HIGH.csv`)에서 종료되고, Cr22.5만 step10 기능주석(BLAST+TAIR10+Araport11) 및 `replot_figures.R`까지 이어져 `10_candidate_annotation.csv`가 최종 송부된 후보 유전자 리스트로 확인됨(2026-07-08 프로젝트 수행자 확인).
- Cr22.5 계통 쪽 일부 스크립트(`alignment/bwamem2-GATK.sh`, `variant_call/Cr22.5/vcf/GATK_filtering.sh`)는 실행 이후 내용이 수정되어 현재 파일 내용이 실제 실행 당시 내용과 다름 — 실행 로그 및 산출물 mtime을 근거로 원래 동작을 확정함.
- Cr1GR1 계통은 SNP hard filtering(GATK SelectVariants/VariantFiltration) 및 vcftools 추가 필터링 단계에 별도 `.sh` 파일이 없으며, 이는 Cr22.5의 동일 코드/파라미터를 대화형으로 재사용해 실행한 것으로 확정됨(프로젝트 수행자 확인 원칙).
- Cr22.5에는 초기/구버전 QTL 분석 스크립트군(`qtl2.R`~`qtl2_v4.R`, `qtl2_downstream*.R`, `get_candidate_genes.R` 등)이 폐기 추정 상태로 남아 있으며, 그 산출물(`Candidate_Genes_List.csv` 등)은 `reanalysis/` 파이프라인으로 대체된 중간 산출물임(최종본 아님).
- 남은 미해결 항목(소스 문서 §4 기준, "확인 필요"):
  1. `alignment/bwamem2-GATK.sh`(Cr22.5)의 원본 매핑 코드 내용 — 확인 필요
  2. Cr22.5 bam 파일이 top-level `bam/`에서 `Cr22.5/bam/`으로 재배치된 경로 — 확인 필요
  3. `vcf_to_qtl2.py`(Cr22.5) 실행 시 정확한 `--vcf` CLI 커맨드 전문 — 확인 필요
  4. `parent_geno/`(chr6 synteny 검증) 산출물을 생성한 정확한 명령어 — 확인 필요
  5. multiqc 정확한 버전 — 확인 필요
  6. rawdata 검증 스크립트(`check_reads.py` 등) 실제 실행 여부/시점 — 확인 필요

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Capsella_rubella_rsDNAseq | 단일 프로젝트지만 alignment 단계부터 Cr1GR1/Cr22.5 두 참조 계통으로 완전 병렬 수행됨; 두 계통이 QTL 재분석 파이프라인 step09 이후 서로 다른 지점(Cr1GR1=step09 종료, Cr22.5=step10 기능주석까지)에서 종료 | bwa-mem2 2.2.1, samtools 1.21, GATK 4.6.2.0, vcftools 0.1.16, fastp v1.0.1, python 3.13.9(scikit-allel 1.3.13), R 4.5.2(qtl2 v0.38) — 두 계통 동일 버전 사용 | Cr22.5 step10: BLAST vs Arabidopsis + TAIR10/Araport11 기능주석 병합; chr6 QTL peak 검증용 부모 계통 synteny 분석(`parent_geno/`, 메인 DAG 곁가지) |
