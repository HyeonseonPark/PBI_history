# Resequencing 변이탐지 DNA-seq (Capsella_rubella_rsDNAseq, Actinia_arguat_rsDNAseq)

두 프로젝트 모두 "참조 게놈 정렬 + GATK 기반 변이 콜링(HaplotypeCaller~GenotypeGVCFs)"이라는 동일한 backbone을 공유하지만, 이후 목적(QTL 매핑 vs 분자마커 개발)에 따라 필터링 전략과 다운스트림 분석이 크게 갈라진다. 아래 Step 1~4는 두 프로젝트가 공유하는 상류 파이프라인이며, Step 5~6은 "변이 필터링"이라는 목적은 공유하되 실제 도구/파라미터가 서로 다르다(각 프로젝트 라벨로 구분 표기). Step 0, 7, 8은 Capsella 전용이며, Actinia의 마커설계 다운스트림(SNP/InDel 특이 마커 탐색~rhAmp 서열 추출)은 프로젝트 고유 확장이므로 [2. 프로젝트별 특이사항]에 기술한다.

## 1. 공통 파이프라인

### 1-1. Capsella_rubella_rsDNAseq DAG

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

### 1-2. Actinia_arguat_rsDNAseq DAG (변이콜링 backbone까지)

12품종(다래) resequencing 기반 분자마커 개발 파이프라인의 상류 구간. Capsella와 달리 참조 계통 병렬 분기는 없으며(단일 참조 게놈 `GWHBJWW00000000`), GenotypeGVCFs 이후 SNP/InDel 마커설계 전용 다운스트림(§2 참고)으로 이어진다.

```
[rawdata] Illumina PE fastq (12품종, 11개 2024-11 + AGB 1개 2025-07 추가)
   ▼
[prinseq.sh] prinseq-lite 0.20.4 QC/트리밍 (dedup 포함, -derep 14)
   ▼
[bwa_mem.sh] bwa mem 0.7.19-r1273 → samtools sort
   ▼
[02.picard.sh] gatk AddOrReplaceReadGroups (RG 태깅. ⚠ MarkDuplicates/dedup 단계 없음 — Capsella와의 핵심 차이)
   ├──▶ [flagstat.sh] samtools flagstat (곁가지, 비블로킹, flagstat.log 0바이트 — 확인 필요)
   ▼
[03.GATK.sh] gatk HaplotypeCaller -ERC GVCF (샘플별 12개)
   ▼
[03.GATK.sh] gatk CombineGVCFs → GenotypeGVCFs
   ▼
[Actinidia_arguta.vcf] (12샘플 통합 VCF)
   ▼
(스크립트 미발견, 수동 실행 추정) vcftools --remove-indels 등으로 SNP/INDEL 분리 + 1~2차 필터링
   ▼
Aarguta_SNP.noMiss.GQ95.homo.vcf.gz / Aarguta_INDEL.noMiss.GQ95.homo.vcf.gz
   ▼
   (이후 SNP/InDel 마커설계 다운스트림 — §2. Actinia_arguat_rsDNAseq 참고)
```

### Step 0. Raw data 검증 (Capsella 전용, 비블로킹)
- **소프트웨어(버전):** md5sum(coreutils), python3, zcat/wc
- **핵심 파라미터:** `xargs -n 2 -P 50` 병렬 50, `ThreadPoolExecutor(max_workers=10)`, batch_size 50
- **Input:** 100 샘플 R1/R2 fastq.gz, md5/summary 텍스트
- **Output:** stdout 리포트만 (파일 산출물 없음)
- **실행 커맨드:**
  ```bash
  # verify_md5.sh
  awk '{print $1, $2}' "$md5_file" | xargs -n 2 -P 50 bash -c 'check_md5 "$@"' _
  # check_reads.py (내부: zcat <sample>.fastq.gz | wc -l, ThreadPoolExecutor(max_workers=10), batch_size=50)
  python3 check_reads.py
  # parse_summary.py
  python3 parse_summary.py <summary.tsv>
  ```

### Step 1. Trimming/QC
- **소프트웨어(버전):**
  - Capsella: fastp **v1.0.1**, multiqc
  - Actinia: prinseq-lite **0.20.4**
- **핵심 파라미터:**
  - Capsella: `--detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --dedup --compression 4 --thread 100`(실제 64 threads로 제한됨)
  - Actinia: `-min_len 50 -min_qual_score 5 -min_qual_mean 15 -derep 14 -trim_qual_left 15 -trim_qual_right 15`
- **Input:**
  - Capsella: `rawdata/*_R1/R2_*.fastq.gz`
  - Actinia: `rawdata/*_R1_001.fastq`, `*_R2_001.fastq` (12품종, AGB만 `gzip -d` 선행)
- **Output:**
  - Capsella: `clean_data/*_clean_R1/R2.fastq.gz`(200개), `qc_reports/*`, `multiqc_result/`
  - Actinia: `*_good_1.fastq`, `*_good_2.fastq` (+ `*_bad_*`, `*_singletons`, `*.log`)
- **실행 커맨드:**
  - Capsella (`fastq.sh`, for 루프 대표 1건):
    ```bash
    fastp --in1 <sample>_R1_*.fastq.gz --in2 <sample>_R2_*.fastq.gz --out1 clean_data/<sample>_clean_R1.fastq.gz --out2 clean_data/<sample>_clean_R2.fastq.gz --detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --dedup --compression 4 --thread 100 --html qc_reports/<sample>_fastp.html --json qc_reports/<sample>_fastp.json
    multiqc qc_reports/ -o multiqc_result
    ```
  - Actinia (`prinseq.sh`, AGB 샘플 대표 1건 — 샘플별 백그라운드 `&` 병렬):
    ```bash
    prinseq-lite -fastq AGB_S0_L009_R1_001.fastq -fastq2 AGB_S0_L009_R2_001.fastq -out_format 3 -out_good AGB_S0_L009_good -out_bad AGB_S0_L009_bad -log AGB_S0_L009.log -min_len 50 -min_qual_score 5 -min_qual_mean 15 -derep 14 -trim_qual_left 15 -trim_qual_right 15
    ```

### Step 2. Alignment + 전처리
- **소프트웨어(버전):**
  - Capsella (Cr1GR1/Cr22.5 병렬): bwa-mem2 **2.2.1**, samtools **1.21**, GATK **4.6.2.0**
  - Actinia: bwa **0.7.19-r1273**(bwa-mem2 아님), samtools **1.21**, GATK **4.6.2.0**(AddOrReplaceReadGroups)
- **핵심 파라미터:**
  - Capsella: `BWA_THREADS=60(Cr1GR1)/80(Cr22.5) SORT_THREADS=20`; RG 태그; MarkDuplicates `--VALIDATION_STRINGENCY SILENT`; HaplotypeCaller `-ERC GVCF --native-pair-hmm-threads 2`
  - Actinia: `bwa mem -t 150`, `samtools sort --threads 150`; AddOrReplaceReadGroups `-RGLB Aarguta -RGPL Illumina -RGPU lane1 -SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT -TMP_DIR ./tmp/`; HaplotypeCaller `-ERC GVCF`(native-pair-hmm-threads 옵션 없음)
- **Input:**
  - Capsella: `clean_data/*_clean_R1/R2.fastq.gz`, 계통별 ragtag reference
  - Actinia: `*_good_1/2.fastq`(Step1 output), `GWHBJWW00000000.genome.fasta`
- **Output:**
  - Capsella: `{Cr1GR1,Cr22.5}/bam/*.dedup.bam(.bai)`, `*.dedup_metrics.txt`, `{Cr1GR1,Cr22.5}/gvcf/*.g.vcf.gz(.tbi)`
  - Actinia: `*_good.sorted.bam`(중간 sam/bam 삭제) → `*_RG.bam(.bai)` → `*_RG.gvcf(.idx)`
- **실행 커맨드:**
  - Cr1GR1 (`bwamem2-GATK-Cr1GR1.sh`, for 루프 대표 1건):
    ```bash
    bwa-mem2 mem -t 60 -R "@RG\tID:<sample>\tSM:<sample>\tPL:ILLUMINA\tLB:<sample>" reference/Cr1GR1.ragtag.fasta clean_data/<sample>_clean_R1.fastq.gz clean_data/<sample>_clean_R2.fastq.gz | samtools sort -@ 20 -o Cr1GR1/bam/<sample>.sorted.bam -
    gatk MarkDuplicates -I Cr1GR1/bam/<sample>.sorted.bam -O Cr1GR1/bam/<sample>.dedup.bam -M Cr1GR1/bam/<sample>.dedup_metrics.txt --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT
    gatk HaplotypeCaller -R reference/Cr1GR1.ragtag.fasta -I Cr1GR1/bam/<sample>.dedup.bam -O Cr1GR1/gvcf/<sample>.g.vcf.gz -ERC GVCF --native-pair-hmm-threads 2
    ```
  - Cr22.5 (`bwamem2-GATK.sh` — 매핑 커맨드 라인은 확인 필요, 스크립트 파일 없음: 실행 후 파일이 수정되어 현재 파일에는 MarkDuplicates/HaplotypeCaller만 남아있음; 로그 실측 `Threads: 80 + 20`으로 BWA_THREADS=80 추정. 현재 파일에 남아있는 부분):
    ```bash
    gatk MarkDuplicates -I bam/<sample>.sorted.bam -O bam/<sample>.dedup.bam -M bam/<sample>.dedup_metrics.txt --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT
    gatk HaplotypeCaller -R reference/Cr22.5.ragtag.fasta -I bam/<sample>.dedup.bam -O gvcf/<sample>.g.vcf.gz -ERC GVCF --native-pair-hmm-threads 2
    ```
  - Actinia (`bwa_mem.sh` + `02.picard.sh` + `03.GATK.sh` HaplotypeCaller 블록, for 루프 대표 1건. ⚠ Capsella와 달리 MarkDuplicates(dedup) 단계가 없고 RG 태깅만 수행):
    ```bash
    # bwa_mem.sh
    bwa mem GWHBJWW00000000.genome.fasta -t 150 <sample>_good_1.fastq <sample>_good_2.fastq > <sample>_good.sam
    samtools view --threads 150 -bSh <sample>_good.sam > <sample>_good.bam
    samtools sort --threads 150 -o <sample>_good.sorted.bam <sample>_good.bam
    # 02.picard.sh
    gatk AddOrReplaceReadGroups -I <sample>_good.sorted.bam -O <sample>_RG.bam -SORT_ORDER coordinate -RGID <sample> -RGLB Aarguta -RGPL Illumina -RGPU lane1 -RGSM <sample> -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT -TMP_DIR ./tmp/
    # 03.GATK.sh (HaplotypeCaller 블록)
    gatk HaplotypeCaller -I <sample>_RG.bam -O <sample>_RG.gvcf -R GWHBJWW00000000.genome.fasta -ERC GVCF
    ```

### Step 3. Alignment 통계 (곁가지, 비블로킹)
- **소프트웨어(버전):** samtools flagstat (1.21, 두 프로젝트 동일)
- **핵심 파라미터:**
  - Capsella: `THREADS=8`, `samtools flagstat -@ 2`
  - Actinia: 옵션 없음(기본), 샘플별 백그라운드 `&` 병렬
- **Input:**
  - Capsella: `{Cr1GR1,Cr22.5}/bam/*.dedup.bam`
  - Actinia: `*_RG.bam`
- **Output:**
  - Capsella: `{Cr1GR1,Cr22.5}/flagstat/*.flagstat`, `*_alignment_stats.tsv`
  - Actinia: `*_RG.flagstat` (`flagstat.log`은 0바이트로 실행 로그 미기록 — 확인 필요)
- **실행 커맨드:**
  - Capsella (`get_alignment_stats.sh`, GENOME in Cr22.5,Cr1GR1 루프 대표 1건, THREADS=8로 병렬 throttle):
    ```bash
    samtools flagstat -@ 2 <ALIGN_DIR>/<GENOME>/bam/<sample>.dedup.bam > <ALIGN_DIR>/<GENOME>/flagstat/<sample>.flagstat
    ```
  - Actinia (`flagstat.sh`):
    ```bash
    samtools flagstat <sample>_RG.bam > <sample>_RG.flagstat
    ```

### Step 4. gVCF 병합 + Genotyping
- **소프트웨어(버전):** GATK **4.6.2.0** (CombineGVCFs, GenotypeGVCFs — 두 프로젝트 동일)
- **핵심 파라미터:** 기본 옵션 (두 프로젝트 동일)
- **Input:**
  - Capsella: 계통별 `gvcf/*.g.vcf.gz`
  - Actinia: `*_RG.gvcf` 12개 (`gvcf.list`로 목록화)
- **Output:**
  - Capsella: `merged.gvcf.gz(.tbi)`, `raw_variants.vcf.gz(.tbi)`
  - Actinia: `combine.vcf(.idx)`, `Actinidia_arguta.vcf(.idx/.stat)` (12샘플 통합)
- **실행 커맨드:**
  - Capsella (`Merging_gvcf.sh`, REF만 상이):
    ```bash
    # Cr1GR1
    gatk CombineGVCFs -R /pbi-acc1/hsPark/Capsella_rubella_rsDNAseq/alignment/reference/Cr1GR1.ragtag.fasta --variant gvcf.list -O merged.gvcf.gz
    gatk GenotypeGVCFs -R /pbi-acc1/hsPark/Capsella_rubella_rsDNAseq/alignment/reference/Cr1GR1.ragtag.fasta -V merged.gvcf.gz -O raw_variants.vcf.gz
    # Cr22.5
    gatk CombineGVCFs -R ../reference/Cr22.5.ragtag.fasta --variant gvcf.list -O merged.gvcf.gz
    gatk GenotypeGVCFs -R ../reference/Cr22.5.ragtag.fasta -V merged.gvcf.gz -O raw_variants.vcf.gz
    ```
  - Actinia (`03.GATK.sh`, CombineGVCFs/GenotypeGVCFs 블록. 파일 하나에 3개 블록이 있으나 실제로는 각각 다른 시점에 주석 처리/해제하며 수동 실행한 것으로 추정 — 확인 필요):
    ```bash
    ls *gvcf | sed 's/\t/\n/g' > gvcf.list
    gatk CombineGVCFs -O combine.vcf -R GWHBJWW00000000.genome.fasta -V gvcf.list
    gatk GenotypeGVCFs -O Actinidia_arguta.vcf -R GWHBJWW00000000.genome.fasta -V combine.vcf
    ```

### Step 5. 변이 필터링 (1차)
Capsella는 GATK 기반 hard filtering(SNP만 선별)을, Actinia는 vcftools 기반 SNP/INDEL 분리를 사용한다 — 목적(노이즈 제거)은 같으나 방법론이 다르다.
- **소프트웨어(버전):**
  - Capsella: GATK **4.6.2.0** (SelectVariants, VariantFiltration)
  - Actinia: VCFtools **0.1.16**
- **핵심 파라미터:**
  - Capsella: `-select-type SNP`; `--filter-expression "QD<2.0||FS>60.0||MQ<40.0||SOR>3.0||MQRankSum<-12.5||ReadPosRankSum<-8.0" --filter-name Basic_Filter` (Cr1GR1은 Cr22.5 코드를 동일 파라미터로 재사용)
  - Actinia: `--remove-indels --recode-INFO-all --recode` (SNP측 확인됨); INDEL측 대응 명령(`--keep-only-indels` 등)은 로그 미발견 — 확인 필요
- **Input:**
  - Capsella: `raw_variants.vcf.gz`
  - Actinia: `Actinidia_arguta.vcf`
- **Output:**
  - Capsella: `raw_snps*.vcf.gz(.tbi)` → `filtered_snps*.vcf.gz(.tbi)`
  - Actinia: `Aarguta_SNP.recode.vcf`(SNP, 92,657,724/108,036,783 사이트 유지), `Aarguta_INDEL.vcf.gz`(INDEL)
- **실행 커맨드:**
  - Capsella (`GATK_filtering.sh` 1~2단계, Cr22.5 파일에 주석 처리되어 남아있음 — 실행 당시 커맨드로 확정; Cr1GR1은 Cr22.5 스크립트와 동일 커맨드 재사용, 대화형 실행 — 별도 .sh 없음):
    ```bash
    gatk SelectVariants -R /pbi-acc1/hsPark/Capsella_rubella_rsDNAseq/alignment/reference/Cr22.5.ragtag.fasta -V raw_variants.vcf.gz -select-type SNP -O raw_snps.vcf.gz
    gatk VariantFiltration -R /pbi-acc1/hsPark/Capsella_rubella_rsDNAseq/alignment/reference/Cr22.5.ragtag.fasta -V raw_snps.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Basic_Filter" -O filtered_snps.vcf.gz
    ```
  - Actinia (스크립트 미발견 — `Aarguta_SNP.log` 로그 근거로 SNP측만 확정, 수동 실행 추정):
    ```bash
    vcftools --vcf Actinidia_arguta.vcf --remove-indels --recode-INFO-all --recode --out Aarguta_SNP
    # INDEL측 대응 명령: 확인 필요 — 스크립트 파일 없음
    ```

### Step 6. 변이 필터링 (2차/추가)
- **소프트웨어(버전):**
  - Capsella: vcftools **0.1.16**
  - Actinia: 확인 필요 — 스크립트 파일 없음 (vcftools/bcftools 조합으로 추정)
- **핵심 파라미터:**
  - Capsella: `--remove-filtered-all --min-alleles 2 --max-alleles 2 --max-missing 0.7 --maf 0.05 --recode --recode-INFO-all` (Cr1GR1은 Cr22.5 코드 재사용)
  - Actinia: 파일명(`noMiss.GQ95.homo`)에서 무결측/GQ≥95/동형접합 필터로 추정되나 정확한 커맨드/파라미터는 확인 필요 — 스크립트 파일 없음
- **Input:**
  - Capsella: `filtered_snps*.vcf.gz`
  - Actinia: `Aarguta_SNP.recode.vcf` / `Aarguta_INDEL.vcf.gz`
- **Output:**
  - Capsella: `clean_snps*.recode.vcf`
  - Actinia: `Aarguta_SNP.noMiss.GQ95.homo.vcf.gz(.tbi)`, `Aarguta_INDEL.noMiss.GQ95.homo.vcf.gz`
- **실행 커맨드:**
  - Capsella (`GATK_filtering.sh` 3단계, Cr22.5 현재 실행 활성 코드; Cr1GR1은 동일 커맨드 재사용, 대화형 실행 — 별도 .sh 없음):
    ```bash
    vcftools --gzvcf filtered_snps.vcf.gz --remove-filtered-all --min-alleles 2 --max-alleles 2 --max-missing 0.7 --maf 0.05 --recode --recode-INFO-all --out clean_snps
    ```
  - Actinia: 확인 필요 — 스크립트 파일 없음. `GATK/Homo_filter.py`가 동형접합 판정 로직(탭구분 GT 3~13열, `/` 또는 `|` 양쪽이 같으면 True)을 갖고 있으나 입출력 형식이 최종 `.vcf.gz`와 맞지 않아 실제 사용 여부 확인 불가.

### Step 7. VCF → R/qtl2 포맷 변환 (Capsella 전용)
- **소프트웨어(버전):** python **3.13.9**, `scikit-allel` **1.3.13** (qtl_analysis conda env)
- **핵심 파라미터:** Chi-square 임계값(분리비 왜곡) 20.0, 결측률 임계값 50%, 인코딩 1=A(Alt)/2=H/3=B(Ref)
- **Input:** `clean_snps*.recode.vcf`, phenotype csv(외부 제공 원시 표현형)
- **Output:** `qtl2*_geno/gmap/pheno.csv`, `qtl2*_control.yaml`
- **실행 커맨드:**
  ```bash
  # Cr1GR1 (argparse 정의 그대로; 정확한 커맨드라인 로그는 없으나 출력 파일명 패턴으로 확정)
  python3 vcf_to_qtl2_1GR1.py --vcf <clean_snps_1GR1.recode.vcf> --pheno <pheno.csv> --out_prefix qtl2_1GR1
  # Cr22.5 — 정확한 --vcf 전문은 확인 필요 (실행 로그 없음, 입력파일명은 동일 디렉토리 vcf_to_rqtl.py 로그로 추정)
  python3 vcf_to_qtl2.py --vcf <clean_snps.recode.vcf> --pheno F2_phenotypes.csv --out_prefix qtl2
  ```

### Step 8. QTL 재분석 파이프라인 (Capsella 전용, 최종 채택)
- **소프트웨어(버전):** R **4.5.2**, R 패키지 `qtl2` **v0.38** (qtl_analysis conda env); Cr22.5 step10 한정 ggplot2
- **핵심 파라미터:** step01: missing rate 20% 개체 제거, 200kb bin thinning; step02: Kosambi cM; step04: HK 회귀(scan1); step05: permutation **n=1000**; step07: Bayes 95% CI; step08: Beavis 보정 PVE; step09: 피크 200kb 이내=HIGH/500kb 이내=MEDIUM; (Cr22.5 step10: BLAST bitscore 기준 best-hit, `dist_to_peak_kb` tier)
- **Input:** `qtl2*_control.yaml`(geno/gmap/pheno), GFF3 주석(공유), (Cr22.5 step10: 프로젝트 외부 `Functional_annotation/` 리소스)
- **Output:** Cr1GR1: `results/09_candidate_genes(_HIGH).csv`(최종); Cr22.5: `results/09_candidate_genes(_HIGH).csv` → step10 `results/10_candidate_annotation.csv`(최종), `10_top100kb_annotated.csv`, `10_gene_map.pdf`
- **실행 커맨드:**
  ```bash
  # 공통(reanalysis/run_all.sh, Cr1GR1/Cr22.5 동일)
  bash run_all.sh   # 또는 bash run_all.sh --from=N
  # 내부적으로 Rscript step01_data_qc.R ~ Rscript step09_candidates.R를 순차 실행, 로그는 logs/stepNN.log
  # Cr22.5 전용 추가 실행 (run_all.sh에 미포함, 별도 실행)
  Rscript step10_candidate_annotation.R
  Rscript replot_figures.R
  ```

## 2. 프로젝트별 특이사항

### Capsella_rubella_rsDNAseq
- 연구목적/샘플 설계: *Capsella rubella* 1GR1 x 22.5 ecotype F2 100개체를 이용한 flowering time, cauline leaf number 형질의 QTL 연구
- 공통 파이프라인과의 차이점: alignment 단계부터 두 참조 계통(`Cr1GR1`, `Cr22.5`)으로 병렬 수행되며, 두 분기는 최종 QTL 재분석 파이프라인(step08) 이후 서로 다른 지점에서 종료된다 — Cr1GR1은 step09(`09_candidate_genes_HIGH.csv`)에서 종료되고, Cr22.5만 step10 기능주석(BLAST+TAIR10+Araport11) 및 `replot_figures.R`까지 이어져 `10_candidate_annotation.csv`가 최종 송부된 후보 유전자 리스트로 확인됨(2026-07-08 프로젝트 수행자 확인).
- 해당 프로젝트 고유 확장 분석: Step 7(VCF→R/qtl2 변환) 및 Step 8(QTL 재분석 파이프라인, step01~09/10)은 Capsella 전용이며, Cr22.5 step10 기능주석(BLAST vs Arabidopsis + TAIR10/Araport11) 및 chr6 QTL peak 검증용 부모 계통 synteny 분석(`parent_geno/`, 메인 DAG 곁가지)이 추가로 존재.
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

### Actinia_arguat_rsDNAseq
- 연구목적/샘플 설계: *Actinidia arguta*(다래) 12품종(오텀센스/칠보/청산/청연/대보/대성/그린볼/그린하트/광산/명주15/남양주11/새한) 구분용 SNP·InDel 분자마커 개발 — IDT rhAmp genotyping assay 설계가 최종 목표. 11품종은 2024-11에, AGB(그린볼) 1품종은 2025-07에 별도 배치로 뒤늦게 처리됨(시료 생산 지연).
- 공통 파이프라인과의 차이점: (1) 참조 계통 병렬 분기 없이 단일 참조 게놈(`GWHBJWW00000000`)만 사용; (2) 정렬 툴이 bwa-mem2가 아닌 bwa mem(0.7.19)이며, GATK MarkDuplicates(dedup) 단계 자체가 없고 AddOrReplaceReadGroups로 RG 태깅만 수행; (3) 변이 필터링이 GATK VariantFiltration(QD/FS/MQ/SOR 기반 hard filter) 방식이 아니라 vcftools 기반 SNP/INDEL 분리 + noMiss/GQ95/homo 필터 방식(정확한 2차 필터 커맨드는 스크립트 부재로 확인 필요).
- 해당 프로젝트 고유 확장 분석 — SNP/InDel 마커 설계 (Step 6 산출물 `Aarguta_{SNP,INDEL}.noMiss.GQ95.homo.vcf.gz` 이후, SNP/INDEL 두 계열이 거의 동일한 구조로 병렬 진행됨):
  - **샘플별 특이 마커 탐색** (`GATK/SNP/bcftools_filter.sh`): 12개 샘플 각각에 대해 해당 샘플만 RR(또는 AA)이고 나머지 11개는 모두 반대 유전형인 사이트를 bcftools 표현식으로 탐색, 샘플별 백그라운드 병렬(`&`) 실행.
    ```bash
    bcftools view -i "(GT[0]=='RR' && GT[1]=='AA' && ... ) || (GT[0]=='AA' && GT[1]=='RR' && ...)" Aarguta_SNP.noMiss.GQ95.homo.vcf.gz -o <sample>_specific.vcf
    ```
    (실제로는 12샘플 전체 조건을 대상 샘플 인덱스에 따라 동적 생성 — 위는 대표 형태)
  - **4샘플(ACB/ACS/ADB/AMJ) 서브셋 다형성 마커 집계** (`convert_to_vcf.sh` → `SNP_combination2.sh`/`SNP_combination3.sh`, INDEL측은 동일 로직을 `GATK/INDEL/SNP_combination3.sh`에 재사용 — 파일명이 "SNP_"이나 내용은 INDEL 처리, 명명 불일치):
    ```bash
    bcftools view -s "ACB_S11_L002,ACS_S19_L003,ADB_S20_L003,AMJ_S12_L002" -o four_samples.subset.vcf.gz -Oz Aarguta_SNP.noMiss.GQ95.homo.vcf.gz
    bcftools index four_samples.subset.vcf.gz
    bcftools query -f '%CHROM\t%POS\t%QUAL[\t%GT]\n' four_samples.subset.vcf.gz | awk '...'  # 유전형 패턴별 최고 QUAL 마커 집계
    ```
    `SNP_combination2.sh`와 `3.sh` 중 실제 최종 채택본은 mtime만으로 특정 불가(둘 다 2025-10-29 동일) — 확인 필요.
  - **12샘플 전체 대상 최소 마커셋 탐색** (`All_combination.py`, SNP/INDEL 각 디렉토리에 동일 로직 존재): strict biallelic(0\|0/1\|1만 허용) 필터 후 greedy set-cover 알고리즘으로 12샘플(66쌍) 전부를 구분하는 최소 SNP/InDel 조합 탐색.
    ```bash
    python3 All_combination.py
    # 내부: bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT]\n' Aarguta_SNP.noMiss.GQ95.homo.vcf.gz
    ```
    출력: `final_fingerprint_matrix.txt`, `backup_rhamp_markers_by_qual.txt` (SNP측). INDEL측은 스크립트 하드코딩과 달리 실제 산출물명이 `INDEL__final_fingerprint_matrix.txt`(언더스코어 2개)/`INDEL_backup_rhamp_markers_by_qual.txt`로 불일치 — 확인 필요.
  - **rhAmp 마커 101bp 서열 추출** (`All_combination_fasta_extract.sh` / INDEL측 `INDEL_fasta_extract.sh`): 변이 위치 기준 양쪽 50bp(`FLANK_SIZE=50`)를 bedtools로 추출 후 `[REF/ALT]` 표기 삽입 — 이 fasta가 rhAmp/IDT 마커 설계 최종 결과물.
    ```bash
    awk -v flank=50 'NR>1 && NF>0 {...; print chr, start, end, marker_id}' backup_rhamp_markers_by_qual.txt | \
      bedtools getfasta -fi ../GWHBJWW00000000.genome.fasta -bed stdin -fo /dev/stdout -name -s | \
      awk '...'  # REF/ALT 주석 삽입
    ```
    출력: `rhamp_marker_sequences_101bp_annotated.fasta`(SNP), `INDEL_marker_sequences_101bp_annotated.fasta`(InDel). InDel측은 마지막 3개 후보를 수동 선별(`INDEL_target_sites.txt`, 스크립트 없음)한 뒤 서열을 추출하며, `Shifted_InDel.bed/fasta`(좌표 보정본, 생성 스크립트 불명)가 추가로 존재.
  - **PCA 분석** (plink → `vcf_pca.R`, SNP/INDEL 각각 별도 수행, plink 실행 스크립트는 미발견·로그로만 확인):
    ```bash
    plink --vcf Aarguta_SNP.noMiss.GQ95.homo.vcf.gz --allow-extra-chr --const-fid --pca --out SNP_pca
    Rscript vcf_pca.R   # ggplot2, PCA_plot.png/pdf 생성
    ```
- 남은 미해결 항목(소스 문서 §4/부록A 기준, "확인 필요"):
  1. Step 8/9(SNP·INDEL 1차/2차 필터링) 실행 스크립트 부재 — `Aarguta_SNP.recode.vcf`, `*.noMiss.GQ95.homo.vcf.gz` 등을 생성한 정확한 커맨드는 로그(`Aarguta_SNP.log`, 1차 필터만) 외에는 확인 불가.
  2. `All_combination.py`(INDEL) 산출물 파일명이 스크립트 하드코딩(`final_fingerprint_matrix.txt`)과 실제 디렉토리(`INDEL__final_fingerprint_matrix.txt`)가 불일치 — 수동 rename 여부 확인 필요.
  3. `INDEL_target_sites.txt`(3개 마커) 수동 선별 기준/과정 — 근거 스크립트 없음.
  4. `Shifted_InDel.bed/fasta` 생성 커맨드 — 확인 필요(좌표 보정 수작업으로 추정).
  5. 프로젝트 전체에 conda env/`environment.yml`/컨테이너 정의 파일이 전혀 없어 실행 환경 재현 불가 — 확인 필요.
  6. `bwa-mem/flagstat.log`가 0바이트 — 실행 커맨드/시각 세부 확인 불가.
  7. `SNP_combination2.sh` vs `SNP_combination3.sh` 중 실제 최종 채택본 불명(둘 다 동일 output을 덮어씀).
  8. 구버전 `SNP_combination.py`(2025-01-25) 및 `SNP_combination_QUAL_sorting.py`의 실제 파이프라인 포함 여부 불명(요구 입력파일 부재).
  9. `rhamp_marker_sequences_101bp.fasta`(주석 없는 버전) 생성 스크립트 불명.
  10. `PCA_plot.*` vs `PCA_SNP_plot.*`, `my_pca.*` vs `SNP_pca.*`(INDEL측도 동일) — 각각 두 세트가 공존하며 어느 쪽이 최종 채택본인지 확인 필요.
  11. `GATK/INDEL/SNP_combination3.sh` 파일명이 SNP측과 동일하게 남아있는 명명 혼선(내용은 INDEL 처리) — 정정 이력 확인 불가.
  12. `INDEL_target_sites.txt`(2025-03-19)가 그 입력인 `INDEL_backup_rhamp_markers_by_qual.txt`(2025-10-30)보다 mtime이 앞서는 날짜 역전 — 원인 확인 필요.

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Capsella_rubella_rsDNAseq | 단일 프로젝트지만 alignment 단계부터 Cr1GR1/Cr22.5 두 참조 계통으로 완전 병렬 수행됨; 두 계통이 QTL 재분석 파이프라인 step09 이후 서로 다른 지점(Cr1GR1=step09 종료, Cr22.5=step10 기능주석까지)에서 종료 | bwa-mem2 2.2.1, samtools 1.21, GATK 4.6.2.0, vcftools 0.1.16, fastp v1.0.1, python 3.13.9(scikit-allel 1.3.13), R 4.5.2(qtl2 v0.38) — 두 계통 동일 버전 사용 | Cr22.5 step10: BLAST vs Arabidopsis + TAIR10/Araport11 기능주석 병합; chr6 QTL peak 검증용 부모 계통 synteny 분석(`parent_geno/`, 메인 DAG 곁가지) |
| Actinia_arguat_rsDNAseq | 단일 참조 게놈(GWHBJWW00000000) 대상, 참조 계통 병렬 분기 없음; MarkDuplicates(dedup) 단계 없이 RG 태깅만 수행; GATK hard filter 대신 vcftools 기반 SNP/INDEL 분리+noMiss/GQ95/homo 필터 사용; QTL 분석 대신 SNP/InDel 조합 기반 rhAmp 분자마커 설계로 종결 | bwa 0.7.19-r1273(bwa-mem2 아님), samtools 1.21, GATK 4.6.2.0, vcftools 0.1.16, bcftools 1.22, bedtools 2.31.1, PLINK 1.9.0-b.7.11, prinseq-lite 0.20.4 | 12품종 전체 구분용 최소 SNP/InDel 조합 탐색(greedy set-cover, `All_combination.py`) 및 샘플별 특이마커 탐색(`bcftools_filter.sh`); rhAmp/IDT genotyping용 101bp 마커 서열 추출(bedtools getfasta); plink PCA 품종 구조 확인 |
