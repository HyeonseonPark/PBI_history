# Multi-omics/Holobiome (Panax_ginseng_holobiome)

## 1. 공통 파이프라인

본 프로젝트는 단일 프로젝트 카테고리이며, 두 개의 오믹스 축(16S 근권 메타지놈, mRNA-seq
전사체)이 각각 독립적으로 처리된 후 MOFA+를 통해 통합된다. 소스 문서 기준으로 실제
채택·완료된 경로만 아래에 표기하며, 실패·미채택 스크립트(구버전 경로)는 표 하단에 별도
기술한다.

```
[16S 축]
[업스트림 pb16S nextflow (ginseng_qiime2), results_ginseng_v2]
   ──{feature-table-tax.biom, best_taxonomy.tsv, metadata.tsv}──▶
[01_load_and_filter.R] ──{ps_filtered.rds, asv_prevalence_abundance.tsv}──▶
   ├──▶ [02_clr_and_tier_classification.R]
   │        ──{otu_clr.rds, ancombc2_result.rds, asv_tier_classification.tsv,
   │           ps_selected.rds, otu_clr_selected.rds}──▶
   │        ├──▶ [03_diversity_baseline.R] (alpha/beta diversity, PERMANOVA)
   │        └──▶ [05_baseline/01_ancombc2_viz.R] ──▶
   │             [05_baseline/02_spieceasi_network.R] ──{spieceasi_*.rds}──▶
   │             [05_baseline/03_sparcc_network.R] ──{sparcc_*.rds, network_comparison.tsv}──▶
   └──▶ (otu_clr_selected.rds → 통합 단계로 전달)

[mRNA-seq 축]
[fastp.sh] ──{clean_data/*_clean_R{1,2}.fastq.gz}──▶
[01_build_hisat2_index.sh] ──{index/Pg_hisat2*.ht2}──▶
[02_align_hisat2.sh] ──{outputs/bam/*.sorted.bam}──▶
[03_stringtie.sh] ──{outputs/stringtie/*/transcripts.gtf}──▶
[04_prepDE.sh] ──{outputs/counts/gene_count_matrix.csv}──▶
[05_edger_normalization.R] ──{dge.rds, logcpm.rds, DE_*.tsv}──▶
   ├──▶ [06_wgcna.R] ──{gene_module_assignment.tsv, module_eigengenes.csv}──▶ (통합 단계로 전달)
   └──▶ [09_go_kegg_enrichment.R] ◀── [08_build_annotation.R] ◀── [07_eggnog_annotation.sh]
                                            (DEG셋 + WGCNA 모듈셋 GO/KEGG enrichment)

[통합]
otu_clr_selected.rds(16S) + logcpm.rds(mRNA-seq) + module_eigengenes.csv(mRNA-seq)
   ──▶ [05_baseline/04_mofa.R] ──{modelA_model.hdf5, modelB_model.hdf5,
        factor_scores.pdf, factor_year_correlation.tsv, top_weights_*.tsv}
```

### Step 0. 원시데이터 QC (mRNA-seq)
- **소프트웨어(버전):** fastp(버전 확인 필요), MultiQC
- **핵심 파라미터:** `--cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --dedup --thread 256`
- **Input:** `*_R1_*/*_R2_*.fastq.gz`(Illumina raw, 30샘플)
- **Output:** `clean_data/*_clean_R{1,2}.fastq.gz`, `qc_reports/*_fastp.{html,json}`, `multiqc_result/`
- **실행 커맨드:**
  ```bash
  for fq1 in *_R1_*.fastq.gz; do
    fq2="${fq1/_R1_/_R2_}"
    fastp --in1 "$fq1" --in2 "$fq2" --out1 "$out1" --out2 "$out2" \
      --detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 4 \
      --cut_mean_quality 20 --dedup --compression 4 --thread 256 \
      --html "$report_html" --json "$report_json"
  done
  multiqc qc_reports/ -o multiqc_result
  ```

### Step 업스트림 16S ASV 생성
- **소프트웨어(버전):** QIIME2 2024.10 amplicon (`ginseng_qiime2` env), PacBio full-length pb16S nextflow
- **핵심 파라미터:** 확인 필요(외부 서버 처리, 상세 파라미터 미확인)
- **Input:** PacBio HiFi 16S raw reads
- **Output:** `feature-table-tax.biom`, `best_taxonomy.tsv`, `metadata.tsv` (rawdata 하위 `results_ginseng_v2/`에 저장)
- **실행 커맨드:** 확인 필요 — 외부 처리/스크립트 없음

### Step 1. 16S ASV 필터링
- **소프트웨어(버전):** R(phyloseq, biomformat, dplyr), `ginseng_r` env
- **핵심 파라미터:** Bacteria/Archaea만 유지, Mitochondria/Chloroplast 제거
- **Input:** `feature-table-tax.biom`, `best_taxonomy.tsv`, `metadata.tsv`
- **Output:** `ps_filtered.rds`, `asv_prevalence_abundance.tsv`, `library_sizes.pdf`
- **실행 커맨드:**
  ```r
  biom_obj <- read_biom("feature-table-tax.biom")
  ps_raw <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE), tax_table(tax_mat), sample_data(meta_df))
  ps_filtered <- ps_raw %>% subset_taxa(
    !is.na(Domain) & Domain %in% c("Bacteria", "Archaea") &
    (is.na(Family) | !Family %in% c("Mitochondria", "Chloroplast")) &
    (is.na(Order)  | !Order  %in% c("Chloroplast"))
  )
  ```

### Step 2. CLR/ANCOM-BC2/Tier 분류
- **소프트웨어(버전):** R(compositions CLR, zCompositions CZM, ANCOMBC ancombc2), `ginseng_r`
- **핵심 파라미터:** `fix_formula="year"`, global+trend+pairwise, BH, n_cl=4; Tier2 게이트 prevalence≥50%(3회 재실행 후 최종 확정)
- **Input:** `ps_filtered.rds`
- **Output:** `otu_clr.rds`, `ancombc2_result.rds`, `asv_tier_classification.tsv`, `ps_selected.rds`, `otu_clr_selected.rds`
- **실행 커맨드:**
  ```r
  otu_czm <- cmultRepl(otu_t, label = 0, method = "CZM", output = "p-counts")
  otu_clr <- clr(otu_czm)
  ancombc2(
    data = ps, assay_name = "counts", fix_formula = "year", rand_formula = NULL,
    group = "year", p_adj_method = "BH", alpha = 0.05, n_cl = 4,
    global = TRUE, pairwise = TRUE, dunnet = FALSE, trend = TRUE,
    trend_control = list(contrast = list(matrix(c(1,0,-1,1), nrow=2, byrow=TRUE), matrix(c(-1,0,1,-1), nrow=2, byrow=TRUE)), node = list(2,2), solver = "ECOS", B = 100),
    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)
  )
  ```

### Step 3. 16S diversity baseline
- **소프트웨어(버전):** R(vegan, phyloseq, ggpubr), `ginseng_r`
- **핵심 파라미터:** rarefy depth=32,957, Aitchison(CLR+Euclidean) PRIMARY, PERMANOVA 999 permutations
- **Input:** `ps_filtered.rds`, `ps_selected.rds`, `otu_clr.rds`, `otu_clr_selected.rds`
- **Output:** `alpha_diversity.csv/.pdf`, `pcoa_aitchison*.pdf`, `permanova_year*.txt` 등
- **실행 커맨드:**
  ```r
  ps_rare <- rarefy_even_depth(ps_all, sample.size = rare_depth, replace = FALSE, rngseed = 42)
  estimate_richness(ps_rare, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"))
  ait_dist <- dist(clr_mat, method = "euclidean")
  adonis2(ait_dist ~ year, data = data.frame(year = year_vec), permutations = 999, method = "euclidean")
  ```

### Step 4. HISAT2 인덱스 구축
- **소프트웨어(버전):** HISAT2(hisat2-build, 버전 확인 필요), 시스템 경로(conda env 미명시)
- **핵심 파라미터:** N_THREADS=38
- **Input:** 참조 게놈 fasta/gtf
- **Output:** `index/Pg_hisat2*.ht2`, `index/splice_sites.txt`, `index/exons.txt`
- **실행 커맨드:**
  ```bash
  python3 /usr/bin/hisat2_extract_splice_sites.py "$GTF" > "$INDEX_DIR/splice_sites.txt"
  python3 /usr/bin/hisat2_extract_exons.py "$GTF" > "$INDEX_DIR/exons.txt"
  hisat2-build -p 38 --ss "$INDEX_DIR/splice_sites.txt" --exon "$INDEX_DIR/exons.txt" "$GENOME" "$INDEX"
  ```

### Step 5. HISAT2 정렬
- **소프트웨어(버전):** HISAT2 + samtools(버전 확인 필요)
- **핵심 파라미터:** `--rna-strandness RF`, N_THREADS=38
- **Input:** `clean_data/*_clean_R{1,2}.fastq.gz`, HISAT2 인덱스
- **Output:** `outputs/bam/*.sorted.bam(+.bai,.flagstat.txt)`
- **실행 커맨드:**
  ```bash
  hisat2 -p 38 --rna-strandness RF -x "$INDEX" -1 "$r1" -2 "$r2" \
    --summary-file "$QC_DIR/<sample>.hisat2_summary.txt" \
    2>>"$LOG" | samtools view -bS -o "$unsorted_bam"
  samtools sort -@ 1 -m 16G -o "$sorted" "$unsorted"
  samtools index "$sorted"
  samtools flagstat "$sorted" > "$QC_DIR/<sample>.flagstat.txt"
  ```

### Step 6. StringTie2 어셈블리/정량
- **소프트웨어(버전):** StringTie2(`stringtie`, 버전 확인 필요)
- **핵심 파라미터:** `-eB --rf -G <GFF3>`, N_THREADS=38
- **Input:** `outputs/bam/*.sorted.bam`
- **Output:** `outputs/stringtie/{sample}/{transcripts.gtf, gene_abundances.tsv}`
- **실행 커맨드:**
  ```bash
  stringtie "$bam" -eB --rf -G "$GFF3" \
    -A "$out_dir/gene_abundances.tsv" -o "$out_dir/transcripts.gtf" -p 38
  ```

### Step 7. 카운트 매트릭스 생성
- **소프트웨어(버전):** `prepDE.py3`(외부 경로 `smJeon` 소유, 버전 확인 필요)
- **핵심 파라미터:** `-l 151`(read length)
- **Input:** StringTie2 transcripts.gtf
- **Output:** `gene_count_matrix.csv`, `transcript_count_matrix.csv` (63,713행)
- **실행 커맨드:**
  ```bash
  python3 /pbi-acc1/smJeon/z_Hisat_stringtie_module/module05_prepDE_after_Stringtie.py3 \
    -i "$SAMPLE_LIST" -g "$OUT_DIR/gene_count_matrix.csv" -t "$OUT_DIR/transcript_count_matrix.csv" -l 151
  ```

### Step 8. edgeR 정규화/DE
- **소프트웨어(버전):** edgeR(TMM, QLF), `ginseng_r`
- **핵심 파라미터:** CPM≥1 in ≥3 samples, TMM, `~year` 모델, FDR<0.05 & |logFC|≥1
- **Input:** `gene_count_matrix.csv`
- **Output:** `dge.rds`, `logcpm.rds`, `DE_{3Y_vs_2Y,4Y_vs_2Y,4Y_vs_3Y}.tsv`
- **실행 커맨드:**
  ```r
  dge <- DGEList(counts = counts_raw, group = year)
  keep <- rowSums(cpm(dge) >= 1) >= 3
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  logcpm <- cpm(dge, log = TRUE, prior.count = 2)
  design <- model.matrix(~ year)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design, robust = TRUE)
  qlf <- glmQLFTest(fit, contrast = makeContrasts(year3Y, levels = design))  # 외 4Y_vs_2Y, 4Y_vs_3Y 대비
  ```

### Step 9. WGCNA 모듈 분석
- **소프트웨어(버전):** WGCNA(R, `blockwiseModules`), `ginseng_r`
- **핵심 파라미터:** top 15,000 variable genes, `networkType="signed"`, power 자동선택(R²≥0.85), minModuleSize=30, mergeCutHeight=0.25
- **Input:** `dge.rds`, `logcpm.rds`
- **Output:** `gene_module_assignment.tsv`, `module_eigengenes.csv`, `hub_genes.tsv` 등 (13개 모듈)
- **실행 커맨드:**
  ```r
  sft <- pickSoftThreshold(datExpr, powerVector = 1:20, networkType = "signed", verbose = 2)
  net <- blockwiseModules(
    datExpr, power = power, TOMType = "signed", networkType = "signed",
    minModuleSize = 30, mergeCutHeight = 0.25, maxBlockSize = 15000,
    numericLabels = FALSE, nThreads = 38, saveTOMs = FALSE, verbose = 3
  )
  MEs <- orderMEs(moduleEigengenes(datExpr, module_colors)$eigengenes)
  hubs <- chooseTopHubInEachModule(datExpr, module_colors, power = power)
  ```

### Step 10. eggNOG 기능주석
- **소프트웨어(버전):** eggNOG-mapper(≥2.1.9) + DIAMOND(sensitive, --iterate), `ginseng_eggnog` env
- **핵심 파라미터:** `--cpu 38 --dbmem --go_evidence non-electronic`, eggNOG DB v5.0.2
- **Input:** 외부 경로 단백질서열(.pep, 프로젝트 외부 참조)
- **Output:** `Pg_eggnog.emapper.annotations`
- **실행 커맨드:**
  ```bash
  conda create -y -n ginseng_eggnog -c bioconda -c conda-forge "eggnog-mapper>=2.1.9"
  emapper.py -i "$PEP" --itype proteins -m diamond --cpu 38 --dbmem \
    --go_evidence non-electronic --override --output Pg_eggnog --output_dir "$OUT_DIR"
  ```

### Step 11. 주석 테이블 구축
- **소프트웨어(버전):** R(GO.db, KEGGREST, AnnotationDbi), `ginseng_r`
- **핵심 파라미터:** KEGG 식물 15종 화이트리스트로 map##### 필터링
- **Input:** `Pg_eggnog.emapper.annotations`, `logcpm.rds`
- **Output:** `gene2go.tsv`, `gene2ko.tsv`, `term2gene_*.tsv`, `term2name_*.tsv`
- **실행 커맨드:**
  ```r
  gene2go <- ann %>% dplyr::select(gene_id, GOs) %>% dplyr::filter(!is.na(GOs)) %>% tidyr::separate_rows(GOs, sep = ",")
  go_onto <- AnnotationDbi::select(GO.db, keys = go_ids, columns = c("ONTOLOGY", "TERM"), keytype = "GOID")
  gene2kegg <- ann %>% dplyr::select(gene_id, KEGG_Pathway) %>% tidyr::separate_rows(KEGG_Pathway, sep = ",") %>% dplyr::filter(str_detect(pathway, "^map[0-9]{5}$"))
  pl <- KEGGREST::keggList("pathway", org)  # PLANT_ORGS 15종 화이트리스트 조회
  ```

### Step 12. GO/KEGG enrichment
- **소프트웨어(버전):** clusterProfiler::enricher/compareCluster, enrichplot, `ginseng_enrich` env
- **핵심 파라미터:** DEG: FDR<0.05 & |logFC|≥1, enrichment p<0.05/q<0.2, minGSSize=5/maxGSSize=500; 9 DEG set + 13 WGCNA 모듈
- **Input:** `term2gene_*`/`term2name_*`, `DE_*.tsv`, `gene_module_assignment.tsv`
- **Output:** `enrichment/{deg,wgcna}/*.tsv`, `summary_top_terms.tsv`
- **실행 커맨드:**
  ```r
  clusterProfiler::enricher(
    gene = genes, universe = universe, TERM2GENE = t2g[[cat]], TERM2NAME = t2n[[cat]],
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 500
  )
  clusterProfiler::compareCluster(
    geneClusters = deg_gene_lists, fun = "enricher", universe = universe,
    TERM2GENE = t2g[[cat]], TERM2NAME = t2n[[cat]], pAdjustMethod = "BH",
    pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 500
  )
  ```

### Step 13. ANCOM-BC2 시각화
- **소프트웨어(버전):** R(phyloseq, pheatmap), `ginseng_r`
- **핵심 파라미터:** q<0.05 유의 ASV, top/bottom 40 표시
- **Input:** `ancombc2_result.rds`, `ps_filtered.rds`
- **Output:** `lollipop_*.pdf`, `fc_heatmap.pdf`
- **실행 커맨드:**
  ```r
  sig3 <- res %>% filter(q_year3Y < 0.05) %>% arrange(lfc_year3Y)
  lollipop_plot(sig3, "lfc_year3Y", "se_year3Y", "ANCOM-BC2: 3Y vs 2Y")
  pheatmap(mat, display_numbers = sig_mat, number_color = "black", fontsize_number = 8,
    color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
    breaks = seq(-max(abs(mat), na.rm=TRUE), max(abs(mat), na.rm=TRUE), length.out=101),
    cluster_cols = FALSE, cellwidth = 40)
  ```

### Step 14. SPIEC-EASI 네트워크
- **소프트웨어(버전):** SpiecEasi(GitHub `zdk123/SpiecEasi`, method="mb") + igraph/ggraph, `ginseng_spieceasi`(R 4.5.3)
- **핵심 파라미터:** pooled: Tier1+Tier2, per-year(n=10): Tier1만(탐색적); `pulsar.params(rep.num=20, seed=42)`
- **Input:** `ps_filtered.rds`, `asv_tier_classification.tsv`
- **Output:** `spieceasi_pooled.rds`, `spieceasi_{2Y,3Y,4Y}.rds`, `spieceasi_hubs.tsv`
- **실행 커맨드:**
  ```r
  spiec.easi(
    data = otu_mat, method = "mb",
    lambda.min.ratio = ifelse(exploratory, 5e-2, 1e-2), nlambda = 20,
    pulsar.params = list(rep.num = 20, seed = 42, ncores = 1L)
  )
  ```

### Step 15. SparCC 네트워크
- **소프트웨어(버전):** SparCC(SpiecEasi 패키지 내 함수 추정) + igraph, `ginseng_spieceasi`
- **핵심 파라미터:** bootstrap=100, SPIEC-EASI 대비 Jaccard 비교
- **Input:** Step14와 동일 입력
- **Output:** `sparcc_pooled_cor.rds`, `sparcc_sig_edges.tsv`, `network_comparison.tsv`
- **실행 커맨드:**
  ```r
  sc_est <- sparcc(otu_mat, iter = 20, inner_iter = 10, th = 0.1)
  sc_boot <- sparccboot(otu_mat, sparcc.params = list(iter = 20, inner_iter = 10, th = 0.1), R = 100, ncpus = 1L)
  pval_res <- pval.sparccboot(sc_boot)
  ```

### Step 16. MOFA+ 통합(최종 통합 스텝)
- **소프트웨어(버전):** MOFA2(R), `ginseng_r`
- **핵심 파라미터:** Model A: CLR(888 ASV)+eigengene(14); Model B: CLR(888)+top3000 genes; 샘플ID 매핑 `{연생}y-{반복}`↔`Pg_{연생}Y_{반복}`
- **Input:** `otu_clr_selected.rds`(16S 축) + `logcpm.rds`, `module_eigengenes.csv`(mRNA-seq 축)
- **Output:** `modelA_model.hdf5`, `modelB_model.hdf5`, `*_factor_year_correlation.tsv`, `*_top_weights_*.tsv`
- **실행 커맨드:**
  ```r
  mofa <- create_mofa(list(microbiome = mic_mat, transcriptome = eigen_mat))  # Model A; Model B는 transcriptome = rna_top3k
  model_opts <- get_default_model_options(mofa); model_opts$num_factors <- 10
  train_opts <- get_default_training_options(mofa)
  train_opts$seed <- 42L; train_opts$maxiter <- 2000L; train_opts$convergence_mode <- "slow"
  mofa <- prepare_mofa(mofa, model_options = model_opts, training_options = train_opts)
  mofa <- run_mofa(mofa, outfile = hdf5_path, use_basilisk = TRUE)
  ```

**미채택/실패 경로(참고용, 소스 문서 확정)**: mRNA-seq 정량 단계에서 `03_featurecounts.R`
(Rsubread)이 2회 실패(zlib.h 컴파일 오류, GTF gene_id 파싱 오류)하여 `04_deseq2_normalization.R`,
구버전 `05_wgcna.R`은 실행되지 않았고 완전 폐기됨 — 실제 채택 경로는 StringTie2/prepDE/edgeR/
WGCNA(06_wgcna.R) 경로임.

## 2. 프로젝트별 특이사항

### Panax_ginseng_holobiome
- (단일 프로젝트 카테고리 — 위 공통 파이프라인이 곧 이 프로젝트의 파이프라인)
- 연구목적/샘플 설계: 인삼특작부 2·3·4년생 인삼의 뿌리 전사체와 근권부 미생물(16S 메타지놈)
  간 상호관계 규명을 통한 연작장해 및 진세노사이드 합성 조절 메커니즘 규명(연차별 10반복,
  총 30샘플).
- 두 오믹스 축의 통합은 MOFA+(Step16) 한 지점에서만 이루어지며, 통합 방식은 factor
  분석(잠재 인자 추출)이지 별도의 통합 매트릭스나 크로스레이어 네트워크 구축이 아님.
  `analysis/03_soil_preprocessing/`, `04_integration/`, `06_quantum_kernel/`,
  `07_rare_imputation/`, `08_crosslayer/`, `09_consortium_opt/`는 모두 빈 디렉토리로
  미착수 상태이며, 소스 문서가 "양자 영감(quantum-inspired) 머신러닝 프레임 적용 시도"라고
  기술한 부분은 계획 단계(RESEARCH_STRUCTURE.md 로드맵)에 머물러 있고, 실제 실행된 것은
  고전적 베이스라인(ANCOM-BC2/SPIEC-EASI/SparCC/MOFA+, Phase1~3)까지임. `envs/ginseng_qml.yml`
  (PennyLane-GPU, Qiskit-Aer-GPU) 및 `scripts/test_gpu_backends.py`가 존재하나 현재 채택된
  DAG 어디에서도 호출되지 않음.
- PromethiOn WGS(영남대 공동연구), 토양 이화학분석, 인삼 뿌리 진세노사이드 HPLC는 모두 예정
  단계이며 현재 데이터 없음.
- 소스 문서에서 "확인 필요"로 표시된 항목(예: fastp 정확 버전 및 실행 conda env, HISAT2/
  samtools/StringTie2 정확 버전, `04_prepDE.sh`가 사용하는 외부 스크립트의 버전/재현성,
  Step1의 GTF와 Step3의 GFF3가 완전히 동일한 좌표계인지 여부, `ginseng_enrich` env의 최초
  생성 커맨드)는 확인 필요 상태로 유지한다.

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Panax_ginseng_holobiome | 16S(근권 메타지놈)와 mRNA-seq(뿌리 전사체) 두 축을 각각 독립 처리 후 MOFA+로 통합하는 유일한 멀티오믹스/홀로바이옴 프로젝트. mRNA-seq 정량 경로에 featureCounts/DESeq2/구버전WGCNA 실패 이력이 있으며 StringTie2/prepDE/edgeR/WGCNA로 전환되어 완료됨 | HISAT2/samtools/StringTie2 정확 버전 확인 필요; R 4.3.3(`ginseng_r`) vs R 4.5.3(`ginseng_spieceasi`) 혼용; eggNOG-mapper ≥2.1.9, eggNOG DB v5.0.2 | ANCOM-BC2/SPIEC-EASI/SparCC 미생물 네트워크 분석, WGCNA 유전자 모듈-형질 상관, MOFA+를 통한 16S-mRNA-seq 잠재 인자 통합; 양자 영감 ML 프레임(Phase6~9)은 계획 단계로 미착수 |
