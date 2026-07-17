# Single-cell RNA-seq (Homo_sapiens_BEST1_scRNAseq, Homo_sapiens_Tcell_scRNAseq)

## 1. 공통 파이프라인

두 프로젝트 모두 "10X Genomics 3' 라이브러리 + Cell Ranger 정렬/카운팅"이라는 상류 단계는 공유하지만, downstream 분석은 서로 다른 언어/툴 생태계(BEST1=Seurat/R, Tcell=Scanpy/Python)로 완전히 분리되어 있고 스텝 경계·개수도 다르다. 이에 따라 Step 1(Cell Ranger)만 두 프로젝트 공통 스텝으로 다루고, 이후 downstream은 프로젝트별로 완전히 별도의 Step 번호 계열(BEST1=Step 2~10, Tcell=Step 11~17)로 기술한다. 두 프로젝트의 DAG도 구조가 크게 달라 아래 1-1/1-2로 분리해 표시한다.

### 1-1. Homo_sapiens_BEST1_scRNAseq DAG

```
[rawdata: Best1_S6_L004_{R1,R2,I1,I2}_001.fastq.gz]
   ──▶ [cellranger.sh: cellranger multi (OCM 4-plex demux, GRCh38-2024-A)]
        ──{per_sample_outs/{2_SC,2_Best1,7_SC,7_Best1}/outs/sample_filtered_feature_bc_matrix.h5}──▶
[01_clustering_annotation.R]
   STEP1 Read10X_h5 + CreateSeuratObject ──▶ rds/00_raw_seu_list.rds
   STEP2 샘플별 QC 필터 ──▶ rds/01_qc_seu_list.rds
   STEP3 DoubletFinder(샘플별) ──▶ rds/02_singlet_seu_list.rds
   STEP4 SCTransform(샘플별, glmGamPoi) ──▶ rds/03_sct_seu_list.rds
        │
        ├──{rds/03_sct_seu_list.rds}──▶ [03_donor_separated_clustering.R] (아래 별도 분기)
        │
   STEP5 Merge + Harmony(group.by.vars="donor_id") ──▶ plots/03_clustering/UMAP_before_harmony.png
   STEP6 UMAP + Clustering(res 0.3/0.5/0.8, default 0.5) ──▶ rds/04_integrated_seu.rds
   STEP7 Module Score 기반 세포유형 Annotation(8 celltype) ──▶ rds/05_annotated_seu.rds
   STEP8-A GBM marker 시각화 ──▶ plots/05_markers/*
        │
        ▼ {rds/05_annotated_seu.rds}
   [02_DEG_within_donor.R] STEP A-F: donor 내 KD vs Ctrl DEG(Wilcoxon)
        + fgsea(Hallmark/KEGG) + cross-donor 일치 DEG
        ──▶ plots/06_DEG/{donor2,donor7,combined}/*

[03_donor_separated_clustering.R]  (input: rds/03_sct_seu_list.rds, Harmony 미사용)
   donor2: merge(2_SC,2_Best1)→PCA→UMAP→Cluster(res 0.3)→annotate ──▶ rds/donor2_annotated.rds
   donor7: merge(7_SC,7_Best1)→PCA→UMAP→Cluster(res 0.5)→annotate ──▶ rds/donor7_annotated.rds
```

### 1-2. Homo_sapiens_Tcell_scRNAseq DAG (2nd_generation)

```
[Batch1 rawdata/ (L009, 2025-12-26)] + [Batch2 rawdata/Additional/ (L009, symlink→L010)]
+ [Batch3 rawdata/Additional_v2/ (L009, symlink→L011)]
   │ (심볼릭 링크로 lane 번호만 다르게 위장 → Cell Ranger가 3개 배치를 동일 샘플의 서로 다른
   │  lane으로 인식하고 병합; run_cellranger_merged_v2.sh)
   ▼
[run_cellranger_merged_v2.sh] cellranger count --fastqs=rawdata,rawdata_Additional_L010,rawdata_Additional_v2_L011
   ──▶ cellranger_output/{C,T}_count_merged_v2/outs/filtered_feature_bc_matrix (★최종, Batch1+2+3★)
        │
        ▼
[01_preprocess_cluster.py] sc.read_10x_mtx(C,T) 병합 → QC(min_genes=400,max_genes=4500,max_mt=15%)
   + Scrublet 이중세포 검출(샘플별, expected_doublet_rate=0.05) → 필터 → normalize_total(1e4)+log1p
   → HVG(2000) → scale → PCA(30) → neighbors(30) → UMAP → Leiden(res 0.3/0.5/0.8)
   ──▶ adata/adata_preprocessed.h5ad
        │
        ├──{adata_preprocessed.h5ad}──▶ [02_Tcell_NgR1_NogoA.py]
        │     CD3D/E/G score(leiden_0.5 클러스터 평균≥0.3)로 T세포 클러스터 판정 → subset
        │     → Hb score(HBB 등)≥0.5 적혈구 오염 제거 → 재HVG/PCA/UMAP/Leiden
        │     → RTN4R(NgR1)/RTN4(NogoA) 발현>0 기준 4-group 분류(RTN4R+RTN4+/+-/-+/--)
        │     → Wilcoxon DEG(Positive vs Double-negative; Treatment vs Control in Positive)
        │     ──▶ adata/adata_Tcell.h5ad, tables/DEG_*.csv, composition_NgR1_NogoA_*.csv,
        │         figures/03_Tcell/*, figures/04_NgR1_NogoA/*
        │
        ├──{adata_preprocessed.h5ad}──▶ [03_add_celltype_annotation.py]
        │     (Step 5와 동일 T세포 재추출 로직 재현, random_state 동일 → 동일 클러스터)
        │     → CELL_TYPE_MAP(leiden_0.5 클러스터 0~9 → Naive T/Treg-like/γδ T/Activated T/Treg 등
        │       10종 수동 마커기반 라벨) 부여
        │     ──▶ adata/adata_Tcell.h5ad(cell_type 컬럼 추가, 덮어씀), tables/composition_celltype_CvsT.csv,
        │         figures/03_Tcell/barplot_celltype_composition_CvsT.png
        │
        ├──{C,T_count_merged_v2 raw filtered_feature_bc_matrix + adata_preprocessed.h5ad QC 바코드}
        │     ──▶ [05_make_cpm_matrix.py] raw count 재로딩→QC-passed 바코드만 유지→CPM(1e6)+log1p
        │     ──▶ cpm/adata_cpm.h5ad, cpm/mtx/{cpm,log1p_cpm}/*
        │
        ▼ (Step 5/6/8 산출물 수치를 하드코딩 리터럴로 집계)
   [04_generate_report.py] ──▶ downstream_analysis/analysis_report.html
        ▼
   [ppt/generate_slides.py] ──▶ CNU_Tcell_NgR1_NogoA_Slides.pptx (그림은 플레이스홀더로만 참조)

[QC 교차검증 사이드 브랜치 — 비블로킹, 메인 파이프라인에 피드백 없음]
   {C,T_count_merged_v2 raw filtered_feature_bc_matrix}
        ├──▶ check_doublet_histogram.py → figures/01_QC/scrublet_histogram_{C,T}.png (Scrublet 결과 재확인용)
        └──▶ debug_df.R → debug_df2.R → debug_df3.R → run_doubletfinder.R (seurat_env conda,
              Seurat+DoubletFinder 기반 교차검증) ──▶ tables/doubletfinder_{C,T}_results.csv,
              doubletfinder_summary.csv, figures/01_QC/doubletfinder_histogram_{C,T}.png
              (Scrublet 대비 15~18% 수준 doublet률 확인되었으나 adata_preprocessed.h5ad에는 미반영 —
               최종 파이프라인은 Scrublet 기반 판정을 그대로 사용)
```

### Step 1. Cell Ranger 정렬/카운팅 (공통 상류, 두 프로젝트 별도 실행)
- **소프트웨어(버전):**
  - BEST1: Cell Ranger 10.0.0 (시스템 PATH 직접 설치, conda/컨테이너 미사용; OCM 기능은 cellranger v9.0+ 필요, 본 프로젝트는 v10.0 사용)
  - Tcell: Cell Ranger 10.0.0 (시스템 바이너리 `/usr/local/bin/cellranger-10.0.0/cellranger`, conda 환경 아님)
- **핵심 파라미터:**
  - BEST1: `--localcores=150 --localmem=500`(GB); config CSV 내 `create-bam=true`, `include-introns=true`, reference=GRCh38-2024-A(10x 프리빌트); OCM 바코드 매핑 2_SC=OB1, 2_Best1=OB2, 7_SC=OB3, 7_Best1=OB4 (`cellranger multi` 1회 실행으로 4-plex 동시 디멀티플렉싱)
  - Tcell: `--localcores=150 --localmem=500 --create-bam=true --include-introns=true`, transcriptome=GRCh38.115(Ensembl); `cellranger count`를 C/T 각각 별도 실행하며, `--fastqs=rawdata,rawdata_Additional_L010,rawdata_Additional_v2_L011`로 3개 시퀀싱 배치(원래 모두 L009로 명명됨)를 심볼릭 링크로 L010/L011 lane인 것처럼 위장해 병합(`run_cellranger_merged_v2.sh`) — Cell Ranger 자체의 멀티플렉싱 기능이 아니라 lane 위장 트릭
- **Input:**
  - BEST1: `rawdata/HN00277842_10X_RawData_Outs/Best1/23KWMHLT4/Best1_S6_L004_{R1,R2,I1,I2}_001.fastq.gz`(4개 샘플이 칩 로딩 단계에서 이미 풀링된 단일 GEX 라이브러리), `reference/refdata-gex-GRCh38-2024-A`
  - Tcell: `rawdata/{C,T}_S0_L009_R{1,2}_001.fastq.gz`(Batch1) + `rawdata_Additional_L010/{C,T}_S0_L010_R{1,2}_001.fastq.gz`(Batch2 symlink) + `rawdata_Additional_v2_L011/{C,T}_S0_L011_R{1,2}_001.fastq.gz`(Batch3 symlink), `reference/GRCh38_cellranger_index`(1차 프로젝트 자산에 대한 심볼릭 링크, 현재 삭제됨 — 확인 필요: 재현 시 참조 인덱스 재확보 필요)
- **Output:**
  - BEST1: `metadata/cellranger_output/BrainOrganoid_OCM_4plex/outs/{raw_feature_bc_matrix*, filtered_feature_bc_matrix*, qc_report.html, per_sample_outs/{2_SC,2_Best1,7_SC,7_Best1}/outs/sample_filtered_feature_bc_matrix.h5}`
  - Tcell: `cellranger_output/{C,T}_count_merged_v2/outs/{filtered_feature_bc_matrix, possorted_genome_bam.bam, web_summary.html, metrics_summary.csv}` (Batch1+2+3 병합, 최종/유일하게 downstream 입력으로 사용됨. 구버전 `{C,T}_count`, `{C,T}_count_merged`는 2026-07-08 인수인계 정리 시 삭제됨)
- **실행 커맨드:**
  ```bash
  # BEST1 (cellranger multi, 1회 실행으로 4샘플 디멀티플렉싱)
  cellranger multi --id="BrainOrganoid_OCM_4plex" \
                    --csv="${CONFIG_CSV}" \
                    --localcores=150 \
                    --localmem=500
  ```
  생성되는 `multi_config.csv` (실제 내용):
  ```csv
  [gene-expression]
  reference,/pbi-acc2/hsPark/Homo_sapiens_BEST1_scRNAseq/reference/refdata-gex-GRCh38-2024-A
  create-bam,true
  include-introns,true

  [libraries]
  fastq_id,fastqs,lanes,feature_types
  Best1,/pbi-acc2/hsPark/Homo_sapiens_BEST1_scRNAseq/rawdata/HN00277842_10X_RawData_Outs/Best1/23KWMHLT4,any,Gene Expression

  [samples]
  sample_id,ocm_barcode_ids
  2_SC,OB1
  2_Best1,OB2
  7_SC,OB3
  7_Best1,OB4
  ```
  ```bash
  # Tcell (run_cellranger_merged_v2.sh — 실제 스크립트에서 발췌, C/T 각각 실행)
  # Step 1: Batch2(L009) → L010 심볼릭 링크
  ln -s "${ADDITIONAL_DIR}/${sample_id}_S0_L009_${read}_001.fastq.gz" \
        "${SYMLINK_L010_DIR}/${sample_id}_S0_L010_${read}_001.fastq.gz"
  # Step 2: Batch3(L009) → L011 심볼릭 링크
  ln -s "${ADDITIONAL_V2_DIR}/${sample_id}_S0_L009_${read}_001.fastq.gz" \
        "${SYMLINK_L011_DIR}/${sample_id}_S0_L011_${read}_001.fastq.gz"
  # Step 3: 병합 count (샘플별)
  cellranger count --id="${sample_id}_count_merged_v2" \
                   --transcriptome="${REF_GENOME_DIR}" \
                   --fastqs="${RAW_DATA_DIR},${SYMLINK_L010_DIR},${SYMLINK_L011_DIR}" \
                   --sample="${sample_id}" \
                   --localcores=150 \
                   --localmem=500 \
                   --create-bam=true \
                   --include-introns=true
  ```

### Step 2. Seurat 객체 생성 및 샘플별 QC 필터링 (BEST1 전용)
- **소프트웨어(버전):** R 4.5.3, Seurat 5.4.0, SeuratObject 5.1.0
- **핵심 파라미터:** `CreateSeuratObject(min.cells=3, min.features=100)`; 샘플별 상이한
  QC 임계값 — 예: `2_SC`(min_feat=150, max_feat=7000, min_umi=300, max_umi=50000, max_mt<25),
  `2_Best1`/`7_SC`(min_feat=200, min_umi=500, max_mt<20), `7_Best1`(min_feat=150, min_umi=300,
  max_mt<25); `percent.mt`는 `^MT-`, `percent.rb`는 `^RP[SL]` 패턴
- **Input:** `metadata/cellranger_output/BrainOrganoid_OCM_4plex/outs/per_sample_outs/{sample}/
  sample_filtered_feature_bc_matrix.h5` (Step 1 산출물)
- **Output:** `rds/00_raw_seu_list.rds`, `rds/01_qc_seu_list.rds`, `plots/01_qc/QC_before_filter.pdf`,
  `plots/01_qc/cell_counts_qc.csv`
- **실행 커맨드:**
  ```r
  mat <- Read10X_h5(file.path(BASE_DATA, sample_id, "sample_filtered_feature_bc_matrix.h5"))
  seu <- CreateSeuratObject(counts = mat, project = sample_id, min.cells = 3, min.features = 100)
  seu$donor_id  <- ifelse(grepl("^2_", sample_id), "donor2", "donor7")
  seu$condition <- ifelse(grepl("_SC$", sample_id), "Control", "Best1_KD")
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

  filtered <- subset(seu,
    subset = nFeature_RNA >= params$min_feat & nFeature_RNA <= params$max_feat &
             nCount_RNA   >= params$min_umi  & nCount_RNA   <= params$max_umi  &
             percent.mt   <  params$max_mt)
  ```

### Step 3. DoubletFinder 이중세포 제거 (샘플별 독립 실행, BEST1 전용)
- **소프트웨어(버전):** R 4.5.3, DoubletFinder 2.0.6, Seurat 5.4.0
- **핵심 파라미터:** 전처리 `NormalizeData → FindVariableFeatures(nfeatures=2000) →
  ScaleData → RunPCA(npcs=20) → FindNeighbors(dims=1:15) → FindClusters(res=0.5)`;
  `paramSweep(PCs=1:15, sct=FALSE)` + `find.pK`로 최적 pK 추정; `doubletFinder(pN=0.25,
  reuse.pANN=FALSE, sct=FALSE)`; 샘플별 예상 doublet rate 하드코딩(`2_SC`=0.008,
  `2_Best1`=0.025, `7_SC`=0.037, `7_Best1`=0.005); `nExp`는 `modelHomotypic`으로 보정
- **Input:** `rds/01_qc_seu_list.rds`
- **Output:** `rds/02_singlet_seu_list.rds`, `plots/02_doublet/{sample}_pK_sweep.png`
- **실행 커맨드:**
  ```r
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>% RunPCA(npcs = 20) %>%
    FindNeighbors(dims = 1:15) %>% FindClusters(resolution = 0.5)

  sweep_res  <- paramSweep(seu, PCs = 1:15, sct = FALSE)
  sweep_stat <- summarizeSweep(sweep_res, GT = FALSE)
  pk_val     <- find.pK(sweep_stat)
  best_pk    <- as.numeric(as.character(pk_val$pK[which.max(pk_val$BCmetric)]))

  n_exp_doublet     <- round(doublet_rate * ncol(seu))
  homotypic_prop    <- modelHomotypic(seu$seurat_clusters)
  n_exp_doublet_adj <- round(n_exp_doublet * (1 - homotypic_prop))

  seu_df <- doubletFinder(seu, PCs = 1:15, pN = 0.25, pK = best_pk,
                          nExp = n_exp_doublet_adj, reuse.pANN = FALSE, sct = FALSE)
  seu <- subset(seu_df, subset = doublet_class == "Singlet")
  ```
  (Seurat v5 비호환 등으로 `doubletFinder()`가 실패하면 `nCount_RNA > mean + 3*SD`
  아웃라이어 기준으로 fallback하도록 `tryCatch` 처리되어 있음 — 실제 실행에서 이 fallback이
  발동했는지는 로그 기준 확인 필요.)

### Step 4. SCTransform 정규화 (샘플별, BEST1 전용)
- **소프트웨어(버전):** R 4.5.3, Seurat 5.4.0(glmGamPoi 백엔드)
- **핵심 파라미터:** `vars.to.regress="percent.mt"`, `method="glmGamPoi"`
- **Input:** `rds/02_singlet_seu_list.rds`
- **Output:** `rds/03_sct_seu_list.rds`
- **실행 커맨드:**
  ```r
  seu_list_sct <- lapply(SAMPLES, function(s) {
    SCTransform(seu_list_singlet[[s]], vars.to.regress = "percent.mt",
                method = "glmGamPoi", verbose = FALSE)
  })
  ```

### Step 5. Merge + Harmony 통합 (donor 통합 경로, BEST1 전용)
- **소프트웨어(버전):** R 4.5.3, Seurat 5.4.0, harmony 1.2.3
- **핵심 파라미터:** `SelectIntegrationFeatures(nfeatures=3000)`, `PrepSCTFindMarkers`,
  `RunPCA(npcs=50)`, Harmony 전 UMAP(`reduction="pca"`)으로 배치효과 확인, Harmony
  `group.by.vars="donor_id", assay.use="SCT", reduction="pca", reduction.save="harmony"`
- **Input:** `rds/03_sct_seu_list.rds` (4개 샘플 전체 병합)
- **Output:** (중간 산출물, RDS 저장 없음) `plots/03_clustering/UMAP_before_harmony.png`
- **실행 커맨드:**
  ```r
  features <- SelectIntegrationFeatures(object.list = seu_list_sct, nfeatures = 3000)
  seu_all <- merge(seu_list_sct[["2_SC"]],
                   y = list(seu_list_sct[["2_Best1"]], seu_list_sct[["7_SC"]], seu_list_sct[["7_Best1"]]),
                   add.cell.ids = SAMPLES)
  VariableFeatures(seu_all) <- features
  seu_all <- PrepSCTFindMarkers(seu_all)
  seu_all <- RunPCA(seu_all, npcs = 50, verbose = FALSE)
  seu_all <- RunUMAP(seu_all, dims = 1:30, reduction = "pca",
                     reduction.name = "umap.before.harmony", verbose = FALSE)
  seu_all <- RunHarmony(seu_all, group.by.vars = "donor_id", assay.use = "SCT",
                        reduction = "pca", reduction.save = "harmony", verbose = FALSE)
  ```

### Step 6. UMAP + 클러스터링 (donor 통합 경로, BEST1 전용)
- **소프트웨어(버전):** R 4.5.3, Seurat 5.4.0, clustree 0.5.1
- **핵심 파라미터:** `RunUMAP(reduction="harmony", dims=1:30)`,
  `FindNeighbors(reduction="harmony", dims=1:30)`, `FindClusters(resolution=c(0.3,0.5,0.8))`,
  기본 해상도는 `SCT_snn_res.0.5` 채택
- **Input:** Step 5 산출 `seu_all`(harmony reduction 포함)
- **Output:** `rds/04_integrated_seu.rds`, `plots/03_clustering/clustree_resolution.pdf`,
  `plots/03_clustering/UMAP_overview.{pdf,png}`
- **실행 커맨드:**
  ```r
  seu_all <- RunUMAP(seu_all, reduction = "harmony", dims = 1:30, verbose = FALSE)
  seu_all <- FindNeighbors(seu_all, reduction = "harmony", dims = 1:30, verbose = FALSE)
  seu_all <- FindClusters(seu_all, resolution = c(0.3, 0.5, 0.8), verbose = FALSE)
  print(clustree(seu_all, prefix = "SCT_snn_res."))

  Idents(seu_all) <- "SCT_snn_res.0.5"
  seu_all$seurat_clusters <- seu_all$SCT_snn_res.0.5
  saveRDS(seu_all, file.path(RDS_DIR, "04_integrated_seu.rds"))
  ```

### Step 7. Module Score 기반 세포유형 Annotation (donor 통합 경로, BEST1 전용)
- **소프트웨어(버전):** R 4.5.3, Seurat 5.4.0(`AddModuleScore`), dplyr 1.2.0
- **핵심 파라미터:** 8개 세포유형 gene signature(GBM_tumor, Astrocyte, Neuron, OPC,
  Oligodendrocyte, Microglia, Endothelial, Proliferating; 유전자 ≥3개 매칭 시에만 스코어링),
  세포별 dominant cell type = 8개 module score 중 최댓값 유형
- **Input:** `rds/04_integrated_seu.rds`
- **Output:** `rds/05_annotated_seu.rds`, `plots/04_annotation/UMAP_celltype_label.{pdf,png}`,
  `plots/04_annotation/cluster_module_scores.csv`,
  `plots/04_annotation/cluster_dominant_celltype.csv`,
  `plots/04_annotation/celltype_composition_by_sample.{pdf,png}`
- **실행 커맨드:**
  ```r
  cell_signatures <- list(
    GBM_tumor = c("EGFR","SOX2","NES","CD44","CHI3L1","OLIG2","GFAP","PTEN","TP53","PDGFRA","MET","FGFR1"),
    Astrocyte = c("GFAP","S100B","ALDH1L1","AQP4","SLC1A2","GJA1","GLUL","ID3","CLU"),
    Neuron    = c("MAP2","RBFOX3","SYP","TUBB3","DCX","NEFM","NEFL","GAD1","GAD2","SLC17A7"),
    OPC       = c("OLIG1","OLIG2","SOX10","PDGFRA","CSPG4","NKX2-2","MYT1","ASCL1"),
    Oligodendrocyte = c("MBP","MOG","MAG","PLP1","CNP","MOBP","CLDN11"),
    Microglia = c("CX3CR1","P2RY12","TMEM119","IBA1","AIF1","CSF1R","FCRLS","HEXB","SALL1"),
    Endothelial = c("PECAM1","VWF","CDH5","ENG","ESAM","KDR","CLDN5"),
    Proliferating = c("MKI67","TOP2A","CDK1","PCNA","CCNB1","CCNA2","UBE2C","CENPF")
  )
  for (ct in names(cell_signatures)) {
    genes_present <- intersect(cell_signatures[[ct]], rownames(seu_all))
    if (length(genes_present) >= 3) {
      seu_all <- AddModuleScore(seu_all, features = list(genes_present),
                                name = paste0("score_", ct), assay = "SCT")
    }
  }
  seu_all$cell_type <- apply(seu_all@meta.data[, score_cols, drop = FALSE], 1,
    function(x) sub("^score_", "", names(which.max(x))))
  saveRDS(seu_all, file.path(RDS_DIR, "05_annotated_seu.rds"))
  ```

### Step 8. GBM 마커 유전자 시각화 (donor 통합 경로, BEST1 전용)
- **소프트웨어(버전):** R 4.5.3, Seurat 5.4.0, ggplot2 4.0.3, patchwork 1.3.2
- **핵심 파라미터:** `GBM_MARKERS = c("EGFR","SOX2","OLIG2","NES","CD44","GFAP")`;
  FeaturePlot(`min.cutoff="q5", max.cutoff="q95", order=TRUE`), 샘플별 split.by 플롯,
  VlnPlot/DotPlot/RidgePlot(assay="SCT")
- **Input:** `rds/05_annotated_seu.rds`(스크립트 내 동일 세션 변수 재사용)
- **Output:** `plots/05_markers/{FeaturePlot_GBM_markers.pdf, FeaturePlot_split_by_sample.pdf,
  VlnPlot_GBM_markers_by_sample.{pdf,png}, DotPlot_GBM_markers_by_cluster.{pdf,png},
  RidgePlot_GBM_markers_by_sample.pdf, GBM_marker_expression_summary.csv}`
- **실행 커맨드:**
  ```r
  GBM_MARKERS <- c("EGFR", "SOX2", "OLIG2", "NES", "CD44", "GFAP")
  for (gene in GBM_MARKERS) {
    FeaturePlot(seu_all, features = gene, reduction = "umap",
                min.cutoff = "q5", max.cutoff = "q95", pt.size = 0.3, order = TRUE)
  }
  VlnPlot(seu_all, features = GBM_MARKERS, group.by = "sample_id", pt.size = 0, assay = "SCT")
  DotPlot(seu_all, features = GBM_MARKERS, group.by = "seurat_clusters", assay = "SCT")
  RidgePlot(seu_all, features = GBM_MARKERS, group.by = "sample_id", assay = "SCT")
  ```
  (`plots/05_markers/CD44_RNA_*`, `CD44_raw_count_by_cluster_sample.csv`는 대화형 R 콘솔에서
  애드혹으로 생성되었으며 스크립트는 보관되지 않았음 — 확인됨, 프로젝트 담당자 응답.)

### Step 9. Donor 분리 클러스터링/annotation (donor별 독립 경로, Harmony 미사용, BEST1 전용)
- **소프트웨어(버전):** R 4.5.3, Seurat 5.4.0, clustree 0.5.1
- **핵심 파라미터:** donor2/donor7 각각 within-donor merge 후 `SelectIntegrationFeatures
  (nfeatures=3000)`, `PrepSCTFindMarkers`, `RunPCA(npcs=50)`, `RunUMAP(reduction="pca",
  dims=1:30)`(Harmony 없음), `FindClusters(resolution=c(0.3,0.5,0.8))`; donor2 기본
  해상도 `SCT_snn_res.0.3`(7 clusters), donor7 기본 해상도 `SCT_snn_res.0.5`(11 clusters);
  세포유형 annotation은 Step 7과 동일한 8종 module score 방식 재사용
- **Input:** `rds/03_sct_seu_list.rds` (donor2: 2_SC+2_Best1, donor7: 7_SC+7_Best1)
- **Output:** `rds/donor{2,7}_clustered.rds`, `rds/donor{2,7}_annotated.rds`,
  `plots/donor{2,7}/{03_clustering,04_annotation,05_markers}/*`
- **실행 커맨드:**
  ```r
  features <- SelectIntegrationFeatures(object.list = seu_list_sct[samples], nfeatures = 3000)
  seu <- merge(seu_list_sct[[samples[1]]], y = seu_list_sct[[samples[2]]], add.cell.ids = samples)
  VariableFeatures(seu) <- features
  seu <- PrepSCTFindMarkers(seu)
  seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = c(0.3, 0.5, 0.8), verbose = FALSE)

  default_res <- if (donor_id == "donor2") "SCT_snn_res.0.3" else "SCT_snn_res.0.5"
  Idents(seu) <- default_res
  seu$seurat_clusters <- seu@meta.data[[default_res]]
  # 이후 Step 7과 동일한 8종 cell_signatures로 AddModuleScore → dominant cell_type 할당
  ```

### Step 10. BEST1 KD DEG 분석 (donor 내 비교) + fgsea (BEST1 전용)
- **소프트웨어(버전):** R 4.5.3, Seurat 5.4.0, msigdbr, fgsea(sessionInfo.txt에 정확한
  버전 별도 기재 없음 — 실행 시점 버전과 현재 환경 버전이 동일함은 확인됨), ggrepel
- **핵심 파라미터:** `FindMarkers(test.use="wilcox", min.pct=0.10, logfc.threshold=0)`;
  유의 DEG 기준 `p_val_adj < 0.05 & |avg_log2FC| > 0.25`; donor2 비교는
  `ident.1="2_Best1", ident.2="2_SC"`, donor7은 `ident.1="7_Best1", ident.2="7_SC"`;
  all-cell 비교 + 4개 세포유형(GBM_tumor, Astrocyte, Neuron, Proliferating)별 비교;
  cross-donor 일치 DEG는 두 donor에서 유의 + 같은 부호(`sign(lfc_d2)==sign(lfc_d7)`)인
  유전자; fgsea는 Hallmark(H) + KEGG_LEGACY(C2) gene set, `minSize=10, maxSize=500,
  nPermSimple=10000, seed=42`
- **Input:** `rds/05_annotated_seu.rds` (donor 통합 경로 Step 7 최종 산출물, RNA assay로
  전환 후 `NormalizeData` 재실행)
- **Output:** `plots/06_DEG/{BEST1_expression_baseline.csv, BEST1_VlnPlot.{pdf,png},
  DEG_summary_table.csv}`, `plots/06_DEG/{donor2,donor7}/{DEG_*.csv, Volcano_*.{pdf,png},
  fgsea_*.csv, fgsea_barplot_*.{pdf,png}}`, `plots/06_DEG/combined/{Consistent_DEG_*.csv,
  LFC_scatter_*.{pdf,png}, Heatmap_consistent_GBM_DEG.{pdf,png}}`
- **실행 커맨드:**
  ```r
  seu <- readRDS(file.path(RDS, "05_annotated_seu.rds"))
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, verbose = FALSE)

  run_deg <- function(obj, donor, ident_kd, ident_ctrl, cell_type_filter = NULL) {
    sub <- obj[, obj$donor_id == donor]
    if (!is.null(cell_type_filter)) sub <- sub[, sub$cell_type %in% cell_type_filter]
    Idents(sub) <- "sample_id"
    FindMarkers(sub, ident.1 = ident_kd, ident.2 = ident_ctrl,
                test.use = "wilcox", min.pct = 0.10, logfc.threshold = 0, verbose = FALSE)
  }
  deg_d2_all <- run_deg(seu, "donor2", "2_Best1", "2_SC")
  deg_d7_all <- run_deg(seu, "donor7", "7_Best1", "7_SC")

  common <- inner_join(d2_sig, d7_sig, by = "gene") %>%
    filter(sign(lfc_d2) == sign(lfc_d7))

  msig_h    <- msigdbr(species = "Homo sapiens", collection = "H")
  msig_kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  ranks <- sort(setNames(deg$avg_log2FC, deg$gene), decreasing = TRUE)
  set.seed(42)
  fgsea(pathways = gsets, stats = ranks, minSize = 10, maxSize = 500, nPermSimple = 10000)
  ```

### Step 11. QC / 정규화 / 클러스터링 (Tcell 전용, Scanpy)
- **소프트웨어(버전):** Python 3.13.0, Scanpy 1.12.1, AnnData 0.11.4, Scrublet 0.2.3, NumPy 2.2.6, Pandas 3.0.0, Matplotlib 3.10.8 (`/root/miniconda3` base 환경, 전용 conda env/environment.yml 없음)
- **핵심 파라미터:** `MIN_GENES=400, MAX_GENES=4500, MAX_PCT_MT=15%`; Scrublet `expected_doublet_rate=0.05, min_counts=2, min_cells=3, n_prin_comps=30`(샘플별 독립 실행); `normalize_total(target_sum=1e4)` + `log1p`; HVG `n_top_genes=2000`; `scale(max_value=10)`; PCA `n_comps=30, svd_solver="arpack"`; `neighbors(n_neighbors=30, n_pcs=30)`; UMAP `random_state=42`; Leiden `resolution=[0.3,0.5,0.8]`(`flavor="igraph", n_iterations=2`); 스레드 204(=256코어×80%, OpenBLAS 64 하드캡)
- **Input:** `cellranger_output/{C,T}_count_merged_v2/outs/filtered_feature_bc_matrix` (Step 1 산출물)
- **Output:** `adata/adata_preprocessed.h5ad`; `figures/01_QC/{QC_metrics_raw.png, QC_violin_postfilter.png, PCA_elbow.png}`; `figures/02_UMAP/{UMAP_overview.png, UMAP_UCB_markers.png, dotplot_leiden0.5_UCB_markers.png}`; `cache/*.h5ad`(sc.read_10x_mtx cache)
- **실행 커맨드:**
  ```bash
  # run_downstream.sh 내부
  python3 01_preprocess_cluster.py 2>&1 | tee logs/01_preprocess.log
  ```
  ```python
  # 01_preprocess_cluster.py 핵심 로직 발췌
  a = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
  adata = ad.concat({"C": a_C, "T": a_T}, label="sample")
  sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True)

  scrub = scr.Scrublet(mat, expected_doublet_rate=0.05)
  scores, doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, n_prin_comps=30, verbose=False)

  keep = ((adata.obs["n_genes_by_counts"] >= 400) & (adata.obs["n_genes_by_counts"] <= 4500) &
          (adata.obs["pct_counts_mt"] <= 15.0) & (~adata.obs["predicted_doublet"]))
  adata = adata[keep].copy()

  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)
  adata.raw = adata
  sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False)
  sc.pp.scale(adata, max_value=10)
  sc.tl.pca(adata, n_comps=30, use_highly_variable=True, svd_solver="arpack")
  sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, use_rep="X_pca")
  sc.tl.umap(adata, random_state=42)
  for res in [0.3, 0.5, 0.8]:
      sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2)
  adata.write_h5ad(f"{ADATA_DIR}/adata_preprocessed.h5ad")
  ```

### Step 12. T세포 서브셋 식별 + NgR1(RTN4R)/NogoA(RTN4) DEG (Tcell 전용)
- **소프트웨어(버전):** Scanpy 1.12.1 (`sc.tl.rank_genes_groups`, Wilcoxon), Step 11과 동일 환경
- **핵심 파라미터:** T세포 판정 `TCELL_SCORE_THRESHOLD=0.3`(CD3D/E/G `score_genes` 평균, `leiden_0.5` 클러스터 단위); 적혈구 오염 제거 `HB_SCORE_THRESHOLD=0.5`(HBB/HBA1/HBA2/HBG1/HBG2); 재처리 `N_HVG_TCELL=2000, N_PCS_TCELL=30, N_NEIGHBORS_T=20`, Leiden `[0.3,0.5,0.8]`; NgR1/NogoA 양성 정의 = log1p 발현 `> 0`; DEG 기준 `padj<0.05, |LFC|>=0.5, pct_nz_group>=0.10`(Wilcoxon, `use_raw=True, pts=True`)
- **Input:** `adata/adata_preprocessed.h5ad` (Step 11 산출물)
- **Output:** `adata/adata_Tcell.h5ad`; `tables/{cluster_Tcell_scores.csv, DEG_NgR1_NogoA_pos_vs_neg_{all,filtered}.csv, DEG_RTN4RposRTN4neg_vs_doubleNeg.csv, DEG_RTN4RnegRTN4pos_vs_doubleNeg.csv, DEG_RTN4RposRTN4pos_vs_doubleNeg.csv, DEG_Treatment_vs_Control_in_NgR1pos.csv, composition_NgR1_NogoA_{counts,pct}.csv}`; `figures/03_Tcell/{UMAP_Tcell_identification.png, UMAP_Tcell_markers.png, dotplot_Tcell_subtype_markers.png, heatmap_UCB_hallmarks_per_cluster.png}`; `figures/04_NgR1_NogoA/*`(UMAP/violin/volcano/heatmap/dotplot/barplot 8종)
- **실행 커맨드:**
  ```bash
  # run_downstream.sh 내부
  python3 02_Tcell_NgR1_NogoA.py 2>&1 | tee logs/02_NgR1_NogoA.log
  ```
  ```python
  # 02_Tcell_NgR1_NogoA.py 핵심 로직 발췌
  sc.tl.score_genes(adata, gene_list=["CD3D","CD3E","CD3G"], score_name="Tcell_score", use_raw=True)
  tcell_clusters = cluster_scores[cluster_scores["Tcell_mean"] >= 0.3].index.tolist()
  adata_t = adata[adata.obs["leiden_0.5"].isin(tcell_clusters)].raw.to_adata()

  sc.tl.score_genes(adata_t, gene_list=hb_genes, score_name="hb_score")
  adata_t = adata_t[adata_t.obs["hb_score"] < 0.5].copy()

  adata_t.obs["NgR1_pos"]  = (expr_RTN4R > 0)
  adata_t.obs["NogoA_pos"] = (expr_RTN4  > 0)
  # → NgR1_NogoA_group ∈ {RTN4R+RTN4+, RTN4R+RTN4-, RTN4R-RTN4+, RTN4R-RTN4-}

  sc.tl.rank_genes_groups(adata_t, groupby="pos_any", groups=["Positive"], reference="Negative",
                          method="wilcoxon", use_raw=True, pts=True, key_added="deg_NgR1_NogoA")
  sig = deg_df[(deg_df["pvals_adj"] < 0.05) & (deg_df["logfoldchanges"].abs() >= 0.5) &
              (deg_df["pct_nz_group"] >= 0.10)]

  # Treatment vs Control, NgR1+/NogoA+ 세포 내부
  sc.tl.rank_genes_groups(pos_adata, groupby="sample", groups=["T"], reference="C",
                          method="wilcoxon", use_raw=True, pts=True)
  adata_t.write_h5ad(f"{ADATA_DIR}/adata_Tcell.h5ad")
  ```

### Step 13. Cell type annotation (Tcell 전용)
- **소프트웨어(버전):** Scanpy 1.12.1, Step 11/12와 동일 환경
- **핵심 파라미터:** Step 12와 동일한 T세포 재추출 로직(`TCELL_SCORE_THRESHOLD=0.3`, `HB_SCORE_THRESHOLD=0.5`, `random_state=42`로 동일 클러스터 재현) 재실행 후, `CELL_TYPE_MAP`으로 `leiden_0.5` 클러스터 0~9를 10개 라벨(0→Naive T, 1/2→Treg-like, 3→γδ T, 4/5→Naive T, 6/7→Activated T, 8→Treg, 9→Treg-like)로 수동 매핑; C/T 조성은 `cell_type` 컬럼과 분리해 별도 테이블로 산출
- **Input:** `adata/adata_preprocessed.h5ad` (Step 11 산출물; run_downstream.sh에는 미포함, 단독 실행)
- **Output:** `adata/adata_Tcell.h5ad`(`cell_type` 컬럼 추가, Step 12 산출물을 덮어씀), `tables/composition_celltype_CvsT.csv`, `figures/03_Tcell/barplot_celltype_composition_CvsT.png`
- **실행 커맨드:**
  ```bash
  python3 03_add_celltype_annotation.py
  ```
  ```python
  # 03_add_celltype_annotation.py 핵심 로직 발췌
  CELL_TYPE_MAP = {
      '0': 'Naive T', '1': 'Treg-like', '2': 'Treg-like', '3': 'γδ T', '4': 'Naive T',
      '5': 'Naive T', '6': 'Activated T', '7': 'Activated T', '8': 'Treg', '9': 'Treg-like',
  }
  adata_t.obs["cell_type"] = adata_t.obs["leiden_0.5"].map(CELL_TYPE_MAP).astype("category")

  counts = adata_t.obs.groupby(["cell_type", "sample"]).size().unstack(fill_value=0)
  comp["C_pct_of_celltype"] = pct_within_celltype["C"].round(1)
  comp.to_csv(f"{TABLE_DIR}/composition_celltype_CvsT.csv")
  adata_t.write_h5ad(f"{ADATA_DIR}/adata_Tcell.h5ad")
  ```

### Step 14. CPM 매트릭스 생성 (Tcell 전용, 병렬 산출물)
- **소프트웨어(버전):** Scanpy 1.12.1, AnnData 0.11.4, `scipy.io.mmwrite`, Step 11과 동일 환경
- **핵심 파라미터:** `normalize_total(target_sum=1e6)`(CPM) 후 `log1p()`; QC-passed 바코드(Step 11 결과)만 유지, raw UMI count는 Cell Ranger filtered matrix에서 재로딩(`backed="r"`로 QC obs만 우선 로드)
- **Input:** `adata/adata_preprocessed.h5ad`(QC-passed 바코드 목록, backed="r") + `cellranger_output/{C,T}_count_merged_v2/outs/filtered_feature_bc_matrix`(raw UMI count 재로딩)
- **Output:** `cpm/adata_cpm.h5ad`(layers: counts/cpm/log1p_cpm), `cpm/mtx/cpm/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}`, `cpm/mtx/log1p_cpm/{...}`
- **실행 커맨드:**
  ```bash
  python3 05_make_cpm_matrix.py 2>&1 | tee logs/05_make_cpm.log
  ```
  ```python
  # 05_make_cpm_matrix.py 핵심 로직 발췌
  qc = sc.read_h5ad(f"{ADATA_DIR}/adata_preprocessed.h5ad", backed="r")
  adata = ad.concat({"C": raw_C, "T": raw_T}, label="_concat_sample")
  adata = adata[qc.obs.index].copy()          # QC-passed 바코드만 유지
  adata.layers["counts"] = adata.X.copy()

  cpm = adata.copy()
  sc.pp.normalize_total(cpm, target_sum=1e6)
  adata.layers["cpm"] = cpm.X.copy()
  log1p_cpm = cpm.copy(); sc.pp.log1p(log1p_cpm)
  adata.layers["log1p_cpm"] = log1p_cpm.X.copy()
  adata.X = adata.layers["cpm"]
  adata.write_h5ad(f"{CPM_DIR}/adata_cpm.h5ad")
  mmwrite(mtx_path, adata.layers["cpm"].T.tocoo())   # genes x cells, 10x 관례
  ```

### Step 15. HTML 종합 리포트 생성 (Tcell 전용)
- **소프트웨어(버전):** 순수 Python 표준 라이브러리(base64, os)만 사용, 외부 패키지 의존 없음
- **핵심 파라미터:** 없음(정적 HTML/CSS 생성); 수치 테이블은 동적 로드가 아닌 Python 코드 내 하드코딩된 리터럴 값(Step 11~14 산출물에서 수동 전사한 것으로 추정 — 향후 산출물 갱신 시 리포트에 자동 반영되지 않음)
- **Input:** `figures/{01_QC,02_UMAP,03_Tcell,04_NgR1_NogoA}` 하위 PNG(base64 임베드)
- **Output:** `downstream_analysis/analysis_report.html`
- **실행 커맨드:**
  ```bash
  python3 04_generate_report.py
  ```

### Step 16. PPT 슬라이드 생성 (Tcell 전용)
- **소프트웨어(버전):** python-pptx (`pptx.Presentation`)
- **핵심 파라미터:** 슬라이드 크기 10.0×5.625 in; 그림 파일은 직접 임베드하지 않고 플레이스홀더 박스 + 경로 텍스트로만 참조(`figure_placeholder()` 함수); 수치는 04_generate_report.py와 동일하게 하드코딩된 Python 리터럴; 디자인 시스템은 저장소 내 참고 문서 기준(폰트 "Pretendard" 지정, 설치 여부 미확인)
- **Input:** Step 15와 동일한 하드코딩 수치/경로 문자열(그림 실체는 임베드하지 않음)
- **Output:** `CNU_Tcell_NgR1_NogoA_Slides.pptx`
- **실행 커맨드:**
  ```bash
  python3 generate_slides.py
  ```

### Step 17. Doublet 검출 교차검증 사이드 브랜치 (Tcell 전용, 비블로킹·메인 파이프라인 미반영)
- **소프트웨어(버전):** Python 쪽 Scrublet 0.2.3(Step 11과 동일); R 쪽 R 4.2.0 + Seurat 5.1.0 + DoubletFinder 2.0.6 + Matrix 1.6.5(conda env `seurat_env`, 직접 조회로 확인)
- **핵심 파라미터:** 기대 multiplet rate `C=0.15, T=0.18`(10x loading guide 기반, 회수 세포수 ~19,943/~23,185 기준); DoubletFinder `pN=0.25`, `pK`는 `paramSweep`+`find.pK`로 자동 산출, `nExp=round(rate×ncol)`; Seurat 5 비호환으로 `doubletFinder()` 원본이 all-NA를 반환하는 문제를 `doubletFinder_fixed()`(pANN을 named numeric vector로 재구현)로 우회
- **Input:** `cellranger_output/{C,T}_count_merged_v2/outs/filtered_feature_bc_matrix` (Step 1과 동일 raw 입력, Step 11 Scrublet과 교차비교 목적)
- **Output:** `figures/01_QC/{scrublet_histogram_{C,T}.png, doubletfinder_histogram_{C,T}.png}`; `tables/{doubletfinder_{C,T}_results.csv, doubletfinder_summary.csv}`
- **실행 커맨드:**
  ```bash
  python3 check_doublet_histogram.py
  # Seurat5 호환성 디버깅 반복 (debug_df.R → debug_df2.R → debug_df3.R) 후 최종본 실행:
  Rscript run_doubletfinder.R
  ```
  ```r
  # run_doubletfinder.R 핵심 로직 발췌
  counts <- Read10X(data.dir = samples[[sid]])
  seu <- CreateSeuratObject(counts = counts, project = sid, min.cells = 3, min.features = 200)
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA(npcs = 30)

  sweep.res  <- paramSweep(seu, PCs = 1:30, sct = FALSE)
  bcmvn      <- find.pK(summarizeSweep(sweep.res, GT = FALSE))
  pK_optimal <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  nExp_poi   <- round(expected_rate[[sid]] * ncol(seu))   # C: 0.15, T: 0.18

  seu <- doubletFinder_fixed(seu, PCs = 1:30, pN = 0.25, pK = pK_optimal, nExp = nExp_poi)
  # 결과: C=15.0%, T와 유사 수준 — Scrublet의 "0 doublets" 판정과 불일치하나 메인 파이프라인에는 미반영
  ```
  (이 결과는 adata_preprocessed.h5ad의 필터링에 피드백되지 않았음이 소스 문서에서 확인됨 —
  최종 파이프라인은 Step 11의 Scrublet 기반 판정을 그대로 사용한다.)

## 2. 프로젝트별 특이사항

### Homo_sapiens_BEST1_scRNAseq
- 연구목적/샘플 설계: Human brain organoid 세포주(#2, #7)에서 BEST1 knockdown 시 단일세포
  수준 전사체 프로파일링. 10X Genomics 4-plex OCM 3' 키트로 하나의 풀링 GEX 라이브러리를
  제작해 마크로젠에서 시퀀싱, OCM 바코드(OB1~OB4)로 4개 샘플(2_SC, 2_Best1, 7_SC, 7_Best1)을
  디멀티플렉싱.
- 공통 파이프라인과의 차이점: Step 1(Cell Ranger) 이후 downstream 전체가 R/Seurat 생태계로
  진행되며(Tcell 프로젝트의 Python/Scanpy와 무관), Seurat 기반 QC → DoubletFinder 이중세포
  제거 → SCTransform 정규화 → (donor 통합 Harmony 경로 / donor 분리 경로) → 클러스터링 →
  module score 기반 세포유형 annotation → BEST1 KD 관련 DEG(도너 내 비교, Wilcoxon)+fgsea
  분석까지 이어지는 프로젝트 고유의 Step 2~10 계열을 갖는다(위 1-1 DAG 참고).
- 프로젝트 내부에 두 개의 병렬 분기가 존재: (1) Step 5~8의 Harmony 기반 donor 통합 분석과
  (2) Step 9의 Harmony 없는 donor별 독립(within-donor) 분석. 두 경로 모두 동일한 `03_sct_seu_list.rds`
  (SCTransform 완료 시점)에서 분기하며, donor 간 GBM 특성 차이가 크다는 판단 하에 별도로
  수행되었다. donor2는 해상도 0.3(소규모, ~3,373 cells 추정), donor7은 해상도 0.5(~4,408 cells
  추정)를 기본값으로 사용.
- DEG 분석(Step 10)은 donor 통합 경로의 최종 annotation 결과(`05_annotated_seu.rds`)를
  입력으로 사용하며, donor 분리 경로(Step 9)의 결과와는 별개로 진행됨 — 즉 donor 내
  KD-vs-Control 비교의 cell_type 라벨은 통합 annotation 기준.
- 소스 문서(`02_pipeline_steps.md`)에서 "확인 필요"로 남긴 항목은 그대로 유지한다:
  `01_clustering_annotation.R` 스크립트 파일 자체의 재수정 시점(2026-07-01)과 최초
  실행판(2026-06-22 산출물 생성 시점) 간의 차이 — 어떤 부분이 수정되었는지(주석/경로
  정리 등 non-functional change로 추정) 소스 diff 확인 불가(버전 관리 시스템 없음, git repo
  아님).
- `plots/05_markers/CD44_RNA_*`, `CD44_raw_count_by_cluster_sample.csv`는 대화형 R 콘솔에서
  애드혹으로 생성되었으며 스크립트는 보관되지 않았음(확인됨, 프로젝트 담당자 응답) — 해당
  산출물에 대응하는 재현 가능한 실행 커맨드는 확인 필요 — 스크립트 파일 없음.
- 별도의 conda `environment.yml`/`requirements.txt`/`config.yaml` 등 환경 설정 파일은
  발견되지 않음 — 필요 시 사용된 R 패키지 버전은 `analysis/sessionInfo.txt`에서 확인 가능
  (R 4.5.3, Seurat 5.4.0, SeuratObject 5.1.0, harmony 1.2.3, DoubletFinder 2.0.6, clustree
  0.5.1, ggplot2 4.0.3, dplyr 1.2.0, patchwork 1.3.2; Ubuntu 24.04.4 LTS).
- `metadata/cellranger.sh`의 정확한 최종 수정일은 파일 자체 mtime 미확인 상태이며 config
  CSV 생성일(2026-06-21)로 추정한 값임 — 확인 필요.

### Homo_sapiens_Tcell_scRNAseq
- 연구목적/샘플 설계: 인간(Homo sapiens) 제대혈(cord blood) T세포에 anti-CD3/anti-CD28
  처리로 activation을 유도한 scRNA-seq 실험. 샘플은 C(Control)/T(Treatment) 2개이며, 동일
  라이브러리를 3회에 걸쳐 재시퀀싱한 배치를 병합해 사용한다(본 문서는 `2nd_generation` 범위만
  다루며, 상위 1차 프로젝트의 `rawdata`/`cellranger_output`/`reference`/
  `run_cellranger_analysis.sh`는 담당자 지시에 따라 스코프 외로 제외됨).
- 공통 파이프라인과의 차이점: Step 1(Cell Ranger) 이후 downstream 전체가 Python/Scanpy
  생태계로 진행되며(BEST1 프로젝트의 R/Seurat과 무관), Harmony 등 배치보정 기법을 전혀
  사용하지 않는다 — 3회의 시퀀싱 배치가 Cell Ranger BAM 레벨에서 이미 병합되었고 동일
  라이브러리의 재시퀀싱이므로 기술적 배치효과가 없다고 판단, C vs T 차이는 순수 생물학적
  차이로 취급함(01_preprocess_cluster.py 주석에 명시). Scanpy 기반 QC(min_genes=400,
  max_genes=4500, max_pct_mt=15%) → Scrublet 이중세포 검출(DoubletFinder 아님) → HVG/PCA/
  UMAP/Leiden 클러스터링 → T세포 서브셋 재추출 → cell type annotation → CPM 매트릭스
  생성까지 이어지는 프로젝트 고유의 Step 11~17 계열을 갖는다(위 1-2 DAG 참고).
- 3-배치 심볼릭 링크 lane 위장 트릭(Cell Ranger 고유 확장): Batch1(L009, 2025-12-26 수신)은
  그대로, Batch2(`rawdata/Additional/`, 원본 파일명도 L009)는 `rawdata_Additional_L010/`
  경로에 L010으로 리네임한 심볼릭 링크를, Batch3(`rawdata/Additional_v2/`, 원본도 L009)는
  `rawdata_Additional_v2_L011/`에 L011로 리네임한 심볼릭 링크를 생성해 Cell Ranger가 3개
  배치를 동일 샘플의 서로 다른 lane으로 인식하고 병합하도록 함. 세대별 산출물
  `C_count`/`T_count`(Batch1만) → `C_count_merged`/`T_count_merged`(Batch1+2) →
  `C_count_merged_v2`/`T_count_merged_v2`(Batch1+2+3, ★최종★)가 순차 생성되었고, 구버전
  2세대(총 ~70G)는 2026-07-08 인수인계 정리 작업에서 삭제, 최종본만 보존됨.
- 프로젝트 고유 확장 분석: NgR1(RTN4R)/NogoA(RTN4) 관련 T세포 서브셋 분석(02_Tcell_NgR1_NogoA.py) —
  CD3D/E/G 기반 T세포 스코어(클러스터 평균≥0.3, leiden_0.5)로 T세포 클러스터를 판정한 뒤
  Hb 스코어(≥0.5)로 적혈구 오염을 제거하고, RTN4R/RTN4 발현(>0)을 기준으로 4개 그룹
  (RTN4R+RTN4+/RTN4R+RTN4-/RTN4R-RTN4+/RTN4R-RTN4-)으로 분류, Wilcoxon 기반 DEG(Positive
  vs Double-negative; Treatment vs Control in Positive)를 수행한다. 이 서브셋 분석은
  BEST1 프로젝트에는 없는 Tcell 프로젝트 고유 확장이다.
- CPM 매트릭스 생성(05_make_cpm_matrix.py)도 Tcell 프로젝트 고유 확장이다 — `adata_preprocessed.h5ad`는
  Script 1에서 `normalize_total(target_sum=1e4)+log1p`로 X를 in-place 덮어쓰므로, raw UMI
  count는 Cell Ranger filtered matrix에서 재로딩한 뒤 QC-passed 바코드만 남기고
  `normalize_total(target_sum=1e6)`(CPM)+log1p 처리하여 별도 저장한다.
- QC 교차검증 사이드 브랜치(메인 파이프라인 비반영): `check_doublet_histogram.py`로 Scrublet
  결과를 재시각화하고, `debug_df.R`→`debug_df2.R`→`debug_df3.R`(Seurat 5 비호환성 디버깅)를
  거쳐 최종 `run_doubletfinder.R`(conda `seurat_env`, R 4.2.0 + Seurat 5.1.0 + DoubletFinder
  2.0.6 + Matrix 1.6.5)로 DoubletFinder 기반 doublet률(C=15.0%, T=18% 기대치 기반)을
  산출했으나, 이 결과가 실제로 `adata_preprocessed.h5ad`의 필터링에 반영되었다는 증거는 없다 —
  최종 파이프라인은 Scrublet 기반 "0 doublets" 판정을 그대로 사용한다(소스 문서에 명시된
  내용, 확인됨).
- 세포유형 annotation(`03_add_celltype_annotation.py`)의 `CELL_TYPE_MAP`은 leiden_0.5
  클러스터 0~9를 Naive T/Treg-like/γδ T/Activated T/Treg 등 10개 라벨로 매핑하는 수동
  마커기반 annotation이며, C/T 조건은 `cell_type`에 인코딩하지 않고 별도 composition
  테이블(`composition_celltype_CvsT.csv`)로 분리해 관리한다(소스 문서 코드 주석에 명시된
  설계 의도).
- HTML 리포트(`04_generate_report.py`)와 PPT(`ppt/generate_slides.py`)의 수치는 동적 로드가
  아닌 Python 코드 내 하드코딩된 리터럴 값으로, Step 11~14 산출물에서 수동 전사한 것으로
  추정된다 — 소스 문서에서도 "확인 필요"로 명확히 결론짓지 못한 부분이며, 향후 산출물이
  갱신되어도 리포트/PPT가 자동 반영되지 않는다는 재현성 caveat이 있다.
- 별도의 conda `environment.yml`/`requirements.txt` 등 환경 설정 파일은 downstream 파이프라인
  자체(01/02/03/05 스크립트, `/root/miniconda3` base 환경)에는 존재하지 않음 — 확인된 실제
  설치 버전은 Python 3.13.0, Scanpy 1.12.1, AnnData 0.11.4, Scrublet 0.2.3, NumPy 2.2.6,
  Pandas 3.0.0, Matplotlib 3.10.8 (소스 문서에서 직접 조회로 확인, `04_generate_report.py`
  기재 버전과 일치). DoubletFinder 교차검증 사이드 브랜치만 별도 conda env(`seurat_env`)를
  사용.
- `reference/GRCh38_cellranger_index` 심볼릭 링크(1차 프로젝트 인덱스를 가리킴)는 현재
  삭제된 상태로, 본 문서 스코프(2nd_generation) 내에서 재현하려면 참조 인덱스를 별도로
  확보해야 함 — 1차 프로젝트 관련 재구성 방법은 스코프 외로 처리(담당자 지시).

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Homo_sapiens_BEST1_scRNAseq | 10X OCM 4-plex 단일 GEX 라이브러리를 cellranger multi 1회 실행으로 4개 샘플 디멀티플렉싱 후, Seurat 파이프라인이 Harmony 기반 donor 통합 경로와 Harmony 미사용 donor 분리 경로로 병행 분기되는 scRNA-seq 프로젝트 | Cell Ranger 10.0.0; R 4.5.3 단일 환경(Seurat 5.4.0, harmony 1.2.3, DoubletFinder 2.0.6); 별도 conda env 미사용 | Module score 기반 8종 세포유형 자동 annotation(GBM_tumor 등 뇌종양 organoid 특이 signature), donor 내 BEST1 KD vs Control DEG(Wilcoxon) + cross-donor 방향 일치 DEG 필터링 + fgsea(Hallmark/KEGG) 통합 분석 |
| Homo_sapiens_Tcell_scRNAseq | 동일 라이브러리를 3회 재시퀀싱한 배치를 심볼릭 링크 lane 위장 트릭으로 병합(cellranger count, Harmony 등 배치보정 미사용) 후, Scanpy 기반 QC/Scrublet 이중세포 검출 → 클러스터링 → T세포 서브셋 재추출 → cell type annotation → CPM 매트릭스 생성으로 이어지는 Python 생태계 파이프라인 | Cell Ranger 10.0.0; Python 3.13.0 단일 환경(Scanpy 1.12.1, AnnData 0.11.4, Scrublet 0.2.3); DoubletFinder 교차검증 사이드 브랜치만 별도 conda env(`seurat_env`, R 4.2.0+Seurat 5.1.0) 사용 | NgR1(RTN4R)/NogoA(RTN4) 발현 기반 4-group T세포 서브셋 분류 + Wilcoxon DEG(Positive vs Double-negative, Treatment vs Control) 분석; QC-passed 바코드 기준 CPM/log1p-CPM 매트릭스 별도 생성(05_make_cpm_matrix.py); Scrublet vs DoubletFinder 이중세포 교차검증 사이드 브랜치(메인 파이프라인 비반영) |
</content>
