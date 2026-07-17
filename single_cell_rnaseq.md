# Single-cell RNA-seq (Homo_sapiens_BEST1_scRNAseq)

## 1. 공통 파이프라인

본 프로젝트는 단일 프로젝트 카테고리이며, Human brain organoid 세포주(#2, #7)에서 BEST1
knockdown 시 단일세포 전사체 변화를 프로파일링한다. 10X Genomics GEM-X OCM(On-Chip
Multiplexing) 4-plex 3' 키트로 하나의 풀링된 GEX 라이브러리를 제작해 시퀀싱했고,
cellranger multi 1회 실행으로 OCM 바코드(OB1~OB4)를 이용해 4개 샘플(2_SC, 2_Best1, 7_SC,
7_Best1)을 동시에 디멀티플렉싱한다. 이후 Seurat 기반 QC → DoubletFinder 이중세포 제거 →
SCTransform 정규화 → (donor 통합 Harmony 경로 / donor 분리 경로) → 클러스터링 → module
score 기반 세포유형 annotation → BEST1 KD 관련 DEG(도너 내 비교, Wilcoxon)+fgsea 분석까지
이어진다. donor 2와 donor 7은 GBM 특성 차이가 크다는 판단 하에, Harmony로 통합한 전체 분석과
Harmony 없이 donor별로 독립 처리한 분석이 병행 수행되었다.

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

### Step 1. CellRanger multi (OCM 4-plex 디멀티플렉싱 + 정렬/카운팅)
- **소프트웨어(버전):** Cell Ranger 10.0.0 (시스템 PATH 직접 설치, conda/컨테이너 미사용;
  OCM 기능은 cellranger v9.0+ 필요, 본 프로젝트는 v10.0 사용)
- **핵심 파라미터:** `--localcores=150 --localmem=500`(GB); config CSV 내
  `create-bam=true`, `include-introns=true`, reference=GRCh38-2024-A(10x 프리빌트);
  OCM 바코드 매핑 2_SC=OB1, 2_Best1=OB2, 7_SC=OB3, 7_Best1=OB4
- **Input:** `rawdata/HN00277842_10X_RawData_Outs/Best1/23KWMHLT4/Best1_S6_L004_{R1,R2,I1,I2}_001.fastq.gz`
  (4개 샘플이 칩 로딩 단계에서 이미 풀링된 단일 GEX 라이브러리), `reference/refdata-gex-GRCh38-2024-A`
- **Output:** `metadata/cellranger_output/BrainOrganoid_OCM_4plex/outs/{raw_feature_bc_matrix*,
  filtered_feature_bc_matrix*, qc_report.html, per_sample_outs/{2_SC,2_Best1,7_SC,7_Best1}/outs/
  sample_filtered_feature_bc_matrix.h5}`
- **실행 커맨드:**
  ```bash
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

### Step 2. Seurat 객체 생성 및 샘플별 QC 필터링
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

### Step 3. DoubletFinder 이중세포 제거 (샘플별 독립 실행)
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

### Step 4. SCTransform 정규화 (샘플별)
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

### Step 5. Merge + Harmony 통합 (donor 통합 경로)
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

### Step 6. UMAP + 클러스터링 (donor 통합 경로)
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

### Step 7. Module Score 기반 세포유형 Annotation (donor 통합 경로)
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

### Step 8. GBM 마커 유전자 시각화 (donor 통합 경로)
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

### Step 9. Donor 분리 클러스터링/annotation (donor별 독립 경로, Harmony 미사용)
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

### Step 10. BEST1 KD DEG 분석 (donor 내 비교) + fgsea
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

## 2. 프로젝트별 특이사항

### Homo_sapiens_BEST1_scRNAseq
- (단일 프로젝트 카테고리 — 위 공통 파이프라인이 곧 이 프로젝트의 파이프라인)
- 연구목적/샘플 설계: Human brain organoid 세포주(#2, #7)에서 BEST1 knockdown 시 단일세포
  수준 전사체 프로파일링. 10X Genomics 4-plex OCM 3' 키트로 하나의 풀링 GEX 라이브러리를
  제작해 마크로젠에서 시퀀싱, OCM 바코드(OB1~OB4)로 4개 샘플(2_SC, 2_Best1, 7_SC, 7_Best1)을
  디멀티플렉싱.
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

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Homo_sapiens_BEST1_scRNAseq | 10X OCM 4-plex 단일 GEX 라이브러리를 cellranger multi 1회 실행으로 4개 샘플 디멀티플렉싱 후, Seurat 파이프라인이 Harmony 기반 donor 통합 경로와 Harmony 미사용 donor 분리 경로로 병행 분기되는 유일한 scRNA-seq 프로젝트 | Cell Ranger 10.0.0; R 4.5.3 단일 환경(Seurat 5.4.0, harmony 1.2.3, DoubletFinder 2.0.6); 별도 conda env 미사용 | Module score 기반 8종 세포유형 자동 annotation(GBM_tumor 등 뇌종양 organoid 특이 signature), donor 내 BEST1 KD vs Control DEG(Wilcoxon) + cross-donor 방향 일치 DEG 필터링 + fgsea(Hallmark/KEGG) 통합 분석 |
