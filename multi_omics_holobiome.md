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

| 스텝 | 소프트웨어(버전) | 핵심 파라미터 | Input | Output |
|---|---|---|---|---|
| 0. 원시데이터 QC (mRNA-seq) | fastp(버전 확인 필요), MultiQC | `--cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --dedup --thread 256` | `*_R1_*/*_R2_*.fastq.gz`(Illumina raw, 30샘플) | `clean_data/*_clean_R{1,2}.fastq.gz`, `qc_reports/*_fastp.{html,json}`, `multiqc_result/` |
| 업스트림 16S ASV 생성 | QIIME2 2024.10 amplicon (`ginseng_qiime2` env), PacBio full-length pb16S nextflow | 확인 필요(외부 서버 처리, 상세 파라미터 미확인) | PacBio HiFi 16S raw reads | `feature-table-tax.biom`, `best_taxonomy.tsv`, `metadata.tsv` (rawdata 하위 `results_ginseng_v2/`에 저장) |
| 1. 16S ASV 필터링 | R(phyloseq, biomformat, dplyr), `ginseng_r` env | Bacteria/Archaea만 유지, Mitochondria/Chloroplast 제거 | `feature-table-tax.biom`, `best_taxonomy.tsv`, `metadata.tsv` | `ps_filtered.rds`, `asv_prevalence_abundance.tsv`, `library_sizes.pdf` |
| 2. CLR/ANCOM-BC2/Tier 분류 | R(compositions CLR, zCompositions CZM, ANCOMBC ancombc2), `ginseng_r` | `fix_formula="year"`, global+trend+pairwise, BH, n_cl=4; Tier2 게이트 prevalence≥50%(3회 재실행 후 최종 확정) | `ps_filtered.rds` | `otu_clr.rds`, `ancombc2_result.rds`, `asv_tier_classification.tsv`, `ps_selected.rds`, `otu_clr_selected.rds` |
| 3. 16S diversity baseline | R(vegan, phyloseq, ggpubr), `ginseng_r` | rarefy depth=32,957, Aitchison(CLR+Euclidean) PRIMARY, PERMANOVA 999 permutations | `ps_filtered.rds`, `ps_selected.rds`, `otu_clr.rds`, `otu_clr_selected.rds` | `alpha_diversity.csv/.pdf`, `pcoa_aitchison*.pdf`, `permanova_year*.txt` 등 |
| 4. HISAT2 인덱스 구축 | HISAT2(hisat2-build, 버전 확인 필요), 시스템 경로(conda env 미명시) | N_THREADS=38 | 참조 게놈 fasta/gtf | `index/Pg_hisat2*.ht2`, `index/splice_sites.txt`, `index/exons.txt` |
| 5. HISAT2 정렬 | HISAT2 + samtools(버전 확인 필요) | `--rna-strandness RF`, N_THREADS=38 | `clean_data/*_clean_R{1,2}.fastq.gz`, HISAT2 인덱스 | `outputs/bam/*.sorted.bam(+.bai,.flagstat.txt)` |
| 6. StringTie2 어셈블리/정량 | StringTie2(`stringtie`, 버전 확인 필요) | `-eB --rf -G <GFF3>`, N_THREADS=38 | `outputs/bam/*.sorted.bam` | `outputs/stringtie/{sample}/{transcripts.gtf, gene_abundances.tsv}` |
| 7. 카운트 매트릭스 생성 | `prepDE.py3`(외부 경로 `smJeon` 소유, 버전 확인 필요) | `-l 151`(read length) | StringTie2 transcripts.gtf | `gene_count_matrix.csv`, `transcript_count_matrix.csv` (63,713행) |
| 8. edgeR 정규화/DE | edgeR(TMM, QLF), `ginseng_r` | CPM≥1 in ≥3 samples, TMM, `~year` 모델, FDR<0.05 & \|logFC\|≥1 | `gene_count_matrix.csv` | `dge.rds`, `logcpm.rds`, `DE_{3Y_vs_2Y,4Y_vs_2Y,4Y_vs_3Y}.tsv` |
| 9. WGCNA 모듈 분석 | WGCNA(R, `blockwiseModules`), `ginseng_r` | top 15,000 variable genes, `networkType="signed"`, power 자동선택(R²≥0.85), minModuleSize=30, mergeCutHeight=0.25 | `dge.rds`, `logcpm.rds` | `gene_module_assignment.tsv`, `module_eigengenes.csv`, `hub_genes.tsv` 등 (13개 모듈) |
| 10. eggNOG 기능주석 | eggNOG-mapper(≥2.1.9) + DIAMOND(sensitive, --iterate), `ginseng_eggnog` env | `--cpu 38 --dbmem --go_evidence non-electronic`, eggNOG DB v5.0.2 | 외부 경로 단백질서열(.pep, 프로젝트 외부 참조) | `Pg_eggnog.emapper.annotations` |
| 11. 주석 테이블 구축 | R(GO.db, KEGGREST, AnnotationDbi), `ginseng_r` | KEGG 식물 15종 화이트리스트로 map##### 필터링 | `Pg_eggnog.emapper.annotations`, `logcpm.rds` | `gene2go.tsv`, `gene2ko.tsv`, `term2gene_*.tsv`, `term2name_*.tsv` |
| 12. GO/KEGG enrichment | clusterProfiler::enricher/compareCluster, enrichplot, `ginseng_enrich` env | DEG: FDR<0.05 & \|logFC\|≥1, enrichment p<0.05/q<0.2, minGSSize=5/maxGSSize=500; 9 DEG set + 13 WGCNA 모듈 | `term2gene_*`/`term2name_*`, `DE_*.tsv`, `gene_module_assignment.tsv` | `enrichment/{deg,wgcna}/*.tsv`, `summary_top_terms.tsv` |
| 13. ANCOM-BC2 시각화 | R(phyloseq, pheatmap), `ginseng_r` | q<0.05 유의 ASV, top/bottom 40 표시 | `ancombc2_result.rds`, `ps_filtered.rds` | `lollipop_*.pdf`, `fc_heatmap.pdf` |
| 14. SPIEC-EASI 네트워크 | SpiecEasi(GitHub `zdk123/SpiecEasi`, method="mb") + igraph/ggraph, `ginseng_spieceasi`(R 4.5.3) | pooled: Tier1+Tier2, per-year(n=10): Tier1만(탐색적); `pulsar.params(rep.num=20, seed=42)` | `ps_filtered.rds`, `asv_tier_classification.tsv` | `spieceasi_pooled.rds`, `spieceasi_{2Y,3Y,4Y}.rds`, `spieceasi_hubs.tsv` |
| 15. SparCC 네트워크 | SparCC(SpiecEasi 패키지 내 함수 추정) + igraph, `ginseng_spieceasi` | bootstrap=100, SPIEC-EASI 대비 Jaccard 비교 | Step14와 동일 입력 | `sparcc_pooled_cor.rds`, `sparcc_sig_edges.tsv`, `network_comparison.tsv` |
| 16. MOFA+ 통합(최종 통합 스텝) | MOFA2(R), `ginseng_r` | Model A: CLR(888 ASV)+eigengene(14); Model B: CLR(888)+top3000 genes; 샘플ID 매핑 `{연생}y-{반복}`↔`Pg_{연생}Y_{반복}` | `otu_clr_selected.rds`(16S 축) + `logcpm.rds`, `module_eigengenes.csv`(mRNA-seq 축) | `modelA_model.hdf5`, `modelB_model.hdf5`, `*_factor_year_correlation.tsv`, `*_top_weights_*.tsv` |

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
