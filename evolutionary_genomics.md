# 진화유전체 분석 (Zoysia_sinica_Genome_Assembly)

## 1. 공통 파이프라인

이 카테고리는 `de_novo_genome_assembly.md`에 정리된 1차 조립·주석 이후 단계로, comparative_genomics.md
(OrthoFinder 오솔로그 클러스터링)의 결과물을 입력으로 사용하는 두 개의 후속 분석 스레드를 다룬다.
스레드 A(유전자가족 진화)는 "단일카피 오솔로그 정렬 → 계통수 추정 → ultrametric 변환 → CAFE5 유전자
가족 확장/축소 분석"과 "동일 정렬/계통수를 재사용한 PAML MCMCTree 분기시간 추정"으로 구성되며, 두
분석 모두 같은 11종(Arabidopsis_thaliana, Brachypodium_distachyon, Oryza_sativa, Oropetium_thomaeum,
Sorghum_bicolor, Setaria_italica, Zoysia_japonica, Zea_mays, Zoysia_matrella, Zoysia_pacifica,
Zoysia_sinica) 단일카피 오솔로그 concatenated 정렬을 공유한다. 스레드 B(PSG, 계통별 양성선택 유전자
탐색)는 OrthoFinder 오솔로그그룹(comparative_genomics.md 참조)에서 유전자별 CDS/단백질 서열을 다시
추출해 별도의 mafft(translatorX)→pal2nal→Gblocks→PAML codeml branch-site 파이프라인을 돌리는 독립
스레드로, 스레드 A와 입력 오솔로그 소스만 공유할 뿐 정렬·트리 산출물을 재사용하지 않는다.

최종 채택된 계통수는 `7_CAFE/Manual_Tree/` 경로다 — 11종 단일카피 오솔로그를 mafft+Gblocks로 직접
재정렬해 RAxML로 새로 추정한 트리이며, GO term 주석·`Zs_Expansion/Contraction` 계통별 추출 파일이
이 경로에만 존재하고 `8_MCMC_Tree/Manual_Tree/`도 동일한 RAxML 트리를 재사용한다(OrthoFinder 자체
종트리를 그대로 쓰는 `OrthoFinder_Tree/` 경로도 남아있으나 비교/검증용 병행 시도로 판단됨 — 확인
필요). Spartina alterniflora를 추가 외군으로 포함한 12종 확장판 CAFE5 실행(`7_CAFE/Spartina/`)도
존재하며, 이는 폐기된 시도가 아니라 별도의 확장 분석이다.

```
[comparative_genomics.md: OrthoFinder 오솔로그그룹 결과]
        │
        ├──▶ (스레드 A: 유전자가족 진화 + 분기시간)
        │      │
        │      ▼
        │  [Step 1. 11종 단일카피 오솔로그 mafft 정렬 + Gblocks 트리밍 (Manual_Tree 경로)]
        │      │         (병행: OrthoFinder_Tree 경로 — OrthoFinder 자체 종트리를 그대로 사용, 상세 불명)
        │      ▼
        │  [Step 2. RAxML 계통수 추정 + ultrametric 변환]
        │      │
        │      ├──▶ [Step 3. CAFE5 유전자가족 확장/축소 분석 (Manual_Tree 채택, OrthoFinder_Tree 병행,
        │      │             Spartina 외군 추가판(P0.05/P0.01)은 후속 확장 실행)]
        │      │
        │      └──▶ [Step 4. PAML MCMCTree 분기시간 추정 (Hessian 계산 → usedata=2 MCMC 1st/2nd round)]
        │
        └──▶ (스레드 B: 계통별 양성선택 유전자 탐색, 독립 진행)
               │
               ▼
           [Step 5. 오솔로그그룹별 CDS 서열 준비 및 ID 정리 (Rename.py, Preparation_CDS_sequence.sh)]
               │
               ▼
           [Step 6. mafft(translatorX) 아미노산 정렬 → pal2nal 코돈 정렬 변환]
               │
               ▼
           [Step 7. Gblocks 코돈 트리밍(-t=c) → PHYLIP 변환/병합(all_combined.phy)]
               │
               ▼
           [Step 8. PAML codeml branch-site 모델(전경) vs null 모델 비교 + LRT 검정
                     (Zs/ZsSa/Zj/Zmt/Zp 5개 전경 계통 반복)]
```

### Step 1. 11종 단일카피 오솔로그 mafft 정렬 + Gblocks 트리밍
- **소프트웨어(버전):** mafft(버전 확인 필요 — 스크립트 파일 없음), Gblocks 0.91b(`.htm` 리포트 헤더에서 확인)
- **핵심 파라미터:** Gblocks — Minimum Number Of Sequences For A Conserved Position: 6 / For A Flanking
  Position: 9 / Maximum Number Of Contiguous Nonconserved Positions: 8 / Minimum Length Of A Block: 10 /
  Allowed Gap Positions: None / Use Similarity Matrices: Yes (11종 기준 기본값과 일치 — `-t=p` 단백질
  모드로 추정되나 명령행 자체는 로그에 남아있지 않음)
- **Input:** `zoysia_orthofinder_gene_family.txt`/`filtered_cafe_input.txt`에서 확인되는 11종 단일카피
  오솔로그(OrthoFinder 결과, comparative_genomics.md 참조)의 concatenated 단백질 서열
  (`Aligned_SingleCopy_11species.fasta`, 헤더가 `SingleCopy_AT`/`SingleCopy_Bdi`/... 형태로 종별 1개씩
  11개 서열)
- **Output:** `Aligned_SingleCopy_11species.fasta-gb`(Gblocks 트리밍 후 정렬, 1943개 block 선택),
  `Aligned_SingleCopy_11species.fasta-gb.htm`(Gblocks 리포트)
- **실행 커맨드:**
  ```bash
  # mafft 정렬 커맨드 자체는 스크립트/로그로 남아있지 않음 — 확인 필요, 추정 금지
  Gblocks Aligned_SingleCopy_11species.fasta-gb  # 실제 인자는 .htm 파라미터로만 역추적 가능, 확인 필요
  ```

### Step 2. RAxML 계통수 추정 + ultrametric 변환
- **소프트웨어(버전):** RAxML 8.2.12 (`RAxML_info.*` 헤더에서 확인)
- **핵심 파라미터:** `raxmlHPC-PTHREADS-AVX -s Aligned_SingleCopy_11species.fasta-gb -n Aligned_SingleCopy_11species.fasta-gb.tree.txt -m PROTGAMMALGX -p 12345 --bootstop-perms=1000 -T 200 -o SingleCopy_AT`
  (`RAxML_info.Aligned_SingleCopy_11species.fasta-gb.tree.txt`의 "RAxML was called as follows" 라인에서
  그대로 확인; LG+GAMMA 모델, ML 추정 아미노산 빈도, Arabidopsis_thaliana(`SingleCopy_AT`)를 outgroup으로
  고정)
- **Input:** `Aligned_SingleCopy_11species.fasta-gb`
- **Output:** `RAxML_bestTree.Aligned_SingleCopy_11species.fasta-gb.tree.txt`(최종 ML 트리, Likelihood
  -730655.898427), `RAxML_bestTree.*.ultrametric.tre`(ultrametric 변환본 — CAFE5 입력용). ultrametric
  트리의 root 분기시각이 160(백만 년)으로 고정되어 있는데, 이는 아래 Step 4 MCMCTree의 root 화석 보정
  구간(151.56–170.06 Mya)과 정확히 겹친다 — 두 분석이 같은 root 나이 가정을 공유하도록 스케일링된
  것으로 보이나, ultrametric 변환에 사용된 스크립트(예: CAFE 튜토리얼의 `make_ultrametric.py` 류) 자체는
  파일로 남아있지 않다 — 확인 필요.
- **실행 커맨드:**
  ```bash
  raxmlHPC-PTHREADS-AVX -s Aligned_SingleCopy_11species.fasta-gb -n Aligned_SingleCopy_11species.fasta-gb.tree.txt -m PROTGAMMALGX -p 12345 --bootstop-perms=1000 -T 200 -o SingleCopy_AT
  # RAxML_bestTree → ultrametric 변환 스크립트/커맨드 자체는 확인 필요 — 스크립트 파일 없음
  ```

### Step 3. CAFE5 유전자가족 확장/축소 분석
- **소프트웨어(버전):** CAFE5(버전 확인 필요 — 실행 커맨드 자체는 `Spartina/` 로그에서만 확인되고
  최종 채택된 Manual_Tree 경로에는 cafe5 로그가 남아있지 않음)
- **핵심 파라미터:** Base 모델(오솔로그그룹별 birth-death 모델), Manual_Tree 경로 Lambda 0.00625
  (최대 가능값과 동일, `Base_results.txt`). `Spartina/` 확장판은
  `cafe5 -i filtered_cafe_input.txt -t SpeciesTree_rooted.txt`(p-value 기본값 0.05) 및
  `--pvalue 0.01` 두 버전을 실행(`cafe5.log`/`cafe5_P0.01.log`에서 커맨드라인 그대로 확인)
- **Input:** `filtered_cafe_input.txt`(오솔로그그룹별 종 카운트 테이블, comparative_genomics.md
  OrthoFinder 결과 기반) + Step 2의 `RAxML_bestTree.*.ultrametric.tre`. Spartina 확장판은
  `Orthogroups.GeneCount.tsv` 기반 재계산된 `filtered_cafe_input.txt`(12종, Spartina alterniflora
  추가) + `SpeciesTree_rooted.txt`
- **Output:** `Base_report.cafe`, `Base_results.txt`, `Base_family_results.txt`(오솔로그그룹별 p-value),
  `Base_clade_results.txt`, `Base_asr.tre`, `Zs_Expansion_ID.txt`/`Zs_Contraction_ID.txt`/
  `Zs_specific_Expansion_ID.txt`와 이를 GO term에 매핑한 `Zs_Expansion_Significant.ZsGOTerm` 등
  Zoysia_sinica 계통 특이 확장/축소 유전자 목록
- **실행 커맨드:**
  ```bash
  # Spartina 확장판(로그에서 커맨드라인 직접 확인):
  cafe5 -i filtered_cafe_input.txt -t SpeciesTree_rooted.txt
  cafe5 -i filtered_cafe_input.txt -t SpeciesTree_rooted.txt --pvalue 0.01
  # Manual_Tree 채택 경로의 실제 cafe5 실행 커맨드 자체는 로그 파일 없음 — 확인 필요
  # (입력 파일 이름과 Base_report.cafe 트리 구조로 미루어 ultrametric 트리 + filtered_cafe_input.txt를
  #  사용했음은 확인되나, 정확한 커맨드라인 인자는 재구성 불가)
  ```

### Step 4. PAML MCMCTree 분기시간 추정
- **소프트웨어(버전):** PAML mcmctree — Hessian 계산 단계 4.10.7(June 2023, `matrix.log` 헤더),
  실제 MCMC 1st/2nd round 단계 4.10.6(November 2022, `Manual_Tree/1st_round/nohup.out` 헤더) — 같은
  분석 파이프라인 내에서 두 세부 버전이 혼용됨
- **핵심 파라미터:** 2단계 표준 PAML 근사 절차 — ① `usedata=3`(Hessian/gradient out.BV 계산, `clock=3`
  correlated rate, `RootAge='<1.0'`, `model=2` WAG+감마, `alpha_gamma=1 1`), ② 산출된 `in.BV`를
  `usedata=2`(normal approximation)로 두 번(1st_round/2nd_round, 서로 다른 랜덤시드 — SeedUsed
  1184356949 vs 543083339) 독립 MCMC 실행하여 수렴 확인. 화석/보정 구간: 벼-포아과 분화
  `B(41.51749, 62, 0.025, 0.025)`, Oropetium-Zoysia류 분화 `B(32.01, 39.48251, 0.025, 0.025)`,
  Zoysia속 내 분화 `B(1.08364, 4.25047, 0.025, 0.025)`, root(Arabidopsis 분기) `B(151.56, 170.05738,
  0.025, 0.025)`. `burnin=5000000, sampfreq=30, nsample=10000000`
- **Input:** Step 1의 Gblocks 트리밍 정렬을 phylip으로 변환한 `Aligned_SingleCopy_11species.fasta-gb.phylip`
  (`Manual_Tree/` 하위), Step 2의 RAxML 트리에 화석 보정 구간을 추가한 `mcmtree.nwk`
- **Output:** `Manual_Tree/2nd_round/out.txt`(최종 채택 posterior) 기준 root(Arabidopsis 분기) 평균
  분기시각 160.53 Mya(95% HPD 151.48–169.97), Oropetium-Zoysia류 분기 33.08 Mya, Zoysia속 내
  Zoysia_pacifica-Zoysia_sinica 분기 2.48 Mya(95% HPD 1.27–3.49), Zoysia_matrella-Zoysia_japonica 분기
  2.48 Mya; `FigTree.tre`(NEXUS, 95% HPD 포함 시각화용 트리), `mcmc.txt`(MCMC 샘플 원본)
- **실행 커맨드:**
  ```bash
  mcmctree mcmctree_hessian_matrix.ctl   # usedata=3, Hessian/gradient(out.BV → in.BV) 계산
  # 이하 1st_round/, 2nd_round/ 각각에서 in.BV 재사용:
  mcmctree mcmctree_hessian_divergence.ctl   # usedata=2, clock=3, burnin=5000000 sampfreq=30 nsample=10000000
  ```

### Step 5. 오솔로그그룹별 CDS 서열 준비 및 ID 정리
- **소프트웨어(버전):** bash/grep/awk 조합 스크립트(`Preparation_CDS_sequence.sh`), Python 3(`Rename.py`)
- **핵심 파라미터:** 없음(단순 텍스트 처리) — `Rename.py`는 헤더에서 `.p` 접미사를 제거하고, `Zm0`로
  시작하는 헤더의 `P`를 `T`로 치환하는 ID 정규화만 수행
- **Input:** `01_protein_sequences/OG*.fa`(오솔로그그룹별 단백질 서열, 최종 채택본은 7종
  구성 — Os/Oro/GWH/Zjn/Zmt/Zp/Zsg. `GWH` 접두사는 Spartina alterniflora 유전자 ID임을
  `01_protein_sequences/*.fa` 헤더(`GWHPCBIM...`)로 확인), `02_dna_sequences/`의 종별 전체 CDS 서열
- **Output:** `01_protein_sequences/OG*_cds.fasta`(단백질 헤더 ID로 매칭된 CDS 서열)
- **실행 커맨드:**
  ```bash
  python Rename.py
  bash Preparation_CDS_sequence.sh
  ```
  (`Preparation_CDS_sequence.log`는 파일만 존재하고 내용은 비어 있음 — 실행 자체는 `nohup.out`
  존재로 정황상 확인되나 로그 상세 확인 불가)
- 최종 채택 구성은 7종(top-level, 2024-07-09~08-02) — Os/Oro/GWH(Spartina)/Zjn/Zmt/Zp/Zsg.

### Step 6. mafft(translatorX) 아미노산 정렬 → pal2nal 코돈 정렬 변환
- **소프트웨어(버전):** translatorX(내부적으로 mafft 호출, 파일명 `*_mafft_translatorx.*`로 확인,
  버전 확인 필요 — translatorX/mafft 실행 스크립트 자체는 파일로 남아있지 않음), pal2nal.pl(버전
  확인 필요)
- **핵심 파라미터:** pal2nal `-output fasta`(`pal2nal.sh`에서 확인) 외 확인 필요
- **Input:** `01_protein_sequences/OG*_cds.fasta`(단백질) + 대응 CDS 서열
- **Output:** `03_alignments_mafft/OG*_mafft_translatorx.aa_ali.fasta`(아미노산 정렬),
  `OG*_mafft_translatorx.nt_ali.fasta`(코돈 정렬), `pal2nal_checks/OG*_pal2nal_mafft_out.fasta`(pal2nal
  변환 결과)
- **실행 커맨드:**
  ```bash
  # translatorX/mafft 정렬 자체의 커맨드라인은 로그 없음 — 확인 필요, 추정 금지
  # pal2nal.sh (실제 스크립트 그대로):
  for aa_file in *_mafft_translatorx.aa_ali.fasta; do
      og_number=$(echo "$aa_file" | grep -o 'OG[0-9]\+')
      nt_file="${og_number}_mafft_translatorx.nt_ali.fasta"
      pal2nal.pl $aa_file $nt_file -output fasta > pal2nal_checks/${og_number}_pal2nal_mafft_out.fasta &
  done
  ```

### Step 7. Gblocks 코돈 트리밍 → PHYLIP 변환/병합
- **소프트웨어(버전):** Gblocks(버전 확인 필요, 코돈 모드), `one_line_fasta.pl` /
  `FASTAtoPHYL.pl`(paml-tutorial 부속 스크립트, `/usr/local/bin/paml-tutorial/positive-selection/00_data/scripts/`
  경로에 설치된 것으로 확인, 버전 없음)
- **핵심 파라미터:** `Gblocks <file> -t=c`(codon 모드, `Gblock.sh`에서 그대로 확인)
- **Input:** `pal2nal_checks/OG*_pal2nal_mafft_out.renamed.fasta`(pal2nal 결과를 ID 정규화한 버전)
- **Output:** `OG*_pal2nal_mafft_out.renamed.fasta-gb`(트리밍된 코돈 정렬, 예: OG0016779는 1707→1503bp,
  88% 유지; OG0016781은 9075→942bp, 10%만 유지 — 오솔로그그룹별 편차가 큼, `Gblock.log`에서 확인),
  개별 PHYLIP 파일을 합친 `all_combined.phy`(265개 오솔로그그룹 병합, PAML codeml 입력)
- **실행 커맨드:**
  ```bash
  # Gblock.sh (실제 스크립트 그대로):
  for x in OG*renamed.fasta; do
      Gblocks ${x} -t=c
  done
  # Final_alignment.sh (실제 스크립트 그대로):
  for y in *gb; do
      /usr/local/bin/paml-tutorial/positive-selection/00_data/scripts/one_line_fasta.pl ${y}
  done
  for z in *one_line.fa; do
      sed 's/[ \t]//g' -i ${z}
      num=$( grep '>' ${z} | wc -l )
      len=$( sed -n '2,2p' ${z} | sed 's/\r//' | sed 's/\n//' | wc -L )
      if [ "$len" -ne 0 ]; then
          perl /usr/local/bin/paml-tutorial/positive-selection/00_data/scripts/FASTAtoPHYL.pl ${z} $num $len
      fi
  done
  for file in *phy; do cat "$file"; echo; done > all_combined.phy
  ```

### Step 8. PAML codeml branch-site 모델 vs null 모델 비교 + LRT 검정
- **소프트웨어(버전):** PAML codeml 4.10.7(`codeml-branchsite.log` 헤더에서 확인), R(dplyr, `LRT_analysis.R`)
- **핵심 파라미터:** branch-site 대안 모델 — `model=2, NSsites=2, CodonFreq=7, estFreq=0, clock=0,
  fix_omega=0, omega=2`(전경 계통 ω 자유 추정); null 모델 — 동일 설정에 `fix_omega=1, omega=1`(전경
  ω=1로 고정, `null_model/codeml-branchsite.ctl`). `ndata=265`(265개 오솔로그그룹 배치 처리),
  `seqtype=1`(코돈). 전경 계통은 트리 파일에 `#1`로 표시 — Zs 계통(`(Zp,Zsg #1)`), ZsSa 계통(Zsg와
  Spartina/GWH를 함께 `#1`로 표시, `(Zp,Zsg #1) ... GWH #1`), 그리고 Zj/Zmt/Zp 계통(2024-08-02 추가
  실행)까지 총 5개 전경 계통 조합을 각각 별도 폴더(`05_branchsite_model_Zs/ZsSa/Zj/Zmt/Zp`)에서 반복
- **Input:** `all_combined.phy`(Step 7 산출물), `SpeciesTree_rooted.txt`(전경 계통 `#1` 표시가 된
  7종 트리, 예: `(Os,(Oro,(((Zmt,Zjn),(Zp,Zsg #1)),GWH)));`)
- **Output:** `branchsite_Zs_MLEs.txt`(사이트클래스별 ω 추정치), `lnL_branchsite_mods.txt` /
  `null_model/lnL_branchsite_mods.txt`(대안/null 모델 log-likelihood), `LRT_analysis.R`로 계산한
  `significant_genes.txt`(LRT 통계량 2×Δ(lnL), df=1 카이제곱 검정, p<0.05), `MLEs_ver2.py`로 추가
  필터링한 `branchsite_Zs_MLEs_filtered.txt`(배경 ω<1이면서 전경 ω≥1인 사이트클래스 존재 조건),
  최종 `Positive_selected_OGs.txt`(양성선택 후보 오솔로그그룹 ID 목록, 예: Zs 계통에서 OG0016779,
  OG0016781, OG0016783 등)
- **실행 커맨드:**
  ```bash
  # codeml-branchsite.ctl (대안 모델, 실제 파일 그대로 — Zs 계통 예시):
  seqfile = all_combined.phy
  treefile = SpeciesTree_rooted.txt
  outfile = out_Zs_unrooted_branchsite.txt
  seqtype = 1
  ndata = 265
  icode = 0
  cleandata = 0
  model = 2
  NSsites = 2
  CodonFreq = 7
  estFreq = 0
  clock = 0
  fix_omega = 0
  omega = 2
  codeml codeml-branchsite.ctl

  # null_model/codeml-branchsite.ctl (fix_omega=1, omega=1 외 동일):
  codeml codeml-branchsite.ctl

  # LRT_analysis.R (실제 스크립트 핵심부):
  lnL_data$LRT_statistic <- 2 * (lnL_data$branchsite_lnL - lnL_data$nullmodel_lnL)
  lnL_data$p_value <- pchisq(lnL_data$LRT_statistic, df = 1, lower.tail = FALSE)
  significant_genes <- lnL_data[lnL_data$p_value < 0.05, ]
  ```

## 2. 프로젝트별 특이사항

### Zoysia_sinica_Genome_Assembly
- (단일 프로젝트 카테고리 — 위 공통 파이프라인이 곧 이 프로젝트의 파이프라인)
- 연구목적/샘플 설계: 들잔디(Zoysia japonica)속 근연 4종(Z. sinica/pacifica/matrella/japonica)을 포함한
  11종(및 확장판에서 Spartina alterniflora 포함 12종) 비교유전체 데이터셋을 이용해, Zoysia_sinica
  계통에서 일어난 유전자가족 확장/축소, 분화시기, 그리고 계통 특이 양성선택 유전자를 규명하려는 분석.
- 프로젝트 내부의 분기:
  - CAFE5 최종 해석(GO term 매핑, 계통 특이 확장/축소 유전자 추출)은 `Manual_Tree`(RAxML로 재구축한
    계통수) 경로 산출물만을 사용했으며, Spartina alterniflora를 추가 외군으로 포함한 12종 확장판
    (p-value 0.05/0.01 두 버전)이 별도 후속 확장으로 존재한다.
  - PSG(양성선택) 분석은 최종적으로 Spartina alterniflora를 포함한 7종 구성으로 수행되었고, 5개
    전경 계통(Zs/ZsSa/Zj/Zmt/Zp) 각각에 대해 독립적으로 branch-site 검정을 실행했다.
  - MCMCTree Hessian 계산 단계와 실제 MCMC 표본추출 단계에서 PAML 세부 버전(4.10.7 vs 4.10.6)이 서로
    다르게 기록되어 있다 — 재현성 측면에서 확인이 필요한 부분.

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Zoysia_sinica_Genome_Assembly | CAFE5는 mafft+Gblocks+RAxML로 재구축한 트리(Manual_Tree)를 채택; PSG는 Spartina 포함 7종 구성으로 최종 5개 전경 계통(Zs/ZsSa/Zj/Zmt/Zp)에 대해 반복 실행 | RAxML 8.2.12, Gblocks 0.91b, PAML(mcmctree 4.10.7/4.10.6 혼용, codeml 4.10.7) — CAFE5/mafft/translatorX/pal2nal 버전은 확인 필요 | Spartina alterniflora를 외군으로 추가한 12종 CAFE5 확장판(P0.05/P0.01); ZsSa(Zoysia_sinica+Spartina 동시 전경) branch-site 모델 |
