# 비교유전체 분석 (Zoysia_sinica_Genome_Assembly)

## 1. 공통 파이프라인

이 카테고리는 `de_novo_genome_assembly.md`에서 다루는 Zoysia_sinica 핵심 조립/주석 파이프라인의
후속 단계 중, "근연종 유전자 재주석 → 오솔로그 클러스터링 → 공선성(Ka/Ks) 분석 → Circos 시각화"로
이어지는 하나의 연속된 분석 흐름을 다룬다. 이 네 개 폴더(`6_OrthoFinder`, `9_Ks_distribution`,
`12_Circos`, `Zoysia_Gene_reannotation`)는 서로 독립적이지 않고, `Zoysia_Gene_reannotation`에서 만든
근연종(Z. japonica/matrella/pacifica) + Z. sinica 자체 재주석 결과가 `6_OrthoFinder`와
`9_Ks_distribution`의 입력(단백질/GFF)으로 그대로 재사용되며, `9_Ks_distribution`의 공선성(MCScanX)
결과 일부가 다시 `12_Circos`의 링크 트랙 입력으로 이어진다. 단, `6_OrthoFinder`의 산출물 자체가
`12_Circos`에 직접 연결되지는 않는다 — `12_Circos`는 Z. sinica 게놈 두 버전(Zsg/Zsg2로 표기된 것으로
추정, 확인 필요) 간의 별도 MCScanX self-synteny 결과를 사용한다.

`6_OrthoFinder`의 오솔로그 그룹 결과는 이 카테고리 밖의 `7_CAFE`(별도 태스크인
`evolutionary_genomics.md`에서 다룸)의 입력(`zoysia_orthofinder_gene_family.txt`)으로도 재사용되는데,
그 파일의 종 구성(11종: Arabidopsis_thaliana/Brachypodium_distachyon/Oryza_sativa/
Oropetium_thomaeum/Sorghum_bicolor/Setaria_italica/Zoysia_japonica/Zea_mays/Zoysia_matrella/
Zoysia_pacifica/Zoysia_sinica)과 생성 시각(2024-06-20)이 `OrthoFinder.log`(2024-06-19 실행,
`Results_Jun19`, 11개 diamond DB 생성)와 정확히 일치하므로, CAFE로 이어지는 것은 이 최초 11종 실행
결과로 판단된다(아래 Step 2 참조).

```
[Zoysia_Gene_reannotation]
  BRAKER2/3 (Zj/Zmt/Zp/Zs 등 근연종 유전자 예측)
        │
        ▼
  AGAT 후처리 (longest isoform 선별 → gene ID rename)
        │  (.longest.rename.gtf/.pep → 아래 두 갈래 입력으로 재사용)
        │
        ├─────────────────────────────────────────────┐
        ▼                                               ▼
[Step 1. Zoysia_Gene_reannotation]              [Step 3. GFF→MCScan bed 변환]
  (위 두 스텝, 아래에 별도 기술)                          │  (9_Ks_distribution)
        │                                               ▼
        ▼                                    [Step 4. MCScanX 공선성 + Ka/Ks 계산]
[Step 2. OrthoFinder 오솔로그 클러스터링]                  │  (11종 전체 blastp all-vs-all,
  (6_OrthoFinder, 11종 Results_Jun19 채택)                 641,046초 소요)
        │                                               ▼
        │  (CAFE 입력으로 전달 — evolutionary_genomics.md 참조)   [Step 5. Ks 분포 / WGD 시각화]
        │                                               │
        └──────────────────(직접 연결 없음)──────────────┘
                                                          │
                                                          ▼
                                          [Step 6. Circos 시각화]
                                            (12_Circos, karyotype+LTR+
                                             genomic-block 링크, 별도 self-synteny 입력)
```

### Step 1. 근연종 유전자 재주석 (Zoysia_Gene_reannotation)
- **소프트웨어(버전):** BRAKER (버전 문자열은 로그에 직접 나타나지 않음 — 확인 필요; `what-to-cite.txt`
  인용 목록으로 실행 모드만 구분 가능), AGAT v1.2.0(각 `*.agat.log` 헤더에서 직접 확인 — 참고: 같은
  프로젝트의 최종 채택 유전자예측(`de_novo_genome_assembly.md` Step 9)은 AGAT v1.4.0을 사용해 버전이
  다르다. 이 폴더의 재주석은 그와 별개의 실행/배치로 보인다)
- **핵심 파라미터:** 종별로 증거 데이터 조합이 다르다 — Zp는 단백질 증거만(BRAKER2 모드,
  `what-to-cite.txt`에 BRAKER3 인용 없음), Zmt/Zj(Whole_SRA)는 단백질+공개 RNA-seq(SRA) 증거
  (BRAKER3 모드), Zs는 단백질+IsoSeq 정렬 bam 증거(BRAKER3 모드)
- **Input:** 각 종의 반복서열 마스킹 게놈(`ZPZ_r1.0.genome.fasta.masked`,
  `ZMW_r1.0.genome.fasta.masked`, `ZJN_r1.1.genome.fasta.masked`,
  `Z.sinica_IO_genome.original.fasta`), 공통 단백질 DB `liliopsida_27sp.ncbi_pep.faa`, (Zmt/Zj)
  SRA RNA-seq accession 목록, (Zs) `isoseq.bam`
- **Output:** 종별 `braker.gtf`(→ `Z*_br*.gtf`로 리네이밍되어 저장), longest-isoform GTF/GFF3/pep
  (`*_br*.longest.gtf`, `.longest.pep`), 최종 gene ID rename 버전(`*.longest.rename.gtf/.pep/.cds`)
  — 이 rename 버전이 `6_OrthoFinder`의 references와 `9_Ks_distribution`에 심볼릭 링크로 재사용됨
- **실행 커맨드:**
  ```bash
  # Zp (Braker2_Zp/braker2_zp.sh) — 단백질 증거만
  braker.pl --genome=ZPZ_r1.0.genome.fasta.masked --prot_seq=liliopsida_27sp.ncbi_pep.faa --threads 20

  # Zmt (Braker_Zmt_Whole_SRA/braker3_Zmt.sh) — 단백질 + 공개 RNA-seq(SRA) 22개
  braker.pl --genome=ZMW_r1.0.genome.fasta.masked --prot_seq=liliopsida_27sp.ncbi_pep.faa \
    --rnaseq_sets_ids=SRR12732216,SRR12732217,SRR12732218,SRR13132433,SRR13132434,SRR13132435,SRR13132436,SRR13132438,SRR13132439,SRR13132440,SRR13132441,SRR13132442,SRR13132443,SRR13132444,SRR13132445,SRR13132437,SRR19710349,SRR19710348,SRR19710343,SRR19710342,SRR19710341,SRR19710340 \
    --threads 20

  # Zs (Braker_Zs/Braker_Zs.sh) — 단백질 + IsoSeq bam (스크립트 원문 그대로, "–-bam" 대시 문자가
  # 일반 하이픈이 아닌 오타로 보이나 원본 그대로 기재)
  #minimap2 -t 128 -ax splice:hq -uf Z.sinica_IO_genome.original.fasta Zs_Pac_hq_transcripts_total.renamed.fastq > isoseq.sam
  #samtools view -bS --threads 128 isoseq.sam -o isoseq.bam
  braker.pl --genome=Z.sinica_IO_genome.original.fasta --prot_seq=../liliopsida_27sp.ncbi_pep.faa –-bam=isoseq.bam --threads=20
  ```
  Zj(japonica)는 `Braker_Zj_original_set/`, `Braker_Zj_Whole_SRA/` 두 폴더 모두 실행 스크립트(`.sh`)가
  남아있지 않다 — 확인 필요 — 스크립트 파일 없음 (단, `braker.log`/`what-to-cite.txt`와 결과 gtf는 존재).
  두 폴더명(original_set vs Whole_SRA)으로 보아 Zmt와 마찬가지로 "원본 단백질 증거만" 버전과
  "전체 공개 RNA-seq 추가" 버전을 각각 실행한 것으로 추정된다.

  AGAT 후처리 3단계(원본 gtf → longest isoform → rename)는 각 단계마다 `*.agat.log`가 남아있으나,
  로그 자체에는 AGAT가 자신을 호출한 커맨드라인이 기록되지 않는다(버전/체크리스트만 출력) — 정확한
  `agat_*.pl` 서브커맨드와 옵션은 확인 필요. 파일명 규칙(`.longest` = longest isoform 선별,
  `.rename` = gene/transcript ID 재명명)으로 미루어 각각 `agat_sp_keep_longest_isoform.pl` 계열과
  ID rename 스크립트(AGAT 또는 커스텀)가 사용된 것으로 추정되나 확정할 수 없다.

### Step 2. OrthoFinder 오솔로그 클러스터링 (6_OrthoFinder)
- **소프트웨어(버전):** OrthoFinder v2.5.5 (모든 `OrthoFinder*.log` 공통 헤더에서 확인)
- **핵심 파라미터:** `-a 200`(알고리즘 스레드 200), `-M msa`(MSA 기반 트리 추론), `-X`(시퀀스 ID
  변환 비활성화); diamond 256 스레드 all-vs-all
- **Input:** `OrthoFinder_input/`의 종별 단백질 FASTA(대부분 `references/`로의 심볼릭 링크 —
  Zjaponica/Zmatrella/Zpacifica/Zsinica는 `Zoysia_Gene_reannotation`의 `*.longest.rename` 계열 단백질을
  가리킴), 총 11개 diamond DB 생성(최초 실행 기준)
- **Output:** `OrthoFinder_input/OrthoFinder/Results_*/`(Orthogroups, Orthologues, 종 트리 등);
  `SingleCopy/`(단일카피 오솔로그 concat FASTA, 종 트리용)
- **최종 채택 실행:** `OrthoFinder.log`(2024-06-19, `Results_Jun19`) — 11종(AT/Bdi/Os/Oro/Sob/Sei/
  Zjn/Zm0/Zmw/Zp/Zsg). `SingleCopy/`(2024-06-20)와 `7_CAFE/zoysia_orthofinder_gene_family.txt`
  (2024-06-20, 동일 11종 헤더)가 이 실행 직후 시각에 생성되어 CAFE로 전달되는 실행으로 판단됨
  (직접 연결을 증명하는 스크립트는 없고 시각/종 구성 일치에 근거한 추정 — 확인 필요). 이후
  Spartina alterniflora를 외군으로 추가한 12종 실행(`OrthoFinder_Spartina.log`, `Results_Jul06`)이
  별도로 존재하며, 이 결과는 `9_Ks_distribution/Spartina/`의 12종 Ks 확장 분석 입력으로 이어진다.
- **실행 커맨드:**
  ```bash
  # OrthoFinder.sh
  orthofinder.py -f OrthoFinder_input/ -a 200 -M msa -X
  ```

### Step 3. GFF→MCScan bed 포맷 변환 (9_Ks_distribution)
- **소프트웨어(버전):** awk 기반 커스텀 변환 스크립트(3종, 별도 버전 없음)
- **핵심 파라미터:** 입력/출력 파일명은 실행 시 대화형으로 입력받음(`read -p`) — 하드코딩된 인자 없음
- **Input:** 종별 GFF3 (예: `Zjaponica.gff3`, `Osativa.gff3` 등)
- **Output:** 4컬럼 bed형 파일(`chr, gene_id, start, end`) — MCScanX 입력용 gff 포맷으로 사용
- **실행 커맨드:**
  ```bash
  # Converting_for_MCScan_gff3.sh — mRNA 라인에서 9번째 컬럼의 Name= 값을 추출(일반 GFF3용)
  awk -F'\t' '
  BEGIN { OFS="\t" }
  $3 == "mRNA" {
      match($9, /Name=([^;]+)/, arr)
      name = arr[1]
      print $1, name, $4, $5
  }
  ' "$input_file" > "$output_file"

  # Converting_for_MCScan_gff3_2.sh — transcript 라인의 9번째 컬럼 전체를 그대로 사용
  awk -F'\t' '
  BEGIN { OFS="\t" }
  $3 == "transcript" { print $1, $9, $4, $5 }
  ' "$input_file" > "$output_file"

  # Converting_for_MCScan_gff3_Os.sh — mRNA 라인에서 첫 ';' 이전 필드의 ID= 접두어 제거(Oryza sativa GFF 전용)
  awk -F'\t' '
  BEGIN { OFS="\t" }
  $3 == "mRNA" {
      split($9, fields, ";")
      sub(/^ID=/, "", fields[1])
      id = fields[1]
      print $1, id, $4, $5
  }
  ' "$input_file" > "$output_file"
  ```

### Step 4. MCScanX 공선성 분석 및 Ka/Ks 계산 (9_Ks_distribution)
- **소프트웨어(버전):** MCScanX(버전 문자열 로그에 없음 — 확인 필요), 내부적으로 ClustalW 2.1을
  호출하여 페어와이즈 정렬 후 Ka/Ks 계산(`MCScanX_11sp_kaks.log`에서 "CLUSTAL 2.1" 확인)
- **핵심 파라미터:** 11종 전체 blastp all-vs-all 결과(`11sp_for_MCScanX.blastp`, 약 3.5GB)를
  두 차례 보정한 `.fixed.blast` → `.fixed2.blast`를 최종 입력으로 사용(중간 보정 내용은 확인 필요 —
  아마 gene ID 불일치 수정)
- **Input:** `11sp_for_MCScanX.fixed2.blast`, `11sp_for_MCScanX.fixed2.gff`(Step 3 산출물 병합)
- **Output:** `11sp_for_MCScanX.fixed2.collinearity`(공선성 블록), `.tandem`(탠덤 중복 쌍),
  `.kaks`(Ka/Ks 값, 108MB), `.html/`(블록별 시각화); Zoysia 전용 실행은 `Zoysia.collinearity` /
  `Zoysia.tandem` 및 duplicate-type 분류 결과 생성
- **실행 커맨드:** 실행 스크립트(`.sh`)는 없고 표준 출력을 리다이렉트한 로그만 남아있어 정확한
  MCScanX 바이너리 호출 커맨드라인 자체는 확인 필요. 로그로부터 재구성한 처리 내용은 다음과 같다.
  ```text
  # MCScanX.log (11종 전체, 641,046초 ≈ 7.4일 소요)
  Reading BLAST file and pre-processing
  28106501 matches imported (17891178 discarded)
  2952311 pairwise comparisons
  73837 alignments generated
  Pairwise collinear blocks written to 11sp_for_MCScanX.fixed2.collinearity
  Tandem pairs written to 11sp_for_MCScanX.fixed2.tandem

  # MCScanX_Zoysia.log (Zoysia 4종 한정, 25,055초 ≈ 7시간)
  4773956 matches imported (2780974 discarded)
  1838465 pairwise comparisons
  21512 alignments generated
  Pairwise collinear blocks written to Zoysia.collinearity
  Tandem pairs written to Zoysia.tandem

  # MCScanX_duplicate_clafssifier_Zoysia.log (중복 유형 분류, 원문 파일명의 오탈자 "clafssifier" 그대로)
  Type of dup  Code  Number
  Singleton    0     12255
  Dispersed    1     32574
  Proximal     2     2392
  Tandem       3     8998
  WGD or segmental  4  113139
  ```
  Ka/Ks 산출 이후 후처리는 pandas 기반 커스텀 스크립트로 수행:
  ```python
  # KaKs_statistics.py — 종쌍(Zj/Zmt/Zp vs Zsg=Z.sinica)별로 필터링해 개별 tsv로 분리
  inputfile = "11sp_for_MCScanX.fixed2.kaks"
  df = pd.read_csv(inputfile, sep='\t', header=None, comment='#', on_bad_lines='skip')
  df[6] = df[4] / df[5]   # Ka/Ks 비율 계산(4=Ka, 5=Ks 컬럼)
  conditions = {
      "Zj_Zsg": (df[1].str.contains('Zj') & df[2].str.contains('Zsg')),
      "Zmt_Zsg": (df[1].str.contains('Zmt') & df[2].str.contains('Zsg')),
      "Zp_Zsg": (df[1].str.contains('Zp') & df[2].str.contains('Zsg')),
      "Zsg_Zj": (df[1].str.contains('Zs') & df[2].str.contains('Zj')),
      "Zsg_Zmt": (df[1].str.contains('Zs') & df[2].str.contains('Zmt')),
      "Zsg_Zp": (df[1].str.contains('Zs') & df[2].str.contains('Zp'))
  }
  # 조건별로 만족하는 행을 Zoysia/{key}.tsv 로 저장

  # KaKs_count.py — 동일 종쌍 조합에서 Ka/Ks > 1(양성선택 후보)인 유전자쌍 개수만 집계·출력
  ```

### Step 5. Ks 분포 / WGD(전유전체중복) 시각화 (9_Ks_distribution)
- **소프트웨어(버전):** Python 3(pandas, matplotlib, seaborn) — 버전 확인 필요
- **핵심 파라미터:** Ks 값 0.0–2.0(또는 0.0–1.0) 구간을 0.01 단위로 바 플롯 + KDE(bw_adjust=0.5)로 이중
  y축에 중첩 표시
- **Input:** `visualization/11sp.kaks`(Step 4 산출물의 사본, 11종), `Spartina/12sp_including_Spartina.kaks`
  (Spartina 포함 12종 버전)
- **Output:** `Speciation_event_based_on_Ks_value.png`(종분화 이벤트 추정), `WGD_event_based_on_Ks_value.png`
  (전유전체중복 피크 추정)
- **실행 커맨드:**
  ```python
  # 최종 채택본: 9_Ks_distribution/Spartina/Ks_Speciation_plot_histo_and_KDE4.py,
  # Ks_WGD_plot_histo_and_KDE2.py 등 (아래는 핵심 로직 발췌, 반복되는 그룹 루프 생략)
  with open('12sp_including_Spartina.kaks', 'r') as file:
      lines = [line.strip() for line in file if not line.startswith('#')]
  df = pd.DataFrame([line.split() for line in lines])
  ks_values = df[6].astype(float).round(2)
  ks_values_in_range = ks_values[(ks_values >= 0.0) & (ks_values <= 2.0)]
  frequency_counts = ks_values_in_range.value_counts().sort_index()
  percentages = (frequency_counts / frequency_counts.sum()) * 100
  ax1.bar(percentages.index, percentages.values, width=0.01, align='center', alpha=0.5)
  sns.kdeplot(ks_values_in_range, bw_adjust=0.5, fill=True, common_norm=False, alpha=0.5, ax=ax2)
  plt.savefig('Speciation_event_based_on_Ks_value.png', dpi=300)
  ```
  최종 채택본은 `9_Ks_distribution/Spartina/Ks_*_plot_histo_and_KDE{,2,3,4}.py` +
  `WGD_Ks_Peak*.py`(2024-07-10~11, `12sp_including_Spartina.kaks` 입력, Spartina 외군 포함)이며,
  실제 출력 파일(`Speciation_event_based_on_Ks_value.png`, `WGD_event_based_on_Ks_value.png`)의
  최종 mtime이 이 폴더와 일치한다(2024-07-11). 다만 같은 폴더의 `Ks_distribution_TEST.png`는
  이보다 9개월 뒤인 2025-04-08에 생성되어 있어, 이것이 그 이후 추가 재검토/최종 확정본인지 여부는
  확인 필요.

### Step 6. Circos 시각화 (12_Circos)
- **소프트웨어(버전):** Circos(버전 확인 필요) — 디렉토리 depth 1에는 wrapper 스크립트나
  `nohup.out`/`.bash_history`가 전혀 없어(확인함) 정확한 실행 커맨드라인은 확인 필요 — 스크립트
  파일 없음. `.conf` 세트와 산출물(`circos_10000.png/svg`)의 존재로 실행 자체는 확인되나 호출 방식은
  config 파일 경로 지정(`circos -conf main2.conf` 형태로 추정)만 가능하다.
- **핵심 파라미터:** `chromosomes_units = 1000000`(Mb 단위), ideogram thickness 180p, 링크 트랙은
  `ribbon=yes, flat=yes`(공선성 블록을 리본으로 표시)
- **Input:** `karyotype.txt`/`new_karyotype.txt`(HiC_scaffold_chr01–20 염색체 좌표), `Masked_Predicted_genes.txt`
  (예측 유전자 tile 트랙), `LTR_elements.txt`/`LTR_Gypsy.txt`/`LTR_Copia.txt`(LTR 레트로트랜스포존 tile
  트랙), `Zsg.Zsg2.anchors.simple.genomic_block.txt`(링크 트랙 — `MCScan/` 하위 폴더에서 별도로 생성된
  Z. sinica 자기 자신 두 버전 간(Zsg vs Zsg2로 표기, 정확한 실체는 확인 필요) MCScanX anchor 결과를
  Python으로 gff 좌표에 매핑한 산출물)
- **Output:** `circos_10000.png` / `circos_10000.svg`(10kb 해상도 최종 원형 플롯)
- **최종 채택 config:** `main2.conf`(LTR_Gypsy/LTR_Copia 2개 tile 트랙, 20개 염색체 전체 링크 색상
  규칙) — `circos_10000.png`/`.svg`의 mtime(2024-07-18)이 이 config와 같은 날짜다. 이후 색상 투명도만
  일부 수정한 `main3.conf`(2024-07-19)가 존재하나, 이를 반영해 재렌더링한 이미지 파일은 폴더에서
  확인되지 않아 최종 채택 여부는 확인 필요.
- **실행 커맨드:**
  ```conf
  # main2.conf 핵심 발췌 — LTR 트랙 2종 + 20개 염색체 링크 색상 규칙(반복되는 chr03~chr20 규칙 생략)
  karyotype = karyotype.txt
  chromosomes_units = 1000000

  <plot>
  type = tile
  file = Masked_Predicted_genes.txt
  r1 = .99r
  r0 = .79r
  </plot>

  <plot>
  type = tile
  file = LTR_Gypsy.txt
  r1 = .59r
  r0 = .41r
  color = dorange
  </plot>

  <plot>
  type = tile
  file = LTR_Copia.txt
  r1 = .77r
  r0 = .60r
  color = dred
  </plot>

  <links>
  <link>
  file   = Zsg.Zsg2.anchors.simple.genomic_block.txt
  ribbon = yes
  flat   = yes
  <rules>
  <rule>
  condition = var(chr1) eq "HiC_scaffold_chr01"
  color     = chr1
  z         = eval(average(var(size1),var(size2)))
  </rule>
  # ... chr02 ~ chr20 동일 패턴 반복 ...
  </rules>
  </link>
  </links>
  ```
  링크 입력(`Zsg.Zsg2.anchors.simple.genomic_block.txt`)은 `MCScan/Genomic_block_MCScan.py`로 생성된다
  (JCVI MCScan 계열 `.anchors` 파일 + gff를 받아 유전자 ID를 좌표로 치환):
  ```python
  # MCScan/Genomic_block_MCScan.py (대화형 입력 3개: gff, anchors, output)
  gff_df = pd.read_csv(gff_file_path, sep='\t', header=None, comment='#')
  anchors_df = pd.read_csv(anchors_file_path, sep='\t', header=None)
  # anchors 파일의 방향(+/-)에 따라 대응 유전자 4개를 gff에서 찾아 좌표(컬럼 0,3,4)로 치환 후 저장
  ```

## 2. 프로젝트별 특이사항

### Zoysia_sinica_Genome_Assembly
- (단일 프로젝트 카테고리 — 위 공통 파이프라인이 곧 이 프로젝트의 파이프라인)
- 연구목적/샘플 설계: 이 카테고리의 네 분석은 Zoysia_sinica 조립·1차 주석 이후, (1) 근연 Zoysia 3종의
  재주석을 통해 비교 가능한 유전자 세트를 확보하고, (2) 오솔로그/공선성 분석으로 유전자 가족 진화와
  전유전체중복(WGD) 흔적을 탐색하며, (3) 그 결과를 Circos로 시각화해 최종 게놈 구조(반복서열 분포 +
  공선성 블록)를 한눈에 보여주는, "비교유전체 근거 자료 생성" 역할을 담당한다.
- 프로젝트 내부의 분기:
  - OrthoFinder는 최초 11종 실행(`Results_Jun19`)이 CAFE로 이어지는 채택 버전이며, Spartina
    alterniflora를 추가한 12종 실행이 별도로 존재해 Ks/WGD 확장 분석(Step 5)의 입력이 된다(Step 2
    참조).
  - Circos 최종 렌더링은 `main2.conf` 시점의 것으로 판단된다(Step 6 참조).
  - `Zoysia_Gene_reannotation`은 Zj/Zmt/Zp 각각에 대해 "원본 단백질 증거만" 버전과 "RNA-seq/IsoSeq
    증거 추가" 버전을 별도로 실행했고(Zp는 예외적으로 단백질 증거만 확인됨), AGAT 버전(v1.2.0)이
    이 프로젝트의 최종 채택 유전자예측(BRAKER3 v3.0.7 / AGAT v1.4.0, `de_novo_genome_assembly.md` Step 9)과
    다르다 — 이는 이 재주석 배치가 그 최종 실행과는 별개의(아마 비교유전체 전용) 부가 실행이었음을
    시사한다.

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Zoysia_sinica_Genome_Assembly | 근연 Zoysia 3종 재주석(BRAKER2/3, 종별 증거 조합 상이) → OrthoFinder(11종 채택 실행) → MCScanX 공선성/Ka·Ks(11종 전체 blastp, 641,046초 소요) → Circos 시각화 순으로 이어지는 단일 스레드 파이프라인. 오솔로그 클러스터링 결과가 `7_CAFE`(별도 카테고리)로 전달되는 유일한 연결점 | OrthoFinder v2.5.5(확인됨), AGAT v1.2.0(이 재주석 배치, 최종 채택 유전자예측의 v1.4.0과 다름), MCScanX/Circos 버전은 로그에 기록되지 않아 확인 필요 | Spartina alterniflora를 외군으로 포함한 12종 Ks/WGD 확장 분석, Z. sinica 게놈 자기 자신 두 버전 간 self-synteny 기반 Circos 시각화 |
