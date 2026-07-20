# 유전체 주석 및 실험 검증 (Zoysia_sinica_Genome_Assembly, Rubus_Pangenome_Project)

## 1. 공통 파이프라인

이 카테고리는 `de_novo_genome_assembly.md`에 문서화된 최종 산출물 — BRAKER3 유전자 예측 +
Manual_curation을 거친 최종 게놈(`Z.sinica_manualcurated.FINAL_merged2.fasta`, 20개 염색체 스캐폴드)
및 그 유전자모델(`Zsg_manualcurated.gff3`, 대표 단백질 `braker_renamed_longest_isoform.pep`) —
을 공통 입력으로 삼아 갈라져 나가는 6개의 독립적인 후속 분석을 다룬다. 각 분석은 서로 의존하지 않고
병렬로 수행되었으며(EnTAP 기능주석, 침수/과습 후보유전자 비교분석, 텔로미어 탐지, 센트로미어 탐지,
프로모터 시스-작용요소 분석, qRT-PCR 프라이머 설계), 우선순위나 순서를 공유하지 않는다. 다만 일부는
서로 다른 시점의 게놈 버전(원본 contig → superscaffold → manual-curation 후 최종본)을 입력으로
사용했다는 차이가 있어, 이 문서에서는 각 분석이 실제로 어떤 게놈 버전을 썼는지 명시한다.

```
[최종 산출물: BRAKER3 + Manual_curation 게놈/유전자모델 (de_novo_genome_assembly.md)]
        │
        ├──▶ [Step 1. 기능 주석: EnTAP (v1.1.0) + AT/Os 상동성 비교]
        │
        ├──▶ [Step 2. 침수(waterlogging) 후보유전자 비교분석: AP2 도메인 hmmsearch
        │     + 아미노산 치환 비교(Biopython) + KaKs]
        │
        ├──▶ [Step 3. 텔로미어 탐지: quarTeT TeloExplorer
        │     (원본 contig / superscaffold / manual-curation 최종본 3개 버전 각각 실행)]
        │
        ├──▶ [Step 4. 센트로미어 탐지: EDTA TE 주석 + quarTeT CentroMiner]
        │
        ├──▶ [Step 5. 프로모터 시스-작용요소 분석: TSS 상류서열 추출 → PlantCARE(외부 웹툴) 제출]
        │
        └──▶ [Step 6. qRT-PCR 프라이머 설계 및 특이도 검증: Primer3 + BLAST 교차검증]
```

### Step 1. 기능 주석 (EnTAP)
- **소프트웨어(버전):** EnTAP v1.1.0 (내부적으로 DIAMOND blastp, EggNOG-mapper 2.1.12 호출)
- **핵심 파라미터:** `--runP`(단백질 입력이므로 TransDecoder 프레임 선별 단계는 건너뜀), 데이터베이스
  `RefSeq.protein.dmnd` + `uniprot_sprot.dmnd`, `e-value=1e-05`, `qcoverage=50`, `tcoverage=50`,
  threads=100
- **Input:** `braker_renamed_longest_isoform.pep` (BRAKER3 최장 이소폼 대표 단백질, 29,551개 서열)
- **Output:** `final_results/annotated.tsv`, `final_results/annotated.faa`,
  `final_results/annotated_gene_ontology_terms.tsv`, `final_results/entap_results.tsv`,
  `final_results/unannotated.tsv`/`.faa`
- **실행 커맨드:** `entap_run.params`/`entap_config.ini`를 인자로 하는 EnTAP 실행이며, 로그에 기록된
  실제 내부 호출은 다음과 같다.
  ```bash
  # DIAMOND 유사서열 검색 (RefSeq, uniprot_sprot 각각 1회)
  diamond blastp -d RefSeq.protein.dmnd -q transcriptomes/braker_renamed_longest_isoform_final.fasta ...
  diamond blastp -d uniprot_sprot.dmnd -q transcriptomes/braker_renamed_longest_isoform_final.fasta ...

  # 유전자 계열/온톨로지 분석 (EggNOG-mapper, 로그에서 그대로 발췌)
  emapper.py --cpu 100 --dmnd_db /pbi-acc1/ENTAP_DB/eggnog4.clustered_proteins.dmnd \
    --no_file_comments -i transcriptomes/braker_renamed_longest_isoform_final.fasta \
    --sensmode more-sensitive --output_dir gene_family/EggNOG \
    --data_dir /pbi-acc1/ENTAP_DB/ -o blastp_braker_renamed_longest_isoform -m diamond
  ```
- **실행 결과 요약(log_file 발췌):** 입력 29,551서열 중 similarity search 정렬 성공 27,686개(93.69%),
  gene family 할당 26,389개(89.30%), GO term 보유 11,960개(40.47%), KEGG pathway 보유 8,034개
  (27.19%). 전체 러닝타임 217분.
- **주의:** 프로젝트 루트의 최상위 `EnTAP/` 폴더(`EnTAP.sh`, `EnTAP.log`)는 위
  `5_Functional_annotation/`과 **다른, 실패한 실행**이다(Zoysia 5개 종 코딩서열 배치 실행 중 첫
  종에서 TransDecoder 에러로 중단되어 Zoysia_sinica 차례까지 도달하지 못함, 2024-02-06). 실제 채택된
  Zoysia_sinica 기능주석 결과는 `5_Functional_annotation/`의 단일-종 단백질 입력 실행(2024-06-11)이다.
- **부가 분석 — 상동성 비교(AT_Similarity / Os_Similarity):** `5_Functional_annotation/AT_Similarity/`,
  `Os_Similarity/`에서 Zsg(Zoysia sinica) 대표 단백질과 Arabidopsis thaliana(TAIR10),
  Oryza sativa(+ O. indica) 단백질 세트 간 상호 blastp(6fmt)를 수행해 `At2Zsg.6fmt`, `Zs2At.6fmt`,
  `Zsg2Os.6fmt` 등의 결과와 `At_Zsg_Pairwise.ID`(1:1 오솔로그 후보 목록)를 생성했다. 다만 이 두
  하위 폴더에는 실행 스크립트 파일이 남아있지 않고 `nohup.out`만 존재하며, 그마저도
  `AT_Similarity/nohup.out`에는 `-max_target_seq`(오타로 인한 blastp 인자 오류) 경고가 섞여 있어
  최종 성공 커맨드라인 전체를 재구성할 수는 없다 — 확인 필요 — 스크립트 파일 없음(정확한 blastp
  커맨드라인).

### Step 2. 침수(waterlogging) 후보유전자 비교분석
- **소프트웨어(버전):** HMMER3(hmmsearch, Pfam AP2 도메인 PF00847.25 모델), Biopython(AlignIO)
- **핵심 파라미터:** 6개 종(Osativa, Othomaeum, Salterniflora, Zjaponica, Zmatrella, Zpacifica,
  Zsinica) 단백질 세트 각각에 대해 AP2 도메인(PF00847) hmmsearch 수행 → 결과에서 ERFVII 전사인자의
  N-degron 특징 모티프인 "MCGG"(`MCGG_motif.fasta`)를 포함하는 서열만 골라 종별
  `AP2_domain_hits_<species>_MCGG.hits`로 필터링
- **Input:** 6개 종 단백질 FASTA(BLAST db로도 인덱싱됨) + `PF00847.hmm`
- **Output:** `AP2_domain_hits_<species>.txt`(hmmsearch domtblout), `AP2_domain_hits_<species>_MCGG.hits`,
  오솔로그 그룹 정렬 FASTA(`OG0000251.fa`, `OG0002467.fa`, `OG0004300.fa`, `OG0018266.fa`),
  `ERF_Group7.KaKs`/`Nitrogen_metabolism.KaKs`/`Suberin_biosynthesis.KaKs`(dN/dS 값),
  `AtHRGs_Results_*.OGs`/`AtHypoxia_signaling_Results_*.OGs`(Arabidopsis 침수 관련 유전자 ID를
  오솔로그 그룹에 매핑한 결과)
- **실행 커맨드:** hmmsearch 자체를 호출한 스크립트/로그는 폴더에 남아있지 않다 — 확인 필요 —
  스크립트 파일 없음. 다만 `AP2_domain_hits_*.txt`가 HMMER3 domtblout 포맷이고 GA/TC cutoff가
  PF00847 모델의 것과 일치해 hmmsearch 실행 자체는 확실하다. 이후 아미노산 치환 비교는 다음
  Biopython 스크립트로 수행되었다(발췌):
  ```python
  # Amino.py — GWH(공개 유전체 accession)와 Zsg(Zoysia sinica)가 동일하되
  # 다른 종들과는 다른 아미노산 위치를 탐색
  alignment = AlignIO.read("OG0001334.fa", "fasta")
  for record in alignment:
      if record.id.startswith("GWH"):
          gwhp_seq = record.seq
      elif record.id.startswith("Zsg"):
          zsg_seq = record.seq
  # ... gwhp_aa == zsg_aa and other_seqs 모두 다른 위치만 출력

  # Amino_acid_change_Only_Zsg_in_Zoysia.py — Zjn/Zp/Zmt는 동일하지만 Zsg만 다른 위치 탐색
  # Amino_acid_change_GWH_Zsg_shared.py — 위 Amino.py 로직을 OG0*.fa 전체에 대해 배치 실행
  ```
- KaKs(`*.KaKs`) 값을 생성한 스크립트는 폴더 내에서 확인되지 않는다 — 확인 필요 — 스크립트 파일 없음.

### Step 3. 텔로미어 탐지 (TeloExplorer)
- **소프트웨어(버전):** quarTeT TeloExplorer 모듈(파일명 규칙상 확인됨; 정확한 quarTeT 버전은
  확인 필요), 시각화는 R `RIdeogram` 패키지
- **핵심 파라미터:** 게놈 버전별로 검출된 텔로미어 반복 단위(monomer)가 서로 다르다 — 아래 3개 실행
  각각 별도 결과를 냈다.
- **Input / Output (3개 게놈 버전, 3회 독립 실행):**
  1. `Z.sinica_IO_genome.original.fasta`(원본 contig 어셈블리, `original_contig` 접두사) → monomer
     `AAACCCT`, 73개 contig 중 양쪽 텔로미어 모두 검출 0개, 한쪽만 검출 15개, 미검출 58개
  2. `Zs_superscaffolds.nogaps.fasta`(manual curation 이전 superscaffold, `quarTeT` 접두사) →
     monomer `AAAATAT`, 20개 스캐폴드 모두 텔로미어 미검출(양쪽 0, 한쪽 0, 미검출 20)
  3. `Manualcurated/Z.sinica_manualcurated.FINAL_merged2.fasta`(manual curation 최종본,
     `Zs_manualcurated` 접두사) → monomer `AAAAAAG`, 20개 염색체 **전부 양쪽 텔로미어 검출 성공**
     (`Zs_manualcurated.telo.info`: `Both telomere found: 20`)
- **실행 커맨드:** 실행 스크립트/로그(nohup.out 등)가 폴더에 남아있지 않다 — 확인 필요 — 스크립트
  파일 없음. `tmp/` 하위의 `*.telo.chr.txt`, `*.telo.label.txt`, `*.telo.genomedrawer.r`,
  `*_telomeric_repeat_windows.tsv` 파일명 규칙으로 quarTeT TeloExplorer 산출물 포맷임만 확인했다.
  시각화는 다음 R 스크립트로 수행:
  ```r
  library(RIdeogram)
  chr <- read.table("tmp/quarTeT.telo.chr.txt", sep = "\t", header = T, stringsAsFactors = F)
  label <- read.table("tmp/quarTeT.telo.label.txt", sep = "\t", header = T, stringsAsFactors = F)
  ideogram(karyotype = chr, label = label, label_type = "marker", output = "quarTeT.telo.svg")
  convertSVG("quarTeT.telo.svg", file = "quarTeT.telo", device = "png")
  ```
- **해석:** manual curation 이후에만 20개 염색체 전부에서 양쪽 텔로미어가 검출된 것은, 수동 큐레이션
  단계(Manual_curation, `de_novo_genome_assembly.md` 참조)가 스캐폴딩 단계에서 누락/절단되었던
  말단부 서열을 복구했음을 시사한다.

### Step 4. 센트로미어 탐지 (CentroMiner)
- **소프트웨어(버전):** EDTA v2.2.2(TE 주석), quarTeT CentroMiner 모듈
- **핵심 파라미터:** EDTA `--sensitive 1 --anno 1 --overwrite 1 --threads 150`; CentroMiner
  `--TE <EDTA TEanno.gff3> --gene <유전자 GFF3> -t 150 --overwrite`
- **Input:** `Z.sinica_manualcurated.FINAL_merged2.fasta`(manual curation 최종 게놈),
  `Zsg_manualcurated.cds`, `Zsg_manualcurated_longest_isoform.bed`(EDTA `--exclude`용), 유전자
  좌표 GFF3(다른 파이프라인 산출물을 심볼릭 링크로 재사용 — 상세 불필요)
- **Output:** `Z.sinica_manualcurated.FINAL_merged2.fasta.mod.EDTA.TEanno.gff3`(TE 40.36%),
  `.mod.MAKER.masked`(마스킹 18.27%), `TandemRepeat/`, `Candidates/`(염색체별 센트로미어 후보 +
  플롯), `Ideogram/`
- **실행 커맨드:**
  ```bash
  #!/bin/bash
  #Building TE annotation using EDTA
  #EDTA.pl --genome Z.sinica_manualcurated.FINAL_merged2.fasta --cds Zsg_manualcurated.cds \
  #  --exclude Zsg_manualcurated_longest_isoform.bed --overwrite 1 --sensitive 1 --anno 1 --threads 150

  #Centromere prediction using CentroMiner
  #conda activate quartet
  quartet.py CentroMiner -i Z.sinica_manualcurated.FINAL_merged2.fasta \
    --TE Z.sinica_manualcurated.FINAL_merged2.fasta.mod.EDTA.TEanno.gff3_20250201_203533 \
    --gene <유전자 GFF3> -t 150 --overwrite
  ```
- **최종 성공 실행:** `CentroMinder.log3`(2025-02-01) — EDTA TE 40.36%, MAKER masking 18.27%까지 정상
  완료. `CentroMinder_quartet.log`는 이 결과(`*.TEanno.gff3_20250201_203533`)를 입력으로 사용해
  CentroMiner 본 단계를 20개 염색체 전부에 대해 45초 만에 1회 성공으로 완료했다.

### Step 5. 프로모터 시스-작용요소 분석 (PlantCARE)
- **소프트웨어(버전):** 커스텀 awk/bash 스크립트로 TSS 상류서열 추출 후 **PlantCARE는 로컬 실행이
  아닌 외부 웹 도구**(http://bioinformatics.psb.ugent.be/webtools/plantcare) 제출용 입력만 준비함
  — 폴더 내에 PlantCARE 자체 실행 파일/로그가 전혀 없고 제출용 FASTA만 존재하는 것으로 이를 확인함
- **핵심 파라미터:** 전사체(transcript) 시작점 기준 + strand는 상류 2000bp~하류 100bp, - strand는
  상류(좌표상 하류) 100bp~2000bp 구간을 추출
- **Input:** 오솔로그 그룹 FASTA(`OG0019154.fa`, 유전자 ID는 "Zsg" 접두사만 사용),
  `Zsg_manualcurated.gff3`
- **Output:** `OG0019154_Zsg_TSS.bed`(→ PlantCARE 제출용 서열의 좌표), 이와 함께 폴더에는
  `OG0000933_Zsg_TSS.fa`, `OG0002182_Zsg_TSS.fa`, `OG0004646_Zsg_TSS.fa` 등 동일한 성격의 산출물이
  더 있으나(2025-03-16, `TSS_extraction.sh`(2025-05-02)보다 이전 날짜) 그것들을 만든 스크립트는
  남아있지 않다 — 확인 필요 — 스크립트 파일 없음(초기 3개 오솔로그 세트의 추출 스크립트)
- **실행 커맨드:**
  ```bash
  #!/bin/bash
  for x in $(cat OG0019154.fa | grep Zsg | cut -c 2-)
  do
      grep ${x} Zsg_manualcurated.gff3 |
      awk '
      $3 == "transcript" {
          if ($7 == "+") {
              start = ($4 - 2000 < 1) ? 1 : $4 - 2000;
              end = $4 + 100;
              print $1"\t"start"\t"end"\t"$9"\t.\t"$7;
          } else if ($7 == "-") {
              start = ($5 - 100 < 1) ? 1 : $5 - 100;
              end = $5 + 2000;
              print $1"\t"start"\t"end"\t"$9"\t.\t"$7;
          }
      }' >> OG0019154_Zsg_TSS.bed
  done
  ```

### Step 6. qRT-PCR 프라이머 설계 및 특이도 검증
- **소프트웨어(버전):** Primer3(`primer3_core`), samtools(faidx), BLAST+(blastn/tblastn),
  MAFFT + EMBOSS `em_cons`(다중 파랄로그 보존서열 추출용), Python(pandas)
- **핵심 파라미터(Primer3):** `PRIMER_OPT_SIZE=20`, `PRIMER_MIN_SIZE=17`, `PRIMER_MAX_SIZE=25`,
  `PRIMER_PRODUCT_SIZE_RANGE=150-210`(v2 배치 스크립트는 `100-200`), `PRIMER_OPT_TM=58`,
  `PRIMER_MISPRIMING_LIBRARY`로 다른 후보/타 종 서열을 off-target 배경으로 지정
- **Input:** 후보 유전자 CDS(`Zsg_braker_renamed_longest_isoform.cds` 등 5종 Zoysia CDS —
  Zj/Zma/Zmt/Zp/Zsg — 를 합친 `Zoysia5sp.cds`를 off-target 배경으로 사용), `SUB1A-1.fasta`
  (Oryza sativa Indica의 Sub1A-1 단백질, ERFVII 침수내성 유전자 참조서열)
- **Output:** `Primer_designed.txt`(Primer3 raw 출력), `Primer_query.fa`(추출된 LEFT/RIGHT 프라이머
  FASTA), `Primer_query2Z*.6fmt`(각 종 CDS DB 대상 특이도 blastn 결과), `filtered_blast_results.6fmt`
- **실행 커맨드(단일 타깃 기본 파이프라인, 발췌):**
  ```bash
  # 1) 프라이머 설계 (qPCR_Search.sh)
  xargs samtools faidx Zoysia5sp.cds < falseID > Offtarget.fa
  # config에 SEQUENCE_ID / SEQUENCE_TEMPLATE / PRIMER_* 파라미터 기록 후
  primer3_core < config > Primer_designed.txt

  # 2) 결과에서 LEFT/RIGHT 프라이머 서열만 FASTA로 추출 (parse_primer_sequences.sh, awk)
  #    Primer_designed.txt → Primer_query.fa

  # 3) 종간 교차서열 확인 (Zoysia2At.tblastn.sh)
  tblastn -db Zsg_braker_renamed_longest_isoform -query Osativa.pep -outfmt 6 \
    -evalue 1e-5 -num_threads 200 -out Os2Zsg.6fmt
  # (Zj/Zma/Zmt/Zp DB 대상으로도 동일하게 반복)

  # 4) 프라이머 특이도 정렬결과를 5' /3' 말단 미스매치 기준으로 정렬 (blastn_sort_by_criteria.sh)
  bash blastn_sort_by_criteria.sh <blast_results.txt> > sorted_results.txt

  # 5) LEFT/RIGHT 모두 동일 subject에 히트하는 프라이머 쌍 필터링 (filter_LEFT_RIGHT_match.py, pandas)
  ```
- **변형 — 다중 파랄로그 보존서열 기반 배치 설계(`Candidates/`):** 여러 파랄로그(중복 유전자)가 있는
  타깃은 MAFFT로 정렬 후 `em_cons`로 보존서열을 얻고, 결과의 축퇴 코드(degenerate code)를 'N'으로
  치환(`VariantsToNNN.py` — 이 스크립트 자체는 `Candidates/` 폴더에 없어 확인 필요 — 스크립트 파일
  없음)한 뒤 `qPCR_Search_v2.sh`로 다건 배치 프라이머 설계(백그라운드 병렬 실행, product size
  100-200bp)를 수행했다:
  ```bash
  # Conserved_templates.sh
  for x in *Zoysia.cds
  do
      mafft --auto ${x} > ${x%.cds}_mafft.fa
      em_cons -sequence ${x%.cds}_mafft.fa -outseq ${x%.cds}_cons.fa -name ${x%.cds}_cons
      python VariantsToNNN.py ${x%.cds}_cons.fa ${x%.cds}_NNN.fa
  done
  ```
- **Zsg_Transcriptome_Target/** 하위 폴더는 위와 동일한 스크립트 세트(`qPCR_Search.sh`,
  `parse_primer_sequences.sh`, `find_min_match_pair.sh`, `blastn_sort_by_criteria.sh`)를 재사용해
  전사체 타깃 특이적으로 재실행한 결과이며, `filter_LEFT_RIGHT_match.py`(pandas 기반 LEFT/RIGHT 쌍
  필터링)가 추가되어 있다.

## 2. 프로젝트별 특이사항

### Zoysia_sinica_Genome_Assembly
- (단일 프로젝트 카테고리 — 위 공통 파이프라인이 곧 이 프로젝트의 파이프라인)
- 연구목적/샘플 설계: 완성된 염색체 수준 게놈에 대해 기능 주석(EnTAP)과 함께 침수/과습 스트레스
  내성 후보유전자(ERFVII/AP2 계열, Sub1A 상동) 탐색, 텔로미어·센트로미어 등 구조적 완결성 검증,
  프로모터 시스-작용요소 및 qRT-PCR 검증용 프라이머까지 마련해 게놈 조립·주석의 품질과 생물학적
  활용 가능성을 동시에 뒷받침하려는 목적으로 판단된다(description.txt/README 부재로 인해 폴더/파일명
  기반 추정 — `de_novo_genome_assembly.md`에서 이미 명시된 한계와 동일).
- 프로젝트 내부의 분기:
  - 실제 채택된 기능주석 결과는 `5_Functional_annotation/`(단일종 단백질 입력)이며, 최상위
    `EnTAP/`은 별도의 실패한 시도다(Step 1 참조).
  - 텔로미어 탐지는 세 게놈 버전(원본 contig/superscaffold/manual-curation 최종본)에 대해 각각
    독립 실행되었고, manual curation 이후에만 완전한 양쪽 텔로미어 검출에 성공했다(Step 3).
  - qRT-PCR 프라이머 설계는 단일 타깃 기본 스크립트 세트 외에, 파랄로그가 여러 개인 타깃을 위한
    MAFFT 보존서열 기반 배치 설계 변형(`Candidates/`)과 전사체 특이적 재실행(`Zsg_Transcriptome_Target/`)
    두 갈래로 확장되었다(Step 6).

### Rubus_Pangenome_Project
- 연구목적/샘플 설계: 4개 근연 Rubus 종(Rco/Rcra/Rpho/Rta) × hap1/hap2, 총 8개
  haplotype-resolved 어셈블리를 대상으로 한 4종 pangenome 비교유전체 연구. 반복서열/유전자 주석
  이후 게놈 완결성 QC(LAI), flavonoid 생합성 경로 효소 동정(KIPEs), 종별 RNA-seq 발현분석
  (Gene_Expression), ellagic acid 생합성 경로 유전자 분석, 공개 유전체 데이터베이스 제출
  (Public_Upload)로 이어진다.
- 공통 파이프라인과의 차이점: 반복서열 마스킹+유전자예측 단계가 두 개의 브랜치
  (Repeat_elements+RepeatMasker_Braker3 vs panEDTA+EDTA_Braker3)로 병행 실행된 뒤, 이후의
  LAI/KIPEs/Gene_Expression/Ellagic acid 분석들이 그 채택된 유전자주석 산출물을 재사용한다.

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Zoysia_sinica_Genome_Assembly | 최종 주석 게놈에서 fan-out하는 6개 독립 분석(EnTAP 기능주석, 침수 후보유전자, 텔로미어, 센트로미어, PlantCARE 프로모터 분석, qRT-PCR 프라이머 설계). 최상위 `EnTAP/`은 별도 실패 실행이며 실제 채택 결과는 `5_Functional_annotation/`. TeloExplorer는 manual curation 이후 버전에서만 완전 검출 | EnTAP v1.1.0, EggNOG-mapper 2.1.12, EDTA v2.2.2, quarTeT(TeloExplorer/CentroMiner, 버전 확인 필요) | AP2 도메인(PF00847)+MCGG 모티프 기반 ERFVII 침수내성 후보유전자 비교(Amino.py류 Biopython 스크립트, KaKs), TSS 상류서열 추출 후 PlantCARE(외부 웹툴) 시스-작용요소 분석, MAFFT 보존서열 기반 다중 파랄로그 qRT-PCR 배치 설계 |
| Rubus_Pangenome_Project | 반복서열 마스킹+유전자예측이 RepeatMasker_Braker3/panEDTA+EDTA_Braker3 두 브랜치로 병행 실행된 뒤 후속 분석이 채택된 산출물을 재사용하는 순차 구조 | (개별 툴 버전은 원본 스크립트에서 확인) | Flavonoid 생합성 경로 효소 동정(KIPEs), Ellagic acid 생합성 경로 유전자 중복/발현 분석, 종별 RNA-seq 발현분석(Gene_Expression), 공개 유전체 데이터베이스 제출(Public_Upload) |
