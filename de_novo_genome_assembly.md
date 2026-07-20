# De novo Genome Assembly (Cannabis_sativa_cv.ChungSam_dnDNAseq, Pichia_Pastoris_dnDNAseq_HanHwa, Abies_koreana_dnDNAseq, Zoysia_macrostachya_dnDNAseq, Zoysia_sinica_Genome_Assembly, Rubus_Pangenome_Project)

## 1. 공통 파이프라인

5개 프로젝트 모두 "원시 long-read(+보조 short-read) → de novo assembly → genome QC" 라는 핵심 골격은
공유하지만, 그 이후 단계(scaffolding 방식, 비교유전체, 유전자주석 수행 여부)는 프로젝트마다 실제로
수행된 범위가 크게 다르다. 특히:

- **Cannabis, Pichia, Zoysia_macrostachya, Zoysia_sinica** 4개는 "assembly → scaffolding → genome QC →
  비교유전체(synteny/dot-plot) → 유전자/반복서열 주석"까지 완결된 파이프라인을 보유한다. 다만
  scaffolding 방식은 프로젝트마다 다르다(Cannabis/Pichia는 참조 유전체 기반 RagTag, Zoysia_macrostachya는
  Omni-C/Hi-C 기반 de novo 스캐폴딩(HapHiC) + 별도 갭 클로징(YagCloser), Zoysia_sinica도 같은
  Omni-C/Hi-C 기반이지만 실제 사용 툴은 다르다 — HapHiC이 아니라 yahs/SALSA2/3D-DNA+Juicer 세 가지를
  병행 시도한 뒤 yahs 계열을 채택했다. 두 Zoysia 프로젝트가 "Hi-C 스캐폴딩을 쓴다"는 점은 같지만
  구체적 툴은 다르므로 혼동하지 않도록 주의).
- **Abies_koreana**는 2026-07 기준 아직 초기 단계 프로젝트로, k-mer profiling → de novo assembly →
  organelle contig 분리(사후 필터링) → 1차 QUAST QC까지만 완료되어 있다. scaffolding, 비교유전체,
  유전자주석 단계는 아직 수행되지 않았으며, BUSCO도 사용되지 않고 QUAST만으로 QC가 이루어진다. 즉
  Abies는 나머지 프로젝트들과 "공통 뼈대"의 앞부분(조립까지)만 공유하고 뒷부분은 아직 해당사항이
  없다 — 이는 프로젝트 진행 단계 차이이지 방법론적 차이가 아니다.
- k-mer 기반 genome survey(FastK+GenomeScope) 단계는 Cannabis/Pichia/Abies 3곳에서 확인되나,
  Zoysia_macrostachya는 이 단계에 해당하는 스크립트/로그가 발견되지 않는다(genome_size=0.323g 값이
  NextDenovo cfg에 직접 기재되어 있으나 그 산출 근거는 확인 필요). Zoysia_sinica도 이 단계에 해당하는
  스크립트/로그가 확인되지 않았다(원본 조립(Original/) 자체가 접근 불가였으므로 genome survey 수행
  여부조차 확인 불가 — 확인 필요).
- organelle contig 분리(Oatk) 단계는 Cannabis와 Abies에만 있다(둘 다 PacBio HiFi 단일 롱리드
  플랫폼을 사용하는 프로젝트라는 공통점이 있음).
- Zoysia_sinica는 description.txt/README가 전혀 없어(확인됨) 이 문서의 연구목적 서술이 폴더/파일명
  기반 추정에 의존하는 유일한 프로젝트이며, 실제 원시 de novo assembly 단계(organelle 분리 → 조립 →
  1차 QC → MAKER → 1차 기능주석 → 1차 RepeatModeler/Masker)가 수행된 `Original/ONT_longread_DNAseq/`
  디렉토리 자체가 조사 시점에 접근 불가(간헐적 NFS 이슈)였다 — 따라서 이 프로젝트는 나머지 4개보다
  검증되지 않은 항목("확인 필요")이 훨씬 많다. 이는 문서화의 결함이 아니라 실제 접근 가능한 증거의
  한계를 정직하게 반영한 것이다.

아래 DAG는 이 5개 프로젝트에서 실제로 확인된 단계만을 표시하며, 특정 프로젝트에만 해당하는 단계는
괄호로 명시한다.

```
[rawdata: long-read (+ short-read)]
        │
        ▼
[Step 1. Genome survey / k-mer profiling: FastK + GenomeScope(2)]  (Cannabis/Pichia/Abies만 해당 — Zoysia_macrostachya/Zoysia_sinica 둘 다 미확인)
        │  (게놈 크기·이형접합률 추정치 → 이후 assembly 파라미터로 재사용)
        ▼
[Step 2. De novo assembly (프로젝트별 도구/전략 상이 — 2절 참조)]
        │  (Zoysia_sinica: 실제 조립이 수행된 Original/ 디렉토리 자체가 접근 불가 — 확인 필요, 추정 금지)
        │
        ├──▶ [Step 3. Organelle contig 분리(Oatk) — Cannabis/Abies만 해당]
        │
        ├──▶ [Step 4. 참조 기반 scaffolding: RagTag — Cannabis/Pichia만 해당]
        │        (Zoysia_macrostachya/Zoysia_sinica도 RagTag를 쓰지만 "1차 scaffolding" 용도가 아님 —
        │         각각 Step 8/Step 5·10 참조)
        │
        ├──▶ [Step 5. Hi-C(Omni-C) 기반 de novo 스캐폴딩 — Zoysia_macrostachya/Zoysia_sinica만 해당
        │     (스캐폴딩 툴 자체는 상이: macrostachya=HapHiC, sinica=yahs/SALSA2/3D-DNA 3종 병행 후 yahs 채택)]
        │            │
        │            ├──▶ [Step 6. Gap closing: YagCloser — Zoysia_macrostachya만 해당]
        │            │
        │            └──▶ [Step 4b. RagTag 참조기반 재정렬/명명(Z. japonica 대비 chr01-20) —
        │                   Zoysia_sinica만 해당, scaffolding 자체가 아니라 renaming/reordering 목적,
        │                   중간 1회]
        │
        ▼
[Step 7. Genome QC: QUAST(5개 공통, Zoysia_sinica는 스캐폴더 3종 비교용으로도 사용) +
 BUSCO(Cannabis/Pichia/Zoysia만, Abies는 미수행) + LAI(Zoysia_macrostachya만)/Merqury(Zoysia_macrostachya/Zoysia_sinica)]
        │
        ├──▶ [Step 8. 비교유전체 분석 (synteny/dot-plot/pan-genome, Abies는 아직 미수행)]
        │
        └──▶ [Step 9. 유전자/반복서열 주석 (de novo 예측 또는 참조 기반 liftover, Abies는 아직 미수행)]
                        │
                        ▼
                [Step 10. 수동 큐레이션(Manual_curation) + 최종 RagTag 재정렬(2차) —
                 Zoysia_sinica만 해당]
                        │
                        ▼
                [프로젝트 고유 확장 분석 / 최종 산출물]
```

### Step 1. Genome survey (k-mer profiling)
- **소프트웨어(버전):** FastK, GenomeScope2 (Cannabis: 버전 확인 필요; Pichia: 사용 확인되었으나 세부 버전은 재현 우선순위 낮음으로 처리됨; Abies: FastK/Fastmerge/Histex/GenomeScope 사용 확인되었으나 FastK가 `-h`/`--version` 플래그를 지원하지 않아 버전 확인 불가) — Zoysia_macrostachya/Zoysia_sinica 둘 다 이 단계에 해당하는 스크립트/로그 없음(확인 필요; Zoysia_sinica는 원본 조립 디렉토리 자체가 접근 불가였으므로 수행 여부조차 확인 불가)
- **핵심 파라미터:** k=31(Cannabis, HiFi reads 대상) / k=21(Pichia, Illumina reads 대상) / k=40(Abies, HiFi reads 대상, `-M64 -T32`, GNU parallel `-j 8`)
- **Input:** 원시 long-read(Cannabis, PacBio HiFi) 또는 QC된 Illumina reads(Pichia) 또는 원시 PacBio Revio HiFi reads(Abies, JSlink 4cell 단독)
- **Output:** k-mer histogram, genome size/heterozygosity 추정치(GenomeScope summary)
- **실행 커맨드:**
  - Cannabis (`<sample>`=Female/Male):
    ```bash
    FastK -k31 -M256 -T125 -t -P/pbi-acc1/hsPark/Cannabis_sativa_cv.ChungSam_dnDNAseq/Genome_survey Chungsamsoo-Tree-<sample>.1Cell.hifi.fastq.gz
    Histex -G Chungsamsoo-Tree-<sample>.1Cell.hifi.hist > Csativa_CS_<sample>_31mer_genomescope.hist
    GeneScopeFK.R -i Csativa_CS_<sample>_31mer_genomescope.hist -o Csativa_CS_<sample>_31mer_genomescope -k 31
    ```
  - Pichia (`<sample>`=BG-10/X-33):
    ```bash
    FastK -k21 -M128 -T200 -v -t -P/pbi-acc1/hsPark/Pichia_Pastoris_dnDNAseq_HanHwa/FastK/ <sample>_good_<R>.fastq
    Fastmerge -ht -T200 -P/pbi-acc1/hsPark/Pichia_Pastoris_dnDNAseq_HanHwa/FastK/ <sample>_merged_k21 <sample>_good_1.ktab <sample>_good_2.ktab
    Histex -G <sample>_merged_k21 > <sample>_merged_k21_genomescope.hist
    GeneScopeFK.R -i <sample>_merged_k21_genomescope.hist -o P.pastoris_<sample>_illumina_merged_k21_genomescope -k 21
    ```
  - Abies (K-mer_distribution_JSLink/FastK.sh, JSlink 4cell 단독):
    ```bash
    ls m84210_2603*.fastq.gz > input_files.txt
    cat input_files.txt | parallel -j 8 "FastK -k40 -M64 -T32 -P/pbi-acc2/hsPark/Abies_koreana_dnDNAseq/K-mer_distribution/tmp -v -t {}"
    Fastmerge -ht -T64 -P/pbi-acc2/hsPark/Abies_koreana_dnDNAseq/K-mer_distribution/tmp/ Abk_merged_k40 <4개 fastq prefix>
    Histex -G Abk_merged_k40 > Abk_merged_k40_genomescope.hist
    GeneScopeFK.R -i Abk_merged_k40_genomescope.hist -o genomescope_out -k 40
    ```
    (JSlink+EYEONCELL 8cell 통합본은 동일 스크립트 패턴을 EYEONCELL fastq.gz 심볼릭링크까지 포함해 실행 중이며, 2026-07-13 기준 아직 완료되지 않음 — 정상 진행중)
  - Zoysia_macrostachya: 확인 필요 — 스크립트 파일 없음
  - Zoysia_sinica: 확인 필요 — 스크립트 파일/로그 없음(원본 조립 디렉토리 접근 불가로 이 단계 수행 여부 자체가 불명)

### Step 2. De novo assembly
- **소프트웨어(버전):** Cannabis: hifiasm 0.25.0-r726 / Pichia: NextDenovo v2.5.2 + NextPolish v1.4.1 / Abies: hifiasm v0.25.0-r726(Cannabis와 동일 버전) / Zoysia_macrostachya: Ratatosk(하이브리드 보정, 버전 확인 필요) + NextDenovo(버전 확인 필요) + NextPolish(버전 확인 필요) + ntLink(스캐폴딩, 버전 확인 필요) / Zoysia_sinica: **확인 필요 — 검증 불가**. 실제 조립이 수행된 `Original/ONT_longread_DNAseq/` 디렉토리가 조사 시점 접근 불가(간헐적 NFS 이슈)였기 때문에 어셈블러/버전/파라미터를 확인할 수 없다. 정황(간접 증거)만 남아있다: `Genome_Assembly/3_Repeat_element_annoation/module09_RepeatModeler.sh`가 `Z.sinica_IO_genome.original.fasta`(contig-level)를 입력으로 사용하고 RepeatMasker 라이브러리 디렉토리명이 `RM_20.FriNov170056382023`(2023-11-17경)인 것으로 미루어, 이 시점 이전에 이미 contig-level de novo assembly가 완료되어 있었을 것으로 강하게 추정되나 어셈블러 자체(예: NextDenovo 등 추정)는 추측하지 않는다.
- **핵심 파라미터:** Cannabis: `-t 100`(스레드) / Pichia: `input_type=corrected, read_type=ont, genome_size=<GenomeScope 추정치>`(NextDenovo), `task=best, sgs_options=-max_depth 100 -bwa`(NextPolish, Illumina 폴리싱만 사용) / Abies: `-t 128`(JSlink 4cell 단독) 또는 `-t 256`(JSlink+EYEONCELL 8cell 통합) / Zoysia_macrostachya: Ratatosk `-c 256`; NextDenovo run.cfg(`read_type=ont, genome_size=0.323g, pa_correction=7, minimap2_options_cns=-t 100 -k 19, nextgraph_options=-a 1`, `job_type=local`); NextPolish(`sgs_options=-max_depth 100 -bwa, polish_options=-10`); ntLink(`k=40 w=500`, 7개 k/w 조합 중 최종 채택) / Zoysia_sinica: 확인 필요(원본 스크립트 접근 불가)
- **Input:** Cannabis: organelle 제거 후 nuclear-only HiFi reads / Pichia: ONT reads를 Illumina로 하이브리드 교정(Ratatosk)한 corrected long reads + Illumina reads(폴리싱용) / Abies: organelle 제거 이전 raw 전체 HiFi reads(JSlink 4개 cell, 통합본은 EYEONCELL 1개 cell 추가) / Zoysia_macrostachya: ONT reads(Z.macrostackya_ONT.fastq.gz)를 Illumina reads(Zma-H1_1/2.fastq.gz)로 Ratatosk 보정한 out_long_reads.fastq / Zoysia_sinica: 확인 필요 — 원본 raw long-read 플랫폼/파일 자체가 미확인(Original/ 접근 불가)
- **Output:** Cannabis: haplotype-resolved contig(hap1/hap2 × Female/Male) / Pichia: 폴리싱된 단일(haploid) genome fasta / Abies: 단일 개체(제주도 표준목)의 hifiasm bp.{p_ctg,hap1.p_ctg,hap2.p_ctg} contig(organelle 포함, 이후 Step3에서 제거) / Zoysia_macrostachya: genome.nextpolish.fasta → ntLink scaffold(k40w500) 결과물인 `Ratatosk_NextDP_NextPsh_ntLinkk40w500.fasta`(1차 QC용 조립본, 이후 Purge_dups/Omni-C 스캐폴딩으로 이어짐) / Zoysia_sinica: `Z.sinica_IO_genome.original.fasta`(contig-level, 존재 자체는 파일로 확인되나 산출 스크립트는 미확인 — 확인 필요). 이 contig 산출물에 대해 Hi-C 매핑(Step 5) 이전인 2023-11경 이미 한 차례 contig-level RepeatModeler/RepeatMasker가 수행된 정황이 있다(상세는 Step 9 하단 참고 노트 참조)
- **실행 커맨드:**
  - Cannabis (`<sample>`=Female/Male; hifiasm.log의 `CMD:` 라인에서 확인):
    ```bash
    hifiasm -o Csativa_nuclear_<sample>.asm -t 100 <sample>_nuclear_reads.fastq.gz
    ```
  - Pichia:
    ```bash
    nextDenovo nextdenovo_run.cfg   # cfg: input_type = corrected, read_type = ont, genome_size = 10923843
    nextPolish nextpolish_run.cfg   # cfg: task = best, sgs_options = -max_depth 100 -bwa
    ```
  - Abies (JSlink 4cell 단독, hifiasm_JSLink/01_hifiasm.sh):
    ```bash
    hifiasm -t 128 -o abies_koreana \
        m84210_260317_114629_s1.hifi_reads.bc2009.fastq.gz \
        m84210_260320_095126_s3.hifi_reads.bc2009.fastq.gz \
        m84210_260320_115450_s4.hifi_reads.bc2009.fastq.gz \
        m84210_260327_044805_s2.hifi_reads.bc2009.fastq.gz \
        2> abies_koreana.hifiasm.log
    awk '/^S/{print ">"$2"\n"$3}' abies_koreana.bp.p_ctg.gfa | fold -w 80 > abies_koreana.bp.p_ctg.fa
    ```
    (JSlink+EYEONCELL 8cell 통합본, hifiasm_JSLink_EYEONCELL/01_hifiasm_combined.sh: 동일 커맨드에 `-t 256`, EYEONCELL BAM→FASTQ 변환분 `CNU_AK_WGS_merged.hifi_reads.fastq.gz`(`samtools fastq -@ 256 <bam> | pigz -p 256 > <out>`로 생성) 1개 fastq 추가)
  - Zoysia_macrostachya (01.Assembly/Ratatosk/Ratatosk.sh + run.cfg + run_NextPolish.cfg + ntLink/ntLink.sh):
    ```bash
    Ratatosk correct -v -c 256 -s short_reads.fastq.gz -l in_long_reads.fastq.gz -o out_long_reads
    nextDenovo run.cfg   # read_type=ont, genome_size=0.323g, pa_correction=7, minimap2_options_cns=-t 100 -k 19
    nextPolish run_NextPolish.cfg   # sgs_options=-max_depth 100 -bwa, polish_options=-10
    ntLink scaffold target=genome.nextpolish.fasta reads=out_long_reads.fastq t=256 overlap=True k=40 w=500
    ```
  - Zoysia_sinica: 확인 필요 — `Original/ONT_longread_DNAseq/` 자체가 접근 불가라 실행 스크립트를 열람할 수 없음. 추측하지 않음.
  - Rubus_Pangenome_Project (Rco/Rcra/Rpho 공통 패턴, `<species>`=Rco/Rcra/Rpho; Rta는 별도
    파이프라인으로 Step 5에서 다룸):
    ```bash
    hifiasm -t 100 -l 3 --h1 <species>_Omni-C_R1.fastq.gz --h2 <species>_Omni-C_R2.fastq.gz \
        -o <species>.asm_Omni-C --telo-m TTTAGGG <species>_DNA.hifi.fastq.gz
    ```
    hifiasm v0.25.0-r726(로그에서 직접 확인), Omni-C reads를 `--h1/--h2`로 투입해 haplotype 분리
    조립, purge level 3(`-l 3`) 사용(2절 참조).

### Step 3. Organelle contig 분리(Oatk) — Cannabis, Abies, Rubus_Pangenome_Project만 해당
- **소프트웨어(버전):** Oatk v1.0(양쪽 공통) + nhmmscan(HMMER 3.4, Abies) + minimap2 2.30-r1287/seqkit v2.3.0(organelle contig 사후 필터링, Abies)
- **핵심 파라미터:** Abies: `-t 128`(단독)/`-t 256`(통합), DB=`acrogymnospermae_{mito,pltd}.fam`(v20230921); 필터링 단계 minimap2 `-x asm5 --secondary=no`, `COV_THRESHOLD=0.90`(organelle 참조에 90% 이상 커버리지로 정렬되는 contig만 제거, NUMT/NUPT로 추정되는 낮은 커버리지 contig는 보존)
- **Input:** Abies: (standalone) 원시 HiFi fastq 4개 또는 (그래프 기반, 보조/검증용) hifiasm p_utg.gfa + Step2 산출 bp.{hap1,hap2,p}.p_ctg.fa
- **Output:** Abies: `{pltd,mito}.ctg.fasta/.gfa/.bed`(organelle 서열) → `organelle_refs.fa`(mito+pltd 병합) → `*.nuclear.fa`(organelle 제거된 최종 핵 게놈 contig, hap1/hap2/p 3종)
- **실행 커맨드:**
  - Abies (standalone, Organelle_JSLink/oatk.sh):
    ```bash
    oatk -t 128 --nhmmscan /usr/local/bin/nhmmscan \
        -m /usr/local/bin/OatkDB/v20230921/acrogymnospermae_mito.fam \
        -p /usr/local/bin/OatkDB/v20230921/acrogymnospermae_pltd.fam \
        -o abies_koreana \
        m84210_260317_114629_s1.hifi_reads.bc2009.fastq.gz m84210_260320_095126_s3.hifi_reads.bc2009.fastq.gz \
        m84210_260320_115450_s4.hifi_reads.bc2009.fastq.gz m84210_260327_044805_s2.hifi_reads.bc2009.fastq.gz
    ```
  - Abies (contig 필터링, hifiasm_JSLink/04_filter_organelle_contigs.sh, `<hap>`=hap1/hap2/p 3종 반복):
    ```bash
    cat abies_koreana.pltd.ctg.fasta abies_koreana.mito.ctg.fasta > organelle_refs.fa
    minimap2 -x asm5 -t 64 --secondary=no organelle_refs.fa abies_koreana.bp.<hap>.p_ctg.fa 2> minimap2_<hap>.log \
        | awk -v thr=0.90 '{aln[$1]+=($4-$3); len[$1]=$2} END{for(c in aln) if(aln[c]/len[c]>=thr) print c}' | sort -u > organelle_contigs_<hap>.txt
    seqkit grep -v -f organelle_contigs_<hap>.txt abies_koreana.bp.<hap>.p_ctg.fa -o abies_koreana.bp.<hap>.p_ctg.nuclear.fa
    ```
  - Cannabis: 확인 필요 — Extract_unmapped_reads.sh(minimap2+samtools 사용 정황만 확인, 스크립트 파일 자체는 확인되지 않음; 상세는 2절 참조)
  - Rubus_Pangenome_Project (4종 공통 패턴, `<species>`=Rco/Rcra/Rpho/Rta):
    ```bash
    oatk -t 100 --nhmmscan /usr/local/bin/nhmmscan \
        -m /usr/local/bin/OatkDB/v20230921/magnoliopsida_mito.fam \
        -p /usr/local/bin/OatkDB/v20230921/magnoliopsida_pltd.fam \
        -o <species>_organelle <species>.fasta.gz
    ```
    oatk v1.0, OatkDB v20230921 magnoliopsida mito/pltd 모델 — 4종 공통(2절 참조).

### Step 4. 참조 기반 scaffolding(RagTag) — Cannabis, Pichia만 해당
- **소프트웨어(버전):** RagTag v2.1.0(Cannabis) / RagTag(Pichia, 버전 확인 필요) — 내부적으로 둘 다 minimap2 사용(Cannabis asm5 프리셋, Pichia v2.24-r1122 asm5 프리셋)
- **핵심 파라미터:** `-t 200`(Cannabis)
- **Input:** 조립 결과(contig-level) + 참조 유전체(Cannabis: 판게놈 유래 AH3Ma/AH3Mb; Pichia: NCBI GS115 GCF_000027005.1)
- **Output:** scaffold fasta(.agp, .stats 포함), 이후 chromosome/pseudochromosome 수준으로 재정리
- **실행 커맨드:**
  - Cannabis (동일 패턴으로 4개 haplotype 조합 반복, ragtag.sh 참조):
    ```bash
    ragtag.py scaffold AH3Mb.softmasked.fasta ../Hifiasm/Male/Csativa_CS_Male.asm.bp.hap1.p_ctg.fa -o CS_Male_AH3Mb_Y -t 200
    ragtag.py scaffold AH3Ma.softmasked.fasta ../Hifiasm/Female/Csativa_CS_Female.asm.bp.hap1.p_ctg.fa -o CS_Female_AH3Ma_hap1 -t 200
    ```
  - Pichia: 확인 필요 — 스크립트 파일 없음

> 참고: Zoysia_macrostachya/Zoysia_sinica 둘 다 RagTag를 사용하지만, 이 프로젝트의 RagTag는 위 Step 4의
> "1차 scaffolding" 용도가 아니다.
> - Zoysia_macrostachya의 RagTag(Pan-genome_analysis/RagTag/)는 자기 게놈의 scaffolding이 아니라 타
>   Zoysia속 게놈과의 비교유전체용 재정렬이다. Zoysia_macrostachya 자신의 scaffolding은 Step 5(Hi-C
>   기반)를 통해 이루어지므로 이 RagTag 사용은 Step 8(비교유전체 분석)에 별도로 기재한다.
> - Zoysia_sinica의 RagTag도 자신의 1차 scaffolding이 아니라는 점은 같지만 용도는 다르다 — Z. japonica
>   기준 chr01~chr20 명명/순서·방향 보정("renaming/reordering")을 위한 것으로, 총 두 차례 실행된다
>   (Hi-C 스캐폴딩 직후 1차, Manual_curation 이후 최종 게놈에 대해 2차). 이는 Zoysia_sinica 자신의 핵심
>   파이프라인 안에 있는 스텝이므로 Step 5의 project-labeled sub-block(1차) 및 Step 10(2차)에 각각
>   기재한다.
> Abies는 아직 scaffolding을 수행하지 않았다.

### Step 5. Hi-C(Omni-C) 기반 de novo 스캐폴딩 — Zoysia_macrostachya/Zoysia_sinica/Rubus_Pangenome_Project만 해당
세 프로젝트 모두 Omni-C(Hi-C)를 이용한 de novo 스캐폴딩이라는 접근 자체는 같지만, 실제 사용한
스캐폴딩 툴은 서로 다르다 — Zoysia_macrostachya는 HapHiC 단일 경로, Zoysia_sinica는 yahs/SALSA2/
3D-DNA+Juicer 세 가지를 병행 시도한 뒤 yahs 계열을 채택했다(둘을 같은 툴을 쓴 것으로 오해하지 말 것).
Rubus_Pangenome_Project는 Rco/Rcra/Rpho 3종이 YaHS를, Rta는 AutoHiC를 사용한다(아래 project-labeled
sub-block 참조).
- **소프트웨어(버전):**
  - Zoysia_macrostachya: purge_dups(버전 확인 필요, 실행 스크립트 부재로 config.json만 확인) → bwa + pairtools + samtools(Omni-C 리드 검증) → HapHiC(`haphic pipeline`, 버전 확인 필요) → Juicebox GUI(수동 curation, juicer_tools 1.9.9/jcuda 0.8)
  - Zoysia_sinica: bwa mem + pairtools + samtools(Omni-C 매핑, 원본 `.sh` 미보존·로그만 확인) → yahs / SALSA-2.3(`/usr/local/bin/SALSA-2.3`, python2 기반) / juicer.sh+3D-DNA(3종을 병행 실행 후 QUAST 비교로 yahs 계열 채택 — Step 7 참조) → (채택된 yahs 계열 스캐폴드에 대해) RagTag 참조기반 재정렬/명명(Z. japonica 대비, 1차, 상세는 아래 실행 커맨드 참조 및 Step 4 상단 참고 블록)
  - Rubus_Pangenome_Project(Rco/Rcra/Rpho): bwa 0.7.19-r1273 + samblaster 0.1.26(Omni-C 정렬/중복마킹) → HapHiC `utils/filter_bam`(BAM 필터링 유틸리티만 사용) → YaHS(no_ec 옵션). Rta는 AutoHiC(Juicer v1.6+BWA 0.7.17-r1188+3D-DNA 180922+딥러닝 모델) 별도 경로.
- **핵심 파라미터:**
  - Zoysia_macrostachya: Purge_dups `config.json`: `isdip=1, ispb=1`, busco `lineage=Poales`; bwa `-5SP -T0 -t 150`; pairtools parse `--min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30`; HapHiC `haphic pipeline purged.fa HiC.filtered.bam 20 --correct_nrounds 2 --max_inflation 5`(20=예상 염색체 수); juicer_tools `-Xmx32G`
  - Zoysia_sinica: bwa mem `-5SP -T0`(재구성값); pairtools parse `--min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30`(재구성값); yahs/SALSA2/3D-DNA 각 툴 자체의 실행 파라미터는 원본 스크립트가 보존되어 있지 않아 확인 필요; RagTag 재정렬은 `reverse_complement_chromosomes.py`(Biopython)에 역상보 대상(13개 염색체: chr01,02,03,06,07,08,09,10,11,14,15,17,19)·순서 리스트가 하드코딩됨
  - Rubus_Pangenome_Project(Rco/Rcra/Rpho): bwa mem `-5SP -t 100`; samtools view `-F 3340`(secondary/supplementary/dup 제외); HapHiC filter_bam `1 --nm 3 --threads 100`; YaHS `--telo-motif AAACCCT --no-contig-ec --no-scaffold-ec --file-type BAM`
- **Input:**
  - Zoysia_macrostachya: `genome.nextpolish.fasta`(Step2 NextPolish 산출물 — 다만 config.json의 `ref` 경로가 Ratatosk 계보/비-Ratatosk 계보 중 어느 쪽인지 파일명만으로 확정 불가, 확인 필요) + Omni-C 원시 리드(`Zma_S222_L002_R1/2_001.fastq.gz`)
  - Zoysia_sinica: `Z.sinica_IO_genome.original.fasta`(Step2 contig 산출물, 확인 필요) + Omni-C paired fastq(파일 자체는 미확인, 로그로만 존재) — 참고로 Hi-C 매핑과 동시에 참조 없는 자체 ragtag pre-scaffold(`Zs_ragtag.scaffold.fasta`)도 별도로 병행 생성되어 yahs/SALSA2/3D-DNA 각각에 "original"과 "ragtag" 두 계열 input으로 함께 투입됨
  - Rubus_Pangenome_Project(Rco/Rcra/Rpho): Step 2 hifiasm 산출 hap1/hap2 p_ctg fasta(purge3) + Omni-C paired fastq
- **Output:**
  - Zoysia_macrostachya: `purged.fa`/`hap.fa`/`dups.bed`(Purge_dups) → `mapped.PT.bam`/`HiC.filtered.bam`(정렬·필터링) → `scaffolds.raw.agp`(HapHiC) → `out_JBAT.FINAL.fa`(1차 Juicebox 수동 curation) → `out_JBAT_SV.FINAL.20chr.fa`(2차 Juicebox 수동 curation, 20개 염색체 확정)
  - Zoysia_sinica: `Zs.genome_mapped.PT.bam`(Omni-C 매핑) → `yahs.out_scaffolds_final.fa`(yahs, 최종 채택) / SALSA `scaffolds/`(SALSA2) / `3D-DNA_scaffold.agp`(3D-DNA+Juicer) → (QUAST 비교로 yahs 계열 확정, Step 7) → `Zs_superscaffolds.nogaps.fasta` 계열 → RagTag 재정렬로 `Z.sinica_IO_genome.original.FINAL.superscaffolds.nogaps_rename_ordered_modified.fasta`(chr01~chr20 명명·역상보 반영, 1차 RagTag 재정렬 완료본)
  - Rubus_Pangenome_Project(Rco/Rcra/Rpho): `<hap>_HiC.filtered.bam`(정렬·필터링) → `<hap>_no_ec_scaffolds_final.fa/.agp`(YaHS)
- **실행 커맨드:**
  - Zoysia_macrostachya (Purge_dups는 config.json만 존재, 실행 스크립트 자체는 확인 필요 — 스크립트 파일 없음)
  - Zoysia_macrostachya (Omni-C 리드 검증, Validating_Omni-C_reads.sh):
    ```bash
    bwa mem -5SP -T0 -t 150 purged.fa Zma_S222_L002_R1_001.fastq.gz Zma_S222_L002_R2_001.fastq.gz -o aligned.sam
    pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 150 --nproc-out 150 --chroms-path Zma_purged.genome aligned.sam > parsed.pairsam
    pairtools sort --nproc 150 parsed.pairsam > sorted.pairsam
    pairtools dedup --nproc-in 150 --nproc-out 150 --mark-dups --output-stats dup_stats.txt --output dedup.pairsam sorted.pairsam
    pairtools split --nproc-in 150 --nproc-out 150 --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam
    samtools sort -@150 -o mapped.PT.bam unsorted.bam && samtools index mapped.PT.bam
    ```
  - Zoysia_macrostachya (HapHiC 스캐폴딩, HapHiC/HapHiC.sh — 전처리(bwa+samblaster+filter_bam) 커맨드는 스크립트 내 주석 처리되어 비활성 상태이며 `HiC.filtered.bam`은 사전 생성된 파일을 그대로 재사용):
    ```bash
    haphic pipeline purged.fa HiC.filtered.bam 20 --correct_nrounds 2 --threads 150 --processes 150 --max_inflation 5
    ```
  - Zoysia_macrostachya (Juicebox 수동 curation): juicer pre/juicer_tools로 `.hic` 파일 생성 후 외부 JuiceBox GUI에서 수동 검토·편집 — GUI 조작 자체는 스크립트화되어 있지 않음(확인 필요)
  - Zoysia_sinica (Omni-C 매핑, 원본 `.sh` 미보존 — 로그와 이후 동일 로직 스크립트(`Contact_map/Generating_Contactmap.sh`)로 재구성한 것으로 완전한 원문은 아님 — 확인 필요):
    ```bash
    bwa mem -5SP -T0 -tN <ref.fasta> R1.fastq R2.fastq \
      | pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 ... \
      | pairtools sort ... \
      | pairtools dedup --mark-dups --output-stats stats.txt \
      | pairtools split --output-pairs mapped.pairs --output-sam - \
      | samtools view -bS | samtools sort -o mapped.PT.bam
    samtools index mapped.PT.bam
    ```
  - Zoysia_sinica (Hi-C 스캐폴딩 3종 병행, 원본 실행 스크립트 미보존 — 디렉토리 구조와 로그로만 확인, 확인 필요):
    - yahs(bin) → `1_Omni-C/yahs/{original,ragtag}/Zs_yahs.log`
    - SALSA2: `RE_sites.py`→`make_links.py`→`fast_scaled_scores.py`→`layout_unitigs.py`→`break_contigs`→`refactor_breaks.py`(iterative, `SALSA_Zs.log`)
    - `juicer.sh`(3D-DNA 파이프라인용) → `1_Omni-C/juicer/{contigs,ragtag}/`
  - Zoysia_sinica (yahs 계열에 대한 1차 RagTag 재정렬/명명, `1_Omni-C/juicer/contigs/aligned/RagTag/Zs_superscaffold2Zj/ordering.py` + `reverse_complement_chromosomes.py`; ragtag 자체 실행 스크립트는 미보존):
    ```bash
    # ragtag scaffold(참조: Z. japonica) 실행 → ragtag.scaffold.fasta (커맨드라인 원문 미보존, 확인 필요)
    python ordering.py            # chr01~chr20 순서 부여
    python reverse_complement_chromosomes.py   # 13개 염색체 역상보(Bio.Seq.reverse_complement)
    ```
  - Rubus_Pangenome_Project(Rco/Rcra/Rpho, `<hap>`=hap1/hap2 각각 반복):
    ```bash
    bwa mem -5SP -t 100 <hap>.p_ctg.fa <species>_Omni-C_R1.fastq.gz <species>_Omni-C_R2.fastq.gz \
        | samblaster | samtools view -F 3340 -o <hap>_HiC.bam
    /usr/local/bin/HapHiC/utils/filter_bam <hap>_HiC.bam 1 --nm 3 --threads 100 \
        | samtools view -b -o <hap>_HiC.filtered.bam
    yahs -o <hap>_no_ec --telo-motif AAACCCT --no-contig-ec --no-scaffold-ec --file-type BAM \
        <hap>.p_ctg.fa <hap>_HiC.filtered.bam
    ```
    Rta(AutoHiC): 실행 커맨드는 2절 Rubus_Pangenome_Project 항목 참조.

### Step 6. Gap closing(YagCloser) — Zoysia_macrostachya만 해당
- **소프트웨어(버전):** minimap2 + samtools(정렬), detgaps, yagcloser.py(`/usr/local/bin/yagcloser/`, 버전 확인 필요)
- **핵심 파라미터:** minimap2 `-x map-ont -t200 --secondary=yes`(스크립트 내 주석 처리 — 실제 실행 시에는 사전 생성된 `aln.s.bam`을 재사용)
- **Input:** `out_JBAT_SV.FINAL.20chr.fa`(Step5 산출물, 20염색체 확정본) + `ONT_corrected.fastq`(원시 ONT 보정 리드)
- **Output:** `out_JBAT_SV.FINAL.20chr.nogaps.fasta`(갭 클로징 완료) → (염색체 정렬/개명 처리, 실행 스크립트 확인 필요) → `Zma_FINAL.20chr.nogaps.sorted.fasta`(이후 모든 하류 스텝의 최종 게놈)
- **실행 커맨드:**
  - Zoysia_macrostachya (03.YagCloser/After_manualcuration/yagcloser.sh):
    ```bash
    samtools index aln.s.bam
    detgaps out_JBAT_SV.FINAL.20chr.fa > out_JBAT_SV.FINAL.20chr.gaps.bed
    yagcloser.py -g out_JBAT_SV.FINAL.20chr.fa -b out_JBAT_SV.FINAL.20chr.gaps.bed -a aln.s.bam -s Zma_superscaffolds -o yagcloser_outdir
    python /usr/local/bin/yagcloser/scripts/update_assembly_edits_and_breaks.py -i out_JBAT_SV.FINAL.20chr.fa -o out_JBAT_SV.FINAL.20chr.nogaps.fasta -e Zma_scaffolds.edits.txt
    ```
  - 이후 염색체 정렬·이름 확정(`chr.id` 매핑 기반으로 추정)하여 `Zma_FINAL.20chr.nogaps.sorted.fasta` 생성 — 해당 정렬 스크립트 자체는 확인 필요 — 스크립트 파일 없음

### Step 7. Genome QC
- **소프트웨어(버전):** BUSCO v5.8.2(Cannabis/Pichia 공통), QUAST v5.2.0(Cannabis/Pichia 공통); Abies: QUAST 5.3.0(fb88221c)만 사용(BUSCO 미수행 — 확인 필요, 프로젝트 초기 단계로 추정); Zoysia_macrostachya: BUSCO(embryophyta/liliopsida/poales_odb10) + QUAST 5.3.0 + LAI(LTR_FINDER_parallel+ltrharvest+LTR_retriever) + Merqury(meryl k=21); Zoysia_sinica: QUAST(버전 확인 필요, 스캐폴더 3종 비교 용도로도 사용) + BUSCO(genome mode, poales_odb10, 버전 확인 필요) + Merqury(meryl+merqury.sh, 버전 확인 필요) — LAI는 사용 정황 없음
- **핵심 파라미터:** BUSCO: `-m genome`, lineage(Cannabis: rosales_odb12/rosids_odb10, Pichia: ascomycota_odb12/saccharomycetes_odb12, Zoysia_macrostachya: embryophyta_odb10/liliopsida_odb10/poales_odb10, Zoysia_sinica: poales_odb10); QUAST: 참조 유전체 지정, 기본 옵션(Abies는 `--large --threads 64 --labels ...` 또는 `--labels JSlink_only,JSlink_EYEONCELL`; Zoysia_sinica는 정확한 커맨드라인 확인 필요, report만 확인); Zoysia_macrostachya LAI: `LTR_retriever -threads 100`; Zoysia_macrostachya/Zoysia_sinica Merqury: `meryl k=21`(macrostachya) / `meryl k=19 memory=800GB threads=130`(sinica)
- **Input:** scaffold(chr 단위) fasta(Cannabis/Pichia/Zoysia_macrostachya); Abies: organelle 제거된 nuclear contig(hap1/hap2/p, scaffolding 이전 단계); Zoysia_sinica: `Zs_superscaffolds.nogaps.fasta`(yahs 계열)와 `Z.sinica_IO_genome.original.FINAL.fasta`(3D-DNA 계열)를 QUAST로 직접 비교, 이후 BUSCO/Merqury는 채택된 yahs 계열(`Zs_superscaffolds.nogaps.fasta`)만 사용
- **Output:** BUSCO short_summary(.txt/.json), QUAST report(.html/.tsv/.pdf), Zoysia_macrostachya는 추가로 LAI score(`*.out.LAI`), Merqury QV/completeness; Zoysia_sinica는 QUAST report에서 yahs 계열(20 contigs, N50 18.8Mb)이 3D-DNA 계열(124 contigs, 동일 N50이나 더 분절적)보다 우수하다는 근거로 스캐폴더를 확정(다만 "왜 yahs를 최종 채택했는지"에 대한 명시적 의사결정 문서는 없고, 이후 모든 하위 스텝이 yahs 계열 파일명을 사용한다는 정황 증거만 있음 — 확인 필요), Merqury QV/completeness.stats
- **실행 커맨드:**
  - Cannabis (`<sample>_<ref>` 조합마다 동일 패턴 반복, quast.log에서 확인):
    ```bash
    busco -m genome -i <sample>.fasta -o <sample>.rosales_odb12 -c 150 --datasets_version odb10 -f -l rosales_odb12
    quast.py -o <sample>_<ref> -r <ref>.softmasked.fasta -t 150 <sample>_<ref>_RagTag.chr.fasta
    ```
  - Pichia (`<lineage>`=ascomycota_odb12/saccharomycetes_odb12, quast.log에서 확인):
    ```bash
    busco -i <sample>_NextDP_NextPsh.fasta -l /pbi-acc1/hsPark/Pichia_Pastoris_dnDNAseq_HanHwa/BUSCO/X-33/busco_downloads/lineages/<lineage> -o <sample>_NextDP_NextPsh_<lineage>_busco --cpu 200 --mode genome
    quast.py --threads 100 --output-dir nextDenovo <sample>_nextND_nextPsh.fasta
    ```
  - Abies (JSlink 단독 vs JSlink+EYEONCELL 통합 비교, hifiasm_JSLink_EYEONCELL/03_quast_compare.sh, `<hap>`=p/hap1/hap2):
    ```bash
    quast.py ../hifiasm/abies_koreana.bp.<hap>.p_ctg.nuclear.fa abies_koreana_combined.bp.<hap>.p_ctg.nuclear.fa \
        --threads 32 --labels "JSlink_only,JSlink_EYEONCELL" --output-dir quast_compare_<hap>
    ```
    (3-way self 비교인 `quast_nuclear/`는 quast.log에서 역추출된 커맨드: `quast.py abies_koreana.bp.p_ctg.nuclear.fa abies_koreana.bp.hap1.p_ctg.nuclear.fa abies_koreana.bp.hap2.p_ctg.nuclear.fa --large --threads 64 --labels p_ctg,hap1,hap2 --output-dir quast_nuclear` — 이 커맨드를 생성한 원본 wrapper 스크립트 자체는 확인 필요)
  - Zoysia_macrostachya (1차 어셈블리 QC, 01.Assembly/QC_results/BUSCO/BUSCO.sh + Quast/Quast.sh):
    ```bash
    busco -m genome -i Ratatosk_NextDP_NextPsh_ntLinkk40w500.fasta -o Ratatosk_NextDP_NextPsh_ntLinkk40w500.embryophyta_odb10 -c 256 --datasets_version odb10 -f -l embryophyta_odb10
    quast.py Ratatosk_NextDP_NextPsh_ntLinkk40w500.fasta -o Ratatosk_NextDP_NextPsh_ntLinkk40w500 -t 256
    ```
  - Zoysia_macrostachya (최종 게놈 LAI/Merqury, 05.Genome_QC/LAI/.../LAI.sh + Mercury/.../mercury.sh):
    ```bash
    LTR_retriever -genome Zma_chromsome_nogap.sorted.fixed.fa -inharvest Zma_chromsome_nogap.sorted.fa.rawLTR.scn -threads 100
    merqury.sh Zma.meryl Zma_FINAL.20chr.nogaps.sorted.fasta Zma_chromosome
    ```
  - Zoysia_sinica (QUAST 스캐폴더 비교, 원본 QUAST 실행 스크립트 미보존 — `quast.log`만 확인, 확인 필요):
    ```bash
    # quast.py <옵션 미확인> Zs_superscaffolds.nogaps.fasta Z.sinica_IO_genome.original.FINAL.fasta -o 2_Genome_QC/quast/...
    ```
  - Zoysia_sinica (BUSCO genome mode, `BUSCO/busco.sh`; 스크립트 자체는 protein 모드 for-loop이고 genome 모드 실행 스크립트 원본은 미보존):
    ```bash
    busco -i Zs_superscaffolds.nogaps.fasta --lineage_dataset poales_odb10 -o Zs_poales_odb10 -m genome --cpu 200   # 확인 필요: 게놈모드 실제 파라미터는 재구성값
    ```
  - Zoysia_sinica (Merqury, `2_Genome_QC/merqury/merqury_QC.sh`):
    ```bash
    meryl k=19 memory=800GB threads=130 count *.fq.gz output Zs.meryl
    merqury.sh Zs.meryl Zs_superscaffolds.nogaps.fasta Zs_superscaffolds.nogaps
    ```

### Step 8. 비교유전체 분석(synteny/dot-plot/pan-genome)
- **소프트웨어(버전):** Cannabis: chromeister + minimap2/syri/plotsr(버전 확인 필요) / Pichia: mummer 4.0.1(nucmer+mummerplot) / Zoysia_macrostachya: minimap2+samtools+SyRI+plotsr(구조변이), RagTag(타 Zoysia속 게놈 비교용 재정렬), OrthoFinder(오솔로그 클러스터링, 버전 확인 필요), CHROMEISTER/MCScan/CandiSSR/CNV_analysis/Purge_dups/Mercury/BUSCO 등(Pan-genome_analysis 하위 병렬 서브모듈군). Zoysia_sinica도 매우 유사한 downstream 비교유전체 스위트(OrthoFinder/CAFE/MCMCTree/Ks-MCScanX/PAML PSG 등)를 보유하지만, macrostachya처럼 하나의 통합 `Pan-genome_analysis/` 디렉토리가 아니라 `5_`~`15_` 번호가 매겨진 개별 폴더로 분리되어 있다 — 상세는 2절 Zoysia_sinica 항목 참조
- **핵심 파라미터:** Cannabis: minimap2 `-ax asm5 --eqx`, syri `-F B` / Pichia: nucmer `--maxmatch -l 20 -g 500 -c 100` / Zoysia_macrostachya: minimap2 `-ax asm5 -t 200 --eqx`, syri `-F B --nc 200`, RagTag `--mm2-params '-x asm5' -w -t 150`
- **Input:** scaffold(chr 단위) fasta 및 비교 대상 유전체(Cannabis: 판게놈 98샘플; Pichia: 공개 Pichia 균주 게놈; Zoysia_macrostachya: `Zma_FINAL.20chr.nogaps.sorted.fasta` + 타 Zoysia속 게놈 ZJN/ZMW/ZPZ/Z.sinica 등, 원본 출처는 확인 필요)
- **Output:** dot-plot 이미지(.png), 구조변이 vcf/summary(Cannabis, Zoysia_macrostachya); Zoysia_macrostachya는 추가로 RagTag scaffold, OrthoFinder Orthogroups.tsv 등(Step 9의 MCMCTree/CAFE5/PSG 입력으로 재사용, 2절 참조)
- **실행 커맨드:**
  - Cannabis (sourmash.sh, plotsr.sh — 다른 genome 쌍은 스크립트 내 대부분 주석 처리되어 미실행 정황; chromeister는 확인 필요 — 스크립트 파일 없음):
    ```bash
    sourmash sketch dna -p scaled=1000,k=31 <genome>.fasta.gz --output <genome>.sig
    minimap2 -ax asm5 -t 256 --eqx CS_Male_AH3Ma_hap2_RagTag.chr.fasta CS_Male_AH3Mb_Y_RagTag.chr.modified.sorted.fasta | samtools sort -O BAM - > D_E.bam
    syri -c D_E.bam -r CS_Male_AH3Ma_hap2_RagTag.chr.fasta -q CS_Male_AH3Mb_Y_RagTag.chr.modified.sorted.fasta -F B --prefix D_E
    ```
  - Pichia (nucmer_plot.sh, `<ref>`은 대상 공개 균주별 .fna 파일마다 반복):
    ```bash
    nucmer --maxmatch -l 20 -g 500 -c 100 -p <ref>_vs_query -t 16 <ref>.fna X-33_flye_GS115_scaffolded.fasta
    mummerplot -t png -p <ref>_vs_query --layout -R <ref>.fna -Q X-33_flye_GS115_scaffolded.fasta -s large --color <ref>_vs_query.delta
    ```
  - Zoysia_macrostachya (구조변이, Pan-genome_analysis/syri/syri.sh):
    ```bash
    minimap2 -ax asm5 -t 200 --eqx Z.sinica_manualcurated.FINAL_merged2.fasta Zma_FINAL.20chr.nogaps.sorted.fasta | samtools sort -O BAM - > curated_Zs2_Zma.bam
    samtools index curated_Zs2_Zma.bam
    syri -c curated_Zs2_Zma.bam -r Z.sinica_manualcurated.FINAL_merged2.fasta -q Zma_FINAL.20chr.nogaps.sorted.fasta -F B --prefix curated_Zs2_Zma --nc 200
    plotsr --sr curated_Zs2_Zmasyri.out --sr Zma_Zjsyri.out --genomes genome.txt -o curated_Zs2-Zma-Zj_plot.png
    ```
    (참고: 이 비교의 참조 서열로 쓰인 `Z.sinica_manualcurated.FINAL_merged2.fasta`는 Zoysia_sinica 프로젝트의 Step 10 수동 큐레이션 최종 산출물과 동일 파일명이다 — 두 프로젝트 문서가 같은 파일을 서로 다른 역할로 참조하는 지점)
  - Zoysia_macrostachya (타 Zoysia속 게놈 비교용 재정렬, Pan-genome_analysis/RagTag/RagTag.sh):
    ```bash
    ragtag.py scaffold Zsg_Zma_ZJN_chr_genome.fa <타종>_r1.0.genome.fasta -o <타종>_ZsgZmaZJN_scaffold_asm5 -w -t 150 --mm2-params '-x asm5'
    ```
  - Zoysia_macrostachya (OrthoFinder/CHROMEISTER/MCScan/CandiSSR/CNV/Purge_dups/Mercury/BUSCO 등): 실행 스크립트가 없고 로그 파일(`orthofinder.log` 등)만 남아있는 경우가 많음 — 정확한 실행 커맨드는 확인 필요(2절 참조)

### Step 9. 유전자/반복서열 주석
- **소프트웨어(버전):** Cannabis: RepeatModeler 2.0.6 + RepeatMasker 4.1.9 → BRAKER3 3.0.8(de novo 유전자 예측, RNA-seq+protein evidence 사용) / Pichia: LiftOff(참조 기반 주석 liftover, GS115 및 CBS7435 기준) / Zoysia_macrostachya: RepeatModeler+RepeatMasker(버전 확인 필요, Cannabis와 동일 계열 툴) → Trimmomatic+prinseq+Hisat2(RNA-seq 전처리) → BRAKER3(실행 스크립트/컨테이너 이미지 부재, 메모 파일만 존재) + GeneMarkS-T(ISOseq 보정) + TSEBRA(통합) / Zoysia_sinica: RepeatModeler+RepeatMasker(최종 스캐폴드 기준, 버전 확인 필요) → BRAKER **v3.0.7**(ETP mode — "braker.pl version 3.0.7" 로그에서 직접 확인됨, AUGUSTUS/GeneMark-ETP/TSEBRA/ProtHint/DIAMOND 동봉) → AGAT **v1.4.0**(rename + longest-isoform 후처리, conda env `agat`) → BUSCO(protein mode, poales_odb10)
- **핵심 파라미터:** Cannabis: RepeatModeler `-LTRStruct`, BRAKER3 `--threads 125`; Pichia: LiftOff `-p 150`(CBS7435 대상은 `-exclude_partial -copies` 추가); Zoysia_macrostachya: RepeatModeler `-LTRStruct -threads 100`, RepeatMasker `-xsmall -gff -pa 100`, Trimmomatic `ILLUMINACLIP:...:2:30:10:2:keepBothReads LEADING:10 TRAILING:10 MINLEN:50 -threads 100`, Hisat2 `-p 200`, TSEBRA `long_reads.cfg`; Zoysia_sinica: RepeatModeler `-LTRStruct -threads 100`, RepeatMasker `-xsmall -gff -pa 100`; BRAKER3 `--gff3 --threads 8`(ETP mode: RNA-seq BAM + protein DB 동시 입력); AGAT rename/longest-isoform 정확한 서브커맨드는 확인 필요; BUSCO(protein) `-m proteins --cpu 200`
- **Input:** Cannabis: scaffolding 이전 contig-level masked genome + mRNA-seq + 판게놈 단백질 세트 / Pichia: scaffold(chr 단위, rename 前) fasta + 참조 gff / Zoysia_macrostachya: `Zma_FINAL.20chr.nogaps.sorted.fasta`(Step6 최종본) + Zma-R1/S1 mRNA-seq + ISOseq reads / Zoysia_sinica: `Zs_superscaffolds.nogaps.fasta`(Step5 산출 최종 스캐폴드, RagTag 재정렬 前 계열 — RepeatMasker/BRAKER 입력은 rename 前 파일명 사용) → RepeatMasker 출력 `.masked` + `Zs_superscaffolds.nogaps.hisat2.sorted.bam`(RNA-seq 정렬, hisat2 사용 추정 — 정렬 스크립트 원본 미보존) + `liliopsida_27sp.ncbi_pep.faa`(근연 27종 NCBI 단백질 DB)
- **Output:** Cannabis: braker.gtf/aa/codingseq(신규 유전자 모델) / Pichia: liftover gff + unmapped_features 목록(참조 유전자 좌표 이전) / Zoysia_macrostachya: `Zma_FINAL.20chr.nogaps.sorted.fasta.masked` → `braker.gtf` → `gmst.global.gtf`(GeneMarkS-T) → `tsebra.gtf`(TSEBRA 통합) → `Zma_braker3_longest.pep`/`Zma_braker3.rename.longest_isoform.gtf(3)`(최종 유전자 모델) / Zoysia_sinica: `Zs_superscaffolds.nogaps.fasta.masked`(RepeatMasker) → `braker.gff3/gtf/aa/codingseq` → AGAT rename/longest-isoform 후 `braker_renamed_longest_isoform.pep/cds` → BUSCO(protein) short_summary(`Zs_braker_poales_odb10/...txt`). 참고: 이 최종 스캐폴드 기준 RepeatModeler/Masker 이전에, Hi-C 매핑 전 contig 단계(Step 2 산출물)에서도 한 차례(2023-11경, RM 디렉토리 `RM_20.FriNov170056382023`) RepeatModeler/RepeatMasker가 수행된 정황이 있으나 그 산출물 실물은 남아있지 않고 스크립트만 확인됨(pre-Hi-C 1차 라운드, 확인 필요·추정 근거는 파일명/디렉토리 타임스탬프뿐)
- **실행 커맨드:**
  - Cannabis (module_RepeatModeler_Masker.sh, `<sample>_<hap>`=Female/Male × hap1/hap2 4종 동일 패턴; braker.log의 `BRAKER CALL:` 라인에서 확인):
    ```bash
    BuildDatabase -name CS_<sample>_<hap>.db Csativa_nuclear_<sample>.asm.bp.<hap>.p_ctg.fa
    RepeatModeler -database CS_<sample>_<hap>.db -LTRStruct Csativa_nuclear_<sample>.asm.bp.<hap>.p_ctg.fa -threads 200
    RepeatMasker -xsmall -gff -pa 200 -lib RM_*/consensi.fa.classified Csativa_nuclear_<sample>.asm.bp.<hap>.p_ctg.fa
    braker.pl --species=Cannabis_sativa --genome=Csativa_nuclear_<sample>.asm.bp.<hap>.p_ctg.fa.masked --rnaseq_sets_ids Csativa_female,Csativa_male --rnaseq_sets_dirs=. --threads 125 --prot_seq=Cannabis_98pangenome.v2.proteins.fasta
    ```
  - Pichia (Liftoff.sh; CBS7435 기준은 GS115 기준 실행 약 6주 후 별도 수행):
    ```bash
    liftoff -g GCF_000027005.1_ASM2700v1_genomic.gff -o <sample>_NextDP_NextPsh_scaffolded.gff -f feature_types.txt -dir itermediate_<sample>_NextDP_NextPsh_scaffolded_files -u <sample>_NextDP_NextPsh_scaffolded_unmapped_features.txt -p 150 <sample>_NextDP_NextPsh_scaffolded.fasta GCF_000027005.1_ASM2700v1_genomic.fna
    liftoff -g GCA_900235035.2_KP_7435-4_genomic.gff -o <sample>_NextDP_NextPsh_scaffolded_CBS7435.gff -f feature_types.txt -dir itermediate_<sample>_NextDP_NextPsh_scaffolded_files -u <sample>_NextDP_NextPsh_scaffolded_CBS7435_unmapped_features.txt -p 150 -exclude_partial -copies <sample>_NextDP_NextPsh_scaffolded.fasta GCA_900235035.2_KP_7435-4_genomic.fna
    ```
  - Zoysia_macrostachya (반복서열, 06.Repeat_element_annotation/module_RepeatModeler_Masker.sh):
    ```bash
    BuildDatabase -name Zma.db Zma_FINAL.20chr.nogaps.sorted.fasta
    RepeatModeler -database Zma.db -LTRStruct Zma_FINAL.20chr.nogaps.sorted.fasta -threads 100
    RepeatMasker -xsmall -gff -pa 100 -lib RM_*/consensi.fa.classified Zma_FINAL.20chr.nogaps.sorted.fasta
    ```
  - Zoysia_macrostachya (RNA-seq 전처리, 07.Gene_annotation/Trimmomatic/trimmomatic.sh + Hisat2/hisat2.sh):
    ```bash
    TrimmomaticPE -threads 100 -summary <sample>.stat <sample>_1.fastq <sample>_2.fastq <sample>_paired_1.fq <sample>_unpaired_1.fq <sample>_paired_2.fq <sample>_unpaired_2.fq ILLUMINACLIP:Trimmomatic_adapters-PE.fa:2:30:10:2:keepBothReads LEADING:10 TRAILING:10 MINLEN:50
    hisat2 -x Zma_FINAL.20chr.nogaps.sorted.fasta.masked.genome -p 200 -1 <sample>_good_paired_1.fq -2 <sample>_good_paired_2.fq -S <sample>.hisat2.sam
    ```
  - Zoysia_macrostachya (BRAKER3/GeneMarkS-T/TSEBRA): 확인 필요 — BRAKER3는 실행 스크립트/컨테이너 이미지 자체가 레포에 없고(`Braker_pipeline_HS.txt` 메모만 존재), GeneMarkS-T/TSEBRA 세부 커맨드도 본 문서 작성 시점에는 확인되지 않음(2절 참조)
  - Zoysia_sinica (반복서열, 최종 스캐폴드 기준, `3_Repeat_element_annoation/module_RepeatModeler_Masker.sh`, `.bash_history`에서 `nohup bash module_RepeatModeler_Masker.sh &`로 백그라운드 실행 확인):
    ```bash
    BuildDatabase -name Zs_superscaffold.db Zs_superscaffolds.nogaps.fasta
    RepeatModeler -database Zs_superscaffold.db -LTRStruct Zs_superscaffolds.nogaps.fasta -threads 100
    RepeatMasker -xsmall -gff -pa 100 -lib RM_*/consensi.fa.classified Zs_superscaffolds.nogaps.fasta
    ```
    (참고: 이보다 앞서 Hi-C 매핑 이전 contig 단계(`Z.sinica_IO_genome.original.fasta`)에 대해서도 동일 툴 체인이 한 차례(2023-11경) 수행된 정황이 있음 — 위 소프트웨어(버전) 항목 참고 노트 참조)
  - Zoysia_sinica (BRAKER3 ETP mode, 로그(`braker.log`/`BRAKER.log`)에서 확인, 원본 실행 스크립트 자체는 4_Gene_annotation에 미보존):
    ```bash
    /opt/BRAKER/scripts/braker.pl --genome=Zs_superscaffolds.nogaps.fasta.masked \
        --bam=Zs_superscaffolds.nogaps.hisat2.sorted.bam \
        --prot_seq=liliopsida_27sp.ncbi_pep.faa --gff3 --threads 8 \
        --workingdir=/pbi-acc1/hsPark/Zoysia_sinica_dnDNAseq_revised/BRAKER
    ```
  - Zoysia_sinica (AGAT rename/longest-isoform 후처리, AGAT v1.4.0, conda env `agat`; 원본 rename 스크립트 미보존, AGAT 로그로만 확인 — 정확한 서브커맨드는 확인 필요):
    ```bash
    # AGAT 기반 gtf/gff rename 및 longest-isoform 추출 (braker.gtf → braker_renamed.gtf → braker_renamed_longest_isoform.gtf)
    ```
  - Zoysia_sinica (BUSCO protein mode, `Gene_busco_poales_odb10.log`; top-level `BUSCO/busco.sh`와 동일 패턴으로 추정):
    ```bash
    busco -i braker_renamed_longest_isoform.pep --lineage_dataset poales_odb10 -m proteins --cpu 200
    ```

### Step 10. 수동 큐레이션(Manual_curation) — Zoysia_sinica만 해당
Step 9의 BRAKER3+AGAT 산출물을 사람이 직접 검토·수정하는 단계로, 다른 4개 프로젝트에는 없는
Zoysia_sinica만의 스텝이다. 이 스텝의 산출물(최종 병합 게놈)은 이후 한 번 더 RagTag로 Z. japonica
기준 재정렬된다(Step 4/Step 5 상단 참고 블록에서 언급한 "2차" RagTag 실행).
- **소프트웨어(버전):** AGAT v1.4.0(conda env `agat`, 로그 `Zsg_braker_rename.agat.log`에서 확인) 기반 gff3 rename → 이후 수기 큐레이션(자동화 스크립트 확인 안 됨, 사람이 직접 gff3를 편집한 것으로 추정) → RagTag(버전 확인 필요, 최종 병합 게놈에 대한 2차 재정렬)
- **핵심 파라미터:** 확인 필요 — Manual_curation 폴더 자체에는 실행 스크립트가 없어(AGAT 로그만 존재) 수동 개입 위주의 스텝이라 파라미터화가 어려움. 2차 RagTag의 실행 파라미터도 미보존(`.err` 로그만 확인)
- **Input:** `braker_renamed_longest_isoform.gff`류(Step 9 산출 계열, 다만 파일명이 `Zsg_braker_rename.*`로 약간 다르게 등장해 정확히 어느 산출물이 직접 입력되었는지는 확인 필요)
- **Output:** `Zsg_manualcurated.gff3`(수동 검토 반영) + `Zsg_unmapped_features.txt`(수동 반영 안 된 피처 45개 이상 목록) + `Z.sinica_manualcurated.FINAL_merged2.fasta`(최종 병합 게놈) → (2차 RagTag 재정렬) → `Manualcurated2Zj/ragtag.scaffold.fasta`(Z. japonica 기준 최종 재정렬본)
- **실행 커맨드:**
  - Zoysia_sinica (AGAT rename, `Zsg_braker_rename.agat.log`에서 버전 확인: "AGAT - Version: v1.4.0"):
    ```bash
    # AGAT 기반 gff3 rename (정확한 서브커맨드는 확인 필요) → Zsg_braker_rename.gff3/gtf
    ```
  - Zoysia_sinica (수동 큐레이션): Manual_curation 폴더 자체에는 실행 스크립트가 없음 — AGAT 로그 외 자동화 스크립트가 확인되지 않아, 사람이 직접 gff3를 편집해 `Zsg_manualcurated.gff3`와 `Z.sinica_manualcurated.FINAL_merged2.fasta`를 만들었을 가능성이 높음(전체 산출물이 2024-12-18 11:30~11:31 한 세션에서 일괄 생성됨) — 정확한 병합 로직은 확인 필요, 추정하지 않음
  - Zoysia_sinica (2차 RagTag 재정렬, `1_Omni-C/juicer/contigs/aligned/RagTag/Manualcurated2Zj/`, ragtag 실행 스크립트 원문 미보존 — `.err` 로그만 확인):
    ```bash
    # ragtag scaffold(참조: Z. japonica pseudomolecule ZJN_r1.1.pseudomol) 대상 Z.sinica_manualcurated.FINAL_merged2.fasta 재정렬 (커맨드라인 원문 확인 필요)
    ```
  - 불확실 사항: 이 2차 RagTag 재정렬 결과가 최종 출판/공개본으로 이어졌는지는 이후 단계가 보존되어 있지 않아 확인 불가. top-level `Zsinica_genome.fasta`(2024-08-14, 321,023,349 bytes)는 같은 시기의 `Z.sinica_IO_genome.original.FINAL.superscaffolds.nogaps_rename_ordered.fasta`(2024-08-07, 317,883,127 bytes)와 파일 크기가 달라 동일 파일이 아니며, 어느 스텝의 산출물을 가공한 것인지 스크립트 근거로 확인하지 못했다 — 추정하지 않고 확인 필요로 남긴다.

## 2. 프로젝트별 특이사항

### Cannabis_sativa_cv.ChungSam_dnDNAseq
- 연구목적/샘플 설계: 한국 재래종 대마(청삼, Chungsam) 자웅이주(암/수) 개체에 대해 PacBio Revio HiFi 리드 기반 de novo 반수체 분리(haplotype-resolved) 전장유전체 조립, 판게놈 비교, 카나비노이드 합성 유전자(THCAS/CBDAS/CBGAS) 분석을 수행.
- 공통 파이프라인과의 차이점:
  - 시퀀싱 플랫폼은 PacBio Revio HiFi 단일 long-read 플랫폼만 사용(하이브리드 short-read 보정 없음). mRNA-seq(Illumina)는 조립이 아니라 BRAKER3 유전자예측의 hint 데이터로만 사용됨.
  - hifiasm 조립이 프로젝트 내에서 두 번 수행됨: 1차(2025-05-18, organelle 포함 전체 reads, 이후 다운스트림에서 재사용되지 않는 폐기 경로로 판단되나 명시적 폐기 근거는 없어 확인 필요) → Oatk(v1.0)로 엽록체/미토콘드리아 게놈 분리 후 organelle-비매핑(nuclear) reads만 추출 → 2차(최종 채택) hifiasm 조립. 이 organelle 분리 단계는 Pichia 프로젝트에는 없음.
  - Female/Male 각 개체를 hap1/hap2로 반수체 분리(haplotype-resolved diploid assembly)하여 총 4개 haplotype 조립본을 생성 — Pichia는 haploid 균주로 haplotype 분리 개념이 없음.
  - RagTag scaffolding 기준 참조는 Sourmash(v4.9.4, k=31 scaled=1000) k-mer 유사도 비교로 판게놈 98개 haplotype 중에서 선정된 AH3Ma(Female 계열)/AH3Mb(Male-Y 계열) — NCBI 참조(GCF_029168945.1)도 별도로 보유하고 있으나 실제 RagTag에 어느 참조가 쓰였는지는 두 참조 체계가 혼용되는 정황이 있어 확인 필요.
  - RepeatModeler/RepeatMasker + BRAKER3(de novo 유전자 예측)는 scaffolding된 chromosome-level fasta가 아니라 scaffolding 이전의 nuclear contig-level fasta를 입력으로 사용한 반면, BUSCO/QUAST/plotsr/chromeister는 RagTag scaffold 결과를 chromosome 단위로 재정리한 파일을 입력으로 사용 — RagTag 이후 "scaffolded 계열"과 "contig 계열"로 다운스트림이 분기됨.
  - Genome survey(FastK+GenomeScope2)의 정확한 도구 버전은 로그에 기록되어 있지 않아 확인 필요. Extract_unmapped_reads.sh(minimap2+samtools) 단계의 minimap2/samtools 버전도 확인 필요.
  - chromosome 단위 fasta(`*_RagTag.chr.fasta`) 생성 스크립트, CHROMEISTER 실행 커맨드/버전, syri/plotsr 버전은 모두 로그/스크립트가 남아있지 않아 확인 필요.
  - BRAKER3 실행이 컨테이너 기반으로 추정되나 정확한 이미지명은 확인 필요.
- (해당 시) 프로젝트 고유 확장 분석:
  - 판게놈(98개 haplotype) 비교: Sourmash 유사도 매트릭스/heatmap, chromeister dot-plot, plotsr+syri 구조변이 시각화(X/Y 성염색체 비교 포함).
  - THCAS/CBDAS/CBGAS 카나비노이드 합성 유전자 분석: (1) BRAKER3 예측 단백질과 THCAS 참조유전자(AB057805.1) 간 tblastn/blastp 비교(메인 조립 파이프라인으로부터 약 6개월 후 수행된 사후 검증, 정확한 BLAST+ 커맨드/버전 확인 필요), (2) 별도 독립적인 범용 targeted gene assembly 파이프라인(BWA-MEM 정렬 → SPAdes de novo assembly → MAFFT 정렬 → 변이 분석, conda env `gene_assembly`) — 이 파이프라인이 사용한 샘플(cBGambient/cBGLimonenon)이 청삼 Female/Male 샘플과 동일한지는 확인 필요.

### Pichia_Pastoris_dnDNAseq_HanHwa
- 연구목적/샘플 설계: 한화솔루션 회사의 Pichia 균주(BG-10, X-33) 2종에 대해 ONT(MinION) 시퀀싱 기반 de novo genome assembly 수행, 한화솔루션에 최종 전달 완료.
- 공통 파이프라인과의 차이점:
  - ONT(MinION) long-read와 Illumina short-read를 병행 사용하는 하이브리드 전략: Illumina reads는 Prinseq으로 QC/trim 후 (a) FastK+GenomeScope 입력, (b) Ratatosk를 이용한 ONT reads의 hybrid error-correction, (c) NextPolish 폴리싱 세 용도로 재사용됨.
  - De novo assembly는 trycycler 기반 multi-assembler(Flye/NextDenovo/miniasm 12개 후보) 합의 조립 경로가 시도되었으나, reconcile 단계에서 pairwise identity 에러(기준 미달)로 실패하여 폐기됨 — 최종적으로는 corrected long reads 전체를 NextDenovo로 raw assembly 후 NextPolish(Illumina 폴리싱)로 마무리하는 단일 경로가 채택됨(X-33은 Flye 전체-리드 조립도 비교 후보로 병행). 채택 사유는 BUSCO/Merqury/Quast QC 지표 중 contiguity가 가장 높았기 때문으로 확인됨.
  - Haploid 단일 genome 조립(Cannabis와 달리 haplotype 분리 개념 없음).
  - RagTag scaffolding 기준 참조는 NCBI GS115(GCF_000027005.1) 단일 참조(Cannabis처럼 유사도 비교를 통한 참조 선정 절차는 없음).
  - 유전자 주석은 BRAKER3류의 de novo 예측이 아니라 LiftOff를 이용한 참조 기반 annotation liftover — GS115 기준 liftover를 먼저 수행하고, 약 6주 후 CBS7435 기준 liftover(`-exclude_partial -copies` 옵션 추가)를 별도로 수행.
  - Scaffold 완료 후 contig 헤더를 rename(`_RagTag` → `_BG-10`/`_X-33`)하여 pseudochromosome fasta를 만드는 단계가 있으나, 이를 수행한 스크립트/커맨드 기록이 없어 확인 필요.
  - RagTag/Quast/Merqury/Hanwha 디렉토리에는 실행 스크립트가 보관되어 있지 않아 CLI 직접 실행으로 추정되며, RagTag의 정확한 버전은 확인 필요.
- (해당 시) 프로젝트 고유 확장 분석:
  - Merqury를 이용한 k-mer 기반 assembly QV/completeness 평가(Cannabis에는 이 QC 도구가 없음).
  - mummer(nucmer+mummerplot) 기반 공개 Pichia 균주(ATCC28485, Phaff 등) 대비 synteny dot-plot 비교.
  - Yagcloser를 이용한 scaffold gap 분석(X-33만 수행, 실제 gap-filling 커맨드는 비활성화되어 있어 분석 전용으로 그침).
  - Genome_editing: gRNA를 BLAST로 게놈에 매핑하고 LiftOff 주석과 교차(intersect)하여 유전체 편집 타겟 위치를 확인(BG-10만 수행).

### Abies_koreana_dnDNAseq
- 연구목적/샘플 설계: 구상나무(제주도 표준목) 전장 유전체 해독. JSLink/EYEONCELL 두 시퀀싱 회사가 각각 PacBio Revio HiFi 4 cell씩(총 8 cell) 생산했으며, 현재 "JSLink 4cell 단독" 조립본과 "JSLink+EYEONCELL 8cell 통합" 조립본 두 갈래를 병행하여 조립/QC를 진행 중.
- 공통 파이프라인과의 차이점:
  - 2026-07-13(문서 작성 기준일) 현재 아직 초기 단계 프로젝트로, Step1(genome survey)·Step2(assembly)·Step3(organelle 분리)·Step7(QUAST QC)까지만 완료되었고 scaffolding, 비교유전체, 유전자주석 단계는 아직 수행되지 않음 — Cannabis/Pichia/Zoysia_macrostachya/Zoysia_sinica와 달리 "완결된" 파이프라인이 아니라 "진행 중" 파이프라인임.
  - 게놈 QC에 BUSCO가 사용되지 않고 QUAST만 사용됨 — 나머지 4개 프로젝트와의 뚜렷한 차이(사유는 확인 필요).
  - 동일 종을 서로 다른 시퀀싱 공급사 조합(JSLink 단독 vs JSLink+EYEONCELL 통합)으로 두 번 독립적으로 조립하고 QUAST로 비교하는 "데이터 소스 조합 검증" 구조 — 다른 프로젝트들에는 없는 패턴.
  - hifiasm 조립 자체는 organelle 제거 이전(raw 전체 reads) 상태로 1회만 수행되고, 이후 별도의 사후 필터링(minimap2+seqkit, coverage threshold 0.90)으로 organelle contig를 제거함 — Cannabis는 Oatk 분리 후 조립을 재수행(2회 조립)하는 방식인 반면, Abies는 조립을 1회만 수행하고 필터링만 사후에 적용하는 점이 다름.
  - `run_pipeline.sh`/`run_remaining.sh` 마스터 오케스트레이션 스크립트가 참조하는 디렉토리명(`Organelle_combined`, `hifiasm_combined`)과 현재 실제 디렉토리명(`Organelle_JSLink_EYEONCELL`, `hifiasm_JSLink_EYEONCELL`)이 rename으로 인해 불일치. 최초 실행(2026-06-08)에서 OatK(통합본) 병렬 스텝이 실패하여 hifiasm(통합본)만 완료되었고, OatK는 이후 별도 세션(`run_remaining.sh`, 2026-06-10~11)에서 재실행되어 성공함.
  - FastK 버전은 `-h`/`--version` 플래그를 지원하지 않아 확인 불가. Ratatosk류 하이브리드 보정은 이 프로젝트에서 사용되지 않음(PacBio HiFi 단일 플랫폼만 사용).
- (해당 시) 프로젝트 고유 확장 분석:
  - JSLink 단독 vs JSLink+EYEONCELL 통합 조립본 간 QUAST 비교(`quast_compare_{hap1,hap2,p}`) — 추가 시퀀싱 투입 효과 검증용.
  - hifiasm 그래프 기반 OatK organelle 추출(-G 모드)과 standalone OatK 결과 간 비교(`03_compare_organelle.sh`, QUAST 기반) — organelle 분리 방법론 자체의 교차검증(실행 완료 여부는 산출물로 확인되나 실행 주체/시점 기록은 확인 필요).
  - 향후 계획: ONT ULK 데이터 및 제주도 취약/건전 10샘플씩 PacBio Kinnex Isoseq 데이터 생산 예정 — 향후 전사체 기반 유전자 주석(Zoysia_macrostachya의 Braker3 경로와 유사한 방식)으로 이어질 것으로 예상되나 현재는 미착수.

### Zoysia_macrostachya_dnDNAseq
- 연구목적/샘플 설계: *Zoysia macrostachya*(들잔디속 근연종)의 de novo 전장유전체 해독. ONT+Illumina 하이브리드 long-read 조립 → Omni-C(Hi-C) 기반 de novo 스캐폴딩 → 수동 curation(Juicebox, 2회) → 갭 클로징(YagCloser) → 폴리싱 → 반복서열/유전자/기능 주석 → 비교유전체(Pan-genome, MCMCTree 분화시기, CAFE5 유전자가족 진화, PSG 양성선택) → 자체 발현정량 및 공개 SRA 염 스트레스 전사체 재분석까지 수행하는, 5개 프로젝트 중 가장 확장된 파이프라인. 20개 염색체 수준 최종 게놈(`Zma_FINAL.20chr.nogaps.sorted.fasta`)을 NCBI에 제출 완료.
- 공통 파이프라인과의 차이점:
  - 참조 유전체 기반 RagTag scaffolding이 아니라, Purge_dups(중복 제거) → Omni-C(Hi-C) 리드 기반 HapHiC 스캐폴딩 → Juicebox 수동 curation(2회, 1차는 초기 스캐폴드 검토·2차는 재정렬/SV 검토로 20개 염색체 확정)이라는 de novo 스캐폴딩 경로를 사용 — Cannabis/Pichia와 근본적으로 다른 접근. RagTag는 Zoysia_macrostachya에서도 사용되나 이는 자기 게놈이 아니라 타 Zoysia속 게놈과의 비교유전체용 재정렬 목적(1절 Step 8 참조). 같은 속의 Zoysia_sinica도 Hi-C 기반 de novo 스캐폴딩을 쓰지만 실제 툴은 다르다(HapHiC 대신 yahs/SALSA2/3D-DNA — 1절 Step 5 참조).
  - 스캐폴딩 이후 YagCloser로 실제 갭 클로징을 수행하여 최종 게놈에 반영 — Pichia도 Yagcloser를 보유하나 X-33에 한해 분석 전용(gap-filling 커맨드 비활성)으로 그친 반면, Zoysia_macrostachya는 실제 갭을 메워 최종 게놈 생성 경로에 포함시킴. Zoysia_sinica는 YagCloser를 사용한 정황이 없다(대신 갭 제거는 스캐폴더 산출물의 "nogaps" 계열 파일로 처리된 것으로 보임 — 확인 필요).
  - 게놈 QC에 LAI(LTR Assembly Index)와 Merqury가 추가로 사용됨 — Cannabis/Pichia에는 없는 QC 지표(Zoysia_sinica는 Merqury는 쓰지만 LAI 사용 정황은 없음).
  - 비교유전체 규모가 압도적으로 큼: OrthoFinder(12종/5-Zoysia 오솔로그 클러스터링) → MCMCTree(PAML, 분화시기 추정) → CAFE5(유전자 가족 확장/축소) → PSG(PAML branch-site 양성선택 유전자, 7개 계통 반복)로 이어지는 다단계 파이프라인 — Cannabis/Pichia/Abies에는 없다. Zoysia_sinica도 매우 유사한 규모의 비교유전체 스위트(OrthoFinder/CAFE/MCMCTree/Ks·MCScanX/PSG 등)를 별도로 갖고 있다(2절 Zoysia_sinica 항목 참조) — 다만 macrostachya는 이를 하나의 `Pan-genome_analysis/` 디렉토리로 통합한 반면 sinica는 번호가 매겨진 개별 폴더로 나뉘어 있다는 조직 방식의 차이가 있다.
  - Pan-genome_analysis 디렉토리 하위에 RagTag/SyRI/plotsr/CHROMEISTER/MCScan/CandiSSR/CNV_analysis/Purge_dups/Mercury/BUSCO 등 다수의 병렬 서브분석이 있어(파일 수 13만 6천여 개), Cannabis/Pichia의 단일 synteny 비교보다 훨씬 포괄적.
  - 자체 조직 발현정량(Hisat2+StringTie+TMM 정규화) 및 공개 SRA 염 스트레스 전사체 재분석(Trimmomatic→Prinseq→Hisat2→edgeR DEG)까지 수행 — 다른 프로젝트들에서는 RNA-seq이 유전자예측용 hint 생성에만 쓰이는 반면, Zoysia_macrostachya는 발현량 정량 및 차등발현유전자(DEG) 분석까지 확장.
  - BRAKER3/EnTAP/OrthoFinder/MCMCTree/CAFE5 등 핵심 스텝 다수가 실행 스크립트/wrapper 없이 로그·설정 파일(ctl, .ini, .log)만 남아있어 재현성이 낮음 — 확인 필요 항목이 Cannabis/Pichia/Abies보다 현저히 많다(다만 Zoysia_sinica는 원본 조립 디렉토리 자체가 접근 불가였던 만큼 확인 필요 항목이 이 프로젝트보다도 더 많다).
  - 06.Repeat_element_annotation의 RepeatModeler 실행은 "before_curation"(1차 manual curation 이전 게놈, 2024-10-26)에 대해서만 확인되고, 최종 20chr.sorted 버전에 대한 재실행 디렉토리는 발견되지 않음(날짜-스크립트 불일치) — 확인 필요.
  - Purge_dups(01.Assembly)의 input이 정확히 어느 NextPolish 산출물 계보(Ratatosk 경유 vs 비-Ratatosk 경유, 파일명이 동일하여 구분 불가)인지 확인 필요.
- (해당 시) 프로젝트 고유 확장 분석:
  - Pan-genome_analysis: RagTag/SyRI/plotsr(구조변이), CHROMEISTER/MCScan(synteny), CandiSSR(SSR 마커), CNV_analysis(카피수변이), Purge_dups/Mercury/BUSCO(QC 반복) 등을 통한 근연 Zoysia종(Z. japonica, Z. matrella, Z. pacifica, Z. sinica 등) 다중 비교.
  - MCMCTree(PAML)를 이용한 종간 분화시기 추정.
  - CAFE5를 이용한 유전자 가족 확장/축소(유전자 가족 진화) 분석.
  - PSG(양성선택 유전자) 분석: translatorX → pal2nal → Gblocks → codeml branch-site model → LRT 유의성 검정, 계통(Zma/Zj/Zmt/Zp/Zs/Zoysia/ZsZma)별 반복 수행 → `Zoysia_PSGs_Orthogroups.tsv`.
  - 자체 조직 발현정량(11.Transcriptome, gene/transcript count matrix + TMM 정규화 TPM) 및 공개 SRA 기반 염 스트레스 전사체 재분석(13.SRA_Salt_Transcriptome, edgeR 기반 DEG 분석).
  - 텔로미어 영역 분석(04.Telomere_region, quarTeT/soft-clip 리드 기반 telomere read 추출·패칭) — 최종 게놈 계보와의 직접 연결점은 확인 필요(독립 QC 서브분석으로 추정).

### Zoysia_sinica_Genome_Assembly
- 연구목적/샘플 설계: **description.txt/README가 전혀 없어(확인됨) 아래 내용은 전부 폴더/파일명 기반 추정이다.** 폴더명(`Zoysia_genus_Project/Zoysia_sinica/Genome_Assembly`)과 파일명 패턴(`Z.sinica_IO_genome...`, Z. japonica 대비 chr01~chr20 재정렬, 근연종 다수 비교)으로 미루어 *Zoysia sinica*(갯벌잔디)의 de novo 전장유전체 조립 및 주석, 그리고 Zoysia_macrostachya와 유사한 규모의 비교유전체 확장 분석을 목표로 하는 프로젝트로 추정된다. 원시 de novo assembly 자체(organelle 분리 → assembly → 1차 QC → MAKER → 1차 기능주석 → 1차 RepeatModeler/Masker)는 `Original/ONT_longread_DNAseq/` 디렉토리에서 수행된 것으로 보이나, 이 디렉토리는 조사 시점에 접근 불가(간헐적 NFS 이슈)였다 — 어셈블러/시퀀싱 플랫폼조차 확정할 수 없으며, `Z.sinica_IO_genome.original.fasta`라는 contig-level 파일의 존재와 그에 대한 2023-11경 RepeatModeler/Masker 실행 정황만으로 "이 시점 이전에 어떤 형태로든 조립이 완료되었다"는 것만 간접적으로 확인된다.
- 공통 파이프라인과의 차이점:
  - Step 1(genome survey)에 해당하는 스크립트/로그가 확인되지 않고, Step 2(de novo assembly) 자체도 원본 디렉토리 접근 불가로 검증 불가 — "확인 필요"가 아니라 "검증 시도 자체가 불가능했음"에 가깝다. 이는 나머지 4개 프로젝트와 근본적으로 다른 지점이다.
  - Hi-C 기반 de novo 스캐폴딩을 쓴다는 점은 Zoysia_macrostachya와 같지만, 실제 스캐폴딩 툴은 yahs/SALSA2/3D-DNA+Juicer 세 가지를 병행 시도한 뒤 yahs 계열을 최종 채택한 것으로, macrostachya의 HapHiC과는 다른 도구다(1절 Step 5 참조). "왜 yahs를 택했는지"에 대한 명시적 의사결정 기록은 없고, 이후 모든 하위 스텝이 yahs 계열 파일명을 사용한다는 정황 증거만 있다.
  - RagTag를 사용하긴 하나 Cannabis/Pichia처럼 1차 scaffolding 목적이 아니라, Z. japonica 기준 chr01~chr20 명명/순서·역상보 보정("renaming/reordering") 목적으로 두 차례(Hi-C 스캐폴딩 직후 1차, Manual_curation 이후 최종 게놈에 대해 2차) 사용된다(1절 Step 4 상단 참고 블록, Step 5, Step 10 참조).
  - 유전자 예측은 BRAKER **v3.0.7**(ETP mode, 로그에서 버전 직접 확인)으로 Zoysia_macrostachya와 같은 BRAKER3 계열이지만, GeneMarkS-T/TSEBRA 통합 대신 AGAT v1.4.0 기반 rename/longest-isoform 후처리를 거친다는 점이 다르다.
  - Zoysia_macrostachya에는 없는 별도의 **수동 큐레이션(Manual_curation)** 스텝이 있다(1절 Step 10) — BRAKER+AGAT 산출물을 사람이 직접 검토해 최종 병합 게놈을 만드는 과정으로, 병합 로직 자체를 설명하는 스크립트는 남아있지 않다(AGAT 로그만 확인).
  - 나머지 4개 프로젝트 전체를 통틀어 가장 "확인 필요"가 많은 프로젝트다: 원본 조립 디렉토리 접근 불가, Step 2~Step 5·Step 10 다수 스텝의 원본 실행 스크립트(.sh) 자체가 보존되지 않고 로그/후속 스크립트로만 재구성됨, 대부분 툴의 정확한 버전 미확인(BRAKER 3.0.7과 AGAT v1.4.0만 로그로 직접 확인됨), 실행 환경(conda/컨테이너)·스케줄러 사용 여부 대부분 미확인, top-level `Zsinica_genome.fasta`의 정확한 출처 미확인. 이는 문서화 결함이 아니라 실제로 접근 가능했던 증거의 한계를 정직하게 반영한 것이다.
- (해당 시) 프로젝트 고유 확장 분석: de novo 조립 이후의 비교유전체·진화유전체·기능주석/실험검증
  후속 분석은 별도 카테고리 문서로 상세 정리했다 — [comparative_genomics.md](comparative_genomics.md)
  (OrthoFinder, Ks/MCScanX, Circos, 근연종 재주석), [evolutionary_genomics.md](evolutionary_genomics.md)
  (CAFE 유전자가족진화, MCMCTree 분화시기, PAML 양성선택), [genome_annotation_and_validation.md](genome_annotation_and_validation.md)
  (EnTAP 기능주석, 침수내성 후보유전자, TeloExplorer/CentroMiner 구조유전체, PlantCARE, qRT-PCR
  검증) 참조.

### Rubus_Pangenome_Project
- 연구목적/샘플 설계: *Rubus*(산딸기속) 4개 근연종(R. coreanus=Rco, R. crataegifolius=Rcra,
  R. phoenicolasius=Rpho, R. takesimensis=Rta) 각각의 haplotype-resolved(hap1/hap2) de novo
  전장유전체 조립 및 판게놈 비교유전체 연구. 각 종은 PacBio Revio HiFi(1 Cell) 롱리드 + Omni-C
  (Hi-C) paired-end + mRNA-seq 원시데이터를 사용한다.
- 공통 파이프라인과의 차이점:
  - k-mer 기반 게놈 서베이 툴이 종마다 다르다 — Rco/Rcra/Rpho는 FastK+GeneScopeFK, Rta는
    jellyfish+GenomeScope v1.0을 사용한다(각 종 실행 스크립트에서 직접 확인).
  - 오가넬(엽록체/미토콘드리아) 분리는 4종 모두 oatk v1.0(OatkDB v20230921 magnoliopsida
    mito/pltd 모델)을 사용하며, 실행 커맨드가 각 종 로그에서 동일한 패턴으로 확인된다.
  - Rco/Rcra/Rpho 3종은 **hifiasm(purge level 3, `-l 3`) 조립 → Omni-C 리드 정렬/필터링 → YaHS
    스캐폴딩(`--no-contig-ec --no-scaffold-ec`, 이하 no_ec)**으로 이어지는 동일한 파이프라인을
    사용한다. Rta만 이 경로 대신 AutoHiC(Juicer+3D-DNA+딥러닝 오류교정/염색체배정 모델 통합
    파이프라인, hap1만 수행)를 사용한다.
- 실행 커맨드(스크립트/로그에서 직접 확인된 것):
  ```bash
  # 오가넬 분리 (4종 공통 패턴, <species>=Rco/Rcra/Rpho/Rta)
  oatk -t 100 --nhmmscan /usr/local/bin/nhmmscan \
      -m /usr/local/bin/OatkDB/v20230921/magnoliopsida_mito.fam \
      -p /usr/local/bin/OatkDB/v20230921/magnoliopsida_pltd.fam \
      -o <species>_organelle <species>.fasta.gz

  # hifiasm 조립, Omni-C reads로 haplotype 분리, purge level 3 (Rco/Rcra/Rpho 공통 패턴)
  hifiasm -t 100 -l 3 --h1 <species>_Omni-C_R1.fastq.gz --h2 <species>_Omni-C_R2.fastq.gz \
      -o <species>.asm_Omni-C --telo-m TTTAGGG <species>_DNA.hifi.fastq.gz

  # Omni-C 리드 정렬 + 필터링 (bwa mem + samblaster + HapHiC filter_bam 유틸리티)
  bwa mem -5SP -t 100 <hap>.p_ctg.fa <species>_Omni-C_R1.fastq.gz <species>_Omni-C_R2.fastq.gz \
      | samblaster | samtools view -F 3340 -o <hap>_HiC.bam
  /usr/local/bin/HapHiC/utils/filter_bam <hap>_HiC.bam 1 --nm 3 --threads 100 \
      | samtools view -b -o <hap>_HiC.filtered.bam

  # YAHS 스캐폴딩 (no_ec: 오류보정 비활성화, hap1/hap2 각각 실행)
  yahs -o <hap>_no_ec --telo-motif AAACCCT --no-contig-ec --no-scaffold-ec --file-type BAM \
      <hap>.p_ctg.fa <hap>_HiC.filtered.bam

  # Rta AutoHiC — Juicer 정렬
  bwa mem -SP5M -t 45 Rta_hap1.fa Rta_Omni-C_R1.fastq.gz Rta_Omni-C_R2.fastq.gz > Rta_Omni-C.fastq.gz.sam
  # Rta AutoHiC — 3D-DNA 스캐폴딩
  run-asm-pipeline.sh -r 2 Rta_hap1.fa merged_nodups.txt
  ```
  Rta의 AutoHiC 최종 산출물(`result/AutoHiC/Rta_hap1_autohic.fasta`)에 대한 QUAST 결과: 5개
  contig, Total length 254,145,080 bp, N50 130,016,229 bp, GC 37.79%(`report.txt`에서 직접 확인).
- (해당 시) 프로젝트 고유 확장 분석: de novo 조립 이후의 비교유전체·진화유전체·기능주석/실험검증
  후속 분석은 별도 카테고리 문서로 정리했다 — [comparative_genomics.md](comparative_genomics.md),
  [evolutionary_genomics.md](evolutionary_genomics.md), [genome_annotation_and_validation.md](genome_annotation_and_validation.md)
  참조.

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Cannabis_sativa_cv.ChungSam_dnDNAseq | PacBio HiFi 단일 플랫폼, organelle(Oatk) 분리 후 hifiasm으로 haplotype-resolved(hap1/hap2) diploid 조립, 판게놈 유래 참조(AH3Ma/AH3Mb)로 RagTag scaffolding, RepeatModeler/RepeatMasker+BRAKER3로 de novo 유전자 예측(scaffolding 이전 contig 계열 입력) | hifiasm 0.25.0-r726, RagTag v2.1.0, BUSCO v5.8.2, QUAST v5.2.0, RepeatModeler 2.0.6/RepeatMasker 4.1.9, BRAKER 3.0.8, sourmash 4.9.4(그 외 FastK/GenomeScope2/minimap2/samtools/syri/plotsr/CHROMEISTER/BLAST+ 버전은 확인 필요) | Sourmash+chromeister 판게놈(98샘플) 비교, plotsr/syri 구조변이 시각화, THCAS/CBDAS/CBGAS 카나비노이드 유전자 분석(BRAKER3 단백질 tblastn 비교 + 독립적 targeted gene assembly 파이프라인) |
| Pichia_Pastoris_dnDNAseq_HanHwa | ONT+Illumina 하이브리드(Ratatosk 보정), haploid 단일 genome 조립(NextDenovo+NextPolish 채택, trycycler multi-assembler 경로는 reconcile 실패로 폐기), NCBI GS115 단일 참조로 RagTag scaffolding, LiftOff로 참조 기반 annotation liftover(de novo 예측 아님) | NextDenovo v2.5.2, NextPolish v1.4.1, BUSCO v5.8.2, QUAST v5.2.0, mummer 4.0.1, minimap2 v2.24-r1122(RagTag 내부)(RagTag 자체 버전은 확인 필요) | Merqury k-mer QV/completeness 평가, mummer 기반 공개 균주 synteny 비교, Yagcloser gap 분석(X-33만, 분석 전용), Genome_editing gRNA BLAST 매핑(BG-10만) |
| Abies_koreana_dnDNAseq | PacBio Revio HiFi(JSLink+EYEONCELL 2개사 8cell), hifiasm으로 raw reads 전체를 1회 조립 후 사후 organelle contig 필터링(Oatk+minimap2/seqkit), scaffolding·비교유전체·유전자주석은 아직 미수행(진행 중 프로젝트), BUSCO 없이 QUAST만으로 QC | hifiasm v0.25.0-r726, Oatk v1.0, minimap2 2.30-r1287, seqkit v2.3.0, QUAST 5.3.0(fb88221c)(FastK 버전은 `-h` 미지원으로 확인 불가) | JSLink 단독 vs JSLink+EYEONCELL 통합 조립 비교(QUAST), hifiasm 그래프 기반 vs standalone OatK organelle 추출 방법론 교차검증 |
| Zoysia_macrostachya_dnDNAseq | ONT+Illumina 하이브리드(Ratatosk), NextDenovo+NextPolish+ntLink 조립, Purge_dups+Omni-C(Hi-C) 기반 HapHiC de novo 스캐폴딩+Juicebox 수동 curation(2회)+YagCloser 갭클로징으로 20염색체 확정, RepeatModeler/RepeatMasker+BRAKER3+TSEBRA 유전자예측, LAI/Merqury 추가 QC | Ratatosk/NextDenovo/NextPolish/ntLink(버전 확인 필요), HapHiC(버전 확인 필요), juicer_tools 1.9.9, yagcloser(버전 확인 필요), RepeatModeler/RepeatMasker(버전 확인 필요), BRAKER3(버전/컨테이너 확인 필요), LTR_retriever, Merqury | OrthoFinder(12종/5-Zoysia)→MCMCTree(분화시기)→CAFE5(유전자가족진화)→PSG(PAML branch-site 양성선택), Pan-genome_analysis(RagTag/SyRI/CHROMEISTER/MCScan/CandiSSR/CNV 등 근연종 다중 비교), 자체 발현정량(TMM 정규화) + 공개 SRA 염스트레스 전사체 재분석(edgeR DEG), 텔로미어 영역 분석 |
| Zoysia_sinica_Genome_Assembly | description.txt 없음(연구목적 전부 폴더/파일명 기반 추정) — 원시 de novo assembly는 접근 불가 디렉토리(`Original/`)에서 수행되어 검증 불가(확인 필요), Omni-C(Hi-C) 스캐폴딩을 yahs/SALSA2/3D-DNA 3종 병행 시도 후 yahs 채택(macrostachya의 HapHiC과 다른 툴), RagTag를 renaming/reordering 목적(Z. japonica 대비 chr01-20, 중간+최종 2회)으로 사용(1차 scaffolding 아님), RepeatModeler/RepeatMasker+BRAKER3(ETP)+AGAT rename/longest-isoform 유전자예측, BUSCO(genome+protein)/Merqury QC, 그리고 다른 4개 프로젝트에는 없는 수동 큐레이션(Manual_curation) 스텝 보유 | yahs/SALSA-2.3/3D-DNA+Juicer(버전 대부분 확인 필요), RagTag(버전 확인 필요), BUSCO/QUAST/Merqury(버전 확인 필요), RepeatModeler/RepeatMasker(버전 확인 필요), **BRAKER v3.0.7**(로그에서 직접 확인), **AGAT v1.4.0**(로그에서 직접 확인) | 비교유전체·진화유전체·기능주석/검증 후속 분석은 별도 문서로 분리 — comparative_genomics.md, evolutionary_genomics.md, genome_annotation_and_validation.md 참조 |
| Rubus_Pangenome_Project | 4개 근연종(Rco/Rcra/Rpho/Rta) 각 hap1/hap2 haplotype-resolved 조립. Rco/Rcra/Rpho는 hifiasm(purge level 3) → Omni-C 리드 정렬/필터링(bwa+samblaster+HapHiC filter_bam) → YaHS(no_ec) 스캐폴딩으로 동일 파이프라인 사용, oatk로 오가넬 분리(4종 공통). Rta만 AutoHiC(Juicer+3D-DNA+딥러닝) 경로 사용 | hifiasm v0.25.0-r726(로그에서 직접 확인, 4종 공통), oatk v1.0(4종 공통), bwa 0.7.19-r1273/samblaster 0.1.26(Rco/Rcra/Rpho), YaHS(Rco/Rcra/Rpho), FastK(Rco/Rcra/Rpho) vs jellyfish+GenomeScope v1.0(Rta), Juicer v1.6+BWA 0.7.17-r1188+3D-DNA 180922(Rta AutoHiC, 로그에서 직접 확인) | de novo 조립 이후 후속 분석은 별도 문서로 분리 — comparative_genomics.md, evolutionary_genomics.md, genome_annotation_and_validation.md 참조 |
