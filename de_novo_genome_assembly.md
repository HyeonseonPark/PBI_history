# De novo Genome Assembly (Cannabis_sativa_cv.ChungSam_dnDNAseq, Pichia_Pastoris_dnDNAseq_HanHwa, Abies_koreana_dnDNAseq, Zoysia_macrostachya_dnDNAseq)

## 1. 공통 파이프라인

4개 프로젝트 모두 "원시 long-read(+보조 short-read) → de novo assembly → genome QC" 라는 핵심 골격은
공유하지만, 그 이후 단계(scaffolding 방식, 비교유전체, 유전자주석 수행 여부)는 프로젝트마다 실제로
수행된 범위가 크게 다르다. 특히:

- **Cannabis, Pichia, Zoysia** 3개는 "assembly → scaffolding → genome QC → 비교유전체(synteny/dot-plot)
  → 유전자/반복서열 주석"까지 완결된 파이프라인을 보유한다. 다만 scaffolding 방식은 셋 다 다르다
  (Cannabis/Pichia는 참조 유전체 기반 RagTag, Zoysia는 Omni-C/Hi-C 기반 de novo 스캐폴딩 + 별도 갭
  클로징).
- **Abies_koreana**는 2026-07 기준 아직 초기 단계 프로젝트로, k-mer profiling → de novo assembly →
  organelle contig 분리(사후 필터링) → 1차 QUAST QC까지만 완료되어 있다. scaffolding, 비교유전체,
  유전자주석 단계는 아직 수행되지 않았으며, BUSCO도 사용되지 않고 QUAST만으로 QC가 이루어진다. 즉
  Abies는 나머지 3개 프로젝트와 "공통 뼈대"의 앞부분(조립까지)만 공유하고 뒷부분은 아직 해당사항이
  없다 — 이는 프로젝트 진행 단계 차이이지 방법론적 차이가 아니다.
- k-mer 기반 genome survey(FastK+GenomeScope) 단계는 Cannabis/Pichia/Abies 3곳에서 확인되나,
  Zoysia는 이 단계에 해당하는 스크립트/로그가 발견되지 않는다(genome_size=0.323g 값이 NextDenovo
  cfg에 직접 기재되어 있으나 그 산출 근거는 확인 필요).
- organelle contig 분리(Oatk) 단계는 Cannabis와 Abies에만 있다(둘 다 PacBio HiFi 단일 롱리드
  플랫폼을 사용하는 프로젝트라는 공통점이 있음).

아래 DAG는 이 4개 프로젝트에서 실제로 확인된 단계만을 표시하며, 특정 프로젝트에만 해당하는 단계는
괄호로 명시한다.

```
[rawdata: long-read (+ short-read)]
        │
        ▼
[Step 1. Genome survey / k-mer profiling: FastK + GenomeScope(2)]  (Cannabis/Pichia/Abies만 해당 — Zoysia는 미확인)
        │  (게놈 크기·이형접합률 추정치 → 이후 assembly 파라미터로 재사용)
        ▼
[Step 2. De novo assembly (프로젝트별 도구/전략 상이 — 2절 참조)]
        │
        ├──▶ [Step 3. Organelle contig 분리(Oatk) — Cannabis/Abies만 해당]
        │
        ├──▶ [Step 4. 참조 기반 scaffolding: RagTag — Cannabis/Pichia만 해당]
        │
        ├──▶ [Step 5. Hi-C(Omni-C) 기반 de novo 스캐폴딩 — Zoysia만 해당]
        │            │
        │            ▼
        │       [Step 6. Gap closing: YagCloser — Zoysia만 해당]
        │
        ▼
[Step 7. Genome QC: QUAST(4개 공통) + BUSCO(Cannabis/Pichia/Zoysia만, Abies는 미수행) + LAI/Merqury(Zoysia만)]
        │
        ├──▶ [Step 8. 비교유전체 분석 (synteny/dot-plot/pan-genome, Abies는 아직 미수행)]
        │
        └──▶ [Step 9. 유전자/반복서열 주석 (de novo 예측 또는 참조 기반 liftover, Abies는 아직 미수행)]
                        │
                        ▼
                [프로젝트 고유 확장 분석 / 최종 산출물]
```

### Step 1. Genome survey (k-mer profiling)
- **소프트웨어(버전):** FastK, GenomeScope2 (Cannabis: 버전 확인 필요; Pichia: 사용 확인되었으나 세부 버전은 재현 우선순위 낮음으로 처리됨; Abies: FastK/Fastmerge/Histex/GenomeScope 사용 확인되었으나 FastK가 `-h`/`--version` 플래그를 지원하지 않아 버전 확인 불가) — Zoysia는 이 단계에 해당하는 스크립트/로그 없음(확인 필요)
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
  - Zoysia: 확인 필요 — 스크립트 파일 없음

### Step 2. De novo assembly
- **소프트웨어(버전):** Cannabis: hifiasm 0.25.0-r726 / Pichia: NextDenovo v2.5.2 + NextPolish v1.4.1 / Abies: hifiasm v0.25.0-r726(Cannabis와 동일 버전) / Zoysia: Ratatosk(하이브리드 보정, 버전 확인 필요) + NextDenovo(버전 확인 필요) + NextPolish(버전 확인 필요) + ntLink(스캐폴딩, 버전 확인 필요)
- **핵심 파라미터:** Cannabis: `-t 100`(스레드) / Pichia: `input_type=corrected, read_type=ont, genome_size=<GenomeScope 추정치>`(NextDenovo), `task=best, sgs_options=-max_depth 100 -bwa`(NextPolish, Illumina 폴리싱만 사용) / Abies: `-t 128`(JSlink 4cell 단독) 또는 `-t 256`(JSlink+EYEONCELL 8cell 통합) / Zoysia: Ratatosk `-c 256`; NextDenovo run.cfg(`read_type=ont, genome_size=0.323g, pa_correction=7, minimap2_options_cns=-t 100 -k 19, nextgraph_options=-a 1`, `job_type=local`); NextPolish(`sgs_options=-max_depth 100 -bwa, polish_options=-10`); ntLink(`k=40 w=500`, 7개 k/w 조합 중 최종 채택)
- **Input:** Cannabis: organelle 제거 후 nuclear-only HiFi reads / Pichia: ONT reads를 Illumina로 하이브리드 교정(Ratatosk)한 corrected long reads + Illumina reads(폴리싱용) / Abies: organelle 제거 이전 raw 전체 HiFi reads(JSlink 4개 cell, 통합본은 EYEONCELL 1개 cell 추가) / Zoysia: ONT reads(Z.macrostackya_ONT.fastq.gz)를 Illumina reads(Zma-H1_1/2.fastq.gz)로 Ratatosk 보정한 out_long_reads.fastq
- **Output:** Cannabis: haplotype-resolved contig(hap1/hap2 × Female/Male) / Pichia: 폴리싱된 단일(haploid) genome fasta / Abies: 단일 개체(제주도 표준목)의 hifiasm bp.{p_ctg,hap1.p_ctg,hap2.p_ctg} contig(organelle 포함, 이후 Step3에서 제거) / Zoysia: genome.nextpolish.fasta → ntLink scaffold(k40w500) 결과물인 `Ratatosk_NextDP_NextPsh_ntLinkk40w500.fasta`(1차 QC용 조립본, 이후 Purge_dups/Omni-C 스캐폴딩으로 이어짐)
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
  - Zoysia (01.Assembly/Ratatosk/Ratatosk.sh + run.cfg + run_NextPolish.cfg + ntLink/ntLink.sh):
    ```bash
    Ratatosk correct -v -c 256 -s short_reads.fastq.gz -l in_long_reads.fastq.gz -o out_long_reads
    nextDenovo run.cfg   # read_type=ont, genome_size=0.323g, pa_correction=7, minimap2_options_cns=-t 100 -k 19
    nextPolish run_NextPolish.cfg   # sgs_options=-max_depth 100 -bwa, polish_options=-10
    ntLink scaffold target=genome.nextpolish.fasta reads=out_long_reads.fastq t=256 overlap=True k=40 w=500
    ```

### Step 3. Organelle contig 분리(Oatk) — Cannabis, Abies만 해당
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

> 참고: Zoysia도 RagTag를 사용하지만(Pan-genome_analysis/RagTag/), 이는 자기 게놈의 scaffolding이
> 아니라 타 Zoysia속 게놈과의 비교유전체용 재정렬이다. Zoysia 자신의 scaffolding은 Step 5(Hi-C 기반)를
> 통해 이루어지므로 RagTag 사용은 Step 8(비교유전체 분석)에 별도로 기재한다. Abies는 아직 scaffolding을
> 수행하지 않았다.

### Step 5. Hi-C(Omni-C) 기반 de novo 스캐폴딩 — Zoysia만 해당
- **소프트웨어(버전):** purge_dups(버전 확인 필요, 실행 스크립트 부재로 config.json만 확인) → bwa + pairtools + samtools(Omni-C 리드 검증) → HapHiC(`haphic pipeline`, 버전 확인 필요) → Juicebox GUI(수동 curation, juicer_tools 1.9.9/jcuda 0.8)
- **핵심 파라미터:** Purge_dups `config.json`: `isdip=1, ispb=1`, busco `lineage=Poales`; bwa `-5SP -T0 -t 150`; pairtools parse `--min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30`; HapHiC `haphic pipeline purged.fa HiC.filtered.bam 20 --correct_nrounds 2 --max_inflation 5`(20=예상 염색체 수); juicer_tools `-Xmx32G`
- **Input:** `genome.nextpolish.fasta`(Step2 NextPolish 산출물 — 다만 config.json의 `ref` 경로가 Ratatosk 계보/비-Ratatosk 계보 중 어느 쪽인지 파일명만으로 확정 불가, 확인 필요) + Omni-C 원시 리드(`Zma_S222_L002_R1/2_001.fastq.gz`)
- **Output:** `purged.fa`/`hap.fa`/`dups.bed`(Purge_dups) → `mapped.PT.bam`/`HiC.filtered.bam`(정렬·필터링) → `scaffolds.raw.agp`(HapHiC) → `out_JBAT.FINAL.fa`(1차 Juicebox 수동 curation) → `out_JBAT_SV.FINAL.20chr.fa`(2차 Juicebox 수동 curation, 20개 염색체 확정)
- **실행 커맨드:**
  - Zoysia (Purge_dups는 config.json만 존재, 실행 스크립트 자체는 확인 필요 — 스크립트 파일 없음)
  - Zoysia (Omni-C 리드 검증, Validating_Omni-C_reads.sh):
    ```bash
    bwa mem -5SP -T0 -t 150 purged.fa Zma_S222_L002_R1_001.fastq.gz Zma_S222_L002_R2_001.fastq.gz -o aligned.sam
    pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 150 --nproc-out 150 --chroms-path Zma_purged.genome aligned.sam > parsed.pairsam
    pairtools sort --nproc 150 parsed.pairsam > sorted.pairsam
    pairtools dedup --nproc-in 150 --nproc-out 150 --mark-dups --output-stats dup_stats.txt --output dedup.pairsam sorted.pairsam
    pairtools split --nproc-in 150 --nproc-out 150 --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam
    samtools sort -@150 -o mapped.PT.bam unsorted.bam && samtools index mapped.PT.bam
    ```
  - Zoysia (HapHiC 스캐폴딩, HapHiC/HapHiC.sh — 전처리(bwa+samblaster+filter_bam) 커맨드는 스크립트 내 주석 처리되어 비활성 상태이며 `HiC.filtered.bam`은 사전 생성된 파일을 그대로 재사용):
    ```bash
    haphic pipeline purged.fa HiC.filtered.bam 20 --correct_nrounds 2 --threads 150 --processes 150 --max_inflation 5
    ```
  - Zoysia (Juicebox 수동 curation): juicer pre/juicer_tools로 `.hic` 파일 생성 후 외부 JuiceBox GUI에서 수동 검토·편집 — GUI 조작 자체는 스크립트화되어 있지 않음(확인 필요)

### Step 6. Gap closing(YagCloser) — Zoysia만 해당
- **소프트웨어(버전):** minimap2 + samtools(정렬), detgaps, yagcloser.py(`/usr/local/bin/yagcloser/`, 버전 확인 필요)
- **핵심 파라미터:** minimap2 `-x map-ont -t200 --secondary=yes`(스크립트 내 주석 처리 — 실제 실행 시에는 사전 생성된 `aln.s.bam`을 재사용)
- **Input:** `out_JBAT_SV.FINAL.20chr.fa`(Step5 산출물, 20염색체 확정본) + `ONT_corrected.fastq`(원시 ONT 보정 리드)
- **Output:** `out_JBAT_SV.FINAL.20chr.nogaps.fasta`(갭 클로징 완료) → (염색체 정렬/개명 처리, 실행 스크립트 확인 필요) → `Zma_FINAL.20chr.nogaps.sorted.fasta`(이후 모든 하류 스텝의 최종 게놈)
- **실행 커맨드:**
  - Zoysia (03.YagCloser/After_manualcuration/yagcloser.sh):
    ```bash
    samtools index aln.s.bam
    detgaps out_JBAT_SV.FINAL.20chr.fa > out_JBAT_SV.FINAL.20chr.gaps.bed
    yagcloser.py -g out_JBAT_SV.FINAL.20chr.fa -b out_JBAT_SV.FINAL.20chr.gaps.bed -a aln.s.bam -s Zma_superscaffolds -o yagcloser_outdir
    python /usr/local/bin/yagcloser/scripts/update_assembly_edits_and_breaks.py -i out_JBAT_SV.FINAL.20chr.fa -o out_JBAT_SV.FINAL.20chr.nogaps.fasta -e Zma_scaffolds.edits.txt
    ```
  - 이후 염색체 정렬·이름 확정(`chr.id` 매핑 기반으로 추정)하여 `Zma_FINAL.20chr.nogaps.sorted.fasta` 생성 — 해당 정렬 스크립트 자체는 확인 필요 — 스크립트 파일 없음

### Step 7. Genome QC
- **소프트웨어(버전):** BUSCO v5.8.2(Cannabis/Pichia 공통), QUAST v5.2.0(Cannabis/Pichia 공통); Abies: QUAST 5.3.0(fb88221c)만 사용(BUSCO 미수행 — 확인 필요, 프로젝트 초기 단계로 추정); Zoysia: BUSCO(embryophyta/liliopsida/poales_odb10) + QUAST 5.3.0 + LAI(LTR_FINDER_parallel+ltrharvest+LTR_retriever) + Merqury(meryl k=21)
- **핵심 파라미터:** BUSCO: `-m genome`, lineage(Cannabis: rosales_odb12/rosids_odb10, Pichia: ascomycota_odb12/saccharomycetes_odb12, Zoysia: embryophyta_odb10/liliopsida_odb10/poales_odb10); QUAST: 참조 유전체 지정, 기본 옵션(Abies는 `--large --threads 64 --labels ...` 또는 `--labels JSlink_only,JSlink_EYEONCELL`); Zoysia LAI: `LTR_retriever -threads 100`; Zoysia Merqury: `meryl k=21`
- **Input:** scaffold(chr 단위) fasta(Cannabis/Pichia/Zoysia); Abies: organelle 제거된 nuclear contig(hap1/hap2/p, scaffolding 이전 단계)
- **Output:** BUSCO short_summary(.txt/.json), QUAST report(.html/.tsv/.pdf), Zoysia는 추가로 LAI score(`*.out.LAI`), Merqury QV/completeness
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
  - Zoysia (1차 어셈블리 QC, 01.Assembly/QC_results/BUSCO/BUSCO.sh + Quast/Quast.sh):
    ```bash
    busco -m genome -i Ratatosk_NextDP_NextPsh_ntLinkk40w500.fasta -o Ratatosk_NextDP_NextPsh_ntLinkk40w500.embryophyta_odb10 -c 256 --datasets_version odb10 -f -l embryophyta_odb10
    quast.py Ratatosk_NextDP_NextPsh_ntLinkk40w500.fasta -o Ratatosk_NextDP_NextPsh_ntLinkk40w500 -t 256
    ```
  - Zoysia (최종 게놈 LAI/Merqury, 05.Genome_QC/LAI/.../LAI.sh + Mercury/.../mercury.sh):
    ```bash
    LTR_retriever -genome Zma_chromsome_nogap.sorted.fixed.fa -inharvest Zma_chromsome_nogap.sorted.fa.rawLTR.scn -threads 100
    merqury.sh Zma.meryl Zma_FINAL.20chr.nogaps.sorted.fasta Zma_chromosome
    ```

### Step 8. 비교유전체 분석(synteny/dot-plot/pan-genome)
- **소프트웨어(버전):** Cannabis: chromeister + minimap2/syri/plotsr(버전 확인 필요) / Pichia: mummer 4.0.1(nucmer+mummerplot) / Zoysia: minimap2+samtools+SyRI+plotsr(구조변이), RagTag(타 Zoysia속 게놈 비교용 재정렬), OrthoFinder(오솔로그 클러스터링, 버전 확인 필요), CHROMEISTER/MCScan/CandiSSR/CNV_analysis/Purge_dups/Mercury/BUSCO 등(Pan-genome_analysis 하위 병렬 서브모듈군)
- **핵심 파라미터:** Cannabis: minimap2 `-ax asm5 --eqx`, syri `-F B` / Pichia: nucmer `--maxmatch -l 20 -g 500 -c 100` / Zoysia: minimap2 `-ax asm5 -t 200 --eqx`, syri `-F B --nc 200`, RagTag `--mm2-params '-x asm5' -w -t 150`
- **Input:** scaffold(chr 단위) fasta 및 비교 대상 유전체(Cannabis: 판게놈 98샘플; Pichia: 공개 Pichia 균주 게놈; Zoysia: `Zma_FINAL.20chr.nogaps.sorted.fasta` + 타 Zoysia속 게놈 ZJN/ZMW/ZPZ/Z.sinica 등, 원본 출처는 확인 필요)
- **Output:** dot-plot 이미지(.png), 구조변이 vcf/summary(Cannabis, Zoysia); Zoysia는 추가로 RagTag scaffold, OrthoFinder Orthogroups.tsv 등(Step 9의 MCMCTree/CAFE5/PSG 입력으로 재사용, 2절 참조)
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
  - Zoysia (구조변이, Pan-genome_analysis/syri/syri.sh):
    ```bash
    minimap2 -ax asm5 -t 200 --eqx Z.sinica_manualcurated.FINAL_merged2.fasta Zma_FINAL.20chr.nogaps.sorted.fasta | samtools sort -O BAM - > curated_Zs2_Zma.bam
    samtools index curated_Zs2_Zma.bam
    syri -c curated_Zs2_Zma.bam -r Z.sinica_manualcurated.FINAL_merged2.fasta -q Zma_FINAL.20chr.nogaps.sorted.fasta -F B --prefix curated_Zs2_Zma --nc 200
    plotsr --sr curated_Zs2_Zmasyri.out --sr Zma_Zjsyri.out --genomes genome.txt -o curated_Zs2-Zma-Zj_plot.png
    ```
  - Zoysia (타 Zoysia속 게놈 비교용 재정렬, Pan-genome_analysis/RagTag/RagTag.sh):
    ```bash
    ragtag.py scaffold Zsg_Zma_ZJN_chr_genome.fa <타종>_r1.0.genome.fasta -o <타종>_ZsgZmaZJN_scaffold_asm5 -w -t 150 --mm2-params '-x asm5'
    ```
  - Zoysia (OrthoFinder/CHROMEISTER/MCScan/CandiSSR/CNV/Purge_dups/Mercury/BUSCO 등): 실행 스크립트가 없고 로그 파일(`orthofinder.log` 등)만 남아있는 경우가 많음 — 정확한 실행 커맨드는 확인 필요(2절 참조)

### Step 9. 유전자/반복서열 주석
- **소프트웨어(버전):** Cannabis: RepeatModeler 2.0.6 + RepeatMasker 4.1.9 → BRAKER3 3.0.8(de novo 유전자 예측, RNA-seq+protein evidence 사용) / Pichia: LiftOff(참조 기반 주석 liftover, GS115 및 CBS7435 기준) / Zoysia: RepeatModeler+RepeatMasker(버전 확인 필요, Cannabis와 동일 계열 툴) → Trimmomatic+prinseq+Hisat2(RNA-seq 전처리) → BRAKER3(실행 스크립트/컨테이너 이미지 부재, 메모 파일만 존재) + GeneMarkS-T(ISOseq 보정) + TSEBRA(통합)
- **핵심 파라미터:** Cannabis: RepeatModeler `-LTRStruct`, BRAKER3 `--threads 125`; Pichia: LiftOff `-p 150`(CBS7435 대상은 `-exclude_partial -copies` 추가); Zoysia: RepeatModeler `-LTRStruct -threads 100`, RepeatMasker `-xsmall -gff -pa 100`, Trimmomatic `ILLUMINACLIP:...:2:30:10:2:keepBothReads LEADING:10 TRAILING:10 MINLEN:50 -threads 100`, Hisat2 `-p 200`, TSEBRA `long_reads.cfg`
- **Input:** Cannabis: scaffolding 이전 contig-level masked genome + mRNA-seq + 판게놈 단백질 세트 / Pichia: scaffold(chr 단위, rename 前) fasta + 참조 gff / Zoysia: `Zma_FINAL.20chr.nogaps.sorted.fasta`(Step6 최종본) + Zma-R1/S1 mRNA-seq + ISOseq reads
- **Output:** Cannabis: braker.gtf/aa/codingseq(신규 유전자 모델) / Pichia: liftover gff + unmapped_features 목록(참조 유전자 좌표 이전) / Zoysia: `Zma_FINAL.20chr.nogaps.sorted.fasta.masked` → `braker.gtf` → `gmst.global.gtf`(GeneMarkS-T) → `tsebra.gtf`(TSEBRA 통합) → `Zma_braker3_longest.pep`/`Zma_braker3.rename.longest_isoform.gtf(3)`(최종 유전자 모델)
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
  - Zoysia (반복서열, 06.Repeat_element_annotation/module_RepeatModeler_Masker.sh):
    ```bash
    BuildDatabase -name Zma.db Zma_FINAL.20chr.nogaps.sorted.fasta
    RepeatModeler -database Zma.db -LTRStruct Zma_FINAL.20chr.nogaps.sorted.fasta -threads 100
    RepeatMasker -xsmall -gff -pa 100 -lib RM_*/consensi.fa.classified Zma_FINAL.20chr.nogaps.sorted.fasta
    ```
  - Zoysia (RNA-seq 전처리, 07.Gene_annotation/Trimmomatic/trimmomatic.sh + Hisat2/hisat2.sh):
    ```bash
    TrimmomaticPE -threads 100 -summary <sample>.stat <sample>_1.fastq <sample>_2.fastq <sample>_paired_1.fq <sample>_unpaired_1.fq <sample>_paired_2.fq <sample>_unpaired_2.fq ILLUMINACLIP:Trimmomatic_adapters-PE.fa:2:30:10:2:keepBothReads LEADING:10 TRAILING:10 MINLEN:50
    hisat2 -x Zma_FINAL.20chr.nogaps.sorted.fasta.masked.genome -p 200 -1 <sample>_good_paired_1.fq -2 <sample>_good_paired_2.fq -S <sample>.hisat2.sam
    ```
  - Zoysia (BRAKER3/GeneMarkS-T/TSEBRA): 확인 필요 — BRAKER3는 실행 스크립트/컨테이너 이미지 자체가 레포에 없고(`Braker_pipeline_HS.txt` 메모만 존재), GeneMarkS-T/TSEBRA 세부 커맨드도 본 문서 작성 시점에는 확인되지 않음(2절 참조)

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
  - 2026-07-13(문서 작성 기준일) 현재 아직 초기 단계 프로젝트로, Step1(genome survey)·Step2(assembly)·Step3(organelle 분리)·Step7(QUAST QC)까지만 완료되었고 scaffolding, 비교유전체, 유전자주석 단계는 아직 수행되지 않음 — Cannabis/Pichia/Zoysia와 달리 "완결된" 파이프라인이 아니라 "진행 중" 파이프라인임.
  - 게놈 QC에 BUSCO가 사용되지 않고 QUAST만 사용됨 — 나머지 3개 프로젝트와의 뚜렷한 차이(사유는 확인 필요).
  - 동일 종을 서로 다른 시퀀싱 공급사 조합(JSLink 단독 vs JSLink+EYEONCELL 통합)으로 두 번 독립적으로 조립하고 QUAST로 비교하는 "데이터 소스 조합 검증" 구조 — 다른 3개 프로젝트에는 없는 패턴.
  - hifiasm 조립 자체는 organelle 제거 이전(raw 전체 reads) 상태로 1회만 수행되고, 이후 별도의 사후 필터링(minimap2+seqkit, coverage threshold 0.90)으로 organelle contig를 제거함 — Cannabis는 Oatk 분리 후 조립을 재수행(2회 조립)하는 방식인 반면, Abies는 조립을 1회만 수행하고 필터링만 사후에 적용하는 점이 다름.
  - `run_pipeline.sh`/`run_remaining.sh` 마스터 오케스트레이션 스크립트가 참조하는 디렉토리명(`Organelle_combined`, `hifiasm_combined`)과 현재 실제 디렉토리명(`Organelle_JSLink_EYEONCELL`, `hifiasm_JSLink_EYEONCELL`)이 rename으로 인해 불일치. 최초 실행(2026-06-08)에서 OatK(통합본) 병렬 스텝이 실패하여 hifiasm(통합본)만 완료되었고, OatK는 이후 별도 세션(`run_remaining.sh`, 2026-06-10~11)에서 재실행되어 성공함.
  - FastK 버전은 `-h`/`--version` 플래그를 지원하지 않아 확인 불가. Ratatosk류 하이브리드 보정은 이 프로젝트에서 사용되지 않음(PacBio HiFi 단일 플랫폼만 사용).
- (해당 시) 프로젝트 고유 확장 분석:
  - JSLink 단독 vs JSLink+EYEONCELL 통합 조립본 간 QUAST 비교(`quast_compare_{hap1,hap2,p}`) — 추가 시퀀싱 투입 효과 검증용.
  - hifiasm 그래프 기반 OatK organelle 추출(-G 모드)과 standalone OatK 결과 간 비교(`03_compare_organelle.sh`, QUAST 기반) — organelle 분리 방법론 자체의 교차검증(실행 완료 여부는 산출물로 확인되나 실행 주체/시점 기록은 확인 필요).
  - 향후 계획: ONT ULK 데이터 및 제주도 취약/건전 10샘플씩 PacBio Kinnex Isoseq 데이터 생산 예정 — 향후 전사체 기반 유전자 주석(Zoysia의 Braker3 경로와 유사한 방식)으로 이어질 것으로 예상되나 현재는 미착수.

### Zoysia_macrostachya_dnDNAseq
- 연구목적/샘플 설계: *Zoysia macrostachya*(들잔디속 근연종)의 de novo 전장유전체 해독. ONT+Illumina 하이브리드 long-read 조립 → Omni-C(Hi-C) 기반 de novo 스캐폴딩 → 수동 curation(Juicebox, 2회) → 갭 클로징(YagCloser) → 폴리싱 → 반복서열/유전자/기능 주석 → 비교유전체(Pan-genome, MCMCTree 분화시기, CAFE5 유전자가족 진화, PSG 양성선택) → 자체 발현정량 및 공개 SRA 염 스트레스 전사체 재분석까지 수행하는, 4개 프로젝트 중 가장 확장된 파이프라인. 20개 염색체 수준 최종 게놈(`Zma_FINAL.20chr.nogaps.sorted.fasta`)을 NCBI에 제출 완료.
- 공통 파이프라인과의 차이점:
  - 참조 유전체 기반 RagTag scaffolding이 아니라, Purge_dups(중복 제거) → Omni-C(Hi-C) 리드 기반 HapHiC 스캐폴딩 → Juicebox 수동 curation(2회, 1차는 초기 스캐폴드 검토·2차는 재정렬/SV 검토로 20개 염색체 확정)이라는 de novo 스캐폴딩 경로를 사용 — Cannabis/Pichia와 근본적으로 다른 접근. RagTag는 Zoysia에서도 사용되나 이는 자기 게놈이 아니라 타 Zoysia속 게놈과의 비교유전체용 재정렬 목적(1절 Step 8 참조).
  - 스캐폴딩 이후 YagCloser로 실제 갭 클로징을 수행하여 최종 게놈에 반영 — Pichia도 Yagcloser를 보유하나 X-33에 한해 분석 전용(gap-filling 커맨드 비활성)으로 그친 반면, Zoysia는 실제 갭을 메워 최종 게놈 생성 경로에 포함시킴.
  - 게놈 QC에 LAI(LTR Assembly Index)와 Merqury가 추가로 사용됨 — Cannabis/Pichia에는 없는 QC 지표.
  - 비교유전체 규모가 압도적으로 큼: OrthoFinder(12종/5-Zoysia 오솔로그 클러스터링) → MCMCTree(PAML, 분화시기 추정) → CAFE5(유전자 가족 확장/축소) → PSG(PAML branch-site 양성선택 유전자, 7개 계통 반복)로 이어지는 다단계 파이프라인 — 다른 3개 프로젝트에는 없음.
  - Pan-genome_analysis 디렉토리 하위에 RagTag/SyRI/plotsr/CHROMEISTER/MCScan/CandiSSR/CNV_analysis/Purge_dups/Mercury/BUSCO 등 다수의 병렬 서브분석이 있어(파일 수 13만 6천여 개), Cannabis/Pichia의 단일 synteny 비교보다 훨씬 포괄적.
  - 자체 조직 발현정량(Hisat2+StringTie+TMM 정규화) 및 공개 SRA 염 스트레스 전사체 재분석(Trimmomatic→Prinseq→Hisat2→edgeR DEG)까지 수행 — 다른 3개 프로젝트에서는 RNA-seq이 유전자예측용 hint 생성에만 쓰이는 반면, Zoysia는 발현량 정량 및 차등발현유전자(DEG) 분석까지 확장.
  - BRAKER3/EnTAP/OrthoFinder/MCMCTree/CAFE5 등 핵심 스텝 다수가 실행 스크립트/wrapper 없이 로그·설정 파일(ctl, .ini, .log)만 남아있어 재현성이 낮음 — 확인 필요 항목이 다른 3개 프로젝트보다 현저히 많음.
  - 06.Repeat_element_annotation의 RepeatModeler 실행은 "before_curation"(1차 manual curation 이전 게놈, 2024-10-26)에 대해서만 확인되고, 최종 20chr.sorted 버전에 대한 재실행 디렉토리는 발견되지 않음(날짜-스크립트 불일치) — 확인 필요.
  - Purge_dups(01.Assembly)의 input이 정확히 어느 NextPolish 산출물 계보(Ratatosk 경유 vs 비-Ratatosk 경유, 파일명이 동일하여 구분 불가)인지 확인 필요.
- (해당 시) 프로젝트 고유 확장 분석:
  - Pan-genome_analysis: RagTag/SyRI/plotsr(구조변이), CHROMEISTER/MCScan(synteny), CandiSSR(SSR 마커), CNV_analysis(카피수변이), Purge_dups/Mercury/BUSCO(QC 반복) 등을 통한 근연 Zoysia종(Z. japonica, Z. matrella, Z. pacifica, Z. sinica 등) 다중 비교.
  - MCMCTree(PAML)를 이용한 종간 분화시기 추정.
  - CAFE5를 이용한 유전자 가족 확장/축소(유전자 가족 진화) 분석.
  - PSG(양성선택 유전자) 분석: translatorX → pal2nal → Gblocks → codeml branch-site model → LRT 유의성 검정, 계통(Zma/Zj/Zmt/Zp/Zs/Zoysia/ZsZma)별 반복 수행 → `Zoysia_PSGs_Orthogroups.tsv`.
  - 자체 조직 발현정량(11.Transcriptome, gene/transcript count matrix + TMM 정규화 TPM) 및 공개 SRA 기반 염 스트레스 전사체 재분석(13.SRA_Salt_Transcriptome, edgeR 기반 DEG 분석).
  - 텔로미어 영역 분석(04.Telomere_region, quarTeT/soft-clip 리드 기반 telomere read 추출·패칭) — 최종 게놈 계보와의 직접 연결점은 확인 필요(독립 QC 서브분석으로 추정).

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Cannabis_sativa_cv.ChungSam_dnDNAseq | PacBio HiFi 단일 플랫폼, organelle(Oatk) 분리 후 hifiasm으로 haplotype-resolved(hap1/hap2) diploid 조립, 판게놈 유래 참조(AH3Ma/AH3Mb)로 RagTag scaffolding, RepeatModeler/RepeatMasker+BRAKER3로 de novo 유전자 예측(scaffolding 이전 contig 계열 입력) | hifiasm 0.25.0-r726, RagTag v2.1.0, BUSCO v5.8.2, QUAST v5.2.0, RepeatModeler 2.0.6/RepeatMasker 4.1.9, BRAKER 3.0.8, sourmash 4.9.4(그 외 FastK/GenomeScope2/minimap2/samtools/syri/plotsr/CHROMEISTER/BLAST+ 버전은 확인 필요) | Sourmash+chromeister 판게놈(98샘플) 비교, plotsr/syri 구조변이 시각화, THCAS/CBDAS/CBGAS 카나비노이드 유전자 분석(BRAKER3 단백질 tblastn 비교 + 독립적 targeted gene assembly 파이프라인) |
| Pichia_Pastoris_dnDNAseq_HanHwa | ONT+Illumina 하이브리드(Ratatosk 보정), haploid 단일 genome 조립(NextDenovo+NextPolish 채택, trycycler multi-assembler 경로는 reconcile 실패로 폐기), NCBI GS115 단일 참조로 RagTag scaffolding, LiftOff로 참조 기반 annotation liftover(de novo 예측 아님) | NextDenovo v2.5.2, NextPolish v1.4.1, BUSCO v5.8.2, QUAST v5.2.0, mummer 4.0.1, minimap2 v2.24-r1122(RagTag 내부)(RagTag 자체 버전은 확인 필요) | Merqury k-mer QV/completeness 평가, mummer 기반 공개 균주 synteny 비교, Yagcloser gap 분석(X-33만, 분석 전용), Genome_editing gRNA BLAST 매핑(BG-10만) |
| Abies_koreana_dnDNAseq | PacBio Revio HiFi(JSLink+EYEONCELL 2개사 8cell), hifiasm으로 raw reads 전체를 1회 조립 후 사후 organelle contig 필터링(Oatk+minimap2/seqkit), scaffolding·비교유전체·유전자주석은 아직 미수행(진행 중 프로젝트), BUSCO 없이 QUAST만으로 QC | hifiasm v0.25.0-r726, Oatk v1.0, minimap2 2.30-r1287, seqkit v2.3.0, QUAST 5.3.0(fb88221c)(FastK 버전은 `-h` 미지원으로 확인 불가) | JSLink 단독 vs JSLink+EYEONCELL 통합 조립 비교(QUAST), hifiasm 그래프 기반 vs standalone OatK organelle 추출 방법론 교차검증 |
| Zoysia_macrostachya_dnDNAseq | ONT+Illumina 하이브리드(Ratatosk), NextDenovo+NextPolish+ntLink 조립, Purge_dups+Omni-C(Hi-C) 기반 HapHiC de novo 스캐폴딩+Juicebox 수동 curation(2회)+YagCloser 갭클로징으로 20염색체 확정, RepeatModeler/RepeatMasker+BRAKER3+TSEBRA 유전자예측, LAI/Merqury 추가 QC | Ratatosk/NextDenovo/NextPolish/ntLink(버전 확인 필요), HapHiC(버전 확인 필요), juicer_tools 1.9.9, yagcloser(버전 확인 필요), RepeatModeler/RepeatMasker(버전 확인 필요), BRAKER3(버전/컨테이너 확인 필요), LTR_retriever, Merqury | OrthoFinder(12종/5-Zoysia)→MCMCTree(분화시기)→CAFE5(유전자가족진화)→PSG(PAML branch-site 양성선택), Pan-genome_analysis(RagTag/SyRI/CHROMEISTER/MCScan/CandiSSR/CNV 등 근연종 다중 비교), 자체 발현정량(TMM 정규화) + 공개 SRA 염스트레스 전사체 재분석(edgeR DEG), 텔로미어 영역 분석 |
