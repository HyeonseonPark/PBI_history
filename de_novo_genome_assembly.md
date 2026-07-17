# De novo Genome Assembly (Cannabis_sativa_cv.ChungSam_dnDNAseq, Pichia_Pastoris_dnDNAseq_HanHwa)

## 1. 공통 파이프라인

두 프로젝트 모두 "원시 long-read(+보조 short-read) → k-mer 기반 genome profiling → de novo assembly →
참조 기반 scaffolding → genome QC → 비교유전체(synteny/dot-plot) 분석 → 유전자 주석" 이라는
동일한 상위 단계 구조를 따른다. 다만 실제 사용 시퀀싱 플랫폼, 조립 도구, 배수성(ploidy)
처리 방식, 주석 방법(de novo 예측 vs liftover)은 프로젝트별로 상이하다 (2절 참조).

```
[rawdata: long-read (+ short-read)]
        │
        ▼
[Step 1. Genome survey / k-mer profiling: FastK + GenomeScope(2)]
        │  (게놈 크기·이형접합률 추정치 → 이후 assembly 파라미터로 재사용)
        ▼
[Step 2. De novo assembly (프로젝트별 도구/전략 상이 — 2절 참조)]
        │
        ▼
[Step 3. 참조 기반 scaffolding: RagTag]
        │
        ├──▶ [Step 4. Genome QC: BUSCO + QUAST]
        │
        ├──▶ [Step 5. 비교유전체 분석 (synteny/dot-plot, 도구는 프로젝트별 상이)]
        │
        └──▶ [Step 6. 유전자/반복서열 주석 (de novo 예측 또는 참조 기반 liftover, 프로젝트별 상이)]
                        │
                        ▼
                [프로젝트 고유 확장 분석 / 최종 산출물]
```

| 스텝 | 소프트웨어(버전) | 핵심 파라미터 | Input | Output |
|---|---|---|---|---|
| 1. Genome survey (k-mer profiling) | FastK, GenomeScope2 (Cannabis: 버전 확인 필요; Pichia: 사용 확인되었으나 세부 버전은 재현 우선순위 낮음으로 처리됨) | k=31(Cannabis, HiFi reads 대상) / k=21(Pichia, Illumina reads 대상) | 원시 long-read(Cannabis, PacBio HiFi) 또는 QC된 Illumina reads(Pichia) | k-mer histogram, genome size/heterozygosity 추정치(GenomeScope summary) |
| 2. De novo assembly | Cannabis: hifiasm 0.25.0-r726 / Pichia: NextDenovo v2.5.2 + NextPolish v1.4.1 | Cannabis: `-t 100`(스레드) / Pichia: `input_type=corrected, read_type=ont, genome_size=<GenomeScope 추정치>`(NextDenovo), `task=best, sgs_options=-max_depth 100 -bwa`(NextPolish, Illumina 폴리싱만 사용) | Cannabis: organelle 제거 후 nuclear-only HiFi reads / Pichia: ONT reads를 Illumina로 하이브리드 교정(Ratatosk)한 corrected long reads + Illumina reads(폴리싱용) | Cannabis: haplotype-resolved contig(hap1/hap2 × Female/Male) / Pichia: 폴리싱된 단일(haploid) genome fasta |
| 3. 참조 기반 scaffolding | RagTag v2.1.0(Cannabis) / RagTag(Pichia, 버전 확인 필요) — 내부적으로 둘 다 minimap2 사용(Cannabis asm5 프리셋, Pichia v2.24-r1122 asm5 프리셋) | `-t 200`(Cannabis) | 조립 결과(contig-level) + 참조 유전체(Cannabis: 판게놈 유래 AH3Ma/AH3Mb; Pichia: NCBI GS115 GCF_000027005.1) | scaffold fasta(.agp, .stats 포함), 이후 chromosome/pseudochromosome 수준으로 재정리 |
| 4. Genome QC | BUSCO v5.8.2(양 프로젝트 공통), QUAST v5.2.0(양 프로젝트 공통) | BUSCO: `-m genome`, lineage(Cannabis: rosales_odb12/rosids_odb10, Pichia: ascomycota_odb12/saccharomycetes_odb12); QUAST: 참조 유전체 지정, 기본 옵션 | scaffold(chr 단위) fasta | BUSCO short_summary(.txt/.json), QUAST report(.html/.tsv/.pdf) |
| 5. 비교유전체 분석(synteny/dot-plot) | Cannabis: chromeister + minimap2/syri/plotsr(버전 확인 필요) / Pichia: mummer 4.0.1(nucmer+mummerplot) | Cannabis: minimap2 `-ax asm5 --eqx`, syri `-F B` / Pichia: nucmer `--maxmatch -l 20 -g 500 -c 100` | scaffold(chr 단위) fasta 및 비교 대상 유전체(Cannabis: 판게놈 98샘플; Pichia: 공개 Pichia 균주 게놈) | dot-plot 이미지(.png), 구조변이 vcf/summary(Cannabis) |
| 6. 유전자/반복서열 주석 | Cannabis: RepeatModeler 2.0.6 + RepeatMasker 4.1.9 → BRAKER3 3.0.8(de novo 유전자 예측, RNA-seq+protein evidence 사용) / Pichia: LiftOff(참조 기반 주석 liftover, GS115 및 CBS7435 기준) | Cannabis: RepeatModeler `-LTRStruct`, BRAKER3 `--threads 125`; Pichia: LiftOff `-p 150`(CBS7435 대상은 `-exclude_partial -copies` 추가) | Cannabis: scaffolding 이전 contig-level masked genome + mRNA-seq + 판게놈 단백질 세트 / Pichia: scaffold(chr 단위, rename 前) fasta + 참조 gff | Cannabis: braker.gtf/aa/codingseq(신규 유전자 모델) / Pichia: liftover gff + unmapped_features 목록(참조 유전자 좌표 이전) |

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

## 3. 요약 비교표

| 프로젝트 | 핵심 차이점 요약 | 소프트웨어 버전 차이 | 고유 확장 분석 |
|---|---|---|---|
| Cannabis_sativa_cv.ChungSam_dnDNAseq | PacBio HiFi 단일 플랫폼, organelle(Oatk) 분리 후 hifiasm으로 haplotype-resolved(hap1/hap2) diploid 조립, 판게놈 유래 참조(AH3Ma/AH3Mb)로 RagTag scaffolding, RepeatModeler/RepeatMasker+BRAKER3로 de novo 유전자 예측(scaffolding 이전 contig 계열 입력) | hifiasm 0.25.0-r726, RagTag v2.1.0, BUSCO v5.8.2, QUAST v5.2.0, RepeatModeler 2.0.6/RepeatMasker 4.1.9, BRAKER 3.0.8, sourmash 4.9.4(그 외 FastK/GenomeScope2/minimap2/samtools/syri/plotsr/CHROMEISTER/BLAST+ 버전은 확인 필요) | Sourmash+chromeister 판게놈(98샘플) 비교, plotsr/syri 구조변이 시각화, THCAS/CBDAS/CBGAS 카나비노이드 유전자 분석(BRAKER3 단백질 tblastn 비교 + 독립적 targeted gene assembly 파이프라인) |
| Pichia_Pastoris_dnDNAseq_HanHwa | ONT+Illumina 하이브리드(Ratatosk 보정), haploid 단일 genome 조립(NextDenovo+NextPolish 채택, trycycler multi-assembler 경로는 reconcile 실패로 폐기), NCBI GS115 단일 참조로 RagTag scaffolding, LiftOff로 참조 기반 annotation liftover(de novo 예측 아님) | NextDenovo v2.5.2, NextPolish v1.4.1, BUSCO v5.8.2, QUAST v5.2.0, mummer 4.0.1, minimap2 v2.24-r1122(RagTag 내부)(RagTag 자체 버전은 확인 필요) | Merqury k-mer QV/completeness 평가, mummer 기반 공개 균주 synteny 비교, Yagcloser gap 분석(X-33만, 분석 전용), Genome_editing gRNA BLAST 매핑(BG-10만) |
