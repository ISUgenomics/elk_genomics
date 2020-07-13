# Need Sequence counts for coverage Calculations

### Pacbio DNA and Illumina DNA -- from minimap2 and hisat2 alignments for pilon
```
FinalGenome.AllCCSReads_sorted.bam: 242770 reads, 82234 filtered, 137964 mapped, 0 proper, 0 stray, Unpaired 100% 13072+/-5312, max 29009
8Bull-OSMn_S4_L006_R1_001.fastq_sorted.reduced.bam: 495708327 reads, 72518901 filtered, 364970856 mapped, 338648670 proper, 2920464 stray, FR 100% 419+/-122, max 1037
2450-OS-Mn_S8_L008_R1_001.fastq_sorted.reduced.bam: 559083608 reads, 87519953 filtered, 402067748 mapped, 371926230 proper, 3524144 stray, FR 100% 391+/-119, max 873
2486-OS-Mn_S1_L002_R1_001.fastq_sorted.reduced.bam: 518187496 reads, 70049044 filtered, 385788749 mapped, 358573832 proper, 2828904 stray, FR 100% 393+/-116, max 811
2510-Os-Mn_S8_L006_R1_001.fastq_sorted.reduced.bam: 545368832 reads, 77651826 filtered, 398885373 mapped, 368147590 proper, 3042110 stray, FR 100% 392+/-115, max 1387
WY2-Pine_S2_L002_R1_001.fastq_sorted.reduced.bam: 545428141 reads, 81259334 filtered, 394059936 mapped, 362692062 proper, 3214010 stray, FR 100% 382+/-118, max 886
8Bull-OSMn_S4_L007_R1_001.fastq_sorted.reduced.bam: 554100989 reads, 88495721 filtered, 398875197 mapped, 370486042 proper, 3487012 stray, FR 100% 397+/-121, max 910
8Bull-OSMn_S4_L008_R1_001.fastq_sorted.reduced.bam: 563674900 reads, 91022079 filtered, 406487642 mapped, 378104030 proper, 3592470 stray, FR 100% 385+/-119, max 900
2758-OS-Mn_S2_L003_R1_001.fastq_sorted.reduced.bam: 547434638 reads, 78039458 filtered, 403149791 mapped, 373317710 proper, 3174874 stray, FR 100% 388+/-117, max 818
2458-0S-Mn_S7_L007_R1_001.fastq_sorted.reduced.bam: 574253897 reads, 88318553 filtered, 410916674 mapped, 378982940 proper, 3445476 stray, FR 100% 388+/-118, max 928
WY1-Pine_S1_L001_R1_001.fastq_sorted.reduced.bam: 566969549 reads, 84789290 filtered, 401595555 mapped, 367892154 proper, 3381686 stray, FR 100% 390+/-119, max 862
WY8-Pine_S5_L005_R1_001.fastq_sorted.reduced.bam: 572953693 reads, 93209447 filtered, 402527807 mapped, 368713328 proper, 3750670 stray, FR 100% 393+/-121, max 964
8Bull-OSMn_S4_L005_R1_001.fastq_sorted.reduced.bam: 582550048 reads, 94801339 filtered, 418630223 mapped, 389176120 proper, 3712806 stray, FR 100% 394+/-121, max 916
WY10-Pine_S6_L006_R1_001.fastq_sorted.reduced.bam: 576931503 reads, 94065151 filtered, 406106721 mapped, 373000486 proper, 3816438 stray, FR 100% 389+/-121, max 928
WY7-Pine_S4_L004_R1_001.fastq_sorted.reduced.bam: 587043778 reads, 99478391 filtered, 409765427 mapped, 376123804 proper, 3952404 stray, FR 100% 393+/-122, max 987
2463-OS-Mn_S3_L004_R1_001.fastq_sorted.reduced.bam: 589807453 reads, 94166068 filtered, 418423520 mapped, 386148712 proper, 3757266 stray, FR 100% 393+/-121, max 897
WY3-Pine_S3_L003_R1_001.fastq_sorted.reduced.bam: 597029323 reads, 96913153 filtered, 418066411 mapped, 383312494 proper, 3894614 stray, FR 100% 381+/-120, max 941


Pacbio CCS mapping % -- (137,964 + 82,234)/242770 = 90.7%

#filtered reads (mapped)
less Pilon_15.o901707|awk 'NR >6 && NR< 23' |awk '{print $4}' |~/common_scripts/summary.sh
#mapped reads
less Pilon_15.o901707|awk 'NR >6 && NR< 23' |awk '{print $6}' |~/common_scripts/summary.sh
#total reads
less Pilon_15.o901707|awk 'NR >6 && NR< 23' |awk '{print $2}' |~/common_scripts/summary.sh

DNA-seq mapping % --(1,392,297,708 + 6,440,317,630) )/8,976,526,175 = 87.26%
```

### Illumina RNA
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/11_AlignSoftmaskedRename

ml samtools; samtools flagstat AllStrandedRNASeq.bam

1825389599 + 0 in total (QC-passed reads + QC-failed reads)
354903545 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
1509460229 + 0 mapped (82.69% : N/A)
1470486054 + 0 paired in sequencing
735243027 + 0 read1
735243027 + 0 read2
1101364886 + 0 properly paired (74.90% : N/A)
1132237898 + 0 with itself and mate mapped
22318786 + 0 singletons (1.52% : N/A)
4989586 + 0 with mate mapped to a different chr
2899793 + 0 with mate mapped to a different chr (mapQ>=5)


#illumina RNA-seq was 82.69% mapped
```
