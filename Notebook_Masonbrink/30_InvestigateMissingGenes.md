#  Need to characterize the genes that are presumed missing from the elk genome


### Andrew's files
```
#Location of files that Andrew created
#/work/gif/severin/Olsen/Elk2020/07_immuneGenes


#Proposed gene deletions -- locations not accurate
|AY644518_TRGJ1| 01| 102107 |
|KT723008_IGHJ2|01| 509129 |
|AC172685_IGHA| 01| 655709|
|IMGT000049_TRDC|01| 3301039|
|D16120_TRGJ2|02| 338 |
|AY149283_IGHJ1| 02 | 783 |
|AY2277782_TRAJ31| 02| 1725 |
|AY644517_TRG| 02 |167787 |
|NW_001494075_IGHJ1|03| 39638 |


```

### Identify Regions in the cattle genome that these genes are located
```
#/work/gif/remkv6/Olsen/Elk/09_MapRnaseq2Cattle

vi Missing_From_Elk.list
################################################################################
AY644518_TRGJ1
KT723008_IGHJ2
AC172685_IGHA
IMGT000049_TRDC
D16120_TRGJ2
AY149283_IGHJ1
AY2277782_TRAJ31
AY644517_TRGC4
NW_001494075_IGHJ1
#################################################################################


ln -s /work/gif/severin/Olsen/Elk2020/07_immuneGenes/bosImmune_2_Elk/bestBlast_output.txt ElkBBH.out
ln -s /work/gif/severin/Olsen/Elk2020/07_immuneGenes/bosImmune_2_Bos/bestBlast_output.txt CowBBH.out


#grep these gene names from the blast output and create a bed file

##################################################################################
IMGT000049_TRDC*01_Bos_taurus_Hereford_F_EX2_g,3231015..3231088_75_nt_1_+1_-1___25_AA_25+0=25_____      NC_037337.1     100.000 25      0       0       1       25      22353349        22353275        1.77e-08
        51.6    25      103308737       0/-2    NC_037337.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 10, ARS-UCD1.2, whole genome shotgun sequence  100
AY644517_TRGC4*02_Bos_taurus_F_EX2A_g,167787..167845_60_nt_1_+1_-1___20_AA_20+0=20_____ NC_037331.1     100.000 20      0       0       1       20      82742556        82742615        7.07e-05        41.2    20      120000601       0/3     NC_037331.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 4, ARS-UCD1.2, whole genome shotgun sequence   100
AY644518_TRGJ1-1*01_Bos_taurus_F_J-REGION_102107..102166_60_nt_3_______19_AA_19+0=19_____       NC_037331.1     100.000 19      0       0       1       19      50155487        50155431        4.40e-05        41.6    19      120000601       0/-3    NC_037331.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 4, ARS-UCD1.2, whole genome shotgun sequence   100
D16120_TRGJ2-1*02_Bos_taurus_(F)_J-REGION_338..396_59_nt_3_______19_AA_19+0=19_partial_in_3'___ NC_037331.1     94.737  19      1       0       1       19      50155487        50155431        8.07e-04        38.1    19      120000601       0/-3    NC_037331.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 4, ARS-UCD1.2, whole genome shotgun sequence   100
NW_001494075_IGHJ1-1*03_Bos_taurus_Hereford_ORF_J-REGION_39812..39874_63_nt_3_______20_AA_20+0=20___rev-compl_  NC_037348.1     100.000 20      0       0       1       20      208561  208502  1.60e-05        42.7    20      69862954        0/-1    NC_037348.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 21, ARS-UCD1.2, whole genome shotgun sequence  100
AY149283_IGHJ1-2*02_Bos_taurus_ORF_J-REGION_783..838_56_nt_2_______18_AA_18+0=18_____   NC_037348.1     94.444  18      1       0       1       18      208381  208328  1.49e-04        40.0    18      69862954
        0/-1    NC_037348.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 21, ARS-UCD1.2, whole genome shotgun sequence  100
NW_001494075_IGHJ1-2*03_Bos_taurus_Hereford_ORF_J-REGION_39638..39693_56_nt_2_______18_AA_18+0=18___rev-compl_  NC_037348.1     100.000 18      0       0       1       18      208381  208328  5.76e-05        41.2    18      69862954        0/-1    NC_037348.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 21, ARS-UCD1.2, whole genome shotgun sequence  100
KT723008_IGHJ2-1*01_Bos_taurus_Holstein_ORF_J-REGION_508948..509010_63_nt_3_______20_AA_20+0=20_____    NC_037348.1     100.000 20      0       0       1       20      208561  208502  1.60e-05        42.7    20      69862954        0/-1    NC_037348.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 21, ARS-UCD1.2, whole genome shotgun sequence  100
KT723008_IGHJ2-2*01_Bos_taurus_Holstein_ORF_J-REGION_509129..509184_56_nt_2_______18_AA_18+0=18_____    NC_037348.1     100.000 18      0       0       1       18      208381  208328  5.76e-05        41.2    18      69862954        0/-1    NC_037348.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 21, ARS-UCD1.2, whole genome shotgun sequence  100

##################################################################################




#Created bed file
cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1<3' |grep -f Missing_From_Elk.list - |awk '{print $2}'  |sed 's/\*/\\*/g' |grep -f - CowBBH.out |awk '{if($9>$10) {print $2,$10,$9,$1} else {print $2,$9,$10,$1} }' |sed 's/\*/\t/g' |sort -k1,1V -k2,3n |tr " " "\t" |cut -f 1,2,3,4 |less

################################################################################
NC_037331.1     50155431        50155487        AY644518_TRGJ1-1
NC_037331.1     50155431        50155487        D16120_TRGJ2-1
NC_037331.1     82742556        82742615        AY644517_TRGC4
NC_037337.1     22353275        22353349        IMGT000049_TRDC
NC_037348.1     208328  208381  AY149283_IGHJ1-2
NC_037348.1     208328  208381  KT723008_IGHJ2-2
NC_037348.1     208328  208381  NW_001494075_IGHJ1-2
NC_037348.1     208502  208561  KT723008_IGHJ2-1
NC_037348.1     208502  208561  NW_001494075_IGHJ1-1

################################################################################


#Merged bed file
cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1<3' |grep -f Missing_From_Elk.list - |awk '{print $2}'  |sed 's/\*/\\*/g' |grep -f - CowBBH.out |awk '{if($9>$10) {print $2,$10,$9,$1} else {print $2,$9,$10,$1} }' |sed 's/\*/\t/g' |sort -k1,1V -k2,3n |tr " " "\t" |cut -f 1,2,3,4 |bedtools merge -i - | less

#############################################################################
NC_037331.1     50155431        50155487  AY644518_TRGJ1-1, D16120_TRGJ2-1
NC_037331.1     82742556        82742615  AY644517_TRGC4
NC_037337.1     22353275        22353349  IMGT000049_TRDC
NC_037348.1     208328  208381  AY149283_IGHJ1-2, KT723008_IGHJ2-2,NW_001494075_IGHJ1-2
NC_037348.1     208502  208561  KT723008_IGHJ2-1, NW_001494075_IGHJ1-1



```


### Bedtools intersect with the genic gff of cattle
```
which mRNA/genes are these proteins associated with?
cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1<3' |grep -f Missing_From_Elk.list - |awk '{print $2}'  |sed 's/\*/\\*/g' |grep -f - CowBBH.out |awk '{if($9>$10) {print $2,$10,$9,$1} else {print $2,$9,$10,$1} }' |sed 's/\*/\t/g' |sort -k1,1V -k2,3n |tr " " "\t" |cut -f 1,2,3,4 |bedtools merge -i - | bedtools intersect -wo -a - -b ../06_Synteny/02_BosTaurus/GCF_002263795.1_ARS-UCD1.2_genomic.gff |awk '$6=="gene"' |less


######################################################################################################################
NC_037331.1     50155431        50155487        NC_037331.1     BestRefSeq%2CGnomon     gene    50141812        50162835        .       -       .       ID=gene-TARP;Dbxref=GeneID:100335800;Name=TARP;description=TCR gamma alternate reading frame protein;gbkey=Gene;gene=TARP;gene_biotype=protein_coding;gene_synonym=TCR,TRGC1,TRGJ1-2    56
NC_037331.1     82742556        82742615        NC_037331.1     Gnomon  gene    82720193        82746245        .       +       .       ID=gene-TRGC4;Dbxref=GeneID:509614;Name=TRGC4;gbkey=Gene;gene=TRGC4;gene_biotype=C_region       59
NC_037337.1     22353275        22353349        NC_037337.1     Gnomon  gene    22350419        22359327        .       -       .       ID=gene-LOC101908015;Dbxref=GeneID:101908015;Name=LOC101908015;gbkey=Gene;gene=LOC101908015;gene_biotype=C_region       74
NC_037337.1     22353275        22353349        NC_037337.1     Gnomon  gene    22263114        23332869        .       -       .       ID=gene-LOC100336282;Dbxref=GeneID:100336282;Name=LOC100336282;gbkey=Gene;gene=LOC100336282;gene_biotype=protein_coding;partial=true    74
######################################################################################################################

So it looks to be 4 genes, not nine.

These two did not intersect with a gene.
NC_037348.1     208328  208381  AY149283_IGHJ1-2, KT723008_IGHJ2-2,NW_001494075_IGHJ1-2
NC_037348.1     208502  208561  KT723008_IGHJ2-1, NW_001494075_IGHJ1-1

Investigated on https://bovinegenome.elsiklab.missouri.edu/ JBROWSE

Interprets to chromosome 21:208238..208561.  
This is clearly a transcribed locus according to isoseq from, and when viewing among the bam dense files for tissues, only appears to be expressed in nasal mucousa

NC_037348.1     208328  208381  AY149283_IGHJ1-2, KT723008_IGHJ2-2,NW_001494075_IGHJ1-2
NC_037348.1     208502  208561  KT723008_IGHJ2-1, NW_001494075_IGHJ1-1



```

### Extract genes and align
```
#essentialy this is just grabbing blast hits only found in cattle, identifying the region they map to in cattle, finding the corresponding gene in cattle, extracting the genomic sequence
cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1<3' |grep -f Missing_From_Elk.list - |awk '{print $2}'  |sed 's/\*/\\*/g' |grep -f - CowBBH.out |awk '{if($9>$10) {print $2,$10,$9,$1} else {print $2,$9,$10,$1} }' |sed 's/\*/\t/g' |sort -k1,1V -k2,3n |tr " " "\t" |cut -f 1,2,3,4 |bedtools merge -i - | bedtools intersect -wo -a - -b ../06_Synteny/02_BosTaurus/GCF_002263795.1_ARS-UCD1.2_genomic.gff |awk '$6=="gene"' |uniq| awk '{print $4":"$7"-"$8}' |while read line; do samtools faidx GCF_002263795.1_ARS-UCD1.2_genomic.fna $line; done >ExtractedMissingGenes.fasta

#Put it together with the elk genome, to make sure spurious reads do not map to the cattle genes.
cat ExtractedMissingGenes.fasta ../08_RenameAgain/01_Move2Box/CervusCanadensisGenome.fa >MissingGenesAddedToElkGenome.fasta


```

### Transfer modified genome to ceres and run expression
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression

#Create a gff that only recognizes these four geneic scaffolds.
ml miniconda
source activate bioawk
bioawk -c fastx '{print $name,"1",length($seq)}' MissingGenesAddedToElkGenome.fasta |head -n 4 |awk '{print $1,"NA","gene",$2,$3,"+",".",".","ID="$1"gene" }' |tr " " "\t" >MissingGenesAddedToElkGenome.gff

#softlink the rnaseq
Elk-spleen_S23_L003_R2_001.fastq -> ../16_RNAseq/lane3/Elk-spleen_S23_L003_R2_001.fastq
Elk-spleen_S23_L003_R1_001.fastq -> ../16_RNAseq/lane3/Elk-spleen_S23_L003_R1_001.fastq
ElkpscapLN_S22_L003_R2_001.fastq -> ../16_RNAseq/lane3/ElkpscapLN_S22_L003_R2_001.fastq
ElkpscapLN_S22_L003_R1_001.fastq -> ../16_RNAseq/lane3/ElkpscapLN_S22_L003_R1_001.fastq
Elk-muscle_S21_L003_R2_001.fastq -> ../16_RNAseq/lane3/Elk-muscle_S21_L003_R2_001.fastq
Elk-muscle_S21_L003_R1_001.fastq -> ../16_RNAseq/lane3/Elk-muscle_S21_L003_R1_001.fastq
Elk-Mes-LN_S24_L003_R2_001.fastq -> ../16_RNAseq/lane3/Elk-Mes-LN_S24_L003_R2_001.fastq
Elk-Mes-LN_S24_L003_R1_001.fastq -> ../16_RNAseq/lane3/Elk-Mes-LN_S24_L003_R1_001.fastq
Elk-lung_S26_L003_R2_001.fastq -> ../16_RNAseq/lane3/Elk-lung_S26_L003_R2_001.fastq
Elk-lung_S26_L003_R1_001.fastq -> ../16_RNAseq/lane3/Elk-lung_S26_L003_R1_001.fastq
Elk-kidney_S25_L003_R2_001.fastq -> ../16_RNAseq/lane3/Elk-kidney_S25_L003_R2_001.fastq
Elk-kidney_S25_L003_R1_001.fastq -> ../16_RNAseq/lane3/Elk-kidney_S25_L003_R1_001.fastq
Elk-spleen_S23_L004_R2_001.fastq -> ../16_RNAseq/lane4/Elk-spleen_S23_L004_R2_001.fastq
Elk-spleen_S23_L004_R1_001.fastq -> ../16_RNAseq/lane4/Elk-spleen_S23_L004_R1_001.fastq
ElkpscapLN_S22_L004_R2_001.fastq -> ../16_RNAseq/lane4/ElkpscapLN_S22_L004_R2_001.fastq
ElkpscapLN_S22_L004_R1_001.fastq -> ../16_RNAseq/lane4/ElkpscapLN_S22_L004_R1_001.fastq
Elk-muscle_S21_L004_R2_001.fastq -> ../16_RNAseq/lane4/Elk-muscle_S21_L004_R2_001.fastq
Elk-muscle_S21_L004_R1_001.fastq -> ../16_RNAseq/lane4/Elk-muscle_S21_L004_R1_001.fastq
Elk-Mes-LN_S24_L004_R2_001.fastq -> ../16_RNAseq/lane4/Elk-Mes-LN_S24_L004_R2_001.fastq
Elk-Mes-LN_S24_L004_R1_001.fastq -> ../16_RNAseq/lane4/Elk-Mes-LN_S24_L004_R1_001.fastq
Elk-lung_S26_L004_R2_001.fastq -> ../16_RNAseq/lane4/Elk-lung_S26_L004_R2_001.fastq
Elk-lung_S26_L004_R1_001.fastq -> ../16_RNAseq/lane4/Elk-lung_S26_L004_R1_001.fastq
Elk-kidney_S25_L004_R2_001.fastq -> ../16_RNAseq/lane4/Elk-kidney_S25_L004_R2_001.fastq
Elk-kidney_S25_L004_R1_001.fastq -> ../16_RNAseq/lane4/Elk-kidney_S25_L004_R1_001.fastq


#map and count
sh runFeatureCounts.sh Elk-kidney_S25_L003_R1_001.fastq Elk-kidney_S25_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-kidney_S25_L004_R1_001.fastq Elk-kidney_S25_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-lung_S26_L003_R1_001.fastq Elk-lung_S26_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-lung_S26_L004_R1_001.fastq Elk-lung_S26_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-Mes-LN_S24_L003_R1_001.fastq Elk-Mes-LN_S24_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-Mes-LN_S24_L004_R1_001.fastq Elk-Mes-LN_S24_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-muscle_S21_L003_R1_001.fastq Elk-muscle_S21_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-muscle_S21_L004_R1_001.fastq Elk-muscle_S21_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh ElkpscapLN_S22_L003_R1_001.fastq ElkpscapLN_S22_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh ElkpscapLN_S22_L004_R1_001.fastq ElkpscapLN_S22_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-spleen_S23_L003_R1_001.fastq Elk-spleen_S23_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-spleen_S23_L004_R1_001.fastq Elk-spleen_S23_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff


##############################################################################################################################################

#runFeatureCounts.sh
#note: I premade the database so it would not be made again and again
###############################################################################
#!/bin/bash

PROC=40
R1_FQ="$1"
R2_FQ="$2"
DBDIR="$3"
GENOME="$4"
GFF="$5"

#module load trimgalore
#module load py-setuptools/35.0.2-py2-hqrsjff
#trim_galore --paired ${R1_FQ} ${R2_FQ}

module load hisat2
#hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p ${PROC} -x ${GENOME%.*} -1 $R1_FQ -2 $R2_FQ -S ${R1_FQ%.*}.sam

module load samtools
samtools view --threads ${PROC} -b -o ${R1_FQ%.*}.bam ${R1_FQ%.*}.sam
samtools sort -m 3G -o ${R1_FQ%.*}_sorted.bam -T ${R1_FQ%.*}_temp --threads ${PROC} ${R1_FQ%.*}.bam

module load GIF2/java/1.8.0_25-b17
module load picard
java -jar /software/7/apps/picard/64/2.9.2/picard.jar CollectAlignmentSummaryMetrics  REFERENCE_SEQUENCE=${DBDIR}/${GENOME} INPUT=${R1_FQ%.*}_sorted.bam OUTPUT=${R1_FQ%.*}.bam_alignment.stats

module load subread
featureCounts -T ${PROC} -p -t gene -g ID -a ${GFF} -o ${R1_FQ%.*}_counts_genes.txt ${R1_FQ%.*}_sorted.bam
##############################################################################################################################################
```

### Alignment stats
```
for f in *stats; do echo "awk '\$1==\"PAIR\"{print FILENAME,\$7}' "$f" >>RNAseqAlignment.txt";done >stats.sh

Elk-kidney_S25_L003_R1_001.bam_alignment.stats 0.923363
Elk-kidney_S25_L004_R1_001.bam_alignment.stats 0.922779
Elk-Mes-LN_S24_L003_R1_001.bam_alignment.stats 0.907378
Elk-Mes-LN_S24_L004_R1_001.bam_alignment.stats 0.906394
Elk-muscle_S21_L003_R1_001.bam_alignment.stats 0.741857
Elk-muscle_S21_L004_R1_001.bam_alignment.stats 0.741641
ElkpscapLN_S22_L003_R1_001.bam_alignment.stats 0.934488
ElkpscapLN_S22_L004_R1_001.bam_alignment.stats 0.934037
Elk-spleen_S23_L003_R1_001.bam_alignment.stats 0.887191
```

# Expression Results
```

for f in *counts_genes.txt; do echo "<(cut -f 7 "$f" |grep -v \"#\" )" ;done |tr "\n" " " |sed 's/^/paste <(cut -f 1 Elk-spleen_S23_L004_R1_001_counts_genes.txt |grep -v "#" ) /g' |less

paste <(cut -f 1 Elk-spleen_S23_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-kidney_S25_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-kidney_S25_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-lung_S26_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-lung_S26_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-Mes-LN_S24_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-Mes-LN_S24_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-muscle_S21_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-muscle_S21_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 ElkpscapLN_S22_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 ElkpscapLN_S22_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-spleen_S23_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-spleen_S23_L004_R1_001_counts_genes.txt |grep -v "#" ) |less

################################################################################
Geneid  Elk-kidney_S25_L003_R1_001_sorted.bam   Elk-kidney_S25_L004_R1_001_sorted.bam   Elk-lung_S26_L003_R1_001_sorted.bam     Elk-lung_S26_L004_R1_001_sorted.bam     Elk-Mes-LN_S24_L003_R1_001_sorted.bam   Elk-Mes-LN_S24_L004_R1_001_sorted.bam   Elk-muscle_S21_L003_R1_001_sorted.bam   Elk-muscle_S21_L004_R1_001_sorted.bam   ElkpscapLN_S22_L003_R1_001_sorted.bam   ElkpscapLN_S22_L004_R1_001_sorted.bam   Elk-spleen_S23_L003_R1_001_sorted.bam   Elk-spleen_S23_L004_R1_001_sorted.bam
NC_037331.1:50141812-50162835gene       0       1       1       0       3       1       0       0       5       4       8       7
NC_037331.1:82720193-82746245gene       1       3       5       2       2       4       0       0       1       0       3       1
NC_037337.1:22350419-22359327gene       0       0       0       0       0       0       0       0       0       0       0       0
NC_037337.1:22263114-23332869gene       59      51      131     115     232     220     9       19      72      66      55      61
################################################################################

```
Only the bottom gene seems to have a significant number of reads, and likely there is a gene in that region conserved between cattle and elk.  I need finer depth windows.

## Investigate missing genes in the pacbio subreads

```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/33_MissingGenesPacbioAlignment

AllFastq.fastq -> ../06_blobtools/AllFastq.fastq
MissingGenesAddedToElkGenome.fasta -> ../32_MissingGeneExpression/MissingGenesAddedToElkGenome.fasta
MissingGenesAddedToElkGenome.gff -> ../32_MissingGeneExpression/MissingGenesAddedToElkGenome.gff

#asm20 allows for a 5% divergence in sequence, since I am mapping elk reads to cattle genes.
minimap2 -H -a -o MissingGenesAddedToElkGenome.sam -t 40 -x asm20 MissingGenesAddedToElkGenome.fasta AllFastq.fastq;ml samtools; samtools sort -o MissingGenesAddedToElkGenome.bam MissingGenesAddedToElkGenome.sam; ml subread; featureCounts -T 40 -p -t gene -g ID -a MissingGenesAddedToElkGenome.gff  -o MissingGenesAddedToElkGenome_counts_genes.txt MissingGenesAddedToElkGenome.bam

```
### Pacbio alignment rate
```

samtools flagstat MissingGenesAddedToElkGenome.sam
12887228 + 0 in total (QC-passed reads + QC-failed reads)
6977634 + 0 secondary
497574 + 0 supplementary
0 + 0 duplicates
10932108 + 0 mapped (84.83% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

this is ~13% lower than what i got using different settings with minimap, assuming it is because I told it -asm20, allowing for divergence.
```

### Pacbio Mapping Results
```
Geneid  Chr     Start   End     Strand  Length  MissingGenesAddedToElkGenome.sam
NC_037331.1:50141812-50162835gene       NC_037331.1:50141812-50162835   1       21024   +       21024   0
NC_037331.1:82720193-82746245gene       NC_037331.1:82720193-82746245   1       26053   +       26053   0
NC_037337.1:22350419-22359327gene       NC_037337.1:22350419-22359327   1       8909    +       8909    0
NC_037337.1:22263114-23332869gene       NC_037337.1:22263114-23332869   1       1069756 +       1069756 32
```
With the bottom gene being 100kb long, it is tough to say whether the gene is missing biologically.  This fourth gene might be missing from our assembly, but be there biologically.  

##  investigate regions these proteins map to in the cattle genome and 100bp surrounding them.
On Nova now
```
#/work/gif/remkv6/Olsen/Elk/09_MapRnaseq2Cattle

cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1<3' |grep -f Missing_From_Elk.list - |awk '{print $2}'  |sed 's/\*/\\*/g' |grep -f - CowBBH.out |awk '{if($9>$10) {print $2,$10,$9,$1} else {print $2,$9,$10,$1} }' |sed 's/\*/\t/g' |sort -k1,1V -k2,3n |tr " " "\t" |cut -f 1,2,3,4 |bedtools merge -i - |awk '{print $1":"$2-100"-"$3+100}' |while read line; do samtools faidx GCF_002263795.1_ARS-UCD1.2_genomic.fna $line; done >ExtractedMissingRegions.fasta

cat ExtractedMissingRegions.fasta ../08_RenameAgain/01_Move2Box/CervusCanadensisGenome.fa >MissingRegionsAddedToElkGenome.fasta
#Transfer to ceres
##############################################################################################
>NC_037331.1:50155331-50155587
CAAGAATGTGGCAGCCACATCAGGGAGAGCCCAGGTTGTTTCTACAAAGCATCTCTCAGT
ATTTTATGACTACAGTAagaatagaaaaatgatactTACCAGGAGGAGTTACTATGATGT
TAGCTCCTTCTCCAAATATCTTCCTCCAGCCTGAGCTGTCACAATGGTAAAAGCTCCTAC
AATAATCCAGCCTTTTTCTACCTTTATAAAGTGACTAGAAGTCAAAAACTTCATCCTGGC
AGTTTGATCTGAATTCT
>NC_037331.1:82742456-82742715
atatacatatatatgtattatacaataacacatgtataaacacacacaatgTCAAATAAT
GGGCACAAATAATTGACtaatacttttttcatttcttgtagTTGTCAGCTCAATTGTCCC
TACCACTGAGTCTCCCAGTGACTGTTTAAACCATGACAGCAGTAAGTTTTTGTAGATGAT
TGCCTTTATGTCATGCTTtagtcttttatttcctttctaatCTATTTATCTTTTGAACTT
CTCTGGCTAAAAGTACAAAA
>NC_037337.1:22353175-22353449
ttgttttcttttgaaatagaATGATTTGGGGCCCCATCTGACAATGACAGAGTTTGCTTT
AACCTCAATTCTGAAGCTGGCCCTTTGGCTTCAACTAACCTTGGGGCTCTTGGCAGgtca
ctggagcttcagcttttgtGCTGTTTTCATATGCCATCGGTTTTGGAGTTGTTTCTGGAA
CAAACAGCAATGGTCAGGGTTACAATGACAACTGAGAGGCTCAAGGGCTTTCAGGTTTTA
GAACTGGAGGTGACTTTAAGTAAGTAAACTAAGGA
>NC_037348.1:208228-208481
ACCCCGTTGGACCAGGACCCTACCTTGCTCAGGGCACGGTGGGGCTGCACGGACCCCAGG
GAACTGCCTGAGGCCCAGGCAGTGCTGCAAGGGAGACACTCCCCAAGGACACGGTGACCG
GGGTGCGCTGGCCCCAGAGATCCATGTCCCAGCAGCATGGCTTCCattctgttgcttctt
ttccACAAAAACACACCCTACCTCCTTGCTGGCCCTGACTCACATCCCAAGTCCACCCAA
GTCCCACTGTCCCC
>NC_037348.1:208402-208661
cttcttttccACAAAAACACACCCTACCTCCTTGCTGGCCCTGACTCACATCCCAAGTCC
ACCCAAGTCCCACTGTCCCCTGACTTGGACACTAGACACCAGATTCAGCTGAGGAGACGG
TGCCCAGGGCAGCCTGGCCACAGAGATGGAAGTCAGCATAGTCACAGTGGGCTGGCCCCG
TGCTGGGAGCACGGAAACCCACCCAGGGTCCCCAGGGAGCCCCCAGGAGCTGCAAGTCTG
TCTTCCTGAAGAcacttgcc
##############################################################################################

```


##### Moved to ceres
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression/01_Region

#grabbed genome from different location, but is same
cat MissinginElk_CattleRegions.fasta  ../../31_CallSNPS/CervusCanadensisGenome.fa > MissingRegionsAddedToElkGenome.fasta


 hisat2-build MissingRegionsAddedToElkGenome.fasta MissingRegionsAddedToElkGenome

 bioawk -c fastx '{print $name,"1",length($seq)}' MissinginElk_CattleRegions.fasta |head -n 4 |awk '{print $1,"NA","gene",$2,$3,"+",".",".","ID="$1"gene" }' |tr " " "\t" >MissingRegionsAddedToElkGenome.gff

paste <(ls -1 *_R1_* )  <(ls -1 *_R2_* ) |while read line; do echo "sh runFeatureCounts.sh "$line" /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression/01_Region MissingRegionsAddedToElkGenome.fasta MissingRegionsAddedToElkGenome.gff"; done >RunPipeline.sh


Region expression

paste <(cut -f 1 Elk-spleen_S23_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-kidney_S25_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-kidney_S25_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-lung_S26_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-lung_S26_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-Mes-LN_S24_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-Mes-LN_S24_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-muscle_S21_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-muscle_S21_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 ElkpscapLN_S22_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 ElkpscapLN_S22_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-spleen_S23_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-spleen_S23_L004_R1_001_counts_genes.txt |grep -v "#" ) |less

#######################################################################################
Geneid  Elk-kidney_S25_L003_R1_001_sorted.bam   Elk-kidney_S25_L004_R1_001_sorted.bam   Elk-lung_S26_L003_R1_001_sorted.bam     Elk-lung_S26_L004_R1_001_sorted.bam     Elk-Mes-LN_S24_L003_R1_001_sorted.bam   Elk-Mes-LN_S24_L004_R1_001_sorted.bam   Elk-muscle_S21_L003_R1_001_sorted.bam   Elk-muscle_S21_L004_R1_001_sorted.bam   ElkpscapLN_S22_L003_R1_001_sorted.bam   ElkpscapLN_S22_L004_R1_001_sorted.bam   Elk-spleen_S23_L003_R1_001_sorted.bam   Elk-spleen_S23_L004_R1_001_sorted.bam
NC_037331.1:50155331-50155587gene       0       0       0       0       0       0       0       0       0       0       0       0
NC_037331.1:82742456-82742715gene       0       0       0       0       0       0       0       0       0       0       0       0
NC_037337.1:22353175-22353449gene       0       0       0       0       0       0       0       0       0       0       0       0
NC_037348.1:208228-208481gene   0       0       0       0       0       0       0       0       0       0       0       0
#######################################################################################
```

### Genomic structure of where missing genes should be.

```
NC_037331.1:50155331-50155587 --> Bt_Chr3  50155331  50155587
NC_037331.1:82742456-82742715 --> Bt_Chr3 82742456  82742715
NC_037337.1:22353175-22353449 --> Bt_Chr10  22353175  22353449
NC_037348.1:208228-208481 --> Bt_Chr21  208228  208481

vi MissingGeneRegionsCattle.bed
##################################
Bt_Chr3  50155331  50155587
Bt_Chr3 82742456  82742715
Bt_Chr10  22353175  22353449
Bt_Chr21  208228  208481
#############################

less SyntenicRibbons.conf |sort -k 4,5 |awk '{print $4,$5,$6,$1,$2,$3}' |tr " " "\t" |bedtools intersect -wo -a - -b MissingGeneRegionsCattle.bed |less
```
These regions were not found to be in synteny.
