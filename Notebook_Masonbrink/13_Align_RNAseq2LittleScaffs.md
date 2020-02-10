## Align RNA-seq to small scaffolds.

Looking to remove redundant small contigs that lack unique rnaseq read mapping.

```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq

#While the little scaffolds were masked prior to mapping, all RNA-seq reads were mapped to just the little scaffolds to ensure completeness of the coding regime

hisat2-build LittleScaffs.fasta.masked maskedLittleScaffs

#align.sh
################################################################################
sh runFeatureCounts.sh  lane3/Elk-kidney_S25_L003_R1_001.fastq lane3/Elk-kidney_S25_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane3/Elk-lung_S26_L003_R1_001.fastq lane3/Elk-lung_S26_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane3/Elk-Mes-LN_S24_L003_R1_001.fastq lane3/Elk-Mes-LN_S24_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane3/Elk-muscle_S21_L003_R1_001.fastq lane3/Elk-muscle_S21_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane3/ElkpscapLN_S22_L003_R1_001.fastq lane3/ElkpscapLN_S22_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane3/Elk-spleen_S23_L003_R1_001.fastq lane3/Elk-spleen_S23_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane3/Undetermined_S0_L003_R1_001.fastq lane3/Undetermined_S0_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane4/Elk-kidney_S25_L004_R1_001.fastq lane4/Elk-kidney_S25_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane4/Elk-lung_S26_L004_R1_001.fastq lane4/Elk-lung_S26_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane4/Elk-Mes-LN_S24_L004_R1_001.fastq lane4/Elk-Mes-LN_S24_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane4/Elk-muscle_S21_L004_R1_001.fastq lane4/Elk-muscle_S21_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane4/ElkpscapLN_S22_L004_R1_001.fastq lane4/ElkpscapLN_S22_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane4/Elk-spleen_S23_L004_R1_001.fastq lane4/Elk-spleen_S23_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
sh runFeatureCounts.sh  lane4/Undetermined_S0_L004_R1_001.fastq lane4/Undetermined_S0_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq maskedLittleScaffs
################################################################################

#runFeatureCounts.sh
################################################################################
#runFeatureCounts.sh
#note: I premade the database so it would not be made again and again
###############################################################################
#!/bin/bash

#Working model.  Things to consider: the script needs different settings for 2 analyses. dna-seq (hisat2), as this is set up for RNA.  free node settings vs regular node settings (procand memory).
#You must provide the following. Note variable DBDIR does not need a "/" at the end.
# sh runFeatureCounts.sh 28 sequence_1.fastq sequence_2.fastq /work/GIF/remkv6/files genome.fa genome.gff3


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

#module load GIF2/java/1.8.0_25-b17
#module load picard
#java -jar /software/7/apps/picard/64/2.9.2/picard.jar CollectAlignmentSummaryMetrics  REFERENCE_SEQUENCE=${DBDIR}/${GENOME} INPUT=${R1_FQ%.*}_sorted.bam OUTPUT=${R1_FQ%.*}.bam_alignment.stats

#module load subread
#featureCounts -T ${PROC} -p -t gene -g ID -a ${GFF} -o ${R1_FQ%.*}_counts_genes.txt ${R1_FQ%.*}_sorted.bam

################################################################################

#merged the bam files, since this is only for determining unique read mapping
module load samtools
samtools merge -@ 40 AllReads_sorted.bam *sorted.bam &

#Create a dummy gff that represents the whole contigs
module load miniconda
source activate bioawk
bioawk -c fastx '{print $name,"rnaseq","gene","1",length($seq),".",".",".","ID="$name}' LittleScaffs.fasta |tr " " "\t" >genome.gff

#run featurecounts to identify unique mapping rnaseq reads
module load subread
featureCounts -T 40 -p -t gene -g ID -a genome.gff -o AllReads_GeneReadCounts.txt AllReads_sorted.bam

```
