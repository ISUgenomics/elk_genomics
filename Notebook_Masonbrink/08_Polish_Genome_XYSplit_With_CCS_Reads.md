# Polish XYSplit Elk genome


### Gather files /structure
```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/11_PolishGenomeWithCCS

ln -s ../04_JuicerElk/01_3DNA/3d-dna/01_PostJBAssembly/FinalAssembly/3d-dna/FinalAssemblyFastaWithY.fasta

#Concatenate CCs reads
cat ../10_CCS_polish/*ccs.fa >AllCCSReads.fa

#What is the coverage of the genome? .5x
bioawk -c fastx '{print $name,length($seq)}' AllCCSReads.fa |awk '{print $2}' |summary.sh
Total:  1,819,991,677
Count:  149,340
Mean:   12,186
Median: 13,197
Min:    117
Max:    48,127
```


#Align the ccs to the genome and run pilon

```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/11_PolishGenomeWithCCS

sh runPilon.sh AllCCSReads.fa /home/rick.masonbrink/elk_bison_genomics/Masonbrink/11_PolishGenomeWithCCS FinalGenome.fa

###############################################################################
#runPilon.sh
###############################################################################
#!/bin/bash

#You must provide the following. Note variable DBDIR does not need a "/" at the end.
# sh runPilon.sh LongReads.fastq /work/GIF/remkv6/files genome.fa ShortReadsR1.fq ShortReadsR2.fq


PBReadsFq="$1"
DIR="$2"
GENOME="$3"
R1_FQ="$4"
R2_FQ="$5"

#module load hisat2
#hisat2-build ${GENOME} ${GENOME%.*}
#hisat2 -p 16 -x ${GENOME%.*} -1 $R1_FQ -2 $R2_FQ -S ${R1_FQ%.*}.sam
#module load samtools
#samtools view --threads 16 -b -o ${GENOME%.*}.${R1_FQ%.*}.bam ${GENOME%.*}.${R1_FQ%.*}.sam
#samtools sort -m 7G -o ${GENOME%.*}.${R1_FQ%.*}_sorted.bam -T Round3PilonPB_temp --threads 16 ${GENOME%.*}.${R1_FQ%.*}.bam
#samtools index ${GENOME%.*}.${R1_FQ%.*}_sorted.bam

module load minimap2
minimap2 -L -ax map-pb ${GENOME} ${PBReadsFq}  >${GENOME%.*}.${PBReadsFq%.*}.sam

module load samtools
samtools view --threads 40 -b -o ${GENOME%.*}.${PBReadsFq%.*}.bam ${GENOME%.*}.${PBReadsFq%.*}.sam
samtools sort -m 3G -o ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam -T Round3PilonPB_temp --threads 40 ${GENOME%.*}.${PBReadsFq%.*}.bam
samtools index ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam

module load pilon
java -Xmx1024g -Djava.io.tmpdir=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/07_blobtools2juiced/TEMP -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome ${GENOME}  --unpaired ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam --output ${GENOME%.*}.Pilon --outdir ${DIR} --changes --fix all --threads 40 --mingap 0 --chunksize 100000


###############################################################################

## This was not sufficient to allow for changes to be made.  Most contigs had a single subread aligning and a minimum depth of 5 is needed to make changes.

```


### Restart mapping with subreads
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/11_PolishGenomeWithCCS

sh runPilon.sh AllSubreads.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/11_PolishGenomeWithCCS FinalGenome.fa

```
