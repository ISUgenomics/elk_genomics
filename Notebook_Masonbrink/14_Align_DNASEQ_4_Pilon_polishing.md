# Now that we have a final genome, we need to do some polishing with dnaseq and long reads


### set up file structure
```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd

ln -s ../../illumina/8_9145_01_8Bull-OSMn_HG73W_1014.tar
ln -s ../../illumina/8_9137_01_2450-OS-Mn_HF557_1013.tar
ln -s ../../illumina/7_9144_01_8Bull-OSMn_HG73W_1014.tar
ln -s ../../illumina/7_9136_01_2458-0S-Mn_HF557_1013.tar
ln -s ../../illumina/6_9143_01_8Bull-OSMn_HG73W_1014.tar
ln -s ../../illumina/6_9135_01_WY10-Pine_HF557_1013.tar
ln -s ../../illumina/5_9142_01_8Bull-OSMn_HG73W_1014.tar
ln -s ../../illumina/5_9134_01_WY8-Pine_HF557_1013.tar
ln -s ../../illumina/4_9141_01_2763-OS-Mn_HG73W_1014.tar
ln -s ../../illumina/4_9133_01_WY7-Pine_HF557_1013.tar
ln -s ../../illumina/3_9140_01_2758-OS-Mn_HG73W_1014.tar
ln -s ../../illumina/3_9132_01_WY3-Pine_HF557_1013.tar
ln -s ../../illumina/2_9139_01_2486-OS-Mn_HG73W_1014.tar
ln -s../../illumina/2_9131_01_WY2-Pine_HF557_1013.tar
ln -s ../../illumina/1_9138_01_2510-Os-Mn_HH7YG_1052.tar
ln -s ../../illumina/1_9130_01_WY1-Pine_HF557_1013.tar

tar -xvf 1_9130_01_WY1-Pine_HF557_1013.tar
tar -xvf 1_9138_01_2510-Os-Mn_HH7YG_1052.tar
tar -xvf 2_9131_01_WY2-Pine_HF557_1013.tar
tar -xvf 2_9139_01_2486-OS-Mn_HG73W_1014.tar
tar -xvf 3_9132_01_WY3-Pine_HF557_1013.tar
tar -xvf 3_9140_01_2758-OS-Mn_HG73W_1014.tar
tar -xvf 4_9133_01_WY7-Pine_HF557_1013.tar
tar -xvf 4_9141_01_2763-OS-Mn_HG73W_1014.tar
tar -xvf 5_9134_01_WY8-Pine_HF557_1013.tar
tar -xvf 5_9142_01_8Bull-OSMn_HG73W_1014.tar
tar -xvf 6_9135_01_WY10-Pine_HF557_1013.tar
tar -xvf 6_9143_01_8Bull-OSMn_HG73W_1014.tar
tar -xvf 7_9136_01_2458-0S-Mn_HF557_1013.tar
tar -xvf 7_9144_01_8Bull-OSMn_HG73W_1014.tar
tar -xvf 8_9137_01_2450-OS-Mn_HF557_1013.tar
tar -xvf 8_9145_01_8Bull-OSMn_HG73W_1014.tar



```
### Set up alignment for DNA-seq
```
 paste <(ls -1 *R1* ) <(ls -1 *R2*) |while read line; do echo "sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa "$line;done >align.sh

sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa 2450-OS-Mn_S8_L008_R1_001.fastq.bz2 2450-OS-Mn_S8_L008_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa 2458-0S-Mn_S7_L007_R1_001.fastq.bz2 2458-0S-Mn_S7_L007_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa 2463-OS-Mn_S3_L004_R1_001.fastq.bz2 2463-OS-Mn_S3_L004_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa 2486-OS-Mn_S1_L002_R1_001.fastq.bz2 2486-OS-Mn_S1_L002_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa 2510-Os-Mn_S8_L006_R1_001.fastq.bz2 2510-Os-Mn_S8_L006_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa 2758-OS-Mn_S2_L003_R1_001.fastq.bz2 2758-OS-Mn_S2_L003_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa 8Bull-OSMn_S4_L005_R1_001.fastq.bz2 8Bull-OSMn_S4_L005_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa 8Bull-OSMn_S4_L006_R1_001.fastq.bz2 8Bull-OSMn_S4_L006_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa 8Bull-OSMn_S4_L007_R1_001.fastq.bz2 8Bull-OSMn_S4_L007_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa 8Bull-OSMn_S4_L008_R1_001.fastq.bz2 8Bull-OSMn_S4_L008_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa WY10-Pine_S6_L006_R1_001.fastq.bz2 WY10-Pine_S6_L006_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa WY1-Pine_S1_L001_R1_001.fastq.bz2 WY1-Pine_S1_L001_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa WY2-Pine_S2_L002_R1_001.fastq.bz2 WY2-Pine_S2_L002_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa WY3-Pine_S3_L003_R1_001.fastq.bz2 WY3-Pine_S3_L003_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa WY7-Pine_S4_L004_R1_001.fastq.bz2 WY7-Pine_S4_L004_R2_001.fastq.bz2
sh runPilon.sh stupid /home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd FinalGenome.fa WY8-Pine_S5_L005_R1_001.fastq.bz2 WY8-Pine_S5_L005_R2_001.fastq.bz2

################################################################################################################################################

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

module load hisat2
#hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p 40 -x ${GENOME%.*} -1 $R1_FQ -2 $R2_FQ -S ${R1_FQ%.*}.sam
module load samtools
samtools view --threads 40 -b -o ${GENOME%.*}.${R1_FQ%.*}.bam ${GENOME%.*}.${R1_FQ%.*}.sam
samtools sort -m 3G -o ${GENOME%.*}.${R1_FQ%.*}_sorted.bam -T Round3PilonPB_temp --threads 40 ${GENOME%.*}.${R1_FQ%.*}.bam
samtools index ${GENOME%.*}.${R1_FQ%.*}_sorted.bam
#samtools merge -@ 40 -m TEMP AllReads.bam *_sorted.bam

#module load minimap2
#minimap2 -L -ax map-pb ${GENOME} ${PBReadsFq}  >${GENOME%.*}.${PBReadsFq%.*}.sam

#module load samtools
#samtools view --threads 40 -b -o ${GENOME%.*}.${PBReadsFq%.*}.bam ${GENOME%.*}.${PBReadsFq%.*}.sam
#samtools sort -m 3G -o ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam -T Round3PilonPB_temp --threads 40 ${GENOME%.*}.${PBReadsFq%.*}.bam
#samtools index ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam

#module load pilon

#java -Xmx120g -Djava.io.tmpdir=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/07_blobtools2juiced/TEMP -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome ${GENOME}  --frags AllReads.bam --output ${GENOME%.*}.Pilon --outdir ${DIR} --changes --fix all --threads 40 --mingap 0 --chunksize 100000

################################################################################################################################################
```


### Alignment stats
```
89.85% overall alignment rate
89.19% overall alignment rate
88.98% overall alignment rate
89.16% overall alignment rate
88.95% overall alignment rate
89.47% overall alignment rate
89.43% overall alignment rate
89.77% overall alignment rate
89.63% overall alignment rate
89.68% overall alignment rate
90.53% overall alignment rate
89.97% overall alignment rate
90.36% overall alignment rate
90.47% overall alignment rate
```

### Fix the conversion and sort
```
#the script failed to run past alignment stage, as it was naming the temp files the same for each script.  
#This is what should have been run

for f in *sam; do mkdir ${f%.*}_temp;done

for f in *sam; do echo "module load samtools; samtools view --threads 20 -b -o "${f%.*}".bam " $f "; samtools sort -m 3G -o "${f%.*}"_sorted.bam -T ${f%.*}_temp --threads 20 "${f%.*}".bam; samtools index "${f%.*}".bam";done >BamNSort.sh

#Sams got converted to bams, so just need sorting.
#Coordinate sort and index
for f in *sam; do echo "module load samtools;samtools sort -m 3G -o "${f%.*}"_sorted.bam -T ${f%.*}_temp --threads 20 "${f%.*}".bam; samtools index "${f%.*}".bam";done >SortNIndex.sh

```

### Need to cap the number of reads aligning to regions of the genome to prevent memory issues.  
```
Found an easy way to cap coverage with http://lindenb.github.io/jvarkit/Biostar154220.html

#just a test on one bam
module load miniconda; source activate my_root; module load java/11.0.2;module load samtools|java -Djava.io.tmpdir=TEMP/ -Xmx120g -jar ../11_PolishGenomeWithCCS/jvarkit/dist/sortsamrefname.jar 8Bull-OSMn_S4_L007_R1_001.fastq.bam |java -jar jvarkit/dist/biostar154220.jar  -n 100  |samtools sort -o 8Bull-OSMn_S4_L007_R1_001.fastq_test.bam -

#adapt to all bam files -- max of 30x coverage for each bam
ls *bai |sed 's/\.bai//g'|while read line; do echo "module load miniconda; source activate my_root; module load java/11.0.2;module load samtools;java -Djava.io.tmpdir="${line%.*}"_TEMP/ -Xmx120g -jar ../11_PolishGenomeWithCCS/jvarkit/dist/sortsamrefname.jar --samoutputformat BAM "$line" |java -jar jvarkit/dist/biostar154220.jar  --samoutputformat BAM -n 30  |samtools sort -o "${line%.*}".reduced.bam - ";done >ReduceBams.sh
```


### test pilon with a part of the reads
```
ls -lrth *bai |awk '{print $9}' |sed 's/\.bai//g' |while read line;do echo "--frags "$line; done|tr "\n" " " |sed 's/$/\n/g'|awk '{print "java -Xmx120g -Djava.io.tmpdir=TEMP2  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.fa "$0" --unpaired  ../11_PolishGenomeWithCCS/FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output FinalGenomePilon.fa --outdir . --changes --fix all --threads 40 --chunksize 100000" }' >PilonTest.sh
#out of heap space, as expected

#Another test
ls -1 *reduced* |while read line;do echo "--frags "$line; done|tr "\n" " " |sed 's/$/\n/g'|awk '{print "java -Xmx1020g -Djava.io.tmpdir=TEMP  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome OurElkMitochondria.fasta "$0" --unpaired  ../11_PolishGenomeWithCCS/FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output FinalGenomePilon.fa --outdir . --changes --fix all --threads 40 --chunksize 100000" }' |less

#grabs all the bams
ls -1 *reduced* |while read line;do echo "--frags "$line; done|tr "\n" " " |sed 's/$/\n/g'|awk '{print "java -Xmx120g -Djava.io.tmpdir=TEMP  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.fa "$0" --unpaired  FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output FinalGenomePilon.fa --outdir . --changes --fix all --threads 40 --chunksize 100000" }' |less
```




#### Split scaffolds for pilon ####
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd/01_SplitForPilon

ln -s ../FinalGenome.fa
fasta-splitter.pl --n-parts 36 FinalGenome.fa
unlink FinalGenome.fa

for f in ../*reduced.bam; do ln -s $f;done
for f in ../*reduced.bam.bai; do ln -s $f;done
ln -s ../../11_PolishGenomeWithCCS/FinalGenome.AllCCSReads_sorted.bam
ln -s ../../11_PolishGenomeWithCCS/FinalGenome.AllCCSReads_sorted.bam.bai


for f in *bam; do mkdir ${f%.*}_temp; done
```

### Create pilon run scripts for each scaffold with reduced bam coverage of PE reads and all CCS reads
```
#use this to make script below
ls -1 *reduced.bam |while read line;do echo "--frags "$line; done|tr "\n" " " |sed 's/$/\n/g'|awk '{print "java -Xmx120g -Djava.io.tmpdir=TEMP  -jar /
software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.fa "$0" --unpaired  FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output FinalGenomePilon.f
a --outdir . --changes --fix all --threads 40 --chunksize 100000" }' |less


#run scripts for pilon
for f in *fa; do echo "module load pilon;java -Xmx120g -Djava.io.tmpdir="${f%.*}"_temp  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome "$f" --frags 2450-OS-Mn_S8_L008_R1_001.fastq_sorted.reduced.bam --frags 2458-0S-Mn_S7_L007_R1_001.fastq_sorted.reduced.bam --frags 2463-OS-Mn_S3_L004_R1_001.fastq_sorted.reduced.bam --frags 2486-OS-Mn_S1_L002_R1_001.fastq_sorted.reduced.bam --frags 2510-Os-Mn_S8_L006_R1_001.fastq_sorted.reduced.bam --frags 2758-OS-Mn_S2_L003_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L005_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L006_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L007_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L008_R1_001.fastq_sorted.reduced.bam --frags WY10-Pine_S6_L006_R1_001.fastq_sorted.reduced.bam --frags WY1-Pine_S1_L001_R1_001.fastq_sorted.reduced.bam --frags WY2-Pine_S2_L002_R1_001.fastq_sorted.reduced.bam --frags WY3-Pine_S3_L003_R1_001.fastq_sorted.reduced.bam --frags WY7-Pine_S4_L004_R1_001.fastq_sorted.reduced.bam --frags WY8-Pine_S5_L005_R1_001.fastq_sorted.reduced.bam  --unpaired  FinalGenome.AllCCSReads_sorted.bam --output "${f%.*}"Pilon.fa --outdir . --changes --fix all --threads 40 --chunksize 100000";done >Pilon.sh


submitted each one to a single short node


######################Results##################################3


#Number of changes that each fasta file/pseudomolecule encountered
[rick.masonbrink@sn-cn-18-1 01_SplitForPilon]$ for f in *.changes; do wc -l $f;done
2369 FinalGenome.part-01Pilon.fa.changes
29155 FinalGenome.part-03Pilon.fa.changes
26624 FinalGenome.part-04Pilon.fa.changes
23091 FinalGenome.part-11Pilon.fa.changes
21005 FinalGenome.part-12Pilon.fa.changes
27877 FinalGenome.part-13Pilon.fa.changes
19986 FinalGenome.part-14Pilon.fa.changes
16510 FinalGenome.part-15Pilon.fa.changes
15092 FinalGenome.part-17Pilon.fa.changes
22802 FinalGenome.part-18Pilon.fa.changes
15887 FinalGenome.part-19Pilon.fa.changes
15517 FinalGenome.part-20Pilon.fa.changes
28002 FinalGenome.part-21Pilon.fa.changes
14326 FinalGenome.part-22Pilon.fa.changes
13529 FinalGenome.part-23Pilon.fa.changes
11402 FinalGenome.part-24Pilon.fa.changes
24279 FinalGenome.part-25Pilon.fa.changes
13618 FinalGenome.part-26Pilon.fa.changes
14757 FinalGenome.part-27Pilon.fa.changes
12722 FinalGenome.part-29Pilon.fa.changes
13944 FinalGenome.part-30Pilon.fa.changes
17720 FinalGenome.part-31Pilon.fa.changes
22096 FinalGenome.part-32Pilon.fa.changes
11842 FinalGenome.part-33Pilon.fa.changes
9409 FinalGenome.part-34Pilon.fa.changes
36215 FinalGenome.part-35Pilon.fa.changes
2406 FinalGenome.part-36Pilon.fa.changes


#how much sequence has changed?
module load miniconda
source activate bioawk

less *Pilon.fa.changes |awk '{print $4}' |sed 's/ //g' |sed 's/^/>fastaName\n/g' |sed 's/\.//g' |bioawk -c fastx '{print length($seq)}' |summary.sh
Total:  17,891,227
Count:  650,590
Mean:   27
Median: 1
Min:    0
Max:    5,789


#How many of the nucleotides were confirmed by pilon?  -- still missing part_16

#total confirmed
cat Pilon_*.o* |grep "Confirmed" |awk '{print $2}' |summary.sh
Total:  2,491,997,436
Count:  25,577
Mean:   97,431
Median: 99,426
Min:    466
Max:    99,991

#Total Genome
 cat Pilon_*.o* |grep "Confirmed" |awk '{print $4}' |summary.sh
Total:  2,512,474,760
Count:  25,577
Mean:   98,231
Median: 99,919
Min:    1,004
Max:    100,000


So far, 99.18% confirmed
```

### Hi-C_scaffold_16 timed out after 48hrs on a high mem node splitting to 100kb chunks, only using one paired end bam and CCS to polish
The process still hung -- see masking below

```
#split and run pilon on the 100kb parts
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd/01_SplitForPilon/01_SplitPilon15


ln -s ../FinalGenome.part-16.fa
for f in ../*reduced.bam; do ln -s $f;done
for f in ../*reduced.bam.bai; do ln -s $f;done
ln -s ../../../11_PolishGenomeWithCCS/FinalGenome.AllCCSReads_sorted.bam
ln -s ../../../11_PolishGenomeWithCCS/FinalGenome.AllCCSReads_sorted.bam.bai


for f in *bam; do mkdir ${f%.*}_temp; done


#use this to make script below
ls -1 *reduced.bam |while read line;do echo "--frags "$line; done|tr "\n" " " |sed 's/$/\n/g'|awk '{print "java -Xmx120g -Djava.io.tmpdir=TEMP  -jar /
software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.fa "$0" --unpaired  FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output FinalGenomePilon.f
a --outdir . --changes --fix all --threads 40 --chunksize 100000" }' |less

java -Xmx120g -Djava.io.tmpdir=TEMP  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.fa --frags 2450-OS-Mn_S8_L008_R1_001.fastq_sorted.reduced.bam --frags 2458-0S-Mn_S7_L007_R1_001.fastq_sorted.reduced.bam --frags 2463-OS-Mn_S3_L004_R1_001.fastq_sorted.reduced.bam --frags 2486-OS-Mn_S1_L002_R1_001.fastq_sorted.reduced.bam --frags 2510-Os-Mn_S8_L006_R1_001.fastq_sorted.reduced.bam --frags 2758-OS-Mn_S2_L003_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L005_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L006_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L007_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L008_R1_001.fastq_sorted.reduced.bam --frags WY10-Pine_S6_L006_R1_001.fastq_sorted.reduced.bam --frags WY1-Pine_S1_L001_R1_001.fastq_sorted.reduced.bam --frags WY2-Pine_S2_L002_R1_001.fastq_sorted.reduced.bam --frags WY3-Pine_S3_L003_R1_001.fastq_sorted.reduced.bam --frags WY7-Pine_S4_L004_R1_001.fastq_sorted.reduced.bam --frags WY8-Pine_S5_L005_R1_001.fastq_sorted.reduced.bam  --unpaired  FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output FinalGenomePilon.fa --outdir . --changes --fix all --threads 40 --chunksize 100000

#run scripts for pilon
for f in *ChromosomeShredder; do echo "module load pilon;java -Xmx20g -Djava.io.tmpdir="$f"_temp  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome "$f" --frags 2450-OS-Mn_S8_L008_R1_001.fastq_sorted.reduced.bam --frags 2458-0S-Mn_S7_L007_R1_001.fastq_sorted.reduced.bam --frags 2463-OS-Mn_S3_L004_R1_001.fastq_sorted.reduced.bam --frags 2486-OS-Mn_S1_L002_R1_001.fastq_sorted.reduced.bam --frags 2510-Os-Mn_S8_L006_R1_001.fastq_sorted.reduced.bam --frags 2758-OS-Mn_S2_L003_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L005_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L006_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L007_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L008_R1_001.fastq_sorted.reduced.bam --frags WY10-Pine_S6_L006_R1_001.fastq_sorted.reduced.bam --frags WY1-Pine_S1_L001_R1_001.fastq_sorted.reduced.bam --frags WY2-Pine_S2_L002_R1_001.fastq_sorted.reduced.bam --frags WY3-Pine_S3_L003_R1_001.fastq_sorted.reduced.bam --frags WY7-Pine_S4_L004_R1_001.fastq_sorted.reduced.bam --frags WY8-Pine_S5_L005_R1_001.fastq_sorted.reduced.bam  --unpaired  FinalGenome.AllCCSReads_sorted.bam --output "${f%.*}"Pilon.fa --outdir . --changes --fix all --threads 1 --chunksize 100000";done >Pilon.sh

#generates a target list
for f in {1..746}; do echo "HiC_scaffold_18: "$f; done |awk '{print $1,(($2*100000)-99999)"-"($2*100000)}' |sed 's/ //g' |while read line; do echo "module load pilon;java -Xmx120g -Djava.io.tmpdir="$line"_temp  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.part-16.fa --frags 2450-OS-Mn_S8_L008_R1_001.fastq_sorted.reduced.bam --frags 2458-0S-Mn_S7_L007_R1_001.fastq_sorted.reduced.bam --frags 2463-OS-Mn_S3_L004_R1_001.fastq_sorted.reduced.bam --frags 2486-OS-Mn_S1_L002_R1_001.fastq_sorted.reduced.bam --frags 2510-Os-Mn_S8_L006_R1_001.fastq_sorted.reduced.bam --frags 2758-OS-Mn_S2_L003_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L005_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L006_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L007_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L008_R1_001.fastq_sorted.reduced.bam --frags WY10-Pine_S6_L006_R1_001.fastq_sorted.reduced.bam --frags WY1-Pine_S1_L001_R1_001.fastq_sorted.reduced.bam --frags WY2-Pine_S2_L002_R1_001.fastq_sorted.reduced.bam --frags WY3-Pine_S3_L003_R1_001.fastq_sorted.reduced.bam --frags WY7-Pine_S4_L004_R1_001.fastq_sorted.reduced.bam --frags WY8-Pine_S5_L005_R1_001.fastq_sorted.reduced.bam  --unpaired  FinalGenome.AllCCSReads_sorted.bam --output "${line%.*}"Pilon.fa --outdir . --changes --fix all --threads 40 --targets "$line" --chunksize 50000";done >Pilon.sh


for f in {1..746}; do echo "HiC_scaffold_18: "$f; done |awk '{print $1,(($2*100000)-99999)"-"($2*100000)}' |sed 's/ //g' |while read line; do mkdir $line"_temp";done

for f in ../*reduced.bam; do ln -s $f;done
for f in ../*reduced.bam.bai; do ln -s $f;done
ln -s ../../../11_PolishGenomeWithCCS/FinalGenome.AllCCSReads_sorted.bam
ln -s ../../../11_PolishGenomeWithCCS/FinalGenome.AllCCSReads_sorted.bam.bai
```

### Try to mask the ribosomal array so Pilon will move
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/14_PilonPairedEnd/01_SplitForPilon/02_MaskPart16

less consensi.fa.classified |grep "Satellite" |sed 's/>//g' |sed 's/#/\t/g' |awk '{print $1}' |while read line; do grep -w $line FinalGenome.part-16.fa.out.gff; done |sort -k 4,5n |awk '{print $1,$4,$5}' |tr " " "\t" |bedtools merge -d 100000|awk '($3-$2)>2000' |bedtools complement -i - -g FinalGenome.part-16.length|bedtools makewindows -b - -w 50000 |sed 's/\t/:/1' |sed 's/\t/-/1' |tr "\n" "," |awk '{print $0}' |while read line; do echo "module load pilon;java -Xmx120g -Djava.io.tmpdir=TEMP  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.part-16.fa  --frags WY7-Pine_S4_L004_R1_001.fastq_sorted.reduced.bam --frags WY8-Pine_S5_L005_R1_001.fastq_sorted.reduced.bam  --unpaired  FinalGenome.AllCCSReads_sorted.bam --output FinalGenome.part-16.Pilon.fa --outdir . --changes --fix all --threads 40 --targets "$line" --chunksize 50000";done >Pilon.sh

#skip all repeats, not just satellites.  have to use high mem node if using all reads, so cut to just one.
less FinalGenome.part-16.fa.out.gff |sort -k 4,5n |awk '{print $1,$4,$5}' |tr " " "\t" |bedtools merge -d 1000|awk '($3-$2)>5000' |bedtools complement -i - -g FinalGenome.part-16.length|bedtools makewindows -b - -w 50000 |sed 's/\t/:/1' |sed 's/\t/-/1' |tr "\n" "," |awk '{print $0}' |while read line; do echo "module load pilon;java -Xmx120g -Djava.io.tmpdir=TEMP  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.part-16.fa --frags 2458-0S-Mn_S7_L007_R1_001.fastq_sorted.reduced.bam --unpaired  FinalGenome.AllCCSReads_sorted.bam --output FinalGenome.part-16.Pilon.fa --outdir . --changes --fix all --threads 40 --targets "$line" --chunksize 50000";done >Pilon.sh


#This is how much of the chromosome is going to be polished.
less FinalGenome.part-16.fa.out.gff |sort -k 4,5n |awk '{print $1,$4,$5}' |tr " " "\t" |bedtools merge -d 1000|awk '($3-$2)>5000' |bedtools complement -i - -g FinalGenome.part-16.length| awk '{print $3-$2}' |summary.sh
Total:  49,328,884
Count:  2,619
Mean:   18,835
Median: 9,104
Min:    1,001
Max:    1,888,303

#How much of the sequence could be confirmed with pilon
bioawk -c fastx '{print $name, length($seq)}' FinalGenome.part-16.Pilon.fa.fasta |awk '{print $2}' |~/common_scripts/summary.sh
Total:  49,260,024
Count:  2,954
Mean:   16,675
Median: 10,450
Min:    373
Max:    50,287


Finished.  REMEMBER NEXT TIME THAT PILON IS 1 COORDINATE BASED, NOT ZERO BASED.
#changed the zero at the start coordinate, and pilon finished.
#so pilon does not concatenate the targeted sections, and leaves them as separate scaffolds.  need to merge in the proper order

#unpolished sections
less Pilon_0.o943881 |grep "changes to" |awk '{print $2}' |sed 's/:/\t/g' |sed 's/-/\t/g' |bedtools complement -i - -g FinalGenome.part-16.length |sed 's/\t/:/1' |sed 's/\t/-/1' |awk 'NR>1' |while read line; do samtools faidx FinalGenome.part-16.fa $line; done |tr "\n" "\t" |sed 's/>/\n>/g' |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' >fastaList2

#polished sections
paste <(less Pilon_0.o943881 |grep "changes to" |awk '{print ">"$2}' ) <(less FinalGenome.part-16.Pilon.fa.fasta |tr "\n" "\t" |sed 's/>/\n>/g' |sed 's/\t//g' |sed 's/pilon/pilon\t/g' |awk 'NR>1' ) |cut -f 1,3 >fastaList1
cat fastaList1 fastaList2|sed 's/:/\t/1' |sed 's/-/\t/1' |sort -k2,3n |awk 'NR>1' |cut -f 4 |tr "\n" "#" |sed 's/#//g' |cat <(echo ">HiC_scaffold_18_pilon") - >Fixed.FinalGenome.part-16.Pilon.fa.fasta

cat ../*Pilon.fa.fasta Fixed.FinalGenome.part-16.Pilon.fa.fasta >FinalGenomePilon.fa

---------------- Information for assembly 'FinalGenomePilon.fa' ----------------


                                         Number of scaffolds        493
                                     Total size of scaffolds 2529349837
                                            Longest scaffold  146388637
                                           Shortest scaffold       1004
                                 Number of scaffolds > 1K nt        493 100.0%
                                Number of scaffolds > 10K nt        177  35.9%
                               Number of scaffolds > 100K nt         35   7.1%
                                 Number of scaffolds > 1M nt         35   7.1%
                                Number of scaffolds > 10M nt         34   6.9%
                                          Mean scaffold size    5130527
                                        Median scaffold size       6008
                                         N50 scaffold length   77654944
                                          L50 scaffold count         13
                                         n90 scaffold length   51438166
                                          L90 scaffold count         29
                                                 scaffold %A      29.18
                                                 scaffold %C      20.79
                                                 scaffold %G      20.79
                                                 scaffold %T      29.16
                                                 scaffold %N       0.08
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      99.8%
              Percentage of assembly in unscaffolded contigs       0.2%
                      Average number of contigs per scaffold       10.8
Average length of break (>25 Ns) between contigs in scaffold        392

                                           Number of contigs       5306
                              Number of contigs in scaffolds       4873
                          Number of contigs not in scaffolds        433
                                       Total size of contigs 2527461845
                                              Longest contig    8486465
                                             Shortest contig        180
                                   Number of contigs > 1K nt       5298  99.8%
                                  Number of contigs > 10K nt       4718  88.9%
                                 Number of contigs > 100K nt       3427  64.6%
                                   Number of contigs > 1M nt        777  14.6%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size     476340
                                          Median contig size     209904
                                           N50 contig length    1121308
                                            L50 contig count        634
                                           n90 contig length     274935
                                            L90 contig count       2339
                                                   contig %A      29.20
                                                   contig %C      20.81
                                                   contig %G      20.81
                                                   contig %T      29.18
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0

```
