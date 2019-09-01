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
### Set up alignment
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
#just a test
module load miniconda; source activate my_root; module load java/11.0.2;module load samtools|java -Djava.io.tmpdir=TEMP/ -Xmx120g -jar ../11_PolishGenomeWithCCS/jvarkit/dist/sortsamrefname.jar 8Bull-OSMn_S4_L007_R1_001.fastq.bam |java -jar jvarkit/dist/biostar154220.jar  -n 100  |samtools sort -o 8Bull-OSMn_S4_L007_R1_001.fastq_test.bam -

#adapt to all bam files
ls *bai |sed 's/\.bai//g'|while read line; do echo "module load miniconda; source activate my_root; module load java/11.0.2;module load samtools;java -Djava.io.tmpdir="${line%.*}"_TEMP/ -Xmx120g -jar ../11_PolishGenomeWithCCS/jvarkit/dist/sortsamrefname.jar --samoutputformat BAM "$line" |java -jar jvarkit/dist/biostar154220.jar  --samoutputformat BAM -n 30  |samtools sort -o "${line%.*}".reduced.bam - ";done >ReduceBams.sh
```


### test pilon with a part of the reads
```
ls -lrth *bai |awk '{print $9}' |sed 's/\.bai//g' |while read line;do echo "--frags "$line; done|tr "\n" " " |sed 's/$/\n/g'|awk '{print "java -Xmx120g -Djava.io.tmpdir=TEMP2  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.fa "$0" --unpaired  ../11_PolishGenomeWithCCS/FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output FinalGenomePilon.fa --outdir . --changes --fix all --threads 40 --chunksize 100000" }' >PilonTest.sh
#out of heap space, as expected

#Another test
ls -1 *reduced* |while read line;do echo "--frags "$line; done|tr "\n" " " |sed 's/$/\n/g'|awk '{print "java -Xmx1020g -Djava.io.tmpdir=TEMP  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome OurElkMitochondria.fasta "$0" --unpaired  ../11_PolishGenomeWithCCS/FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output FinalGenomePilon.fa --outdir . --changes --fix all --threads 40 --chunksize 100000" }' |less

ls -1 *reduced* |while read line;do echo "--frags "$line; done|tr "\n" " " |sed 's/$/\n/g'|awk '{print "java -Xmx120g -Djava.io.tmpdir=TEMP  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.fa "$0" --unpaired  FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output FinalGenomePilon.fa --outdir . --changes --fix all --threads 40 --chunksize 100000" }' |less





#### Split scaffolds for pilon ####

ln -s ../FinalGenome.fa
fasta-splitter.pl --n-parts 36 FinalGenome.fa
unlink FinalGenome.fa

for f in ../*reduced.bam; do ln -s $f;done
for f in ../*reduced.bam.bai; do ln -s $f;done
ln -s ../../11_PolishGenomeWithCCS/FinalGenome.AllCCSReads_sorted.bam
ln -s ../../11_PolishGenomeWithCCS/FinalGenome.AllCCSReads_sorted.bam.bai


for f in *bam; do mkdir ${f%.*}_temp; done


#use this to make script below
ls -1 *reduced.bam |while read line;do echo "--frags "$line; done|tr "\n" " " |sed 's/$/\n/g'|awk '{print "java -Xmx120g -Djava.io.tmpdir=TEMP  -jar /
software/7/apps/pilon/1.23/pilon-1.23.jar --genome FinalGenome.fa "$0" --unpaired  FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output FinalGenomePilon.f
a --outdir . --changes --fix all --threads 40 --chunksize 100000" }' |less


#run scripts for pilon
for f in *fa; do echo "module load pilon;java -Xmx120g -Djava.io.tmpdir="${f%.*}_temp  -jar /software/7/apps/pilon/1.23/pilon-1.23.jar --genome "$f" --frags 2450-OS-Mn_S8_L008_R1_001.fastq_sorted.reduced.bam --frags 2458-0S-Mn_S7_L007_R1_001.fastq_sorted.reduced.bam --frags 2463-OS-Mn_S3_L004_R1_001.fastq_sorted.reduced.bam --frags 2486-OS-Mn_S1_L002_R1_001.fastq_sorted.reduced.bam --frags 2510-Os-Mn_S8_L006_R1_001.fastq_sorted.reduced.bam --frags 2758-OS-Mn_S2_L003_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L005_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L006_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L007_R1_001.fastq_sorted.reduced.bam --frags 8Bull-OSMn_S4_L008_R1_001.fastq_sorted.reduced.bam --frags WY10-Pine_S6_L006_R1_001.fastq_sorted.reduced.bam --frags WY1-Pine_S1_L001_R1_001.fastq_sorted.reduced.bam --frags WY2-Pine_S2_L002_R1_001.fastq_sorted.reduced.bam --frags WY3-Pine_S3_L003_R1_001.fastq_sorted.reduced.bam --frags WY7-Pine_S4_L004_R1_001.fastq_sorted.reduced.bam --frags WY8-Pine_S5_L005_R1_001.fastq_sorted.reduced.bam  --unpaired  FinalAssemblyFastaWithY.AllCCSReads_sorted.bam --output "${f%.*}Pilon.fa --outdir . --changes --fix all --threads 40 --chunksize 100000";done >Pilon.sh


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


```
