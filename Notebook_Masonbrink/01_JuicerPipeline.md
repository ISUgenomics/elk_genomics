# Need to improve the elk assembly using HiC reads and the juicer pipeline.

# What was the previous mapping percentage for the hic reads
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/01_ExtractFastq

samtools flagstat -@ 40 lib_002.sorted.md.bam
333824140 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
135601503 + 0 duplicates
289958389 + 0 mapped (86.86% : N/A)
333824140 + 0 paired in sequencing
166912070 + 0 read1
166912070 + 0 read2
264469578 + 0 properly paired (79.22% : N/A)
264469578 + 0 with itself and mate mapped
25488811 + 0 singletons (7.64% : N/A)
125097186 + 0 with mate mapped to a different chr
75127210 + 0 with mate mapped to a different chr (mapQ>=5)
[rick.masonbrink@sn-cn-15-2 01_ExtractFastq]$ samtools flagstat -@ 40 lib_001.sorted.md.bam
302679286 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
121184776 + 0 duplicates
251726759 + 0 mapped (83.17% : N/A)
302679286 + 0 paired in sequencing
151339643 + 0 read1
151339643 + 0 read2
215556276 + 0 properly paired (71.22% : N/A)
215556276 + 0 with itself and mate mapped
36170483 + 0 singletons (11.95% : N/A)
131992084 + 0 with mate mapped to a different chr
87874904 + 0 with mate mapped to a different chr (mapQ>=5)
[rick.masonbrink@sn-cn-15-2 01_ExtractFastq]$ samtools flagstat -@ 40 lib_003.sorted.md.bam
367950046 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
151707012 + 0 duplicates
321002811 + 0 mapped (87.24% : N/A)
367950046 + 0 paired in sequencing
183975023 + 0 read1
183975023 + 0 read2
292808190 + 0 properly paired (79.58% : N/A)
292808190 + 0 with itself and mate mapped
28194621 + 0 singletons (7.66% : N/A)
138772472 + 0 with mate mapped to a different chr
83321293 + 0 with mate mapped to a different chr (mapQ>=5)
[rick.masonbrink@sn-cn-15-2 01_ExtractFastq]$ ^C
[rick.masonbrink@sn-cn-15-2 01_ExtractFastq]$


```

### Extract fastq
```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/01_ExtractFastq

for f in *bam ; do echo "module load picard/64/2.9.2; java -jar /software/7/apps/picard/64/2.9.2/picard.jar SamToFastq I="$f" FASTQ="${f%.*}"R1.fq F2="${f%.*}"R2.fq FU="${f%.*}"_unpaired.fq" ;done >picard.sh
#######################################################################################
 module load picard/64/2.9.2; java -jar /software/7/apps/picard/64/2.9.2/picard.jar SamToFastq I=lib_001.sorted.md.bam FASTQ=lib_001.sorted.mdR1.fq F2=lib_001.sorted.mdR2.fq FU=lib_001.sorted.md_unpaired.fq
module load picard/64/2.9.2; java -jar /software/7/apps/picard/64/2.9.2/picard.jar SamToFastq I=lib_002.sorted.md.bam FASTQ=lib_002.sorted.mdR1.fq F2=lib_002.sorted.mdR2.fq FU=lib_002.sorted.md_unpaired.fq
module load picard/64/2.9.2; java -jar /software/7/apps/picard/64/2.9.2/picard.jar SamToFastq I=lib_003.sorted.md.bam FASTQ=lib_003.sorted.mdR1.fq F2=lib_003.sorted.mdR2.fq FU=lib_003.sorted.md_unpaired.fq
#######################################################################################

```


###  Repeatmask the genome
```

#note that this needed at least 96 hours to complete, Better to split into two nodes for the next round with this large elk genome
sh ~/common_scripts/runRepeatModeler.sh north_american_elk_15Jun2018_oY8t2.fasta


##################################################################################
#!/bin/bash


if [ $# -lt 1 ] ; then
        echo "usage: runRepeatModeler <genome.fasta>"
        echo ""
        echo "To build custom repeat library and mask the repeats of the genome"
        echo ""
exit 0
fi
GENOME="$1"

module purge
module load parallel
module load repeatmasker/4.0.7
module load repeatmodeler/1.0.8
module load perl/5.24.1
DATABASE="$(basename ${GENOME%.*}).DB"
BuildDatabase -name ${DATABASE} -engine ncbi ${GENOME}
RepeatModeler -database ${DATABASE}  -engine ncbi -pa 40
ln -s $(find $(pwd) -name "consensi.fa.classified")
RepeatMasker -pa 40  -gff -lib consensi.fa.classified ${GENOME}
##################################################################################

```

### Set up Juicer
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk

#chromosome sizes file
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' ../north_american_elk_15Jun2018_oY8t2.fasta |tr "\n" "\t" |sed 's/>/\n/g' >chrom.sizes

mkdir HIC_tmp
mkdir fastq

mkdir references
cd references
module load bwa;bwa index ../north_american_elk_15Jun2018_oY8t2.fasta
cd ..

mkdir restriction
cd restriction
vi generate_site_positions.py
module load juicer
python generate_site_positions.py DpnII north_american_elk_15Jun2018_oY8t2.fasta ../north_american_elk_15Jun2018_oY8t2.fasta
cd ..

```
