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
 module load picard/64/2.9.2; java -Xmx200G -Xms50G -jar /software/7/apps/picard/64/2.9.2/picard.jar SamToFastq I=lib_001.sorted.md.bam FASTQ=lib_001.sorted.mdR1.fq F2=lib_001.sorted.mdR2.fq FU=lib_001.sorted.md_unpaired.fq
module load picard/64/2.9.2; java -Xmx200G -Xms50G -jar /software/7/apps/picard/64/2.9.2/picard.jar SamToFastq I=lib_002.sorted.md.bam FASTQ=lib_002.sorted.mdR1.fq F2=lib_002.sorted.mdR2.fq FU=lib_002.sorted.md_unpaired.fq
module load picard/64/2.9.2; java -Xmx200G -Xms50G -jar /software/7/apps/picard/64/2.9.2/picard.jar SamToFastq I=lib_003.sorted.md.bam FASTQ=lib_003.sorted.mdR1.fq F2=lib_003.sorted.mdR2.fq FU=lib_003.sorted.md_unpaired.fq
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

#need this guy for sorting
mkdir HIC_tmp

#softlink fastq
mkdir fastq
cd fastq
for f in ../../01_ExtractFastq/*mdR*fq; do ln -s $f;done


#softlink and index references
mkdir references
cd references
ln -s ../north_american_elk_15Jun2018_oY8t2.fasta
module load bwa;bwa index ../north_american_elk_15Jun2018_oY8t2.fasta
cd ..


#get your proper restriction enzyme for a fragment delimited map
mkdir restriction
cd restriction
vi generate_site_positions.py
module load juicer
python generate_site_positions.py DpnII north_american_elk_15Jun2018_oY8t2.fasta ../north_american_elk_15Jun2018_oY8t2.fasta
cd ..

#get modified juicer script
cp ../02_TestJuicer/new_Juicer.sh .

#softlink and rename according to "*_R*", the fastq files to the appropriate format for juicer
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/fastq
for f in ../../01_ExtractFastq/*fq; do echo "ln -s "$f" "${f%.*}".fastq" ;done |sed 's|../../01_ExtractFastq/||2' |sed 's/.mdR/_R/2' |less
###############################################################################
ln -s ../../01_ExtractFastq/lib_001.sorted.mdR1.fq lib_001.sorted_R1.fastq
ln -s ../../01_ExtractFastq/lib_001.sorted.mdR2.fq lib_001.sorted_R2.fastq
ln -s ../../01_ExtractFastq/lib_002.sorted.mdR1.fq lib_002.sorted_R1.fastq
ln -s ../../01_ExtractFastq/lib_002.sorted.mdR2.fq lib_002.sorted_R2.fastq
ln -s ../../01_ExtractFastq/lib_003.sorted.mdR1.fq lib_003.sorted_R1.fastq
ln -s ../../01_ExtractFastq/lib_003.sorted.mdR2.fq lib_003.sorted_R2.fastq
###############################################################################

new_Juicer.sh -y restriction_sites/MaskedMisAssFixed.Pilon.fasta_DpnII.txt -q short -Q 48:00:00 -z references/MaskedMisAssFixed.Pilon.fasta -p chrom.sizes -t 40

#job1.sh just prints the version numberss of the software -- skipped

#modify and split job2.sh
sed -i 's/96:00:00/48:00:00/g' job2.sh
#had time to wait for the masked genome, so just ran as is without split
sbatch job2.sh


#job5.sh needs few mods, time, and path to countligations.sh
sed -i 's/96:00:00/48:00:00/g' job5.sh

#count countligations.sh was in path, removing faulty path
#/gpfs0/juicer//scripts/countligations.sh
countligations.sh



```
