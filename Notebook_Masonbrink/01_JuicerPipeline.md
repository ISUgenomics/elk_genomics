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
ln -s ../../03_repeatmodeler/north_american_elk_15Jun2018_oY8t2.fasta.masked
module load bwa;bwa index ../north_american_elk_15Jun2018_oY8t2.fasta.masked
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
```
### modify and split job2.sh
```
sed -i 's/96:00:00/48:00:00/g' job2.sh
#had time to wait for the masked genome, so just ran as is without split
sbatch job2.sh

```

### job5.sh needs few mods, time, and path to countligations.sh
```
sed -i 's/96:00:00/48:00:00/g' job5.sh

#count countligations.sh was in path, removing faulty path
#/gpfs0/juicer//scripts/countligations.sh
countligations.sh
```

### job6.sh  and job7.sh
```
#this is just aligning the reads
for f in *fastq ;do echo "module load bwa; module load juicer; bwa mem -t 40 ../references/north_american_elk_15Jun2018_oY8t2.fasta.masked "$f" > "$f".sam";done >align.sh
```

### job8.sh
```
#sort the sam files
for f in *sam; do echo " sort --parallel=40 -S 110G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -k1,1 "$f" > "${f%.*}"_sort.sam";done  >sortsam.sh

#add read ends to the sam file read names
for f in splits/*_R1*_sort.sam ;do echo "awk 'NF >= 11{\$1 = \$1\"/1\";print}' "$f" > "${f%.*}"1.sam";done >addReadEnds.sh
for f in splits/*_R2*_sort.sam ;do echo "awk 'NF >= 11{\$1 = \$1\"/2\";print}' "$f" > "${f%.*}"1.sam";done >>addReadEnds.sh

#merge the sam files
 paste <(ls -1 *R1*.fastq_sort1.sam ) <(ls -1 *R2*.fastq_sort1.sam ) |while read a b; do echo "sort --parallel=40 -S 110G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -k1,1 -m " $a" "$b" > "${a%.*}merge.sam; done > merge.sh


 #Create the empty files and find the chimeric read alignments
for f in *sort1merge.sam;do echo "touch "${f%.*}".fastq_abnorm.sam"; done >createSamSubsets.sh
for f in *sort1merge.sam;do echo "touch "${f%.*}".fastq_unmapped.sam"; done >>createSamSubsets.sh
for f in *sort1merge.sam;do echo "touch "${f%.*}".fastq_norm.txt"; done >>createSamSubsets.sh
sh createSamSubsets.sh

#call chimeric reads in teh dataset
paste <(ls -1 *norm.txt ) <(ls -1 *abnorm.sam ) <(ls -1 *unmapped.sam ) <(ls -1 *1merge.sam)|while read a b c d; do echo "awk -v \"fname1\"="$a" -v \"fname2\"="$b" -v \"fname3\"="$c"  -f /software/7/apps/juicer/1.6.2/scripts/chimeric_blacklist.awk "$d;done >chimericReadAlignments
sh chimericReadAlignments

# something with the restriction fragment file, runs super fast
for f in *norm.txt; do echo "/software/7/apps/juicer/1.6.2/scripts/fragment.pl "$f" "${f%.*}".frag.txt ../restriction_sites/north_american_elk_15Jun2018_oY8t2.fasta_DpnII.txt" ;done |sed 's/_norm//2' >frags.sh

# sort by chromosome, fragment, strand, and position
for f in *frag.txt; do echo "sort -S 110G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n "$f" > "${f%.*}".fastq.sort.txt" ;done >sortByChromFragStrandPos.sh

#merge the sort.txt files into a single file
echo "sort --parallel=40 -S110G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/splits/*.sort.txt > /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/aligned/merged_sort.txt" >mergeAllSort.sh

awk -v queue=long -v groupname=a1560892488 -v debugdir=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/debug -v dir=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/aligned -v topDir=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk -v  juicedir=/software/7/apps/juicer/1.6.2/scripts/ -v site=DpnII -v genomeID=Elk -v genomePath=chrom.sizes -v user=rick.masonbrink -v guardjid="a1560892488_dedup_guard" -f /software/7/apps/juicer/1.6.2/scripts/split_rmdups.awk /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/aligned/merged_sort.txt
```

### 3ddna assembly and installation
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA
module load miniconda
conda create 3d-dna
conda create -n 3d-dna
source activate 3d-dna
git clone https://github.com/lastz/lastz.git
cd lastz/
# changed the install directory to one up from current
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA
cd src/
make
make install
make test

pip install numpy
pip install numpy --user
pip install scipy --user
pip install matplotlib --user

#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA
git clone https://github.com/theaidenlab/3d-dna.git
PATH=$PATH:/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA/lastz/
PATH=$PATH:/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA/3d-dna

module load miniconda
source activate 3d-dna
module load java
module load parallel
cd /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA/3d-dna
bash run-asm-pipeline.sh /home/rick.masonbrink/elk_bison_genomics/Masonbrink/north_american_elk_15Jun2018_oY8t2.fasta /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/aligned/merged_nodups.txt


#java error, had a couple reads that had coordinates off of the end of the chromosome.  removed.
#java error says this scaffold is the culprit. moving to new file
grep "ScoY8t2_23286\;HRSCAF\=23550" ../aligned/merged_nodups.txt >scaffScoY8t2_23286

#aiden form says reads can be mapping off of the scaffolds, these need removed.  some random error with bwa?
#are there any?
grep "ScoY8t2_23286\;HRSCAF\=23550" ../aligned/merged_nodups.txt >scaffScoY8t2_23286
awk '{if($2=="ScoY8t2_23286;HRSCAF=23550" && $3>115005250){next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7>115005250){next } else {print $0}}' scaffScoY8t2_23286 >fixedscafScoY8t2_23286

#fixed count
(bioawk) [rick.masonbrink@sn-cn-20-3 01_3DNA]$ wc -l fixedscafScoY8t2_23286
6700402 fixedscafScoY8t2_23286
#old count
(bioawk) [rick.masonbrink@sn-cn-20-3 01_3DNA]$ wc -l scaffScoY8t2_23286
6700404 scaffScoY8t2_23286

#apply this to merged_nodups.txt
awk '{if($2=="ScoY8t2_23286;HRSCAF=23550" && $3>115005250){next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7>115005250){next } else {print $0}}' merged_nodups.txt >fixedmerged_nodups.txt


#rerun 3d-dna

module load miniconda
source activate 3d-dna
module load java
module load parallel
cd /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA/3d-dna
bash run-asm-pipeline.sh /home/rick.masonbrink/elk_bison_genomics/Masonbrink/north_american_elk_15Jun2018_oY8t2.fasta /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/aligned/fixedmerged_nodups.txt

```
