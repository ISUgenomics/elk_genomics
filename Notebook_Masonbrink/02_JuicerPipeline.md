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

samtools flagstat -@ 40 lib_003.sorted.md.bam
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

# something with the restriction fragment file
for f in *norm.txt; do echo "/software/7/apps/juicer/1.6.2/scripts/fragment.pl "$f" "${f%.*}".frag.txt ../restriction_sites/north_american_elk_15Jun2018_oY8t2.fasta_DpnII.txt" ;done |sed 's/_norm//2' >frags.sh

###rerunning the above step
for f in *norm.txt; do echo "/software/7/apps/juicer/1.6.2/scripts/fragment.pl "$f" "${f%.*}".frag.txt ../restriction_sites/north_american_elk_15Jun2018_oY8t2.fasta.masked_DpnII.txt" ;done |sed 's/_norm//2' >frags.sh

# sort by chromosome, fragment, strand, and position
for f in *frag.txt; do echo "sort -S 110G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n "$f" > "${f%.*}".fastq.sort.txt" ;done >sortByChromFragStrandPos.sh
 ###

#merge the sort.txt files into a single file
echo "sort --parallel=40 -S110G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/splits/*.sort.txt > /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/aligned/merged_sort.txt" >mergeAllSort.sh

awk -v queue=long -v groupname=a1560892488 -v debugdir=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/debug -v dir=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/aligned -v topDir=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk -v  juicedir=/software/7/apps/juicer/1.6.2 -v site=DpnII -v genomeID=Elk -v genomePath=chrom.sizes -v user=rick.masonbrink -v guardjid="a1560892488_dedup_guard" -f split_rmdups.awk /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/aligned/merged_sort.txt
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


#no go.
Try again
awk '{if($2=="ScoY8t2_23286;HRSCAF=23550") {print $3} else if($6=="ScoY8t2_23286;HRSCAF=23550"){print $7 } else {next}}' merged_nodups.txt |sort |uniq -c |sort -k2,2nr |tail
      8 1362
      1 1361
      2 1360
     39 1359
      1 1356
      1 1355
      1 1354
      1 1353
      3 1352
      1 1351
[rick.masonbrink@sn-cn-22-0 aligned]$ awk '{if($2=="ScoY8t2_23286;HRSCAF=23550") {print $3} else if($6=="ScoY8t2_23286;HRSCAF=23550"){print $7 } else {next}}' merged_nodups.txt |sort |uniq -c |sort -k2,2nr |head
      1 115005254
      1 115005252
      2 115005250
      2 115005238
      1 115005236
      1 115005235
      1 115005221
      1 115005216
      1 115005211
      1 115005209
 awk '{if($2=="ScoY8t2_23286;HRSCAF=23550") {print $3} else if($6=="ScoY8t2_23286;HRSCAF=23550"){print $7 } else {next}}' merged_nodups.txt |sort |uniq -c |sort -k1,1nr |head
  37378 9290251
  34265 41362900
  29150 19388140
  25922 9288883
  21454 5228452
  19680 19388117
  19124 19388038
  15615 9290155
  15155 9289093
  11243 26246534


#stats from another large scaffold
awk '{if($2=="ScoY8t2_329;HRSCAF=386") {print $3} else if($6=="ScoY8t2_329;HRSCAF=386"){print $7 }
 else {next}}' merged_nodups.txt |sort |uniq -c |sort -k1,1nr |head
   2832 65883156
   1112 14005317
   1077 8335593
    970 39151896
    957 65883139
    796 8335576
    329 47970340
    317 47970348
    288 47970357
    281 6407397


It is likely that ScoY8t2_23286;HRSCAF=23550 has too much depth.  I will remove some of the entries, only those found in column 2 and 3.  Leaving col 6 and 7 alone.  
awk '{if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==9290251) {next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==41362900){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==19388140){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==9288883){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==5228452){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==19388117){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==19388038){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==9290155){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==9289093){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==26246534){next} else {print}}' merged_nodups.txt >fixedmerged_nodups.txt
#this did not remove enough reads, adding col 6 and 7 to the removal
awk '{if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==9290251) {next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==41362900){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==19388140){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==9288883){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==5228452){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==19388117){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==19388038){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==9290155){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==9289093){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $3==26246534){next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7==9290251) {next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7==41362900){next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7==19388140){next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7==9288883){next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7==5228452){next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7==19388117){next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7==19388038){next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7==9290155){next} else if($6=="ScoY8t2_23286;HRSCAF=23550" && $7==9289093){next} else if($2=="ScoY8t2_23286;HRSCAF=23550" && $7==26246534){next} else {print}}' merged_nodups.txt >fixedmerged_nodups.txt

#this was not complete also.  I checked to see if the merged_nodups.txt file had lines with 16 fields.  One line did not have 16 fields
awk 'NR!=72954653' fixedmerged_nodups.txt >fixed2merged_nodups.txt

#this fixed the issue, but I still have issues with the presentation of all intrachromosomal contacts in JBAT, likely due to a mismatched restriction fragment file.  7/12/19

#this issue is because the mapping information does not have a mapping site for restriction fragments in some chromosomes.  Just changing all values to dummy values, since I have only interest in assembly and 3d-dna is not using the fragment map
awk '{$8="1";$4="0";print $0}' fixed2merged_nodups.txt >fixed3merged_nodups.txt
```

### Run 3d-dna
```
module load miniconda
source activate 3d-dna
module load java
module load parallel
cd /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA/3d-dna
bash run-asm-pipeline.sh -i 100 /home/rick.masonbrink/elk_bison_genomics/Masonbrink/north_american_elk_15Jun2018_oY8t2.fasta /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/aligned/fixed3merged_nodups.txt
```


### Fix Scaffolds in Juicebox

```
I essentially tried to place all pseudomolecules in the same space by placing all unincorporated contigs at the end of each chromosome.  When all contigs greater than 10kb were placed, I edited each chromosome manually.  I followed the rule of only incorporating scaffolds to their site of highest signal, when signal was even between something repetitive and something non repetitive, I placed the contig in a the nonrepetitive region.   I kept repeats together when possible and only edited obvious misjoins, ~400 scaffolds.

#exported assembly
Final.assembly



```

### Seal and finalize assembly
```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA/3d-dna/01_PostJBAssembly/FinalAssembly/3d-dna

###Mostly Ns in this one
###ln -s ../../3d-dna/north_american_elk_15Jun2018_oY8t2.FINAL.FINAL.fasta
ln -s ../../3d-dna/north_american_elk_15Jun2018_oY8t2.FINAL.mnd.txt
ln -s ../../../north_american_elk_15Jun2018_oY8t2.FINAL.fasta
samtools faidx north_american_elk_15Jun2018_oY8t2.FINAL.fasta


module load miniconda
conda init bash
source activate 3d-dna
module load java
module load parallel
cd /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA/3d-dna/01_PostJBAssembly/FinalAssembly/3d-dna
sh run-asm-pipeline-post-review.sh -r Final.assembly -g 100 -s seal -i 100 --build-gapped-map north_american_elk_15Jun2018_oY8t2.FINAL.fasta north_american_elk_15Jun2018_oY8t2.FINAL.mnd.txt


#Fasta did not generate on the initial attempt, while the hic and assembly file were correct.
```

### Figured out that I merged the Y and the X chromosomes, so I had to work it out in juicebox.
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA/3d-dna/01_PostJBAssembly/FinalAssembly/01_XYSplit3d-dna

Copied in directory from ../3d-dna so I would not have to recreate the fasta database thingy.



#run 3d-dna with seal.  #turns out that this was a poor choice that led to more splits. #also does not generate the fasta, but the assembly and hic files were of appropriate size.

#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --partition=long
#SBATCH -t 168:00:00
#SBATCH -J finalizeAssembly1_0
#SBATCH -o finalizeAssembly1_0.o%j
#SBATCH -e finalizeAssembly1_0.e%j
#SBATCH --mail-user=remkv6@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module load miniconda
conda init bash
source activate 3d-dna
module load java
module load parallel
cd /home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk/01_3DNA/3d-dna/01_PostJBAssembly/FinalAssembly/01_XYSplit3d-dna
sh run-asm-pipeline-post-review.sh -r RemoveY13.assembly -g 100 -s seal -i 100  FirstScaffoldsStep0.FINAL.fasta north_american_elk_15Jun2018_oY8t2.mnd.txt
scontrol show job $SLURM_JOB_ID


#Perhaps run it like this #Still did not get a proper fasta
sh run-asm-pipeline-post-review.sh -r RemoveY13.assembly -g 100 -s finalize -i 100  FirstScaffoldsStep0.final.fasta north_american_elk_15Jun2018_oY8t2.mnd.txt


#Trying this
sh run-asm-pipeline-post-review.sh -r RemoveY13.assembly -g 100 -s finalize -i 100  ../3d-dna/FirstScaffoldsStep0.FINAL.fasta north_american_elk_15Jun2018_oY
#no effect.  

#had to generate the fasta manually

bash finalize/construct-fasta-from-asm.sh FirstScaffoldsStep0.final.final.cprops FirstScaffoldsStep0.final.final.asm ../3d-dna/FirstScaffoldsStep0.FINAL.fasta >test &

#this left some gaps in the sequence, so I had to make it one line fasta tab, and then correct it as below.
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < test >test2
tr "\t" "\n" < test2 | fold -w 60 >test3
mv test3 FirstScaffoldsStep0.final.FINAL.fasta


#is it the correct size?
bioawk -c fastx '{print $name,length($seq)}' FirstScaffoldsStep0.final.FINAL.fasta |awk '{print $2}' |summary.sh
Total:  2,559,916,354
Count:  22,564
Mean:   113,451
Median: 387
Min:    0
Max:    145,554,050


bioawk -c fastx '{print $name,length($seq)}' FirstScaffoldsStep0.FINAL.fasta |awk '{print $2}' |summary.sh
Total:  2,560,751,668
Count:  22,577
Mean:   113,423
Median: 387
Min:    119
Max:    153,425,878

#This is a negligable difference that will be fixed by pilon
```
