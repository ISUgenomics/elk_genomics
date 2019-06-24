### Testing juicer pipeline on ceres with small genome

## Use SCN as example 129mb genome size

Copied directory contents from condo do ceres, requires folders: fastq, references, restriction.  Also need files for fastq, reference +index, restriction enzyme file.

```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer
module load juicer

#generate the job scripts
new_Juicer.sh -y restriction_sites/MisAssFixed.Pilon.fasta_MboI.txt -q short -Q 48:00:00 -z references/MaskedMisAssFixed.Pilon.fasta -p chrom.sizes -t 40
```
### run the jobs generated
```
#job1.sh just has scripts to get version numbers of the software, skipping.

#job3.sh is just splitting of the fastq files
sed -i 's/96:00:00/48:00:00/g' job2.sh
split -l 11 job3.sh
for f in *sub; do sbatch $f;done
```

### job5.sh count ligations

```
sed -i 's/96:00:00/48:00:00/g' job5.sh
#count countligations.sh was in path, removing faulty path
#/gpfs0/juicer//scripts/countligations.sh
countligations.sh
```

### align reads to reference (job6 and job7)
```
job6.sh and job7.sh are aligning both ends of the read pair as single end.

#GRAB the respective line from each file that needs to be ran (alignment)

#job6.sh
###############################################################################
bwa mem -t 40 references/MaskedMisAssFixed.Pilon.fasta /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R1*.fastq > /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R1*.fastq.sam
###############################################################################

#job7.sh
###############################################################################
 bwa mem -t 40 references/MaskedMisAssFixed.Pilon.fasta /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R2*.fastq > /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R2*.fastq.sam
###############################################################################

#job6 and 7 are just separate the read ends, which is not necessary here.
#make sure you are in the splits directory or you will get the folder name attached in the submission scripts  
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits
for f in *fastq ;do echo "module load bwa; module load juicer; bwa mem -t 40 references/MaskedMisAssFixed.Pilon.fasta "$f" > "$f".sam";done >align.sh

```

### job 8.sh
```
There are many steps here.  
#sort bam files by read name
#########################################################################################
sort --parallel=40 -S 14G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -k1,1 /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R1*.fastq.sam > /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R1*.fastq_sort.sam
#########################################################################################

for f in splits/*sam; do echo " sort --parallel=40 -S 100G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -k1,1 "$f" > "${f%.*}"_sort.sam";done  |less

# remove sam header and add read end indicatory
################################################################################
awk 'NF >= 11{\$1 = \$1"/1";print}' /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R1*.fastq_sort.sam > /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R1*.fastq_sort1.sam
awk 'NF >= 11{\$1 = \$1"/2";print}' /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R2*.fastq_sort.sam > /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R2*.fastq_sort1.sam
################################################################################

for f in splits/*_R1*_sort.sam ;do echo "awk 'NF >= 11{\$1 = \$1"/1";print}' "$f" > "${f%.*}"1.sam";done >addReadEnds.sh
for f in splits/*_R2*_sort.sam ;do echo "awk 'NF >= 11{\$1 = \$1"/2";print}' "$f" > "${f%.*}"1.sam";done >>addReadEnds.sh


# merge bam files, juicer wants to merge these and overwrite the original sam.  I added a merged at the end to differentiate.
################################################################################
sort --parallel=8 -S 14G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -k1,1 -m /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R1*.fastq_sort1.sam /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/*_R2*.fastq_sort1.sam > /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq.sam
##################################################################################

paste <(ls -1 splits/*R1*fastq_sort1.sam ) <(ls -1 splits/*R2*fastq_sort1.sam ) |while read a b; do echo "sort --parallel=40 -S 110G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -k1,1 -m " $a" "$b" > "${a%.*}merge.sam; done > merge.sh


# Identify chimeric read alignments
# call chimeric_blacklist.awk to deal with chimeric reads; sorted file is sorted by read name at this point
###########################################################################################################################################################################################################################
touch /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq_abnorm.sam /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq_unmapped.sam
awk -v "fname1"=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq_norm.txt -v "fname2"=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq_abnorm.sam -v "fname3"=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq_unmapped.sam -f /gpfs0/juicer//scripts/chimeric_blacklist.awk /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq.sam
###########################################################################################################################################################################################################################

#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits

for f in *sort1merge.sam;do echo "touch "${f%.*}".fastq_abnorm.sam"; done >createSamSubsets.sh
for f in *sort1merge.sam;do echo "touch "${f%.*}".fastq_unmapped.sam"; done >>createSamSubsets.sh
for f in *sort1merge.sam;do echo "touch "${f%.*}".fastq_norm.txt"; done >>createSamSubsets.sh

#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits
paste <(ls -1 *norm.txt ) <(ls -1 *abnorm.sam ) <(ls -1 *unmapped.sam ) <(ls -1 *1merge.sam)|while read a b c d; do echo "awk -v \"fname1\"="$a" -v \"fname2\"="$b" -v \"fname3\"="$c"
 -f /software/7/apps/juicer/1.6.2/scripts/chimeric_blacklist.awk "$f;done >chimericReadAlignments


#do something with the fragment delimitation
##############################################################################################
/software/7/apps/juicer/1.6.2/scripts/fragment.pl /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq_norm.txt /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq.frag.txt restriction_sites/MisAssFixed.Pilon.fasta_MboI.txt

##############################################################################################

for f in *norm.txt; do echo "/software/7/apps/juicer/1.6.2/scripts/fragment.pl "$f" "${f%.*}".frag.txt restriction_sites/MisAssFixed.Pilon.fasta_DpnII.txt" ;done |sed 's/_norm//2' >frags.sh

# sort by chromosome, fragment, strand, and position
#####################################################################################################################################################################################################################
sort -S 2G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq.frag.txt > /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/splits/**.fastq.sort.txt
######################################################################################################################################################################################################################

for f in *frag.txt; do echo "sort -S 110G -T /home/rick.masonbrink/elk_bison_genomics/Masonbrink/02_TestJuicer/HIC_tmp -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n "$f" > "${f%.*}".fastq.sort.txt" ;done >sortByChromFragStrandPos.sh


```
