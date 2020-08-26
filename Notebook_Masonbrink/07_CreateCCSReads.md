# Create ccs reads for polishing

Pilon has a hard time polishing the genome with subreads being so error prone.  Need to create CCS reads from the subreads.
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/10_CCS_polish

#softlink bax.h5 files
for f in ../../pacbio/*/*/Analysis_Results/*bax.h5; do ln -s $f;done
for f in ../../pacbio/*/Analysis_Results/*bax.h5; do ln -s $f;done

#convert old pacbio format to new bam format
ls -1 *bax.h5 |sed 's/\./\t/1' |awk '{print $1}' |sort|uniq|while read line;do echo "module load miniconda; source activate my_root; bax2bam "$line".1.bax.h5 "$line".2.bax.h5 "$line".3.bax.h5";done >bax2bam.sh


# make ccs
for f in *subreads.bam; do echo "sh PolishSubreads.sh "$f ; done >MakeCcs.sh

###############################################################################
#!/bin/bash
sr=$1
ccs=$(echo $1 |sed 's/.subreads/.ccs/g')
module load miniconda
source activate my_root
module load bamtools
ccs --polish --numThreads 40 --logLevel "FATAL" --maxLength=100000 --reportFile ${ccs%.*}_report.txt --minPasses 1 $sr $ccs
bamtools convert -format fasta -in $ccs > ${ccs%.*}.fa
###############################################################################
```
###  CCS Stats
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/10_CCS_polish
cat *ccs.fa |awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' |cut -f 2 |~/common_scripts/summary.sh
Total:  1,819,991,677
Count:  149,340
Mean:   12,186
Median: 13,197
Min:    117
Max:    48,127

```
### Subread stats
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/10_CCS_polish
for f in *subreads.bam; do echo "ml samtools; samtools view "$f" |awk '{print length(\$10)}' >"${f%.*}".lengths";done  >samtoolsLengths.sh

cat *lengths |~/common_scripts/summary.sh
Total:  52,353,735,044
Count:  5,702,794
Mean:   9,180
Median: 8,093
Min:    1

```
