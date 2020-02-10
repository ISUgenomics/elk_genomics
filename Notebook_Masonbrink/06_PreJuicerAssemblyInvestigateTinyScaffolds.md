# Juicebox doesnt know where to put tiny scaffolds, are there duplicate seqs in the genome already?
```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/05_InvestigateTinyScaffolds

#separate out the little scaffs
 bioawk -c fastx '{print $name,length($seq)}' Modnorth_american_elk_15Jun2018_oY8t2.fasta |awk '$2<1000 {print $1}' |while read line; do samtools faidx Modnorth_american_elk_15Jun2018_oY8t2.fasta $line; done>LittleScaffs.fasta
grep -c ">" LittleScaffs.fasta
19190


#separate out the big scaffs
 bioawk -c fastx '{print $name,length($seq)}' Modnorth_american_elk_15Jun2018_oY8t2.fasta |awk '$2>1000 {print $1}' |while read line; do samtools faidx Modnorth_american_elk_15Jun2018_oY8t2.fasta $line; done>BigScaffs.fasta
 grep -c ">" BigScaffs.fasta
4111


#make the blast database
module load blast+
makeblastdb -in BigScaffs.fasta -dbtype nucl
blastn -db BigScaffs.fasta -query LittleScaffs.fasta -num_threads 40 -outfmt 6 -out Small2Big.blastout

# get the best hit for each.
sort -u -k1,1nr Small2Big.blastout >TopHit.blastout


## get scaffold length associated with blast hits
awk '{print $1}' TopHit.blastout |while read line; do samtools faidx Modnorth_american_elk_15Jun2018_oY8t2.fasta $line; done |bioawk -c fastx '{print $name,length($seq)}' |paste - TopHit.blastout >ScaffLen2TopHit.blastout &


## How many are 80% contained and higher than 90% identical
less ScaffLen2TopHit.blastout |awk '$6>(.8*$2) && $5>.9' |wc
  18795  263130 2323854

#How many are 90% contained and and 90% identical
less ScaffLen2TopHit.blastout |awk '$6>(.9*$2) && $5>.9' |wc
  18579  260106 2296849
```
We are likely to remove 18,579 scaffolds from the current assembly of 23302 scaffolds, leaving 4,723 scaffolds.
