# Juicebox doesnt know where to put tiny scaffolds, are there duplicate seqs in the genome already?
```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/05_InvestigateTinyScaffolds

#separate out the little scaffs
 bioawk -c fastx '{print $name,length($seq)}' Modnorth_american_elk_15Jun2018_oY8t2.fasta |awk '$2<1000 {print $1}' |while read line; do samtools faidx Modnorth_american_elk_15Jun2018_oY8t2.fasta $line; done>LittleScaffs.fasta

#separate out the big scaffs
 bioawk -c fastx '{print $name,length($seq)}' Modnorth_american_elk_15Jun2018_oY8t2.fasta |awk '$2>1000 {print $1}' |while read line; do samtools faidx Modnorth_american_elk_15Jun2018_oY8t2.fasta $line; done>BigScaffs.fasta

#make the blast database
module load blast+
makeblastdb -in BigScaffs.fasta -dbtype nucl
blastn -db BigScaffs.fasta -query LittleScaffs.fasta -num_threads 40 -outfmt 6 -out Small2Big.blastout

```
