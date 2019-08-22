#  How many of those little scaffolds remaining are actually uncollapsed haplotigs and are nearly completely contained in the larger scaffolds anyway?

### Set up fasta subsets
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/08_InvestigateTinyScaffolds4Juiced

ln -s ../04_JuicerElk/01_3DNA/3d-dna/01_PostJBAssembly/FinalAssembly/3d-dna/FirstScaffoldsStep0.FINAL.fasta

#Top 34 scaffolds vs all others

bioawk -c fastx '{print $name,length($seq)}' FirstScaffoldsStep0.FINAL.fasta |sort -k2,2n |awk 'nr <35' |while read line;do samtools faidx FirstScaffoldsStep0.FINAL.fasta ;done >BigScaffs.fasta
bioawk -c fastx '{print $name,length($seq)}' FirstScaffoldsStep0.FINAL.fasta |sort -k2,2n |awk 'nr >35' |while read line;do samtools faidx FirstScaffoldsStep0.FINAL.fasta;done >littleScaffs.fasta
```
### set up blast
```
module load blast+
makeblastdb -in BigScaffs.fast -dbtype nucl

blastn -db BigScaffs.fasta -query littleScaffs.fasta -num_threads 40 -culling_limit 5  -outfmt 6 -out little2big.blastout
```

### Filter and analyze
```
How many little contigs are there?
grep -c ">" littleScaffs.fasta
22543


sort -k1,1 -u little2big.blastout >TopHitUniq.blastout
bioawk -c fastx '{print $name,length($seq)}' FirstScaffoldsStep0.FINAL.fasta >Scaffold.lengths

awk '{print $1}' TopHitUniq.blastout |while read line; do grep -w $line Scaffold.lengths; done |paste - TopHitUniq.blastout >LengthNBlast.out &
# A single blast hit from query scaffold to larger scaffolds that contained more than 90% of the scaffold length into a large scaffold.

90%
less LengthNBlast.out |awk '($6/$2)>.9' |wc
 16036  224504 1723107

 80%?
 less LengthNBlast.out |awk '($6/$2)>.8' |wc
   20446  286244 2209194

How many did not have a perfect blast hit that enclosed the entire scaffold within a pseudomolecule?
less LengthNBlast.out |awk '($6/$2)>.8 {print $1}' |cat - <(grep ">" littleScaffs.fasta|sed 's/>//g' ) |sort|uniq -c |awk '$1==1'|wc
   1937    3874   50529
less LengthNBlast.out |awk '($6/$2)>.8 {print $1}' |cat - <(grep ">" littleScaffs.fasta|sed 's/>//g' ) |sort|uniq -c |awk '$1==1 {print $2}'>littleScaffsWOPerfectHit.list

#filter 90% identity
grep -f littleScaffsWOPerfectHit.list little2big.blastout|awk '$3>90' >littleScaffsWOPerfectHitFiltered.blastout

#create bed file and merge coordinates from hits
less littleScaffsWOPerfectHitFiltered.blastout |awk '{if($7>$8){print $1,$8,$7}else {print $1,$7,$8}}' |sort -k1,1V -k2,2n |tr " " "\t" |bedtools merge -d 300 >littleScaffsWOPerfectHitFilteredHitsMerged.bed

less littleScaffsWOPerfectHitFilteredHitsMergedLengths.tab|awk '$2!=$5' |awk '{print $1,$2,$5,$4,$5-$4}' |awk '$5>200' >List2Refine.bedlike

#how many left?
wc List2Refine.bedlike
 1909  9545 69077 List2Refine.bedlike


Remove those that do not have a single read of pacbio coverage. how many are left?
less  List2Refine.bedlike |awk '$5<(.9*$2) {print $1}' |sort |uniq|grep -w -f - ../07_blobtools2juiced/blobplot_out.FirstScaffoldsStep0.FINAL.AllFastq_sorted.bam.cov |awk '$2!=0' |wc
    473    1419   13075

```

###  manual checks to see if these scaffolds can fit anywhere with a consensus
```
#HiC_scaffold_14 is split among many spots and is ~3.6kb --high coverage, but obviously a repeat
#HiC_scaffold_19 is ~10kb and split among many spots.  --high coverage, obviously a repeat
#HiC_scaffold_21 is ~3-4kb and split among many spots.  --high coverage, obviously a repeat


```
