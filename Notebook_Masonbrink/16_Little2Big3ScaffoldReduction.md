# Finished polishing the genome, so lets check to see if any of the smaller scaffolds will now align to the larger ones

```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/26_RefineLittle2Big3
ln -s ../14_PilonPairedEnd/01_SplitForPilon/02_MaskPart16/FinalGenomePilon.fa

module load samtools

samtools faidx FinalGenomePilon.fa

bioawk -c fastx '{print $name,length($seq)}' FinalGenomePilon.fa |sort -k2,2nr |awk 'NR>35 {print $1}' |while read line; do samtools faidx FinalGenomePilon.fa $line; done >LittleScaffs.fasta

#separate big scaffs from little
bioawk -c fastx '{print $name,length($seq)}' FinalGenomePilon.fa |sort -k2,2nr |awk 'NR<36 {print $1}' |while read line; do samtools faidx FinalGenomePilon.fa $line; done >BigScaffs.fasta


#blast little to big
module load blast+
makeblastdb -in BigScaffs.fasta -dbtype nucl
echo "module load blast+;blastn -db BigScaffs.fasta -query LittleScaffs.fasta -outfmt 6 -num_threads 40 -dust no -max_target_seqs 5 -out Small2Big.blastout" >blast.sh

#get little scaff lengths
bioawk -c fastx '{print $name,length($seq)}' LittleScaffs.fasta >LittleScaffLengths.txt

#get the covered regions of the short scaffolds that are in the large scaffolds
cat *blastout|awk '$12>300'|awk '{print $1,$7,$8}' |awk '{if($2>$3){print $1,$3,$2}else {print $1,$2,$3}}' |sort -k1,1V -k2,2n|tr " " "\t" |bedtools merge -d 1000 >LittleScaffoldRepeatBlastRepMasker.bed

bioawk -c fastx '{print $name,length($seq)}' LittleScaffs.fasta >LittleScaffLengths.txt
less LittleScaffoldRepeatBlastRepMasker.bed |awk '{print $1}' |while read line; do grep -w $line LittleScaffLengths.txt;done |paste - LittleScaffoldRepeatBlastRepMasker.bed >ScaffLengthLittleScaffoldRepeatBlastRepMasker.bedlike

awk -v lengths=0 '{if(name==$1) {count=($5-$4)+count;lengths=$2;name=$1} else if($1!=name) {print name,lengths,count;name=$1;lengths=$2;count=0;count=count+($5-$4)} }' ScaffLengthLittleScaffoldRepeatBlastRepMasker.bedlike >List2refine.bedlike


#How many have at least 90% of their length contained in a pseudomolecule
less List2refine.bedlike  |awk '$3<($2*.9)' |wc
    151     451    5086

```
