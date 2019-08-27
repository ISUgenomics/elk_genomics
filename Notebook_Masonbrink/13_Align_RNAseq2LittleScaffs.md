## Align RNA-seq to small scaffolds.

I will likely remove scaffolds that do not have a unique rnaseq reads mapping.

```
16_RNAseq/

hisat2-build FinalAssemblyFastaWithY.fasta FinalAssemblyFastaWithY


module load samtools
samtools merge -@ 40 AllReads_sorted.bam *sorted.bam &



cp ~/common_scripts/runFeatureCounts.sh


module load miniconda
source activate bioawk
bioawk -c fastx '{print $name,"rnaseq","gene","1",length($seq),".",".",".","ID="$name}' LittleScaffs.fasta |tr " " "\t" >genome.gff



paste <(ls -1 *R1*fastq) <(ls -1 *R2*fastq) |while read line; do echo "sh runFeatureCounts.sh "$line" /home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseqLittleScaffs.fasta.masked genome.gff" ;done >featureCounts.sh

 less AllReads_GeneReadCounts.txt |awk '$7!=0 {print $1}' |while read line; do grep -w $line ../13_Little2Big2/goThroughTheseByHand.bedlike;done >RNAseqMappingOnly.bedlike

```
