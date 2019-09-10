## Need to see if these scaffolds are just haplotypes and investigate for elimination


```
############################
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/13_Little2Big2

ln -s ../04_JuicerElk/01_3DNA/3d-dna/01_PostJBAssembly/FinalAssembly/3d-dna/FinalAssemblyFastaWithY.fasta

module load samtools
samtools faidx FinalAssemblyFastaWithY.fasta &



module load miniconda
source activate bioawk
#separate little scaffs from big
bioawk -c fastx '{print $name,length($seq)}' FinalAssemblyFastaWithY.fasta |sort -k2,2nr |awk 'NR>35 {print $1}' |while read line; do samtools faidx FinalAssemblyFastaWithY.fasta $line; done >LittleScaffs.fasta

#separate big scaffs from little
bioawk -c fastx '{print $name,length($seq)}' FinalAssemblyFastaWithY.fasta |sort -k2,2nr |awk 'NR<36 {print $1}' |while read line; do samtools faidx FinalAssemblyFastaWithY.fasta $line; done >BigScaffs.fasta


#blast little to big
module load blast+
makeblastdb -in BigScaffs.fasta -dbtype nucl
echo "module load blast+;blastn -db BigScaffs.fasta -query LittleScaffs.fasta -outfmt 6 -num_threads 40 -dust no -max_target_seqs 5 -out Small2Big.blastout" >blast.sh

#get little scaff lengths
bioawk -c fastx '{print $name,length($seq)}' LittleScaffs.fasta >LittleScaffLengths.txt

fasta-splitter --n-parts 7 LittleScaffs.fasta
#not sure if I used this one bc the blast file was too large
less Small2Big.blastout|awk '{print $1}'|while read line; do grep -w $line LittleScaffLengths.txt; done >LengthsOrderedByBlast
```

##RepeatMask the little scaffolds, mostly just to get labeled coordinates that would not show up in a blast with dust on (i.e. tandem repeats)
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/13_Little2Big2/01_RepeatMasker

ln -s ../LittleScaffs.fasta
ln -s ../../03_repeatmodeler/consensi.fa.classified


sh ./runRepeatModeler.sh  LittleScaffs.fasta

#runRepeatModeler.sh  -- modified to only mask
###############################################################################################################################################
#!/bin/bash
# runs repeat masking for the genome after constructing custom repeat library
# uses repeat modeler for building custom db and RepeatMasking for masking
# run it as:
# runRepeatModeler.sh Genome.fasta
# based on Rick's guide https://intranet.gif.biotech.iastate.edu/doku.php/people:remkv6:genome738polished_repeatmodeler_--de_novo_repeat_identification

if [ $# -lt 1 ] ; then
        echo "usage: runRepeatModeler <genome.fasta>"
        echo ""
        echo "To build custom repeat library and mask the repeats of the genome"
        echo ""
exit 0
fi


GENOME="$1"
#module use /shared/software/GIF/modules/
module purge
#module load parallel

module load repeatmasker/4.0.7
#module load repeatmodeler/1.0.8
#module load perl/5.24.1
DATABASE="$(basename ${GENOME%.*}).DB"
#BuildDatabase -name ${DATABASE} -engine ncbi ${GENOME}
#RepeatModeler -database ${DATABASE}  -engine ncbi -pa 40
#ln -s $(find $(pwd) -name "consensi.fa.classified")
RepeatMasker -pa 40  -gff -lib consensi.fa.classified ${GENOME}
##########################################################################################################################################################

==================================================
file name: LittleScaffs.fasta
sequences:         22522
total length:   35608595 bp  (35583022 bp excl N/X-runs)
GC level:         53.07 %
bases masked:   28508858 bp ( 80.06 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
SINEs:             1356       210530 bp    0.59 %
      ALUs            0            0 bp    0.00 %
      MIRs          432        53217 bp    0.15 %

LINEs:            13883      5784020 bp   16.24 %
      LINE1        7880      3636397 bp   10.21 %
      LINE2         255        36822 bp    0.10 %
      L3/CR1          0            0 bp    0.00 %

LTR elements:      2102       990609 bp    2.78 %
      ERVL          207        64002 bp    0.18 %
      ERVL-MaLRs    595       106737 bp    0.30 %
      ERV_classI   1177       703332 bp    1.98 %
      ERV_classII   123       116538 bp    0.33 %

DNA elements:       911       185676 bp    0.52 %
     hAT-Charlie    312        52276 bp    0.15 %
     TcMar-Tigger   127        34210 bp    0.10 %

Unclassified:       400        90323 bp    0.25 %

Total interspersed repeats:  7261158 bp   20.39 %


Small RNA:          940       158444 bp    0.44 %

Satellites:       18687     21103509 bp   59.27 %
Simple repeats:    2401       138771 bp    0.39 %
Low complexity:     318        15981 bp    0.04 %
==================================================

```




# Figure out which scaffolds deserve elimination.  

Scaffolds less than 1kb were eliminated, scaffolds that were at least 80% covered in the large scaffolds were removed, and scaffolds that did not map a unique rnaseq read were removed.
```
cat *blastout|awk '$12>300'|awk '{print $1,$7,$8}' |cat - <(awk '{print $1,$4,$5}' 01_RepeatMasker/LittleScaffs.fasta.out.gff) |awk '{if($2>$3){print $1,$3,$2}else {print $1,$2,$3}}' |sort -k1,1V -k2,2n|tr " " "\t" |bedtools merge -d 1000 >LittleScaffoldRepeatBlastRepMasker.bed


less LittleScaffoldRepeatBlastRepMasker.bed |awk '{print $1}' |while read line; do grep -w $line LittleScaffLengths.txt;done |paste - LittleScaffoldRepeatBlastRepMasker.bed >ScaffLengthLittleScaffoldRepeatBlastRepMasker.bedlike



#How many are left after removing those without rnaseq reads mapping
awk '{print $1}' goThroughTheseByHand.bedlike|sort|uniq|while read line; do grep -w $line ../16_RNAseq/AllReads_GeneReadCounts.txt;done |awk '$7>0' |wc
    901    6307   48165



#get the total length of blast+repeatmasker coordinates for each scaffold and eliminate scaffolds with 80% coverage.
awk -v lengths=0 '{if(name==$1) {count=($5-$4)+count;lengths=$2;name=$1} else if($1!=name) {print name,lengths,count;name=$1;lengths=$2;count=0;count=count+($5-$4)} }' ScaffLengthLittleScaffoldRepeatBlastRepMasker.bedlike >List2refine.bedlike

awk '$7>0 {print $1}' ../16_RNAseq/AllReads_GeneReadCounts.txt |while read line; do grep -w $line List2refine.bedlike;done |awk '($3)<(.9*$2)' |awk '$2>1000' >goThroughTheseByHand.bedlike

awk '{print $1}' goThroughTheseByHand.bedlike|sort|uniq|wc
    458     458    8152

#Manually deleted HiC_scaffold_58 -- the concatenated mitochondria

less goThroughTheseByHand.bedlike|awk '{print $1}' |while read line; do samtools faidx LittleScaffs.fasta $line; done >LittleScaffoldsRemaining.fasta

cat ../12_C.elaphusHippelaphus2/OurElkMitochondria.fasta LittleScaffoldsRemaining.fasta BigScaffs.fasta >FinalGenome.fa


---------------- Information for assembly '../FinalGenome.fa' ----------------


                                         Number of scaffolds        493
                                     Total size of scaffolds 2529467193
                                            Longest scaffold  146207705
                                           Shortest scaffold       1004
                                 Number of scaffolds > 1K nt        493 100.0%
                                Number of scaffolds > 10K nt        177  35.9%
                               Number of scaffolds > 100K nt         35   7.1%
                                 Number of scaffolds > 1M nt         35   7.1%
                                Number of scaffolds > 10M nt         34   6.9%
                                          Mean scaffold size    5130765
                                        Median scaffold size       6000
                                         N50 scaffold length   77665974
                                          L50 scaffold count         13
                                                 scaffold %A      29.18
                                                 scaffold %C      20.79
                                                 scaffold %G      20.79
                                                 scaffold %T      29.16
                                                 scaffold %N       0.09
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      99.9%
              Percentage of assembly in unscaffolded contigs       0.1%
                      Average number of contigs per scaffold       16.3
Average length of break (>25 Ns) between contigs in scaffold        298

                                           Number of contigs       8016
                              Number of contigs in scaffolds       7594
                          Number of contigs not in scaffolds        422
                                       Total size of contigs 2527225143
                                              Longest contig    4878448
                                             Shortest contig        180
                                   Number of contigs > 1K nt       8004  99.9%
                                  Number of contigs > 10K nt       6781  84.6%
                                 Number of contigs > 100K nt       4321  53.9%
                                   Number of contigs > 1M nt        653   8.1%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size     315273
                                          Median contig size     121491
                                           N50 contig length     826534
                                            L50 contig count        887
                                                   contig %A      29.20
                                                   contig %C      20.80
                                                   contig %G      20.81
                                                   contig %T      29.18
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
```
