# Need final counts of repetitive elements in the genome using repeatmasker


### setup
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/22_RepeatMasker2

#softlink consensii file from first repeatmasker/modeler run
ln -s ../03_repeatmodeler/consensi.fa.classified

ln -s ../27_RenameScaffolds/02_minimap/FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa

sh runRepeatModeler.sh FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa

#runRepeatModeler.sh
###########################################################################################################################################################################
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
module purge

ml repeatmasker/4.1.0
#module load repeatmodeler/1.0.8
#DATABASE="$(basename ${GENOME%.*}).DB"
#BuildDatabase -name ${DATABASE} -engine ncbi ${GENOME}
#RepeatModeler -database ${DATABASE}  -engine ncbi -pa 16
#ln -s $(find $(pwd) -name "consensi.fa.classified")
RepeatMasker -pa 40 -gff -lib consensi.fa.classified ${GENOME}
###########################################################################################################################################################################
```

### Results
```
==================================================
file name: FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa
sequences:           185
total length: 2526613007 bp  (2524723970 bp excl N/X-runs)
GC level:         41.62 %
bases masked:  950095574 bp ( 37.60 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements       3074508    871921964 bp   34.51 %
   SINEs:           398464     59471874 bp    2.35 %
   Penelope              0            0 bp    0.00 %
   LINEs:           2378908    727023500 bp   28.77 %
    CRE/SLACS            0            0 bp    0.00 %
     L2/CR1/Rex      84489     13232972 bp    0.52 %
     R1/LOA/Jockey       0            0 bp    0.00 %
     R2/R4/NeSL          0            0 bp    0.00 %
     RTE/Bov-B      1041248    326959230 bp   12.94 %
     L1/CIN4        1220946    373830151 bp   14.80 %
   LTR elements:    297136     85426590 bp    3.38 %
     BEL/Pao             0            0 bp    0.00 %
     Ty1/Copia           0            0 bp    0.00 %
     Gypsy/DIRS1        59       175262 bp    0.01 %
       Retroviral   297077     85251328 bp    3.37 %

DNA transposons     237111     38319877 bp    1.52 %
   hobo-Activator   125493     19695825 bp    0.78 %
   Tc1-IS630-Pogo    51864     13034360 bp    0.52 %
   En-Spm                0            0 bp    0.00 %
   MuDR-IS905            0            0 bp    0.00 %
   PiggyBac              0            0 bp    0.00 %
   Tourist/Harbinger     0            0 bp    0.00 %
   Other (Mirage,        0            0 bp    0.00 %
    P-element, Transib)

Rolling-circles          0            0 bp    0.00 %

Unclassified:        57641     13210554 bp    0.52 %

Total interspersed repeats:   923452395 bp   36.55 %


Small RNA:          249385     40566928 bp    1.61 %

Satellites:           2946      2527124 bp    0.10 %
Simple repeats:     480790     19602725 bp    0.78 %
Low complexity:      87206      4259264 bp    0.17 %
==================================================
```


## Combine this value with edta using bedtools merge
```
#/work/gif/remkv6/Olsen/Elk/08_RenameAgain

#they have to have overlap to be merged.
cat  Repeatmasker.gff3 HeaderDashedtidyAllAnnoEDTA.gff |sort -k1,1V -k4,5n |bedtools merge -i - >edtaRepmaskerMerge.gff3

#check to make sure there is not stranded information
awk '$2>$3' edtaRepmaskerMerge.gff3|wc
      0       0       0

#total repeats repeatmasker gff + edta gff       
awk '{print $3-$2}' edtaRepmaskerMerge.gff3|summary.sh
Total:  960,683,749
Count:  2,828,366
Mean:   339
Median: 197
Min:    6
Max:    84,693

# get bp of repeat per chromosome
awk '{print $1}' edtaRepmaskerMerge.gff3 |sort|uniq|while read line; do echo "awk '\$1==\""$line"\"{print \$3-\$2}' edtaRepmaskerMerge.gff3|summary.sh >>RepeatsPerChromosomeSummaries.txt"; done >test.sh



```
