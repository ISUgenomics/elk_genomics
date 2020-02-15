#  EDTA seems like the better program for repeat prediction than repeatmodeler

This will be used to filter gene models
```
#/work/GIF/remkv6/Elk/30_EDTA

ln -s ../../24_mikado/FinalGenomePilonReducedSoftMaskedRecode.fa

ml miniconda3; source activate EDTA; cd /work/GIF/remkv6/Elk/30_EDTA/EDTA; ./EDTA.pl --genome FinalGenomePilonReducedSoftMaskedRecode.fa --threads 16 --overwrite 1 --anno 1 --sensitive 1
```

### Split and run each chromosome separately
```

fasta-splitter.pl --n-parts 36 FinalGenomePilonReducedSoftMaskedRecode.fa
for f in *part*; do echo "ml miniconda3; source activate EDTA; cd /work/GIF/remkv6/Elk/30_EDTA/EDTA; ./EDTA.pl --genome "$f" --threads 16 --overwrite 1 --anno 1 --sensitive 1";done >>edta.sh


awk '{if($4>$5) {print $1,$5,$4} else {print $1,$4,$5}}' *anno.gff |awk 'substr($1,1,1)!="#"' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools merge -i - |awk '{print $3-$2}' |summary.sh
Total:  636,339,795
Count:  1,630,819
Mean:   390
Median: 234
Min:    2
Max:    103,396
cat *anno.gff >AllEDTARepeatAnnotations.gff


for f in *sum; do awk 'NR>5 && NR<25' $f;done |grep -v -i -e "-" -e "=" -e"total" -e "DNA" -e "LTR" >>SUM.db

awk '{print $1}' SUM.db |sort|uniq|sed 's/ //g' |sed '/^$/d' |awk '{print $0}' |while read line; do grep  "$line" SUM.db;done |awk '{print $1}' |sort|uniq -c |sort -k1,1nr |less
sort -n | awk ' $1 ~ /^[0-9]*(\.[0-9]*)?$/ { a[c++] = $1; sum += $1; } END { printf ("Total:\t""%'"'"'d\n", sum)}



awk '{print $1}' SUM.db |sort|uniq|sed 's/ //g' |sed '/^$/d' |awk '{print $0}' |while read line; do grep  "$line" SUM.db;done |awk '{print $1}' |sort |uniq >>sum1
awk '{print $1}' SUM.db |sort|uniq|sed 's/ //g' |sed '/^$/d' |awk '{print $0}' |while read line; do grep  "$line" SUM.db;done |awk '{print $1}' |sort|uniq |while read line; do grep $line SUM.db | awk '{print $2}' |sort -n | awk ' $1 ~ /^[0-9]*(\.[0-9]*)?$/ { a[c++] = $1; sum += $1; } END { printf ("%'"'"'d\n", sum)}' ;done >>sum2
[remkv6@condo2017 EDTA]$ awk '{print $1}' SUM.db |sort|uniq|sed 's/ //g' |sed '/^$/d' |awk '{print $0}' |while read line; do grep  "$line" SUM.db;done |awk '{print $1}' |sort|uniq |while read line; do grep $line SUM.db | awk '{print $3}' |sort -n | awk ' $1 ~ /^[0-9]*(\.[0-9]*)?$/ { a[c++] = $1; sum += $1; } END { printf ("%'"'"'d\n", sum)}' ;done  >>sum3
 paste sum1 sum2 sum3 >EDTA.sum

2526613007
 ###############################################################################
Class                  Count        bpMasked    %masked
 Copia   950     375,186
 DTA     75,223  13,578,319
 DTC     1,060,611       221,454,042
 DTH     128,084 20,954,105
 DTM     633,261 133,727,925
 DTT     101,761 16,676,500
 Gypsy   62,324  26,383,325
 Helitron        127,314 45,359,518
 Repeat  0       0
 unknown 629,492 174,212,026

 #################################################################################

 #create the gff conglomerate
cat *TEanno.gff >AllAnno.gff
cat *TElib.fa >AllTElib.fa
cat *MAKER.masked >FinalGenomePilonReducedSoftMaskedRecodeEDTAMasked.fa


#How many Genes have a cds that overlaps with a repeat?
bedtools intersect -wo -b AllAnno.gff -a <(awk '$3=="CDS"' ../../24_mikado/01_mikado2/mikado.loci.gff3 ) |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |sed 's/\./\t/2' |awk '{print $1}' |sort|uniq|wc
   3428    3428   47452


#Genes to remove
bedtools intersect -wo -b AllAnno.gff -a <(awk '$3=="CDS"' ../../24_mikado/01_mikado2/mikado.loci.gff3 ) |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |sed 's/\./\t/2' |awk '{print $1}' |sort|uniq >GenesToRemove.list
```
