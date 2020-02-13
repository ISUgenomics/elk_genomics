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
```
