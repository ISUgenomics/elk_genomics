#  EDTA seems like the better program for repeat prediction than repeatmodeler

This will be used to filter gene models
```
#/work/GIF/remkv6/Elk/30_EDTA

ln -s ../../24_mikado/FinalGenomePilonReducedSoftMaskedRecode.fa

ml miniconda3; source activate EDTA; cd /work/GIF/remkv6/Elk/30_EDTA/EDTA; ./EDTA.pl --genome FinalGenomePilonReducedSoftMaskedRecode.fa --threads 16 --overwrite 1 --anno 1 --sensitive 1

```
