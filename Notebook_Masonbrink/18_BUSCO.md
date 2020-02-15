# Assess genome completeness using BUSCO, benchmarking universal single copy orthologs

### Busco genome mode
```
run_BUSCO.py -i FinalGenomePilonReduced.fa  -l /home/rick.masonbrink/elk_bison_genomics/Masonbrink/25_BUSCO/mammalia_odb9 -o test1 -m geno -c 39  -f --long --limit 6  --augustus_parameters '--AUGUSTUS_CONFIG_PATH=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/25_BUSCO/config'

# Busco failed to use augustus in the second round, so this is a rudimentary measure.  However, a second run on condo with the same version ran augustus, and got a lower score.  So perhaps this is the best I can get.

INFO    C:92.0%[S:89.7%,D:2.3%],F:4.0%,M:4.0%,n:4104
INFO    3776 Complete BUSCOs (C)
INFO    3680 Complete and single-copy BUSCOs (S)
INFO    96 Complete and duplicated BUSCOs (D)
INFO    166 Fragmented BUSCOs (F)
INFO    162 Missing BUSCOs (M)
INFO    4104 Total BUSCO groups searched


```

### BUSCO genome mode Full run

```
#/work/GIF/remkv6/Elk/31_busco

#version 3.0.2
ln -s ../../24_mikado/FinalGenomePilonReducedSoftMaskedRecode.fa

ml miniconda2; source activate busco; export AUGUSTUS_CONFIG_PATH=/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/05_pseudomolecule/config; export ; run_BUSCO.py -i FinalGenomePilonReducedSoftMaskedRecode.fa -l /work/GIF/remkv6/Elk/31_busco/mammalia_odb9 -o GenomeBUSCO -m geno -c 15 -s CervusCanadensis3 -f

INFO    Results:
INFO    C:91.3%[S:89.4%,D:1.9%],F:4.5%,M:4.2%,n:4104
INFO    3748 Complete BUSCOs (C)
INFO    3670 Complete and single-copy BUSCOs (S)
INFO    78 Complete and duplicated BUSCOs (D)
INFO    185 Fragmented BUSCOs (F)
INFO    171 Missing BUSCOs (M)
INFO    4104 Total BUSCO groups searched

```
### BUSCO protein mode Full run
```
#/work/GIF/remkv6/Elk/31_busco/01_peptideBusco

#version 3.0.2

ln -s ../../24_mikado/01_mikado2/PrimaryIsoformsMikado.proteins.fasta

ml miniconda2; source activate busco; export AUGUSTUS_CONFIG_PATH=/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/05_pseudomolecule/config; export ; run_BUSCO.py -i PrimaryIsoformsMikado.proteins.fasta -l /work/GIF/remkv6/Elk/31_busco/mammalia_odb9 -o GenomeBUSCO -m prot -c 15 -s CervusCanadensis3 -f

INFO    Results:
INFO    C:79.4%[S:77.0%,D:2.4%],F:5.8%,M:14.8%,n:4104
INFO    3259 Complete BUSCOs (C)
INFO    3159 Complete and single-copy BUSCOs (S)
INFO    100 Complete and duplicated BUSCOs (D)
INFO    239 Fragmented BUSCOs (F)
INFO    606 Missing BUSCOs (M)
INFO    4104 Total BUSCO groups searched


```

### BUSCO 4
```
ml miniconda3; source activate busco4; export BUSCO_CONFIG_FILE=/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/05_pseudomolecule/config; busco -i FinalGenomePilonReducedSoftMaskedRecode.fa --autolineage-euk -o GenomeBUSCO -m geno -c 15 -f --augustus_species CervusCanadensis3

```
