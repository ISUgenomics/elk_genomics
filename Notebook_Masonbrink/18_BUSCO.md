# Assess genome completeness using BUSCO, benchmarking universal single copy orthologs

### Busco3 genome mode
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

### BUSCO3 genome mode Full run

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
### BUSCO3 protein mode Full run
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

### BUSCO 4 Genome mode
```
/work/GIF/remkv6/Elk/31_busco/02_Busco4


cp -rf /work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/05_pseudomolecule/config/spec
ies/CervusCanadensis3/ ~/.conda/pkgs/augustus-3.2.3-boost1.60_0/config/species/.

ml miniconda3; source activate busco4; export AUGUSTUS_CONFIG_PATH=/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/05_pseudomolecule/config;export BUSCO_CONFIG_FILE=/home/remkv6/.conda/envs/busco4/config/config.ini ; busco -i FinalGenomePilonReducedSoftMaskedRecode.fa --autolineage-euk -o GenomeBUSCO4 -m geno -c 15 -f --augustus_species CervusCanadensis3

```

### BUSCO 4 Protein Mode Final Primary Annotations
```

ln -s ../../../24_mikado/01_mikado2/FinalGenePrediction.proteins.fasta
ml cdbfasta
cdbfasta FinalGenePrediction.proteins.fasta

#get only primary isoforms to prevent artifactual duplicated buscos
less FinalGenePrediction.proteins.fasta |awk '{print $1}' |grep ">" |sed 's/>//g' |grep "\.1" |cdbyank FinalGenePrediction.proteins.fasta.cidx >PrimaryIsoformsFinalGenePrediction.proteins.fasta

ml miniconda3; source activate busco4; export BUSCO_CONFIG_FILE=/home/remkv6/.conda/envs/busco4/config/config.ini ; busco -i FinalGenomePilonReducedSoftMaskedRecode.fa --autolineage-euk -o GenomeBUSCO4 -m prot -c 15 -f --augustus_species CervusCanadensis3


INFO:   Running BUSCO using lineage dataset eukaryota_odb10 (eukaryota, 2019-11-20)
INFO:   ***** Run HMMER on gene sequences *****
INFO:   Running 255 job(s) on hmmsearch
INFO:   [hmmsearch]     26 of 255 task(s) completed
INFO:   [hmmsearch]     51 of 255 task(s) completed
INFO:   [hmmsearch]     77 of 255 task(s) completed
INFO:   [hmmsearch]     102 of 255 task(s) completed
INFO:   [hmmsearch]     128 of 255 task(s) completed
INFO:   [hmmsearch]     153 of 255 task(s) completed
INFO:   [hmmsearch]     179 of 255 task(s) completed
INFO:   [hmmsearch]     204 of 255 task(s) completed
INFO:   [hmmsearch]     230 of 255 task(s) completed
INFO:   [hmmsearch]     255 of 255 task(s) completed
INFO:   Results:        C:77.2%[S:68.2%,D:9.0%],F:5.5%,M:17.3%,n:255

```

### BUSCO 4 on Protein mode on Initial Primary annotations
```
#/work/GIF/remkv6/Elk/31_busco/02_Busco4/02_ProteinBusco4MikadoRound1NoFilter


ln -s ../../../24_mikado/01_mikado2/FinalGenePrediction.proteins.fasta
ml cdbfasta
cdbfasta FinalGenePrediction.proteins.fasta

#get only primary isoforms to prevent artifactual duplicated buscos
less FinalGenePrediction.proteins.fasta |awk '{print $1}' |grep ">" |sed 's/>//g' |grep "\.1" |cdbyank FinalGenePrediction.proteins.fasta.cidx >PrimaryIsoformsFinalGenePrediction.proteins.fasta


ml miniconda3; source activate busco4; export BUSCO_CONFIG_FILE=/home/remkv6/.conda/envs/busco4/config/config.ini ; busco -i PrimaryIsoformsMikado.proteins.fasta --auto-lineage-euk -o ProtBUSCO4 -m prot -c 15 -f --augustus_species CervusCanadensis3


INFO:   Running BUSCO using lineage dataset eukaryota_odb10 (eukaryota, 2019-11-20)
INFO:   ***** Run HMMER on gene sequences *****
INFO:   Running 255 job(s) on hmmsearch
INFO:   [hmmsearch]     26 of 255 task(s) completed
INFO:   [hmmsearch]     51 of 255 task(s) completed
INFO:   [hmmsearch]     77 of 255 task(s) completed
INFO:   [hmmsearch]     102 of 255 task(s) completed
INFO:   [hmmsearch]     128 of 255 task(s) completed
INFO:   [hmmsearch]     153 of 255 task(s) completed
INFO:   [hmmsearch]     179 of 255 task(s) completed
INFO:   [hmmsearch]     204 of 255 task(s) completed
INFO:   [hmmsearch]     230 of 255 task(s) completed
INFO:   [hmmsearch]     255 of 255 task(s) completed
INFO:   Results:        C:91.8%[S:80.0%,D:11.8%],F:6.3%,M:1.9%,n:255

```

### BUSCO4 on Red Deer proteins
```
#/work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4

# transferrred from ceres
GCA_002197005.1_CerEla1.0_protein.faa
cp buscoprot_0.sub  RedDeerbuscoprot_0.sub

ml miniconda3; source activate busco4; export BUSCO_CONFIG_FILE=/home/remkv6/.conda/envs/busco4/config/config.ini ; busco -i GCA_002197005.1_CerEla1.0_protein.faa --auto-lineage-euk -o RedDeerProtBUSCO4 -m prot -c 15 -f --augustus_species CervusCanadensis3


```
### Busco4 on cattle proteins
```
#/work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4
wget ftp://ftp.ensembl.org/pub/release-99/fasta/bos_taurus/pep/Bos_taurus.ARS-UCD1.2.pep.all.fa.gz
gunzip Bos_taurus.ARS-UCD1.2.pep.all.fa.gz
cp buscoprot_0.sub  Cowbuscoprot_0.sub

ml miniconda3; source activate busco4; export BUSCO_CONFIG_FILE=/home/remkv6/.conda/envs/busco4/config/config.ini ; busco -i Bos_taurus.ARS-UCD1.2.pep.all.fa --auto-lineage-euk -o Bos_taurus.ProtBUSCO4 -m prot -c 15 -f --augustus_species CervusCanadensis3
```

### BUSCO4 on mikado/braker merged elk proteins
```
#/work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4

#copied from nova
gffmergeElkGenesVHEJ_proteins.fasta

cp buscoprot_0.sub  gffmergeElkbuscoprot_0.sub

ml miniconda3; source activate busco4; export BUSCO_CONFIG_FILE=/home/remkv6/.conda/envs/busco4/config/config.ini ; busco -i gffmergeElkGenesVHEJ_proteins.fasta --auto-lineage-euk -o mikadobrakermergeProtBUSCO4 -m prot -c 15 -f --augustus_species CervusCanadensis3
```

### BUSCO4 on mikado/braker proteins reduced by bos/red deer gffs
```
#/work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4

#copied from nova
Bos_ReddeerReductionVHEJ_proteins.fasta

cp buscoprot_0.sub Bos_ReddeerReductionbuscoprot_0.sub

```

### Investigate BUSCO proteins by aligning proteins to genome
```

#Eukaryota dataset
#/work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4/01_AlignBuscoProteinsEuk
grep -c ">" /work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4/busco_downloads/lineages/eukaryota_odb10/ancestral
255
ln -s /work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4/busco_downloads/lineages/eukaryota_odb10/ancestral
ln -s ../../FinalGenomePilonReducedSoftMaskedRecode.fa

echo "ml miniconda3; source activate Genomethreader;gth  -genomic FinalGenomePilonReducedSoftMaskedRecode.fa -protein ancestral -skipalignmentout -gff3out -o EukAncestralBusco.gff3 -force" >gth.sh


#cetartiodactyla dataset
#/work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4/02_AlignBuscoProteinsCer
grep -c ">" /work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4/busco_downloads/lineages/cetartiodactyla_odb10/ancestral
13335
ln -s /work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4/busco_downloads/lineages/cetartiodactyla_odb10/ancestral
ln -s ../../FinalGenomePilonReducedSoftMaskedRecode.fa

ml miniconda3; source activate Genomethreader;gth  -genomic FinalGenomePilonReducedSoftMaskedRecode.fa -protein ancestral -skipalignmentout -gff3out -o CetAncestralBusco.gff3 -force

```
### Run Busco on High Confidence genes
```
#/work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4


ml miniconda3; source activate busco4; export BUSCO_CONFIG_FILE=/home/remkv6/.conda/envs/busco4/config/config.ini ; busco -i ProperCodingHighConfidenceBos_ReddeerReductionVHEJ_proteins.fasta --auto-lineage-euk -o ProperCodingHighConfidenceBos_ReddeerReductionVHEJ_proteinProtBUSCO4 -m prot -c 15 -f --augustus_species CervusCanadensis3


--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:97.7%[S:41.2%,D:56.5%],F:1.2%,M:1.1%,n:255     |
|249    Complete BUSCOs (C)                       |
|105    Complete and single-copy BUSCOs (S)       |
|144    Complete and duplicated BUSCOs (D)        |
|3      Fragmented BUSCOs (F)                     |
|3      Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------

--------------------------------------------------
|Results from dataset cetartiodactyla_odb10       |
--------------------------------------------------
|C:92.1%[S:40.2%,D:51.9%],F:1.7%,M:6.2%,n:13335   |
|12285  Complete BUSCOs (C)                       |
|5367   Complete and single-copy BUSCOs (S)       |
|6918   Complete and duplicated BUSCOs (D)        |
|227    Fragmented BUSCOs (F)                     |
|823    Missing BUSCOs (M)                        |
|13335  Total BUSCO groups searched               |
--------------------------------------------------

```
