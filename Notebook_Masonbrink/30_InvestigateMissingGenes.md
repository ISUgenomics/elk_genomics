#  Need to characterize the genes that are presumed missing from the elk genome


### Andrew's files
```
#Location of files that Andrew created
#/work/gif/severin/Olsen/Elk2020/07_immuneGenes


#Proposed gene deletions -- locations not accurate
|AY644518_TRGJ1| 01| 102107 |
|KT723008_IGHJ2|01| 509129 |
|AC172685_IGHA| 01| 655709|
|IMGT000049_TRDC|01| 3301039|
|D16120_TRGJ2|02| 338 |
|AY149283_IGHJ1| 02 | 783 |
|AY2277782_TRAJ31| 02| 1725 |
|AY644517_TRG| 02 |167787 |
|NW_001494075_IGHJ1|03| 39638 |


```

### Identify Regions in the cattle genome that these genes are located
```
#/work/gif/remkv6/Olsen/Elk/09_MapRnaseq2Cattle

vi Missing_From_Elk.list
################################################################################
AY644518_TRGJ1
KT723008_IGHJ2
AC172685_IGHA
IMGT000049_TRDC
D16120_TRGJ2
AY149283_IGHJ1
AY2277782_TRAJ31
AY644517_TRGC4
NW_001494075_IGHJ1
#################################################################################


ln -s /work/gif/severin/Olsen/Elk2020/07_immuneGenes/bosImmune_2_Elk/bestBlast_output.txt ElkBBH.out
ln -s /work/gif/severin/Olsen/Elk2020/07_immuneGenes/bosImmune_2_Bos/bestBlast_output.txt CowBBH.out


#grep these gene names from the blast output and create a bed file

##################################################################################
IMGT000049_TRDC*01_Bos_taurus_Hereford_F_EX2_g,3231015..3231088_75_nt_1_+1_-1___25_AA_25+0=25_____      NC_037337.1     100.000 25      0       0       1       25      22353349        22353275        1.77e-08
        51.6    25      103308737       0/-2    NC_037337.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 10, ARS-UCD1.2, whole genome shotgun sequence  100
AY644517_TRGC4*02_Bos_taurus_F_EX2A_g,167787..167845_60_nt_1_+1_-1___20_AA_20+0=20_____ NC_037331.1     100.000 20      0       0       1       20      82742556        82742615        7.07e-05        41.2    20      120000601       0/3     NC_037331.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 4, ARS-UCD1.2, whole genome shotgun sequence   100
AY644518_TRGJ1-1*01_Bos_taurus_F_J-REGION_102107..102166_60_nt_3_______19_AA_19+0=19_____       NC_037331.1     100.000 19      0       0       1       19      50155487        50155431        4.40e-05        41.6    19      120000601       0/-3    NC_037331.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 4, ARS-UCD1.2, whole genome shotgun sequence   100
D16120_TRGJ2-1*02_Bos_taurus_(F)_J-REGION_338..396_59_nt_3_______19_AA_19+0=19_partial_in_3'___ NC_037331.1     94.737  19      1       0       1       19      50155487        50155431        8.07e-04        38.1    19      120000601       0/-3    NC_037331.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 4, ARS-UCD1.2, whole genome shotgun sequence   100
NW_001494075_IGHJ1-1*03_Bos_taurus_Hereford_ORF_J-REGION_39812..39874_63_nt_3_______20_AA_20+0=20___rev-compl_  NC_037348.1     100.000 20      0       0       1       20      208561  208502  1.60e-05        42.7    20      69862954        0/-1    NC_037348.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 21, ARS-UCD1.2, whole genome shotgun sequence  100
AY149283_IGHJ1-2*02_Bos_taurus_ORF_J-REGION_783..838_56_nt_2_______18_AA_18+0=18_____   NC_037348.1     94.444  18      1       0       1       18      208381  208328  1.49e-04        40.0    18      69862954
        0/-1    NC_037348.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 21, ARS-UCD1.2, whole genome shotgun sequence  100
NW_001494075_IGHJ1-2*03_Bos_taurus_Hereford_ORF_J-REGION_39638..39693_56_nt_2_______18_AA_18+0=18___rev-compl_  NC_037348.1     100.000 18      0       0       1       18      208381  208328  5.76e-05        41.2    18      69862954        0/-1    NC_037348.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 21, ARS-UCD1.2, whole genome shotgun sequence  100
KT723008_IGHJ2-1*01_Bos_taurus_Holstein_ORF_J-REGION_508948..509010_63_nt_3_______20_AA_20+0=20_____    NC_037348.1     100.000 20      0       0       1       20      208561  208502  1.60e-05        42.7    20      69862954        0/-1    NC_037348.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 21, ARS-UCD1.2, whole genome shotgun sequence  100
KT723008_IGHJ2-2*01_Bos_taurus_Holstein_ORF_J-REGION_509129..509184_56_nt_2_______18_AA_18+0=18_____    NC_037348.1     100.000 18      0       0       1       18      208381  208328  5.76e-05        41.2    18      69862954        0/-1    NC_037348.1 Bos taurus isolate L1 Dominette 01449 registration number 42190680 breed Hereford chromosome 21, ARS-UCD1.2, whole genome shotgun sequence  100

##################################################################################



This bed file was incomplete in that it included a gene found in elk, but not cattle.  I removed this
cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1==1 {print $2}'  |sed 's/\*/\\\*/1' |awk '{print $1}' |grep -f - ElkBBH.out |less
AY227782_TRAJ11*02_Bos_taurus_F_J-REGION_37141..37200_60_nt_3_______19_AA_19+0=19_____  Chromosome_06   89.474  19      2       0       1       19      20872867        20872811        4.51e-04     38.5    19      96780817        0/-1    Chromosome_06   100
KT723008_IGHA*01_Bos_taurus_Holstein_F_M_g,655709..655893_186_nt_1_+1_____62_AA_62+0=62_____    Chromosome_17   93.548  62      4       0       1       62      7599393 7599208 3.98e-11    61.2     62      65378136        0/-1    Chromosome_17   100
```
### Genes missing in Elk Restart
```
/work/gif/remkv6/Olsen/Elk/09_MapRnaseq2Cattle

This grabs only those in Cow, and cuts the names after the asterisk.  
cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1==2 {print $2}'  |sed 's/\*/\\\*/1' |awk '{print $1}' |grep -f - CowBBH.out |awk '{if($10>$9) {print $2,$9,$10,$1}else {print $2,$10,$9,$1}}' |sort -k1,1 -k2,3n |tr " " "\t" |uniq|sed 's/\*/\t/g' |awk '{print $1,$2,$3,$4}' |less


NC_037331.1 50155431 50155487 AC172685_TRGJ2-1
NC_037331.1 50155431 50155487 AY644518_TRGJ1-1
NC_037331.1 50155431 50155487 D16118_TRGJ2-1
NC_037331.1 50155431 50155487 D16120_TRGJ2-1
NC_037331.1 82699801 82699857 D13648_TRGJ3-1
NC_037331.1 82707936 82707992 AY644517_TRGC3
NC_037331.1 82742556 82742615 AY644517_TRGC4
NC_037337.1 22272459 22272521 IMGT000049_TRAJ2
NC_037337.1 22273005 22273061 IMGT000049_TRAJ3
NC_037337.1 22277164 22277220 IMGT000049_TRAJ5
NC_037337.1 22278352 22278411 IMGT000049_TRAJ6
NC_037337.1 22280478 22280534 IMGT000049_TRAJ8-1
NC_037337.1 22283266 22283322 IMGT000049_TRAJ8-2
NC_037337.1 22284736 22284792 IMGT000049_TRAJ11
NC_037337.1 22285297 22285353 IMGT000049_TRAJ12
NC_037337.1 22290540 22290599 IMGT000049_TRAJ17
NC_037337.1 22292039 22292095 IMGT000049_TRAJ19
NC_037337.1 22297656 22297712 AY227782_TRAJ25
NC_037337.1 22297656 22297712 IMGT000049_TRAJ25
NC_037337.1 22300122 22300178 IMGT000049_TRAJ27
NC_037337.1 22301761 22301817 IMGT000049_TRAJ29
NC_037337.1 22304918 22304974 AY227782_TRAJ31
NC_037337.1 22307467 22307520 IMGT000049_TRAJ33
NC_037337.1 22308157 22308210 IMGT000049_TRAJ34
NC_037337.1 22309187 22309240 IMGT000049_TRAJ35
NC_037337.1 22312866 22312925 IMGT000049_TRAJ38
NC_037337.1 22315617 22315670 IMGT000049_TRAJ40
NC_037337.1 22317784 22317846 IMGT000049_TRAJ42
NC_037337.1 22321231 22321290 IMGT000049_TRAJ46
NC_037337.1 22323995 22324054 IMGT000049_TRAJ48
NC_037337.1 22324885 22324938 IMGT000049_TRAJ49
NC_037337.1 22334740 22334796 IMGT000049_TRAJ56
NC_037337.1 22335369 22335428 IMGT000049_TRAJ57
NC_037337.1 22353275 22353349 IMGT000049_TRDC
NC_037348.1 59222 59287 KT723008_IGHD
NC_037348.1 107140 107256 KT723008_IGHD1-3
NC_037348.1 208328 208381 AY149283_IGHJ1-2
NC_037348.1 208328 208381 KT723008_IGHJ1-2
NC_037348.1 208328 208381 KT723008_IGHJ2-2
NC_037348.1 208328 208381 NW_001494075_IGHJ1-2
NC_037348.1 208502 208561 KT723008_IGHJ2-1
NC_037348.1 208502 208561 NW_001494075_IGHJ1-1




#What are the fasta sequences to the hits in cattle so I can map them to elk?
cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1==2 {print $2}'  |sed 's/\*/\\\*/1' |awk '{print $1}' |grep -f - CowBBH.out |awk '{if($10>$9) {print $2,$9,$10,$1}else {print $2,$10,$9,$1}}' |sort -k1,1 -k2,3n |tr " " "\t" |uniq|sed 's/\*/\t/g' |awk '{print $1,$2,$3}' |tr " " "\t" |bedtools merge -i - |bedtools getfasta -bed - -fi GCF_002263795.1_ARS-UCD1.2_genomic.fna -fo MissingGenesInElkCattleNucSequence.fasta

>NC_037331.1:50155431-50155487
GGAGGAGTTACTATGATGTTAGCTCCTTCTCCAAATATCTTCCTCCAGCCTGAGCT
>NC_037331.1:82699801-82699857
GTGCAGGCTGGAAGAAGATATTTGGAAAAGCAACTGAGCTCATAGTAGCTCCTCCT
>NC_037331.1:82707936-82707992
TTGTCTCTTCAGTTGTCACTGCTACTAAACCTCCAAATGATGGTTTGAAGGATAAA
>NC_037331.1:82742556-82742615
TTGTCAGCTCAATTGTCCCTACCACTGAGTCTCCCAGTGACTGTTTAAACCATGACAGC
>NC_037337.1:22272459-22272521
GATATTATGGACACTTGGGTTCCTTTCCCAAAAGTGAGTTTACCGATCGCTCCTCCTGTATT
>NC_037337.1:22273005-22273061
GGATGGACACTCAGTCTGGTCCCTGCTCCGAAGTTTAACTTATCAGTGCCGAATCT
>NC_037337.1:22277164-22277220
GGATGCACTCGGAGTCTTGTTCCACTCCCAAAAGTGAGTGTTCTGCTGCCTGTGTT
>NC_037337.1:22278352-22278411
GGGTGAACAACAAGCTTGGTCCCTGTTCCAAATATAAATCCATACTTTAGTCTTGATAC
>NC_037337.1:22280478-22280534
GGGTTGATCAAAAGTTGGGTGCCAGTTCCAAATACGAATTTCTGATAACCTGTGCT
>NC_037337.1:22283266-22283322
GGGTTGATTAAAAGTTGGGTGCCAGTTCCAAATGTGAGTTTCTGATAACCTGTGTT
>NC_037337.1:22284736-22284792
GGAAAGACAAGAAGCACAGTGCCCTTTCCAAAAGTAAGTGTGTTGTATCCTGAATT
>NC_037337.1:22285297-22285353
GGCCTGACCAGCAGTCTAGTCCCACTTCCGAAGGTCCATGTATAGCCTCCATCCTT
>NC_037337.1:22290540-22290599
AGCTTGACTATCAGCCTGGTTCCTTCTCCAAAGGTTAGCTTGGTCCCTGCAGTGGTGGT
>NC_037337.1:22292039-22292095
GGATTGTCATTGTGTTTGGATCCTTTTCCAAAGGTGAACTTGCGAAAACTCTGATA
>NC_037337.1:22297656-22297712
GGCTTGACAAGCAGCCTTGTCCCCTTCCCAAAGACGAGGGAGAAGCCTTGTCCTTC
>NC_037337.1:22300122-22300178
GGCTTCACAGTGAGCGCAGTCCCGTCCCCAAAGGTTAATTTGCCTGTACCGGTGTT
>NC_037337.1:22301761-22301817
GGAGTCACAGCAAGTCTTGTGCCTTTTCCAAAGACAAGTTGCCTGTTTCCTGAATC
>NC_037337.1:22304918-22304974
GGCTTTACCACCACCTGggttccagttccaaagaagatTCTGGAATTGTCAGTCCC
>NC_037337.1:22307467-22307520
GGCTTTATAATTAGTTTGGTCCCAGAGCCCCAGATCCACTGATAGTTGCCATC
>NC_037337.1:22308157-22308210
GGAAAAACTTGTAATCTGGTTCCAGCCCCAAAGATGAGTTTGTCCCTGTTGGA
>NC_037337.1:22309187-22309240
GGCATAACGATCACTTGAGTGCCAGACCCCAAATGCAGCACATTCCCAAAGTT
>NC_037337.1:22312866-22312925
GGATTTACTGCCAGACTTGTCCCCAGTCCCCAAATCAGCTTACGGTTGTTGCCAGTATT
>NC_037337.1:22315617-22315670
GTTAAAACCTGCAGCCTAGTGCCCGCTCCAAAGACGTATTTGTAGTTTCCTGA
>NC_037337.1:22317784-22317846
GGTTTAACAGAGACCATAGTGCCTTTTCCAAAGATGAGCTGTCCTTGGCCGCTGCCAACATT
>NC_037337.1:22321231-22321290
GGCCTCACTGCTAAACGCGTCCCGGTCCCAAAAGTCAGCCTGTCTCCGCTGCCACTCTT
>NC_037337.1:22323995-22324054
GCTGTGATGGTGAGTCTAGTTCCTGTTCCAAAGTTTAATTGACTGCCTTGGTAGTTAGA
>NC_037337.1:22324885-22324938
gGAATGACTGTCAAACTTGTCCCTCTCCCAAAATAGTTCTGGCCGTAGCTGTT
>NC_037337.1:22334740-22334796
GGTCTAACATTCAGGATTGTTCCTTTTCCAAATGTCAGCTTATTATTGGAGTATCA
>NC_037337.1:22335369-22335428
GGGTTTACTGTCAGTTTCGTTCCCTTTCCAAAGAGGAGTCTTTCAGATCCGCCCTGAGT
>NC_037337.1:22353275-22353349
TGGGGCTCTTGGCAGgtcactggagcttcagcttttgtGCTGTTTTCATATGCCATCGGTTTTGGAGTTGTTTC
>NC_037348.1:59222-59287
CTGCAGGGTTCCTTTGGACACGCCAGTGTCTGGGCTTCTGTCTGTCTGTGGCTGGTGGTCAGTGC
>NC_037348.1:107140-107256
GGGAGCCCAGTGCCCCTGGGGATGCATTGGCCAGCTCCACACCTGTGTGCGGGTCAGACTTTATGTCAGGGCCTGAGTCACTGTGGGTGTAGCAGTAACCATCATCACGATAGTCT
>NC_037348.1:208328-208381
CCCAAGGACACGGTGACCGGGGTGCGCTGGCCCCAGAGATCCATGTCCCAGCA
>NC_037348.1:208502-208561
GATTCAGCTGAGGAGACGGTGCCCAGGGCAGCCTGGCCACAGAGATGGAAGTCAGCATA

```



## Error in getting the genes restart is above.
```
###########################Created bed file
######################################cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1<3' |grep -f Missing_From_Elk.list - |awk '{print $2}'  |sed 's/\*/\\*/g' |grep -f - CowBBH.out |awk '{if($9>$10) {print $2,$10,$9,$1} else {print $2,$9,$10,$1} }' |sed 's/\*/\t/g' |sort -k1,1V -k2,3n |tr " " "\t" |cut -f 1,2,3,4 |less

################################################################################
#####NC_037331.1     50155431        50155487        AY644518_TRGJ1-1
#####NC_037331.1     50155431        50155487        D16120_TRGJ2-1
#####NC_037331.1     82742556        82742615        AY644517_TRGC4
#####NC_037337.1     22353275        22353349        IMGT000049_TRDC
#####NC_037348.1     208328  208381  AY149283_IGHJ1-2
#####NC_037348.1     208328  208381  KT723008_IGHJ2-2
#####NC_037348.1     208328  208381  NW_001494075_IGHJ1-2
#####NC_037348.1     208502  208561  KT723008_IGHJ2-1
#####NC_037348.1     208502  208561  NW_001494075_IGHJ1-1
################################################################################


#Merged bed file
#####cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1<3' |grep -f Missing_From_Elk.list - |awk '{print $2}'  |sed 's/\*/\\*/g' |grep -f - CowBBH.out |awk '{if($9>$10) {print $2,$10,$9,$1} else {print $2,$9,$10,$1} }' |sed 's/\*/\t/g' |sort -k1,1V -k2,3n |tr " " "\t" |cut -f 1,2,3,4 |bedtools merge -i - | less

#############################################################################
#####NC_037331.1     50155431        50155487  AY644518_TRGJ1-1, D16120_TRGJ2-1
#####NC_037331.1     82742556        82742615  AY644517_TRGC4
#####NC_037337.1     22353275        22353349  IMGT000049_TRDC
#####NC_037348.1     208328  208381  AY149283_IGHJ1-2, KT723008_IGHJ2-2,NW_001494075_IGHJ1-2
#####NC_037348.1     208502  208561  KT723008_IGHJ2-1, NW_001494075_IGHJ1-1
```


### Bedtools intersect with the genic gff of cattle
```
which mRNA/genes are these proteins associated with?
cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1<3' |grep -f Missing_From_Elk.list - |awk '{print $2}'  |sed 's/\*/\\*/g' |grep -f - CowBBH.out |awk '{if($9>$10) {print $2,$10,$9,$1} else {print $2,$9,$10,$1} }' |sed 's/\*/\t/g' |sort -k1,1V -k2,3n |tr " " "\t" |cut -f 1,2,3,4 |bedtools merge -i - | bedtools intersect -wo -a - -b ../06_Synteny/02_BosTaurus/GCF_002263795.1_ARS-UCD1.2_genomic.gff |awk '$6=="gene"' |less


######################################################################################################################
NC_037331.1     50155431        50155487        NC_037331.1     BestRefSeq%2CGnomon     gene    50141812        50162835        .       -       .       ID=gene-TARP;Dbxref=GeneID:100335800;Name=TARP;description=TCR gamma alternate reading frame protein;gbkey=Gene;gene=TARP;gene_biotype=protein_coding;gene_synonym=TCR,TRGC1,TRGJ1-2    56
NC_037331.1     82742556        82742615        NC_037331.1     Gnomon  gene    82720193        82746245        .       +       .       ID=gene-TRGC4;Dbxref=GeneID:509614;Name=TRGC4;gbkey=Gene;gene=TRGC4;gene_biotype=C_region       59
NC_037337.1     22353275        22353349        NC_037337.1     Gnomon  gene    22350419        22359327        .       -       .       ID=gene-LOC101908015;Dbxref=GeneID:101908015;Name=LOC101908015;gbkey=Gene;gene=LOC101908015;gene_biotype=C_region       74
NC_037337.1     22353275        22353349        NC_037337.1     Gnomon  gene    22263114        23332869        .       -       .       ID=gene-LOC100336282;Dbxref=GeneID:100336282;Name=LOC100336282;gbkey=Gene;gene=LOC100336282;gene_biotype=protein_coding;partial=true    74
######################################################################################################################

So it looks to be 4 genes, not nine.

These two did not intersect with a gene.
NC_037348.1     208328  208381  AY149283_IGHJ1-2, KT723008_IGHJ2-2,NW_001494075_IGHJ1-2
NC_037348.1     208502  208561  KT723008_IGHJ2-1, NW_001494075_IGHJ1-1

Investigated on https://bovinegenome.elsiklab.missouri.edu/ JBROWSE

Interprets to chromosome 21:208238..208561.  
This is clearly a transcribed locus according to isoseq from, and when viewing among the bam dense files for tissues, only appears to be expressed in nasal mucosa

NC_037348.1     208328  208381  AY149283_IGHJ1-2, KT723008_IGHJ2-2,NW_001494075_IGHJ1-2
NC_037348.1     208502  208561  KT723008_IGHJ2-1, NW_001494075_IGHJ1-1



```

### Extract genes and align
```
#essentialy this is just grabbing blast hits only found in cattle, identifying the region they map to in cattle, finding the corresponding gene in cattle, extracting the genomic sequence
cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1<3' |grep -f Missing_From_Elk.list - |awk '{print $2}'  |sed 's/\*/\\*/g' |grep -f - CowBBH.out |awk '{if($9>$10) {print $2,$10,$9,$1} else {print $2,$9,$10,$1} }' |sed 's/\*/\t/g' |sort -k1,1V -k2,3n |tr " " "\t" |cut -f 1,2,3,4 |bedtools merge -i - | bedtools intersect -wo -a - -b ../06_Synteny/02_BosTaurus/GCF_002263795.1_ARS-UCD1.2_genomic.gff |awk '$6=="gene"' |uniq| awk '{print $4":"$7"-"$8}' |while read line; do samtools faidx GCF_002263795.1_ARS-UCD1.2_genomic.fna $line; done >ExtractedMissingGenes.fasta

#Put it together with the elk genome, to make sure spurious reads do not map to the cattle genes.
cat ExtractedMissingGenes.fasta ../08_RenameAgain/01_Move2Box/CervusCanadensisGenome.fa >MissingGenesAddedToElkGenome.fasta


```

### Transfer modified genome to ceres and run expression
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression

#Create a gff that only recognizes these four geneic scaffolds.
ml miniconda
source activate bioawk
bioawk -c fastx '{print $name,"1",length($seq)}' MissingGenesAddedToElkGenome.fasta |head -n 4 |awk '{print $1,"NA","gene",$2,$3,"+",".",".","ID="$1"gene" }' |tr " " "\t" >MissingGenesAddedToElkGenome.gff

#softlink the rnaseq
Elk-spleen_S23_L003_R2_001.fastq -> ../16_RNAseq/lane3/Elk-spleen_S23_L003_R2_001.fastq
Elk-spleen_S23_L003_R1_001.fastq -> ../16_RNAseq/lane3/Elk-spleen_S23_L003_R1_001.fastq
ElkpscapLN_S22_L003_R2_001.fastq -> ../16_RNAseq/lane3/ElkpscapLN_S22_L003_R2_001.fastq
ElkpscapLN_S22_L003_R1_001.fastq -> ../16_RNAseq/lane3/ElkpscapLN_S22_L003_R1_001.fastq
Elk-muscle_S21_L003_R2_001.fastq -> ../16_RNAseq/lane3/Elk-muscle_S21_L003_R2_001.fastq
Elk-muscle_S21_L003_R1_001.fastq -> ../16_RNAseq/lane3/Elk-muscle_S21_L003_R1_001.fastq
Elk-Mes-LN_S24_L003_R2_001.fastq -> ../16_RNAseq/lane3/Elk-Mes-LN_S24_L003_R2_001.fastq
Elk-Mes-LN_S24_L003_R1_001.fastq -> ../16_RNAseq/lane3/Elk-Mes-LN_S24_L003_R1_001.fastq
Elk-lung_S26_L003_R2_001.fastq -> ../16_RNAseq/lane3/Elk-lung_S26_L003_R2_001.fastq
Elk-lung_S26_L003_R1_001.fastq -> ../16_RNAseq/lane3/Elk-lung_S26_L003_R1_001.fastq
Elk-kidney_S25_L003_R2_001.fastq -> ../16_RNAseq/lane3/Elk-kidney_S25_L003_R2_001.fastq
Elk-kidney_S25_L003_R1_001.fastq -> ../16_RNAseq/lane3/Elk-kidney_S25_L003_R1_001.fastq
Elk-spleen_S23_L004_R2_001.fastq -> ../16_RNAseq/lane4/Elk-spleen_S23_L004_R2_001.fastq
Elk-spleen_S23_L004_R1_001.fastq -> ../16_RNAseq/lane4/Elk-spleen_S23_L004_R1_001.fastq
ElkpscapLN_S22_L004_R2_001.fastq -> ../16_RNAseq/lane4/ElkpscapLN_S22_L004_R2_001.fastq
ElkpscapLN_S22_L004_R1_001.fastq -> ../16_RNAseq/lane4/ElkpscapLN_S22_L004_R1_001.fastq
Elk-muscle_S21_L004_R2_001.fastq -> ../16_RNAseq/lane4/Elk-muscle_S21_L004_R2_001.fastq
Elk-muscle_S21_L004_R1_001.fastq -> ../16_RNAseq/lane4/Elk-muscle_S21_L004_R1_001.fastq
Elk-Mes-LN_S24_L004_R2_001.fastq -> ../16_RNAseq/lane4/Elk-Mes-LN_S24_L004_R2_001.fastq
Elk-Mes-LN_S24_L004_R1_001.fastq -> ../16_RNAseq/lane4/Elk-Mes-LN_S24_L004_R1_001.fastq
Elk-lung_S26_L004_R2_001.fastq -> ../16_RNAseq/lane4/Elk-lung_S26_L004_R2_001.fastq
Elk-lung_S26_L004_R1_001.fastq -> ../16_RNAseq/lane4/Elk-lung_S26_L004_R1_001.fastq
Elk-kidney_S25_L004_R2_001.fastq -> ../16_RNAseq/lane4/Elk-kidney_S25_L004_R2_001.fastq
Elk-kidney_S25_L004_R1_001.fastq -> ../16_RNAseq/lane4/Elk-kidney_S25_L004_R1_001.fastq


#map and count
sh runFeatureCounts.sh Elk-kidney_S25_L003_R1_001.fastq Elk-kidney_S25_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-kidney_S25_L004_R1_001.fastq Elk-kidney_S25_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-lung_S26_L003_R1_001.fastq Elk-lung_S26_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-lung_S26_L004_R1_001.fastq Elk-lung_S26_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-Mes-LN_S24_L003_R1_001.fastq Elk-Mes-LN_S24_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-Mes-LN_S24_L004_R1_001.fastq Elk-Mes-LN_S24_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-muscle_S21_L003_R1_001.fastq Elk-muscle_S21_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-muscle_S21_L004_R1_001.fastq Elk-muscle_S21_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh ElkpscapLN_S22_L003_R1_001.fastq ElkpscapLN_S22_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh ElkpscapLN_S22_L004_R1_001.fastq ElkpscapLN_S22_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-spleen_S23_L003_R1_001.fastq Elk-spleen_S23_L003_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff
sh runFeatureCounts.sh Elk-spleen_S23_L004_R1_001.fastq Elk-spleen_S23_L004_R2_001.fastq /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression MissingGenesAddedToElkGenome.fasta MissingGenesAddedToElkGenome.gff


##############################################################################################################################################

#runFeatureCounts.sh
#note: I premade the database so it would not be made again and again
###############################################################################
#!/bin/bash

PROC=40
R1_FQ="$1"
R2_FQ="$2"
DBDIR="$3"
GENOME="$4"
GFF="$5"

#module load trimgalore
#module load py-setuptools/35.0.2-py2-hqrsjff
#trim_galore --paired ${R1_FQ} ${R2_FQ}

module load hisat2
#hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p ${PROC} -x ${GENOME%.*} -1 $R1_FQ -2 $R2_FQ -S ${R1_FQ%.*}.sam

module load samtools
samtools view --threads ${PROC} -b -o ${R1_FQ%.*}.bam ${R1_FQ%.*}.sam
samtools sort -m 3G -o ${R1_FQ%.*}_sorted.bam -T ${R1_FQ%.*}_temp --threads ${PROC} ${R1_FQ%.*}.bam

module load GIF2/java/1.8.0_25-b17
module load picard
java -jar /software/7/apps/picard/64/2.9.2/picard.jar CollectAlignmentSummaryMetrics  REFERENCE_SEQUENCE=${DBDIR}/${GENOME} INPUT=${R1_FQ%.*}_sorted.bam OUTPUT=${R1_FQ%.*}.bam_alignment.stats

module load subread
featureCounts -T ${PROC} -p -t gene -g ID -a ${GFF} -o ${R1_FQ%.*}_counts_genes.txt ${R1_FQ%.*}_sorted.bam
##############################################################################################################################################
```

### Alignment stats
```
for f in *stats; do echo "awk '\$1==\"PAIR\"{print FILENAME,\$7}' "$f" >>RNAseqAlignment.txt";done >stats.sh

Elk-kidney_S25_L003_R1_001.bam_alignment.stats 0.923363
Elk-kidney_S25_L004_R1_001.bam_alignment.stats 0.922779
Elk-Mes-LN_S24_L003_R1_001.bam_alignment.stats 0.907378
Elk-Mes-LN_S24_L004_R1_001.bam_alignment.stats 0.906394
Elk-muscle_S21_L003_R1_001.bam_alignment.stats 0.741857
Elk-muscle_S21_L004_R1_001.bam_alignment.stats 0.741641
ElkpscapLN_S22_L003_R1_001.bam_alignment.stats 0.934488
ElkpscapLN_S22_L004_R1_001.bam_alignment.stats 0.934037
Elk-spleen_S23_L003_R1_001.bam_alignment.stats 0.887191
```

# Expression Results
```

for f in *counts_genes.txt; do echo "<(cut -f 7 "$f" |grep -v \"#\" )" ;done |tr "\n" " " |sed 's/^/paste <(cut -f 1 Elk-spleen_S23_L004_R1_001_counts_genes.txt |grep -v "#" ) /g' |less

paste <(cut -f 1 Elk-spleen_S23_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-kidney_S25_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-kidney_S25_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-lung_S26_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-lung_S26_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-Mes-LN_S24_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-Mes-LN_S24_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-muscle_S21_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-muscle_S21_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 ElkpscapLN_S22_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 ElkpscapLN_S22_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-spleen_S23_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-spleen_S23_L004_R1_001_counts_genes.txt |grep -v "#" ) |less

################################################################################
Geneid  Elk-kidney_S25_L003_R1_001_sorted.bam   Elk-kidney_S25_L004_R1_001_sorted.bam   Elk-lung_S26_L003_R1_001_sorted.bam     Elk-lung_S26_L004_R1_001_sorted.bam     Elk-Mes-LN_S24_L003_R1_001_sorted.bam   Elk-Mes-LN_S24_L004_R1_001_sorted.bam   Elk-muscle_S21_L003_R1_001_sorted.bam   Elk-muscle_S21_L004_R1_001_sorted.bam   ElkpscapLN_S22_L003_R1_001_sorted.bam   ElkpscapLN_S22_L004_R1_001_sorted.bam   Elk-spleen_S23_L003_R1_001_sorted.bam   Elk-spleen_S23_L004_R1_001_sorted.bam
NC_037331.1:50141812-50162835gene       0       1       1       0       3       1       0       0       5       4       8       7
NC_037331.1:82720193-82746245gene       1       3       5       2       2       4       0       0       1       0       3       1
NC_037337.1:22350419-22359327gene       0       0       0       0       0       0       0       0       0       0       0       0
NC_037337.1:22263114-23332869gene       59      51      131     115     232     220     9       19      72      66      55      61
################################################################################

```
Only the bottom gene seems to have a significant number of reads, and likely there is a gene in that region conserved between cattle and elk.  I need finer depth windows.

## Investigate missing genes in the pacbio subreads

```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/33_MissingGenesPacbioAlignment

AllFastq.fastq -> ../06_blobtools/AllFastq.fastq
MissingGenesAddedToElkGenome.fasta -> ../32_MissingGeneExpression/MissingGenesAddedToElkGenome.fasta
MissingGenesAddedToElkGenome.gff -> ../32_MissingGeneExpression/MissingGenesAddedToElkGenome.gff

#asm20 allows for a 5% divergence in sequence, since I am mapping elk reads to cattle genes.
minimap2 -H -a -o MissingGenesAddedToElkGenome.sam -t 40 -x asm20 MissingGenesAddedToElkGenome.fasta AllFastq.fastq;ml samtools; samtools sort -o MissingGenesAddedToElkGenome.bam MissingGenesAddedToElkGenome.sam; ml subread; featureCounts -T 40 -p -t gene -g ID -a MissingGenesAddedToElkGenome.gff  -o MissingGenesAddedToElkGenome_counts_genes.txt MissingGenesAddedToElkGenome.bam

```
### Pacbio alignment rate
```

samtools flagstat MissingGenesAddedToElkGenome.sam
12887228 + 0 in total (QC-passed reads + QC-failed reads)
6977634 + 0 secondary
497574 + 0 supplementary
0 + 0 duplicates
10932108 + 0 mapped (84.83% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

this is ~13% lower than what i got using different settings with minimap, assuming it is because I told it -asm20, allowing for divergence.
```

### Pacbio Mapping Results
```
Geneid  Chr     Start   End     Strand  Length  MissingGenesAddedToElkGenome.sam
NC_037331.1:50141812-50162835gene       NC_037331.1:50141812-50162835   1       21024   +       21024   0
NC_037331.1:82720193-82746245gene       NC_037331.1:82720193-82746245   1       26053   +       26053   0
NC_037337.1:22350419-22359327gene       NC_037337.1:22350419-22359327   1       8909    +       8909    0
NC_037337.1:22263114-23332869gene       NC_037337.1:22263114-23332869   1       1069756 +       1069756 32
```
With the bottom gene being 100kb long, it is tough to say whether the gene is missing biologically.  This fourth gene might be missing from our assembly, but be there biologically.  

##  investigate regions these proteins map to in the cattle genome and 100bp surrounding them.
On Nova now
```
#/work/gif/remkv6/Olsen/Elk/09_MapRnaseq2Cattle

cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1<3' |grep -f Missing_From_Elk.list - |awk '{print $2}'  |sed 's/\*/\\*/g' |grep -f - CowBBH.out |awk '{if($9>$10) {print $2,$10,$9,$1} else {print $2,$9,$10,$1} }' |sed 's/\*/\t/g' |sort -k1,1V -k2,3n |tr " " "\t" |cut -f 1,2,3,4 |bedtools merge -i - |awk '{print $1":"$2-100"-"$3+100}' |while read line; do samtools faidx GCF_002263795.1_ARS-UCD1.2_genomic.fna $line; done >ExtractedMissingRegions.fasta

cat ExtractedMissingRegions.fasta ../08_RenameAgain/01_Move2Box/CervusCanadensisGenome.fa >MissingRegionsAddedToElkGenome.fasta
#Transfer to ceres
##############################################################################################
>NC_037331.1:50155331-50155587
CAAGAATGTGGCAGCCACATCAGGGAGAGCCCAGGTTGTTTCTACAAAGCATCTCTCAGT
ATTTTATGACTACAGTAagaatagaaaaatgatactTACCAGGAGGAGTTACTATGATGT
TAGCTCCTTCTCCAAATATCTTCCTCCAGCCTGAGCTGTCACAATGGTAAAAGCTCCTAC
AATAATCCAGCCTTTTTCTACCTTTATAAAGTGACTAGAAGTCAAAAACTTCATCCTGGC
AGTTTGATCTGAATTCT
>NC_037331.1:82742456-82742715
atatacatatatatgtattatacaataacacatgtataaacacacacaatgTCAAATAAT
GGGCACAAATAATTGACtaatacttttttcatttcttgtagTTGTCAGCTCAATTGTCCC
TACCACTGAGTCTCCCAGTGACTGTTTAAACCATGACAGCAGTAAGTTTTTGTAGATGAT
TGCCTTTATGTCATGCTTtagtcttttatttcctttctaatCTATTTATCTTTTGAACTT
CTCTGGCTAAAAGTACAAAA
>NC_037337.1:22353175-22353449
ttgttttcttttgaaatagaATGATTTGGGGCCCCATCTGACAATGACAGAGTTTGCTTT
AACCTCAATTCTGAAGCTGGCCCTTTGGCTTCAACTAACCTTGGGGCTCTTGGCAGgtca
ctggagcttcagcttttgtGCTGTTTTCATATGCCATCGGTTTTGGAGTTGTTTCTGGAA
CAAACAGCAATGGTCAGGGTTACAATGACAACTGAGAGGCTCAAGGGCTTTCAGGTTTTA
GAACTGGAGGTGACTTTAAGTAAGTAAACTAAGGA
>NC_037348.1:208228-208481
ACCCCGTTGGACCAGGACCCTACCTTGCTCAGGGCACGGTGGGGCTGCACGGACCCCAGG
GAACTGCCTGAGGCCCAGGCAGTGCTGCAAGGGAGACACTCCCCAAGGACACGGTGACCG
GGGTGCGCTGGCCCCAGAGATCCATGTCCCAGCAGCATGGCTTCCattctgttgcttctt
ttccACAAAAACACACCCTACCTCCTTGCTGGCCCTGACTCACATCCCAAGTCCACCCAA
GTCCCACTGTCCCC
>NC_037348.1:208402-208661
cttcttttccACAAAAACACACCCTACCTCCTTGCTGGCCCTGACTCACATCCCAAGTCC
ACCCAAGTCCCACTGTCCCCTGACTTGGACACTAGACACCAGATTCAGCTGAGGAGACGG
TGCCCAGGGCAGCCTGGCCACAGAGATGGAAGTCAGCATAGTCACAGTGGGCTGGCCCCG
TGCTGGGAGCACGGAAACCCACCCAGGGTCCCCAGGGAGCCCCCAGGAGCTGCAAGTCTG
TCTTCCTGAAGAcacttgcc
##############################################################################################

```


##### Moved to ceres
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression/01_Region

#grabbed genome from different location, but is same
cat MissinginElk_CattleRegions.fasta  ../../31_CallSNPS/CervusCanadensisGenome.fa > MissingRegionsAddedToElkGenome.fasta


 hisat2-build MissingRegionsAddedToElkGenome.fasta MissingRegionsAddedToElkGenome

 bioawk -c fastx '{print $name,"1",length($seq)}' MissinginElk_CattleRegions.fasta |head -n 4 |awk '{print $1,"NA","gene",$2,$3,"+",".",".","ID="$1"gene" }' |tr " " "\t" >MissingRegionsAddedToElkGenome.gff

paste <(ls -1 *_R1_* )  <(ls -1 *_R2_* ) |while read line; do echo "sh runFeatureCounts.sh "$line" /home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression/01_Region MissingRegionsAddedToElkGenome.fasta MissingRegionsAddedToElkGenome.gff"; done >RunPipeline.sh


Region expression

paste <(cut -f 1 Elk-spleen_S23_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-kidney_S25_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-kidney_S25_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-lung_S26_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-lung_S26_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-Mes-LN_S24_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-Mes-LN_S24_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-muscle_S21_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-muscle_S21_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 ElkpscapLN_S22_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 ElkpscapLN_S22_L004_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-spleen_S23_L003_R1_001_counts_genes.txt |grep -v "#" ) <(cut -f 7 Elk-spleen_S23_L004_R1_001_counts_genes.txt |grep -v "#" ) |less

#######################################################################################
Geneid  Elk-kidney_S25_L003_R1_001_sorted.bam   Elk-kidney_S25_L004_R1_001_sorted.bam   Elk-lung_S26_L003_R1_001_sorted.bam     Elk-lung_S26_L004_R1_001_sorted.bam     Elk-Mes-LN_S24_L003_R1_001_sorted.bam   Elk-Mes-LN_S24_L004_R1_001_sorted.bam   Elk-muscle_S21_L003_R1_001_sorted.bam   Elk-muscle_S21_L004_R1_001_sorted.bam   ElkpscapLN_S22_L003_R1_001_sorted.bam   ElkpscapLN_S22_L004_R1_001_sorted.bam   Elk-spleen_S23_L003_R1_001_sorted.bam   Elk-spleen_S23_L004_R1_001_sorted.bam
NC_037331.1:50155331-50155587gene       0       0       0       0       0       0       0       0       0       0       0       0
NC_037331.1:82742456-82742715gene       0       0       0       0       0       0       0       0       0       0       0       0
NC_037337.1:22353175-22353449gene       0       0       0       0       0       0       0       0       0       0       0       0
NC_037348.1:208228-208481gene   0       0       0       0       0       0       0       0       0       0       0       0
#######################################################################################
```

### Genomic structure of where missing genes should be.

```
NC_037331.1:50155331-50155587 --> Bt_Chr3  50155331  50155587
NC_037331.1:82742456-82742715 --> Bt_Chr3 82742456  82742715
NC_037337.1:22353175-22353449 --> Bt_Chr10  22353175  22353449
NC_037348.1:208228-208481 --> Bt_Chr21  208228  208481

vi MissingGeneRegionsCattle.bed
##################################
Bt_Chr3  50155331  50155587
Bt_Chr3 82742456  82742715
Bt_Chr10  22353175  22353449
Bt_Chr21  208228  208481
#############################

less SyntenicRibbons.conf |sort -k 4,5 |awk '{print $4,$5,$6,$1,$2,$3}' |tr " " "\t" |bedtools intersect -wo -a - -b MissingGeneRegionsCattle.bed |less
```
These regions were not found to be in synteny.

### Find regions next to these areas
```
#Create chrom sizes file
less karyotype.conf |awk '{print $4,$6}' |grep "Bt" |tr " " "\t" >Btaurus.ChromSizes.txt
less karyotype.conf |awk '{print $1,$3}' |grep "C" |tr " " "\t" >Ccanadensis.ChromSizes.txt


less SyntenicRibbons.conf |sort -k 4,5 |awk '{if($5>$6){print $4,$6,$5,$1,$2,$3}else {print  $4,$5,$6,$1,$2,$3}}' |tr " " "\t" |bedtools slop -b 10000000 -i - -g Btaurus.ChromSizes.txt |bedtools intersect -wo -a - -b MissingGeneRegionsCattle.bed |less -S
#flanked syntenic coordinate,
Bt_Chr_10       4362722 24905244        Chromosome_06   13394204        13947631        Bt_Chr_10       22353175        22353449        274
Bt_Chr_21       0       16821217        Chromosome_17   62603736        63147732        Bt_Chr_21       208228  208481  253
Bt_Chr_3        47499490        68618306        Chromosome_02   54362766        55451623        Bt_Chr_3        50155331        50155587        256
Bt_Chr_3        60212107        84185631        Chromosome_02   66520471        70230081        Bt_Chr_3        82742456        82742715        259
Bt_Chr_3        67186963        87922758        Chromosome_02   73036876        73736481        Bt_Chr_3        82742456        82742715        259
Bt_Chr_3        77268009        97479276        Chromosome_02   83109250        83348104        Bt_Chr_3        82742456        82742715        259
Bt_Chr_3        81445021        101562854       Chromosome_02   86622338        86740167        Bt_Chr_3        82742456        82742715        259



#find the corresponding regions in elk, though none are syntenic

less SyntenicRibbons.conf |sort -k 4,5 |awk '{if($5>$6){print $4,$6,$5,$1,$2,$3}else {print  $4,$5,$6,$1,$2,$3}}' |tr " " "\t" |bedtools flank -b 10000000 -i - -g Btaurus.ChromSizes.txt |bedtools intersect -wo -a - -b MissingGeneRegionsCattle.bed |cut -f 4,5,6 |bedtools slop -b 10000000 -i - -g Ccanadensis_RenamedSizes.txt |sort -k1,1 -k2,2n |bedtools getfasta -fi ../../../../08_RenameAgain/FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa  -bed - -fo test
new_Assemblathon.pl test
############################################################################################

--------------- Information for assembly 'test' ----------------


                                         Number of scaffolds          7
                                     Total size of scaffolds  139182582
                                            Longest scaffold   23709610
                                           Shortest scaffold   12774400
                                 Number of scaffolds > 1K nt          7 100.0%
                                Number of scaffolds > 10K nt          7 100.0%
                               Number of scaffolds > 100K nt          7 100.0%
                                 Number of scaffolds > 1M nt          7 100.0%
                                Number of scaffolds > 10M nt          7 100.0%
                                          Mean scaffold size   19883226
                                        Median scaffold size   20553427
                                         N50 scaffold length   20553427
                                          L50 scaffold count          4
                                         n90 scaffold length   20117829
                                          L90 scaffold count          6
                                                 scaffold %A      29.82
                                                 scaffold %C      20.17
                                                 scaffold %G      20.16
                                                 scaffold %T      29.78
                                                 scaffold %N       0.06
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

#################################################################################################


0.0006 x 139182582 = ~83,510 N's in the 139Mb section that these genes could reside



#How much of these regions are covered by repeats?
less SyntenicRibbons.conf |sort -k 4,5 |awk '{if($5>$6){print $4,$6,$5,$1,$2,$3}else {print  $4,$5,$6,$1,$2,$3}}' |tr " " "\t" |bedtools flank -b 10000000 -i - -g Btaurus.ChromSizes.txt |bedtools intersect -wo -a - -b MissingGeneRegionsCattle.bed |cut -f 4,5,6 |bedtools slop -b 10000000 -i - -g Ccanadensis_RenamedSizes.txt |sort -k1,1 -k2,2n |bedtools coverage -a - -b ../../../../08_RenameAgain/edtaRepmaskerMerge.gff3 |less

Chromosome_02   44362766        65451623        22599   7947143 21088857        0.3768409
Chromosome_02   56520471        80230081        26142   9158586 23709610        0.3862816
Chromosome_02   63036876        83736481        23754   7946625 20699605        0.3839023
Chromosome_02   73109250        93348104        23719   7705751 20238854        0.3807405
Chromosome_02   76622338        96740167        23642   7568086 20117829        0.3761880
Chromosome_06   3394204 23947631        23371   7414734 20553427        0.3607541
Chromosome_17   52603736        65378136        13703   5695192 12774400        0.4458285


36.1-44.6 % repetitive

```


### Map extracted regions from cattle to elk
```

ln -s ../../08_RenameAgain/FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa
makeblastdb -in FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa  -dbtype nucl -out FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES
blastn -query NucleotidePositionsMissingGenes.bed -db FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES -outfmt 6 -num_threads 6 -out GeneFragments.blastout

Only 4/5 genes had a blast hit to the C. canadensis genome
NC_037331.1:50155331-50155587   Chromosome_03   94.024  251     13      2       1       250     47250736        47250985        2.93e-103       379
NC_037331.1:82742456-82742715   Chromosome_03   93.130  262     13      4       1       260     78406493        78406751        2.96e-103       379
NC_037337.1:22353175-22353449   Chromosome_06   90.114  263     23      1       1       263     20940274        20940533        5.36e-91        339
NC_037348.1:208402-208661       Chromosome_17   87.854  247     20      7       12      257     7755424 7755661 8.60e-74        281

## extract regions for blastx to nr
less GeneFragments.blastout |awk '{print $2,$9,$10}' |tr " " "\t" |bedtools getfasta -fi FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa -bed - -fo NuclRegions.fasta
>Chromosome_03:47250736-47250985
AGGACTGTGGCAGCCAGATCAGGGAGACCCCAAGTTGTTTCTACAAAGCATCTAACAGTATTTTATGACTACAGTAAAAATAGAAAAATGATACTTACCAGGAGGAGTTACTATGATGATAGCTCCTTCTCCAAATATCTTCTTCCAGCCTGAGCTGTCACGATGGTAAAAGCTCCTACAATAATCCAGCTTCTTCTACCTTTATAAAGTGACTAGAAGTCAAAAAACTTCATACTGGCAGTTTGATCT
>Chromosome_03:78406493-78406751
tatacatatatatatatagcaataacacatatataaacatataCAATGCCAAATAATGGGCACAAATAATTGACTAATACTTTTTTTCATTTCTTGTAGTTGTCAGCTCAATTATCCCTACCACTGAGTCTCCCAGTGACTGTTTAAACCATGAAAGTGGTAAGTGTTTGTATATGATTGCCTTTATGTCATGCTTCAGTATTTTATTTCCTTTCTAATCTATTTATCTTTTGCACTTCTCTGGCTAAAAGTACAAAA
>Chromosome_06:20940274-20940533
TGTTTTCTTTTGAACTAGAATGATTTTGGGCCCCGTCTAACAATGACAGAGTTTGCTTTAACCTCAATTCTGAAGTTGGCCCTTTGGCTTCAACTAACCTTGGGGCTCATAGCAGGTCTCTGAAACTGAAGTTTGTTTGCTGTTTTGAGGTTCCATCGGTTTTGGAGTTCCTGGAACAAACAGCAGTGGTCAGGGTTACAATGACAGCTGAGAGGCTCAAGGGCTTTCAGGTTTTAGAATTGGAGGTGACTTTAAGTAA
>Chromosome_17:7755424-7755661
AAAAACACACCCTGCCTCCTTTCGGTCCTGACTCCCATCCCAAGTTCACCCACATCCCGCTGTCCCCGGACTCTGGACACCAGATTCACCTGAGGAGACGGTGACCAGGGCACCCTGGCCCAAGATGGAAGTCAGCATAGTCACAGTGGGCTGGCCCTGTACTGGGGGGCATGGAAACCCACCCAGCATCCCCAGGGAGCCCCCAGGAGCTGCAAGTCTGTCTTTCTGAAGACACTT

#blastx to nr
>Chromosome_03:47250736-47250985
T-cell receptor gamma chain V-J-C region (clone SFTG4) - sheep (fragment) [Ovis aries] 	Ovis aries 	52.0 	52.0 	33% 	2e-05 	82.14% 	246 	A43546
>Chromosome_03:78406493-78406751
No hit
>Chromosome_06:20940274-20940533
hypothetical protein Celaphus_00006467 [Cervus elaphus hippelaphus] 	Cervus elaphus hippelaphus 	53.9 	53.9 	30% 	2e-06 	92.31% 	168 	OWK09794.1
>Chromosome_17:7755424-7755661
hypothetical protein E2562_002192 [Oryza meyeriana var. granulata] 	Oryza meyeriana var. granulata 	43.1 	43.1 	65% 	0.025 	38.46% 	333 	KAF0922982.1
```

Determine if these regions have any expression
```
less GeneFragments.blastout |awk '{print $2,"NA","gene",$9,$10,"+",".",".","ID="$1"gene" }' |less
Chromosome_03 NA gene 47250736 47250985 + . . ID=NC_037331.1:50155331-50155587gene
Chromosome_03 NA gene 78406493 78406751 + . . ID=NC_037331.1:82742456-82742715gene
Chromosome_06 NA gene 20940274 20940533 + . . ID=NC_037337.1:22353175-22353449gene
Chromosome_17 NA gene 7755424 7755661 + . . ID=NC_037348.1:208402-208661gene


```



### Add 100bp to each gene 5'
```
less GeneFragments.blastout |awk '{print $2,$9-100,$10}' |tr " " "\t" |bedtools getfasta -fi FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa -bed - -fo NuclRegions.fasta
##################################
>Chromosome_03:47250636-47250985
ACCGCTTACAAGCAGGGACAAGAATAAGCCATTTACAAAGTGGTTTAGTTTCCATATAGAGTGACCAGTAGCTATAAAACTAAGTTTAGACGCTGAGACCAGGACTGTGGCAGCCAGATCAGGGAGACCCCAAGTTGTTTCTACAAAGCATCTAACAGTATTTTATGACTACAGTAAAAATAGAAAAATGATACTTACCAGGAGGAGTTACTATGATGATAGCTCCTTCTCCAAATATCTTCTTCCAGCCTGAGCTGTCACGATGGTAAAAGCTCCTACAATAATCCAGCTTCTTCTACCTTTATAAAGTGACTAGAAGTCAAAAAACTTCATACTGGCAGTTTGATCT
>Chromosome_03:78406393-78406751
CAAATGAAATCTATGAATTACATGGAATCAGACCAACTAGTTTTGGAACTGTAGGGCTCATATTTTTAATAACTGTCTCTCTCAGCTGATTatctatctatatacatatatatatatagcaataacacatatataaacatataCAATGCCAAATAATGGGCACAAATAATTGACTAATACTTTTTTTCATTTCTTGTAGTTGTCAGCTCAATTATCCCTACCACTGAGTCTCCCAGTGACTGTTTAAACCATGAAAGTGGTAAGTGTTTGTATATGATTGCCTTTATGTCATGCTTCAGTATTTTATTTCCTTTCTAATCTATTTATCTTTTGCACTTCTCTGGCTAAAAGTACAAAA
>Chromosome_06:20940174-20940533
TCACAGTCAGCAAAGCAAGCTTTGGAAAAACCCATTCCTATAGTCAGGAAGACTCTGCTTGAAGAGTTAGAGGGAGCCAACCGGACTCAAGTATGGTTTTTGTTTTCTTTTGAACTAGAATGATTTTGGGCCCCGTCTAACAATGACAGAGTTTGCTTTAACCTCAATTCTGAAGTTGGCCCTTTGGCTTCAACTAACCTTGGGGCTCATAGCAGGTCTCTGAAACTGAAGTTTGTTTGCTGTTTTGAGGTTCCATCGGTTTTGGAGTTCCTGGAACAAACAGCAGTGGTCAGGGTTACAATGACAGCTGAGAGGCTCAAGGGCTTTCAGGTTTTAGAATTGGAGGTGACTTTAAGTAA
>Chromosome_17:7755324-7755661
GCAAGGGAGACGCCCCCTGAGGACACGGTGACCGGGGTGTGCTGGCCCCAGAGCTCCAAGTCCCAGCAACACAGCTTCTCCTCTGCTGCTGCTTATATGCAAAAACACACCCTGCCTCCTTTCGGTCCTGACTCCCATCCCAAGTTCACCCACATCCCGCTGTCCCCGGACTCTGGACACCAGATTCACCTGAGGAGACGGTGACCAGGGCACCCTGGCCCAAGATGGAAGTCAGCATAGTCACAGTGGGCTGGCCCTGTACTGGGGGGCATGGAAACCCACCCAGCATCCCCAGGGAGCCCCCAGGAGCTGCAAGTCTGTCTTTCTGAAGACACTT
```



Blast output from above to NR with 100bp added to 5' end, making sure I didnt cut off the start codon
```
>Chromosome_03:47250636..47250985
T-cell receptor gamma chain V-J-C region (clone SFTG4) - sheep (fragment) [Ovis aries] 	Ovis aries 	52.0 	52.0 	24% 	5e-05 	82.14% 	246 	A43546
>Chromosome_03:78406393..78406751
No Hit
>Chromosome_06:20940174..20940533
hypothetical protein Celaphus_00006467 [Cervus elaphus hippelaphus] 	Cervus elaphus hippelaphus 	53.9 	53.9 	21% 	5e-06 	92.31% 	168 	OWK09794.1
>Chromosome_17:7755324-7755661
No Hit
```
24 and 21% coverage is low.


### Are these regions expressed?
```
 for f in *Frag*txt; do echo "<( cut -f 7 "$f"|grep -v \"#\" ) "; done |tr "\n" " " |less
paste <( cut -f 1,7 Elk-kidney_S25_L003_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 Elk-kidney_S25_L004_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 Elk-lung_S26_L003_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 Elk-lung_S26_L004_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 Elk-Mes-LN_S24_L003_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 Elk-Mes-LN_S24_L004_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 Elk-muscle_S21_L003_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 Elk-muscle_S21_L004_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 ElkpscapLN_S22_L003_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 ElkpscapLN_S22_L004_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 Elk-spleen_S23_L003_R1_001Fragments_counts_genes.txt|grep -v "#" )  <( cut -f 7 Elk-spleen_S23_L004_R1_001Fragments_counts_genes.txt|grep -v "#" ) |less
#################################################################################################################################################
Geneid  Elk-kidney_S25_L003_R1_001_sorted.bam   Elk-kidney_S25_L004_R1_001_sorted.bam   Elk-lung_S26_L003_R1_001_sorted.bam     Elk-lung_S26_L004_R1_001_sorted.bam     Elk-Mes-LN_S24_L003_R1_001_sorted.bam   Elk-Mes-LN_S24_L004_R1_001_sorted.bam   Elk-muscle_S21_L003_R1_001_sorted.bam   Elk-muscle_S21_L004_R1_001_sorted.bam   ElkpscapLN_S22_L003_R1_001_sorted.bam   ElkpscapLN_S22_L004_R1_001_sorted.bam   Elk-spleen_S23_L003_R1_001_sorted.bam   Elk-spleen_S23_L004_R1_001_sorted.bam
NC_037331.1:50155331-50155587gene       0       0       0       0       15      4       0       0       15      23      35      27
NC_037331.1:82742456-82742715gene       0       0       0       0       0       2       0       0       2       1       1       0
NC_037337.1:22353175-22353449gene       8       9       84      90      147     129     1       1       303     286     203     222
NC_037348.1:208402-208661gene   0       0       0       1       29      56      0       0       15      10      7       14
#################################################################################################################################################

```
