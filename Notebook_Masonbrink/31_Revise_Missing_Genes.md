# More detail is needed on these missing immunoglobulin genes in elk. Restart


### Identify which genes are not found in Elk
```
/work/gif/remkv6/Olsen/Elk/09_MapRnaseq2Cattle

This grabs only those protein names in Cattle (missing from elk), and cuts the names after the asterisk.  
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
```

### What are the nucleotide fasta sequences to the hits in cattle that are "missing" in Elk, so I can map them to elk?
```
cat ElkBBH.out CowBBH.out CowBBH.out |cut -f 1 |sort |uniq -c |sort -rn |awk '$1==2 {print $2}'  |sed 's/\*/\\\*/1' |awk '{print $1}' |grep -f - CowBBH.out |awk '{if($10>$9) {print $2,$9,$10,$1}else {print $2,$10,$9,$1}}' |sort -k1,1 -k2,3n |tr " " "\t" |uniq|sed 's/\*/\t/g' |awk '{print $1,$2,$3}' |tr " " "\t" |bedtools merge -i - |bedtools getfasta -bed - -fi GCF_002263795.1_ARS-UCD1.2_genomic.fna -fo MissingGenesInElkCattleNucSequence.fasta

vi queryDirect.fasta
#####################################################
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
################################################################################
How many genees
grep -c ">" queryDirect.fasta
34
```

### Blast cattle nucleotide sequences to Elk genome
```
ml gcc; ml blast-plus/2.11.0-py3-4pqzweg
blastn -query queryDirect.fasta -db FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES -outfmt 6 -num_threads 6 -out queryDirectGeneFragments.blastout

#Blast output
####################################################################################################################################################
NC_037331.1:82699801-82699857   Chromosome_03   95.833  48      2       0       3       50      78366825        78366872        1.57e-13        78.7
NC_037331.1:82707936-82707992   Chromosome_03   96.429  56      2       0       1       56      78375135        78375190        5.59e-18        93.5
NC_037331.1:82742556-82742615   Chromosome_03   98.182  55      1       0       1       55      78406593        78406647        4.76e-19        97.1
NC_037337.1:22272459-22272521   Chromosome_06   93.548  62      4       0       1       62      20860857        20860918        6.71e-18        93.5
NC_037337.1:22277164-22277220   Chromosome_06   92.727  55      4       0       1       55      20865405        20865459        4.35e-14        80.5
NC_037337.1:22280478-22280534   Chromosome_06   98.148  54      1       0       1       54      20868786        20868839        1.55e-18        95.3
NC_037337.1:22280478-22280534   Chromosome_06   98.148  54      1       0       1       54      20871340        20871393        1.55e-18        95.3
NC_037337.1:22292039-22292095   Chromosome_06   94.643  56      3       0       1       56      20879847        20879902        2.60e-16        87.9
NC_037337.1:22297656-22297712   Chromosome_06   98.182  55      0       1       1       55      20885452        20885505        1.55e-18        95.3
NC_037337.1:22301761-22301817   Chromosome_06   96.429  56      2       0       1       56      20889617        20889672        5.59e-18        93.5
NC_037337.1:22304918-22304974   Chromosome_06   94.643  56      3       0       1       56      20892783        20892838        2.60e-16        87.9
NC_037337.1:22308157-22308210   Chromosome_06   97.674  43      1       0       1       43      20896025        20896067        1.82e-12        75.0
NC_037337.1:22309187-22309240   Chromosome_06   90.196  51      5       0       1       51      20897062        20897112        3.05e-10        67.6
NC_037337.1:22312866-22312925   Chromosome_06   94.828  58      3       0       1       58      20900746        20900803        2.21e-17        91.6
NC_037337.1:22323995-22324054   Chromosome_06   94.915  59      3       0       1       59      20911916        20911974        6.15e-18        93.5
NC_037337.1:22335369-22335428   Chromosome_06   98.182  55      1       0       5       59      20923357        20923411        4.76e-19        97.1
NC_037348.1:107140-107256       Chromosome_17   87.156  109     12      2       1       108     7781035 7781142 2.12e-26        122
####################################################################################################################################################
```
### Coordinates of the cattle gene nucleotide sequences found in the elk genome, not found by tblastn
```
less queryDirectGeneFragments.blastout |awk '{if($9>$10){print $2,$10,$9} else {print $2,$9,$10}} ' |less
###################################################################################################################
Chromosome_03 78366825 78366872
Chromosome_03 78375135 78375190
Chromosome_03 78406593 78406647
Chromosome_06 20860857 20860918
Chromosome_06 20865405 20865459
Chromosome_06 20868786 20868839
Chromosome_06 20871340 20871393
Chromosome_06 20879847 20879902
Chromosome_06 20885452 20885505
Chromosome_06 20889617 20889672
Chromosome_06 20892783 20892838
Chromosome_06 20896025 20896067
Chromosome_06 20897062 20897112
Chromosome_06 20900746 20900803
Chromosome_06 20911916 20911974
Chromosome_06 20923357 20923411
Chromosome_17 7781035 7781142
###########################################################################################################################
Chromosome 3, 10, and 21 of Bt, so these are not involved in fusion/fissions
```


#### Genes not in Elk, but in Cattle
```
Those that were not found by blasting nucleotide sequences.  
less queryDirectGeneFragments.blastout |awk '{print $1}' |grep -v -f - <(grep ">" queryDirect.fasta) |less
########################################
>NC_037337.1:22273005-22273061
>NC_037337.1:22278352-22278411
>NC_037337.1:22283266-22283322
>NC_037337.1:22284736-22284792
>NC_037337.1:22285297-22285353
>NC_037337.1:22290540-22290599
>NC_037337.1:22300122-22300178
>NC_037337.1:22307467-22307520
>NC_037337.1:22315617-22315670
>NC_037337.1:22317784-22317846
>NC_037337.1:22321231-22321290
>NC_037337.1:22324885-22324938
>NC_037337.1:22334740-22334796
>NC_037337.1:22353275-22353349
>NC_037348.1:59222-59287
>NC_037348.1:208328-208381
>NC_037348.1:208502-208561
########################################

#Add annotations. Manually investigated for gene names and to associate with the IMGT database proteins, some of which are almost identical.
ProtNames.list
#######################################
TRAJ3
TRAJ6
TRAJ8-2
TRAJ11
TRAJ12
TRAJ17
TRAJ27
TRAJ33
TRAJ40
TRAJ42
TRAJ46
TRAJ49
TRAJ56
TRDC
IGHD
IGHJ1-2
IGHJ2-2
IGHJ1-2
IGHJ2-1
IGHJ1-1
#######################################
```

### Are there just certain isoforms missing for each gene type.  
```
Verify which are missing and which have other variants present.  These are present in the elk genome
tr "," "\n" <ProtNames.list |sed 's/>//g' | while read line; do grep "$line" ElkBBH.out ;done |less#######################################################################
AY227782_TRAJ8-2*02_Bos_taurus_F_J-REGION_38611..38670_60_nt_3_______19_AA_19+0=19_____ Chromosome_06   89.474  19      2       0       1       19      20868841        20868785        5.27e-04        38.5    19
      96780817        0/-1    Chromosome_06   100
AY227782_TRAJ11*02_Bos_taurus_F_J-REGION_37141..37200_60_nt_3_______19_AA_19+0=19_____  Chromosome_06   89.474  19      2       0       1       19      20872867        20872811        4.51e-04        38.5    19
            96780817        0/-1    Chromosome_06   100      
IMGT000049_TRDC*01_Bos_taurus_Hereford_F_EX1_n,3230147..3230424_279_nt_1_+1_-1__93_AA_93+0=93_____      Chromosome_06   91.209  91      8       0       2       92      20941304        20941032        3.82e-50
                    174     93      96780817        0/-3    Chromosome_06   98
IMGT000049_TRDC*01_Bos_taurus_Hereford_F_EX3_g,3231497..3231609_114_nt_1_+1_-1___38_AA_38+0=38_____     Chromosome_06   100.000 38      0       0       1       38      20940006        20939893        1.46e-15
                    72.4    38      96780817        0/-2    Chromosome_06   100
IMGT000049_TRDC*01_Bos_taurus_Hereford_F_EX4UTR_g,3233021..3233919_900_nt_1_+1_____300_AA_300+0=300_____        Chromosome_06   74.026  154     27      3       45      190     20938355        20937909        8.69e-95        194     300     96780817        0/-3    Chromosome_06   96
KT723008_IGHD*02_Bos_taurus_Holstein_ORF_CH1_n,525572..525894_324_nt_1_+1_-1___108_AA_108+0=108_____    Chromosome_17   89.130  92      9       1       3       93      7747488 7747213 1.49e-37        138     108     65378136        0/-1    Chromosome_17   84
KT723008_IGHD*02_Bos_taurus_Holstein_ORF_CHS_g,531422..531606_186_nt_1_+1_____62_AA_62+0=62_____        Chromosome_17   77.419  62      14      0       1       62      7731939 7731754 2.60e-22        92.8    62
                  65378136        0/-1    Chromosome_17   100
KT723008_IGHD*02_Bos_taurus_Holstein_ORF_CH2_g,527840..528162_324_nt_1_+1_-1___108_AA_108+0=108_____    Chromosome_17   83.333  108     18      0       1       108     7735874 7735551 7.91e-44        156     108     65378136        0/-2    Chromosome_17   100
KT723008_IGHD*02_Bos_taurus_Holstein_ORF_CH3_g,528268..528587_321_nt_1_+1_-1___107_AA_107+0=107_____    Chromosome_17   85.047  107     16      0       1       107     7735429 7735109 1.19e-52        181     107     65378136        0/-3    Chromosome_17   100
KT723008_IGHD*02_Bos_taurus_Holstein_ORF_M1_g,532785..532927_144_nt_1_+1_____48_AA_48+0=48_____ Chromosome_17   91.667  48      4       0       1       48      7730800 7730657 1.43e-10        58.9    48      65378136        0/-3    Chromosome_17   100

#######################################################################
Yes,  there were 4 of the above 20 genes that had another isoform variant present in the elk genome.
```

### Compile Annotations for these "missing" genes lacking a nucleotide and peptide locus in the elk genome.
```
>NC_037337.1:22273005-22273061 IMGT000049_TRAJ3 -- No other copies present
>NC_037337.1:22278352-22278411 IMGT000049_TRAJ6 -- No other copies present
>NC_037337.1:22283266-22283322 IMGT000049_TRAJ8-2 -- AY227782_TRAJ8-2*02_Bos_taurus_F_J-REGION_38611..38670_60_nt_3_______19_AA_19+0=19_____
>NC_037337.1:22284736-22284792 IMGT000049_TRAJ11 -- AY227782_TRAJ11*02_Bos_taurus_F_J-REGION_37141..37200_60_nt_3_______19
>NC_037337.1:22285297-22285353 IMGT000049_TRAJ12 -- No other copies present
>NC_037337.1:22290540-22290599 IMGT000049_TRAJ17 -- No other copies present
>NC_037337.1:22300122-22300178 IMGT000049_TRAJ27 -- No other cppies present
>NC_037337.1:22307467-22307520 IMGT000049_TRAJ33 -- No other copies present
>NC_037337.1:22315617-22315670 IMGT000049_TRAJ40 -- No other copies present
>NC_037337.1:22317784-22317846 IMGT000049_TRAJ42 -- No other copies present
>NC_037337.1:22321231-22321290 IMGT000049_TRAJ46 -- No other copies present
>NC_037337.1:22324885-22324938 IMGT000049_TRAJ49 -- No other copies present
>NC_037337.1:22334740-22334796 IMGT000049_TRAJ56 -- No other copies present
>NC_037337.1:22353275-22353349 MGT000049_TRDC -- multiple other isoforms present
>NC_037348.1:59222-59287 KT723008_IGHD -- multiple other isoforms present
>NC_037348.1:208328-208381 AY149283_IGHJ1-2,KT723008_IGHJ2-2,NW_001494075_IGHJ1-2
>NC_037348.1:208502-208561 KT723008_IGHJ2-1,NW_001494075_IGHJ1-1 -- No other copies present
```

### Was the coding portion of this gene deleted and surrounding sequences persist? Extract 20bp borders
```
#genes in the elk genome via blastn

#genes not in the elk genome via direct nucleotide sequence, extracting cattle sequences plus 20bp on each side.
less queryDirectGeneFragments.blastout |awk '{print $1}' |grep -v -f - <(grep ">" queryDirect.fasta) |sed 's/>//g' |sed 's/:/\t/g' |sed 's/-/\t/g' |awk '{print $1":"$2-20"-"$3+20}' |while read line ;do samtools faidx  GCF_002263795.1_ARS-UCD1.2_genomic.fna  $line; done >Cattle20bpGeneExtractsMissing.fasta

#see if I can BLAST those genes using their peripheral sequences to get a hit
blastn -query Cattle20bpGeneExtractsMissing.fasta  -db FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES -outfmt 6 -num_threads 4 -out Cattle20bpGeneExtractsMissing.blastout
################################################################################
NC_037331.1:50155411-50155507   Chromosome_03   96.907  97      3       0       1       97      47250816        47250912        9.81e-39        163
NC_037337.1:22278332-22278431   Chromosome_06   92.000  100     6       1       1       100     20866573        20866670        1.72e-31        139
NC_037337.1:22283246-22283342   Chromosome_06   94.792  96      5       0       1       96      20871319        20871414        7.63e-35        150
NC_037337.1:22283246-22283342   Chromosome_06   93.750  96      6       0       1       96      20868765        20868860        3.55e-33        145
NC_037337.1:22284716-22284812   Chromosome_06   94.253  87      5       0       8       94      20872798        20872884        7.69e-30        134
NC_037337.1:22285277-22285373   Chromosome_06   98.592  71      1       0       1       71      20873358        20873428        1.29e-27        126
NC_037337.1:22300102-22300198   Chromosome_06   95.876  97      4       0       1       97      20887918        20888014        4.56e-37        158
NC_037337.1:22307447-22307540   Chromosome_06   97.872  94      2       0       1       94      20895321        20895414        9.39e-39        163
NC_037337.1:22315597-22315690   Chromosome_06   94.681  94      5       0       1       94      20903523        20903616        9.45e-34        147
NC_037337.1:22321211-22321310   Chromosome_06   96.000  100     4       0       1       100     20909128        20909227        1.02e-38        163
NC_037337.1:22353255-22353369   Chromosome_06   83.478  115     16      1       1       115     20940354        20940465        7.58e-21        104
NC_037348.1:208482-208581       Chromosome_17   92.391  92      5       2       7       98      7755492 7755581 1.04e-28        130
################################################################################

#How many genes have direct deletions causing problems with tblastn and blastn direct
less Cattle20bpGeneExtractsMissing.blastout |awk '{print $1}' |sort|uniq|wc
     11      11     326
#Bed file of "missing" genes with direct deletions    
less Cattle20bpGeneExtractsMissing.blastout|awk '{if($9>$10){print $2,$10,$9} else {print $2,$9,$10}} ' |less
##################################
Chromosome_03 47250816 47250912
Chromosome_06 20866573 20866670
Chromosome_06 20871319 20871414
Chromosome_06 20868765 20868860
Chromosome_06 20872798 20872884
Chromosome_06 20873358 20873428
Chromosome_06 20887918 20888014
Chromosome_06 20895321 20895414
Chromosome_06 20903523 20903616
Chromosome_06 20909128 20909227
Chromosome_06 20940354 20940465
Chromosome_17 7755492 7755581
##################################



#How many are still not found
less Cattle20bpGeneExtractsMissing.blastout |awk '{print $1}' |sort|uniq|grep -v -w -f - <(grep ">" Cattle20bpGeneExtractsMissing.fasta ) |wc
      7       7     207
Not found
##############################
>NC_037337.1:22272985-22273081
>NC_037337.1:22290520-22290619
>NC_037337.1:22317764-22317866
>NC_037337.1:22324865-22324958
>NC_037337.1:22334720-22334816
>NC_037348.1:59202-59307
>NC_037348.1:208308-208401
###############################

#Grab those 7 fastas with another 800bp on each side for mapping rnaseq
less Cattle20bpGeneExtractsMissing.blastout|awk '{print $1}' |grep -v -f - <(grep ">"  Cattle20bpGeneExtractsMissing.fasta) |sed 's/>//g' |sed 's/:/\t/g' |sed 's/-/\t/g' |awk '
{print $1":"$2-80"-"$3+80}' |while read line ;do samtools faidx  GCF_002263795.1_ARS-UCD1.2_genomic.fna  $line; done >Cattle100bpGeneExtractsMissing.fasta
```

# Create a cattle+elk Fasta and gff to measure read counts
```

Missing.bed
########################################################
#bed file of those found with direct blast
Chromosome_03 78366825 78366872
Chromosome_03 78375135 78375190
Chromosome_03 78406593 78406647
Chromosome_06 20860857 20860918
Chromosome_06 20865405 20865459
Chromosome_06 20868786 20868839
Chromosome_06 20871340 20871393
Chromosome_06 20879847 20879902
Chromosome_06 20885452 20885505
Chromosome_06 20889617 20889672
Chromosome_06 20892783 20892838
Chromosome_06 20896025 20896067
Chromosome_06 20897062 20897112
Chromosome_06 20900746 20900803
Chromosome_06 20911916 20911974
Chromosome_06 20923357 20923411
Chromosome_17 7781035 7781142
#bed file of those found with blast +borders
Chromosome_03 47250816 47250912
Chromosome_06 20866573 20866670
Chromosome_06 20871319 20871414
Chromosome_06 20868765 20868860
Chromosome_06 20872798 20872884
Chromosome_06 20873358 20873428
Chromosome_06 20887918 20888014
Chromosome_06 20895321 20895414
Chromosome_06 20903523 20903616
Chromosome_06 20909128 20909227
Chromosome_06 20940354 20940465
Chromosome_17 7755492 7755581
#those missing completely
NC_037337.1:22272905-22273161   1       257
NC_037337.1:22290440-22290699   1       260
NC_037337.1:22317684-22317946   1       263
NC_037337.1:22324785-22325038   1       254
NC_037337.1:22334640-22334896   1       257
NC_037348.1:59122-59387 1       266
NC_037348.1:208228-208481       1       254
##########################################################



#Fasta sequences of those that could not be found with tblastn, blastn, or blastn +20bp borders. Extracted with an additional 80bp+20bp+seq+2-bp+80bp so these regions might map whole reads to compete with the genome

>NC_037337.1:22272985-22273081
>NC_037337.1:22290520-22290619
>NC_037337.1:22317764-22317866
>NC_037337.1:22324865-22324958
>NC_037337.1:22334720-22334816
>NC_037348.1:59202-59307
>NC_037348.1:208308-208401

#gets the sequence lengths of the added cattle regions
bioawk -c fastx '{print $name,"1",length($seq)}' Cattle100bpGeneExtractsMissing.fasta |
cat Cattle100bpGeneExtractsMissing.fasta FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa >AddCattleFinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa

#putative gff for expression
less Missing.bed |tr "\t" " " |awk '{print $1,"NA","gene",$2,$3,"-",".",".","ID="$1":"$2"-"$3}' |tr " " "\t" >AddCattleFinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.gff

```

#transfer to ceres for expression
```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/32_MissingGeneExpression
sed -i 's/MissingGenesAddedToElkGenome.fasta/AddCattleFinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa/g' AlignAndCount*sub
sed -i 's/Fragments.gff/AddCattleFinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.gff/g' AlignAndCount*sub

ml hisat2;hisat2-build AddCattleFinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa

```
