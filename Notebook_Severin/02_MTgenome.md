# Finding the Mitochondrial genome

* 2019 June 13
* /project/elk_bison_genomics/Severin/02_MTgenome

Always good to know which scaffold is the mitochondrial genome or to determine if it is missing or present in multiple copies.

I searched NCBI for  ```"(mitochondrial "Cervus canadensis"[porgn:__txid1574408])"```

  * https://www.ncbi.nlm.nih.gov/nuccore

I then downloaded the the fasta sequences and placed them in a file named Cercan_mitogenes.fasta


## Blast mitogenes against the Elk genome assembly

```
blastn  -db ../01_phiXcontamination/elk_blastDB -query Cercan_mitogenes.fasta -num_threads 20 -outfmt 6 | sort -k 7n > mitogenes_2_Elk.blastnout
```

There are hits to 28 scaffolds.

```
more mitogenes_2_Elk.blastnout | awk '{print $1,$2}' | sort | uniq | awk '{print $2}' | sort | uniq -c | wc
     28      56     918
```

11 scaffolds have at least 10 hits from different mitochondrial IDs

```
more mitogenes_2_Elk.blastnout | awk '{print $1,$2}' | sort | uniq | awk '{print $2}' | sort | uniq -c | awk '$1>10'
    119 ScoY8t2_1063;HRSCAF=1174
     68 ScoY8t2_1222;HRSCAF=1339
     64 ScoY8t2_203;HRSCAF=236
     29 ScoY8t2_23257;HRSCAF=23437
     71 ScoY8t2_23276;HRSCAF=23464
    114 ScoY8t2_23283;HRSCAF=23527
     15 ScoY8t2_23300;HRSCAF=23700
    123 ScoY8t2_238;HRSCAF=281
     35 ScoY8t2_282;HRSCAF=329
     25 ScoY8t2_288;HRSCAF=336
     44 ScoY8t2_401;HRSCAF=471

more mitogenes_2_Elk.blastnout | awk '{print $1,$2}' | sort | uniq | awk '{print $2}' | sort | uniq -c | awk '$1>10' | wc
     11      22     361
```

Assumption: Expected mitochondrial genome size is around 16kb?  

* [Roe Deer](https://www.tandfonline.com/doi/full/10.1080/23802359.2017.1365645)
* [white tail deer](https://academic.oup.com/jmammal/article/97/1/234/2459714)


## Blast a related full mitochondrial genome to see if we can get a clearer picture of which scaffolds may be Mitochondrial.

#### Download white tail deer MT genome and unzip
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Odocoileus_virginianus_texanus/CHR_MT/9880_ref_Ovir.te_1.0_chrMT.fa.gz
gunzip 9880_ref_Ovir.te_1.0_chrMT.fa.gz
```

## Blast mitogenes against the Elk genome assembly

```
blastn  -db ../01_phiXcontamination/elk_blastDB -query 9880_ref_Ovir.te_1.0_chrMT.fa -num_threads 20 -outfmt 6 | sort -k 7n > OvirMT_2_Elk.blastnout
```

#### Best candidates for the MT genome

It looks like there are several scaffolds that could be the MT genome and it appears that it may have rolling assembled so we will need to manually fix that.  The best candidates so far are **ScoY8t2_238;HRSCAF=281** and **ScoY8t2_401;HRSCAF=471**

```
scaffold                length
ScoY8t2_401;HRSCAF=471 48531  
ScoY8t2_238;HRSCAF=281 42018

```

```
more OvirMT_2_Elk.blastnout  | more| sort -k 2rn | sort -k 2,2 -k 7,7n | awk '$9<1000000'
ref|NC_015247.1|        ScoY8t2_184;HRSCAF=213  79.326  445     68      15      28      456     60921   61357   2.85e-75        292
ref|NC_015247.1|        ScoY8t2_184;HRSCAF=213  74.391  2054    419     60      5288    7317    46646   48616   0.0     784
ref|NC_015247.1|        ScoY8t2_238;HRSCAF=281  87.535  11280   1366    32      1       11266   28100   39353   0.0     13051
ref|NC_015247.1|        ScoY8t2_238;HRSCAF=281  87.591  11282   1356    35      1       11266   11670   22923   0.0     13084
ref|NC_015247.1|        ScoY8t2_238;HRSCAF=281  90.501  5169    462     24      11329   16477   22940   28099   0.0     6802
ref|NC_015247.1|        ScoY8t2_238;HRSCAF=281  91.943  2656    186     23      11329   13965   39370   42016   0.0     3694
ref|NC_015247.1|        ScoY8t2_238;HRSCAF=281  87.309  654     83      0       15824   16477   11016   11669   0.0     749
ref|NC_015247.1|        ScoY8t2_281;HRSCAF=328  74.540  2066    419     65      5288    7326    58068   56083   0.0     806
ref|NC_015247.1|        ScoY8t2_298;HRSCAF=349  74.008  2066    420     71      5289    7326    27351   25375   0.0     736
ref|NC_015247.1|        ScoY8t2_300;HRSCAF=351  74.554  2075    407     69      5288    7326    53944   51955   0.0     798
ref|NC_015247.1|        ScoY8t2_384;HRSCAF=451  74.867  2065    406     69      5288    7326    4600    2623    0.0     837
ref|NC_015247.1|        ScoY8t2_400;HRSCAF=470  74.201  2066    413     70      5289    7326    103381  105354  0.0     756
ref|NC_015247.1|        ScoY8t2_401;HRSCAF=471  86.952  7388    938     20      1       7379    20299   12929   0.0     8303
ref|NC_015247.1|        ScoY8t2_401;HRSCAF=471  95.053  283     14      0       11003   11285   25714   25432   2.07e-121       446
ref|NC_015247.1|        ScoY8t2_401;HRSCAF=471  89.762  5167    480     26      11329   16477   25435   20300   0.0     6571
```

There is also something about the gene located between bases 5288 and 7326 in the white tail deer MT genome as that 2000 base sequence is found in several small scaffolds.  
