# Phi X contamination

Contamination from PhiX is a real problem. All assemblies should be checked.

## Download PhiX sequence from NCBI

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna

seqlen.awk NC_001422.fna
gi|9626372|ref|NC_001422.1| 5386

```

## Make a Blast database for the Elk genome

```
module load blast+/2.6.0
makeblastdb -in north_american_elk_15Jun2018_oY8t2.fasta -input_type fasta -dbtype nucl -parse_seqids -out elk_blastDB
makeblastdb -in north_american_elk_15Jun2018_oY8t2.fasta -input_type fasta -dbtype prot -parse_seqids -out elk_blastDB

#Adding sequences from FASTA; added 23302 sequences in 60.8183 seconds.
```

## Perform blast to find potential contamination in the genome

```
blastn  -db elk_blastDB -query NC_001422.fna -outfmt 6 | sort -k 7n > phiX_2_Elk.blastnout

tblastx -db elk_blastDB -query NC_001422.fna -outfmt 6 | sort -k 7n > phiX_2_Elk.tblastxout
```

#### Header information for blast output
```
 1.	 qseqid	 query (e.g., gene) sequence id
 2.	 sseqid	 subject (e.g., reference genome) sequence id
 3.	 pident	 percentage of identical matches
 4.	 length	 alignment length
 5.	 mismatch	 number of mismatches
 6.	 gapopen	 number of gap openings
 7.	 qstart	 start of alignment in query
 8.	 qend	 end of alignment in query
 9.	 sstart	 start of alignment in subject
 10.	 send	 end of alignment in subject
 11.	 evalue	 expect value
 12.	 bitscore	 bit score
 ```


 # No significant hits were identified to indicate contamination with PhiX

 ```
 gi|9626372|ref|NC_001422.1|     ScoY8t2_23300;HRSCAF=23700      40.476  42      25      0       833     958     10871943        10872068        9.5     38.2
 gi|9626372|ref|NC_001422.1|     ScoY8t2_288;HRSCAF=336  30.952  42      29      0       1605    1730    68865940        68866065        9.8     28.1
 gi|9626372|ref|NC_001422.1|     ScoY8t2_288;HRSCAF=336  46.667  15      8       0       1843    1887    68866162        68866206        9.8     25.4
 gi|9626372|ref|NC_001422.1|     ScoY8t2_23282;HRSCAF=23520      39.130  46      28      0       2362    2499    36233823        36233686        3.7     39.6
 gi|9626372|ref|NC_001422.1|     ScoY8t2_5897;HRSCAF=6038        34.483  29      19      0       2686    2772    37840227        37840141        9.8     30.9
 gi|9626372|ref|NC_001422.1|     ScoY8t2_5897;HRSCAF=6038        39.286  28      17      0       2788    2871    37840098        37840015        9.8     22.6
 gi|9626372|ref|NC_001422.1|     ScoY8t2_117;HRSCAF=135  30.357  56      39      0       3086    2919    21508566        21508399        9.5     38.2
 gi|9626372|ref|NC_001422.1|     ScoY8t2_288;HRSCAF=336  32.558  43      29      0       3193    3065    80424996        80425124        0.31    33.6
 gi|9626372|ref|NC_001422.1|     ScoY8t2_23240;HRSCAF=23404      30.645  62      43      0       3253    3068    13822539        13822354        9.5     38.2
 gi|9626372|ref|NC_001422.1|     ScoY8t2_288;HRSCAF=336  50.000  36      18      0       3265    3158    23063935        23063828        2.4     38.6
 gi|9626372|ref|NC_001422.1|     ScoY8t2_288;HRSCAF=336  34.286  35      23      0       3271    3167    80424882        80424986        0.31    25.8
 gi|9626372|ref|NC_001422.1|     ScoY8t2_117;HRSCAF=135  30.882  68      47      0       3381    3584    7671150 7670947 2.7     40.0
 gi|9626372|ref|NC_001422.1|     ScoY8t2_23291;HRSCAF=23585      41.667  48      28      0       3473    3330    49271374        49271517        9.5     38.2
 gi|9626372|ref|NC_001422.1|     ScoY8t2_288;HRSCAF=336  45.714  35      19      0       3830    3726    78159550        78159446        2.4     39.6
 gi|9626372|ref|NC_001422.1|     ScoY8t2_247;HRSCAF=290  36.585  41      26      0       3948    4070    3226986 3226864 5.0     39.1
 gi|9626372|ref|NC_001422.1|     ScoY8t2_23291;HRSCAF=23585      43.902  41      23      0       4420    4298    14246463        14246585        0.21    43.7
 gi|9626372|ref|NC_001422.1|     ScoY8t2_5897;HRSCAF=6038        34.884  43      28      0       5369    5241    51644473        51644601        9.5     38.2
 ```
