# Raw data exploration and QC

* /project/elk_bison_genomics/Severin/00_rawdata
* 2019 June 12

```
mkdir 00_rawdata
cp ../final.genome.scf.fasta.bz2 .
module load bzip2
bunzip2 final.genome.scf.fasta.bz2
```


## Assembly statistics of the genome
These are the summary statistics of the genome to get a sense of the continuity of the assembly based on scaffolds and contigs. Dr Zimin made this assembly. I will need to ask him for methods.

```
perl  ~/isu_gif_vrsc/isugif/utilities/utilities/new_Assemblathon.pl final.genome.scf.fasta
---------------- Information for assembly 'final.genome.scf.fasta' ----------------


                                         Number of scaffolds      29125
                                     Total size of scaffolds 2559821282
                                            Longest scaffold    7681499
                                           Shortest scaffold        300
                                 Number of scaffolds > 1K nt       9934  34.1%
                                Number of scaffolds > 10K nt       6499  22.3%
                               Number of scaffolds > 100K nt       3320  11.4%
                                 Number of scaffolds > 1M nt        786   2.7%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size      87891
                                        Median scaffold size        441
                                         N50 scaffold length    1224689
                                          L50 scaffold count        617
                                                 scaffold %A      29.11
                                                 scaffold %C      20.87
                                                 scaffold %G      20.87
                                                 scaffold %T      29.08
                                                 scaffold %N       0.06
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      37.2%
              Percentage of assembly in unscaffolded contigs      62.8%
                      Average number of contigs per scaffold        1.0
Average length of break (>25 Ns) between contigs in scaffold       1393

                                           Number of contigs      30285
                              Number of contigs in scaffolds       2104
                          Number of contigs not in scaffolds      28181
                                       Total size of contigs 2558204714
                                              Longest contig    7681499
                                             Shortest contig        300
                                   Number of contigs > 1K nt      11091  36.6%
                                  Number of contigs > 10K nt       7561  25.0%
                                 Number of contigs > 100K nt       3869  12.8%
                                   Number of contigs > 1M nt        705   2.3%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size      84471
                                          Median contig size        457
                                           N50 contig length     939833
                                            L50 contig count        760
                                                   contig %A      29.13
                                                   contig %C      20.88
                                                   contig %G      20.89
                                                   contig %T      29.10
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
```


## untar the dovetail folders that were provided

* /home/andrew.severin/elk_bison_genomics/Severin/00_rawdata/dovetail
* Need to check to see if we have the raw data as well from the Hi-C Illumina runs

```
ln -s ../../../dovetail/elk_dovetail_delivery_WIL462_KF0fHeNqkz.tar
ln -s ../../../dovetail/elk_dovetail_delivery_WIL462_SQND8mmmSa.tar
tar -xvf elk_dovetail_delivery_WIL462_KF0fHeNqkz.tar &
tar -xvf elk_dovetail_delivery_WIL462_SQND8mmmSa.tar &

gunzip elk_dovetail_delivery_WIL462_KF0fHeNqkz/north_american_elk_11Jun2018_1UW0F.fasta.gz
gunzip elk_dovetail_delivery_WIL462_SQND8mmmSa/north_american_elk_15Jun2018_oY8t2.fasta.gz

```


## Assembly statistics after Dovetail

* /home/andrew.severin/elk_bison_genomics/Severin/00_rawdata/dovetail/delivery_WIL462_KF0fHeNqkz
* north_american_elk_11Jun2018_1UW0F.fasta

```
perl  ~/isu_gif_vrsc/isugif/utilities/utilities/new_Assemblathon.pl north_american_elk_11Jun2018_1UW0F.fasta
---------------- Information for assembly 'north_american_elk_11Jun2018_1UW0F.fasta' ----------------


                                         Number of scaffolds      23758
                                     Total size of scaffolds 2560425402
                                            Longest scaffold   55855843
                                           Shortest scaffold        300
                                 Number of scaffolds > 1K nt       4567  19.2%
                                Number of scaffolds > 10K nt       2015   8.5%
                               Number of scaffolds > 100K nt        525   2.2%
                                 Number of scaffolds > 1M nt        244   1.0%
                                Number of scaffolds > 10M nt         80   0.3%
                                          Mean scaffold size     107771
                                        Median scaffold size        394
                                         N50 scaffold length   17659217
                                          L50 scaffold count         40
                                                 scaffold %A      29.09
                                                 scaffold %C      20.86
                                                 scaffold %G      20.88
                                                 scaffold %T      29.08
                                                 scaffold %N       0.09
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      96.6%
              Percentage of assembly in unscaffolded contigs       3.4%
                      Average number of contigs per scaffold        1.3
Average length of break (>25 Ns) between contigs in scaffold        307

                                           Number of contigs      30975
                              Number of contigs in scaffolds       7715
                          Number of contigs not in scaffolds      23260
                                       Total size of contigs 2558203134
                                              Longest contig    5383423
                                             Shortest contig        180
                                   Number of contigs > 1K nt      11777  38.0%
                                  Number of contigs > 10K nt       8293  26.8%
                                 Number of contigs > 100K nt       4375  14.1%
                                   Number of contigs > 1M nt        635   2.1%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size      82589
                                          Median contig size        468
                                           N50 contig length     782806
                                            L50 contig count        928
                                                   contig %A      29.12
                                                   contig %C      20.88
                                                   contig %G      20.90
                                                   contig %T      29.11
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
```

I suspect that this is the chicago only assembly as it didn't really pull that much together.


* /home/andrew.severin/elk_bison_genomics/Severin/00_rawdata/dovetail/delivery_WIL462_SQND8mmmSa


```
perl  ~/isu_gif_vrsc/isugif/utilities/utilities/new_Assemblathon.pl

---------------- Information for assembly 'north_american_elk_15Jun2018_oY8t2.fasta' ----------------


                                         Number of scaffolds      23302
                                     Total size of scaffolds 2560463354
                                            Longest scaffold  128024478
                                           Shortest scaffold        300
                                 Number of scaffolds > 1K nt       4111  17.6%
                                Number of scaffolds > 10K nt       1559   6.7%
                               Number of scaffolds > 100K nt         93   0.4%
                                 Number of scaffolds > 1M nt         40   0.2%
                                Number of scaffolds > 10M nt         37   0.2%
                                          Mean scaffold size     109882
                                        Median scaffold size        391
                                         N50 scaffold length   76232311
                                          L50 scaffold count         14
                                                 scaffold %A      29.11
                                                 scaffold %C      20.86
                                                 scaffold %G      20.88
                                                 scaffold %T      29.07
                                                 scaffold %N       0.09
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      97.8%
              Percentage of assembly in unscaffolded contigs       2.2%
                      Average number of contigs per scaffold        1.3
Average length of break (>25 Ns) between contigs in scaffold        297

                                           Number of contigs      30899
                              Number of contigs in scaffolds       7739
                          Number of contigs not in scaffolds      23160
                                       Total size of contigs 2558203086
                                              Longest contig    5383423
                                             Shortest contig        180
                                   Number of contigs > 1K nt      11703  37.9%
                                  Number of contigs > 10K nt       8220  26.6%
                                 Number of contigs > 100K nt       4325  14.0%
                                   Number of contigs > 1M nt        642   2.1%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size      82792
                                          Median contig size        467
                                           N50 contig length     794616
                                            L50 contig count        921
                                                   contig %A      29.13
                                                   contig %C      20.88
                                                   contig %G      20.89
                                                   contig %T      29.10
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
```

N50's go from 1.2M to 17.6M to 76M.  So there was a massive improvement in the genome assembly.  There must be a lot of small unplaced contigs or potentially contamination that can be removed.


This assembly is the current best assembly **north_american_elk_15Jun2018_oY8t2.fasta**
