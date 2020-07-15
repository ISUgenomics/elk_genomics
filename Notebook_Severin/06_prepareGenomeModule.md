# Explore the Gene model statistics

* 07/15/2020
* /work/gif/severin/Olsen/06_prepareGenomeModules

Set up required variables for genomeModules to work

```
export GENMODgit="/work/gif/severin/isugif/genomeModules/"
export GENMOD="/work/gif/severin/isugif/"
```

Soft link the genome and gff3 file

```
ln -s /work/gif/archiveResults/Olsen/Elk/Cercan.fasta
ln -s /work/gif/archiveResults/Olsen/Elk/Cercan.gff3
```

Decided to just use the fasta2gff.pl script directly outside the container as the container is dated and giving me trouble

```
perl /work/gif/severin/isugif/genomeModules//bin/gff2fasta.pl Cercan.fasta Cercan.gff3 Cercan
```

Output

```
ls -lha
total 805M
drwxr-sr-x. 2 severin its-hpc-nova-gif   11 Jul 15 12:00 .
drwxr-sr-x. 4 severin its-hpc-nova-gif    4 Jul 15 11:33 ..
-rw-r--r--. 1 severin its-hpc-nova-gif  63M Jul 15 11:58 Cercan.cdna.fasta
-rw-r--r--. 1 severin its-hpc-nova-gif  63M Jul 15 11:58 Cercan.cds.fasta
-rw-r--r--. 1 severin its-hpc-nova-gif  74M Jul 15 11:58 Cercan.exon.fasta
lrwxrwxrwx. 1 severin its-hpc-nova-gif   47 Jul 15 11:36 Cercan.fasta -> /work/gif/archiveResults/Olsen/Elk/Cercan.fasta
-rw-r--r--. 1 severin its-hpc-nova-gif  48K Jul 15 1  1:57 Cercan.fasta.index
-rw-r--r--. 1 severin its-hpc-nova-gif 858M Jul 15 11:58 Cercan.gene.fasta
lrwxrwxrwx. 1 severin its-hpc-nova-gif   46 Jul 15 11:36 Cercan.gff3 -> /work/gif/archiveResults/Olsen/Elk/Cercan.gff3
-rw-r--r--. 1 severin its-hpc-nova-gif  22M Jul 15 11:58 Cercan.pep.fasta
-rw-r--r--. 1 severin its-hpc-nova-gif  56M Jul 15 11:58 Cercan.upstream3000.fasta
```


## Assembly Stats on the genes (start - stop)

```
assemblyStats Cercan.gene.fasta
singularity loaded
Number of Scaffolds:                 18960
Total Nucleotide content             884234743
Longest Scaffolds:                    2126538 1821443 1767020 1281368 1232833    >Cercan_00011563 >Cercan_00007552 >Cercan_00006234 >Cercan_00002203 >Cercan_00011057
Shortest Scaffolds:                  140 144 150 156 159         >Cercan_00006611 >Cercan_00002873 >Cercan_00010533 >Cercan_00015951 >Cercan_00011846
Mean Scaffold Size                   46636
Median Scaffold length               17744.5
N50 Scaffold length                  123631
L50 Scaffold length                  1723
N90 Scaffold length                  23381
L90 Scaffold length                  8183

                                      #Scaffs   % Scaffolds      Nucleotides     % Nucleotide Content
Number of Scaffolds [0-1K) nt         1494       7.879 %         917860                  0.103 %
Number of Scaffolds [1K-10K) nt       5453       28.76 %         26396320        2.985 %
Number of Scaffolds [10K-100K) nt     9733       51.33 %         352715761       39.88 %
Number of Scaffolds [100K-1M) nt      2267       11.95 %         486856990       55.05 %
Number of Scaffolds [1M-10M) nt       13         0.068 %         17347812        1.961 %
Number of Scaffolds > 10M nt          0          0.0 %   0       0.0 %

```
There are 18,960 gene models the longest where the stop - start position of the gene model in the genome is 2,126,538 bases (Cercan_00011563) and is annotated as a dystrophin gene.  This happens to be the largest gene in Humans as well.


## Assembly stats on the peptides

```
assemblyStats Cercan.pep.fasta
singularity loaded
Number of Scaffolds:                 33432
Total Nucleotide content             21343973
Longest Scaffolds:                    8822 8799 8785 7665 7613   >Cercan_00018031-R2 >Cercan_00018031-R1 >Cercan_00018031-R3 >Cercan_00014077-R5 >Cercan_00014077-R1
Shortest Scaffolds:                  18 19 22 22 35      >Cercan_00006611-R1 >Cercan_00011769-R1|CM008041_1_cds_OWJ99395_1_18858_mrna1 >Cercan_00009704-R1|CM008031_1_cds_OWK02783_1_15739_mrna1 >Cercan_00016988-R1 >Cercan_00012317-R1|CM008008_1_cds_OWK18193_1_373_mrna2
Mean Scaffold Size                   638
Median Scaffold length               479.0
N50 Scaffold length                  856
L50 Scaffold length                  7403
N90 Scaffold length                  329
L90 Scaffold length                  23127

                                      #Scaffs   % Scaffolds      Nucleotides     % Nucleotide Content
Number of Scaffolds [0-1K) nt         27943      83.58 %         12440699                58.28 %
Number of Scaffolds [1K-10K) nt       5489       16.41 %         8903274         41.71 %
Number of Scaffolds [10K-100K) nt     0          0.0 %   0       0.0 %
Number of Scaffolds [100K-1M) nt      0          0.0 %   0       0.0 %
Number of Scaffolds [1M-10M) nt       0          0.0 %   0       0.0 %
Number of Scaffolds > 10M nt          0          0.0 %   0       0.0 %
```
There are 33,432 proteins (from the isoforms) shortest proteins are 18-35 Amino Acids long and  the longest is 8822 AA (spectrin repeat containing nuclear envelope protein 1 (SYNE1)) [ie longest protein in AA content]

## Mikado gff stats

```
mikado util stats Cercan.gff3 --tab-stats Cercan.gff3.table > Cercan.gff.stats
```

|Stat|Total|Average|Mode|Min|1%|5%|10%|25%|Median|75%|90%|95%|99%|Max|  
|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|
|Number of genes|18079|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA
|Number of genes (coding)|18003|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA
|Number of monoexonic genes|1412|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA
|Transcripts per gene|33534|1.85|1|1|1|1|1|1|1|2|3|4|7|14
|Coding transcripts per gene|33401|1.85|1|0|1|1|1|1|1|2|3|4|7|14
|CDNA lengths|75396610|2,248.36|678|57|249|402|558|962|1,674|2,874|4,471|5,943|9,933|31,719
|CDNA lengths (mRNAs)|75199014|2,251.40|678|57|249|402|561|963|1,677|2,875|4,476|5,951|9,945|31,719
|CDS lengths|63981389|1,907.96|0|0|216|351|483|858|1,431|2,391|3,756|5,056|8,469|26,466
|CDS lengths (mRNAs)|NA|1,915.55|678|54|231|363|492|864|1,437|2,397|3,762|5,064|8,472|26,466
|CDS/cDNA ratio|NA|89.57|100.0|5|31|40|53|91|100|100|100|100|100|100
|Monoexonic transcripts|1928|942.94|225|57|201|240|290|438|756|1,129|1,767|2,321|3,939|14,659
|MonoCDS transcripts|2348|844.23|225|54|178|225|263|390|686|1,080|1,572|1,980|3,269|6,948
|Exons per transcript|400695|11.95|2|1|1|1|2|4|9|16|25|32|51|146
|Exons per transcript (mRNAs)|3345|11.97|2|1|1|1|2|4|9|16|25|32|51|146
|Exon lengths|NA|188.16|126|1|8|42|59|88|126|177|272|494|1,546|18,639
|Exon lengths (mRNAs)|NA|188.12|126|1|8|42|59|88|126|177|271|493|1,546|18,639
|Intron lengths|NA|4,157.15|88|1|76|93|133|406|1,221|3,240|8,286|15,612|55,221|490,824
|Intron lengths (mRNAs)|NA|4,154.18|88|1|76|93|133|405|1,221|3,239|8,278|15,596|55,188|490,824
|CDS exons per transcript|2313|11.68|2|0|1|1|2|4|9|16|25|31|51|145
|CDS exons per transcript (mRNAs)|2313|11.73|2|1|1|1|2|4|9|16|25|31|51|145
|CDS exon lengths|63981389|163.37|126|1|6|39|57|86|124|171|242|371|1,068|12,624
|CDS Intron lengths|1417693176|3,957.51|87|0|75|92|130|399|1,201|3,159|7,940|14,598|51,348|490,823
|5'UTR exon number|33401|0.46|0|0|0|0|0|0|0|1|2|2|4|28
|3'UTR exon number|33401|0.30|0|0|0|0|0|0|0|1|1|1|2|25
|5'UTR length|3055629|91.48|0|0|0|0|0|0|0|39|298|517|1,171|12,476
|3'UTR length|8161996|244.36|0|0|0|0|0|0|0|22|834|1,526|3,300|18,537
|Stop distance from junction|NA|11.07|0|0|0|0|0|0|0|0|0|0|234|6,905
|Intergenic distances|NA|93,327.44|61;96|-1,240,468|-76,699|-1,269|93|3,852|19,357|70,010|213,609|422,160|1,344,865|8,641,691
|Intergenic distances (coding)|NA|93,813.30|61;96|-1,240,468|-76,627|-1,252|96|3,862|19,465|70,470|215,214|428,042|1,351,464|8,641,691


Mikado ignores any genes that do not have transcripts (ie mRNA lines in the gff3 file) which is why it only reports the 18,079 out of the 18,960 gene models.
