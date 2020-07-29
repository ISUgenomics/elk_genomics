# Results

Collecting the results and final files we will need for the paper and Jbrowse.

* Nova: /work/gif/archiveResults/Olsen/Elk

## Final Genome and GFF files are located here:


* /work/gif/remkv6/Olsen/Elk/05a_RenameScaffsNGenes/GTfixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3
* /work/gif/remkv6/Olsen/Elk/05a_RenameScaffsNGenes/FinalGenomePilonReducedSoftMaskedRecode.fa

```
cp /work/gif/remkv6/Olsen/Elk/05a_RenameScaffsNGenes/GTfixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 Cercan.gff3
cp /work/gif/remkv6/Olsen/Elk/05a_RenameScaffsNGenes/FinalGenomePilonReducedSoftMaskedRecode.fa Cercan.fasta
```



## Assembly statistics on Genome assembly

```
assemblyStats Cercan.fasta > Cercan.stats
loading singularity
Number of Scaffolds:                 185
Total Nucleotide content             2526613007
Longest Scaffolds:                    146388637 127605827 114865875 114606702 105318381          >Chromosome_X >Chromosome_05 >Chromosome_20 >Chromosome_18 >Chromosome_09
Shortest Scaffolds:                  1015 1024 1061 1088 1138    >Scaffold_091 >Scaffold_090 >Scaffold_082 >Scaffold_072 >Scaffold_150
Mean Scaffold Size                   13657367
Median Scaffold length               9169.0
N50 Scaffold length                  77654944
L50 Scaffold length                  13
N90 Scaffold length                  51438166
L90 Scaffold length                  29

                                      #Scaffs   % Scaffolds      Nucleotides     % Nucleotide Content
Number of Scaffolds [0-1K) nt         0          0.0 %           0               0.0 %
Number of Scaffolds [1K-10K) nt       94         50.81 %         393769          0.015 %
Number of Scaffolds [10K-100K) nt     56         30.27 %         1472118         0.058 %
Number of Scaffolds [100K-1M) nt      0          0.0 %           0               0.0 %
Number of Scaffolds [1M-10M) nt       1          0.540 %         7618728         0.301 %
Number of Scaffolds > 10M nt          34         18.37 %         2517128392      99.62 %
```

## Gene model Statistics

```
more Cercan.gff3 | awk '$3=="gene"' | wc
  18960  170640 1644756
more Cercan.gff3 | awk '$3=="mRNA"' | wc
  33433 1417723 18296056
```
