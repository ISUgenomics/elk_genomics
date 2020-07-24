# Results


## QC
A lack of contamination was confirmed using both the initial Masurca assembly and final pseudomolecule assembly.  We aligned the Pacbio subreads to the genome, aligned contigs to the NT database with BLAST, and compiled this information with Blobtools.  All putative contaminant contigs were ruled out via manual investigation.  
Additional steps were taken to ensure the completeness of the genome, via mapping the input assembly reads to the genome.  90.7% of PacBio CCS reads aligned to the genome, while 87.3% of Illumina DNA-seq aligned.
### PhiX
No significant hits were identified to indicate contamination with PhiX

### BUSCO
Busco 4.0.2 was used to estimate the completeness of the genome assembly and annotation. BUSCO in genome mode predicted C:64.4%[S:62.0%,D:2.4%],F:3.1%,M:32.5%,n:255 and C:88.1%[S:86.0%,D:2.1%],F:2.1%,M:9.8%,n:13335 on eukaryota_odb10 and cetartiodactyla_odb10, respectively.  In protein mode, a remarkable score improvement is seen in both eukaryota and cetartiodactyla at C:97.7%[S:41.2%,D:56.5%],F:1.2%,M:1.1%,n:255 and C:92.1%[S:40.2%,D:51.9%],F:1.7%,M:6.2%,n:13335.
s

## Assembly
Without a genome sequence, identifying Elk genes that respond to Brucella infection would be impossible. Here we present the first pseudomolecule assembly of the rocky mountain elk, generated with xxx of sequencing data and xx fold coverage of the genome. An initial assembly of 23,302 contigs was generated using MaSuRCA with an N50 of xx and L90 of xxx contigs. This assembly was assembled into 35 psuedomolecules using 150 gbp of HiC sequencing and the Juicer, Juice Box, and the 3D-DNA pipeline. Assembly quality was verified via mapping illumina and Pacbio CCS reads to the assembly, with alignment percentages at 90.7% and 87.3% respectively.   


## Mitochondrial Genome

The mitochondrial genome is 16,429 bases long with 39 gene models, including tRNAs.  Three large mitochondrial integrations were identified in the genome ranging from 93-99% identical on Chromosome_02, Chromosome_12, and Chromosome_23.


## Gene Prediction
For gene annotation we generated 1.5 billion paired end reads of sequencing from six tissues, including kidney, lung, Mes, muscle, scap, and spleen. After masking repeat sequences using Repeatmodeler and Repeatmaker, we performed five de novo transcript/gene predictions with RNA-seq. The best transcripts were discerned using Mikado, followed by clustering with Cufflinks using Bos Taurus mRNAs to determine gene loci. Using this approach 18,960 genes were predicted to encode 33,433 mRNAs.

## Repeat prediction
Two repeat predictions were used on the genome, annotating common areas, but also distinct.  With EDTA, DNA transposons comprised the largest percentage of repeats in the genome, at 16%. Both Copia and Gypsy retrotransposons  were present at 1% of the genome, while Helitrons were nearly twice as large at 2% of the genome. While EDTA relies heavily on homology and transposon structure, Repeatmodeler/Repeatmasker uses homology and copy number to determine repeats in the genome.  Repeatmodeler  assessed xx percent of the genome as repetitive, with xx contributing the most to the genome's repetitiveness.

## Functional Annotation
16,853 of the 18,960 genes were functionally annotated with one of Interproscan or Blast to NR, NT, and Uniprot.  
## Synteny
All elk chromosomes were syntenic with all C. elaphus and B. Taurus chromosomes, with the exception of the Y chromosome (Fig circos). As has been seen in previous Cervus assemblies,  multiple pairs of chromosomes are tandemly fused in B. Taurus and vise-versa.  Previous reports indicate that 8 chromosomes of red deer were fused in B. Taurus, yet we identified twelve chromosomes of elk are tandemly fused into six Bos taurus chromosomes.
eight acrocentric deer chromosomes appear to be tandemly fused into four acrocentric bovine ones (for a typical example, see Fig. S3), and conversely, one acrocen-tric deer chromosome appears to be split to two acrocentric bovine ones (Fig. S1). One metacentric deer chromosome appears to be centrically split to two acrocentric ones in the B. taurus complement (Fig. 3). The structural relation between deer chromosomes 19 and 31 vs. bovine 1 (i.e., Ce19 and 31 vs. Bt1) as well as for Ce26 and 28 vs. Bt9 seemed more complex. In the case of “Ce19 and 31 vs. Bt1”, a tandem fusion of the two acrocentric C. elaphus chromo-somes and a translocation are combined (Fig. 4), whereas the tandem fusion of acrocentric chromosomes Ce26 and Ce28 in acrocentric Bt9 displayed a paracentric inversion in relation of Ce28 and Bt9 (Fig. S2). All structural relations of Ce and Bt chromosomes are summarized in Table 2. These observations are in complete accordance with the compara-tive cytogenetic analyses by Bonnet et al. (2001) between cattle and sika de
er (C. nippon)

Two chromosomal translocations were inferred between the two Cervus species, both having strong HiC support in Elk.  Chromosome_08 and Chromosome_33 of elk, comprised large portions of CM…40.1 and CM…15.1. With the majority of Chromosome_08 homologous to CM…15.1, a xx MB chunk seems to have moved between red deer’s CM…40.1 and CM…15.1. Another smaller chromosome translocation of XX MB has occurred between CM…29.1 and CM…10.1, attributed to Chromosome_22 and Chromosome_03 in elk.  Interestingly, both of these translocations are between chromosomes in elk that are fused chromosomes in Bos Taurus, NC…29 and NC…32. Because the B. taurus assembly was used to orient and join scaffolds in the red deer assembly, it is likely that these translocations are due to misassembly.   


 the larger translocation may be due to different outcomes of the fission/fusion that occurred between these chromosomes in most recent common Cervus ancestor.
