Rick E Masonbrink, Darrel O. Bayles, David Alt, Aleksey Zimin, Fred Tatum, Paola Boggiatto, Jenny Wilson-Welder, Jeffrey Williams, Hand Edwards, Mary Wood, Andrew Severin, Steven Olsen

# Results


## QC
A lack of contamination was confirmed using both the initial Masurca assembly and final pseudomolecule assembly.  We aligned the Pacbio subreads to the genome, aligned contigs to the NT database with BLAST, and compiled this information with Blobtools(Supplemental).  All putative contaminant contigs were ruled out via manual investigation.  
Additional steps were taken to ensure the completeness of the genome, via mapping the input assembly reads to the genome.  90.7% of PacBio CCS reads aligned to the genome, while 87.3% of Illumina DNA-seq aligned (Table SeqNCoverage).
### PhiX
No significant hits were identified to indicate contamination with PhiX

### BUSCO
Busco 4.0.2 was used to estimate the completeness of the genome assembly and annotation. BUSCO in genome mode predicted C:64.4%[S:62.0%,D:2.4%],F:3.1%,M:32.5%,n:255 and C:88.1%[S:86.0%,D:2.1%],F:2.1%,M:9.8%,n:13,335 on eukaryota_odb10 and cetartiodactyla_odb10, respectively.  In protein mode, a remarkable score improvement is seen in both eukaryota and cetartiodactyla at C:97.7%[S:41.2%,D:56.5%],F:1.2%,M:1.1%,n:255 and C:92.1%[S:40.2%,D:51.9%],F:1.7%,M:6.2%,n:13,335.

## Assembly
To improve conservation of the Rocky Mountain Elk (C. canadensis), a genome sequence is needed to improve …breeding and identify and assess …  herd diversity.  Here we present the first pseudomolecule assembly of the C. canadensis, generated with 1.7 trillion base pairs of sequencing at a 686-fold coverage of the genome. An initial assembly of 23,302 contigs was generated using MaSuRCA with an N50 of xx and L90 of xxx contigs. This assembly was scaffolded into 35 psuedomolecules leveraging 150 gbp of HiC sequencing with the Juicer/Juice Box/3D-DNA pipeline. A high quality assembly was ascertained via the mapping of Illumina and Pacbio CCS reads to the assembly, with alignment percentages at 90.7 and 87.3% respectively(Table SeqNcoverage).   


## Mitochondrial Genome
The mitochondrial genome is 16,429 bases long with 39 gene models, including tRNAs.  Three large mitochondrial integrations were identified in the genome ranging from 93-99% identical on Chromosome_02, Chromosome_12, and Chromosome_23.

## Gene Prediction
For gene annotation we generated 1.5 billion paired end reads of sequencing from six tissues, including kidney, lung, mesenteric lymph node, muscle, prescapular lymph node, and spleen. After masking repeat sequences using Repeatmodeler and Repeatmaker, we performed five de novo transcript/gene predictions with RNA-seq (InputTranscriptomes). The best transcripts were discerned using Mikado, followed by clustering with Cufflinks using Bos Taurus mRNAs to cluster transcripts into gene loci. Using this approach 18,013 genes were predicted to encode 33,433 mRNAs (InputTranscriptomes).

## Repeat prediction
Two repeat predictions were used on the genome, annotating common areas, but also distinct (Repeats).  While EDTA relies heavily on homology and transposon structure, Repeatmodeler/Repeatmasker uses homology and copy number to determine repeats in the genome.  With EDTA, 25.8% of the genome was marked repetitive, with DNA transposons comprised the largest percentage of repeats in the genome, at 16%. In contrast, Repeatmasker assessed 36.5% of the genome as an interspersed repeat, with 28.8% of the genome being comprised LINE retrotransposons. We merged these repeat annotations with Bedtools  to reveal that 38% of the genome is repetitive. This large repetitive content as well as assembly methods are likely responsible for the differences in chromosome sizes between red deer and elk, as repetitive content in the red deer was estimated at 22.7% and gaps comprised more than 1.5Gbp of the assembly.  

## Functional Annotation
17,938 of the 18,013 genes or 99.6% were functionally annotated with one of Interproscan or BLAST to NR, NT, and Uniprot (FunctionalAnnotationsTable).  
## Synteny
All elk chromosomes were syntenic with all C. elaphus and B. Taurus chromosomes, though chromosome Y lacked the genes required for gene-based synteny (Fig circos). As has been seen in previous Cervus assemblies,  multiple pairs of chromosomes are tandemly fused in B. Taurus and vise-versa.  We confirmed previous reports chromosome fusions and fissions indicated that 12 cervus chromosomes fused into six in B. taurus, as well as 4 chromosomes in B. Taurus are fused into two cervus chromosomes (table).  

Two inter-chromosomal translocations were inferred between the two Cervus species, both having strong HiC support in Elk (table, plot). Chromosome_15 and Chromosome_24 of elk, comprised large portions of Ce_Chr_33 and a minor portion of C. elaphus Ce_Chr_8. With the majority of Chromosome_24 homologous to C. elaphus Ce_Chr_8, a 17 MB region of Ce_Chr_33 may have been falsely attached to Ce_Chr_8 in C. elaphus. Another smaller chromosome translocation of 13.6 MB occurred between Ce_Chr_22 and Ce_Chr_3  of C. elaphus, attributed to chromosomes 21 and 25 in C. canadensis.  A small region of 22 was likely falsely attached to Ce_Chr_3 in C. elaphus.  Interestingly, both of these translocations are between chromosomes in elk that are fused chromosomes in Bos Taurus, Bt_Chr_2 and Bt_Chr_5 (table). While it is possible that these translocations occurred since the divergence of these two species, because the B. taurus assembly was used to orient and join scaffolds in the C. elaphus genome assembly, it is likely that these translocations are misassemblies in the C. elaphus genome.
