# Methods


## Sample Collection and Sequencing

## Genome QC

The phi-X174 nucleotide sequence obtained from NCBI (NC_001422.1) was BLASTed (blastn and tblastx version 2.6.0) against the assembled genome with default settings.



### Genome Assembly



1,004,453,472 HI-C reads were used for scaffolding the initial assembly with the Juicer 1.5.6, 3D-DNA 180922, and JuiceBox 1.9.8 pipeline. Reads were extracted from bam files with Picard 2.9.2. The initial assembly was masked using Repeatmodeler 4.0.7 and Repeatmasker 1.0.8, prior to the alignment of HI-C reads with BWA mem 0.7.17.  Alignments were processed using Juicer, 3D-DNA, and Juicebox. For the juicebox assembly strategy, all contigs greater than 10kb were placed manually, scaffolds were only incorporated in regions with the highest HI-C signal, scaffolds were placed in non-repetitive regions when HI-C signal was even between a repetitive and non-repetitive regions, repeats were clustered whenever possible, and only obvious mis-joins were edited, a total of ~400 scaffolds. The initial Juicebox scaffolding created 34 pseudomolecules, which was then compared to the Cervus elaphus genome assembly to reveal the merger of the X and Y chromosomes.  Using BLASTn of the Cervus elaphus genome sequence (GCA_002197005.1) to create coordinates, allowed the correct separation the X and Y chromosome via the heatmap in Juicebox.  The 3D-DNA assembly finished with 22,557 scaffolds.  
The contigs that could not be integrated into the pseudomolecules, 22,522 contigs, were examined  and eliminated based on repetitiveness, heterozygous haplotigs, RNA-seq mapping potential, and contig size. Bedtools 2.25.0 was used to merge coordinates from mapping these contigs to the pseudomolecules with BLAST+ 2.9 (score >300) and Repeatmasker 1.0.8 masking coordinates. 22,065 contigs were eliminated that were less than 1kb, had at least 90% query coverage, and lacked a single unique mapping RNA-seq read, leaving 457 contigs, 35 scaffolds, and the mitochondrial genome.  
The assembly was polished with Pilon 1.23 using CCS PacBio reads and paired end Illumina DNA-seq.  CCS PacBio reads were created from the PacBio subreads using bax2bam (https://github.com/PacificBiosciences/bax2bam) and Bamtools 2.5.1 and then aligned using minimap 2.6. Paired end reads were aligned using Hisat2 2.0.5, followed by bam conversion and sorting with Samtools 1.9. Due to uneven and excessive coverage in repetitive regions, paired end alignments were set at a max coverage of 30x using jvarkit (http://lindenb.github.io/jvarkit/). Due to the excessive repetitiveness of Scaffold 18, 214/264MB on scaffold_18 were not polished.
After polishing, another round of small contig elimination was performed by merging repeatmasker coordinates and coordinates from BLAST+ 2.9 (score >300) to the pseudomolecules with Bedtools 2.25.0. If 90% of query length was repetitive and contained within the pseudomolecules, it was eliminated, leaving 151 contigs, 35 pseudomolecules, and the mitochondrial contig.

### Identification of Mitochondria
BLAST+ 2.9 was used to identify the mitochondrial genome in the assembly, by querying the mitochondrial scaffold of the Cervus elaphus GCA_002197005.1. Scaffold HiC_scaffold_58 was identified as the Mitochondrial genome.  Due to the long read technology, this scaffold was a rolling assembly that contained 3 copies of the MT genome attached end to end.  The scaffold was manually corrected with Samtools 1.9 and renamed to XXX. The mitochondrial genome was annotated for genes using Mitos2 webserver using RefSeq 89 Metazoa, genetic code of 2, and default settings.  

### Gene model prediction
753,228,475 rnaseq reads aligned to the genome using hisat2 2.0.5 followed by bam conversion and sorting with Samtools 1.9. RNA-seq read counts were obtained using Subread 1.5.2. The alignments were assembled into genome-guided transcriptomes using Trinity 2.8.4, Strawberry 1.1.1, Stringtie 1.3.3b, and Class2 2.1.7.  The RNA-seq alignments were also used for a gene prediction via Braker2 2.1.4 with Augustus 3.3.3 on a genome softmasked by Repeatmasker 1.0.8 with a Repeatmodeler 4.0.7 library. High confidence exon splicing junctions were identified using Portcullis 1.1.2. Each of these assemblies was then supplied to Mikado 2.0rc6 to pick consensus transcripts, while utilizing Cervus-specific proteins from Uniprot (downloaded 12-28-19).  This mikado prediction was filtered for transposable elements using bedtools 2.25.0 intersect and filtered for pseudogenes via removing genes with 5 or fewer mapping RNA-seq reads.  Using bedtools2.25.0  intersect these filtered Mikado genes were used to find corresponding Braker2 2.1.4 genes. Both of these predictions, together with a Genomethreader 1.7.1 alignment of Uniprot proteins from the Pecora infraorder (downloaded 02-07-20) were used for a final round of Mikado gene prediction.     


### Functional Gene Annotation
Predicted transcript and protein fasta files were generated using Cufflinks gffread (2.2.1). Functional gene annotations were compiled using interproscan 5.27-66.0 and BLAST searches to NCBI NT and NR databases downloaded on 10-23-19, as well as swissprot/uniprot databases downloaded on 12/09/2019.

### Repeat prediction
A final version of predicted repeats was obtained using EDTA 1.7.9.


### BUSCO
Universal single copy orthologs were assessed using BUSCO 4.0, with the eukaryota_odb10 dataset.
