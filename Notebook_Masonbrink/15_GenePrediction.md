# Gene prediction for Cervus canadensis

# Gather transcriptomic files
```
#Downloaded EST's for cervus mrna from NCBI's nucleotide database on 9/9/2019
#desktop/condo/ceres -->
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/CervusNucleotideESTsFromNCBI.fasta

#Downloaded Cervus Nippon Elaphus cds from NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/197/005/GCA_002197005.1_CerEla1.0/GCA_002197005.1_CerEla1.0_cds_from_genomic.fna.gz

#Downloaded Cervus Nippon Elaphus predicted proteins from NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/197/005/GCA_002197005.1_CerEla1.0/GCA_002197005.1_CerEla1.0_protein.faa.gz


#Downloaded all proteins from Trembl and swissprot, 22,346 proteins downloaded on 9/9/2019
#Desktop/condo/ceres -->
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/uniprot-cervus.fasta

#12 paired end rnaseq samples
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-spleen_S23_L004_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-spleen_S23_L004_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-spleen_S23_L003_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-spleen_S23_L003_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/ElkpscapLN_S22_L004_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/ElkpscapLN_S22_L004_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/ElkpscapLN_S22_L003_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/ElkpscapLN_S22_L003_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-muscle_S21_L004_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-muscle_S21_L004_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-muscle_S21_L003_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-muscle_S21_L003_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-Mes-LN_S24_L004_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-Mes-LN_S24_L004_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-Mes-LN_S24_L003_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-Mes-LN_S24_L003_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-lung_S26_L004_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-lung_S26_L004_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-lung_S26_L003_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-lung_S26_L003_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-kidney_S25_L004_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-kidney_S25_L004_R1_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-kidney_S25_L003_R2_001.fastq
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/16_RNAseq/Elk-kidney_S25_L003_R1_001.fastq
```

### Braker gene prediction for mikado
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker

# Grabbing previously extracted rnaseq
for f in  ../16_RNAseq/*fastq; do ln -s $f;done

ln -s ../16_RNAseq/GCA_002197005.1_CerEla1.0_protein.faa.gz
ln -s ../16_RNAseq/uniprot-cervus.fasta
ln -s ../16_RNAseq/GCA_002197005.1_CerEla1.0_cds_from_genomic.fna

mkdir 01_AlignRNA




module load braker/2.1.2

braker --species=CervusCanadensis --gff3 --cores 40 --genome=SCNgenome.fasta --bam=SCNgenome.consensus_isoforms_sorted.bam,SCNgenome.H.glycinesEST_sorted.bam,AllRNASEQ_sorted.bam

```

### Portcullis for gene prediciton
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/18_PortCullis
#uses lots of memory, so beware
module load miniconda
source activate portcullis
portcullis full --threads 9 --verbose --use_csi --output portcullis_out --orientation FR SCNgenome.fasta NonRiboReads.bam
```

### Stringtie transcript assembly for gene prediction
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/19_Stringtie

module load stringtie
stringtie NonRiboReads.bam -j 5 -p 16 -v -o NonRrnaRNASEQ_stringtie.gtf

```

### run class
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/20_Class2

module load miniconda
source activate Class2
run_class.pl -a AllRNASEQ_sorted.bam -o AllRNASEQClass2.gtf -p 40 --verbose

```

### Trinity assembly for gene prediction
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/21_Trinity

sh runTrinity.sh NonRiboReads.bam


#runTrinity.sh
################################################################################
#!/bin/bash

module load trinityrnaseq/2.8.4

bam="$1"
out=$(basename ${bam%.*} |cut -f 1 -d "_")
Trinity \
   --genome_guided_bam ${bam} \
   --max_memory 120G \
   --genome_guided_max_intron 30000 \
   --full_cleanup \
--CPU 40
################################################################################

#finished in 32hrs
#how many transcripts did we get?
grep -c ">" Trinity-GG.fasta
110683

#map trinity transcripts back to genome
module load gmap_gsnap/2017-03-17
sh runGmap.sh SCNgenome /work/GIF/remkv6/Baum/04_Dovetail2Restart/12_Trinity/trinity_out_dir/ ../SCNgenome.fasta Trinity-GG.fasta
###############################################################################
#!/bin/bash

#Makes a database and searches your sequences.
#sh runGmap.sh <database name> <folder of database file ending with a "/"> <Fasta file> <query file>

#examples
#sh run_gmap.sh red_abalone_02Jun2017_5fUJu /work/GIF/remkv6/Serb/03_DavideGMAP/ red_abalone_02Jun2017_5fUJu.fasta DavideQuerydna.fasta
#sh run_gmap.sh  m.yessoensisGenome /work/GIF/remkv6/Serb/03_DavideGMAP DavideQuerydna.fasta
#sh run_gmap.sh Crassostreagigasgenome /work/GIF/remkv6/Serb/03_DavideGMAP Crassostreagigasgenome.fa DavideQuerydna.fasta


module load gmap_gsnap/2017-03-17
dbname=$1
dbloc=$2
dbfasta=$3
query=$4
gmap_build -d $dbname  -D $dbloc $dbfasta
gmap -D $dbloc -d $dbname -B 5 -t 40  --input-buffer-size=1000000 --output-buffer-size=1000000 -f gff3_gene  $query >${dbname%.*}.${query%.*}.gff3
################################################################################
```

### RepeatMasker of new genome
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/22_RepeatMasker2

ln -s ../03_repeatmodeler/consensi.fa.classified

###############################################################################
#!/bin/bash
# runs repeat masking for the genome after constructing custom repeat library
# uses repeat modeler for building custom db and RepeatMasking for masking
# run it as:
# runRepeatModeler.sh Genome.fasta
# based on Rick's guide https://intranet.gif.biotech.iastate.edu/doku.php/people:remkv6:genome738polished_repeatmodeler_--de_novo_repeat_identification

if [ $# -lt 1 ] ; then
        echo "usage: runRepeatModeler <genome.fasta>"
        echo ""
        echo "To build custom repeat library and mask the repeats of the genome"
        echo ""
exit 0
fi


GENOME="$1"
module purge
module load parallel

module load repeatmasker/4.0.7
#module load repeatmodeler/1.0.8
#DATABASE="$(basename ${GENOME%.*}).DB"
#BuildDatabase -name ${DATABASE} -engine ncbi ${GENOME}
#RepeatModeler -database ${DATABASE}  -engine ncbi -pa 16
#ln -s $(find $(pwd) -name "consensi.fa.classified")
RepeatMasker -pa 40 -gff -lib consensi.fa.classified ${GENOME}
#################################################################################

```
