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

### RepeatMasker of new genome
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/02_RepeatMasker
ln -s ../../26_RefineLittle2Big3/FinalGenomePilonReduced.fa
ln -s ../../03_repeatmodeler/consensi.fa.classified

sh runRepeatModeler.sh FinalGenomePilonReduced.fa

#runRepeatModeler
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

#rename the scaffolds and get rid of special characters
sed 's/_//g' FinalGenomePilonReduced.fa.masked |sed 's/pilon//g' >FinalGenomePilonReducedRenamedMasked.fa
```



### Perform alignments of expression data
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/06_AlignRNAMasked

for f in  ../*fastq; do ln -s $f;done
ln -s ../02_RepeatMasker/FinalGenomePilonReducedRenamedMasked.fa

sh runHISAT2.sh Elk-kidney_S25_L003_R1_001.fastq Elk-kidney_S25_L003_R2_001.fastq
sh runHISAT2.sh Elk-kidney_S25_L004_R1_001.fastq Elk-kidney_S25_L004_R2_001.fastq
sh runHISAT2.sh Elk-lung_S26_L003_R1_001.fastq Elk-lung_S26_L003_R2_001.fastq
sh runHISAT2.sh Elk-lung_S26_L004_R1_001.fastq Elk-lung_S26_L004_R2_001.fastq
sh runHISAT2.sh Elk-Mes-LN_S24_L003_R1_001.fastq Elk-Mes-LN_S24_L003_R2_001.fastq
sh runHISAT2.sh Elk-Mes-LN_S24_L004_R1_001.fastq Elk-Mes-LN_S24_L004_R2_001.fastq
sh runHISAT2.sh Elk-muscle_S21_L003_R1_001.fastq Elk-muscle_S21_L003_R2_001.fastq
sh runHISAT2.sh Elk-muscle_S21_L004_R1_001.fastq Elk-muscle_S21_L004_R2_001.fastq
sh runHISAT2.sh ElkpscapLN_S22_L003_R1_001.fastq ElkpscapLN_S22_L003_R2_001.fastq
sh runHISAT2.sh ElkpscapLN_S22_L004_R1_001.fastq ElkpscapLN_S22_L004_R2_001.fastq
sh runHISAT2.sh Elk-spleen_S23_L003_R1_001.fastq Elk-spleen_S23_L003_R2_001.fastq
sh runHISAT2.sh Elk-spleen_S23_L004_R1_001.fastq Elk-spleen_S23_L004_R2_001.fastq
sh runHISAT2.sh Undetermined_S0_L003_R1_001.fastq Undetermined_S0_L003_R2_001.fastq
sh runHISAT2.sh Undetermined_S0_L004_R1_001.fastq Undetermined_S0_L004_R2_001.fastq

#runHISAT2.sh
################################################################################
#!/bin/bash

module load hisat2
module load samtools
DBDIR="/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/06_AlignRNAMasked"
GENOME="FinalGenomePilonReducedRenamedMasked.fa"

p=40
R1_FQ="$1"
R2_FQ="$2"


hisat2 \
  -p ${p} \
  -x ${DBDIR}/${GENOME} \
  -1 ${R1_FQ} \
  -2 ${R2_FQ}  \
  -S  ${R1_FQ}.sam &> ${R1_FQ}.log
samtools view --threads 40 -b -o ${R1_FQ}.bam ${R1_FQ}.sam
mkdir ${R1_FQ}_temp
samtools sort  -o ${R1_FQ}_sorted.bam -T ${R1_FQ}_temp --threads 40 ${R1_FQ}.bam
samtools index ${R1_FQ}_sorted.bam
################################################################################
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
for f in ../../16_RNAseq/*fastq; do ln -s $f;done


module load braker/2.1.2

#copy the augustus config folder so it is writeable
cp -rf /software/7/apps/augustus/3.3.2/config/ .

 echo "module load miniconda;source activate Braker; braker --species=CervusCanadensis -gff3 ----useexisting cores 40 --genome=FinalGenomePilonReducedRenamedMasked.fa --bam=Elk-kidney_S25_L003_R1_001.fastq_sorted.bam,Elk-kidney_S25_L004_R1_001.fastq_sorted.bam,Elk-lung_S26_L003_R1_001.fastq_sorted.bam,Elk-lung_S26_L004_R1_001.fastq_sorted.bam,Elk-Mes-LN_S24_L003_R1_001.fastq_sorted.bam,Elk-Mes-LN_S24_L004_R1_001.fastq_sorted.bam,Elk-muscle_S21_L003_R1_001.fastq_sorted.bam,Elk-muscle_S21_L004_R1_001.fastq_sorted.bam,ElkpscapLN_S22_L003_R1_001.fastq_sorted.bam,ElkpscapLN_S22_L004_R1_001.fastq_sorted.bam,Elk-spleen_S23_L003_R1_001.fastq_sorted.bam,Elk-spleen_S23_L004_R1_001.fastq_sorted.bam,Undetermined_S0_L003_R1_001.fastq_sorted.bamm --AUGUSTUS_CONFIG_PATH=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/05_BrakerRun/config/" >braker.sh
```
Lots of problems with braker installations/file formatting issues.  A single bam file will progress through the entire prediction, but all files do not.  I will try all rna-seq alignments only and eliminate the transcript alignments.  


### Portcullis for gene prediciton
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/18_PortCullis
#uses lots of memory, so beware
module load miniconda
source activate portcullis
echo "ml miniconda;source activate portcullis; portcullis full --threads 36 --verbose --use_csi --output portcullis_out --orientation FR FinalGenomePilonReduced.fa AllRNASEQ_sorted.bam" >portcullis.sh
```

### Stringtie transcript assembly for gene prediction
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/19_Stringtie

module load stringtie
ln -s ../17_Braker/01_AlignRNA/AllRNASEQ_sorted.bam
ln -s ../17_Braker/01_AlignRNA/AllRNASEQ_sorted.bam.bai
ln -s ../17_Braker/01_AlignRNA/FinalGenomePilonReduced.fa

stringtie AllRNASEQ_sorted.bam -j 5 -p 40 -v -o FinalGenomePilonReduced.gtf

```

### run class
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/20_Class2

module load miniconda
source activate Class2
echo "ml miniconda; source actiate Class2; run_class.pl -a AllRNASEQ_sorted.bam -o AllRNASEQ_Class2.gtf -p 40 --verbose" > class2.sh
```

### Trinity assembly for gene prediction
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/21_Trinity

ln -s ../17_Braker/01_AlignRNA/AllRNASEQ_sorted.bam
ln -s ../17_Braker/01_AlignRNA/AllRNASEQ_sorted.bam.bai
echo "module load trinityrnaseq/2.8.4; sh runTrinity.sh AllRNASEQ_sorted.bam" >trinity.sh


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

### run strawberry
```
#/project/elk_bison_genomics/Masonbrink/23_strawberry/strawberry

ln -s ../17_Braker/01_AlignRNA/AllRNASEQ_sorted.bam
ln -s ../17_Braker/01_AlignRNA/AllRNASEQ_sorted.bam.bai
ln -s ../17_Braker/01_AlignRNA/FinalGenomePilonReduced.fa

ml minioonda
source activate strawberry
strawberry  -T strawberry.log  --no-quant -p 40 -v AllRNASEQ_sorted.bam

echo "ml miniconda; source activate strawberry; strawberry  -T strawberry.log  --no-quant -p 40 -v AllRNASEQ_sorted.bam" >strawberry.sh


```
