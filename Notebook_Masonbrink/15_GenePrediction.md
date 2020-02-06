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


#Found a faster way to softmask the genome
module load repeatmasker/4.0.7
less FinalGenomePilonReduced.fa.out.gff |awk '{if($5>$4) {print $1,$4,$5} else {print $1,$5,$4}}' |tr " " "\t">
repeats.bed
 bedtools maskfasta -soft -fi FinalGenomePilonReduced.fa -bed repeats.bed -fo FinalGenomePilonReducedSoftMasked.fa


```



##### Perform alignments of expression data -- needed to realign later
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

# Merge all of the bam files

module load samtools; samtools merge -@ 40 AllRNASEQ.bam Elk-kidney_S25_L003_R1_001.fastq_sorted.bam Elk-kidney_S25_L004_R1_001.fastq_sorted.bam Elk-lung_S26_L003_R1_001.fastq_sorted.bam Elk-lung_S26_L004_R1_001.fastq_sorted.bam Elk-Mes-LN_S24_L003_R1_001.fastq_sorted.bam Elk-Mes-LN_S24_L004_R1_001.fastq_sorted.bam Elk-muscle_S21_L003_R1_001.fastq_sorted.bam Elk-muscle_S21_L004_R1_001.fastq_sorted.bam ElkpscapLN_S22_L003_R1_001.fastq_sorted.bam ElkpscapLN_S22_L004_R1_001.fastq_sorted.bam Elk-spleen_S23_L003_R1_001.fastq_sorted.bam Elk-spleen_S23_L004_R1_001.fastq_sorted.bam Undetermined_S0_L003_R1_001.fastq_sorted.bam Undetermined_S0_L004_R1_001.fastq_sorted.bam
mkdir TEMP
samtools sort -@ 40 -T TEMP -o AllRNASEQ_sorted.bam AllRNASEQ.bam;samtools index AllRNASEQ_sorted.bam


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
Lots of problems with braker installations/file formatting issues.  A single bam file will progress through the entire prediction, but all files do not.  I will try all rna-seq alignments only and eliminate the transcript alignments.  Still no go, genemark fails on the ET_C_2 step with empty stop codon files.


### Align RNA-seq reades with softmasked genome, recoded to only numbers
```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/11_AlignSoftmaskedRename

#genome with original scaffold names
ln -s ../02_RepeatMasker/FinalGenomePilonReducedSoftMasked.fa

#fix genome to have only numbers for names
less FinalGenomePilonReducedSoftMasked.fa |sed 's/HiC_scaffold_//g' |sed 's/_pilon//g' >FinalGenomePilonReducedSoftMaskedRecode.fa

module load hisat2
hisat2-build FinalGenomePilonReducedSoftMaskedRecode.fa FinalGenomePilonReducedSoftMaskedRecode

for f in ../../16_RNAseq/*fastq; do ln -s $f;done
vi runHISAT2.sh
#################################################################################
#!/bin/bash

module load hisat2
module load samtools
DBDIR="/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/11_AlignSoftmaskedRename"
GENOME="FinalGenomePilonReducedSoftMaskedRecode"

p=40
R1_FQ="$1"
R2_FQ="$2"


hisat2 \
  -p ${p} \
  -x ${DBDIR}/${GENOME} \
  -1 ${R1_FQ} \
  -2 ${R2_FQ}  \
  --rna-strandness RF \
  -S  ${R1_FQ}.sam &> ${R1_FQ}.log
samtools view --threads 40 -b -o ${R1_FQ}.bam ${R1_FQ}.sam
mkdir ${R1_FQ}_temp
samtools sort -m 3G -o ${R1_FQ}_sorted.bam -T ${R1_FQ}_temp --threads 40 ${R1_FQ}.bam
samtools index ${R1_FQ}_sorted.bam
samtools view -h -F 16 ${R1_FQ}.sam |samtools view -b -o Forward${R1_FQ}.bam
samtools view -h -f 16 ${R1_FQ}.sam | samtools view -b -o Reverse${R1_FQ}.bam
samtools sort -m 3G -n -o Forward${R1_FQ}_sorted.bam -T ${R1_FQ}_temp --threads 40  Forward${R1_FQ}.bam
samtools sort -m 3G -n -o Reverse${R1_FQ}_sorted.bam -T ${R1_FQ}_temp --threads 40  Reverse${R1_FQ}.bam
samtools index Forward${R1_FQ}_sorted.bam
samtools index Reverse${R1_FQ}_sorted.bam

##################################################################################

paste <(ls -1 *R1* ) <(ls -1 *R2*) |awk '{print "sh runHISAT2.sh "$0 }' >hisat2.sh


# Merge the RNAseq alignments

samtools merge -@ 20 AllStrandedRNASeq Elk-kidney_S25_L003_R1_001.fastq_sorted.bam Elk-kidney_S25_L004_R1_001.fastq_sorted.bam Elk-lung_S26_L003_R1_001.fastq_sorted.bam Elk-lung_S26_L004_R1_001.fastq_sorted.bam Elk-Mes-LN_S24_L003_R1_001.fastq_sorted.bam Elk-Mes-LN_S24_L004_R1_001.fastq_sorted.bam Elk-muscle_S21_L003_R1_001.fastq_sorted.bam Elk-muscle_S21_L004_R1_001.fastq_sorted.bam ElkpscapLN_S22_L003_R1_001.fastq_sorted.bam ElkpscapLN_S22_L004_R1_001.fastq_sorted.bam Elk-spleen_S23_L003_R1_001.fastq_sorted.bam Elk-spleen_S23_L004_R1_001.fastq_sorted.bam
samtools index AllStrandedRNASeq.bam

```

### Run Braker
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/09_BrakerRunMasked


#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH -t 168:00:00
#SBATCH --partition=mem
#SBATCH -J braker_0
#SBATCH -o braker_0.o%j
#SBATCH -e braker_0.e%j
#SBATCH --mail-user=remkv6@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module purge;module load miniconda;source activate Braker;braker.pl --species=CervusCanadensis --nocleanup --gff3  --cores 40 --genome=FinalGenomePilonReducedSoftMaskedRenamed.fa --bam=Elk-kidney_S25_L003_R1_001.fastq_sorted.bam,Elk-kidney_S25_L004_R1_001.fastq_sorted.bam,Elk-lung_S26_L003_R1_001.fastq_sorted.bam,Elk-lung_S26_L004_R1_001.fastq_sorted.bam,Elk-Mes-LN_S24_L003_R1_001.fastq_sorted.bam,Elk-Mes-LN_S24_L004_R1_001.fastq_sorted.bam,Elk-muscle_S21_L003_R1_001.fastq_sorted.bam,Elk-muscle_S21_L004_R1_001.fastq_sorted.bam,ElkpscapLN_S22_L003_R1_001.fastq_sorted.bam,ElkpscapLN_S22_L004_R1_001.fastq_sorted.bam,Elk-spleen_S23_L003_R1_001.fastq_sorted.bam,Elk-spleen_S23_L004_R1_001.fastq_sorted.bam,Undetermined_S0_L003_R1_001.fastq_sorted.bam --AUGUSTUS_CONFIG_PATH=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/05_BrakerRun/config/ --AUGUSTUS_BIN_PATH=/home/rick.masonbrink/.conda/envs/Braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/05_BrakerRun/Augustus/scripts/ --GENEMARK_PATH=/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/09_BrakerRunMasked/gm_et_linux_64/
scontrol show job $SLURM_JOB_ID


```




### Portcullis for gene prediciton
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/18_PortCullis
#uses lots of memory, so beware

ln -s ../17_Braker/11_AlignSoftmaskedRename/AllStrandedRNASeq.bam
ln -s ../17_Braker/11_AlignSoftmaskedRename/AllStrandedRNASeq.bam.bai
ln -s ../17_Braker/11_AlignSoftmaskedRename/FinalGenomePilonReducedSoftMaskedRecode.fa

ml miniconda;source activate portcullis; portcullis full --threads 36 --verbose --output portcullis_out --orientation RF --strandedness firststrand FinalGenomePilonReducedSoftMaskedRecode.fa AllStrandedRNASeq.bam
```

### Stringtie transcript assembly for gene prediction
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/19_Stringtie

ln -s ../17_Braker/11_AlignSoftmaskedRename/FinalGenomePilonReducedSoftMaskedRecode.fa
ln -s ../17_Braker/11_AlignSoftmaskedRename/AllStrandedRNASeq.bam
ln -s ../17_Braker/11_AlignSoftmaskedRename/AllStrandedRNASeq.bam.bai

ml stringtie; stringtie AllStrandedRNASeq.bam -j 5 --rf -p 40 -v -o FinalGenomePilonReducedSoftMaskedRecode.gtf

```

### run class
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/20_Class2

ln -s ../17_Braker/11_AlignSoftmaskedRename/AllStrandedRNASeq.bam
ln -s ../17_Braker/11_AlignSoftmaskedRename/AllStrandedRNASeq.bam.bai

ml miniconda; source activate Class2; run_class.pl -a AllStrandedRNASeq.bam -o AllStrandedRNASeqClass2.gtf -p 40 --verbose
```

### Trinity assembly for gene prediction
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/21_Trinity

ln -s ../17_Braker/11_AlignSoftmaskedRename/AllStrandedRNASeq.bam
ln -s ../17_Braker/11_AlignSoftmaskedRename/AllStrandedRNASeq.bam.bai
module load python_3/3.6.6
module load trinityrnaseq/2.8.4;module load salmon/0.10.1; sh runTrinity.sh AllStrandedRNASeq.bam



#runTrinity.sh
################################################################################
#!/bin/bash

module load jellyfish2
module load samtools
module load salmon
module load bowtie2
module load java_8_sdk/1.8.0_121
module load python_2/2.7.14
module load trinityrnaseq/2.8.4


bam="$1"
out=$(basename ${bam%.*} |cut -f 1 -d "_")
Trinity \
   --genome_guided_bam ${bam} \
   --SS_lib_type RF \
   --max_memory 120G \
   --genome_guided_max_intron 30000 \
   --full_cleanup \
--CPU 40

################################################################################


##  A couple of processes died due to memory shortage,  ran separately
module load jellyfish2
module load samtools
module load salmon
module load bowtie2
module load java_8_sdk/1.8.0_121
module load python_2/2.7.14
module load trinityrnaseq/2.8.4


/software/7/apps/trinityrnaseq/2.8.4/util/support_scripts/../../Trinity --single "Dir_AllStrandedRNASeq.bam.+.sam.minC1.gff/13/74/55883420_56479821.trinity.reads" --output "Dir_AllStrandedRNASeq.bam.+.sam.minC1.gff/13/74/55883420_56479821.trinity.reads.out" --CPU 10 --max_memory 120G --SS_lib_type F --seqType fa --trinity_complete --full_cleanup
/software/7/apps/trinityrnaseq/2.8.4/util/support_scripts/../../Trinity --single "Dir_AllStrandedRNASeq.bam.+.sam.minC1.gff/14/31/24198809_24586881.trinity.reads" --output "Dir_AllStrandedRNASeq.bam.+.sam.minC1.gff/14/31/24198809_24586881.trinity.reads.out" --CPU 10 --max_memory 120G --SS_lib_type F --seqType fa --trinity_complete --full_cleanup



#map trinity transcripts back to genome
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/21_Trinity/trinity_out_dir

module load gmap_gsnap/2017-03-17
sh runGmap.sh FinalGenomePilonReducedSoftMaskedRecode /home/rick.masonbrink/elk_bison_genomics/Masonbrink/21_Trinity/OldTrinityRun/ ../../23_strawberry/FinalGenomePilonReducedSoftMaskedRecode.fa Trinity-GG.fasta

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
gmap -D $dbloc -d $dbname -B 5 -t 16  --input-buffer-size=1000000 --output-buffer-size=1000000 -f gff3_match_cdna  $query >${dbname%.*}.${query%.*}.gff3
################################################################################



```

### run strawberry
```
#/project/elk_bison_genomics/Masonbrink/23_strawberry/strawberry

ln -s ../17_Braker/11_AlignSoftmaskedRename/AllStrandedRNASeq.bam
ln -s ../17_Braker/11_AlignSoftmaskedRename/AllStrandedRNASeq.bam.bai
ln -s ../17_Braker/11_AlignSoftmaskedRename/FinalGenomePilonReducedSoftMaskedRecode.fa

ml miniconda; source activate strawberry; strawberry  -T strawberry.log  --rf --no-quant -p 40 -v AllStrandedRNASeq.bam
```

### Run Mikado
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/24_mikado

#Transcriptome Assemblies
ln -s ../23_strawberry/strawberry_assembled.gtf
ln -s ../20_Class2/AllStrandedRNASeqClass2.gtf
ln -s ../18_PortCullis/portcullis_out/2-junc/portcullis_all.junctions.bed
ln -s ../17_Braker/11_AlignSoftmaskedRename/braker/augustus.hints.gff3
ln -s  ../21_Trinity/trinity_out_dir/FinalGenomePilonReducedSoftMaskedRecode.Trinity-GG.gff3

#for mikado blasting
ln -s ../17_Braker/uniprot-cervus.fasta

#genome
ln -s ../23_strawberry/FinalGenomePilonReducedSoftMaskedRecode.fa


# config script for mikado
################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH -p medium
#SBATCH -t 168:00:00
#SBATCH -J mikado_0
#SBATCH -o mikado_0.o%j
#SBATCH -e mikado_0.e%j
#SBATCH --mail-user=remkv6@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
cd /home/rick.masonbrink/elk_bison_genomics/Masonbrink/24_mikado
ml miniconda/3.6

source activate mikado2


#!/bin/bash
#setup variables
genome="FinalGenomePilonReducedSoftMaskedRecode.fa"
bam="AllStrandedRNASeq.bam"
list="list.txt"
#run splice junction prediction
junctions="portcullis_all.junctions.bed"
#configure
#mikado configure \
   --list $list \
   --reference $genome \
   --mode permissive \
   --scoring mammalian.yaml \
   --junctions $junctions \
     configuration.yaml
#prepare
mikado prepare \
   --json-conf configuration.yaml
#blast db
makeblastdb \
   -in uniprot-cervus.fasta \
   -dbtype prot \
   -parse_seqids
#blast
blastx \
  -max_target_seqs 5 \
   -num_threads 40 \
   -query mikado_prepared.fasta \
   -outfmt 5 \
   -db uniprot-cervus.fasta \
   -evalue 0.000001 2> blast.log | sed '/^$/d' > mikado.blast.xml
blastxml=mikado.blast.xml
#transdecoder
TransDecoder.LongOrfs \
   -t mikado_prepared.fasta
TransDecoder.Predict \
   -t mikado_prepared.fasta \
   --cpu 40
orfs=$(find $(pwd) -name "mikado_prepared.fasta.transdecoder.bed")
#serialise

mikado serialise \
   --start-method spawn \
   --procs 40 \
   --blast_targets uniprot-cervus.fasta \
   --json-conf configuration.yaml \
   --xml ${blastxml} \
   --orfs ${orfs}
#pick
mikado pick \
   --start-method spawn \
   --procs 40 \
   --json-conf configuration.yaml \
   --subloci_out mikado.subloci.gff3

#serialize and pick had to be run a earlier version of mikado, as the current had issues with installation
###################################################################################
```
