#  Need better datasets to filter genes with, how do cattle genes look aligned to the genome?


### Align Bos taurus CDS to Elk genome
```
#Need to use the longest models from each gene prediction,  braker gets preference


#/work/gif/remkv6/Olsen/Elk/01_cattleAlignment

wget ftp://ftp.ensembl.org/pub/release-99/fasta/bos_taurus/cds/Bos_taurus.ARS-UCD1.2.cds.all.fa.gz
gunzip Bos_taurus.ARS-UCD1.2.cds.all.fa.gz
grep ">" Bos_taurus.ARS-UCD1.2.cds.all.fa |awk '{print $1}' |grep "\.1" |wc
  17660   17660  388520


echo "sh runGmap.sh FinalGenomePilonReducedSoftMaskedRecode /work/gif/remkv6/Olsen/Elk/01_cattleAlignment/ FinalGenomePilonReducedSoftMaskedRecode.fa Bos_taurus.ARS-UCD1.2.cds.all.fa" >gmap.sh

#runGmap.sh
################################################################################
module load gmap-gsnap/2017-06-16-4oy56bt
dbname=$1
dbloc=$2
dbfasta=$3
query=$4
#gmap_build -d $dbname  -D $dbloc $dbfasta
gmap -D $dbloc -d $dbname -B 5 -t 35  --input-buffer-size=1000000 --output-buffer-size=1000000 -f gff3_gene  $query >${dbname%.*}.${query%.*}.gff3
#################################################################################

#### Quality metrics

#/work/gif/remkv6/Olsen/Elk/01_cattleAlignment
#04/13/20
#How many of the CDS's aligned?

#how many didnt align?
less gmap_0.e2868817 |grep -i "No Paths" |grep  "\.1" |wc
    187     935    7480

17,660 - 187 = 17,473 CDS aligned

sed -i 's/Parent=/Parent=BOS/g' FinalGenomePilonReducedSoftMaskedRecode.Bos_taurus.ARS-UCD1.2.cds.all.gff3
sed -i 's/ID=/ID=BOS/g' FinalGenomePilonReducedSoftMaskedRecode.Bos_taurus.ARS-UCD1.2.cds.all.gff3
sed -i 's/Name=/Name=BOS/g' FinalGenomePilonReducedSoftMaskedRecode.Bos_taurus.ARS-UCD1.2.cds.all.gff3

```

#### transcript alignment of cervus (Red deer)
```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/03_CervusTranscripts

sh runGmap.sh FinalGenomePilonReducedRenamedMasked /home/rick.masonbrink/elk_bison_genomics/Masonbrink/17_Braker/07_CervusTranscriptsMasked/ FinalGenomePilonReducedRenamedMasked.fa transcripts.fasta

runGmap.sh
################################################################################
module load gmap_gsnap/2017-03-17
dbname=$1
dbloc=$2
dbfasta=$3
query=$4
#gmap_build -d $dbname  -D $dbloc $dbfasta
gmap -D $dbloc -d $dbname -B 5 -t 40  --input-buffer-size=1000000 --output-buffer-size=1000000 -f gff3_gene  $query >${dbname%.*}.${query%.*}.gff3
#module load samtools
#samtools view --threads 40 -b -o ${dbname%.*}.${query%.*}.bam ${dbname%.*}.${query%.*}.sam
#mkdir ${dbname%.*}.${query%.*}_temp
#samtools sort  -o ${dbname%.*}.${query%.*}_sorted.bam -T ${dbname%.*}.${query%.*}_temp --threads 40 ${dbname%.*}.${query%.*}.bam
#samtools index ${dbname%.*}.${query%.*}_sorted.bam
##############################################################################

#input transcripts
grep -c ">" transcripts.fasta
    22259

less gmap_0.e946284|grep "No paths" |grep "\.1" |wc
    361    1805   17731

#how many of the 22,259 get genes in the
less bosgffmergeElkGenes.gff |awk '$3=="locus"'|grep "ENSBT" |wc
      22079  198711 4644131

#Transfer to Nova
###############
FinalGenomePilonReducedRenamedMasked.transcripts.gff3
############
#/work/gif/remkv6/Olsen/Elk/02_mergeMikadoBraker      

sed -i 's/HiCscaffold//g' FinalGenomePilonReducedRenamedMasked.transcripts.gff3
sed -i 's/Parent=/Parent=RED/g' FinalGenomePilonReducedRenamedMasked.transcripts.gff3
sed -i 's/ID=/ID=RED/g' FinalGenomePilonReducedRenamedMasked.transcripts.gff3


```
