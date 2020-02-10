# Identify regions of X and Y homology between Cervus elaphus and our Elk

### Need blast to blast to related elk genome for to have synteny help inform Juicebox assembly.

```
#/home/rick.masonbrink/elk_bison_genomics/Masonbrink/12_C.elaphusHippelaphus2

ln -s ../09_CervusElaphusHippelaphus/^CA_002197005.1_CerEla1.0_genomic.fna
ln -s ../04_JuicerElk/01_3DNA/3d-dna/01_PostJBAssembly/FinalAssembly/3d-dna/FinalAssemblyFastaWithY.fasta

module load miniconda
source activate bioawk
bioawk -c fastx '{print $name,length($seq)}' FinalAssemblyFastaWithY.fasta |sort -k2,2nr >ScaffLengths.txt



makeblastdb -in FinalAssemblyFastaWithY.fasta -dbtype nucl -out FinalAssemblyFastaWithY.fasta

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -zxvf taxdb.tar.gz


fasta-splitter.pl --n-parts 36  GCA_002197005.1_CerEla1.0_genomic.fna
for f in *part*; do echo "sh runMegaBlast.sh "$f;done >blast.sh

#runMegaBlast.sh
###############################################################################
#!/bin/bash
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
#tar -zxvf taxdb.tar.gz

module load blast+
FASTA="$1"
blastn \
-task megablast \
-query ${FASTA} \
-db FinalAssemblyFastaWithY.fasta \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 40 \
-evalue 1e-5 \
-out ${FASTA%.**}.vs.nt.cul5.1e5.megablast.out
#################################################################################


for f in *sub;do echo "sbatch "$f;done |sed '0~15 s/$/\nsleep 15s/g' >submit.sh
cat *megablast.out >AllBlasts.out

```


```
less AllBlasts.blast.out |grep "chromosome Y" |awk '{if($10>$11){print $1,$11,$10,$5":"$12"-"$13}else {print $1,$10,$11,$5":"$12"-"$13}}' |sort -k1,1V -k2,2n |tr " " "\t" >Ychromosome.bed
[rick.masonbrink@sn-cn-18-2 09_CervusElaphusHippelaphus]$ less AllBlasts.blast.out |grep "chromosome X" |awk '{if($10>$11){print $1,$11,$10,$5":"$12"-"$13}else {print $1,$10,$11,$5":"$12"-"$13}}' |sort -k1,1V -k2,2n |tr " " "\t" >Xchromosome.bed

```

### Juicebox edits
```
Have best genome I could get out of juicebox including the splitting of the x and y, a bit of syntenic restructuring in areas with clear bowties but unclear resolutions. Placed an extra 20 scaffolds
Ran the new assembly back through seal as before in 01_JuicerPipeline.
```
