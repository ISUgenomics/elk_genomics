# Create functional annotations for elk gene predictions


### Interproscan
```
/work/GIF/remkv6/Elk/25_Interpro
#version 5.27-66.0


fasta-splitter.pl --n-parts 99 ../Bos_ReddeerReductionVHEJ_proteins.fasta
#remove the characters from the sequences that interproscan flags and stumbles upon
for f in *fasta; do sed -i -e '/^[^>]/s/[^GPAVLIMCFYWHKRQNEDST]//g' $f;done



for f in *part*; do mkdir $f.TEMP; done
for f in *part*.fasta; do echo "module use /work/GIF/software/modules;ml GIF/interproscan/5.27-66.0; interproscan.sh -cpu 16 -f TSV,GFF3 -goterms -i "$f" -pa -T "$f".TEMP";done >interproscan.sh
```


### Proteins to NR
```
#/work/GIF/remkv6/Elk/26_Prots2Nr
#downloaded December 8 2019

fasta-splitter.pl --n-parts 99 ../Bos_ReddeerReductionVHEJ_proteins.fasta
for f in *fasta; do sed -i -e '/^[^>]/s/[^GPAVLIMCFYWHKRQNEDST]//g' $f;done
for f in *part*; do echo "sh runBlastP2NR.sh "$f;done >blast.sh

#runBlastP2NR.sh
###########################################################################
#!/bin/bash
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
#tar -zxvf taxdb.tar.gz

module load blast-plus
FASTA="$1"
blastp \
-query ${FASTA} \
-db /work/GIF/databases/ncbi_nr/nr \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 16 \
-evalue 1e-5 \
-out ${FASTA%.**}.vs.nr.cul5.1e5.blastp.out
###########################################################################

```

### Proteins to uniprot/swissprot
```
#/work/GIF/remkv6/Elk/27_ProtsUniprot
#Downloaded December 9 2019

fasta-splitter.pl --n-parts 99 ../Bos_ReddeerReductionVHEJ_proteins.fasta
for f in *fasta; do sed -i -e '/^[^>]/s/[^GPAVLIMCFYWHKRQNEDST]//g' $f;done

for f in *part*; do echo "sh runBlastP2uniprot_swissprot.sh "$f;done >blast.sh

#runBlastP2uniprot_swissprot.sh
##########################################################################
#!/bin/bash
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
#tar -zxvf taxdb.tar.gz

module load blast-plus
FASTA="$1"
blastp \
-query ${FASTA} \
-db /work/GIF/databases/uniprot_sprot/uniprot_sprot.fasta \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 16 \
-evalue 1e-5 \
-out ${FASTA%.**}.vs.uniprot_sprot.cul5.1e5.blastp.out
##########################################################################


```

### Transcripts to NT
```
#Downloaded 10-23-19
#/work/GIF/remkv6/Elk/28_Transcrips2Nt
ln -s ~/common_scripts/runBlastN2NT.sh

 fasta-splitter.pl --n-parts 99 ../Bos_ReddeerReductionVHEJ_transcripts.fasta

for f in *part*; do echo "sh runBlastP2uniprot_swissprot.sh "$f;done >blast.sh

#runBlastP2uniprot_swissprot.sh
##########################################################################
#!/bin/bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -zxvf taxdb.tar.gz

module load blast-plus
FASTA="$1"
blastn \
-query ${FASTA} \
-db /work/GIF/databases/ncbi_nt/nt \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 16 \
-evalue 1e-5 \
-out ${FASTA%.**}.vs.nt.cul5.1e5.blastn.out
##########################################################################

```

### Transcripts to uniprot/swissprot
```
#/work/GIF/remkv6/Elk/29_TransUniprot

ln -s ../24_mikado/01_mikado2/FinalGenePrediction.transcripts.fasta
fasta-splitter.pl --n-parts 16 FinalGenePrediction.transcripts.fasta

for f in *part*; do echo "sh runBlastX2uniprot_swissprot.sh "$f;done >blast.sh

#runBlastX2uniprot_swissprot.sh
##########################################################################
#!/bin/bash
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
#tar -zxvf taxdb.tar.gz

module load blast-plus
FASTA="$1"
blastx \
-query ${FASTA} \
-db /work/GIF/databases/uniprot_sprot/uniprot_sprot.fasta \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 16 \
-evalue 1e-5 \
-out ${FASTA%.**}.vs.uniprot_sprot.cul5.1e5.blastx.out
##########################################################################

```
