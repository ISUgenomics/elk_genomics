#  There appears to be some contamination in the genome assembly, run blobtools

```
#/work/GIF/remkv6/Elk/01_Blobtools


sh runMegablast.sh north_american_elk_15Jun2018_oY8t2.fasta
################################################################################
#!/bin/bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -zxvf taxdb.tar.gz

module load blast-plus
FASTA="$1"
blastn \
-task megablast \
-query ${FASTA} \
-db /work/GIF/GIF3/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/nt/nt \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 16 \
-evalue 1e-5 \
-out ${FASTA%.**}.vs.nt.cul5.1e5.megablast.out
################################################################################


#how many are not cattle related?
less north_american_elk_15Jun2018_oY8t2.vs.nt.cul5.1e5.megablast.out |sort -k1,1 -u |grep -v "Ovis canadensis" |grep -v "Bos taurus" |grep -v "Muntiacus reevesi" |grep -v "Muntiacus vaginalis" |grep -v "Ovis aries" |grep -v "Bison bison" |less

```
