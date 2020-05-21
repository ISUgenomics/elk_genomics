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

## Compile all annotations for gff


### NR Blast compilation
```
#/work/GIF/remkv6/Elk/26_Prots2Nr
cat *blastp.out |sort -k1,1V -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tNRBlastp:/g' >NRBlastp.tab
```
### NT Blastn annotation
```
#/work/GIF/remkv6/Elk/28_Transcrips2Nt
cat *blastn.out |sort -k1,1V -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tNTBlastn:/g' >NTBlastn.tab
```
### Uniprot Blastp
```
#/work/GIF/remkv6/Elk/27_ProtsUniprot
cat *blastp.out |sort -u -k1,1V |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tUPBlastp:/g'> UPBlastp.tab
```
### Uniprot Blastx
```
/work/GIF/remkv6/Elk/29_TransUniprot
cat *blastx.out |sort -k1,1 -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tUPBlastx:/g' >UPBlastx.tab
```

### Interproscan annotation
```
#figuring this out via semicolon separated fields, col two vs number of fields, with usable info.
cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |awk '$2=="TIGRFAM"' |awk -F";" '{print NF}' |sort|uniq|less

#this does not work out of the box, you will have to make sure how many columns each database retreives after you change semicolons to tabs.  Then make sure you are getting informative details rather than garbage info.
cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="SFLD" || $2=="MobiDBLite") {print $1"\t"$10"_"$9} else if($2=="Coils") {print $1"\t"$9} else if($2=="PRINTS"){print $1"\t"$8} else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==12) {print $1"\t"$10}else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==13) {print $1"\t"$9"_"$11} else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==14) {print $1"\t"$9"_"$11} else if($2=="Pfam" &&NF==11) {print $1"\t"$9"_"$10} else if($2=="Pfam" &&NF==12) {print $1"\t"$9"_"$10"_"$12} else if($2=="Pfam" &&NF==13) {print $1"\t"$8"_"$10"_"$12} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap")&& NF==11 ) {print $1,$9"_"$10} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap" )&& NF==12 ) {print $1,$9"_"$10"_"$12} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap")&& NF==13 ) {print $1,$8"_"$10"_"$11"_"$13} else if ($2=="CDD" && NF==13) {print $1"\t"$9"_"$10} else if ($2=="CDD" && NF==14) {print $1"\t"$9"_"$10"_"$12} else if ($2=="CDD" && NF==15) {print $1"\t"$8"_"$10"_"$11"_"$13} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") && NF==10) {print $1"\t"$9} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") && NF==11) {print $1"\t"$9"_"$11} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") &&  NF==12) {print $1"\t"$8"_"$10"_"$12} else if(($2=="PIRSF" ||$2=="Gene3D" )&& NF==10) {print $1"\t"$9} else if (($2=="PIRSF"||$2=="Gene3D" ) && NF==11) {print $1"\t"$9"_"$11} else if ($2=="CDD" && NF==11) {print $1"\t"$9"_"$10} else if ($2=="CDD" && NF==12) {print $1"\t"$9"_"$10"_"$13} else {print $0}}' |sort|uniq >Interpro.tab

#none should be longer than 2 tab delimited columns
 awk -F"\t" 'NF>2' Interpro.tab |wc
 0

 #how many interpro annotations did I start with?
 cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |wc
   8667  127093 1605177

#how many did I get tab formatted? Did I get them all?
cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="SFLD" || $2=="MobiDBLite") {print $1"\t"$10"_"$9} else if($2=="Coils") {print $1"\t"$9} else if($2=="PRINTS"){print $1"\t"$8} else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==12) {print $1"\t"$10}else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==13) {print $1"\t"$9"_"$11} else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==14) {print $1"\t"$9"_"$11} else if($2=="Pfam" &&NF==11) {print $1"\t"$9"_"$10} else if($2=="Pfam" &&NF==12) {print $1"\t"$9"_"$10"_"$12} else if($2=="Pfam" &&NF==13) {print $1"\t"$8"_"$10"_"$12} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap")&& NF==11 ) {print $1,$9"_"$10} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap" )&& NF==12 ) {print $1,$9"_"$10"_"$12} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap")&& NF==13 ) {print $1,$8"_"$10"_"$11"_"$13} else if ($2=="CDD" && NF==13) {print $1"\t"$9"_"$10} else if ($2=="CDD" && NF==14) {print $1"\t"$9"_"$10"_"$12} else if ($2=="CDD" && NF==15) {print $1"\t"$8"_"$10"_"$11"_"$13} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") && NF==10) {print $1"\t"$9} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") && NF==11) {print $1"\t"$9"_"$11} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") &&  NF==12) {print $1"\t"$8"_"$10"_"$12} else if(($2=="PIRSF" ||$2=="Gene3D" )&& NF==10) {print $1"\t"$9} else if (($2=="PIRSF"||$2=="Gene3D" ) && NF==11) {print $1"\t"$9"_"$11} else if ($2=="CDD" && NF==11) {print $1"\t"$9"_"$10} else if ($2=="CDD" && NF==12) {print $1"\t"$9"_"$10"_"$13} else {print $0}}' |wc
   8667   26420  662489

#I got them all, but what is the total after sort|uniq?
wc Interpro.tab
  5525  16809 429072 Interpro.tab

  cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") && NF==10) {print $1"\t"$9} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") && NF==11) {print $1"\t"$9"_"$11} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") &&  NF==12) {print $1"\t"$8"_"$10"_"$12} else {next}}' >PantherSmartSuperfamily.tab
  cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="CDD" && NF==11) {print $1"\t"$9"_"$10} else if($2=="CDD" && NF==12) {print $1"\t"$9"_"$10"_"$13} if ($2=="CDD" && NF==13) {print $1"\t"$8"_"$10"_"$11"_"$13} else {next}}' >CDD.tab
  cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if(($2=="PIRSF" ||$2=="Gene3D" )&& NF==10) {print $1"\t"$9} else if (($2=="PIRSF"||$2=="Gene3D" ) && NF==11) {print $1"\t"$9"_"$11} else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==12) {print $1"\t"$10}else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==13) {print $1"\t"$9"_"$11} else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==14) {print $1"\t"$9"_"$11} else {next}}' > Gene3dPirsf.tab
  cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="SFLD" || $2=="MobiDBLite") {print $1"\t"$10"_"$9} else if($2=="Coils") {print $1"\t"$9} else {next}}' >CoilsSfldMobidblite.tab
  cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="Pfam" &&NF==11) {print $1"\t"$9"_"$10} else if($2=="Pfam" &&NF==12) {print $1"\t"$9"_"$10"_"$12} else if($2=="Pfam" && NF==13) {print $1"\t"$8"_"$10"_"$12} else {next}}' >Pfam.tab
  cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if(($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap")&& NF==11 ) {print $1,$9"_"$10} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap" )&& NF==12 ) {print $1,$9"_"$10"_"$12} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap")&& NF==13 ) {print $1,$8"_"$10"_"$11"_"$13} else {next}}' >PrintsProdomPrositepatternsPrositeprofilesTigrfamHamap.tab

  cat *tab |sort|uniq >Interpro.tab1

#How many are not exact duplicate annotations?
wc Interpro.tab1
5334  16905 440155 Interpro.tab1

  

```
### Combine all annotationsxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```
#Fix the oddball annotation of transcript in $3
 awk -F"\t" '{if($3=="transcript"){print $1"\t"$2"\tmRNA\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9} else {print $0}}' OrderedSCNGenePredictions.gff3 >OddTranscriptFixOrderedSCNGenePredictions.gff3


#get all the blasts together
awk -F"\t" '{arr[$1]=arr[$1] "\t" $0}END{for(i in arr)print i,arr[i]}' *tab |cut -f 1,3,5,7,9,11  |sed 's/\t/#/1'|sed 's/\t/|/g' |sed 's/=//g'|sed 's/#/\tNote=/1' >CombineAnnot.tab1

#add in the interpro annotations
awk '{arr[$1]=arr[$1] "\t" $0}END{for(i in arr)print i,arr[i]}' *tab1 |sed 's/ ;/\t/1' |cut -f 2,3,5 |sed 's/\t/;/2' >CombeinAnnotIPRs.tab2


less ../07_NewGenes/OddTranscriptFixOrderedSCNGenePredictions.gff3 |awk '$3=="mRNA"' |sed 's/ID=/ID=\t/g' |sed 's/;/\t;/1' >OddTranscriptFixOrderedSCNGenePredictionsgff3.GREPMOD

#find the mrnas lacking annotations
awk '{print $1}' CombeinAnnotIPRs.tab2 |cat - <(cut -f 10 OddTranscriptFixOrderedSCNGenePredictionsgff3.GREPMOD) |sort|uniq -c |awk '$1==1 {print $2}' |cat - CombeinAnnotIPRs.tab2 >AllGenesWWOannot.tab3

#get the mrnas in the gff in the proper order and paste
paste <(sort -k10,10 OddTranscriptFixOrderedSCNGenePredictionsgff3.GREPMOD|sed 's/;Name=.*//g') <(sort -k1,1 AllGenesWWOannot.tab3) |awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9$10$11";"$13}' >AllGenesAnnotated.gff


#combine with exons,genes, utrs, cds, etc annotations
awk '$3!="mRNA"' ../07_NewGenes/OddTranscriptFixOrderedSCNGenePredictions.gff3 |sed 's/Name=.*//g'|grep -v "#" |cat - AllGenesAnnotated.gff >AllTypesWWOAnnotationsDisorganized.gff

#get the proper order for the gff, so it will show up in jbrowse
 perl gff3sort/gff3sort.pl --precise --chr_order natural AllTypesWWOAnnotationsDisorganized.gff > SCNgenomeFunctionalGeneAnnotations.gff3
```
