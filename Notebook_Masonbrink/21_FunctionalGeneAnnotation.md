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


    cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") && NF==10) {print $1"\t"$9} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") && NF==11) {print $1"\t"$9","$11} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") &&  NF==12) {print $1"\t"$8","$10","$12} else {next}}' >PantherSmartSuperfamily.tab

    cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="CDD" && NF==11) {print $1"\t"$9","$10} else if($2=="CDD" && NF==12) {print $1"\t"$9","$10","$13} if ($2=="CDD" && NF==13) {print $1"\t"$8","$10","$11","$13} else {next}}' >CDD.tab

    cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if(($2=="PIRSF" ||$2=="Gene3D" )&& NF==10) {print $1"\t"$9} else if (($2=="PIRSF"||$2=="Gene3D" ) && NF==11) {print $1"\t"$9","$11} else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==12) {print $1"\t"$10}else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==13) {print $1"\t"$9","$11} else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==14) {print $1"\t"$9","$11} else {next}}' > Gene3dPirsf.tab

    cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="SFLD" || $2=="MobiDBLite") {print $1"\t"$10","$9} else if($2=="Coils") {print $1"\t"$9} else {next}}' >CoilsSfldMobidblite.tab
    cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="Pfam" &&NF==11) {print $1"\t"$9","$10} else if($2=="Pfam" &&NF==12) {print $1"\t"$9","$10","$12} else if($2=="Pfam" && NF==13) {print $1"\t"$8","$10","$12} else {next}}' >Pfam.tab

    cat *fasta.gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if(($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap")&& NF==11 ) {print $1"\t"$9","$10} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap" )&& NF==12 ) {print $1"\t"$9","$10","$12} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap")&& NF==13 ) {print $1"\t"$8","$10","$11","$13} else {next}}' >PrintsProdomPrositepatternsPrositeprofilesTigrfamHamap.tab

  cat *tab |sort|uniq >Interpro.tab1

#How many are not exact duplicate annotations?
wc Interpro.tab1
5334  16905 440155 Interpro.tab1



```
### Combine all annotations
```
#/work/GIF/remkv6/Elk/32_CombineFunctionalAnnotations

#get interproscan tab file into the righ format
cp ../25_Interpro/Interpro.tab1 .
sed -i 's/Name=//g' Interpro.tab1

cat *tab |sort|uniq|awk -F"\t" '{arr[$1]=arr[$1] "," $2}END{for(i in arr)print i,arr[i]}'   |sed 's/,/\t/1' >Blasts.tab

cat Blasts.tab Interpro.tab1 |awk -F"\t" '{arr[$1]=arr[$1] "," $2}END{for(i in arr)print i,arr[i]}' |sed 's/ ,$//g' |sed 's/  / /g' |sed 's/ /\t/1' |sed 's/\t,/\t/g' >ImperfectListAnnotations.tab

#get all the mRNAs, even those without annotations
awk '$3=="mRNA" ||$3=="transcript"' Bos_ReddeerReductionVHEJ.gff |sed 's/geneID=.*locus=/locus=/g' |sed 's/Parent=.*locus=/locus=/g'|sed 's/;/\t/g' |sed 's/ID=//g' |cut -f 9 |cat ImperfectListAnnotations.tab - |sort -k1,1 -u >AllGenesAddedAllAnnotations.tab

paste <(awk '$3=="mRNA" ||$3=="transcript"' Bos_ReddeerReductionVHEJ.gff |sed 's/ID=/ID=\t/1' |sed 's/;/\t;/1' |tr " " "\t" |sort -k10,10 |sed 's/geneID=.*locus=/Parent=/g' )  <(sort -k1,1 AllGenesAddedAllAnnotations.tab) |awk -F"\t"  '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9$10$11";Note="$13}' >AnnotatedmRNAs.gff3

cat <(awk '$3!="mRNA" && $3!="transcript"' Bos_ReddeerReductionVHEJ.gff) AnnotatedmRNAs.gff3 |sort -k1,1V -k4,5nr >AnnotatedGeneModels.gff

#get the proper order for the gff, so it will show up in jbrowse
sed 's/=/#/g' OrderedAnnotatedGeneModels.gffReorder.gff |sed 's/ID#/ID=/g' |sed 's/Parent#/Parent=/g' |sed 's/transcripts#/transcripts=/g' |sed 's/geneid#/geneid=/g' |sed 's/Note#/Note=/g' |sed 's/genes#/genes=/g' >EqualsFixed.gff

 perl gff3sort/gff3sort.pl --precise --chr_order natural EqualisFixed.gff > OrderedAnnotatedGeneModels.gff3

#How many were annotated?
less OrderedAnnotatedGeneModels.gff3 |awk '$3=="mRNA"' |sed 's/Note=./\t/' |awk -F"\t" '{print NF}' |sort|uniq -c| less

110262 10 Annotated
 17339 9 Not annotated


#reorder and confirm the ontology of the gff for mikado util grep
sed -i 's/\tgeme\t/\tgene\t/g'  OrderedAnnotatedGeneModels.gff3
gffread -F -O -E OrderedAnnotatedGeneModels.gff3 -o OrderedAnnotatedGeneModels.gffReorder.gff

#just get teh high confidence genes
mikado util grep ../02_mergeMikadoBraker/HighConfidence4MikadoGrep.list OrderedAnnotatedGeneModels.gffReorder.gff   HighConfidenceOrderedAnnotatedGeneModels.gff3
#get rid of genes annotated as TE's
mikado util grep -v AllUniqueTransposonGenesmRNAs.list  HighConfidenceOrderedAnnotatedGeneModels.gff3 NOTEHighConfidenceOrderedAnnotatedGeneModels.gff3

perl gff3sort/gff3sort.pl --precise --chr_order natural GTReOrderedAnnotatedGeneModels.gffReorder.gff > OrderedGTReOrderedAnnota
tedGeneModels.gffReorder.gff
gffread -F -O -E GTReOrderedAnnotatedGeneModels.gffReorder.gff -o FixedGTReOrderedAnnotatedGeneModels.gffReorder.gff
mikado util grep ../02_mergeMikadoBraker/HighConfidence4MikadoGrep.list FixedGTReOrderedAnnotatedGeneModels.gffReorder.gff >HighConfidencetest.gff3
gffread -F -O -E HighConfidencetest.gff3 -o GTHighConfidencetest.gff3
perl gff3sort/gff3sort.pl --precise --chr_order natural GTHighConfidencetest.gff3 > OrderedGTHighConfidencetest.gff3
sh ~/common_scripts/runTabix.sh OrderedGTHighConfidencetest.gff3
#WORKS! now remove TEs from high confidence
mikado util grep  -v AllUniqueTransposonGenesmRNAs.list HighConfidencetest.gff3 >NOTEHighConfidencetest.gff3
gffread -F -O -E NOTEHighConfidencetest.gff3 -o GTNOTEHighConfidencetest.gff3
perl gff3sort/gff3sort.pl --precise --chr_order natural GTNOTEHighConfidencetest.gff3 > OrderedGTNOTEHighConfidencetest.gff3
sh ~/common_scripts/runTabix.sh OrderedGTNOTEHighConfidencetest.gff3
#works!, this is the non transposable element, high confidence gene model gff


#create for genes that overlap repeats
RepetitiveMrnasOrderedAnnotatedGeneModels.list
mikado util grep  RepetitiveMrnasOrderedAnnotatedGeneModels.list FixedGTReOrderedAnnotatedGeneModels.gffReorder.gff >RepetitiveGenes.gff3
gffread -F -O -E RepetitiveGenes.gff3 -o GTRepetitiveGenes.gff3
perl gff3sort/gff3sort.pl --precise --chr_order natural GTRepetitiveGenes.gff3 > OrderedGTRepetitiveGenes.gff3
sh ~/common_scripts/runTabix.sh OrderedGTRepetitiveGenes.gff3


#creat the EDTA track

#modify EDTA gff to get to display in jbrowse
less AllAnno.gff |sed 's/;Method=.*//g' |sed  's/ID=/ID=\t/1' |sed 's/Parent=/ID=\t/1' |awk '{print $1"\t"$2"\tmRNA\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9NR$10";Note="$3}' |sed 's/Parent=/ID=/g' >AllAnnoMod.gff
gt gff3 -force -tidy -sortlines -o AllAnnoEDTA.gff AllAnnoMod.gff
sh ~/common_scripts/runTabix.sh AllAnnoEDTA.gff

#create repeatmodeler track



#All genes and transcripts assembled
gffread -F -O -E OrderedAnnotatedGeneModels_sorted.gff -o GTOrderedAnnotatedGeneModels_sorted.gff
perl gff3sort/gff3sort.pl --precise --chr_order natural GTOrderedAnnotatedGeneModels_sorted.gff > OrderedGTOrderedAnnotatedGeneModels_sorted.gff
sh ~/common_scripts/runTabix.sh OrderedGTOrderedAnnotatedGeneModels_sorted.gff

#Low confidence genes
mikado util grep  -v ../02_mergeMikadoBraker/HighConfidence4MikadoGrep.list FixedGTReOrderedAnnotatedGeneModels.gffReorder.gff >LOWConfidencetest.gff3
gffread -F -O -E LOWConfidencetest.gff3 -o GTLOWConfidencetest.gff3
perl gff3sort/gff3sort.pl --precise --chr_order natural GTLOWConfidencetest.gff3 > OrderedGTLOWConfidencetest.gff3
sh ~/common_scripts/runTabix.sh OrderedGTLOWConfidencetest.gff3


#gene that could possibly be added to high confidence, due to large numbers of annotations from multiple sources. Though most of these may already be represented by another better assembled transcript...
less GTLOWConfidencetest.gff3 |awk '$3=="mRNA"' |grep -v -i -e "integrase" -e "transposon" -e "helitron" -e "maverick" |grep -v -e "LINE" -e "SINE" |grep -v -i -e "hypothetical" -e "transposable" |awk -F";" 'NF>4' |grep -v "null)" |grep -v "BAC" >LowConfidenceOrderedAnnotatedGeneModelsThatCouldBeAdded2Highconfidence.gff3

wc LowConfidenceOrderedAnnotatedGeneModelsThatCouldBeAdded2Highconfidence.gff3
   3228   84232 1073751 LowConfidenceOrderedAnnotatedGeneModelsThatCouldBeAdded2Highconfidence.gff3

```



###  Files in proper format, filter for number
```



#Genes lacking Annotations
less OrderedAnnotatedGeneModels.gff3 | awk '$3=="mRNA"' |sed 's/mobidb-lite,signature_desc#consensus disorder prediction//g' |sed 's/(null)//g' |grep -v  "Note=." |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1,2 |sed 's/Parent=//g' >NoAnnotationGenesmRNAs.list
wc NoAnnotationGenesmRNAs.list
17574  35148 452703 NoAnnotationGenesmRNAs.list



#transposon genes to be removed (determined through annotation)
less OrderedAnnotatedGeneModels.gff3 |awk '$3=="mRNA"' |sed 's/mobidb-lite,signature_desc#consensus disorder prediction//g' |sed 's/(null)//g' |grep  -e "LINE" -e "SINE" |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1,2 |sed 's/Parent=//g' |sed 's/ID=//g'  >TransposonGenesmRNAs.list

#add in lines and sines with different statement
less OrderedAnnotatedGeneModels.gff3 |awk '$3=="mRNA"' |sed 's/mobidb-lite,signature_desc#consensus disorder prediction//g' |sed 's/(null)//g' |grep  -e "LINE" -e "SINE" |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1,2 |sed 's/Parent=//g' |sed 's/ID=//g'  >TransposonGenesmRNAs.list
 less OrderedAnnotatedGeneModels.gff3 |awk '$3=="mRNA"' |sed 's/mobidb-lite,signature_desc#consensus disorder prediction//g' |sed 's/(null)//g' |grep  -i -e "transpos" -e "integrase" -e "helitron" -e "maverick" |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1,2 |sed 's/Parent=//g' |sed 's/ID=//g' |cat - TransposonGenesmRNAs.list |sort|uniq >AllUniqueTransposonGenesmRNAs.list

 wc AllUniqueTransposonGenesmRNAs.list
  6735  13470 205784 AllUniqueTransposonGenesmRNAs.list





# Repeat GeneFiltering

#/work/gif/remkv6/Olsen/Elk/05_EvaluatePrediction

awk '$3=="mRNA"' OrderedAnnotatedGeneModels.gff3 |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |sed 's/Parent=//g' |cut -f 1,2 |sed  's/geneRLOC/RLOC/g' >OrderedAnnotatedGeneModels.list

cat <( bedtools intersect   -wo -a OrderedAnnotatedGeneModels.gff3 -b ../05_EvaluatePrediction/AllAnno.gff |awk '$3=="CDS"' |awk '{print $9}'|sort|uniq -c|awk '$1>1 {print $2}'|sed 's/Parent=//g') <(bedtools intersect -f .2  -wo -a gttest.gff3 -b ../05_EvaluatePrediction/AllAnno.gff |awk '$3=="mRNA"' |awk '{print $9}'|sort|uniq |sed 's/ID=//g'|sed 's/;/\t/g' |cut -f 1)|sort|uniq |awk '{print $1"\tRepeats"}' |cat OrderedAnnotatedGeneModels.list - |awk -F"\t" '{arr[$1]=arr[$1] "\t" $2}END{for(i in arr)print i,arr[i]}' |awk '$3=="Repeats"' |cut -f 1,2 >RepetitiveMrnasOrderedAnnotatedGeneModels.list

wc RepetitiveMrnasOrderedAnnotatedGeneModels.list
 23740  47480 496636 RepetitiveMrnasgttest.list


# survey the low confidence genes that have functional annotations other than hypothetical, BAC library, or transposon
less GTLOWConfidencetest.gff3 |awk '$3=="mRNA"' |grep -v -i -e "integrase" -e "transposon" -e "helitron" -e "maverick" |grep -v -e "LINE" -e "SINE" |grep -v -i -e "hypothetical" -e "transposable" |awk -F";" 'NF>4' |grep -v "null)" |grep -v "BAC" >LowConfidenceOrderedAnnotatedGeneModelsThatCouldBeAdded2Highconfidence.gff3

cat LowConfidenceOrderedAnnotatedGeneModelsThatCouldBeAdded2Highconfidence.gff3 GTNOTEHighConfidencetest.gff3 >HCLCtest.gff3


```


### Rename Scaffolds, genes and transcripts
```
#/work/gif/remkv6/Olsen/Elk/05a_RenameScaffsNGenes

#maker has issues with the gff format, so it will not rename the genes/mRNAs, need to figure out how to get the format right.

#gets rid of geneID and transcripts associations that are not always correct
sed 's/geneID=.*;//g' GTOrderedGTNOTEHighConfidencetest.gff3 |sed 's/transcripts=.*;//g' |sed 's/transcripts=.*//g
' >GeneidTranscriptsRemovedGTOrderedGTNOTEHighConfidencetest.gff3

awk '$3=="gene"' GeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 >genes.gff
less GeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |awk '$3=="mRNA"' |grep -v "Parent" >mRNAsNoParent.gff


cat <(less GeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |awk '$3=="mRNA"' |grep  "Parent") <(awk '$3!="mRNA"' GeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3)|less


#adds parent= to each mRNA lacking it.
bedtools intersect -wo -f .8 -a mRNAsNoParent.gff -b genes.gff |sort -u -k9,9|cut -f 1,2,3,4,5,6,7,8,9,18|sed 's/;/\t/1' |awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9";"$11";Name="$9";"$10}' |sed 's/ID=/Parent=/2' |sed 's/;;/;/g'  >415Fix.gff

#concatenates and sorts the gff again
cat <(less GeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |awk '$3=="mRNA"' |grep  "Parent" |sed 's/;/\t/1' |awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9";"$9";"$10}' |sed 's/ID=/Name=/2' |less) <(awk '$3!="mRNA"' GeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3) 415Fix.gff |awk -F"\t" '{if($3=="CDS") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="substr($9,8,length($9))"-cds"NR";"$9} else if ($3=="exon") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="substr($9,8,length($9))"-exon"NR";"$9} else {print $0}} ' |sort -k1,1V -k4,5n |sed 's/ mikado\./ mikado_/g'  >RevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3



paste <(cut -f 1-8 RevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 ) <(cut -f 9 RevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |sed 's/\./_/g' ) |sed 's/,/-/g' |sed 's/Name=ID=/Name=/g' |sed 's/Name=Parent=/Name=/g' |sed 's/ID=/Parent=/2'  |sort -k1,1V -k4,5n >NoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3


#run through genometools to get rid of some errors
gt gff3  -fixregionboundaries -sort -o GTNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 NoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3


gff3_QC -f ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa  -noncg -g GTNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3  -o error.txt -s statistic.txt

gff3_fix -qc_r error.txt -g GTNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 -og CorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3


grep -v  "Parent=RLOC" NoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |awk '$3=="mRNA"' >messeUpmRNAs.gff

bedtools intersect -wo -f .8 -a <(sed 's/Parent=.*;//g' messeUpmRNAs.gff |sed 's/Parent=.*//g' ) -b genes.gff |sort -u -k9,9|cut -f 1,2,3,4,5,6,7,8,9,18|sed 's/ID=/Parent=/2' |sed 's/\t//9' >fixedmesseUpmRNAs.gff

cat <(grep  "Parent=RLOC" NoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |awk '$3=="mRNA"' )  <(awk '$3!="mRNA"' NoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 ) fixedmesseUpmRNAs.gff |sort -k1,1V -k4,5n >fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3


ml maker
singularity shell "/opt/rit/singularity/images/maker/2.31.10_3.1/maker.simg"
maker_map_ids --prefix Cercan_  --iterate 1 --justify 8 fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 >test
map_gff_ids test fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3


cp ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa .
#copied from ceres and manipulated to fix spaces
sed 's/ //g' ScaffoldNames.map |sed 's/S/\tS/g' |sed 's/Chr/\tChr/g' >FixedScaffoldNames.map

map_fasta_ids FixedScaffoldNames.map FinalGenomePilonReducedSoftMaskedRecode.fa
map_data_ids FixedScaffoldNames.map fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3

samtools faidx FinalGenomePilonReducedSoftMaskedRecode.fa

```

#### Gene functions are not displaying, add note to gene features
```
#my laptop
#/home/remkv6/Documents/1ElkUSDAGenome/RenamedJBrowseFinal



awk '$3=="gene"' fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 >GenesOnly.gff3

#check to make sure we got all the genes.
awk '$3=="mRNA" {print $9}' fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |grep -v "Note" |sed 's/Parent=//g' |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 2 |sort|uniq|cat  annotatedGenes.list -|sort|uniq|wc
wc GenesOnly.gff3

awk '{print $1}'  annotatedGenes.list |grep -v - GenesOnly.gff3 >UnannotatedGenes.gff3
awk '{print $1}'  annotatedGenes.list |while read line; do grep  $line GenesOnly.gff3;done >AnnotatedGenes.gff3

paste AnnotatedGenes.gff3 annotatedGenes.list |cut -f 1-9,11|sed 's/\t/;/9' |cat  - UnannotatedGenes.gff3 |cat - <(awk '$3!="gene"' fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3) |sort -k1,1V -k4,5nr >GenePredictionsDisplayFixed.gff3

#transfer to nova to get sorted properly
#/work/gif/remkv6/Olsen/Elk/05a_RenameScaffsNGenes

```

### Statistics of functional gene predictions
```
#/home/remkv6/Documents/1ElkUSDAGenome/RenamedJBrowseFinal


```
