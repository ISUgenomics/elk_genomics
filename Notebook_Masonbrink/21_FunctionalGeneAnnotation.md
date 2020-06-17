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

```
