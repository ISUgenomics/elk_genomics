#  Need to prepare synteny plots for comparisons to C. elaphus and B. taurus
```
#/work/gif/remkv6/Olsen/Elk/06_Synteny
```
## get proteins and transcripts from our elk
```
#/work/gif/remkv6/Olsen/Elk/05a_RenameScaffsNGenes
sed  's/Name=;//g' fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 >NameRemovedfixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3

ml cufflinks
gffread NameRemovedfixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3  -g FinalGenomePilonReducedSoftMaskedRecode.fa -t mRNA -x  fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest_transcripts.fasta -y  fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest_proteins.fasta
```

## C. elaphus synteny
```
#/work/gif/remkv6/Olsen/Elk/06_Synteny/01_RedDeer

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/197/005/GCA_002197005.1_CerEla1.0/GCA_002197005.1_CerEla1.0_cds_from_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/197/005/GCA_002197005.1_CerEla1.0/GCA_002197005.1_CerEla1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/197/005/GCA_002197005.1_CerEla1.0GCA_002197005.1_CerEla1.0_protein.faa.gz

```
### BlastP ortholog calling by 50% length and 70% identity
```
##/work/gif/remkv6/Olsen/Elk/06_Synteny/01_RedDeer

mkdir 01_BlastpOrthology
cd 01_BlastpOrthology

#make blast database
ml gcc/7.3.0-xegsmw4; ml blast-plus/2.7.1-py2-ybam4tg;makeblastdb -in GCA_002197005.1_CerEla1.0_protein.faa -dbtype prot

#run the blast
ml gcc/7.3.0-xegsmw4; ml blast-plus/2.7.1-py2-ybam4tg;blastp -db GCA_002197005.1_CerEla1.0_protein.faa -query fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest_proteins.fasta -outfmt 6 -num_threads 4 -out Elk2RedDeer.blastout

#Get only the best blast hits
sort -u -k1,1 Elk2RedDeer.blastout >BestHitOnlyElk2RedDeer.blastout

#filter the best blast hits for greater than 50% query coverage and greater than 70% identity
awk '{print $1}' BestHitOnlyElk2RedDeer.blastout |cdbyank fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfid
encetest_proteins.fasta.cidx |bioawk -c fastx '{print $name,length($seq)}' |paste BestHitOnlyElk2RedDeer.blastout  - |awk '($4/$14)>.5 && $3>70' |awk '{print $1,$2}'  >PairwiseOrthology.list

#leave it in mcscanx format
awk '{print $1}' BestHitOnlyElk2RedDeer.blastout |cdbyank fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest_proteins.fasta.cidx |bioawk -c fastx '{print $name,length($seq)}' |paste BestHitOnlyElk2RedDeer.blastout  - |awk -F"\t" '($4/$14)>.5 && $3>70'  |cut -f 1-12 >xyz.blast

#elaphus genome has proteins and mRNA's with different names :/ rename
singularity shell "/opt/rit/singularity/images/maker/2.31.10_3.1/maker.simg"

awk '$3=="CDS"' GCA_002197005.1_CerEla1.0_genomic.gff |cut -f 9 |awk -F";" '{print $1,$2}' |sed 's/ID=cds-//g' |sed 's/Parent=//g' >makerIndexForProteinNames.map
map_data_ids -col 2 makerIndexForProteinNames.map 01_BlastpOrthology/xyz.blast

```

### mcscanX of red deer vs elk
```
##/work/gif/remkv6/Olsen/Elk/06_Synteny/01_RedDeer
mkdir 02_mCscanX
cd 02_mCscanX/

ln -s ../01_BlastpOrthology/PairwiseOrthology.list
ln -s ../GCA_002197005.1_CerEla1.0_protein.faa
ln -s ../01_BlastpOrthology/fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest_proteins.fasta
ln -s ../GCA_002197005.1_CerEla1.0_genomic.fna
ln -s ../../../05a_RenameScaffsNGenes/FinalGenomePilonReducedSoftMaskedRecode.fa
cp ../../../05a_RenameScaffsNGenes/fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3.gz .
gunzip fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3.gz


## Create the gff files for mcscanx

awk '$3=="mRNA" {print $1,$4,$5,$9}' fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print $1,$4,$2,$3}' |tr " " "\t" >CanadensisMCFormat.gff
awk '$3=="mRNA" {print $1,$4,$5,$9}' GCA_002197005.1_CerEla1.0_genomic.gff |sed 's/ID=//g' |sed 's/;/\t/g' |awk  '{print $1,$4,$2,$3}' |tr " " "\t" >ElaphusMCFormat.gff

#get only the pseudomolecules for the control file

cat <( bioawk -c fastx '{print $name,length($seq)}' GCA_002197005.1_CerEla1.0_genomic.fna |awk '$2>500000 {print $1}' |sort|uniq|tr "\n" "," ) <( bioawk -c fastx '{print $name,length($seq)}' FinalGenomePilonReducedSoftMaskedRecode.fa |awk '$2>500000 {print $1}' |sort|uniq|tr "\n" "," ) |less

control_file.ctl
##################################################################
3000
3000
Ce_Chr_1,Ce_Chr_2,Ce_Chr_3,Ce_Chr_4,Ce_Chr_5,Ce_Chr_6,Ce_Chr_7,Ce_Chr_8,Ce_Chr_9,Ce_Chr_10,Ce_Chr_11,Ce_Chr_12,Ce_Chr_13,Ce_Chr_14,Ce_Chr_15,Ce_Chr_16,Ce_Chr_17,Ce_Chr_18,Ce_Chr_19,Ce_Chr_20,Ce_Chr_21,Ce_Chr_22,Ce_Chr_23,Ce_Chr_24,Ce_Chr_25,Ce_Chr_26,Ce_Chr_27,Ce_Chr_28,Ce_Chr_29,Ce_Chr_30,Ce_Chr_31,Ce_Chr_32,Ce_Chr_33,Ce_Chr_X,Ce_Chr_Y
Chromosome_01,Chromosome_02,Chromosome_03,Chromosome_04,Chromosome_05,Chromosome_06,Chromosome_07,Chromosome_08,Chromosome_09,Chromosome_10,Chromosome_11,Chromosome_12,Chromosome_13,Chromosome_14,Chromosome_15,Chromosome_16,Chromosome_17,Chromosome_18,Chromosome_19,Chromosome_20,Chromosome_21,Chromosome_22,Chromosome_23,Chromosome_24,Chromosome_25,Chromosome_26,Chromosome_27,Chromosome_28,Chromosome_29,Chromosome_30,Chromosome_31,Chromosome_32,Chromosome_33,Chromosome_X,Chromosome_Y
##################################################################

cat CanadensisMCFormat.gff ElaphusMCFormat.gff > xyz.gff

#transferred to condo, having trouble with mcscanx installation
#/work/GIF/remkv6/Elk/33_MCSCANX

#produce the xyz.collinearity file.  
MCScanX xyz


#plot all chromosomes vs all chromosomes
java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o elkvsReddeer
cp control_file.ctl Circlecontrol_file.ctl

#Circlecontrol_file.ctl
#################################################################################
3000
CM008008.1,CM008009.1,CM008010.1,CM008011.1,CM008012.1,CM008013.1,CM008014.1,CM008015.1,CM008016.1,CM008017.1,CM008018.1,CM008019.1,CM008020.1,CM008021.1,CM008022.1,CM008023.1,CM008024.1,CM008025.1,CM008026.1,CM008027.1,CM008028.1,CM008029.1,CM008030.1,CM008031.1,CM008032.1,CM008033.1,CM008034.1,CM008035.1,CM008036.1,CM008037.1,CM008038.1,CM008039.1,CM008040.1,CM008041.1,CM008042.1,Chromosome_01,Chromosome_02,Chromosome_03,Chromosome_04,Chromosome_05,Chromosome_06,Chromosome_07,Chromosome_08,Chromosome_09,Chromosome_10,Chromosome_11,Chromosome_12,Chromosome_13,Chromosome_14,Chromosome_15,Chromosome_16,Chromosome_17,Chromosome_18,Chromosome_19,Chromosome_20,Chromosome_21,Chromosome_22,Chromosome_23,Chromosome_24,Chromosome_25,Chromosome_26,Chromosome_27,Chromosome_28,Chromosome_29,Chromosome_30,Chromosome_31,Chromosome_32,Chromosome_33,Chromosome_X,Chromosome_Y
#################################################################################


#made some manual modifcations to the circlecontrol_file.ctl  to get just chromosomes with large rearrangments.  

java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c RearrangementsCirclecontrol_file.ctl -o RearrangeCircleV1.png




#create a plot of only X and Y chromosomes

vi ChromosomeX.ctl
###################################################################
3000
CM008041.1,CM008042.1,Chromosome_X,Chromosome_Y
################################################################


java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c ChromosomeX.ctl -o ChromosomeXCircleV1.png
```

## B. taurus synteny
```
/work/gif/remkv6/Olsen/Elk/06_Synteny/02_BosTaurus
#downloaded proteins manually and uploaded via globus to 01_BlastpOrthology
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz

ml gcc/7.3.0-xegsmw4
ml blast-plus/2.7.1-py2-ybam4tg
ml gcc/7.3.0-xegsmw4; ml blast-plus/2.7.1-py2-ybam4tg;blastp -db GCA_002197005.1_CerEla1.0_protein.faa -query fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest_proteins.fasta -outfmt 6 -num_threads 4 -out Elk2RedDeer.blastout

mkdir 01_BlastpOrthology
cd 01_BlastpOrthology


makeblastdb -in GCF_002263795.1_ARS-UCD1.2_protein.faa -dbtype prot -out GCF_002263795.1_ARS-UCD1.2_protein.blastdb
echo "ml gcc/7.3.0-xegsmw4; ml blast-plus/2.7.1-py2-ybam4tg; blastp -db GCF_002263795.1_ARS-UCD1.2_protein.blastdb -query fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest_proteins.fasta -outfmt 6 -num_threads 36 -out Canadensis2Elaphus.blastout" >blastp.sh


#get just the best blast hit, (yes correcting the naming too)
sort -u -k1,1 Canadensis2Elaphus.blastout >BestHitOnlyElk2Taurus


#filter the best blast hits for greater than 50% query coverage and greater than 70% identity
awk '{print $1}' BestHitOnlyElk2RedDeer.blastout |cdbyank fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest_proteins.fasta.cidx |bioawk -c fastx '{print $name,length($seq)}' |paste BestHitOnlyElk2RedDeer.blastout  - |awk '($4/$14)>.5 && $3>70' |awk '{print $1,$2}'  >PairwiseOrthology.list

#get blast in correct format for mcscanx
awk '{print $1}' BestHitOnlyElk2Taurus |cdbyank fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest_proteins.fasta.cidx |bioawk -c fastx '{print $name,length($seq)}' |paste BestHitOnlyElk2Taurus  - |awk -F"\t" '($4/$14)>.5 && $3>70'  |cut -f 1-12 >xyz.blast

ln -s ../../01_RedDeer/02_mCscanX/CanadensisMCFormat.gff
awk '$3=="mRNA" {print $1,$4,$5,$9}' ../GCF_002263795.1_ARS-UCD1.2_genomic.gff |sed 's/ID=rna-//g' |sed 's/;/\t/g' |awk  '{print $1,$4,$2,$3}' |tr " " "\t" >TaurusMCFormat.gff
 cat CanadensisMCFormat.gff TaurusMCFormat.gff >xyz.gff


#make the control file

#softlink genomes
ln -s ../../../05a_RenameScaffsNGenes/FinalGenomePilonReducedSoftMaskedRecode.fa
ln -s ../GCF_002263795.1_ARS-UCD1.2_genomic.fna


cat <( bioawk -c fastx '{print $name,length($seq)}' GCF_002263795.1_ARS-UCD1.2_genomic.fna |awk '$2>5000000 {print $1}' |sort|uniq|tr "\n" "," ) <( bioawk -c fastx '{print $name,length($seq)}' FinalGenomePilonReducedSoftMaskedRecode.fa |awk '$2>5000000 {print $1}' |sort|uniq|tr "\n" "," ) |less

vi control_file.ctl
##################################################################
3000
3000
NC_037328.1,NC_037329.1,NC_037330.1,NC_037331.1,NC_037332.1,NC_037333.1,NC_037334.1,NC_037335.1,NC_037336.1,NC_037337.1,NC_037338.1,NC_037339.1,NC_037340.1,NC_037341.1,NC_037342.1,NC_037343.1,NC_037344.1,NC_037345.1,NC_037346.1,NC_037347.1,NC_037348.1,NC_037349.1,NC_037350.1,NC_037351.1,NC_037352.1,NC_037353.1,NC_037354.1,NC_037355.1,NC_037356.1,NC_037357.1,NW_020192292.1
Chromosome_01,Chromosome_02,Chromosome_03,Chromosome_04,Chromosome_05,Chromosome_06,Chromosome_07,Chromosome_08,Chromosome_09,Chromosome_10,Chromosome_11,Chromosome_12,Chromosome_13,Chromosome_14,Chromosome_15,Chromosome_16,Chromosome_17,Chromosome_18,Chromosome_19,Chromosome_20,Chromosome_21,Chromosome_22,Chromosome_23,Chromosome_24,Chromosome_25,Chromosome_26,Chromosome_27,Chromosome_28,Chromosome_29,Chromosome_30,Chromosome_31,Chromosome_32,Chromosome_33,Chromosome_X,Chromosome_Y
##################################################################

#obtain xyz collinearity file
ml mcscanx/20170403
MCScanX xyz


#plot all chromosomes vs all chromosomes
java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o elkvscow.png


cp control_file.ctl Circlecontrol_file.ctl

vi Circlecontrol_file.ctl
###############################################################################
3000
NC_037328.1,NC_037329.1,NC_037330.1,NC_037331.1,NC_037332.1,NC_037333.1,NC_037334.1,NC_037335.1,NC_037336.1,NC_037337.1,NC_037338.1,NC_037339.1,NC_037340.1,NC_037341.1,NC_037342.1,NC_037343.1,NC_037344.1,NC_037345.1,NC_037346.1,NC_037347.1,NC_037348.1,NC_037349.1,NC_037350.1,NC_037351.1,NC_037352.1,NC_037353.1,NC_037354.1,NC_037355.1,NC_037356.1,NC_037357.1,NW_020192292.1,Chromosome_01,Chromosome_02,Chromosome_03,Chromosome_04,Chromosome_05,Chromosome_06,Chromosome_07,Chromosome_08,Chromosome_09,Chromosome_10,Chromosome_11,Chromosome_12,Chromosome_13,Chromosome_14,Chromosome_15,Chromosome_16,Chromosome_17,Chromosome_18,Chromosome_19,Chromosome_20,Chromosome_21,Chromosome_22,Chromosome_23,Chromosome_24,Chromosome_25,Chromosome_26,Chromosome_27,Chromosome_28,Chromosome_29,Chromosome_30,Chromosome_31,Chromosome_32,Chromosome_33,Chromosome_X,Chromosome_Y
####################################################################################

java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c Circlecontrol_file.ctl -o CircleV1.png





#Found that Bos taurus have different names for the mRNA and peptides, so addressing
awk '$3=="CDS" {print $9}' GCF_002263795.1_ARS-UCD1.2_genomic.gff |sed 's/;/\t/g' |awk '{print $1,$2}' |sed 's/ID=cds-//g' |sed 's/Parent=rna-//g' |sort |uniq >makerIndexForProteinNames.map


singularity shell "/opt/rit/singularity/images/maker/2.31.10_3.1/maker.simg"

map_data_ids -col 2 makerIndexForProteinNames.map xyz.blast



## Create rearrangments/fission/fusion only circle plot

vi RearrangedCircolecontrol_file.ctl
###################################################################
3000
NC_037328.1,NC_037329.1,NC_037332.1,NC_037335.1,NC_037336.1,NC_037344.1,NC_037346.1,NC_037350.1,NC_037352.1,Chromosome_03,Chromosome_05,Chromosome_06,Chromosome_08,Chromosome_15,Chromosome_17,Chromosome_19,Chromosome_22,Chromosome_26,Chromosome_28,Chromosome_31,Chromosome_33
####################################################################

java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c RearrangedCirclecontrol_file.ctl -o RearrangedCircleV1.png


vi ChromosomeX.ctl
#################################################################
3000
NC_037357.1,NW_020192292.1,Chromosome_X,Chromosome_Y
####################################################################

java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c ChromosomeX.ctl -o ChromosomeXCircleV1.png

```


# iadhore

MCScanX did not provide high enough quality figures for publication, so performing synteny with iadhore
### C.elaphus iadhore
```
#/work/GIF/remkv6/Elk/34_iadhoreSynteny

#genomes, gffs, proteins, and pairwise orthology files linked
ln -s ../../33_MCSCANX/fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3
ln -s ../../33_MCSCANX/fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest_proteins.fasta
ln -s ../../33_MCSCANX/FinalGenomePilonReducedSoftMaskedRecode.fa
ln -s ../../33_MCSCANX/GCA_002197005.1_CerEla1.0_genomic.fna
ln -s ../../33_MCSCANX/GCA_002197005.1_CerEla1.0_genomic.gff
ln -s ../../33_MCSCANX/GCA_002197005.1_CerEla1.0_protein.faa
#copied PairwiseOrthology.list from nova

#change names of proteins to transcript names
ml maker/2.31.10_3.1
 singularity shell "/opt/rit/singularity/images/maker/2.31.10_3.1/maker.simg"
 map_data_ids -col 2 makerIndexForProteinNames.map PairwiseOrthology.list

 sed -i 's/ /\t/g' *list

#set up query files
mkdir query
cd query/

awk '$3=="mRNA"' ../fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print $9$7,$1}' |grep "R1" |awk '{print >> $2 ".lst"; close ($2)}'
ls *lst >input.txt
sed -i 's/ .*//g' *.lst
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini


#set up subject files
cd ../
mkdir subject
cd subject/
awk 'substr($1,1,2)=="CM"' ../../../33_MCSCANX/xyz.gff |awk '{if($4>$3) {print $2"+",$1} else {print $2"-",$1}}' |grep "RA" |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1-2 -d"." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
sed -i 's/rna-gnl|WGS:MKHE|//g' *lst


cd ../
cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini

ml GIF2/iAdHoRe/3.0.01
i-adhore iadhore.ini
#301 multiplicons
```


### Bos taurus Synteny with iadhore
```
/work/GIF/remkv6/Elk/34_iadhoreSynteny/02_Bostaurus

ln -s ../../33_MCSCANX/01_Btaurus/xyz.gff
ln -s ../../33_MCSCANX/01_Btaurus/xyz.blast
cp -rf ../01_Celaphus/query/ .

mkdir subject
cd subject/

awk 'substr($1,1,2)!="Ch" && substr($1,1,2)!="Sc"'  ../xyz.gff|awk '{if($4>$3) {print $2"+",$1} else {print $2"-",$1}}' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1-2 -d"." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini


less xyz.blast |sort -u -k1,1 |cut -f 1,2 >TopHitPairwiseOrthology.list

# add this to the iadhore file
genome=Ccanadensis

#scaffolds and locations
genome=Celaphus
#scaffolds and locations

blast_table=TopHitPairwiseOrthology.list
prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=false
output_path= output
alignment_method=gg2
gap_size=5
cluster_gap=15
level_2_only=true
q_value=.01


i-adhore iadhore.ini

389 multiplicons
```

### Create circos plot for c. elaphus vs c. canadensis
```
/work/GIF/remkv6/Elk/34_iadhoreSynteny/01_Celaphus/01_Circos

#get genoems and gffs
ln -s ../../../33_MCSCANX/FinalGenomePilonReducedSoftMaskedRecode.fa
ln -s ../../../33_MCSCANX/GCA_002197005.1_CerEla1.0_genomic.fna
ln -s ../../../33_MCSCANX/GCA_002197005.1_CerEla1.0_genomic.gff
ln -s ../../../33_MCSCANX/fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3

#create grepmod files for whole word matches on mRNA name.
 awk '$3=="mRNA"' GCA_002197005.1_CerEla1.0_genomic.gff |sed 's/ID=rna-gnl|WGS:MKHE|//g' | sed 's/\;.*//g' >Celaphus.Grepmod

 awk '$3=="mRNA"' fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |sed 's/ID=//g' |sed 's/\;.*//g' >Ccanadensis.Grepmod



#make the syntenic ribbons

 ln -s ../output/segments.txt


 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Ccanadensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line Ccanadensis.Grepmod; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Ccanadensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line Ccanadensis.Grepmod; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Ccanadensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -m 1 -w $line Celaphus.Grepmod; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Ccanadensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -m 1 -w $line Celaphus.Grepmod; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8


less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '$5!=$1' |awk '{if($5=="Ccanadensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' | paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf

less SyntenicRibbons.conf |awk '{print $3-$2}' |summary.sh


# how much gene based synteny in C. canadensis
[remkv6@condo041 01_Circos]$ less SyntenicRibbons.conf |awk '{print $3-$2}' |summary.sh
Total:  417,249,141
Count:  301
Mean:   1,386,209
Median: 706,604
Min:    11,960
Max:    10,098,677

#how much gene based synteny in C. elaphus
[remkv6@condo041 01_Circos]$ less SyntenicRibbons.conf |awk '{print $6-$5}' |summary.sh
Total:  555,237,827
Count:  301
Mean:   1,844,643
Median: 887,828
Min:    11,451
Max:    13,546,564


#Create karyotypes


bioawk -c fastx '{print $name,length($seq)}' GCA_002197005.1_CerEla1.0_genomic.fna |awk '{print "chr","-",$1,$1,"0",$2,"blue"}'  >CelaphusKaryotype.conf

bioawk -c fastx '{print $name,length($seq)}' FinalGenomePilonReducedSoftMaskedRecode.fa |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >CcanadensisKaryotype.conf

awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' CelaphusKaryotype.conf >>tmpKaryotype.conf1";done >CelaphusKaryotype.sh
 sh CelaphusKaryotype.sh

awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' CcanadensisKaryotype.conf >>tmpKaryotype.conf2";done >CcanadensisKaryotype.sh
sh CcanadensisKaryotype.sh
cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf


sort tmpKaryotype.conf2 |uniq|awk '{print $3}' |tr "\n" "," |sed 's/.$//' |awk '{print "circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - "$0" -static_rx "$0 }' |less


circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Chromosome_01,Chromosome_02,Chromosome_03,Chromosome_04,Chromosome_05,Chromosome_06,Chromosome_07,Chromosome_08,Chromosome_09,Chromosome_10,Chromosome_11,Chromosome_12,Chromosome_13,Chromosome_14,Chromosome_15,Chromosome_16,Chromosome_17,Chromosome_18,Chromosome_19,Chromosome_20,Chromosome_21,Chromosome_22,Chromosome_23,Chromosome_24,Chromosome_25,Chromosome_27,Chromosome_28,Chromosome_29,Chromosome_30,Chromosome_31,Chromosome_32,Chromosome_33,Chromosome_X -static_rx Chromosome_01,Chromosome_02,Chromosome_03,Chromosome_04,Chromosome_05,Chromosome_06,Chromosome_07,Chromosome_08,Chromosome_09,Chromosome_10,Chromosome_11,Chromosome_12,Chromosome_13,Chromosome_14,Chromosome_15,Chromosome_16,Chromosome_17,Chromosome_18,Chromosome_19,Chromosome_20,Chromosome_21,Chromosome_22,Chromosome_23,Chromosome_24,Chromosome_25,Chromosome_27,Chromosome_28,Chromosome_29,Chromosome_30,Chromosome_31,Chromosome_32,Chromosome_33,Chromosome_X



chromosomes_order = CM008022.1,Chromosome_05,CM008018.1,Chromosome_20,CM008027.1,Chromosome_04,CM008016.1,Chromosome_09,CM008011.1,Chromosome_11,CM008012.1,Chromosome_15,Chromosome_08,CM008039.1,CM008029.1,CM008021.1,Chromosome_12,CM008009.1,Chromosome_23,CM008026.1,Chromosome_19,CM008030.1,Chromosome_10,CM008041.1,Chromosome_X,CM008008.1,Chromosome_01,CM008017.1,Chromosome_02,CM008019.1,Chromosome_14,Chromosome_22,CM008034.1,CM008010.1,Chromosome_03,CM008020.1,Chromosome_18,Chromosome_07,CM008014.1,CM008024.1,CM008013.1,Chromosome_06,Chromosome_17,CM008025.1,Chromosome_13,Chromosome_24,CM008031.1,CM008036.1,CM008023.1,Chromosome_16,Chromosome_33,Chromosome_29,CM008037.1,Chromosome_30,CM008040.1,Chromosome_21,CM008028.1,Chromosome_28,CM008035.1,Chromosome_25,CM008032.1,Chromosome_31,CM008038.1,Chromosome_27,CM008015.1,Chromosome_32
```


### Circos with Bos taurus
```
/work/GIF/remkv6/Elk/34_iadhoreSynteny/02_Bostaurus/01_Circos

ln -s ../../../33_MCSCANX/01_Btaurus/GCF_002263795.1_ARS-UCD1.2_genomic.gff
ln -s ../../../33_MCSCANX/fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3
ln -s ../../../33_MCSCANX/FinalGenomePilonReducedSoftMaskedRecode.fa

#copied bos taurus genome from nova
GCF_002263795.1_ARS-UCD1.2_genomic.fna

ln -s ../output/segments.txt


#Create grepmod gff files
awk '$3=="mRNA"' GCF_002263795.1_ARS-UCD1.2_genomic.gff |sed 's/ID=rna-//g' |sed 's/;.*//g' >Btaurus.Grepmod
ln -s ../../01_Celaphus/01_Circos/Ccanadensis.Grepmod

less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Ccanadensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line Ccanadensis.Grepmod; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Ccanadensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line Ccanadensis.Grepmod; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Ccanadensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -m 1 -w $line Btaurus.Grepmod; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Ccanadensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -m 1 -w $line Btaurus.Grepmod; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8


less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '$5!=$1' |awk '{if($5=="Ccanadensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' | paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf


#How much gene based synteny in c. canadensis
awk '{print $3-$2}' SyntenicRibbons.conf |summary.sh
Total:  208,429,454
Count:  389
Mean:   535,808
Median: 217,865
Min:    8,957
Max:    6,789,973
#How much gene based synteny in b. taurus
awk '{print $6-$5}' SyntenicRibbons.conf |summary.sh
Total:  204,471,024
Count:  389
Mean:   525,632
Median: 203,727
Min:    5,134
Max:    7,186,933



#Create karyotypes


bioawk -c fastx '{print $name,length($seq)}' GCF_002263795.1_ARS-UCD1.2_genomic.fna |awk '{print "chr","-",$1,$1,"0",$2,"blue"}'  >BtaurusKaryotype.conf

bioawk -c fastx '{print $name,length($seq)}' FinalGenomePilonReducedSoftMaskedRecode.fa |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >CcanadensisKaryotype.conf

awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' BtaurusKaryotype.conf >>tmpKaryotype.conf1";done >BtaurusKaryotype.sh
 sh BtaurusKaryotype.sh

awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' CcanadensisKaryotype.conf >>tmpKaryotype.conf2";done >CcanadensisKaryotype.sh
sh CcanadensisKaryotype.sh
cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf

#organize ribbons so they dont overlap a whole bunch
orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Chromosome_01,Chromosome_02,Chromosome_03,Chromosome_04,Chromosome_05,Chromosome_06,Chromosome_07,Chromosome_08,Chromosome_09,Chromosome_10,Chromosome_11,Chromosome_12,Chromosome_13,Chromosome_14,Chromosome_15,Chromosome_16,Chromosome_17,Chromosome_18,Chromosome_19,Chromosome_20,Chromosome_21,Chromosome_22,Chromosome_23,Chromosome_24,Chromosome_25,Chromosome_26,Chromosome_27,Chromosome_28,Chromosome_29,Chromosome_30,Chromosome_31,Chromosome_32,Chromosome_33,Chromosome_X -static_rx Chromosome_01,Chromosome_02,Chromosome_03,Chromosome_04,Chromosome_05,Chromosome_06,Chromosome_07,Chromosome_08,Chromosome_09,Chromosome_10,Chromosome_11,Chromosome_12,Chromosome_13,Chromosome_14,Chromosome_15,Chromosome_16,Chromosome_17,Chromosome_18,Chromosome_19,Chromosome_20,Chromosome_21,Chromosome_22,Chromosome_23,Chromosome_24,Chromosome_25,Chromosome_26,Chromosome_27,Chromosome_28,Chromosome_29,Chromosome_30,Chromosome_31,Chromosome_32,Chromosome_33,Chromosome_X

chromosomes_order = Chromosome_05,NC_037344.1,NC_037346.1,NC_037345.1,Chromosome_04,NC_037329.1,Chromosome_08,Chromosome_33,NC_037332.1,Chromosome_22,Chromosome_03,NC_037330.1,Chromosome_20,NC_037334.1,Chromosome_09,NC_037328.1,Chromosome_19,Chromosome_31,NC_037342.1,Chromosome_01,NC_037338.1,Chromosome_11,NC_037331.1,Chromosome_18,NC_037352.1,Chromosome_10,Chromosome_15,NC_037355.1,NC_037353.1,NC_037333.1,Chromosome_06,Chromosome_17,NC_037356.1,Chromosome_02,NC_037335.1,Chromosome_29,Chromosome_16,NC_037340.1,Chromosome_23,NC_037337.1,Chromosome_12,NC_037349.1,Chromosome_24,NC_037348.1,Chromosome_13,NC_037357.1,Chromosome_X,NC_037341.1,Chromosome_21,NC_037350.1,Chromosome_07,NC_037343.1,Chromosome_14,NC_037347.1,Chromosome_25,NC_037336.1,Chromosome_28,Chromosome_26,NC_037339.1,Chromosome_30,NC_037351.1,Chromosome_27,NC_037354.1,Chromosome_32

```
