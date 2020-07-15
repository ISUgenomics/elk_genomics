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
CM008008.1,CM008009.1,CM008010.1,CM008011.1,CM008012.1,CM008013.1,CM008014.1,CM008015.1,CM008016.1,CM008017.1,CM008018.1,CM008019.1,CM008020.1,CM008021.1,CM008022.1,CM008023.1,CM008024.1,CM008025.1,CM008026.1,CM008027.1,CM008028.1,CM008029.1,CM008030.1,CM008031.1,CM008032.1,CM008033.1,CM008034.1,CM008035.1,CM008036.1,CM008037.1,CM008038.1,CM008039.1,CM008040.1,CM008041.1,CM008042.1
Chromosome_01,Chromosome_02,Chromosome_03,Chromosome_04,Chromosome_05,Chromosome_06,Chromosome_07,Chromosome_08,Chromosome_09,Chromosome_10,Chromosome_11,Chromosome_12,Chromosome_13,Chromosome_14,Chromosome_15,Chromosome_16,Chromosome_17,Chromosome_18,Chromosome_19,Chromosome_20,Chromosome_21,Chromosome_22,Chromosome_23,Chromosome_24,Chromosome_25,Chromosome_26,Chromosome_27,Chromosome_28,Chromosome_29,Chromosome_30,Chromosome_31,Chromosome_32,Chromosome_33,Chromosome_X,Chromosome_Y
##################################################################

cat CanadensisMCFormat.gff ElaphusMCFormat.gff > xyz.gff

#transferred to condo, having trouble with mcscanx installation
#/work/GIF/remkv6/Elk/33_MCSCANX

MCScanX xyz


java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o elkvsReddeer
cp control_file.ctl Circlecontrol_file.ctl

#Circlecontrol_file.ctl
#################################################################################
3000
CM008008.1,CM008009.1,CM008010.1,CM008011.1,CM008012.1,CM008013.1,CM008014.1,CM008015.1,CM008016.1,CM008017.1,CM008018.1,CM008019.1,CM008020.1,CM008021.1,CM008022.1,CM008023.1,CM008024.1,CM008025.1,CM008026.1,CM008027.1,CM008028.1,CM008029.1,CM008030.1,CM008031.1,CM008032.1,CM008033.1,CM008034.1,CM008035.1,CM008036.1,CM008037.1,CM008038.1,CM008039.1,CM008040.1,CM008041.1,CM008042.1,Chromosome_01,Chromosome_02,Chromosome_03,Chromosome_04,Chromosome_05,Chromosome_06,Chromosome_07,Chromosome_08,Chromosome_09,Chromosome_10,Chromosome_11,Chromosome_12,Chromosome_13,Chromosome_14,Chromosome_15,Chromosome_16,Chromosome_17,Chromosome_18,Chromosome_19,Chromosome_20,Chromosome_21,Chromosome_22,Chromosome_23,Chromosome_24,Chromosome_25,Chromosome_26,Chromosome_27,Chromosome_28,Chromosome_29,Chromosome_30,Chromosome_31,Chromosome_32,Chromosome_33,Chromosome_X,Chromosome_Y
#################################################################################


java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c RearrangementsCirclecontrol_file.ctl -o RearrangeCircleV1.png

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

```
