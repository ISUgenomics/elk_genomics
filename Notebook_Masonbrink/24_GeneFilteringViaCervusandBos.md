#  Need to reduce gene models and identify those that are high confidence

```
#/work/gif/remkv6/Olsen/Elk/02_mergeMikadoBraker

#red deer modified transcript alignment to the masked genome
FinalGenomePilonReducedRenamedMasked.transcripts.gff3


ml cufflinks
gffread -VHEJ -K -M -Q -g ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa ../05_EvaluatePrediction/Fixedaugustus.hints.gff3 ../05_EvaluatePrediction/ExpressedIntactORFmikado.loci.gff3 -o gffmergeElkGenes.gff

#get proteins from this merge to run busco
 gffread gffmergeElkGenes.gff  -g ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa -VHEJ -t mRNA -x gffmergeElkGenesVHEJ_transcripts.fasta -y gffmergeElkGenesVHEJ_proteins.fasta


#merge the genes that -K -M -Q  with the red deer and bos taurus alignments
gffread -K -VHEJ -M -Q -d dupinfo -g../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa gffmergeElkGenes.gff FinalGenomePilonReducedSoftMaskedRecode.Bos_taurus.ARS-UCD1.2.cds.all.gff3 FinalGenomePilonReducedRenamedMasked.transcripts.gff3 -o Bos_ReddeerReductionTEST.gff

#input merged gene prediction
awk '$3=="locus"'  gffmergeElkGenes.gff|wc
  92756  834804 11089322


#total loci from merging elk to bos and red deer
awk '$3=="locus"'  Bos_ReddeerReduction.gff|wc
  96744  870696 13605644

#How many of the original merged braker/mikado set?
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -e "mRNA" -e  "mikado"  -e "MIKADO" |wc
  89762  807858 12573055
#~3k reduction

#How many loci are there that have either RED or BOS(ENS) in the locus tag?
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -e "RED" -e "ENS" |wc
  23829  214461 5592002

#How many have RED or BOS(ENS) also have a transcript in mikado or braker "mikado"  "MIKADO" and "mRNA"   
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -e "RED" -e "ENS" |grep -e "mikado" -e "mRNA" -e "MIKADO" |wc
  16847  151623 4559413

### How well did each prediction do?
  Braker was much better
#### How many BOS(ENS) or Elk(RED) loci were also in braker
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -e "RED" -e "ENS" |grep "mRNA" |wc
  14766  132894 4039738


#### How many BOS(ENS) or Elk(RED) loci were also in mikado
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -e "RED" -e "ENS" |grep -e "mikado"  -e "MIKADO" |wc
   8451   76059 2585545

# Get ready for busco4 evaluation of proteins
gffread Bos_ReddeerReduction.gff  -g ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa -VHEJ -t mRNA -x  Bos_ReddeerReductionVHEJ_transcripts.fasta -y  Bos_ReddeerReductionVHEJ_proteins.fasta



```

### Filter out all other genes that are not Elk

```
#How many loci are either deer or cattle
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -c -e "RED" -e "ENS"
23829

How many of the cattle/deer loci do not have a gene assembled from the rnaseq in elk?
less -S Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep "," |grep -v -e "mikado" -e "MIKADO" -e "mRNA" |wc
   1643   14787  338059

#Dual support that has transcriptional support
less -S Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep "," |grep -e "mRNA" -e "mikado" -e "MIKADO"|wc
     18563  167067 4783092


#How many of the braker genes were merged via mikado/red deer/cattle
bedtools intersect -wo -b <(awk '$3=="locus"' Bos_ReddeerReductionTEST.gff|grep  "," ) -a ../05_EvaluatePrediction/Fixedaugustus.hints.gff3 |awk '$3=="gene"' |grep -e "ENS" -e "RED"|cut -f 9 |sort|uniq|wc
  36736   36736  471500


#How many genes are there after collapsing via mikado/red deer/cattle?
bedtools intersect -wo -b <(awk '$3=="locus"' Bos_ReddeerReductionTEST.gff|grep  "," ) -a ../05_EvaluatePrediction/Fixedaugustus.hints.gff3 |awk '$3=="gene"' |grep -e "ENS" -e "RED"|cut -f 10-18 |sort|uniq|wc
 15585  140265 4195710


#How many braker genes were left over after removing those that had an intersect with a gffmerge -MKQ models (mikado, red deer, cattle )
bedtools intersect -v -wo -b <(awk '$3=="locus"' Bos_ReddeerReductionTEST.gff|grep  "," ) -a ../05_EvaluatePrediction/Fixedaugustus.hints.gff3 |awk '$3=="gene"' |sort|uniq|wc
  26509  238581 1478904



```

#Which are busco genes and should not be eliminated?
```
#make gene (mRNA)list of complete buscos
awk '{print $3}' Cetfull_table.tsv |awk 'NR>3'|sort|uniq|cat <(awk '{print $3}' Eukfull_table.tsv) -|sort|uniq|sed '/^$/d' |grep -v "version" >BuscoGenes.list
wc BuscoGenes.list
 25255  25255 493072 BuscoGenes.list


#make it so grep for exact matches of gene names
###less Bos_ReddeerReduction.gff |awk '$3=="locus"' |sed 's/,/\t/g' |sed 's/=/=\t/g' |sed 's/;/\t/g' >Bos_ReddeerReductionGREPMOD.gff
###### this was killed
### grep -w -f BuscoGenes.list Bos_ReddeerReductionGREPMOD.gff >BuscoBos_ReddeerReductionGREPMOD.gff &


#Test how we can extract genes/mrnas
less Bos_ReddeerReduction.gff|awk '$3=="locus"' |cut -f 9 |sed 's/ID=//g'|sed 's/;/\t/1' |sed 's/transcripts=/\t/g' |cut -f 1,3 |sed 's/,/\t/g' |awk -F"\t" 'NF>2' |awk '{print $2,$1}' |tr " " "\t" >Bos_ReddeerReductionDualEvidLoci.list

#All braker
less Bos_ReddeerReduction.gff|awk '$3=="locus"' |cut -f 9 |sed 's/ID=//g'|sed 's/;/\t/1' |sed 's/transcripts=/\t/g' |cut -f 1,3 |grep "mRNA" |sed 's/,/\t/g'  |awk '{print $2,$1}' |tr " " "\t" >Bos_ReddeerReductionBrakerEvidLoci.list

 wc -l Bos_ReddeerReductionBrakerEvidLoci.list                                                                         
 56180

#which of these are busco genes?
cat BuscoGenes.list Bos_ReddeerReductionBrakerEvidLoci.list| sort -k1,1 -u |awk 'NF==2' |cut -f 2 |sort|uniq|wc
49697   49697  695758
cat BuscoGenes.list Bos_ReddeerReductionDualEvidLoci.list| sort -k1,1 -u |awk 'NF==2' |cut -f 2 |sort|uniq|wc
12629   12629  176806

cat Bos_ReddeerReductionDualEvidLoci.list Bos_ReddeerReductionBrakerEvidLoci.list |sort|uniq|cat - BuscoGenes.list|sort -k1,1 -u |awk 'NF==2'|wc
  59959  119918 1618471

Combine dual evidence,braker, Cet BUSCOs, and euk BUSCOs

#number of genes with busco complete?
cat <(less Eukfull_table.tsv |awk -F"\t" 'NF>3' |cut -f 1,3 |sort -k1,1 -u |grep -v "#" |cut -f 2 ) <(less Cetfull_table.tsv|awk -F"\t" 'NF>3' |cut -f 1,3 |sort -k1,1 -u |grep -v "#" |cut -f 2)|sort|uniq|wc
  12577   12577  238575

#write the above to file
cat <(less Eukfull_table.tsv |awk -F"\t" 'NF>3' |cut -f 1,3 |grep -v "#" |cut -f 2 ) <(less Cetfull_table.tsv|awk -F"\t" 'NF>3' |cut -f 1,3  |grep -v "#" |cut -f 2)|sort|uniq >BuscoCompleteMrnaListEukCet.list
wc -l BuscoCompleteMrnaListEukCet.list
25253 BuscoCompleteMrnaListEukCet.list



cat <(cut -f 1 Bos_ReddeerReductionDualEvidLoci.list) BuscoCompleteMrnaListEukCet.list |sort|uniq >BuscoDualEvidencemRNA.list


#whole genome database to grep
less Bos_ReddeerReduction.gff |awk '$3=="mRNA" || $3=="transcript"' |cut -f 9 |sed 's/;/\t/1' | sed 's/locus=/locus=\t/g' |sed 's/ID=//g' |cut -f 1,3 |sort|uniq>AllmRNAGenesBos_ReddeerReduction.list
(mikado) [remkv6@nova031 02_mergeMikadoBraker]$ wc -l AllmRNAGenesBos_ReddeerReduction.list
154729 AllmRNAGenesBos_ReddeerReduction.list

awk -F"\t" '{arr[$1]=arr[$1] "\t" $2}END{for(i in arr)print i,arr[i]}'  <(awk '{print $1"\tHighQuality"}' BuscoDualEvidencemRNA.list) AllmRNAGenesBos_ReddeerReduction.list |awk '$2=="HighQuality" {print $1"\t"$3}' >HighConfidence4MikadoGrep.list

mikado util grep HighConfidence4MikadoGrep.list Bos_ReddeerReduction.gff HighConfidenceBos_ReddeerReduction.gff

awk '$3=="locus"' HighConfidenceBos_ReddeerReduction.gff |wc
20646  185814 5165459
gffread -VHEJ HighConfidenceBos_ReddeerReduction.gff  -g ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa  -t mRNA -x  HighConfidenceBos_ReddeerReductionVHEJ_transcripts.fasta -y  HighConfidenceBos_ReddeerReductionVHEJ_proteins.fasta

#remove those that do not produce correctly coding proteins
grep ">" HighConfidenceBos_ReddeerReductionVHEJ_proteins.fasta|awk '{print $1}' |sed 's/>//g' |grep -v -w -f HighConfidence4MikadoGrep.list >ProperCodingHighConfidence4MikadoGrep.list
awk -F"\t" '{arr[$1]=arr[$1] "\t" $2}END{for(i in arr)print i,arr[i]}'  <(awk '{print $1"\tHighQuality"}' ProperCodingHighConfidence4MikadoGrep.list) AllmRNAGenesBos_ReddeerReduction.list |awk '$2=="HighQuality" {print $1"\t"$3}' >GenesmRNAProperCodingHighConfidence4MikadoGrep.list
mikado util grep  GenesmRNAProperCodingHighConfidence4MikadoGrep.list Bos_ReddeerReduction.gff ProperCodingHighConfidenceBos_ReddeerReduction.gff

gffread -VHEJ ProperCodingHighConfidenceBos_ReddeerReduction.gff  -g ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa  -t mRNA -x  ProperCodingHighConfidenceBos_ReddeerReductionVHEJ_transcripts.fasta -y  ProperCodingHighConfidenceBos_ReddeerReductionVHEJ_proteins.fasta

awk '$3=="mRNA" ||$3=="transcript"' ProperCodingHighConfidenceBos_ReddeerReduction.gff |wc
  33605  302445 4803183
  grep -c ">" ProperCodingHighConfidenceBos_ReddeerReductionVHEJ_transcripts.fasta
  33605

  awk '$3=="locus"' ProperCodingHighConfidenceBos_ReddeerReduction.gff |wc
    17873  160857 4515020

```
# Summary
Essentially I ran a gffread -VHEJKMQ  to collapse transcripts between gene models from braker and mikado, and alignments from by cattle and red deer alignments.  I retained every gene that had dual evidence from two of the above four sources. I combined these mRNA/genes with all identified complete busco mRNAs (both eukaryotes and even-toed ungulates (cet...?)). There are 26,509 braker genes that did not have dual evidence.


## Single source data exploration
```
mikado util grep -v GenesmRNAProperCodingHighConfidence4MikadoGrep.list Bos_ReddeerReduction.gff LowConfidenceBos_ReddeerReduction.gff

awk '$3=="locus"' LowConfidenceBos_ReddeerReduction.gff |wc
  93151  838359 12926698

gffread -VHEJ LowConfidenceBos_ReddeerReduction.gff  -g ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa  -t mRNA -x  LowConfidenceBos_ReddeerReductionVHEJ_transcripts.fasta -y  LowConfidenceBos_ReddeerReductionVHEJ_proteins.fasta

grep -c ">" LowConfidenceBos_ReddeerReductionVHEJ_proteins.fasta
93996

#number of loci in the LC dataset?  
awk '$3=="locus"' LowConfidenceBos_ReddeerReduction.gff |wc
  93151  838359 12926698

#number of loci that have a 10% overlap with an EDTA repeat, While these are confidently repeats, the size of gene content required needed
bedtools intersect -wo -f .1 -a LowConfidenceBos_ReddeerReduction.gff -b ../05_EvaluatePrediction/AllAnno.gff |awk '$3=="locus"' |awk '{print $9}' |sed 's/;/\t/g' |cut -f 1 |sed 's/ID=//g' |sort|uniq|wc
26816   26816  375424

#number of loci in the LC dataset?  
awk '$3=="locus"' LowConfidenceBos_ReddeerReduction.gff |wc
  93151  838359 12926698

#number of loci with any repeat overlap
bedtools intersect -wo -a LowConfidenceBos_ReddeerReduction.gff -b ../05_EvaluatePrediction/AllAnno.gff |awk '$3=="locus"' |awk '{print $9}' |sed 's/;/\t/g' |cut -f 1 |sed 's/ID=//g' |sort|uniq|wc
  49465   49465  692510

#How many Low confidence mrnas overlap with with Expressed mikado loci?
bedtools intersect -f .3 -wo -a LowConfidenceBos_ReddeerReduction.gff -b ../05_EvaluatePrediction/Expressedmikado.loci.gff3 |awk '$3=="mRNA"' |cut -f 9 |sed 's/locus=/locus=\t/g' |cut -f 2 |sort|uniq|wc
  73232   73232 1025248

#how many are likely silenced?
#comparison with existing annotation
bedtools intersect -f .3 -v -wo -a LowConfidenceBos_ReddeerReduction.gff -b ../05_EvaluatePrediction/Expressedmikado.loci.gff3 |awk '$3=="mRNA"' |wc
 23309  209781 2890600

#How many loci are expressed and not associated with repeats? ## Old calculation, waiting for new expression values
 bedtools intersect -f .3 -v -wo -a LowConfidenceBos_ReddeerReduction.gff -b ../05_EvaluatePrediction/Expressedmikado.loci.gff3 |awk '$3=="mRNA"' |cut -f 1-9 |bedtools intersect -wo -v -a - -b ../05_EvaluatePrediction/AllAnno.gff | awk '{print $9}' |sed 's/locus=/\t/g' |cut -f 2 |sort|uniq|wc
   13982   13982  195748

#grab the zero expression genes
awk -F"\t" '{arr[$1]=arr[$1] "\t" $2}END{for(i in arr)print i,arr[i]}' AllmRNAGenesBos_ReddeerReduction.list  <(awk '{print $1"\tNOUREADS"}'   AllmRNAGenesBos_ReddeerReduction.list) |awk '$3=="NOUREADS"' |cut -f 1,2 >NoExpressionList4MikadoGrep.list

 bedtools intersect -f .3 -v -wo -a LowConfidenceBos_ReddeerReduction.gff -b ../05_EvaluatePrediction/Expressedmikado.loci.gff3 |awk '$3=="mRNA"' |cut -f 1-9 |bedtools intersect -wo -v -a - -b ../05_EvaluatePrediction/AllAnno.gff | awk '{print $9}' |sed 's/locus=/\t/g' |cut -f 2 |sort|uniq|wc

mikado util grep -v NoExpressionList4MikadoGrep.list Bos_ReddeerReduction.gff AllExpressed.Bos_ReddeerReduction.gff
mikado util NoExpressionList4MikadoGrep.list Bos_ReddeerReduction.gff NoExpression.Bos_ReddeerReduction.gff

```

### Run Busco on High Confidence genes
```
#/work/GIF/remkv6/Elk/31_busco/02_Busco4/01_ProteinBusco4


ml miniconda3; source activate busco4; export BUSCO_CONFIG_FILE=/home/remkv6/.conda/envs/busco4/config/config.ini ; busco -i ProperCodingHighConfidenceBos_ReddeerReductionVHEJ_proteins.fasta --auto-lineage-euk -o ProperCodingHighConfidenceBos_ReddeerReductionVHEJ_proteinProtBUSCO4 -m prot -c 15 -f --augustus_species CervusCanadensis3


--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:97.7%[S:41.2%,D:56.5%],F:1.2%,M:1.1%,n:255     |
|249    Complete BUSCOs (C)                       |
|105    Complete and single-copy BUSCOs (S)       |
|144    Complete and duplicated BUSCOs (D)        |
|3      Fragmented BUSCOs (F)                     |
|3      Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------

--------------------------------------------------
|Results from dataset cetartiodactyla_odb10       |
--------------------------------------------------
|C:92.1%[S:40.2%,D:51.9%],F:1.7%,M:6.2%,n:13335   |
|12285  Complete BUSCOs (C)                       |
|5367   Complete and single-copy BUSCOs (S)       |
|6918   Complete and duplicated BUSCOs (D)        |
|227    Fragmented BUSCOs (F)                     |
|823    Missing BUSCOs (M)                        |
|13335  Total BUSCO groups searched               |
--------------------------------------------------

```
