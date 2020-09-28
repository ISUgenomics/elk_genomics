# Need to identify the issues with interproscan annotations
THere were only ~200 interpro annotations in the final gff, so something broke.


It appears to be an issue with how I had to name the gene names for the interproscan to run.  I had to remove all of the periods from the gene names.  


### get necessary files
```
#/work/GIF/remkv6/Elk/35_InvestigateInterproFlawedAnnots


# copied geneRenaming.map from nova /work/gif/remkv6/Olsen/Elk/08_RenameAgain
#this had gene names that were removed later (low confidence), so need to cull those
less OrderedfixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 |awk '$3=="mRNA"' |sed 's/;/\t/g' |cut -f 9 |sort|uniq >mRNAPresentLists

#it has the correct number of mRNA's
wc -l mRNAPresentLists
33433 mRNAPresentLists


#copied file to condo at /work/GIF/remkv6/Elk/35_InvestigateInterproFlawedAnnots
#cull gene names
less mRNAPresentLists |sed 's/ID=//g' |sed 's/|/\t/g' |cut -f 1 |while read line; do echo "awk '\$2==\""$line"\"' geneRenaming.map >>NewGeneNamesOnlyPresentInAnnotation.map";done >NewRenamer.sh
sh NewRenamer.sh
Creates # NewGeneNamesOnlyPresentInAnnotation.map

wc NewGeneNamesOnlyPresentInAnnotation.map
  33433   66866 1257830 NewGeneNamesOnlyPresentInAnnotation.map


#interproscan file with all annotations
ln -s ../25_Interpro/interproAnnot.tab1

#blast for all genes
ln -s ../32_CombineFunctionalAnnotations/Blasts.tab
```


### Rename genes in the interproannotation and get counts
```
#because the names are different, I need to sort out how to make changes to each gene name without affecting the annotation definition.

# four sets of annotations: mikado, mRNA, ENS, RED
awk 'substr($1,1,3)=="mik"' interproAnnot.tab1 |sed 's/\./_/1' |sed 's/\./_/1' >mikMatched.interproAnnot.tab
awk 'substr($1,1,3)=="RED"' interproAnnot.tab1 |sed 's/\./_/1' |sed 's/\./_/1' |sed 's/\./_/1'| >redMatched.interproAnnot.tab1
awk 'substr($1,1,3)=="ENS"' interproAnnot.tab1 |sed 's/\./_/1' |sed 's/\./_/1'| >ENSMatched.interproAnnot.tab1
awk 'substr($1,1,3)=="mrn"' interproAnnot.tab1  >mrnaMatched.interproAnnot.tab1

#concatenate them all
cat mrnaMatched.interproAnnot.tab1 redMatched.interproAnnot.tab1 mikMatched.interproAnnot.tab1 ENSMatched.interproAnnot.tab1 >CATAllModNamesinterproAnnot.tab1




less NewGeneNamesOnlyPresentInAnnotation.map |grep "-" |cut  -f 1 |while read line ; do echo "awk '\$1==\""$line"\"' CATAllModNamesinterproAnnot.tab1 >>AnnotatedmRNAsInFinalAnnot.tab1";done >grabAnnots.sh
sh grabAnnots.sh

wc AnnotatedmRNAsInFinalAnnot.tab1
  32917    65834 27530067 AnnotatedmRNAsInFinalAnnot.tab1


#change the naming scheme to Cercan
map_data_ids NewGeneNamesOnlyPresentInAnnotation.map CATAllModNamesinterproAnnot.tab1

#number of mRNA's annotated
less CATAllModNamesinterproAnnot.tab1  |grep -c "Cercan"
  32917
#number of genes annotated
less AnnotatedmRNAsInFinalAnnot.tab1|sed 's/-/\t/1' |cut -f 1|sort|uniq|wc
   17607   17607  281712

```

### get the total of all functionally annotated genes and mRNA's
```
#Need to perform the same operations on the blast results to get accurate counts of annotated genes and mRNA's

awk 'substr($1,1,3)=="ENS"' Blasts.tab |sed 's/\./_/1' |sed 's/\./_/1' >ENSMatched.Blasts.tab
awk 'substr($1,1,3)=="RED"' Blasts.tab |sed 's/\./_/1' |sed 's/\./_/1' |sed 's/\./_/1' >REDMatched.Blasts.tab
awk 'substr($1,1,3)=="mik"' Blasts.tab |sed 's/\./_/1'|sed 's/\./_/1'   >mikMatched.Blasts.tab
awk 'substr($1,1,3)=="mRN"' Blasts.tab    >mRNAMatched.Blasts.tab
cat mikMatched.Blasts.tab mRNAMatched.Blasts.tab ENSMatched.Blasts.tab REDMatched.Blasts.tab >AllMatchingBlasts.tab


#grab only the ids that are in the final annotation.
less NewGeneNamesOnlyPresentInAnnotation.map |grep "-" |cut  -f 1 |while read line ; do echo "awk '\$1==\""$line"\"' AllMatchingBlasts.tab >>AnnotatedmRNAsInFinalBLAST.tab1";done >grabAnnotsBlast.sh

sh grabAnnotsBlast.sh &

wc AnnotatedmRNAsInFinalBLAST.tab1
   33089  1227941 15321603 AnnotatedmRNAsInFinalBLAST.tab1

map_data_ids -delimit " \t" NewGeneNamesOnlyPresentInAnnotation.map AnnotatedmRNAsInFinalBLAST.tab1

grep -c "Cercan" AnnotatedmRNAsInFinalBLAST.tab1
33089

cat AnnotatedmRNAsInFinalBLAST.tab1 AnnotatedmRNAsInFinalAnnot.tab1 |awk '{print $1}' |sort|uniq|wc
  33355   33355  633809

```
