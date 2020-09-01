# Renaming of scaffolds and genes for Cervus canadensis

### Gather all pertinent files
```
#/work/gif/remkv6/Olsen/Elk/08_RenameAgain

#Copied from Ceres 22_repeatmasker
FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa

#Copied from Ceres 22_repeatmasker
FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa.out.gff


fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3


#EDTA transposon annotator track
cp ../05a_RenameScaffsNGenes/01_OtherGFFs/tidyAllAnnoEDTA.gff.gz .

#Expressed genes annotated as transposable elements
cp ../05_EvaluatePrediction/TransposableElementsOrderedAnnotatedGeneModels_sorted.gff .

#Genes of low confidence
cp ../05a_RenameScaffsNGenes/01_OtherGFFs/tidyOrderedGTLOWConfidencetest.gff3.gz .

#Genes that overlap with repeats
cp ../05a_RenameScaffsNGenes/01_OtherGFFs/tidyOrderedGTRepetitiveGenes.gff3.gz .

#Final gene predictions
cp ../05a_RenameScaffsNGenes/fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3.gz .


for f in *gz; do gunzip $f;done
```

### Create map of name changes "old \t new"
```
#/work/gif/remkv6/Olsen/Elk/08_RenameAgain

#Create a map names to change
#needs modified to work.  A names file will have to be ran twice, since the style of naming is to be kept consistent
Old2NewNamesByChrSize.map

awk '{print "Chromosome_"$2"\tCHUhromosome_"$1}' Old2NewNamesByChrSize.map >FixedOld2NewNamesByChrSize.map

#for run 1
awk '{print "Chromosome_"$2"\tC22hromosome_"$1}' Old2NewNamesByChrSize.map >Fixed1Old2NewNamesByChrSize.map
#forgot to pad with a zero.  vi used

#for run 2
awk '{print "C22hromosome_"$1"\tChromosome_"$1}' Old2NewNamesByChrSize.map >Fixed2Old2NewNamesByChrSize.map
#forgot to pad with a zero. vi used
```

## Rename the chromosomes

### Expressed Transposable elements -- Rename scaffolds in files
```
#/work/gif/remkv6/Olsen/Elk/08_RenameAgain

#expressed transposable elements
## this was missing the chromosome at the beginning of chromosome name. Unfortunately it adds it to scaffolds too
sed  -i 's/^/Chromosome_/g'  TransposableElementsOrderedAnnotatedGeneModels_sorted.gff

#Modified the map file so it did not have padding, which is more difficult to add to the gff file.
cp Fixed1Old2NewNamesByChrSize.map OnlyTransposableElementsOrderedAnnotatedGeneModels_sorted.map


#rename the first set of names
map_data_ids OnlyTransposableElementsOrderedAnnotatedGeneModels_sorted.map TransposableElementsOrderedAnnotatedGeneModels_sorted.gff
#fix those huge chromosome numbers that should be scaffolds

sed  -i 's/Chromosome_/Scaffold_/g'  TransposableElementsOrderedAnnotatedGeneModels_sorted.gff

map_data_ids Fixed2Old2NewNamesByChrSize.map TransposableElementsOrderedAnnotatedGeneModels_sorted.gff

perl ../02_mergeMikadoBraker/gff3sort/gff3sort.pl --precise --chr_order natural TransposableElementsOrderedAnnotatedGeneModels_sorted.gff >OrderedTransposableElementsOrderedAnnotatedGeneModels_sorted.gff3
bgzip OrderedTransposableElementsOrderedAnnotatedGeneModels_sorted.gff3
tabix -p gff OrderedTransposableElementsOrderedAnnotatedGeneModels_sorted.gff3.gz
```


### rename the fasta scaffolds of the genome
```
#/work/gif/remkv6/Olsen/Elk/08_RenameAgain

map_fasta_ids Fixed1Old2NewNamesByChrSize.map FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa

map_fasta_ids Fixed2Old2NewNamesByChrSize.map FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa

samtools faidx FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa
```

### rename the repeatmasker gff scaffolds
```
#/work/gif/remkv6/Olsen/Elk/08_RenameAgain

map_data_ids Fixed1Old2NewNamesByChrSize.map FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa.out.gff
map_data_ids Fixed2Old2NewNamesByChrSize.map FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa.out.gff

#convert to gff3
less FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa.out.gff |awk  '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"="$10,$11,$12}' |sed 's/"//g' |grep -v "#" |sort -k1,1V -k4,5n >Repeatmasker.gff3

bgzip Repeatmasker.gff3
tabix -p gff Repeatmasker.gff3.gz

```

### EDTA TRANSPOSON ANNOTATOR  scaffold rename
```

#this one needs to have column 2 changed also, as it has indices in the header that would be incorrect and may mislead

map_data_ids Fixed1Old2NewNamesByChrSize.map tidyAllAnnoEDTA.gff
map_data_ids Fixed2Old2NewNamesByChrSize.map tidyAllAnnoEDTA.gff

#just getting rid of the header
grep -v "#" tidyAllAnnoEDTA.gff  |sort -k1,1V -k4,5n >HeaderDashedtidyAllAnnoEDTA.gff
bgzip HeaderDashedtidyAllAnnoEDTA.gff
tabix -p gff HeaderDashedtidyAllAnnoEDTA.gff.gz
```

### Low confidence genes rename scaffolds
```
#this one needs to have column 2 changed also, as it has indices in the header that would be incorrect and may mislead

map_data_ids Fixed1Old2NewNamesByChrSize.map tidyOrderedGTLOWConfidencetest.gff3
map_data_ids Fixed2Old2NewNamesByChrSize.map tidyOrderedGTLOWConfidencetest.gff3

grep -v "#" tidyOrderedGTLOWConfidencetest.gff3  >HeaderDashedtidyOrderedGTLOWConfidencetest.gff3

perl ../02_mergeMikadoBraker/gff3sort/gff3sort.pl --precise --chr_order natural HeaderDashedtidyOrderedGTLOWConfidencetest.gff3 >OrderedHeaderDashedtidyOrderedGTLOWConfidencetest.gff3
bgzip OrderedHeaderDashedtidyOrderedGTLOWConfidencetest.gff3
tabix -p gff OrderedHeaderDashedtidyOrderedGTLOWConfidencetest.gff3.gz
```


### rename scaffolds in genes with repetitive overlap
```
#this one needs to have column 2 changed also, as it has indices in the header that would be incorrect and may mislead

map_data_ids Fixed1Old2NewNamesByChrSize.map tidyOrderedGTRepetitiveGenes.gff3
map_data_ids Fixed2Old2NewNamesByChrSize.map tidyOrderedGTRepetitiveGenes.gff3

grep -v "#" tidyOrderedGTRepetitiveGenes.gff3 >HeaderDashedtidyOrderedGTRepetitiveGenes.gff3

perl ../02_mergeMikadoBraker/gff3sort/gff3sort.pl --precise --chr_order natural HeaderDashedtidyOrderedGTRepetitiveGenes.gff3 >OrderedHeaderDashedtidyOrderedGTLOWConfidencetest.gff3
bgzip HeaderDashedtidyOrderedGTRepetitiveGenes.gff3
tabix -p gff HeaderDashedtidyOrderedGTRepetitiveGenes.gff3.gz
```


### Rename the scaffolds in the main gene gff, High confidence genes
```
map_data_ids Fixed1Old2NewNamesByChrSize.map fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3
map_data_ids Fixed2Old2NewNamesByChrSize.map fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3

perl ../02_mergeMikadoBraker/gff3sort/gff3sort.pl --precise --chr_order natural fixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3 >OrderedfixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3
bgzip OrderedfixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3
tabix -p gff OrderedfixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3.gz
```


### move folders and rename files to transfer
```
mkdir 01_Move2Box
 cd 01_Move2Box/
cp ../HeaderDashedtidyAllAnnoEDTA.gff.gz EDTATransposons.gff.gz
cp ../HeaderDashedtidyAllAnnoEDTA.gff.gz.tbi EDTATransposons.gff.gz.tbi
cp ../HeaderDashedtidyOrderedGTRepetitiveGenes.gff3.gz GenesWithRepetitiveOverlap.gff3.gz
cp ../HeaderDashedtidyOrderedGTRepetitiveGenes.gff3.gz.tbi GenesWithRepetitiveOverlap.gff3.gz.tbi
cp ../ OrderedHeaderDashedtidyOrderedGTLOWConfidencetest.gff3.gz LowConfidenceGenes.gff3.gz
cp ../OrderedHeaderDashedtidyOrderedGTLOWConfidencetest.gff3.gz LowConfidenceGenes.gff3.gz
cp ../OrderedHeaderDashedtidyOrderedGTLOWConfidencetest.gff3.gz.tbi LowConfidenceGenes.gff3.gz.tbi
cp ../OrderedTransposableElementsOrderedAnnotatedGeneModels_sorted.gff3.gz ExpressedTransposableElements.gff3.gz
cp ../OrderedTransposableElementsOrderedAnnotatedGeneModels_sorted.gff3.gz.tbi ExpressedTransposableElements.gff3.gz.tbi
cp ../OrderedfixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3.gz HighConfidenceGenes.gff3.gz
cp ../OrderedfixedNoHashCorrectedNoPeriodsRevisedGeneidTranscriptsRemovedOrderedGTNOTEHighConfidencetest.gff3.gz.tbi HighConfidenceGenes.gff3.gz.tbi
cp ../Repeatmasker.gff3.gz Repeatmasker.gff3.gz
cp ../Repeatmasker.gff3.gz.tbi Repeatmasker.gff3.gz.tbi
cp ../FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa CervusCanadensisGenome.fa
cp ../FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa.fai CervusCanadensisGenome.fa.fai

rclone copy 01_Move2Box/ OlsenElkBox:GIF_shared/Projects/SteveOlsen/ELK/
```
