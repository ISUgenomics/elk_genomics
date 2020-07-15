# need to rename the gffs that display in jbrowse to match genome now. 

```
#/work/gif/remkv6/Olsen/Elk/05a_RenameScaffsNGenes/01_OtherGFFs

for f in ../../05_EvaluatePrediction/FinalGFFs/* ; do cp -rf $f .; done

rm *tbi
rm FinalGenomePilonReducedSoftMaskedRecode.*
gunzip OrderedGTLOWConfidencetest.gff3.gz
gunzip OrderedGTRepetitiveGenes.gff3.gz
gunzip OrderedTransposableElementsOrderedAnnotatedGeneModels_sorted.gff.gz
gunzip AllAnnoEDTA.gff.gz

singularity shell "/opt/rit/singularity/images/maker/2.31.10_3.1/maker.simg"
for f in *gff; do map_data_ids FixedScaffoldNames.map $f;done
for f in *gff3; do map_data_ids FixedScaffoldNames.map $f;done

gt gff3 -tidy -sortlines <(awk 'substr($1,1,1)!="#"' AllAnnoEDTA.gff) >tidyAllAnnoEDTA.gff



```
