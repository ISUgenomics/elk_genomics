# Need to rename scaffolds for distribution

### run minimap
```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/27_RenameScaffolds/02_minimap

#run minimap alignment
module load minimap2/2.6; minimap2 -t 40 -ax asm5 FinalGenomePilonReducedSoftMasked.fa GCA_002197005.1_CerEla1.0_genomic.fna > OurElkvsExisting.sam

#filter down the hits
less -S OurElkvsExisting.sam |awk '$1!="#" && substr($1,1,1)!="@"{print $1,$3}' |sort |uniq -c |sort -k1,1nr >chr2chr.list &

#filter down the good hits
less chr2chr.list |awk '$1>100' |sort -k3,3V |less
```
Looked up the genome assembly for Cervus elaphus hippelaphus, and used vlookup to assign names correctly.  



### info table from excel
```
Chr	AlignmentCounts	Scaffold Name	Our Scaffold Name
Chr_1	7203	CM008008.1	HiC_scaffold_11_pilon
Chr_2	4599	CM008009.1	HiC_scaffold_33_pilon
Chr_3	1115	CM008010.1	HiC_scaffold_9_pilon
Chr_3	5519	CM008010.1	HiC_scaffold_10_pilon
Chr_4	5545	CM008011.1	HiC_scaffold_22_pilon
Chr_5	11230	CM008012.1	HiC_scaffold_21_pilon
Chr_6	5445	CM008013.1	HiC_scaffold_8_pilon
Chr_7	4790	CM008014.1	HiC_scaffold_4_pilon
Chr_8	3470	CM008015.1	HiC_scaffold_15_pilon
Chr_9	9579	CM008016.1	HiC_scaffold_35_pilon
Chr_10	3687	CM008017.1	HiC_scaffold_5_pilon
Chr_11	9373	CM008018.1	HiC_scaffold_14_pilon
Chr_12	8727	CM008019.1	HiC_scaffold_23_pilon
Chr_13	6263	CM008020.1	HiC_scaffold_6_pilon
Chr_14	7283	CM008021.1	HiC_scaffold_26_pilon
Chr_15	8296	CM008022.1	HiC_scaffold_25_pilon
Chr_16	4243	CM008023.1	HiC_scaffold_29_pilon
Chr_17	6159	CM008024.1	HiC_scaffold_20_pilon
Chr_18	10947	CM008025.1	HiC_scaffold_17_pilon
Chr_19	8843	CM008026.1	HiC_scaffold_12_pilon
Chr_20	9480	CM008027.1	HiC_scaffold_16_pilon
Chr_21	7693	CM008028.1	HiC_scaffold_19_pilon
Chr_22	4692	CM008029.1	HiC_scaffold_9_pilon
Chr_23	7036	CM008030.1	HiC_scaffold_34_pilon
Chr_24	5203	CM008031.1	HiC_scaffold_31_pilon
Chr_25	7295	CM008032.1	HiC_scaffold_28_pilon
Chr_26	3733	CM008033.1	HiC_scaffold_30_pilon
Chr_27	6158	CM008034.1	HiC_scaffold_13_pilon
Chr_28	6320	CM008035.1	HiC_scaffold_24_pilon
Chr_29	5820	CM008036.1	HiC_scaffold_18_pilon
Chr_30	8738	CM008037.1	HiC_scaffold_3_pilon
Chr_31	6145	CM008038.1	HiC_scaffold_27_pilon
Chr_32	4623	CM008039.1	HiC_scaffold_32_pilon
Chr_33	7482	CM008040.1	HiC_scaffold_7_pilon
Chr_33	1380	CM008040.1	HiC_scaffold_15_pilon
Chr_X	17548	CM008041.1	HiC_scaffold_1_pilon
Chr_Y	312	CM008042.1	HiC_scaffold_2_pilon
```
Decided to go with a padded numbering scheme, pasted from excel to scaffoldRename.list
```
#create naming scheeme from list above
awk '{print "sed -i s/"$1"/"$2"/g FinalGenomePilonReducedSoftMaskedFINALSCAFFNAMES.fa" }' scaffoldRename.list  |sed "s|s/|'s/|g" |sed "s|/g|/g'|g" >ScaffoldRenamer.sh

#submitted to node, as it was taking quite a while to loop through this file 180 some times.  

```
