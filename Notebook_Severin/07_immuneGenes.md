# Immune gene identification

* Sep 2, 2020

I am going to take a database of immune genes and blast them to Elk and bovine genomes. Then see if what may be missing in one or the other.

## Objectives

* BLAST immune genes to Elk
* BLAST immune genes to Bovine
* Create Venn diagram of overlap of presence and absence
* Take closer look at what is missing in Elk vs Bovine

## Transfer to Condo

Nova is oversubscribed so I am going to perform these analyses on Condo in ptmp

* /ptmp/GIF/severin/Olsen

```bash
rsync -avz -e ssh severin@novadtn.its.iastate.edu:/work/gif/remkv6/Olsen/Elk/08_RenameAgain/01_Move2Box .
```

## Download database of Immune genes

* http://www.imgt.org/download/GENE-DB/README.txt

```
wget http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP
```

Replace spaces and `|` with `_`

```
perl -i -pe 's/[\| ]/_/g' IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP
```

Move it to an easier name

```
mv IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP IMGTGENEDB.fasta
```

## Blast Immune genes against Elk genome

```
ln -s 01_Move2Box/CervusCanadensisGenome.fa
nextflow run isugifnf/blast --query IMGTGENEDB.fasta --genome CervusCanadensisGenome.fa --chunkSize 1000 --app tblastn -profile condo --options '-evalue 1e-3'
```

## BLAST just the bos_taurus immune genes to Elk

```
fastaMulti2singleLine.pl IMGTGENEDB.fasta | paste - - | awk '{print $1,$2}' | grep -i bos_taurus | tr ' ' '\n' > IMGTGENEDB_bos.fasta
```

```
module load gcc/7.3.0-xegsmw4 nextflow
nextflow run isugifnf/blast --query IMGTGENEDB_bos.fasta --dbDir "./" --dbName "CervusCanadensisGenome" --chunkSize 100 --app tblastn -profile condo --outdir bosImmune_2_Elk
N E X T F L O W  ~  version 20.07.1
Launching `isugifnf/blast` [marvelous_church] - revision: 273a14a8df [master]
executor >  slurm (9)
executor >  slurm (9)
[a8/7cb1ae] process > software_check [100%] 1 of 1 ✔
[07/47b137] process > runBlast (5)   [100%] 8 of 8 ✔
```

grab the best hits

* condo: /ptmp/GIF/severin/Olsen/bosImmune_2_Elk

```
cd bosImmune_2_Elk
cat blast_output_combined.txt | firstInstanceOf.awk > bestBlast_output.txt
```

## BLAST just the bos_taurus immune genes to Bos taurus

Get the genome we used in the paper

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz
mv GCF_002263795.1_ARS-UCD1.2_genomic.fna bostaurus.fasta
```

BLAST

```
nextflow run isugifnf/blast --query IMGTGENEDB_bos.fasta  --genome bostaurus.fasta --chunkSize 100 --app tblastn -profile condo --outdir bosImmune_2_Bos --threads 8
```

grab the best hits

* condo: /ptmp/GIF/severin/Olsen/bosImmune_2_Bos

```
cd bosImmune_2_Bos
cat blast_output_combined.txt | firstInstanceOf.awk > bestBlast_output.txt
```

## Best hits only found in bos taurus and not in Elk

```
cat bosImmune_2_Bos/bestBlast_output.txt  bosImmune_2_Bos/bestBlast_output.txt bosImmune_2_Elk/bestBlast_output.txt | cut -f 1 | sort |uniq -c | sort -rn | awk '$1<3'
     2 NW_001494075_IGHJ1-2*03_Bos_taurus_Hereford_ORF_J-REGION_39638..39693_56_nt_2_______18_AA_18+0=18___rev-compl_
     2 NW_001494075_IGHJ1-1*03_Bos_taurus_Hereford_ORF_J-REGION_39812..39874_63_nt_3_______20_AA_20+0=20___rev-compl_
     2 KT723008_IGHJ2-2*01_Bos_taurus_Holstein_ORF_J-REGION_509129..509184_56_nt_2_______18_AA_18+0=18_____
     2 KT723008_IGHJ2-1*01_Bos_taurus_Holstein_ORF_J-REGION_508948..509010_63_nt_3_______20_AA_20+0=20_____
     2 KT723008_IGHJ1-2*01_Bos_taurus_Holstein_ORF_J-REGION_320315..320370_56_nt_2_______18_AA_18+0=18_____
     2 KT723008_IGHD1-3*01_Bos_taurus_Holstein_ORF_D-REGION_422641..422759_119_nt_1_______39_AA_39+0=39_____
     2 KT723008_IGHD*02_Bos_taurus_Holstein_ORF_H2_g,527527..527591_66_nt_1_+1_-1___22_AA_22+0=22_____
     2 IMGT000049_TRDC*01_Bos_taurus_Hereford_F_EX2_g,3231015..3231088_75_nt_1_+1_-1___25_AA_25+0=25_____
     2 IMGT000049_TRAJ8-2*01_Bos_taurus_Hereford_F_J-REGION_3301039..3301098_60_nt_3_______19_AA_19+0=19_____
     2 IMGT000049_TRAJ8-1*01_Bos_taurus_Hereford_F_J-REGION_3303827..3303886_60_nt_3_______19_AA_19+0=19_____
     2 IMGT000049_TRAJ6*01_Bos_taurus_Hereford_F_J-REGION_3305951..3306012_62_nt_2_______20_AA_20+0=20_____
     2 IMGT000049_TRAJ57*01_Bos_taurus_Hereford_F_J-REGION_3248933..3248995_63_nt_3_______20_AA_20+0=20_____
     2 IMGT000049_TRAJ56*01_Bos_taurus_Hereford_F_J-REGION_3249567..3249624_58_nt_1_______19_AA_19+0=19_____
     2 IMGT000049_TRAJ5*01_Bos_taurus_Hereford_F_J-REGION_3307141..3307200_60_nt_3_______19_AA_19+0=19_____
     2 IMGT000049_TRAJ49*01_Bos_taurus_Hereford_F_J-REGION_3259424..3259479_56_nt_2_______18_AA_18+0=18_____
     2 IMGT000049_TRAJ48*01_Bos_taurus_Hereford_F_J-REGION_3260307..3260369_63_nt_3_______20_AA_20+0=20_____
     2 IMGT000049_TRAJ46*01_Bos_taurus_Hereford_F_J-REGION_3263071..3263133_63_nt_3_______20_AA_20+0=20_____
     2 IMGT000049_TRAJ42*01_Bos_taurus_Hereford_F_J-REGION_3266515..3266580_66_nt_3_______21_AA_21+0=21_____
     2 IMGT000049_TRAJ40*01_Bos_taurus_Hereford_F_J-REGION_3268691..3268747_57_nt_3_______18_AA_18+0=18_____
     2 IMGT000049_TRAJ38*01_Bos_taurus_Hereford_F_J-REGION_3271437..3271498_62_nt_2_______20_AA_20+0=20_____
     2 IMGT000049_TRAJ35*01_Bos_taurus_Hereford_ORF_J-REGION_3275122..3275177_56_nt_2_______18_AA_18+0=18_____
     2 IMGT000049_TRAJ34*01_Bos_taurus_Hereford_F_J-REGION_3276151..3276207_57_nt_3_______18_AA_18+0=18_____
     2 IMGT000049_TRAJ33*01_Bos_taurus_Hereford_F_J-REGION_3276841..3276897_57_nt_3_______18_AA_18+0=18_____
     2 IMGT000049_TRAJ3*01_Bos_taurus_Hereford_F_J-REGION_3311301..3311359_59_nt_2_______19_AA_19+0=19_____
     2 IMGT000049_TRAJ29*01_Bos_taurus_Hereford_F_J-REGION_3282544..3282603_60_nt_3_______19_AA_19+0=19_____
     2 IMGT000049_TRAJ27*01_Bos_taurus_Hereford_F_J-REGION_3284184..3284242_59_nt_2_______19_AA_19+0=19_____
     2 IMGT000049_TRAJ25*01_Bos_taurus_Hereford_F_J-REGION_3286649..3286708_60_nt_3_______19_AA_19+0=19_____
     2 IMGT000049_TRAJ2*01_Bos_taurus_Hereford_F_J-REGION_3311840..3311905_66_nt_3_______21_AA_21+0=21_____
     2 IMGT000049_TRAJ19*01_Bos_taurus_Hereford_P_J-REGION_3292268..3292325_58_nt_1_______19_AA_19+0=19_____
     2 IMGT000049_TRAJ17*01_Bos_taurus_Hereford_F_J-REGION_3293762..3293824_63_nt_3_______20_AA_20+0=20_____
     2 IMGT000049_TRAJ12*01_Bos_taurus_Hereford_F_J-REGION_3299008..3299067_60_nt_3_______19_AA_19+0=19_____
     2 IMGT000049_TRAJ11*01_Bos_taurus_Hereford_F_J-REGION_3299569..3299628_60_nt_3_______19_AA_19+0=19_____
     2 D16120_TRGJ2-1*02_Bos_taurus_(F)_J-REGION_338..396_59_nt_3_______19_AA_19+0=19_partial_in_3'___
     2 D16118_TRGJ2-1*03_Bos_taurus_(F)_J-REGION_343..399_57_nt_1_______19_AA_19+0=19_partial_in_3'___
     2 D13648_TRGJ3-1*02_Bos_taurus_(F)_J-REGION_348..405_58_nt_2_______19_AA_19+0=19_partial_in_3'___
     2 AY644518_TRGJ1-1*01_Bos_taurus_F_J-REGION_102107..102166_60_nt_3_______19_AA_19+0=19_____
     2 AY644517_TRGC4*02_Bos_taurus_F_EX2A_g,167787..167845_60_nt_1_+1_-1___20_AA_20+0=20_____
     2 AY644517_TRGC3*02_Bos_taurus_F_EX2A_g,133153..133208_57_nt_1_+1_-1___19_AA_19+0=19_____
     2 AY227782_TRAJ31*02_Bos_taurus_F_J-REGION_17025..17082_58_nt_1_______19_AA_19+0=19_____
     2 AY227782_TRAJ25*02_Bos_taurus_F_J-REGION_24285..24344_60_nt_3_______19_AA_19+0=19_____
     2 AY149283_IGHJ1-2*02_Bos_taurus_ORF_J-REGION_783..838_56_nt_2_______18_AA_18+0=18_____
     2 AC172685_TRGJ2-1*01_Bos_taurus_F_J-REGION_236573..236632_60_nt_3_______19_AA_19+0=19_____
     1 KT723008_IGHA*01_Bos_taurus_Holstein_F_M_g,655709..655893_186_nt_1_+1_____62_AA_62+0=62_____
     1 AY227782_TRAJ11*02_Bos_taurus_F_J-REGION_37141..37200_60_nt_3_______19_AA_19+0=19_____
```


So It appears the blast file contains all the different regions of the genes so we really are seeing 9 genes

| gene | chromosome | start position in Bos|
| -- | -- | --|
|AY644518_TRGJ1| 01| 102107 |
|KT723008_IGHJ2|01| 509129 |
|AC172685_IGHA| 01| 655709|
|IMGT000049_TRDC|01| 3301039|
|D16120_TRGJ2|02| 338 |
|AY149283_IGHJ1| 02 | 783 |
|AY2277782_TRAJ31| 02| 1725 |
|AY644517_TRG| 02 |167787 |
|NW_001494075_IGHJ1|03| 39638 |


And for some reason 2 bos immune genes did not have a blast hit to the bos genome, which is strange.  

It is notable that the immune genes lost in Elk happen to occur on Elk chromosomes that underwent a large scale chromosomal rearrangment and fusion. This chromosomal assembly of the Elk genome will provide an excellent resource for further investigation of how these genes may be involved in the lack of adaptive cellular immune responses to existing Brucella vaccines.                              

## Literature

Looks like there is some interesting aspects to what is present in cattle based on this [paper: Immunoglobulin genes and diversity: what we have learned from domestic animals](https://jasbsci.biomedcentral.com/articles/10.1186/2049-1891-3-18 )

It makes me wonder if there immune genes that are only found in Elk and not in cow.

I was hoping for a result that would just jump out and be aha but looks like I would need to do a much deeper literature review.  Time to bring in the domain experts.
