# Background

This file really is to sort out the data that we have and where it is located.

## File structure when starting of project

```
├── dovetail
│   ├── elk_dovetail_delivery_WIL462_KF0fHeNqkz.tar
│   ├── elk_dovetail_delivery_WIL462_SQND8mmmSa.tar
│   ├── elk_dovetail_local.md5sum
│   └── elk_dovetail.md5sum
├── elk_sequencing.md5sum
├── final.genome.scf.fasta.bz2
├── illumina
│   ├── 1_9130_01_WY1-Pine_HF557_1013.tar
│   ├── 1_9138_01_2510-Os-Mn_HH7YG_1052.tar
│   ├── 2_9131_01_WY2-Pine_HF557_1013.tar
│   ├── 2_9139_01_2486-OS-Mn_HG73W_1014.tar
│   ├── 3_9132_01_WY3-Pine_HF557_1013.tar
│   ├── 3_9140_01_2758-OS-Mn_HG73W_1014.tar
│   ├── 4_9133_01_WY7-Pine_HF557_1013.tar
│   ├── 4_9141_01_2763-OS-Mn_HG73W_1014.tar
│   ├── 5_9134_01_WY8-Pine_HF557_1013.tar
│   ├── 5_9142_01_8Bull-OSMn_HG73W_1014.tar
│   ├── 6_9135_01_WY10-Pine_HF557_1013.tar
│   ├── 6_9143_01_8Bull-OSMn_HG73W_1014.tar
│   ├── 7_9136_01_2458-0S-Mn_HF557_1013.tar
│   ├── 7_9144_01_8Bull-OSMn_HG73W_1014.tar
│   ├── 8_9137_01_2450-OS-Mn_HF557_1013.tar
│   ├── 8_9145_01_8Bull-OSMn_HG73W_1014.tar
│   ├── download_161102_bz2.md5sum
│   └── download_170119_bz2.md5sum
├── mRNA
│   ├── 3_4_22061_Elk-muscle_HVW3J_1422.tar
│   └── elk_mRNA.md5sum
└── pacbio
    ├── 10_9649_Bull-Elk1_PacBio_1060_1.tar.bz2
    ├── 10_9649_Bull-Elk1_PacBio_1060.tar.bz2
    ├── 10_9650_Bull-Elk2_PacBio_1061_1.tar.bz2
    ├── 10_9650_Bull-Elk2_PacBio_1061_2.tar.bz2
    ├── 10_9650_Bull-Elk2_PacBio_1061.tar.bz2
    ├── 10_9651_Bull-Elk3_PacBio_1067.tar.bz2
    ├── 10_9651_Bull-Elk3_PacBio_1078.tar.bz2
    ├── 10_9652_Bull-Elk4_PacBio_1079.tar.bz2
    └── elk_pacbio_bz2.md5sum

```

## md5sum

* /home/andrew.severin/elk_bison_genomics/dovetail
* 2019 June 12

I am going to assume that the md5sums that Darrell generated were to confirm they were completely downloaded but I will need to verify this.  The one I did below matches.
```
md5sum *.tar > elk_dovetail_local.md5sum &
```


## What is inside the tar files


```
tar -tf elk_dovetail_delivery_WIL462_KF0fHeNqkz.tar
delivery_WIL462_KF0fHeNqkz/
delivery_WIL462_KF0fHeNqkz/lib_001.sorted.md.bam
delivery_WIL462_KF0fHeNqkz/north_american_elk_11Jun2018_1UW0F.fasta.gz
delivery_WIL462_KF0fHeNqkz/north_american_elk_11Jun2018_1UW0F.input_breaks.txt
delivery_WIL462_KF0fHeNqkz/north_american_elk_11Jun2018_1UW0F.table.txt
delivery_WIL462_KF0fHeNqkz/lib_002.sorted.md.bam
delivery_WIL462_KF0fHeNqkz/north_american_elk_11Jun2018_1UW0F.report.pdf
delivery_WIL462_KF0fHeNqkz/md5sums.txt
delivery_WIL462_KF0fHeNqkz/delivery_WIL462_KF0fHeNqkz.sh
delivery_WIL462_KF0fHeNqkz/manifest.txt
delivery_WIL462_KF0fHeNqkz/north_american_elk_11Jun2018_1UW0F.gapclosing_table.txt
delivery_WIL462_KF0fHeNqkz/lib_003.sorted.md.bam

tar -tf elk_dovetail_delivery_WIL462_SQND8mmmSa.tar
delivery_WIL462_SQND8mmmSa/
delivery_WIL462_SQND8mmmSa/north_american_elk_15Jun2018_oY8t2.table.txt
delivery_WIL462_SQND8mmmSa/lib_001.sorted.md.bam
delivery_WIL462_SQND8mmmSa/delivery_WIL462_SQND8mmmSa.sh
delivery_WIL462_SQND8mmmSa/lib_002.sorted.md.bam
delivery_WIL462_SQND8mmmSa/north_american_elk_15Jun2018_oY8t2.gapclosing_table.txt
delivery_WIL462_SQND8mmmSa/md5sums.txt
delivery_WIL462_SQND8mmmSa/manifest.txt
delivery_WIL462_SQND8mmmSa/north_american_elk_15Jun2018_oY8t2.report.pdf
delivery_WIL462_SQND8mmmSa/north_american_elk_15Jun2018_oY8t2.input_breaks.txt
delivery_WIL462_SQND8mmmSa/lib_003.sorted.md.bam
delivery_WIL462_SQND8mmmSa/north_american_elk_15Jun2018_oY8t2.fasta.gz
```


## Note from Darrell

The Dovetail Genomics assembly improvement was performed in two separate assembly phases.  Each phase used a different type of proximity ligation data named Chicago library reads (three libraries) and Dovetail HiC reads (three libraries), respectively.  Both of these newly generated sequencing read types served as sequential inputs into the Dovetail Genomics' HiRise assembly pipeline that is specifically designed to use proximity ligation data.

Two helpful metrics for evaluating the overall state of an assembly are the L90 and N90 metrics.  The L90 is the smallest number of scaffolds that comprise 90% of the total assembly.  The N90 is the is the size of the smallest scaffold such that summing the lengths of all scaffolds that size and larger equal 90% of the total assembly.  

ARS supplied Dovetail Genomics with the current assembly that consisted of a total length of 2,559.8 Mb and an L90 of 2,500 scaffolds with and N90 of 0.198 Mb.  The first round of Dovetail assembly used the HiRise pipeline and the three Chicago libraries to yield a new assembly of 2,560.4 Mb, a L90 of 175 scaffolds, and a N50 of 2.292 Mb.  This first-round Dovetail assembly served as the input into the second round of assembly.  The second round of assembly used the HiRise pipeline and the three HiC libraries.  The second round of assembly yielded a final assembly of 2,560.5 Mb with a L90 of 31 scaffolds, and a N90 of 43.374 Mb.  Additionally, the longest scaffold in the original assembly was 7,681,499 bp in length, but after the two rounds of assembly improvement, the longest scaffold was 128,024,478 bp long.
