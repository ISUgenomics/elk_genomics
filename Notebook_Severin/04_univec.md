# UniVec contamination Check

* /project/elk_bison_genomics/Severin/04_univec
* 2019 June 17

## UniVec sequence database of known vector contaminants

Here is a description taken directly from the README.uv file.

```UniVec is a non-redundant database of sequences commonly attached to
cDNA or genomic DNA during the cloning process.  It was developed by
staff at the National Center for Biotechnology Information, part of
the National Library of Medicine at the National Institutes of
Health. UniVec_Core is a subset of the full UniVec database. ```

The sequences for Univec are found on [NCBI ftp site](ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec)

## Make a Blast database for UniVec

```
module load blast+/2.6.0
makeblastdb -in UniVec -input_type fasta -dbtype nucl -parse_seqids -out UniVec

```


blastn  -db UniVecDB -query north_american_elk_15Jun2018_oY8t2.fasta  -num_threads 20 -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 > univec_blastOut


## Load a version of python 2 and install jcvi module locally

```
module load python_2/2.7.14
pip install --user jcvi
```

## Run the jcvi vecscreen of the univec datbase

```
module load blast/2.7.1
python -m jcvi.apps.vecscreen mask --db=UniVec north_american_elk_15Jun2018_oY8t2.fasta
```


## Univec output

There are 45 scaffolds that have the potential of contamination from one of the vectors in univec

```
more north_american_elk_15Jun2018_oY8t2.UniVec_Core.blast | awk   '{print $1}' | sort | uniq | wc
     45      45    1131
```

Sorting by lowest evalue we see that there are only 2 regions on scaffold ScoY8t2_23241;HRSCAF=23412 between 113529490 and 113529692 that may warrant further investigation. Though even longest alignment for these blast results are 123bp.

```
sort -gk 11 north_american_elk_15Jun2018_oY8t2.UniVec_Core.blast | head
ScoY8t2_23241;HRSCAF=23412      gnl|uv|KF947528.1:4787-5103     93.496  123     5       1       113529490       113529609       135     13      2.47e-44        185
ScoY8t2_23241;HRSCAF=23412      gnl|uv|AF416744.1:635-951       88.288  111     4       2       113529582       113529692       260     159     8.79e-28        130
ScoY8t2_23302;HRSCAF=23755      gnl|uv|JX069762.1:1-134 100.000 34      0       0       37766328        37766361        17      50      7.22e-09        67.7
ScoY8t2_23302;HRSCAF=23755      gnl|uv|JX069764.1:1-100 100.000 34      0       0       37766328        37766361        9       42      7.22e-09        67.7
ScoY8t2_23302;HRSCAF=23755      gnl|uv|JX069764.1:9741-9848-49  100.000 34      0       0       37766328        37766361        117     150     7.22e-09        67.7
ScoY8t2_23302;HRSCAF=23755      gnl|uv|NGB00363.1:1-34  100.000 34      0       0       37766328        37766361        34      1       7.22e-09        67.7
ScoY8t2_23302;HRSCAF=23755      gnl|uv|NGB00847.1:1-64  100.000 34      0       0       37766328        37766361        64      31      7.22e-09        67.7
ScoY8t2_23302;HRSCAF=23755      gnl|uv|NGB00848.1:1-64  100.000 34      0       0       37766328        37766361        64      31      7.22e-09        67.7
ScoY8t2_23302;HRSCAF=23755      gnl|uv|NGB00849.1:1-64  100.000 34      0       0       37766328        37766361        64      31      7.22e-09        67.7
ScoY8t2_23302;HRSCAF=23755      gnl|uv|NGB00850.1:1-64  100.000 34      0       0       37766328        37766361        64      31      7.22e-09        67.7
```

Only 2 scaffolds in this list had a length less than 1000 bp.

```
more *blast| sort -k 1,1 | awk '{print $1}' | sort | uniq | xargs -I xx awk '$1=="'xx'"' ../00_rawdata/north_american_elk_15Jun2018_oY8t2.len | awk '$2<1000'
ScoY8t2_20888;HRSCAF=21037 319
ScoY8t2_4955;HRSCAF=5094 784

```

```
grep "ScoY8t2_4955;HRSCAF=5094" *blast | sort -gk 11 | head -n 4
ScoY8t2_4955;HRSCAF=5094        gnl|uv|JX069762.1:1-134 100.000 28      0       0       1       28      44      17      2.54e-05        55.9
ScoY8t2_4955;HRSCAF=5094        gnl|uv|JX069762.1:9774-9902-49  100.000 28      0       0       1       28      173     146     2.54e-05        55.9
ScoY8t2_4955;HRSCAF=5094        gnl|uv|JX069764.1:1-100 100.000 28      0       0       1       28      36      9       2.54e-05        55.9
ScoY8t2_4955;HRSCAF=5094        gnl|uv|JX069764.1:9741-9848-49  100.000 28      0       0       1       28      144     117     2.54e-05        55.9

grep "ScoY8t2_20888;HRSCAF=21037" *blast | sort -gk 11 | head -n 4
ScoY8t2_20888;HRSCAF=21037      gnl|uv|JX069762.1:1-134 100.000 32      0       0       1       32      48      17      1.10e-07        63.8
ScoY8t2_20888;HRSCAF=21037      gnl|uv|JX069762.1:9774-9902-49  100.000 32      0       0       1       32      177     146     1.10e-07        63.8
ScoY8t2_20888;HRSCAF=21037      gnl|uv|JX069764.1:1-100 100.000 32      0       0       1       32      40      9       1.10e-07        63.8
ScoY8t2_20888;HRSCAF=21037      gnl|uv|JX069764.1:9741-9848-49  100.000 32      0       0       1       32      148     117     1.10e-07        63.8

```


### Note E coli file wasn't found or downloaded so it did not run contamination check against E coli genome.

```
Traceback (most recent call last):
  File "/software/7/apps/python_2/2.7.14/lib/python2.7/runpy.py", line 174, in _run_module_as_main
    "__main__", fname, loader, pkg_name)
  File "/software/7/apps/python_2/2.7.14/lib/python2.7/runpy.py", line 72, in _run_code
    exec code in run_globals
  File "/home/andrew.severin/.local/lib/python2.7/site-packages/jcvi/apps/vecscreen.py", line 108, in <module>
    main()
  File "/home/andrew.severin/.local/lib/python2.7/site-packages/jcvi/apps/vecscreen.py", line 25, in main
    p.dispatch(globals())
  File "/home/andrew.severin/.local/lib/python2.7/site-packages/jcvi/apps/base.py", line 96, in dispatch
    globals[action](sys.argv[2:])
  File "/home/andrew.severin/.local/lib/python2.7/site-packages/jcvi/apps/vecscreen.py", line 52, in mask
    assert op.exists(ecolifile)
AssertionError

[1]+  Exit 1                  python -m jcvi.apps.vecscreen mask --db=UniVec north_american_elk_15Jun2018_oY8t2.fasta
```
