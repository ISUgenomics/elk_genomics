### Testing juicer pipeline on ceres with small genome

## Use SCN as example 129mb genome size
```
/home/rick.masonbrink/elk_bison_genomics/Masonbrink/04_JuicerElk
ln -s ../north_american_elk_15Jun2018_oY8t2.fasta

mkd
mkdir fastq


bioawk -c fastx '{print $name, length($seq)}'
```
