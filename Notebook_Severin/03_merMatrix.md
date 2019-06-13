# Kmer Matrix

One way to look for contamination is to look at the frequency of 3-mers in each scaffold and look for scaffolds that cluster differently from the rest.  This is a function I requested of the jellyfish author but have yet to find a good place for publishing it. The program can be found in the scripts folder in this repo.


## Create a softlink to the mer_matrix program
```
ln -s ~/isugif/privateGIF/mer_matrix-0.0.1/bin/mer_matrix
```

## Run Mer matrix to create a 64 column matrix of 3mers

```
module load jellyfish2/2.2.9
#was getting an error that it couldn't find libjellyfish-2.0.so.2 so I added it to the LD_LIBRARY_PATH variable
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/software/7/apps/jellyfish2/2.2.9/lib/"
./mer_matrix ../00_rawdata/north_american_elk_15Jun2018_oY8t2.fasta >elk.3mer

```


## Explore output in R

Transferred the output to my laptop

* /Users/severin/Desktop/Projects/elk_genomics/03_merMatrix

#### Load the 3mer data and calculate PCA of vectors
```
elkmer<-read.table('elk.3mer')
rownames(elkmer)<-elkmer[,1]
elkmer<-elkmer[,2:dim(elkmer)[2]]
elkmer.rs<-rowSums(elkmer)
elkmer.norm<-10*sweep(elkmer,1,elkmer.rs,FUN="/")

elkpca<-prcomp((elkmer.norm))

```

#### Load length data and normalize

```
elklen<-read.table('north_american_elk_15Jun2018_oY8t2.len')

```


```
library(ggplot2)

# Here I am creating a new x11 window then plotting the elkpca x values for PC1 and PC6 then changing the point size from 2 to the log of the scaffold size with an open circle shape.  I had to sort the elklen in the same order as the elkpca data as well.
x11()
ggplot(as.data.frame(elkpca$x), aes(PC1, PC4),size=2)+ geom_point(aes(size = log(elklen[rownames(elkpca$x),])),shape=1)



```

#### TODO

* See if I can color the top 37 scaffolds a different color in this plot to highlight the glob of points that definitely correspond to the Elk genome.
* blast some of the scaffolds in different clusters to see if they are contamination quick check, can do full check with blobtools
