### Need some of the files generated by blobtools.Run Blobtools

```
cd 15_Blobtools3/

ln -s ../11_PolishGenomeWithCCS/FinalAssemblyFastaWithY.FinalAssemblyFastaWithY_sorted.bam
ln -s ../11_PolishGenomeWithCCS/FinalAssemblyFastaWithY.FinalAssemblyFastaWithY_sorted.bam.bai
ln -s ../11_PolishGenomeWithCCS/FinalAssemblyFastaWithY.fasta
ln -s ../12_C.elaphusHippelaphus2/AllBlasts.out
sed 's|N/A|9860|g' AllBlasts.out >TaxonedAllBlasts.out


```
