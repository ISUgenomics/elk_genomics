#  Need to reduce gene models and identify those that are high confidence

```
#/work/gif/remkv6/Olsen/Elk/02_mergeMikadoBraker

#red deer modified transcript alignment to the masked genome
FinalGenomePilonReducedRenamedMasked.transcripts.gff3


ml cufflinks
gffread -K -M -Q -g ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa ../05_EvaluatePrediction/Fixedaugustus.hints.gff3 ../05_EvaluatePrediction/ExpressedIntactORFmikado.loci.gff3 -o gffmergeElkGenes.gff

#get proteins from this merge to run busco
 gffread gffmergeElkGenes.gff  -g ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa -VHEJ -t mRNA -x gffmergeElkGenesVHEJ_transcripts.fasta -y gffmergeElkGenesVHEJ_proteins.fasta


#merge the genes that -K -M -Q  with the red deer and bos taurus alignments
gffread -K -M -Q -g../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa gffmergeElkGenes.gff FinalGenomePilonReducedSoftMaskedRecode.Bos_taurus.ARS-UCD1.2.cds.all.gff3 FinalGenomePilonReducedRenamedMasked.transcripts.gff3 -o Bos_ReddeerReduction.gff

#input merged gene prediction
awk '$3=="locus"'  gffmergeElkGenes.gff|wc
  92965  836685 11113531

#total loci from merging elk to bos and red deer
awk '$3=="locus"'  Bos_ReddeerReduction.gff|wc
  96744  870696 13605644

#How many of the original merged braker/mikado set?
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -e "mRNA" -e  "mikado"  -e "MIKADO" |wc
  89762  807858 12573055
#~3k reduction


#How many loci are there that have either RED or BOS(ENS) in the locus tag?
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -e "RED" -e "ENS" |wc
  23829  214461 5592002

#How many have RED or BOS(ENS) also have a transcript in mikado or braker "mikado"  "MIKADO" and "mRNA"   
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -e "RED" -e "ENS" |grep -e "mikado" -e "mRNA" -e "MIKADO" |wc
  16847  151623 4559413

### How well did each prediction do?
  Braker was much better
#### How many BOS(ENS) or Elk(RED) loci were also in braker
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -e "RED" -e "ENS" |grep "mRNA" |wc
  14766  132894 4039738


#### How many BOS(ENS) or Elk(RED) loci were also in mikado
less Bos_ReddeerReduction.gff |awk '$3=="locus"' |grep -e "RED" -e "ENS" |grep -e "mikado"  -e "MIKADO" |wc
   8451   76059 2585545

# Get ready for busco4 evaluation of proteins
gffread Bos_ReddeerReduction.gff  -g ../05_EvaluatePrediction/FinalGenomePilonReducedSoftMaskedRecode.fa -VHEJ -t mRNA -x  Bos_ReddeerReductionVHEJ_transcripts.fasta -y  Bos_ReddeerReductionVHEJ_proteins.fasta



```


### Busco
