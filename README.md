# Population genomics of the threatened New Zealand storm petrel
This repository contains details for the data processing and analyses used to investigate the population structure and genetic diversity of the New Zealand storm petrel (NZSP).
### Contents:
* [DNA sequence processing with UNEAK](#dna-sequence-processing-with-uneak)
* [DNA sequence processing with Stacks](#dna-sequence-processing-with-stacks)
* [Population structure analysis](#population-structure-analysis)
# DNA sequence processing with UNEAK
## Filtering with KGD 
All KGD R scripts are provided in the [KGD R package](https://github.com/AgResearch/KGD). Some modifications were made to the R scripts, as outlined below:
### run_kgd.R
Added an argument to writeVCF() to only include SNPs with Hardy_weinberg disequilibrium > -0.05
```r
#Original:
writeVCF(outname="GHW05", ep=.001)
#Modified:
writeVCF(outname="GHW05X", ep=.001, snpsubset=which(HWdis.sep > -0.05))
```
### GBS-Chip-Gmatrix.R
Changed relatedness threshold so that all relatedness values would be included in HighRelatednessHWdgm.csv, rather than only the top most related pairs
```r
#Original:
uhirel <- which(GGBS5 > hirel.thresh & upper.tri(GGBS5), arr.ind = TRUE)
#Modified:
uhirel <- which(GGBS5 > -5 & upper.tri(GGBS5), arr.ind = TRUE)
```
## Allelic depth filter
Allelic depth filter of 8 was applied to the VCF (GHW05.vcf) produced by KGD
```
module load BCFtools
bgzip GHW05.vcf
bcftools index GHW05.vcf.gz 
bcftools filter -S . -i 'FMT/AD[*:*]>7' -O v -o GHW05depth8.vcf GHW05full.vcf.gz #Genotypes are only included if at least one of the AD fields is greater than 7 (depth of 8)
```
## Filtering with Plink in R
Files were filtered using Plink commands in R
### Preparing VCF file for Plink
Currently, each SNP is given its own chromosome. This will clean up the sample names and make all SNPs on chromosome 1.
```
sed -i 's/_merged_2_0_X4//g' GHW05.vcf  #Tidy up the sample names
head -n 15 GHW05.vcf > GHW05head.vcf #Take just the header
awk 'BEGIN {OFS="\t"}; NR>22 {$1=1; print}' GHW05.vcf > GHW05edit.vcf #Replace all first column chrosome values with 1 and add tab spaces, also removes header
cat GHW05head.vcf GHW05edit.vcf > GHW05full.vcf #append files together to create the full, cleaned up VCF file
```
### Applying SNP and sample missingness filters
```r
#Set working directory in R and PATH in console
shell("plink")
shell("plink  --allow-no-sex  --out NZSP  --recode  --vcf GHW05full.vcf") #load VCF file into Plink

#Make dataset for population structure analysis
shell("plink --file NZSP --out NZSP25 --geno 0.25 --make-bed") #Filter for maximum threshold of 25% missingness per variant
shell("plink --bfile NZSP25 --out NZSP25_50 --mind 0.5 --make-bed") #Filter for maximum threshold of 50% missingness per sample
shell("plink --bfile NZSP25_50 --rel-cut 0.3 --out NZSP25_50_rel03") #Remove close relatives with threshold of 0.3

#Make dataset for parantage and effective population size analysis
shell("plink --file NZSP --out NZSP25 --geno 0.15 --make-bed") #Filter maximum threshold of 15% missingness per variant
shell("plink --bfile NZSP25 --out NZSP15_50 --mind 0.5 --make-bed") #Filter maximum threshold of 50% missingness per sample
```
# DNA sequence processing with Stacks
# Population structure analysis
