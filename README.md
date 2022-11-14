# Population genomics of the threatened New Zealand storm petrel
This repository contains details for the data processing and analyses used to investigate the population structure and genetic diversity of the New Zealand storm petrel (NZSP).
### Contents:
* [DNA sequence processing with UNEAK](#dna-sequence-processing-with-uneak)
* [DNA sequence processing with Stacks](#dna-sequence-processing-with-stacks)
* [Population structure analysis](#population-structure-analysis)
* [Minor allele frequency distribution plots](#minor-allele-frequency-distribution-plots)
* [Kinship analysis](#kinship-analysis)
* [Effective population size](#effective-population-size)
# DNA sequence processing with UNEAK
## Filtering with KGD 
All KGD R scripts are provided in the [KGD R package](https://github.com/AgResearch/KGD). Some modifications were made to lines in the R scripts, as outlined below:
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

#Make BED file dataset for population structure analysis
shell("plink --file NZSP --out NZSP25 --geno 0.25 --make-bed") #Filter for maximum threshold of 25% missingness per variant
shell("plink --bfile NZSP25 --out NZSP25_50 --mind 0.5 --make-bed") #Filter for maximum threshold of 50% missingness per sample
shell("plink --bfile NZSP25_50 --rel-cut 0.3 --out NZSP25_50_rel03") #Remove close relatives with threshold of 0.3

#Make BED file dataset for parentage and effective population size analysis
shell("plink --file NZSP --out NZSP25 --geno 0.15 --make-bed") #Filter maximum threshold of 15% missingness per variant
shell("plink --bfile NZSP25 --out NZSP15_50 --mind 0.5 --make-bed") #Filter maximum threshold of 50% missingness per sample
```
# DNA sequence processing with Stacks
## Adapter trimming
Trim adapters off of both raw DNA sequence files using cutadapt
```
module purge
module load cutadapt
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a CGAGATCGGAAGAGCGGACTTTAAGC -o NZSP_1_fastq.gz SQ1792_HN7WGDRXY_s_1_fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a CGAGATCGGAAGAGCGGACTTTAAGC -o NZSP_2_fastq.gz SQ1792_HN7WGDRXY_s_2_fastq.gz
```
Check the trim worked using FastQC quality checking
```
module purge
module load FastQC
fastqc NZSP_1_fastq.gz 
fastqc NZSP_2_fastq.gz 
```
## Assembly with Stacks
Demulitplex and clean the DNA sequences with process_radtags
```
cat NZSP_1_fastq.gz NZSP_2_fastq.gz > NZSP_fastq.gz #Combine both files of DNA sequences

module purge
module load Stacks
process_radtags -1 ./NZSP_fastq.gz -2 ./NZSP_2_fastq.gz --renz-1 pstI --renz-2 mspI -c -q -b NZSP_Barcodes.txt -t 65 -o ./radtags_out2 #t= truncate to length 65
```
### Kmer filtering with Stacks
1. Use R to print out lines of code to be used with Stacks for each sample
```r
#Print out lines of code in R for each sample
df<-as.data.frame(All_sheep) #Create dataframe of each sample name from barcode file
ls<-df[!duplicated(df$V1), ] #Remove duplicate samples
ls<-as.list(ls)

for(i in ls) {
    print(paste("./kmer_filter -1 ./radtags_out/",i,".1.fq.gz -2 ./radtags_out/",i,".2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant", sep=""), quote = FALSE)
}
```
2. Use printed lines of code to run a kmer filter on each sample
<details><summary>Kmer filter code</summary>
  <p>
    
```
module load Stacks
kmer_filter -1 ./radtags_out/K64513.1.fq.gz -2 ./radtags_out/K64513.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64514.1.fq.gz -2 ./radtags_out/K64514.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64515.1.fq.gz -2 ./radtags_out/K64515.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64516.1.fq.gz -2 ./radtags_out/K64516.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64517.1.fq.gz -2 ./radtags_out/K64517.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64518.1.fq.gz -2 ./radtags_out/K64518.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64519.1.fq.gz -2 ./radtags_out/K64519.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64520.1.fq.gz -2 ./radtags_out/K64520.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64521.1.fq.gz -2 ./radtags_out/K64521.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64522.1.fq.gz -2 ./radtags_out/K64522.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64523.1.fq.gz -2 ./radtags_out/K64523.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64524.1.fq.gz -2 ./radtags_out/K64524.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64525.1.fq.gz -2 ./radtags_out/K64525.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64526.1.fq.gz -2 ./radtags_out/K64526.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64527.1.fq.gz -2 ./radtags_out/K64527.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64528.1.fq.gz -2 ./radtags_out/K64528.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64529.1.fq.gz -2 ./radtags_out/K64529.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64530.1.fq.gz -2 ./radtags_out/K64530.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64531.1.fq.gz -2 ./radtags_out/K64531.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64532.1.fq.gz -2 ./radtags_out/K64532.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64533.1.fq.gz -2 ./radtags_out/K64533.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64534.1.fq.gz -2 ./radtags_out/K64534.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64535.1.fq.gz -2 ./radtags_out/K64535.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64536.1.fq.gz -2 ./radtags_out/K64536.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64537.1.fq.gz -2 ./radtags_out/K64537.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64538.1.fq.gz -2 ./radtags_out/K64538.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64539.1.fq.gz -2 ./radtags_out/K64539.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64540.1.fq.gz -2 ./radtags_out/K64540.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64541.1.fq.gz -2 ./radtags_out/K64541.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64542.1.fq.gz -2 ./radtags_out/K64542.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64543.1.fq.gz -2 ./radtags_out/K64543.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64544.1.fq.gz -2 ./radtags_out/K64544.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64545.1.fq.gz -2 ./radtags_out/K64545.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64546.1.fq.gz -2 ./radtags_out/K64546.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64547.1.fq.gz -2 ./radtags_out/K64547.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64548.1.fq.gz -2 ./radtags_out/K64548.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64549.1.fq.gz -2 ./radtags_out/K64549.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64550.1.fq.gz -2 ./radtags_out/K64550.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64551.1.fq.gz -2 ./radtags_out/K64551.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64552.1.fq.gz -2 ./radtags_out/K64552.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64553.1.fq.gz -2 ./radtags_out/K64553.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64554.1.fq.gz -2 ./radtags_out/K64554.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64555.1.fq.gz -2 ./radtags_out/K64555.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64556.1.fq.gz -2 ./radtags_out/K64556.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64557.1.fq.gz -2 ./radtags_out/K64557.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64558.1.fq.gz -2 ./radtags_out/K64558.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64559.1.fq.gz -2 ./radtags_out/K64559.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64560.1.fq.gz -2 ./radtags_out/K64560.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64561.1.fq.gz -2 ./radtags_out/K64561.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64562.1.fq.gz -2 ./radtags_out/K64562.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64563.1.fq.gz -2 ./radtags_out/K64563.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64564.1.fq.gz -2 ./radtags_out/K64564.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64565.1.fq.gz -2 ./radtags_out/K64565.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64566.1.fq.gz -2 ./radtags_out/K64566.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64567.1.fq.gz -2 ./radtags_out/K64567.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64568.1.fq.gz -2 ./radtags_out/K64568.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64569.1.fq.gz -2 ./radtags_out/K64569.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64570.1.fq.gz -2 ./radtags_out/K64570.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64571.1.fq.gz -2 ./radtags_out/K64571.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64572.1.fq.gz -2 ./radtags_out/K64572.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64573.1.fq.gz -2 ./radtags_out/K64573.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64574.1.fq.gz -2 ./radtags_out/K64574.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64575.1.fq.gz -2 ./radtags_out/K64575.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64576.1.fq.gz -2 ./radtags_out/K64576.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64577.1.fq.gz -2 ./radtags_out/K64577.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64578.1.fq.gz -2 ./radtags_out/K64578.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64579.1.fq.gz -2 ./radtags_out/K64579.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64580.1.fq.gz -2 ./radtags_out/K64580.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64581.1.fq.gz -2 ./radtags_out/K64581.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64582.1.fq.gz -2 ./radtags_out/K64582.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64583.1.fq.gz -2 ./radtags_out/K64583.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64584.1.fq.gz -2 ./radtags_out/K64584.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64585.1.fq.gz -2 ./radtags_out/K64585.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64586.1.fq.gz -2 ./radtags_out/K64586.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64587.1.fq.gz -2 ./radtags_out/K64587.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64588.1.fq.gz -2 ./radtags_out/K64588.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64589.1.fq.gz -2 ./radtags_out/K64589.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64590.1.fq.gz -2 ./radtags_out/K64590.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64591.1.fq.gz -2 ./radtags_out/K64591.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64592.1.fq.gz -2 ./radtags_out/K64592.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64593.1.fq.gz -2 ./radtags_out/K64593.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64594.1.fq.gz -2 ./radtags_out/K64594.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64595.1.fq.gz -2 ./radtags_out/K64595.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64596.1.fq.gz -2 ./radtags_out/K64596.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64597.1.fq.gz -2 ./radtags_out/K64597.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64598.1.fq.gz -2 ./radtags_out/K64598.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64599.1.fq.gz -2 ./radtags_out/K64599.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64600.1.fq.gz -2 ./radtags_out/K64600.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64601.1.fq.gz -2 ./radtags_out/K64601.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64602.1.fq.gz -2 ./radtags_out/K64602.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64603.1.fq.gz -2 ./radtags_out/K64603.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64604.1.fq.gz -2 ./radtags_out/K64604.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
kmer_filter -1 ./radtags_out/K64605.1.fq.gz -2 ./radtags_out/K64605.2.fq.gz -i gzfastq -o ./kmerfil_out/ --rare --abundant
```
</p>
</details>
 
3. Clean up file names for next steps
```
for filename in *.fq; do      [ -f "$filename" ] || continue;     mv "$filename" "${filename//fil./}";  done #remove fil. from file names
```
4. Optional quality check of sequences 
```
module purge
module load FastQC
fastqc *1.fq -o fastqc
fastqc *2.fq -o fastqc

module purge 
module load MultiQC
multiqc . #Summarises all FastQC reports for each sample into single report
```
5. Run denovo Stacks wrapper for de novo assembly of sequences. Samples are listed and categorised into populations depending on location of origin in the popmap file
```
module purge 
module load Stacks
denovo_map.pl -T 8 -M 2 -o ./denovo_M2/ --samples ./kmerfil_out3 --popmap ./Stacks_pop.txt 
```
6. Make blacklist of loci with two or more SNPs. The number of SNPs for each loci was counted and then counted loci were filtered in excel to create a blacklist of loci with two or more SNPs
```
awk -F"\t" '{print $1}' Test.txt | uniq -c > Newtest.txt #Counts duplicate loci in column 1
```
7. Rerun last step of denovo_map wrapper (populations) with blacklist. 

Blacklist contains all loci with two or more SNPs. A duplicate sample and three samples of unknown origin and date were also removed from the sample list, and an argument was added to create a VCF of the final Stacks dataset
```
module purge 
module load Stacks
populations -P ./denovo_M2 -t 8 -M ./stacks_popF.txt -O ./pop_M2F_S3 -B Blacklist_2.txt --vcf --plink #Create population statistics and VCF of SNPs for specified samples and loci
```
**KGD filtering was applied using modifications descibed above for [DNA sequence processing with UNEAK](#dna-sequence-processing-with-uneak)**

One additional modification was made to run_KGD.R to change the format of the input file to a VCF rather than a UNEAK dataset:
```r
#Original:
gform <- "uneak"
#Modified:
gform <- "VCF"
```
**Allelic depth and SNP/sample missingness filters were also applied as described above for the [DNA sequence processing with UNEAK](#dna-sequence-processing-with-uneak) filtering steps**
# Population structure analysis
## PCA
### Conduct PCA in Plink for the UNEAK and Stacks datasets
Conduct PCA using Plink in R with the UNEAK and Stacks BED files to create eigenvector and eigenvalue files. An additional PCA plot was created using the UNEAK dataset with a relatedness cutoff of 0.3.
```r
#Set working directory in R and PATH in console
shell("plink")
shell("plink  --allow-no-sex --make-rel  --bfile NZSP25_50_rel03 --out NZSP25_50_rel03  --recode --pca 100")
shell("plink  --allow-no-sex --make-rel  --bfile Stack10_50 --out Stacks_10_50  --recode --pca 100")
```
### Create PCA plots in R for both datasets
```r
file<-"NZSP25_50_rel03" #UNEAK or Stacks dataset name
in1<-read.table(paste("C:\\path\\to\\file\\",file,".eigenvec", sep=""))
plot(in1$V3,in1$V4)

#Requires sample key file with ID column and location column for each sample
names<-read.table("C:\\path\\to\\file\\samples_25_rel03.txt")
head(names)
colnames(in1)[2] <- "ID"
colnames(names)[1] <- "ID"
allt<-merge(in1,names)
length(allt[,1]) #86
head(allt)

Eigenvalue<-scan(paste("C:\\path\\to\\file\\",file,".eigenval", sep=""))
Eigenvalue
pve<-data.frame(PC = 1:100, pve = Eigenvalue/sum(Eigenvalue)*100)
head(pve)

#Define colours 
G_colour<-"#E1BE6A" #Gold
N_colour<-"#62A4DC" #Blue

#Extract Far North samples
allt_North<-subset(allt, V2 =="Far_North") 
length(allt_North[,1]) #27

#Plot just the Far North samples with increased point size (cex) and in triangle shape (pch)
plot(allt_North$V3,allt_North$V4,pch=17,xlab = (paste0("PC1 (", signif(pve$pve[1], digits = 3), "%)")), ylab = (paste0("PC2 (", signif(pve$pve[2], digits = 3), "%)")), col=N_colour, ylim = range(allt$V4), xlim = range(allt$V3), cex=1.2)
legend("topright", c("Far North", "Hauraki Gulf"), 
       col=c(N_colour, G_colour), pch=c(17, 19))

#Plot Hauraki samples as well with circle shape and increased size
allt_Hauturu<-subset(allt, V2 =="Hauraki_Gulf")
length(allt_Hauturu[,1]) 
points(allt_Hauturu$V3,allt_Hauturu$V4,col=G_colour, pch=19, cex=1.2)

#Optional - add labels to plot points with the last three digits of the sample ID
#Uses only the last three digits; 4 to 6 (identifier)
#pos = position of label; 1 under, 2 left, 3 top and 4 right
text(allt$V3, allt$V4, labels = substring(names$ID, 4, 6), cex = 0.8, pos = 1)
```
## FastStructure
Conduct fastStructure analysis on UNEAK and Stacks BED files produced by Plink filtering
```
#Move filtered BED files (bed, bim and fam) file into directory
module purge
module load fastStructure
#Test K=1 to 5
for k in {1..10}
do
  structure.py -K $k --input=NZSP8_g1_m5 --output=NZSP1_structure
done

chooseK.py --input=NZSP1_structure #Identify the best K 
distruct.py -K 1 --input=NZSP1_structure  --popfile=Structure_stacks_g1.txt --output=NZSPg1_distruct1.svg # Make fastStructure plot with the best K
```
# Minor allele frequency distribution plots
Create minor allele frequency (MAF) distribution plots for UNEAK and Stacks datasets

Packages: 
* ggplot2
```r
#load in Gulf and North sample IDs for the dataset (UNEAK or Stacks)
Gulf_Samples<-read.table("C:\\path\\to\\file\\Gulf_Samples.txt")
North_Samples<-read.table("C:\\path\\to\\file\\North_Samples.txt")

#Generate 100 random subsamples of 20 samples for each location 
for(i in 1:100) {
  a<-as.vector(Gulf_Samples$V1)
  V1<-a[sample(1:length(a), 20)] #Select 20 samples randomly
  Subsamp<-as.data.frame(V1)
  Subsamp$V2<-Subsamp$V1 #Make ID and family ID column (identical) for use in plink
  write.table(Subsamp, file=paste("C:\\path\\to\\directory\\Subsamp_G_",i,".txt", sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote = FALSE)
}
for(i in 1:100) {
  a<-as.vector(North_Samples$V1)
  V1<-a[sample(1:length(a), 20)] #Select 20 samples randomly 
  Subsamp<-as.data.frame(V1)
  Subsamp$V2<-Subsamp$V1 #Make ID and family ID column (identical) for use in plink
  write.table(Subsamp, file=paste("C:\\path\\to\\directory\\Subsamp_N_",i,".txt", sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote = FALSE)
}

#Set working directory in R and PATH in console 
#Move BED files into working directory
shell("cd")
shell("plink")
#Generate MAF statistics for the 100 subsamples using the filtered UNEAK or Stacks BED files 
for (i in 1:100) {
  shell(paste("plink --bfile NZSP25_50 --out NZSP_G",i," --keep Subsamp_G_",i,".txt --make-bed", sep="")) #Keep only the Gulf subsample, remove all other samples
  shell(paste("plink --bfile NZSP_G",i," --freq --out NZSP_G",i, sep="")) #Generate MAFs stats
}
for (i in 1:100) {
  shell(paste("plink --bfile NZSP25_50 --out NZSP_N",i," --keep Subsamp_N_",i,".txt --make-bed", sep="")) #Keep only the North subsample, remove all other samples
  shell(paste("plink --bfile NZSP_N",i," --freq --out NZSP_N",i, sep="")) #Generate MAFs stats
}

#Generate histogram counts for each subsample
CountsN<-data.frame() #create empty data frame
CountsG<-data.frame() #create empty data frame

for(i in 1:100) {
  MAF<-read.table(paste("C:\\path\\to\\directory\\NZSP_N",i,".frq", sep=""), header = TRUE)
  hist_info <- hist(MAF$MAF)
  CountsN<-rbind(CountsN, hist_info$counts)
}
colnames(CountsN)<-c("1","2","3","4","5","6","7","8","9","10")

for(i in 1:100) {
  MAF<-read.table(paste("C:\\path\\to\\directory\\NZSP_G",i,".frq", sep=""), header = TRUE)
  hist_info <- hist(MAF$MAF)
  CountsG<-rbind(CountsG, hist_info$counts)
}
colnames(CountsG)<-c("1","2","3","4","5","6","7","8","9","10")

#Calculate standard deviation and mean
Sd<- sapply(CountsG, FUN = sd) #calculate standard deviation for all columns 
sd<- as.data.frame(Sd)
mean<-sapply(CountsG,FUN=mean) #calculate the mean for all columns
mean<-as.data.frame(mean)

range<-c("0-0.05", "0.05-0.10", "0.10-0.15", "0.15-0.20", "0.20-0.25", "0.25-0.30", "0.30-0.35", "0.35-0.40","0.40-0.45", "0.45-0.50")
data<-cbind(mean, range)
data<-as.data.frame(data)
data$mean<-as.numeric(data$mean)

#Create histogram plots for each location using mean and standard deviation
library(ggplot2)
ggplot(data) + 
  geom_bar(aes(x=range, y=mean), fill="#E1BE6A", stat="identity", width=0.97) + 
  geom_errorbar( aes(x=range, ymin=mean-sd$Sd, ymax=mean+sd$Sd), width=0.2, size=1) + 
  labs(x="MAF", y="Frequency") + ggtitle("Hauraki Gulf")

ggplot(data) + 
  geom_bar(aes(x=range, y=mean), fill="#62A4DC", stat="identity", width=0.97) +
  geom_errorbar( aes(x=range, ymin=mean-sd$Sd, ymax=mean+sd$Sd), width=0.2, size=1) +
  labs(x="MAF", y="Frequency") + ggtitle("Far North")
```
# Kinship analysis
## Sequoia
Kinship assignment with Sequoia for the UNEAK and Stacks smaller datasets - sequoia R package only works using R v4.1.1

Packages:
* plyr
* sequoia v2.3.5 (removed from CRAN repository; archived package)
```r
install.packages("plyr")
install.packages("C:/path/to/Downloads/sequoia_2.3.5.tar.gz", repos = NULL, type = "source") 
library(sequoia)

#Set working directory in R and set PATH in console
#Move UNEAK or Stacks BED files
shell("cd")
shell("plink")
shell("plink --bfile NZSP15_50 --recodeA --out Sequoia_25_50 ") #Make .raw file using Plink and Stacks or UNEAK BED files 
GenoM <- GenoConvert(InFile = "C:\\path\\to\\file\\Sequoia_15_50.raw", InFormat = "raw")
GetMaybeRel(GenoM = GenoM)
```
## Plink and KGD relatedness estimates
The UNEAK and Stacks filtered BED files were also used to generate relatedness estimates for each pairing of individuals in Plink
```r
shell("plink --make-rel square --bfile NZSP25_50 --make-bed --out NZSP25_50_square") #Produce relatedness estimates in square matrix

#Put the Plink relatedness estimates in a similar format as KGD relatedness estimates:
Square <- read.table("C:\\path\\to\\file\\NZSP25_50_square.rel", strip.white=TRUE) #Plink relatedness estimates
RelID <- read.table("C:\\path\\to\\file\\NZSP25_50_square.rel.id", strip.white=TRUE) #Plink relatedness IDs)
Plink_vector <- as.vector(Square) 
len<-length(RelID$V1)
Col2 <- rep(RelID$V1, len) #Replicate each sample in order by length of ID list
head (Col2)
Col1 <- rep(RelID$V1, each=len) #Repeat each sample individually 
head(Col1)
Plinkrel<-cbind(Col1, Col2) 
Plinkrel<-data.frame(Plinkrel)
Plinkmerge<-paste(Plinkrel$Col1, Plinkrel$Col2, sep="_") #Merge sample pairs into one column, seperated by an underscore
Plink<-cbind(Plinkmerge, Plink_vector) #Merge with relatedness estimates
Plink <- as.data.frame(Plink)
```
The relatedness estimates from Plink and KGD were also plotted for comparison. KGD pairings were organised in Excel (Microsoft). All self-pairs in KGD were assigned to a relatedness of 1 and subsequently removed from the plot.
```r
KGDrel<-read.table("C:\\path\\to\\file\\Relatedness_KGD.txt", header=TRUE)
KGDrel<-KGDrel[!(KGDrel$KGDrel=="1"),] #Remove self pairs
KGDrel<-as.data.frame(KGDrel)
KGDmerge<-paste(KGDrel$Indiv1, KGDrel$Indiv2, sep="_")
KGD<-cbind(KGDmerge, KGDrel$KGDrel)
colnames(Plink)[1] <- "Sample Pair"
colnames(Plink)[2] <- "Plinkrel"
colnames(KGD)[1] <- "Sample Pair"
colnames(KGD)[2] <- "KGDrel"
KGD <- as.data.frame(KGD)
Plink <- as.data.frame(Plink)
KGD_Plink<-merge(KGD, Plink, by="Sample Pair")
colnames(KGD_Plink)[2]<- "KGDrel"

KGD_Plink$KGDrel<-as.numeric(KGD_Plink$KGDrel)
KGD_Plink$Plinkrel<-as.numeric(KGD_Plink$Plinkrel)

plot(KGD_Plink$KGDrel, KGD_Plink$Plinkrel,xlab = "KGD rel", ylab = "Plink rel", col="cornflowerblue")
```
## Sample homozygosity
Expected and observed homozygosity values for the UNEAK dataset were plotted to investigate inbreeding  
```r
#Set working directory in R and PATH in console
shell("plink --bfile NZSP25_50 --het") #Calculates observed and expected autosomal homozygous genotype counts for each sample in dataset
#Import sample IDs, observed and expected homozygosity values into R
plot(Homozygosity$O.HOM, Homozygosity$E.HOM, xlab = "O(HOM)", ylab = "E(HOM)") 
text(Homozygosity$O.HOM, Homozygosity$E.HOM, labels = substring(Homozygosity$ID, 4, 6), cex = 0.8, pos = 1) #Add sample ID identifiers to plot points
```
# Effective population size
VCF files were created in Plink using the smaller filtered Stacks and UNEAK BED files for conversion into GENEPOP files using PGDspider and then Ne estimation with NeEstimator
```r
#Set working directory in R and PATH in console
shell("plink --bfile NZSP25_50 --recode vcf --out NZSP25_50") #Recode UNEAK BED files into VCF
shell("plink --bfile Stacks8_10_50 --recode vcf --out Stacks8_10_50")#Recode Stacks BED files into VCF
```
# DNA Sexing
A binomial test in R was used to investigate if sex biases in sample groups were statistically significant
```r
#example of binomial test for one sample group
binom.test(x=11, n=29, alternative = 'two.sided') 
```
