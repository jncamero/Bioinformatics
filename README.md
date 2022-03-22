# Bioinformatics
DNA and RNA sequence analyses

# 1 RNAseq.sh
This is a pipeline for analyzing RNAseq data, via BASH one liners, making use of parellization where possible, starting with reads in fasta files, and ending with a table of variant counts. 
It was developed to run on Ubuntu 18.04.6 LTS OS.  

**Required:**
Samtools,
STAR,
GATK,
Picard Tools 

**Optional:**
bcftools


# 2 Microhaplotyping.sh
Calling phased variants present in targetted sequencing reads

**Required**
Bedtools,
Samtools

1 Microhaplotyping pipeline from targetted sequencing reads
This pipeline calls microhaplotypes when multiple variants are present in the same sequencing read.  
It uses a VCF to identity polymorphic loci, then calls these bp positions directly from the individual BAM files
and returns them phased in a string
The output (i.e. .haps files) can then be passed to the R script MH.processing.R to create a matrix indicating the frequency of
Each mulitallelic microhaplotype at each targetted sequencing locus

The microhaplotyping pipeline takes as input:
1. A vcf file indicating polymorphic loci
2. Individual BAM files from which the VCF was derived
3. A bed file indicating the genomic interval where you want to call a microhaplotype

Directions: 
1. set the directory to the folder containing 1,2 and 3 at the top of MH.sh
2. set the variables "vcf" and "bed" to their respective files at the top of MH.sh
3. set the "proc" variable to the number of files you want to process simultaneously
(On a 12 core linux machine with 65 GB of RAM I find setting proc equal to 35-40 has good performance) 
4. make file executable, and make sure you have read a write permissions to the directory
sudo chmod +x Microhaplotyping.sh
5. Run

nohup ./Microhaplotyping.sh &
