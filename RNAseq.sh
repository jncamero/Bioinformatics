#############################
############Running in-line
#############################
#Set directory to alignment files
cd /path/to/bam/

refdir=/path/to/reference/directory/
refas=/path/to/reference/fasta
seqdict='Ref.dict'

#Create Ref.dict file
java -jar $pic/picard.jar CreateSequenceDictionary R=$refas O=$seqdict
#Adding STAR, GATK to PATH 
export GATK_LOCAL_JAR=/path/to/gatk.jar
export PATH=$PATH:/path/to/picard.jar
export PATH=$PATH:/path/to/gatk-4.1.9.0

#INDEXING
STAR --runMode genomeGenerate --genomeDir $ref --genomeFastaFiles $refas --genomeSAindexNbases 13 --runThreadN 20

#1: Alignment of raw data with STAR
#For paired reads here, I substitute the "R1" ending for "R2" when I specify R2 in the call to the STAR aligner
for i in *_R1*.fastq.gz; do
STAR --runMode alignReads --genomeLoad  LoadAndKeep --readFilesCommand zcat --outSAMtype BAM Unsorted --genomeDir $ref --readFilesIn $i ${i%_R1_001.fastq.gz}_R2_001.fastq.gz --runThreadN 10 --outFileNamePrefix ${i%_R1_001.fastq.gz}
done

#2: 2-pass Alignment
mkdir 2pass.sj 
mkdir align-2
mv *SJ.out.tab 2pass.sj

cd 2pass.sj 
for i in *SJ.out.tab ; do
cd $zz
STAR --runMode genomeGenerate --genomeDir $ref --genomeFastaFiles $refas --genomeSAindexNbases 13 --runThreadN 20

#CREATE INDEX TO ADD EVALUATION OF SPLICE JUNCTIONS TO REGULAR REFERENCE ALIGNMENT
cd 2pass.sj
a=$(echo ${i} | sed 's/SJ\.out\.tab/''/g')

STAR --runMode genomeGenerate --genomeDir $ref --genomeFastaFiles $refas \
--sjdbFileChrStartEnd $i --sjdbOverhang 75 --genomeSAindexNbases 13 --runThreadN 20 

cd ../align-2
#FINAL ALIGNMENT
#Variable for file name

STAR --genomeDir $ref --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix $a --readFilesIn $zz/$a'_R1_001.fastq.gz' $zz/$a'_R2_001.fastq.gz' --runThreadN 20

done

#################################################
#################################################
#If necessary,combine all bam files for each entry (i.e. if data spread across multiple files)
#https://gatkforums.broadinstitute.org/gatk/discussion/6057/i-have-multiple-read-groups-for-1-sample-how-should-i-pre-process-them

java -jar picard.jar MergeSamFiles \
I=Sample1.lane1.Aligned.out.bam \
I=Sample1.lane2.Aligned.out.bam \
I=Sample1.lane3.Aligned.out.bam \
I=Sample1.lane4.Aligned.out.bam \
O=Sample1.out.bam &

#Validate SAM files if we want
java -jar picard.jar ValidateSamFile \
I=Some.bam.file.bam \
O=cheeks \
MODE=SUMMARY

##############################################################
##############################################################
#2. Marking duplicates, and index

mkdir dup
for line in FG*;
do

#echo $line
x=$(echo ${line} | sed s/'out.'*'bam'/'tot.md.bam'/g)
#echo $x
java -jar picard.jar MarkDuplicates \
I=${line} \
O=$x \
M=$x'.txt'
echo $x
done
 
#3. Filter based on CIGAR string quality
#set output directory
out3dir=/set/out/directory/
ls|grep md.bam$ > files

for line in $(cat files); do
x=$(echo ${line} | sed s/'md.bam'/'cig.bam'/g)
echo $x
gatk SplitNCigarReads \
-R $refas \
-I ${line}  \
-O $outdir$x
done

#If necessary, change read groups (in parallel)
ls|grep .bam>names.bam
cat names.bam|parallel '
pic=/usr/local/bin
  line=$(echo {});
  #Renaming out file suffix from bam to nom.bam
  out=$(echo $line| sed 's/bam/nom\.bam/g');
  #Giving all reads in bam file a single read group identifier
  sn=$(echo $line|sed 's/\.bam//g');
           
java -jar $pic/picard.jar AddOrReplaceReadGroups \
       I=$line \
       O=$out \
       RGID=$sn \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=$sn'

#INDEX Files: PARALLEL 
ls|grep .*sorted.nom.bam| parallel '
samtools index {}'

#4 Merge cleaned bam files, and index
samtools merge -r Merged.bam *.nom.bam 


#5. Haplotype Caller; call individaual GVCF files 
ls|grep .*sorted.nom.bam| parallel '
ref=home/jncameron/Ref;
gee=/usr/local/bin/gatk-4.1.9.0; 
line=$(echo {})
out=$(echo $line|sed 's/\.sorted\.nom\.bam/\.vcf\.gz/g')
 $gee/gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $refas \
   -I {} \
   -O $out \
   -ploidy 4 \
   -ERC GVCF &
'
#6. COMBINE GVCFs
nohup gatk CombineGVCFs \
   -R $refas \
   --variant Sample1.vcf.gz \
   --variant Sample2..vcf.gz \
   --variant Sample3..vcf.gz \
   --variant Sample4..vcf.gz \
   -O combined.g.vcf.gz &
 
#7 Genotyping group
nohup gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R $refas \
   -V combined.g.vcf.gz \
   -O combined.genotype.vcf.gz &

#8 VARIANT FILTERING
   gatk VariantFiltration \
   -R $refas \
   -V combined.genotype.vcf.gz \
   -O output-filtered.vcf.gz \
   --filter-name "my_filter3" \
   --filter-expression "DP < 20"
  
#EVENTLENGTH (length of the event)
# TRANSITION (1 for a bi-allelic transition (SNP), 0 for bi-allelic transversion (SNP), -1 for INDELs and multi-allelics)
# HET (count of het genotypes)
# HOM-REF (count of homozygous reference genotypes)
# HOM-VAR (count of homozygous variant genotypes)
# NO-CALL (count of no-call genotypes)
# TYPE (type of variant, possible values are NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, and MIXED
# VAR (count of non-reference genotypes)
# NSAMPLES (number of samples)
# NCALLED (number of called samples)
# MULTI-ALLELIC (is this variant multi-allelic? true/false)
 
 
#9 Create table of reads from
     gatk VariantsToTable \
     -V output-filtered.vcf.gz \
     -F CHROM -F POS -F TYPE -F DP -F TRANSITION -F NCALLED -F EVENTLENGTH -GF AD \
     --split-multi-allelic \
     -O output.table
