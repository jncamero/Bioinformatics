#!/bin/bash
#Set directory to location of your bam files
cd /Directory/with/BAM/files;

#Set location to bed file for targetted read locations
bed='/bed/file'

#Set location to vcf filed called on all bam files
vcf='/vcf/file'

#Set number of processes to run in parallel
proc=35

task(){
cat $bed|while read tline; do
echo $tline
start=$(echo $tline|cut -f2 -d" ")
end=$(echo $tline|cut -f3 -d" ")
chr=$(echo $tline|cut -f1 -d" ")
prob=$(echo $tline|cut -f4 -d" ")
#zjk=$(($zjk+1))
#echo $chr $zjk;

#echo $start $end $chr $prob
#start=$(head -4 $bed |tail -1|cut -f2)
#end=$(head -4 $bed |tail -1|cut -f3)
#chr=$(head -4 $bed |tail -1|cut -f1)

#Position of SNPs in VCF
snps=$(bcftools view --no-header -t $chr:$start-$end $vcf |cut -f2);

#Filter based VARIANT POSITIONS in VCF.
one=$(echo $snps|tr -cd ' \t' | wc -c);
len=$(echo $(($one+1)));
first=$(echo $snps|cut -f1 -d" ");
lst=$(echo $snps|cut -f$len -d" ");
#formatted SNP positions with comma separation
snpform=$(echo $snps|sed -e 's/ /,/g');
nos=$(samtools view -f"0x2" -F"0x800,0x100,0x400,0x200,0x4" $1 "$chr:$first-$lst"|wc -l);

#If there are reads for the region, AND there are variants present in VCF...GO
if (( nos > 0 && ${#snps} > 0 ));then
samtools view -f"0x2" -F"0x800,0x100,0x400,0x200,0x4" $1 "$chr:$first-$lst"|while read line; do
xx=$();

for (( i=1; i<=$len; i++ )); do
#echo $i;
	#Change this for loop
	#Change this for loo
     refpt=$(echo $line|cut -f4 -d" ");

a=$(echo $snps|cut -f$i -d" ");
	#Need to add 1...
#Scaled position of SNP
b=$(($a-$refpt+1));
xx=$xx,$b;
done
#cut first comma off of string
hapind=$(echo $xx| cut -c 2-);
check1=$(echo $line|cut -f4 -d" ");
check2int=$(echo $line|cut -f10 -d" "|wc -c);
check2=$(($check2int+1+$check1));

#if the targetted variants are contained in the read, get haplotype
	#if [ "$check1" -le "$first" ] && [ "$check2" -ge "$lst" ];
if (("check1"<="first" && "check2">="lst"));
then
sequence=$(echo $line|cut -f10 -d" ");
hap=$(echo $sequence|cut -c$hapind);
#echo $hap;
if [ ! -z "$hap" ]; then       
echo $chr $snpform $hap $prob>> $1.haps;
fi
#echo $chr $snpform $hap $prob
#echo $hap| sort| uniq -c |sort -nr
#echo $check1 $check2;
fi

done

fi #END IF

done
}

export -f task;

#Run function in parallel on your bam files

parallel --jobs $proc 'task' ::: *bam; 
