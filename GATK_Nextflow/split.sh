#! /usr/bin/env bash

#Splitting multisample vcf file to individual sample vcf files
mkdir SplitResults
echo Splitting
for file in *.vcf*; do
	echo $file 
	for sample in `bcftools query -l $file`; do
		echo $sample    
		bcftools view -c1 -Oz -s $sample -o SplitResults/${file/.vcf*/.$sample.vcf.gz} $file
	done
done

#Uncompressing individual sample vcf files
echo Gunzipping
cd SplitResults
for file in *.vcf*; do
	echo $file
	gunzip $file
done

#Extracting SNPs from vcf files
echo Selecting Variants
mkdir SNP_vcfs
for vcf in *.vcf; do
echo $vcf
gatk SelectVariants --select-type SNP -V $vcf -O SNP_vcfs/SNPs_$vcf
done

#Converting vcf files to nexus matrix files
echo Generating nexus matrix files
cd SNP_vcfs
mkdir SNPs_nex 
for file in *.vcf; do
python ../../vcf2phylip.py -i $file -p -n --output-folder SNPs_nex
done

#Replacing asterisks with N symbol
cd SNPs_nex
echo Replacing incompatible symbols with N
for file in *min1.nexus; do
sed 's/*/N/g' $file > phylo_$file 
done

#The output files with the 'phylo_' are then to be fed into BEAUTi to allow for parameter configuration after which the resulting XML file is passed into BEAST.
#FigTree can be used for tree visualization.
