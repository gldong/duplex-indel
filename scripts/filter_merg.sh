#!/bin/bash

#-----------------------------------------------
# Filter: read merging
#-----------------------------------------------
# Filter out variants that are not verified in unmerged bam (pile up of unmerged bam at variant candidate sites)

# Uncomment the line below for cluster computing, otherwise make sure you have these tools installed
# module load gcc/14.2.0 R/4.4.0 samtools/1.21 htslib/1.21 bcftools/1.21

SAMPLE_ID=$1
bulk_bam=$2
base_dir=$3
ref=$4
script_dir=.
barcodes_dir=../references/barcodes
output_suffix=filtered_merg

echo ${SAMPLE_ID} ${bulk_bam}

# Extract variant site with barcode for pileup
input_f=${base_dir}/${SAMPLE_ID}.indel_calls.vcf
output_f=${base_dir}/indel_candidates.txt
Rscript ${script_dir}/filter_merg_extract_candidates.R ${input_f} ${output_f}

# Extract reads from barcode bam and pileup at each candidate site
rm -rf ${base_dir}/indel_candidates_unmerged # remove any old files
mkdir -p ${base_dir}/indel_candidates_unmerged
bam_file=${base_dir}/${SAMPLE_ID}.unmerged.mem.bam
candidate_file=${base_dir}/indel_candidates.txt
if [ $(($(wc -l ${candidate_file} | awk '{print $1}')*3)) -ne $(ls ${base_dir}/indel_candidates_unmerged | wc -l) ]; then
while read -r line; do
	arr=(${line//\t/ })
	mut_type=${arr[0]}
	chr=${arr[1]}
	pos=${arr[2]}
	barcode=${arr[5]}
	region=${chr}:${pos}-${pos}
	tagfile=${barcodes_dir}/META-CS_fwdrev_barcodes_${barcode}.txt
	out_bam=${base_dir}/indel_candidates_unmerged/${mut_type}_${chr}_${pos}.bam
	out_vcf=${base_dir}/indel_candidates_unmerged/${mut_type}_${chr}_${pos}.vcf
	samtools view -b -h --tag-file BC:${tagfile} ${bam_file} ${region} > ${out_bam} && samtools index ${out_bam}
	../pileup -nL1 -r ${region} -P20 -C -q20,30 -Q20,30 -s4 -c -f ${ref} ${bulk_bam} ${out_bam} > ${out_vcf}
	sleep .1s
done < ${candidate_file}
fi

# Filter out indel candidates that cannot be verified in unmerged bam and output filtered calls
Rscript ${script_dir}/filter_merg_filter_candidates.R ${input_f} ${base_dir} ${SAMPLE_ID} ${output_suffix}






