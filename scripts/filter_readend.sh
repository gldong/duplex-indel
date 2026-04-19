#!/bin/bash

#-----------------------------------------------
# Filter: near read ends
#-----------------------------------------------
# Filter out variants within <int> bp from the end of the read

# Uncomment the line below for cluster computing, otherwise make sure you have these tools installed
# module load gcc/14.2.0 htslib/1.21 bcftools/1.21

# Parse arguments
SAMPLE_ID=$1
SAMPLE_DIR=$2
CUTOFF=$3
OUTPUT_SUFFIX=filtered_readend

# Extract ds variant pos on the read and filter variant if pos <= CUTOFF
cat ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.vcf | awk '!/^#/ && $8 ~ /VT=DS/ { split($10,BC,":"); printf("%s\t%s\t%s\n", $1, $2, BC[1])}' > ${SAMPLE_DIR}/ds_region_with_BC_bam.txt
echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
if [ -s ${SAMPLE_DIR}/ds_region_with_BC_bam.txt ]; then # if no ds calls, write empty VCF
  while IFS=$'\t' read -r -a array
  do
    bcftools view -H -v indels -r ${array[0]}:${array[1]} --regions-overlap 0 -s ${array[2]}.bam ${SAMPLE_DIR}/${SAMPLE_ID}.vcf.gz | \
    awk -v cutoff="$CUTOFF" ' BEGIN{FS="\t";OFS="\t"} {\
      split($10,fields,":");
      split(fields[5],ALEN,",");
      if(ALEN > cutoff) print $1,$2,$3,$4,$5,ALEN[2],".","."
    }' >> ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
  done < ${SAMPLE_DIR}/ds_region_with_BC_bam.txt
fi

# Extract ss variant pos on the read and filter variant if pos <= CUTOFF
cat ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.vcf | awk '!/^#/ && $8 ~ /VT=SS/ { split($10,BC,":"); printf("%s\t%s\t%s\n", $1, $2, BC[1])}' > ${SAMPLE_DIR}/ss_region_with_BC_bam.txt
echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ss.${OUTPUT_SUFFIX}.vcf
if [ -s ${SAMPLE_DIR}/ss_region_with_BC_bam.txt ]; then # if no ss calls, write empty VCF
  while IFS=$'\t' read -r -a array
  do
    bcftools view -H -v indels -r ${array[0]}:${array[1]} --regions-overlap 0 -s ${array[2]}.bam ${SAMPLE_DIR}/${SAMPLE_ID}.vcf.gz | \
    awk -v cutoff="$CUTOFF" ' BEGIN{FS="\t";OFS="\t"} {\
      split($10,fields,":");
      split(fields[5],ALEN,",");
      if(ALEN > cutoff) print $1,$2,$3,$4,$5,ALEN[2],".","."
    }' >> ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ss.${OUTPUT_SUFFIX}.vcf
  done < ${SAMPLE_DIR}/ss_region_with_BC_bam.txt
fi



