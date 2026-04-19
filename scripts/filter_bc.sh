#!/bin/bash

#-----------------------------------------------
# Filter: unique barcode pair
#-----------------------------------------------
# Filter out variants covered by multiple barcode pairs

# Uncomment the line below for cluster computing, otherwise make sure you have these tools installed
# module load gcc/14.2.0 htslib/1.21 bcftools/1.21 

# Parse arguments
SAMPLE_ID=$1
SAMPLE_DIR=$2
OUTPUT_SUFFIX=filtered_bc

# Filter ds calls
cat ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.vcf | awk '!/^#/ && $8 ~ /VT=DS/ { split($10,BC,":"); printf("%s\t%s\t%s\n", $1, $2, BC[1])}' > ${SAMPLE_DIR}/ds_region_with_BC_bam.txt
if [ ! -s ${SAMPLE_DIR}/ds_region_with_BC_bam.txt ]; then # if no ds calls, write empty VCF
  echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf; 
else
  bcftools view -H -v indels -R ${SAMPLE_DIR}/ds_region_with_BC_bam.txt --regions-overlap 0 ${SAMPLE_DIR}/${SAMPLE_ID}.vcf.gz | \
  if [ $(wc -l - | cut -d ' ' -f1) -ne $(wc -l ${SAMPLE_DIR}/ds_region_with_BC_bam.txt | cut -d ' ' -f1) ]; then 
    echo ${SAMPLE_ID} "ERROR: mismatching double-stranded variant numbers between indel_calls and intermediate vcf"; continue;
  else
    bcftools view -H -v indels -R ${SAMPLE_DIR}/ds_region_with_BC_bam.txt --regions-overlap 0 ${SAMPLE_DIR}/${SAMPLE_ID}.vcf.gz | \
    awk ' BEGIN{FS="\t";OFS="\t"} {\
      if(NR==1) print "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
      flt=0; \
      for(i=11;i<=NF;i+=2) {\
        split($i,mem,":"); \
        split(mem[2],memADF,","); split(mem[3],memADR,","); \
        if(memADF[1]+memADF[2]+memADR[1]+memADR[2] > 0) flt++; \
      } \
      if(flt == 1) print $1,$2,$3,$4,$5,".",".","."
    }' > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
  fi
fi

# Filter ss calls
cat ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.vcf | awk '!/^#/ && $8 ~ /VT=SS/ { split($10,BC,":"); printf("%s\t%s\t%s\n", $1, $2, BC[1])}' > ${SAMPLE_DIR}/ss_region_with_BC_bam.txt
if [ ! -s ${SAMPLE_DIR}/ss_region_with_BC_bam.txt ]; then # if no ss calls, write empty VCF
  echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ss.${OUTPUT_SUFFIX}.vcf; 
else
  bcftools view -H -v indels -R ${SAMPLE_DIR}/ss_region_with_BC_bam.txt --regions-overlap 0 ${SAMPLE_DIR}/${SAMPLE_ID}.vcf.gz | \
  if [ $(wc -l - | cut -d ' ' -f1) -ne $(wc -l ${SAMPLE_DIR}/ss_region_with_BC_bam.txt | cut -d ' ' -f1) ]; then 
    echo ${SAMPLE_ID} "ERROR: mismatching single-stranded variant numbers between indel_calls and intermediate vcf"; continue;
  else
    bcftools view -H -v indels -R ${SAMPLE_DIR}/ss_region_with_BC_bam.txt --regions-overlap 0 ${SAMPLE_DIR}/${SAMPLE_ID}.vcf.gz | \
    awk ' BEGIN{FS="\t";OFS="\t"} {\
      if(NR==1) print "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
      flt=0; \
      for(i=11;i<=NF;i+=2) {\
        split($i,mem,":"); \
        split(mem[2],memADF,","); split(mem[3],memADR,","); \
        if(memADF[1]+memADF[2]+memADR[1]+memADR[2] > 0) flt++; \
      } \
      if(flt == 1) print $1,$2,$3,$4,$5,".",".","."
    }' > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ss.${OUTPUT_SUFFIX}.vcf
  fi
fi



