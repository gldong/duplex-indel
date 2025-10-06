#!/bin/bash

#-----------------------------------------------
# Filter: multiple Tn5 sites
#-----------------------------------------------
# Filter out variants with multiple Tn5 sites (under the same barcode pair)

SAMPLE_DIR=$1
OUTPUT_SUFFIX=filtered_tn5

# Filter ds calls
# extract variants with a unique pair of Tn5 start|end sites
echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${SAMPLE_DIR}/${SAMPLE}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
cat ${SAMPLE_DIR}/${SAMPLE}.indel_calls.txt | grep ^NV | \
awk ' BEGIN{FS="\t";OFS="\t"} {\
  split($9,BC,":"); \
  if(BC[4]==1) print $2,$3,".",$4,$5,".",".","."
}' >> ${SAMPLE_DIR}/${SAMPLE}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf

# Filter ss calls
# extract variants with a unique pair of Tn5 start|end sites
echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${SAMPLE_DIR}/${SAMPLE}.indel_calls.ss.${OUTPUT_SUFFIX}.vcf
cat ${SAMPLE_DIR}/${SAMPLE}.indel_calls.txt | grep ^DV | \
awk ' BEGIN{FS="\t";OFS="\t"} {\
  split($9,BC,":"); \
  if(BC[4]==1) print $2,$3,".",$4,$5,".",".","."
}' >> ${SAMPLE_DIR}/${SAMPLE}.indel_calls.ss.${OUTPUT_SUFFIX}.vcf




