#!/bin/bash

#-----------------------------------------------
# Filter: unique Tn5 site
#-----------------------------------------------
# Filter out variants with multiple Tn5 sites (under the same barcode pair)

SAMPLE_ID=$1
SAMPLE_DIR=$2
OUTPUT_SUFFIX=filtered_tn5

# Filter ds calls
# extract variants with a unique pair of Tn5 start|end sites
echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
cat ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.vcf | \
awk ' BEGIN{FS="\t";OFS="\t"} !/^#/ && $8 ~ /VT=DS/ {\
  split($10,BC,":"); \
  if(BC[4]==1) print $1,$2,".",$4,$5,".",".","."
}' >> ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf

# Filter ss calls
# extract variants with a unique pair of Tn5 start|end sites
echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ss.${OUTPUT_SUFFIX}.vcf
cat ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.vcf | \
awk ' BEGIN{FS="\t";OFS="\t"} !/^#/ && $8 ~ /VT=SS/ {\
  split($10,BC,":"); \
  if(BC[4]==1) print $1,$2,".",$4,$5,".",".","."
}' >> ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ss.${OUTPUT_SUFFIX}.vcf




