#!/bin/bash

set -Eeuo pipefail

#-----------------------------------------------
# Argument parser
#-----------------------------------------------
usage() {
  echo "Usage: $0 -s <sample_ID> -c <config_file> -d <sample_dir> [-m] [-t] [-b] [-e <int>]"
  echo
  echo "Options:"
  echo "  -s <sample_ID>    Sample ID [required]"
  echo "  -c <config_file>  Config file [required]"
  echo "  -d <sample_dir>   Output directory for this sample [required]"
  echo "  -m                Enable read merging filter"
  echo "  -t                Enable unique Tn5 site filter"
  echo "  -b                Enable unique barcode pair filter"
  echo "  -e <int>          Enable read end filter within the last <int> bp from the end of the read"
  echo "  -h                Show help message and exit"
}

# Print help message if no arguments are provided
if [ $# -eq 0 ]; then
  usage
  exit 1
fi

# Default
FILTER_MERG=0
FILTER_TN5=0
FILTER_BC=0
FILTER_READEND=0

while getopts "s:c:d:mtbe:h" opt; do
  case $opt in
    s) SAMPLE_ID="$OPTARG" ;;
    c) CONFIG_FILE="$OPTARG" ;;
    d) SAMPLE_DIR="$OPTARG" ;;
    m) FILTER_MERG=1 ;;
    t) FILTER_TN5=1 ;;
    b) FILTER_BC=1 ;;
    e) FILTER_READEND="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

# Shift parsed options away
shift $((OPTIND - 1))

# Check required arguments
if [ -z "${SAMPLE_ID:-}" ]; then
  echo "Error: -s <sample_ID> is required." >&2
  exit 1
fi
if [ -z "${CONFIG_FILE:-}" ]; then
  echo "Error: -c <config_file> is required." >&2
  exit 1
fi
if [ -z "${SAMPLE_DIR:-}" ]; then
  echo "Error: -d <sample_dir> is required." >&2
  exit 1
fi

echo "SAMPLE_ID: ${SAMPLE_ID}"
echo "Applying filters for: "
if [ ${FILTER_MERG} ]; then echo "--read merging"; fi
if [ ${FILTER_TN5} ]; then echo "--unique Tn5 site"; fi
if [ ${FILTER_BC} ]; then echo "--unique barcode pair"; fi
if [ ${FILTER_READEND} ]; then echo "--read end"; fi

SCRIPT_DIR=scripts


#-----------------------------------------------
# Filter: read merging confounding
#-----------------------------------------------
# Filter out variants that are not verified in unmerged bam (pile up of unmerged bam at variant candidate sites)

echo "Apply read merging filter..."
# Extract bulk bam and reference file from config file
SAMPLE_BULK=$(cat ${CONFIG_FILE} | grep bulk_bam: | sed 's/bulk_bam: //' | sed 's/"//g')
REF_FILE=$(cat ${CONFIG_FILE} | grep ref_fasta | sed 's/ref_fasta: //' | sed 's/"//g')
bash ${SCRIPT_DIR}/filter_merg.sh ${SAMPLE_ID} ${SAMPLE_BULK} ${SAMPLE_DIR} ${REF_FILE}


#-----------------------------------------------
# Filter: unique Tn5 site
#-----------------------------------------------
# Filter out variants with multiple Tn5 sites (under the same barcode pair)

echo "Apply unique Tn5 site filter..."
bash ${SCRIPT_DIR}/filter_tn5.sh ${SAMPLE_DIR}


#-----------------------------------------------
# Filter: unique barcode pair
#-----------------------------------------------
# Filter out variants covered by multiple barcode pairs

echo "Apply unique barcode pair filter..."
bash ${SCRIPT_DIR}/filter_bc.sh ${SAMPLE_ID} ${SAMPLE_DIR}


#-----------------------------------------------
# Filter: near read end
#-----------------------------------------------
# Filter out variants within <int> bp from the end of the read

echo "Apply read end filter..."
bash ${SCRIPT_DIR}/filter_readend.sh ${SAMPLE_ID} ${SAMPLE_DIR} ${FILTER_READEND}


#-----------------------------------------------
# Combine filters
#-----------------------------------------------
# Combine applied filters above

echo "Combine filters..."

OUTPUT_SUFFIX=filtered
FILE1=""
FILE2=""
if [ ${FILTER_MERG} ]; then 
  OUTPUT_SUFFIX+="_merg"
  if [ -z ${FILE1} ]; then FILE1=${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.filtered_merg.vcf; fi
fi
if [ ${FILTER_TN5} ]; then 
  OUTPUT_SUFFIX+="_tn5"
  if [ -z ${FILE1} ]; then 
    FILE1=${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.filtered_tn5.vcf
  else
    FILE2=${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.filtered_tn5.vcf
    grep -f ${FILE1} ${FILE2} > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
    FILE1=${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
    FILE2=""
  fi
fi
if [ ${FILTER_BC} ]; then 
  OUTPUT_SUFFIX+="_bc"
  if [ -z ${FILE1} ]; then 
    FILE1=${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.filtered_tn5.vcf
  else
    FILE2=${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.filtered_tn5.vcf
    grep -f ${FILE1} ${FILE2} > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
    FILE1=${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
    FILE2=""
  fi
fi
if [ ${FILTER_READEND} ]; then 
  OUTPUT_SUFFIX+="_readend"
  if [ -z ${FILE1} ]; then 
    FILE1=${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.filtered_tn5.vcf
  else
    FILE2=${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.filtered_tn5.vcf
    grep -f ${FILE1} ${FILE2} > ${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
    FILE1=${SAMPLE_DIR}/${SAMPLE_ID}.indel_calls.ds.${OUTPUT_SUFFIX}.vcf
    FILE2=""
  fi
fi

echo "Final filtered output: ${FILE1}"



