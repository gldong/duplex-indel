#### Filter indel calls using unmerged bam ####
# For ds calls, keep 1/1 in unmerged bam; for ss calls, keep 1/0 or 0/1 (i.e. filter out 0/0 or ./.)
# - input: indel_calls
# - output: filtered_calls

options(stringsAsFactors = FALSE)
library(stringr)

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
input_f <- args[1]
base_dir <- args[2]
SAMPLE_ID <- args[3]
output_suffix <- args[4]

# Extract original ds and ss calls
indel_calls <- read.csv(input_f, header=FALSE, skip=13, sep=c("\t", " "), na.strings="")
ds_calls <- indel_calls[indel_calls$V1 == "NV",c(2:5)]
ss_calls <- indel_calls[indel_calls$V1 == "DV",c(2:5)]

# Apply filter to ds and ss calls
filtered_ds_calls <- c()
if (nrow(ds_calls)!=0) {
	for (j in 1:nrow(ds_calls)){
		vcf <- read.table(sprintf("%s/%s/indel_candidates_unmerged/NV_%s_%s.vcf",base_dir,SAMPLE_ID,ds_calls$V2[j],ds_calls$V3[j]), header=F, skip=97, sep="\t")
		if (grepl("1/1",vcf$V11)){ #if indel also exists in unmerged bam
			filtered_ds_calls <- rbind(filtered_ds_calls, ds_calls[j,])
		}
	}
}
filtered_ss_calls <- c()
if (nrow(ss_calls)!=0) {
	for (j in 1:nrow(ss_calls)){
		vcf <- read.table(sprintf("%s/%s/indel_candidates_unmerged/DV_%s_%s.vcf",base_dir,SAMPLE_ID,ss_calls$V2[j],ss_calls$V3[j]), header=F, skip=97, sep="\t")
		if (grepl("1/0",vcf$V11)|grepl("0/1",vcf$V11)) {
			filtered_ss_calls <- rbind(filtered_ss_calls, ss_calls[j,])
		}
	}
}

# Reformat for VCF output
if (!is.null(filtered_ds_calls)){
	if (nrow(filtered_ds_calls)!=0) {
		filtered_ds_calls <- cbind(filtered_ds_calls, ".",".",".",".")
		filtered_ds_calls <- filtered_ds_calls[,c(1,2,5,3,4,6:8)]
	}
}
if (!is.null(filtered_ss_calls)){
	if (nrow(filtered_ss_calls)!=0) {
		filtered_ss_calls <- cbind(filtered_ss_calls, ".",".",".",".")
		filtered_ss_calls <- filtered_ss_calls[,c(1,2,5,3,4,6:8)]
	}
}

# Write output
new_header <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
vcf_head <- paste0("##fileformat=VCFv4.0\n#",
					 paste(new_header,collapse="\t"),"\n")
out <- sprintf("%s/%s.indel_calls.ds.%s.vcf",base_dir,SAMPLE_ID,output_suffix)
cat(vcf_head, file=out)
write.table(filtered_ds_calls, file=out, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
out <- sprintf("%s/%s.indel_calls.ss.%s.vcf",base_dir,SAMPLE_ID,output_suffix)
cat(vcf_head, file=out)
write.table(filtered_ss_calls, file=out, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
