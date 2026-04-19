#### Extract indel candidates from indel_calls.vcf for pileup in unmerged bam ####
options(stringsAsFactors = FALSE)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
input_f <- args[1]
output_f <- args[2]

indel_calls <- read.table(pipe(paste("grep -v '^#'", input_f)), header=FALSE, sep="\t", na.strings="", quote="")
ds_calls <- indel_calls[grepl("VT=DS", indel_calls$V8), c(1, 2, 4, 5, 10)]
ss_calls <- indel_calls[grepl("VT=SS", indel_calls$V8), c(1, 2, 4, 5, 10)]
if (nrow(ds_calls)!=0){
	ds_calls$var_type <- "DS"
	ds_calls$barcode <- str_split(str_split(ds_calls$V10, ":", simplify=T)[,1], ".mem.", simplify=T)[,2]
	ds_calls <- ds_calls[, c("var_type", "V1", "V2", "V4", "V5", "barcode")]
} else { 
	ds_calls <- data.frame() 
}
if (nrow(ss_calls)!=0){
	ss_calls$var_type <- "SS"
	ss_calls$barcode <- str_split(str_split(ss_calls$V10, ":", simplify=T)[,1], ".mem.", simplify=T)[,2]
	ss_calls <- ss_calls[, c("var_type", "V1", "V2", "V4", "V5", "barcode")]
} else { 
	ss_calls <- data.frame() 
}
df <- rbind(ds_calls, ss_calls)
write.table(df, output_f, quote=F, sep="\t", row.names = F, col.names = F)
