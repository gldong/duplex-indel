#### Extract indel candidates from indel_calls.txt for pileup in unmerged bam ####
options(stringsAsFactors = FALSE)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
input_f <- args[1]
output_f <- args[2]

indel_calls <- read.csv(input_f, header=FALSE, skip=13, sep=c("\t", " "), na.strings="")
ds_calls <- indel_calls[indel_calls$V1 == "NV",c(1:5,9)]
ss_calls <- indel_calls[indel_calls$V1 == "DV",c(1:5,9)]
if (nrow(ds_calls)!=0){
	ds_calls$barcode <- str_split(str_split(ds_calls$V9, ":", simplify=T)[,1], ".mem.", simplify=T)[,2]
	ds_calls$V9 <- NULL
}
if (nrow(ss_calls)!=0){
	ss_calls$barcode <- str_split(str_split(ss_calls$V9, ":", simplify=T)[,1], ".mem.", simplify=T)[,2]
	ss_calls$V9 <- NULL
}
df <- rbind(ds_calls, ss_calls)
write.table(df, output_f, quote=F, sep="\t", row.names = F, col.names = F)
