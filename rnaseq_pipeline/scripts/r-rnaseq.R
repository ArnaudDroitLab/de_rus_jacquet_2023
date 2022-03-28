library(DESeq2)
library(readr)
library(rnaseq)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

indir = args[1]
prefix = args[2]
anno = args[3]

print(indir)
print(prefix)
print(anno)

fnames = get_filenames(indir)

txi = import_kallisto(fnames, anno=anno, ignoreTxVersion=TRUE)
saveRDS(txi, paste(prefix, "txi.rds", sep="/"))

tpm = get_tpm_anno_df(txi)
write_csv(tpm, paste(prefix, "tpm.csv", sep='/'))


