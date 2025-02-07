# This script is used to generate fasta file for all LNCs in the LNC_name.xlsx file



rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyr)
library(dplyr)
# use remotes to install rutil
# if (!('remotes' %in% .packages(T))) install.packages('remotes');
# remotes::install_github('dongzhuoer/rutil');
library(rutil)
if (!require("rtracklayer")) {
  BiocManager::install("rtracklayer")
}else{
  library("rtracklayer")
}
library(readxl)

# input fasta file
fasta_df <- read_fasta("./data/transcripts_latest.fa")

# get transcripts table
trans <- separate(fasta_df, name, into = c('trans_id', 'other'), ' ') %>%
  separate(other, into = c('other2', 'gene_symbol'), '=')
trans <- trans %>% mutate(length = str_length(trans$seq)) %>%
  dplyr::select(c('trans_id', 'gene_symbol', 'seq', 'length'))

# import gtf file
gtf_file <- import('./data/gencode.v22.annotation.gtf')
gtf_df <- as.data.frame(gtf_file)

# import lncDATA
alllnc <- read_xlsx('./data/WMDS_latest/LNC_name.xlsx', col_names = F)

# id translation
gtf_all_trans <- gtf_df[gtf_df$type == 'transcript', ]
all_gene_list <- list()
index <- 0
all_gene_with_V <- gtf_all_trans$gene_id
for (gene in all_gene_with_V) {
  index <- index + 1
  tmp <- str_split(gene, '\\.')[[1]][1]
  all_gene_list[index] <- tmp
}
# get all lnc
k_all <- all_gene_list %in% alllnc$...1
gtf_all_lnc <- gtf_all_trans[k_all, ]

# select fasta file
k <- trans$trans_id %in% gtf_all_lnc$transcript_id
fasta_df_driver <- trans[k, c(1, 3)]
colnames(fasta_df_driver) <- c('name', 'seq')
write_fasta(fasta_df_driver, './data/allLNC_fasta.fa')

