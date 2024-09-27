# This script is used to calculate the TS scores
# of pan-cancer lncRNAs and non-driver lncRNAs in 14 cancer types.
# pandri_125710.csv is a file that contains mulitple groups of lncRNAs



rm(list = ls())
source("./code/utils.R")
# prepare the input file
pandri_TSscore_file(
  pan_cancer_lnc = "./data/WMDS_latest/gene_14.xlsx",
  all_lnc = "./data/WMDS_latest/LNC_name.xlsx",
  method = c(1, 2, 5, 7, 10),
  output_path = "./data/pandri_TS_scores/"
)
# calc
# pandri
# the input fpkm file is domwloaded from UCSC Xena
lnc_cancer_TSscores(
  input_file_path = "./data/TCGA_14cancer/fpkm/",
  cancer_list = colnames(read.csv("./data/pandri_TS_scores/pandri_125710.csv"))[-1],
  pan_cancer_lnc = "./data/pandri_TS_scores/pandri_125710.csv",
  output_path = "./data/pandri_TS_scores/",
  csv_method = T
)
# non-drivers
lnc_cancer_TSscores(
  input_file_path = "./data/TCGA_14cancer/fpkm/",
  cancer_list = colnames(read.csv("./data/pandri_TS_scores/nondri_all.csv"))[-1],
  pan_cancer_lnc = "./data/pandri_TS_scores/nondri_all.csv",
  output_path = "./data/pandri_TS_scores/",
  csv_method = T
)
