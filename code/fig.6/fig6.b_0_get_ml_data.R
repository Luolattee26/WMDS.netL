# This code is to obtain joint data from GTEx as well as TCGA
# for machine learning model construction



# init code
rm(list = ls())
gc()
source("./code/utils_TCGAmodel.R")
cancer <- "THCA"
tissue <- "Thyroid"


# # This gtex data also can be downloaded from Xena browser
# # https://xenabrowser.net/datapages/?cohort=GTEX&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# # read gtex data
# # old code to get GTEx Rdata form gtex_RSEM_gene_fpkm file
# data1 <- read.table('./data/GTEx/gtex_RSEM_gene_fpkm', header = T, sep = '\t', check.names = F)
# GTEx_data <- data1
# colnames(GTEx_data)[1] <- 'Ensembl_ID'
# # log2 transform
# # from log2(fpkm+0.001) to log2(fpkm+1)
# tmp <- GTEx_data[2:dim(GTEx_data)[2]]
# tmp <- apply(tmp, 2, function(x) log2(2 ^ as.numeric(x) - 0.001 + 1))
# GTEx_data[2:dim(GTEx_data)[2]] <- tmp
# save(GTEx_data, file = './data/GTEx/GTEx_form_Xena_in_log2(fpkm+1).Rdata')

# read GTEx Rdata (in log2(fpkm+1))
path_tcga <- paste("./data/TCGA_14cancer/fpkm/TCGA-", cancer, ".htseq_fpkm.tsv", sep = "")
load("./data/GTEx/GTEx_form_Xena_in_log2(fpkm+1).Rdata")
raw_data <- read.table(path_tcga, header = T, sep = "\t", check.names = F)
cox_data <- read.csv("./output/pan_cancer_cox/cox_p.csv", row.names = 1)

# ID trans for merge
GTEx_data$ID <- lapply(GTEx_data$Ensembl_ID, function(x) {
  str_split(x, "\\.")[[1]][1]
})
raw_data$ID <- lapply(raw_data$Ensembl_ID, function(x) {
  str_split(x, "\\.")[[1]][1]
})
GTEx_data$ID <- as.character(GTEx_data$ID)
GTEx_data <- GTEx_data[!duplicated(GTEx_data$ID), ]
raw_data$ID <- as.character(raw_data$ID)

# get tissue sample from GTEx data
gtex_phenotype <- read.table("./data/GTEx/GTEX_phenotype", header = T, sep = "\t", check.names = F)
tissues <- c(tissue)
lst <- get_tissue_sample(
  gtex_phenotype = gtex_phenotype,
  tissues = tissues
)
k <- colnames(GTEx_data) %in% lst[[tissue]]
table(k)
tmp <- cbind(GTEx_data$Ensembl_ID, GTEx_data[k], GTEx_data$ID)
colnames(tmp)[ncol(tmp)] <- "ID"
colnames(tmp)[1] <- "Ensembl_ID"

# get merged data(TCGA + GTEx)
merged_data <- merge(raw_data, tmp, by = "ID", all = FALSE)
cols_to_delete <- c("ID", "Ensembl_ID.y")
merged_data <- subset(merged_data, select = -which(names(merged_data) %in% cols_to_delete))

# get_samples
data <- get_samples(
  exp_data = merged_data,
  cancer_type = cancer,
  threshold = 0.05,
  cox_data = cox_data,
  path = "./output/expData_MLmodel/combine_GTEx/useTopFeatures/",
  method = "best",
  useTop = 5,
  padj = F
)
dim(data)

# PCA to check batch effect
tmp <- t(data)
tmp <- data.frame(tmp)
# get normal sample
tmp <- subset(tmp, tmp$status == 1)
runPCA(tmp)

# combat to remove batch effect
tmp_combat <- t(data)
tmp_combat <- data.frame(tmp_combat)
corrected_data <- runComBat(tmp_combat,
                            cancer_type = cancer,
                            path = "./output/expData_MLmodel/combine_GTEx/useTopFeatures/"
)
print("here is the overview of samples:")
print(table(unlist(corrected_data["status", ])))

# check again
tmp <- t(corrected_data)
tmp <- data.frame(tmp)
# get normal sample
tmp <- subset(tmp, tmp$status == 1)
runPCA(tmp)
