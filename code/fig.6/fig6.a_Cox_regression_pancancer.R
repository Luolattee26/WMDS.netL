# This code is used to perform the pan cancer cox regression analysis



# init code
rm(list = ls())
source("./code/utils_TCGAmodel.R")


# data input
# TCGA exp data
file_path <- "./data/TCGA_14cancer/fpkm/"
fileName <- dir(file_path)
cancer_list_all <- fileName2cancerList(fileName)
filePath <- sapply(fileName, function(x) {
  paste(file_path, x, sep = "")
})
exp_data <- lapply(filePath, function(x) {
  read.table(x, header = T, sep = "\t", check.names = F)
})
names(exp_data) <- cancer_list_all

# survival data
# the survival data is from TCGA (Xena platform)
# all cancer
file_path <- "./data/survival_data/"
fileName <- dir(file_path)
cancer_list <- fileName2cancerList(fileName)
filePath <- sapply(fileName, function(x) {
  paste(file_path, x, sep = "")
})
survival_data <- lapply(filePath, function(x) {
  read.table(x, header = T, sep = "\t")
})
names(survival_data) <- cancer_list_all

# input driver lnc data
# here use pandri-7
pan_7_driver <- read.csv("./data/pancancer_driver/pan_dri_7.csv")
lncRNA_list <- unlist(pan_7_driver$lnc)
all_lnc <- read_excel("./data/WMDS_latest/LNC_name.xlsx", col_names = F)
# get driver lnc list with version
# just use the exp_data[[1]] to transfer the ID
# all data is same
ID_noV <- lapply(exp_data[[1]]$Ensembl_ID, function(x) {
  str_split(x, "\\.")[[1]][1]
})
k <- ID_noV %in% lncRNA_list
lncRNA_list <- exp_data[[1]]$Ensembl_ID[k]

# pan cancer cox using pan cancer survival data
pan_cox_result <- c()
KM_group_threshold <- 0.5
for (cancer_index in 1:length(cancer_list_all)) {
  # get single cancer data
  cancer_name <- names(survival_data)[cancer_index]
  one_sur_data <- survival_data[[cancer_index]]
  cancer_exp <- exp_data[[cancer_name]]
  fixed_samples <- seprate_samples(cancer_exp)
  tmp_tumor_df <- fixed_samples[[2]]
  rownames(tmp_tumor_df) <- cancer_exp[, 1]
  # FPKM to TPM
  tmp_tumor_df <- as.matrix(tmp_tumor_df)
  # select survival data
  key <- one_sur_data$sample %in% colnames(tmp_tumor_df)
  table(key)
  one_sur_data <- one_sur_data[key, ]
  k1 <- one_sur_data$OS.time >= 30
  table(k1)
  k2 <- !(is.na(one_sur_data$OS.time) | is.na(one_sur_data$OS))
  table(k2)
  one_sur_data <- one_sur_data[k1 & k2, ]
  # select the overlap data
  k <- intersect(colnames(tmp_tumor_df), one_sur_data$sample)
  k_tmp <- one_sur_data$sample %in% k
  one_sur_data <- one_sur_data[k_tmp, ]
  tmp_tumor_df <- tmp_tumor_df[, k]
  one_sur_data <- one_sur_data[match(colnames(tmp_tumor_df), one_sur_data$sample), ]
  # get KM group
  tmp_exp_group <- add_exp_feature(tmp_tumor_df, lncRNA_list, KM_group_threshold)
  # mix survival with exp data
  tmp_survival_mixed_expression <- mix_exp(one_sur_data, tmp_exp_group)
  # cox
  tmp_all_cox_result <- data.frame()
  for (i in 5:ncol(tmp_survival_mixed_expression)) {
    cox_result <- coxph(Surv(OS.time, OS) ~ tmp_survival_mixed_expression[, i], data = tmp_survival_mixed_expression)
    tmp <- summary(cox_result)
    pvalue_wald <- tmp$waldtest[3]
    beta <- tmp[["coefficients"]][1]
    HR <- tmp[["coefficients"]][2]
    se_beta <- tmp[["coefficients"]][3]
    pvalue_likelihood <- tmp$logtest[3]
    df <- data.frame(
      gene_ID = colnames(tmp_survival_mixed_expression)[i],
      pvalue_wald = tmp$waldtest[3],
      beta = tmp[["coefficients"]][1],
      HR = tmp[["coefficients"]][2],
      se_beta = tmp[["coefficients"]][3],
      pvalue_likelihood = tmp$logtest[3],
      n = tmp$n
    )
    rownames(df) <- df[, 1]
    df <- df[, -1]
    tmp_all_cox_result <- rbind(tmp_all_cox_result, df)
  }
  pan_cox_result[[cancer_index]] <- tmp_all_cox_result
}
names(pan_cox_result) <- cancer_list

# combine the pan cancer data into dataframe
# and save it
pan_cox_result_mixed <- data.frame(a = seq(1, length(lncRNA_list)))
pan_cox_p_mixed <- data.frame(a = seq(1, length(lncRNA_list)))
for (cancer_index in 1:length(cancer_list)) {
  cancer_name <- names(pan_cox_result)[cancer_index]
  one_cox_data <- pan_cox_result[[cancer_index]]
  pan_cox_result_mixed <- cbind(pan_cox_result_mixed, one_cox_data$HR)
  pan_cox_p_mixed <- cbind(pan_cox_p_mixed, one_cox_data$pvalue_wald)
}
pan_cox_result_mixed <- pan_cox_result_mixed[, -1]
colnames(pan_cox_result_mixed) <- names(pan_cox_result)
rownames(pan_cox_result_mixed) <- rownames(one_cox_data)
symbol_list_pan <- unlist(get_symbol_list(rownames(pan_cox_result_mixed), all_lnc))
rownames(pan_cox_result_mixed) <- symbol_list_pan
pan_cox_p_mixed <- pan_cox_p_mixed[, -1]
colnames(pan_cox_p_mixed) <- names(pan_cox_result)
rownames(pan_cox_p_mixed) <- rownames(one_cox_data)
symbol_list_pan <- unlist(get_symbol_list(rownames(pan_cox_p_mixed), all_lnc))
rownames(pan_cox_p_mixed) <- symbol_list_pan
# set path
filePath_cox <- "./output/pan_cancer_cox"
fileName1 <- paste(filePath_cox, "/", "cox_HR.csv", sep = "")
fileName2 <- paste(filePath_cox, "/", "cox_p.csv", sep = "")
write.csv(pan_cox_result_mixed, fileName1)
write.csv(pan_cox_p_mixed, fileName2)
