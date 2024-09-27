# This script is used to find the lncRNAs that
# are significantly associated with survival in multiple cancer types.



# init code
rm(list = ls())
source("./code/utils_TCGAmodel.R")


# read cox result
cox_data <- read.csv("./output/pan_cancer_cox/cox_p.csv", row.names = 1)
# sig count for each cancer
sig_count <- cox_sig_count(
  cox_data = cox_data,
  threshold = 0.05,
  padj = F,
  padj.method = "BH"
)
# how many times being sig. for each lnc
sig_count_forLNC <- cox_sig_count_forLNC(
  cox_data = cox_data,
  threshold = 0.05,
  padj = F,
  padj.method = "BH"
)

# select better lnc
threshold <- 5
lst <- sig_count_forLNC[sig_count_forLNC >= threshold]
all_lnc <- read_excel("./data/WMDS_latest/LNC_name.xlsx", col_names = F)
lst_ID <- get_ID_list(names(lst), all_lnc)


# save result
save(list = c("lst", "lst_ID"), file = "./data/important_lnc_in_cox.Rdata")


sig_count_forLNC <- sig_count_forLNC[sig_count_forLNC >= 5]
