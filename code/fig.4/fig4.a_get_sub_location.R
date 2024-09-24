# This code is used to get the subcellular localization of lncRNA
# in HepG2 and K562 cell lines.



# init code
rm(list = ls())
source("./code/utils_TCGAmodel.R")


# get data list
# These data are downloaded from ENCODE project
data_list <- get_heatmap_data_list("./data/ENCODE_SUBCELLULAR_LOCALIZATION/transcript_quantifications/")

# gtf
# The gtf file is downloaded from GENCODE project
gtf_file <- import("./data/gencode.v22.annotation.gtf")
gtf_df <- as.data.frame(gtf_file)
gtf_all_trans <- gtf_df[gtf_df$type == "transcript", ]
# delete the transcript which is not lncRNA
Retained_lst <- c(
  "non_coding", "3prime_overlapping_ncrna",
  "antisense", "lincRNA", "retained_intron",
  "sense_intronic",
  "sense_overlapping", "macro_lncRNA"
)
k <- gtf_all_trans$transcript_type %in% Retained_lst
table(k)
gtf_all_trans <- gtf_all_trans[k, ]

# read lnc files
# input pandri
file_path <- "./data/pancancer_driver/"
index <- 0
pan_list <- list()
pancancer_driver_list <- c(1, 2, 5, 7)
for (i in pancancer_driver_list) {
  index <- index + 1
  fileName <- paste(file_path, "pan_dri_", i, ".csv", sep = "")
  pan_list[[index]] <- read.csv(fileName, header = T, sep = ",", check.names = F)
}
names(pan_list) <- pancancer_driver_list
# all lncRNA
all_LNC <- read_xlsx("./data/WMDS_latest/LNC_name.xlsx", col_names = F)

# ID trans
all_ID_v <- gtf_all_trans$gene_id
all_ID <- lapply(all_ID_v, function(x) {
  return(str_split_1(x, "\\.")[1])
})
gtf_all_trans$ID <- all_ID

# final result
result <- data.frame()
# non-drivers
k_non <- all_LNC$...1 %in% pan_list[["1"]]$lnc
table(k_non)
all_LNC <- all_LNC[!k_non, ]
k <- gtf_all_trans$ID %in% all_LNC$...1
table(k)
gtf_tmp <- gtf_all_trans[k, ]
res_tmp <- Get_SUBCELLULAR_LOCALIZATION_exp(data_list, gtf_tmp)
result <- rbind(result, res_tmp)
colnames(result) <- names(data_list)
# all drivers
k <- gtf_all_trans$ID %in% pan_list[["1"]]$lnc
table(k)
gtf_tmp <- gtf_all_trans[k, ]
res_tmp <- Get_SUBCELLULAR_LOCALIZATION_exp(data_list, gtf_tmp)
result <- rbind(result, res_tmp)
colnames(result) <- names(data_list)
# more than 2 drivers
k <- gtf_all_trans$ID %in% pan_list[["2"]]$lnc
table(k)
gtf_tmp <- gtf_all_trans[k, ]
res_tmp <- Get_SUBCELLULAR_LOCALIZATION_exp(data_list, gtf_tmp)
result <- rbind(result, res_tmp)
colnames(result) <- names(data_list)
# more than 5 drivers
k <- gtf_all_trans$ID %in% pan_list[["5"]]$lnc
table(k)
gtf_tmp <- gtf_all_trans[k, ]
res_tmp <- Get_SUBCELLULAR_LOCALIZATION_exp(data_list, gtf_tmp)
result <- rbind(result, res_tmp)
colnames(result) <- names(data_list)
# more than 7 drivers
k <- gtf_all_trans$ID %in% pan_list[["7"]]$lnc
table(k)
gtf_tmp <- gtf_all_trans[k, ]
res_tmp <- Get_SUBCELLULAR_LOCALIZATION_exp(data_list, gtf_tmp)
result <- rbind(result, res_tmp)
colnames(result) <- names(data_list)
rownames(result) <- c(
  "all non-driver", "all drivers", "more than 2",
  "more than 5", "more than 7"
)
df_HepG2 <- result[, 1:4]
df_K562 <- result[, -c(1:4)]


# save
save(
  list = c("result", "df_HepG2", "df_K562"),
  file = "./data/ENCODE_SUBCELLULAR_LOCALIZATION/transcript_quantifications/data_cleaned.rdata"
)
