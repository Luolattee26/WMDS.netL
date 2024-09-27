# This code is used to plot the proportion of
# different types of lncRNAs in different groups.



# init code
rm(list = ls())
source("./code/utils_TCGAmodel.R")


# gtf
# the gtf file is downloaded from gencode
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
count_df <- gtf_tmp %>%
  group_by(transcript_type) %>%
  summarise(count = n())
count_df$sample <- rep("Non-drivers", dim(count_df)[1])
result <- bind_rows(result, count_df)
# all drivers
k <- gtf_all_trans$ID %in% pan_list[["1"]]$lnc
table(k)
gtf_tmp <- gtf_all_trans[k, ]
count_df <- gtf_tmp %>%
  group_by(transcript_type) %>%
  summarise(count = n())
count_df$sample <- rep("All drivers", dim(count_df)[1])
result <- bind_rows(result, count_df)
# more than 2 drivers
k <- gtf_all_trans$ID %in% pan_list[["2"]]$lnc
table(k)
gtf_tmp <- gtf_all_trans[k, ]
count_df <- gtf_tmp %>%
  group_by(transcript_type) %>%
  summarise(count = n())
count_df$sample <- rep(">= 2", dim(count_df)[1])
result <- bind_rows(result, count_df)
# more than 5 drivers
k <- gtf_all_trans$ID %in% pan_list[["5"]]$lnc
table(k)
gtf_tmp <- gtf_all_trans[k, ]
count_df <- gtf_tmp %>%
  group_by(transcript_type) %>%
  summarise(count = n())
count_df$sample <- rep(">= 5", dim(count_df)[1])
result <- bind_rows(result, count_df)
# more than 7 drivers
k <- gtf_all_trans$ID %in% pan_list[["7"]]$lnc
table(k)
gtf_tmp <- gtf_all_trans[k, ]
count_df <- gtf_tmp %>%
  group_by(transcript_type) %>%
  summarise(count = n())
count_df$sample <- rep(">= 7", dim(count_df)[1])
result <- bind_rows(result, count_df)

# plot
result$sample <- factor(result$sample, levels = c(
  "Non-drivers",
  "All drivers",
  ">= 2",
  ">= 5",
  ">= 7"
))
p <- ggplot(result, aes(sample, fill = transcript_type)) +
  geom_bar(aes(sample, count, fill = transcript_type),
           width = 0.8,
           stat = "identity", position = "fill", color = "#303030"
  ) +
  scale_fill_manual(values = c(
    "#e3a6a5", "#ecdd9e", "#71aacd",
    "#3e97ac", "#806da5", "#af5366",
    "#de9758", "#518b83", "#e9738d"
  )) +
  labs(
    fill = "Transcript type",
    tittle = "Transcript type of lncRNAs"
  ) +
  ylab("Proportion") +
  xlab("") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 50),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 50),
    axis.title.y = element_text(size = 40, vjust = 3),
    axis.title.x = element_text(size = 50, vjust = 1),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 38)
  ) +
  coord_flip()
# save
ggsave("./output/Transcript_type.tiff",
       height = 7.4, width = 21, dpi = 600
)
