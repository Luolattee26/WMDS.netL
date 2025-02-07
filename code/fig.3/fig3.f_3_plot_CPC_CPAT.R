# This code is used to plot the coding ability of drivers and non-drivers.



# init code
rm(list = ls())


# gtf
gtf_file <- import("./data/gencode.v22.annotation.gtf")
gtf_df <- as.data.frame(gtf_file)
gtf_all_trans <- gtf_df[gtf_df$type == "transcript", ]
Retained_lst <- c(
  "non_coding", "3prime_overlapping_ncrna",
  "antisense", "lincRNA", "retained_intron",
  "sense_intronic",
  "sense_overlapping", "macro_lncRNA"
)
k <- gtf_all_trans$transcript_type %in% Retained_lst
table(k)
gtf_all_trans <- gtf_all_trans[k, ]

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
all_LNC <- read_xlsx("./data/WMDS_latest/LNC_name.xlsx", col_names = F)

# input coding data
cpc <- read.csv("./data/coding/Coding_CPC.txt", sep = "\t")
cpat <- read.csv("./data/coding/cpat.ORF_prob.best.tsv", sep = "\t")

# ID translation
all_ID_V <- gtf_all_trans$gene_id
all_ID <- lapply(all_ID_V, function(x) {
  str_split(x, "\\.")[[1]][1]
})

# pan anal cpc
df_cpc_long <- data.frame()
index <- 0
for (pan in pan_list) {
  index <- index + 1
  lnc_list <- pan$lnc
  k <- all_ID %in% lnc_list
  tmp <- gtf_all_trans[k, ]
  trans_list <- tmp$transcript_id
  k <- cpc$X.ID %in% trans_list
  data <- cpc[k, ]
  data$group <- rep(paste("pandri", names(pan_list)[index], sep = ""), dim(data)[1])
  df_cpc_long <- rbind(df_cpc_long, data)
}
df_cpc_long <- df_cpc_long[, c(1, 7, 8, 9)]

# pan anal cpat
df_cpat_long <- data.frame()
index <- 0
for (pan in pan_list) {
  index <- index + 1
  lnc_list <- pan$lnc
  k <- all_ID %in% lnc_list
  tmp <- gtf_all_trans[k, ]
  trans_list <- tmp$transcript_id
  k <- cpat$seq_ID %in% trans_list
  data <- cpat[k, ]
  data$group <- rep(paste("pandri", names(pan_list)[index], sep = ""), dim(data)[1])
  df_cpat_long <- rbind(df_cpat_long, data)
}
df_cpat_long <- df_cpat_long[, c(1, 11, 12)]

# non-driver
k <- all_LNC$...1 %in% pan_list[[1]]$lnc
table(k)
non_driver <- all_LNC[!k, ]
# cpc
lnc_list <- non_driver$...1
k <- all_ID %in% lnc_list
tmp <- gtf_all_trans[k, ]
trans_list <- tmp$transcript_id
k <- cpc$X.ID %in% trans_list
data <- cpc[k, ]
data$group <- rep("non-driver", dim(data)[1])
data <- data[, c(1, 7, 8, 9)]
df_cpc_long <- rbind(df_cpc_long, data)
# cpat
lnc_list <- non_driver$...1
k <- all_ID %in% lnc_list
tmp <- gtf_all_trans[k, ]
trans_list <- tmp$transcript_id
k <- cpat$seq_ID %in% trans_list
data <- cpat[k, ]
data$group <- rep("non-driver", dim(data)[1])
data <- data[, c(1, 11, 12)]
df_cpat_long <- rbind(df_cpat_long, data)

# plot
# cpc
p <- ggplot(df_cpc_long, aes(
  x = reorder(group, -as.numeric(coding_probability)), y = as.numeric(coding_probability),
  fill = group
)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = list(
      c("non-driver", "pandri1"), c("non-driver", "pandri2"), c("non-driver", "pandri7"),
      c("non-driver", "pandri5")
    ),
    label = "p.format", method = "wilcox.test"
  ) +
  labs(x = "Driver group", y = "coding") +
  theme_minimal()
p
# cpat
p <- ggplot(df_cpat_long, aes(x = reorder(group, -as.numeric(Coding_prob)), y = as.numeric(Coding_prob), fill = group)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = list(
      c("non-driver", "pandri1"), c("non-driver", "pandri2"), c("non-driver", "pandri7"),
      c("non-driver", "pandri5")
    ),
    label = "p.format", method = "wilcox.test"
  ) +
  labs(x = "Driver group", y = "coding") +
  theme_minimal()
p

# coding percent
threshold_coding <- 0.364
label <- names(table(df_cpc_long$group))
df_cpat_long$label <- ifelse(df_cpat_long$Coding_prob > threshold_coding, "coding", "noncoding")
cpc_percent <- list()
cpat_percent <- list()
index <- 0
for (group in names(table(df_cpc_long$group))) {
  print(group)
  index <- index + 1
  tmp1 <- df_cpc_long[df_cpc_long$group == group, ]
  per1 <- round(table(tmp1$label)["coding"] / dim(tmp1)[1] * 100, 2)
  cpc_percent[index] <- per1

  tmp2 <- df_cpat_long[df_cpat_long$group == group, ]
  per2 <- round(table(tmp2$label)["coding"] / dim(tmp2)[1] * 100, 2)
  cpat_percent[index] <- per2
}
df_per <- cbind(cpc_percent, cpat_percent, names(table(df_cpc_long$group)))
colnames(df_per) <- c("cpc", "cpat", "group")
df_per <- as.data.frame(df_per)
df_per$mean_per <- round((as.numeric(df_per$cpc) + as.numeric(df_per$cpat)) / 2, 2)
df_per$cpc <- unlist(df_per$cpc)
df_per$cpat <- unlist(df_per$cpat)
df_per$group <- unlist(df_per$group)
df_per <- df_per[order(df_per$mean_per), ]
rank <- df_per$group
df_per <- tidyr::pivot_longer(df_per, cols = c(cpc, cpat, mean_per), names_to = "score_group", values_to = "per")
df_per$group <- factor(df_per$group, levels = rank)
# remove mean per
df_per <- subset(df_per, score_group != "mean_per")

# mixed plot
p <- ggplot(df_per, aes(x = group, y = per, fill = score_group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +
  labs(x = "drivers group", y = "Conding percent", fill = "Algorithms", tittle = "Coding ability of drivers") +
  scale_fill_manual(values = c("#ffc24b", "#1d3557"), labels = c("CPAT", "CPC")) +
  theme_classic() +
  scale_x_discrete(labels = c(
    "Non-drivers", "All drivers", ">= 2", ">= 5",
    ">= 7"
  )) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 55),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 30, vjust = 0.6, hjust = 0.5, size = 53, face = "bold"),
    axis.text.y = element_text(size = 59),
    axis.title.y = element_text(size = 60, vjust = 0.5, hjust = 0.5),
    legend.text = element_text(size = 50, vjust = 1),
    legend.title = element_text(size = 52, vjust = 1),
    title = element_text(size = 55)
  )
# ggtitle('Coding ability of drivers')
p
# save
ggsave("./output/coding.tiff",
       height = 9, width = 21, dpi = 600
)
