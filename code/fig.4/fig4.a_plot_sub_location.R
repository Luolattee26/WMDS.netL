# This code is used to plot the subcellular localization of lncRNA
# in HepG2 and K562 cell lines with stack bar plot.



# init code
rm(list = ls())
source("./code/utils_TCGAmodel.R")


# read data
load("./data/ENCODE_SUBCELLULAR_LOCALIZATION/transcript_quantifications/data_cleaned.rdata")
rownames(df_HepG2) <- c(
  "Non-drivers", "All drivers", ">= 2",
  ">= 5", ">= 7"
)
colnames(df_HepG2) <- lapply(colnames(df_HepG2), function(x) {
  return(str_split_1(x, "_")[2])
})
colnames(df_HepG2)[2] <- "insoluble-\ncytoplasmic"
rownames(df_K562) <- c(
  "Non-drivers", "All drivers", ">= 2",
  ">= 5", ">= 7"
)
colnames(df_K562) <- lapply(colnames(df_K562), function(x) {
  return(str_split_1(x, "_")[2])
})
colnames(df_K562)[2] <- "insoluble-\ncytoplasmic"
# data transform
df_long_hepg2 <- df_HepG2 %>%
  rownames_to_column(var = "sample") %>%
  gather(key = "type", value = "value", -sample)
df_long_K562 <- df_K562 %>%
  rownames_to_column(var = "sample") %>%
  gather(key = "type", value = "value", -sample)

# plot
# HepG2
df_long_hepg2$sample <- factor(df_long_hepg2$sample, levels = c(
  "Non-drivers",
  "All drivers",
  ">= 2",
  ">= 5",
  ">= 7"
))
p <- ggplot(df_long_hepg2, aes(sample, fill = type)) +
  geom_bar(aes(sample, value, fill = type),
           width = 0.8,
           stat = "identity", position = "fill", color = "#303030"
  ) +
  scale_fill_manual(values = c(
    "#e3a6a5", "#ecdd9e", "#71aacd",
    "#3e97ac", "#806da5", "#af5366",
    "#de9758", "#518b83", "#e9738d"
  )) +
  labs(
    fill = "Location",
    tittle = "Transcript type of lncRNAs"
  ) +
  ylab("Proportion") +
  xlab("") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 50),
    axis.text.x = element_text(size = 42, angle = 30, vjust = 0.6),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 42),
    axis.title.y = element_text(size = 42, vjust = 1, hjust = 0.5),
    axis.title.x = element_text(size = 42, vjust = 1),
    legend.text = element_text(size = 32),
    legend.title = element_text(size = 34),
    legend.position = "right",
    title = element_text(size = 50),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.5, "cm")
  ) +
  ggtitle("HepG2") +
  guides(fill = guide_legend(nrow = 4))
# save
ggsave("./output/RCI/pancancer/HepG2_ENCODE.tiff",
       height = 11.5, width = 10, dpi = 600
)
# K562
df_long_K562$sample <- factor(df_long_K562$sample, levels = c(
  "Non-drivers",
  "All drivers",
  ">= 2",
  ">= 5",
  ">= 7"
))
p <- ggplot(df_long_K562, aes(sample, fill = type)) +
  geom_bar(aes(sample, value, fill = type),
           width = 0.8,
           stat = "identity", position = "fill", color = "#303030"
  ) +
  scale_fill_manual(values = c(
    "#e3a6a5", "#ecdd9e", "#71aacd",
    "#3e97ac", "#806da5", "#af5366",
    "#de9758", "#518b83", "#e9738d"
  )) +
  labs(
    fill = "Location",
    tittle = "Transcript type of lncRNAs"
  ) +
  ylab("Proportion") +
  xlab("") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 50),
    axis.text.x = element_text(size = 42, angle = 30, vjust = 0.6),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 42),
    axis.title.y = element_text(size = 42, vjust = 1, hjust = 0.5),
    axis.title.x = element_text(size = 42, vjust = 1),
    legend.text = element_text(size = 32),
    legend.title = element_text(size = 34),
    legend.position = "right",
    title = element_text(size = 50),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.5, "cm")
  ) +
  ggtitle("K562") +
  guides(fill = guide_legend(nrow = 4))
# save
ggsave("./output/RCI/pancancer/K562_ENCODE.tiff",
       height = 11.5, width = 10, dpi = 600
)
