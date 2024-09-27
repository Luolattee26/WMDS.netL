# init code
rm(list = ls())
source("./code/utils.R")
output_path <- "./output/"


# data import and prepare
TS_df <- pandriTS_data_df(TS_result_path = "./data/pandri_TS_scores/")
phastCon_df <- phastCon_data_df(input_path = "./data/Conservation_Scores/")
colnames(TS_df) <- c("Score", "group")
colnames(phastCon_df) <- c("Score", "group")
TS_df <- dplyr::select(TS_df, group, Score)
phastCon_df <- dplyr::select(phastCon_df, group, Score)
TS_df$group <- factor(TS_df$group, levels = c(
  "all non-driver", "all drivers", "more than 2",
  "more than 5", "more than 7"
))
phastCon_df$group <- factor(phastCon_df$group, levels = c(
  "all non-driver", "all drivers", "more than 2",
  "more than 5", "more than 7"
))
phastCon_df$Score <- as.numeric(phastCon_df$Score)
TS_df$Score <- as.numeric(TS_df$Score)

# statistic test
stat.test.phastCon <- phastCon_df %>%
  wilcox_test(
    Score ~ group,
    p.adjust.method = "BH",
    detailed = T
  )
stat.test.TS <- TS_df %>%
  wilcox_test(
    Score ~ group,
    p.adjust.method = "BH",
    detailed = T
  )

# plot of phastCon
colors <- c(
  "#eb4b3a", "#48bad0", "#1a9781",
  "#355783", "#ef9a80"
)
p <- ggplot(phastCon_df) +
  geom_boxplot(aes(group, Score, color = group),
               outlier.shape = 21
  ) +
  scale_color_manual(values = c(
    "#eb4b3a", "#48bad0", "#1a9781",
    "#355783", "#ef9a80"
  )) +
  scale_x_discrete(labels = c("Non-drivers", "All drivers", ">= 2", ">= 5", ">= 7")) +
  xlab("") +
  ylab("Sequence conservation") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 30, vjust = 0.8, face = "bold",
      color = colors, size = 13, hjust = 0.65
    ),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 21)
  ) +
  # ggtitle('Sequence conservation') +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 1), limits = c(0, 1.1))
# add stat.test result
x_value <- c(2, 3, 4, 5, 3, 4, 5, 4, 5, 5)
y_value_data <- unlist(lapply(levels(phastCon_df$group), function(x) {
  return(max(phastCon_df[phastCon_df$group == x, 2]))
}))[2:5]
names(y_value_data) <- levels(phastCon_df$group)[2:5]
y_value <- rep(y_value_data, 1:4) + 0.01
y_value <- y_value + c(0.015, 0.015 * 1:2, 0.015 * 1:3, 0.015 * 1:4)
color_value <- c(
  colors[1], colors[1], colors[1], colors[1],
  colors[2], colors[2], colors[2],
  colors[3], colors[3],
  colors[4]
)
for (i in 1:nrow(stat.test.phastCon)) {
  if (stat.test.phastCon$p.adj.signif[i] != "ns") {
    y_tmp <- y_value[i]
    p <- p + annotate(
      geom = "text",
      label = stat.test.phastCon$p.adj.signif[i],
      x = x_value[i],
      y = y_tmp,
      color = color_value[i],
      size = 6
    )
  }
}
ggsave(paste(output_path, "/", "fancy_conservation.tiff", sep = ""),
       height = 7.4, width = 8, dpi = 600, units = "cm"
)

# plot of TS
colors <- c(
  "#eb4b3a", "#48bad0", "#1a9781",
  "#355783", "#ef9a80"
)
p <- ggplot(TS_df) +
  geom_boxplot(aes(group, Score, color = group),
               outlier.shape = 21
  ) +
  scale_color_manual(values = c(
    "#eb4b3a", "#48bad0", "#1a9781",
    "#355783", "#ef9a80"
  )) +
  scale_x_discrete(labels = c("Non-drivers", "All drivers", ">= 2", ">= 5", ">= 7")) +
  xlab("") +
  ylab("Tissue specificity") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 30, vjust = 0.8, face = "bold",
      color = colors, size = 13, hjust = 0.65
    ),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 21)
  ) +
  # ggtitle('Tissue specific') +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 1), limits = c(0, 1))
# add stat.test result
x_value <- c(2, 3, 4, 5, 3, 4, 5, 4, 5, 5)
y_value_data <- unlist(lapply(levels(TS_df$group), function(x) {
  return(max(TS_df[TS_df$group == x, 2]))
}))[2:5]
names(y_value_data) <- levels(TS_df$group)[2:5]
y_value <- rep(y_value_data, 4:1) + 0.0005
y_value <- y_value + c(0.015 * 1:4, 0.015 * 1:3, 0.015 * 1:2, 0.015)
# # use y_value_new to limit the height
# y_value_new <- c(0.55, 0.55, 0.55, 0.55)
# names(y_value_new) <- levels(TS_df$group)[1:4]
# y_value <- rep(y_value_new, 4:1) + 0.0005
# y_value <- y_value + c(0.015 * 1:4, 0.015 * 1:3, 0.015 * 1:2, 0.015)
color_value <- c(
  colors[1], colors[1], colors[1], colors[1],
  colors[2], colors[2], colors[2],
  colors[3], colors[3],
  colors[4]
)
for (i in 1:nrow(stat.test.TS)) {
  if (stat.test.TS$p.adj.signif[i] != "ns") {
    y_tmp <- y_value[i]
    p <- p + annotate(
      geom = "text",
      label = stat.test.TS$p.adj.signif[i],
      x = x_value[i],
      y = y_tmp,
      color = color_value[i],
      size = 6
    )
  }
}
ggsave(paste(output_path, "/", "fancy_TS.tiff", sep = ""),
       height = 7.4, width = 8, dpi = 600, units = "cm"
)
