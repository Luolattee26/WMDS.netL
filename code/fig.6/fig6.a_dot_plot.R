# init code
rm(list = ls())
source("./code/utils.R")
source("./code/utils_TCGAmodel.R")
library(ggplot2)
library(reshape2)


# data input
cox_data <- read.csv("./output/pan_cancer_cox/cox_p.csv", row.names = 1)
HR_data <- read.csv("./output/pan_cancer_cox/cox_HR.csv", row.names = 1)
# long data prepare
cox_data <- t(as.matrix(cox_data))
HR_data <- t(as.matrix(HR_data))
HR_data_long <- melt(HR_data, varnames = c("type", "gene"),
                     value.name = "HR")
cox_data_long <- melt(cox_data, varnames = c("type", "gene"),
                      value.name = "pvalue")
mixed_HRcox <- cbind(HR_data_long, cox_data_long$pvalue)
colnames(mixed_HRcox) <- c("type", "gene", "HR", "pvalue")
# whether to use full list of drivers
load("./data/important_lnc_in_cox.Rdata")
useFull <- FALSE
if (!useFull) {
  k <- mixed_HRcox$gene %in% names(lst)
  table(k)
  mixed_HRcox <- mixed_HRcox[k, ]
}

# plot
p1 <- ggplot(mixed_HRcox, aes(x = type, y = gene)) +
  geom_point(aes(size = -log10(pvalue + 0.0001), fill = HR),
             shape = 21,
             color = "black"
  ) +
  scale_fill_gradient2(
    name = "HR",
    limit = c(0, 4),
    breaks = c(0, 1, 4),
    low = "#444283",
    high = "#943934",
    mid = "white",
    midpoint = 1,
    guide = "colourbar"
  ) +
  scale_size_continuous(
    name = "-Log10 p",
    limits = c(-0.001, 10),
    breaks = c(0, 2, 4, 8, 10)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.y = element_text(angle = 0, size = 8),
    axis.text.x = element_text(angle = 45, size = 14, hjust = 0, vjust = 0.5)
  ) +
  coord_flip() +
  scale_y_discrete(position = "right")


# save
ggsave(
  filename = "./output/pan_cancer_cox/dot_plot_sig.tiff",
  width = 20, height = 11,
  dpi = 600, plot = p1, units = "cm"
)
