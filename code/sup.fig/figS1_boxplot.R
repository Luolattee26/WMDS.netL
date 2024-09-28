# This code is used to plot the expression
# of drivers and non-drivers in 14 cancers
# The out put is a boxplot with breaked y-axis



# init code
rm(list = ls())
source("./code/utils_TCGAmodel.R")
driver_type <- 7
output_path <- "./output/exp_plot/pandri7/"
# "the pattern supported:"
# "1--drivers in tumor vs. non-drivers in tumor--Fig.4C"
# "2--drivers in normal vs. non-drivers in normal--Fig.S1A"
# "3--drivers in tumor vs. drivers in normal--Fig.S1B"
# "4--non-drivers in tumor vs. non-drivers in normal--Fig.S1C"
pattern <- 2


# input driver lnc data
driver <- paste("./data/pancancer_driver/pan_dri_", 7, ".csv", sep = "")
pan_7_driver <- read.csv(driver)
lncRNA_list <- unlist(pan_7_driver$lnc)

# prepare violin data
# This data is get from UCSC-Xena platform
data <- violin_data("./data/TCGA_14cancer/fpkm/",
                    driver_list = lncRNA_list,
                    pattern = pattern
)
colnames(data) <- c("expr", "group", "project")

# order data
data_new <- data %>%
  group_by(project) %>%
  mutate(median = median(expr), group_max = max(expr)) %>%
  arrange(desc(median))
# factor
data_new$project <- factor(data_new$project, levels = unique(data_new$project))
data_new1 <- data_new %>%
  group_by(project) %>%
  mutate(group_number = length(unique(group)))
# stat. test
stat_data <- data_new1[which(data_new1$group_number == 2), ] %>%
  group_by(project) %>%
  wilcox_test(expr ~ group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p")
data_new1 <- unique(data_new1[, c(3, 5, 6)])
index <- which(data_new1$group_number == 2)
bracket_data <- data.frame(
  x = index - 0.3,
  xend = index + 0.3,
  y = as.numeric(unlist(data_new1[index, 2])) + 0.5,
  annotation = ""
)

# new plot
p <- ggplot() +
  geom_boxplot(data = data_new, aes(x = project, y = expr, fill = group),
               outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(name = NULL,
                    values = c("drivers" = "#d6503a",
                               "non-drivers" = "#5488ef")) +
  theme_bw() +
  ylab("Log2(TPM + 1)") +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 14),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # legend.position = c(0.99, 0.99),
    legend.position = "none",
    legend.justification = c(1, 1)
  ) +
  geom_signif(
    stat = "identity",
    position = "identity",
    data = bracket_data,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = y,
      annotation = annotation
    )
  ) +
  annotate(geom = "text", x = index, y = bracket_data$y + 0.5,
           label = as.character(stat_data$p.signif)) +
  scale_y_continuous(limits = c(0, 20)) +
  scale_y_break(breaks = c(3, 14), scales = 0.3, space = 0.3,
                ticklabels = c(14, 17)) +
  scale_x_discrete(labels = sort(as.character(unique(data_new$project)))) +
  ggtitle("Expression of >=7 drivers and non-drivers in tumor samples")
# save
ggsave(paste(output_path, "/", "pan_cancer_boxplot_breaked_pattern_",
             pattern, '.tiff', sep = ""),
       height = 9, width = 21, dpi = 600, units = "cm"
)
