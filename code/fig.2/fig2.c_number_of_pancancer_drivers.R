# This code was used to count how many drivers were
# identified in several cancers in the paper,
# termed as "pancancer drivers"
# The final output will be a bar plot



# init code
rm(list = ls())
source("./code/utils.R")


# read drivers
cancer_list <- colnames(read_excel("./data/WMDS_latest/cancer_list.xlsx"))
drivers_14 <- read_excel("./data/WMDS_latest/gene_14.xlsx",
                         col_names = cancer_list)

# calculation
drivers_list <- na.omit(unlist(drivers_14))
element_counts <- data.frame(table(drivers_list))

# data transform
result <- list()
for (i in 1:14) {
  Freq <- dim(subset(element_counts, Freq >= i))[1]
  label <- paste(">=", i, sep = "")
  result[[label]] <- Freq
}
count_df <- data.frame(
  count = unlist(result),
  group = names(result)
)
count_df$group <- factor(count_df$group, levels = lapply(1:14, function(x) {
  paste(">=", x, sep = "")
}))

# plot
p <- ggplot(count_df) +
  geom_bar(aes(group, count),
           fill = "#80b1d3",
           stat = "identity"
  ) +
  xlab("Number of cancer types") +
  ylab("drivers") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 55),
    axis.title.x = element_text(size = 48, vjust = 1,
                                hjust = 0.5),
    axis.text.x = element_text(angle = 30, vjust = 0.5,
                               hjust = 0.5, size = 40, face = "bold"),
    axis.text.y = element_text(size = 52),
    axis.title.y = element_text(size = 53, vjust = 1, hjust = 0.5),
    title = element_text(size = 55)
  ) +
  geom_text(aes(group, count, label = count), vjust = 0,
            position = position_dodge(width = 0.7), size = 15)
p <- p + scale_y_continuous(limits = c(0, 727))
# save
ggsave("./output/drivers_pancancer_overlap.tiff",
       height = 9.5, width = 21, dpi = 600, plot = p
)
