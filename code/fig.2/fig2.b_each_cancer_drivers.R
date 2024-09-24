# This code was used to count how many drivers were
# identified for each of the 14 cancers in the paper
# The final output will be a bar plot



# init code
rm(list = ls())
source("./code/utils.R")


# read drivers
cancer_list <- colnames(read_excel("./data/WMDS_latest/cancer_list.xlsx"))
drivers_14 <- read_excel("./data/WMDS_latest/gene_14.xlsx",
                         col_names = cancer_list)
# get driver number
driver_count <- unlist(lapply(cancer_list, function(x) {
  return(length(unlist(na.omit(drivers_14[x]))))
}))
count_df <- data.frame(
  count = driver_count,
  cancer = cancer_list
)

# plot
p <- ggplot(count_df) +
  geom_bar(aes(cancer, count),
           fill = "#fb8073",
           stat = "identity"
  ) +
  xlab("cancers") +
  ylab("drivers") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 55),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 30, vjust = 0.5,
                               hjust = 0.5, size = 40, face = "bold"),
    axis.text.y = element_text(size = 44),
    axis.title.y = element_text(size = 53, vjust = 1, hjust = 0.5),
    title = element_text(size = 55)
  ) +
  geom_text(aes(cancer, count, label = count), vjust = 0,
            position = position_dodge(width = 0.7), size = 15)
p <- p + scale_y_continuous(limits = c(0, 220))



# save
ggsave("./output/drivers_14cancers.tiff",
       height = 9.5, width = 21, dpi = 600, plot = p
)
