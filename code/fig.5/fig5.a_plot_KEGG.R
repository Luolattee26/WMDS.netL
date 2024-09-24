# This code is to plot the enriched KEGG pathways



# init code
rm(list = ls())
source("./code/utils.R")


# read enrichment data
load("./output/enrichment/important_drivers_RRAselected_enrichment.Rdata")
# data prepare
tmp_data <- kegg[1:25, ]
data <- data.frame()
for (type in names(table(tmp_data$category))) {
  data <- rbind(
    data,
    c(type, NA, NA),
    tmp_data[tmp_data$category == type, c(4, 11, 8)]
  )
}
colnames(data) <- c("Pathway", "Number", "P.adj")
data$Number <- as.numeric(data$Number)
data$P.adj <- as.numeric(data$P.adj)
data$P.adj <- format(data$P.adj, scientific = TRUE, digits = 2)


# plot
# basic plot
ggplot(data) +
  geom_bar(aes(Pathway, Number), stat = "identity") +
  coord_flip()
# add group info
na_index <- which(is.na(data$Number))
group <- c(
  NA, rep(data$Pathway[na_index[1]], na_index[2] - na_index[1] - 1),
  NA, rep(data$Pathway[na_index[2]], na_index[3] - na_index[2] - 1),
  NA, rep(data$Pathway[na_index[3]], na_index[4] - na_index[3] - 1),
  NA, rep(data$Pathway[na_index[4]], na_index[5] - na_index[4] - 1),
  NA, rep(data$Pathway[na_index[5]], nrow(data) - na_index[5])
)
data$Group <- group
data$Pathway <- factor(data$Pathway, levels = rev(data$Pathway))
# add color and face
colors <- c(
  "black", rep("#9dd1c9", na_index[2] - na_index[1] - 1),
  "black", rep("#f2b06f", na_index[3] - na_index[2] - 1),
  "black", rep("#bebbd7", na_index[4] - na_index[3] - 1),
  "black", rep("#eb8776", na_index[5] - na_index[4] - 1),
  "black", rep("#88afcf", nrow(data) - na_index[5])
)
face <- c(
  "bold", rep(NULL, na_index[2] - na_index[1] - 1),
  "bold", rep(NULL, na_index[3] - na_index[2] - 1),
  "bold", rep(NULL, na_index[4] - na_index[3] - 1),
  "bold", rep(NULL, na_index[5] - na_index[4] - 1),
  "bold", rep(NULL, nrow(data) - na_index[5])
)
# final plot
p <- ggplot(data, aes(Pathway, Number)) +
  geom_bar(aes(fill = Group), stat = "identity") +
  geom_text(aes(label = P.adj, y = Number - 5), size = 10) +
  scale_fill_manual(values = c("#9dd1c9", "#f2b06f", "#bebbd7", "#eb8776", "#88afcf")) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 85, 20)) +
  ylab("Number") +
  xlab("") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_line(linetype = "dashed"),
    panel.grid.major.y = element_line(linetype = "dashed"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(
      face = "bold",
      color = rev(colors),
      hjust = 0,
      size = 30, lineheight = 2
    ),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    plot.title = element_text(hjust = 0.5, size = 30),
    text = element_text(family = "Times")
  ) +
  ggtitle("KEGG annotation")


# save
ggsave("./output/enrichment/kegg_mixed_barplot.tiff",
       height = 14.5, width = 21, dpi = 1000
)
