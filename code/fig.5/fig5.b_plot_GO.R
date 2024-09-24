# This code is to plot the enriched GO terms



# init code
rm(list = ls())
source("./code/utils.R")


# read enrichment data
load("./output/enrichment/important_drivers_RRAselected_enrichment.Rdata")
# data prepare
data <- rbind(
  c("BP", NA, NA),
  na.omit(goBP)[1:10, c(3, 10, 7)]
)
data <- rbind(
  data,
  c("MF", NA, NA),
  na.omit(goMF)[1:5, c(3, 10, 7)]
)
data <- rbind(
  data,
  c("CC", NA, NA),
  na.omit(goCC)[1:5, c(3, 10, 7)]
)
colnames(data) <- c("Pathway", "Number", "P.adj")
data$Number <- as.numeric(data$Number)
data$P.adj <- as.numeric(data$P.adj)
# re-write
data[15, 1] <- "DNA-binding transcription activator activity(II-specific)"
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
  NA, rep(data$Pathway[na_index[3]], nrow(data) - na_index[3])
)
data$Group <- group
data$Pathway <- factor(data$Pathway, levels = rev(data$Pathway))
# add color and face
colors <- c(
  "black", rep("#9dd1c9", na_index[2] - na_index[1] - 1),
  "black", rep("#f2b06f", na_index[3] - na_index[2] - 1),
  "black", rep("#bebbd7", nrow(data) - na_index[3])
)
face <- c(
  "bold", rep(NULL, na_index[2] - na_index[1] - 1),
  "bold", rep(NULL, na_index[3] - na_index[2] - 1),
  "bold", rep(NULL, nrow(data) - na_index[3])
)
# final plot
p <- ggplot(data, aes(Pathway, Number)) +
  geom_bar(aes(fill = Group), stat = "identity") +
  geom_text(aes(label = P.adj, y = Number - 5), size = 10) +
  scale_fill_manual(values = c("#9dd1c9", "#bebbd7", "#f2b06f")) +
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
  ggtitle("GO annotation")
# p + scale_x_discrete(labels = function(x) str_replace_all(x, ",.*", ""))


# save
ggsave("./output/enrichment/go_mixed_barplot.tiff",
       height = 14.5, width = 21, dpi = 1000
)
