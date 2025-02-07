# This code is used to generate the figure 4b in the manuscript,
# to show the interactions of lncRNAs in different driver groups.



# init code
rm(list = ls())
source("./code/utils_TCGAmodel.R")


# interaction data
# These data are from starBase
miRNA2lnc <- read.csv("./data/starBase_data/starBase_miRNA_lnc.txt", sep = "\t")
RBP <- read.csv("./data/starBase_data/starBase_RBP_lnc.txt", sep = "\t")

# input lncDATA
file_path <- "./data/pancancer_driver/"
index <- 0
pan_list <- list()
pancancer_driver_list <- c(1, 2, 5, 7)
for (i in pancancer_driver_list) {
  index <- index + 1
  fileName <- paste(file_path, "pan_dri_", i, ".csv", sep = "")
  pan_list[[index]] <- read.csv(fileName, header = T, sep = ",",
                                check.names = F)
}
names(pan_list) <- pancancer_driver_list
all_lnc <- read_excel("./data/WMDS_latest/LNC_name.xlsx", col_names = F)
all_lnc_list <- all_lnc$...1

# count
# miRNA
mir_count <- list()
RBP_count <- list()
fractionMir <- list()
fractionRbp <- list()
index <- 0
for (pan in pan_list) {
  index <- index + 1
  k1 <- miRNA2lnc$geneID %in% pan$lnc
  mirTmp <- miRNA2lnc[k1, ]
  fractionMir[index] <- round(length(levels(factor(mirTmp$geneID))) / dim(pan)[1] * 100, 2)
  mir_count[index] <- dim(mirTmp)[1]
  k2 <- RBP$geneID %in% pan$lnc
  RbpTmp <- RBP[k2, ]
  fractionRbp[index] <- round(length(levels(factor(RbpTmp$geneID))) / dim(pan)[1] * 100, 2)
  RBP_count[index] <- dim(RbpTmp)[1]
}
df <- data.frame(
  group = names(pan_list), mirCount = log10(unlist(mir_count)), rbpCount = log10(unlist(RBP_count)),
  mirPercent = unlist(fractionMir), rbpPercent = unlist(fractionRbp)
)
# add non-driver
k <- all_lnc_list %in% pan_list[[1]]$lnc
all_non <- all_lnc_list[!k]
k1 <- miRNA2lnc$geneID %in% all_non
mirTmp <- miRNA2lnc[k1, ]
f1 <- round(length(levels(factor(mirTmp$geneID))) / length(all_non) * 100, 2)
c1 <- dim(mirTmp)[1]
k2 <- RBP$geneID %in% all_non
RbpTmp <- RBP[k1, ]
f2 <- round(length(levels(factor(RbpTmp$geneID))) / length(all_non) * 100, 2)
c2 <- dim(RbpTmp)[1]
tmp <- list("non", c1, c2, f1, f2)
df <- rbind(df, tmp)

# plot
df_noCount <- df[, c(1, 4, 5)]
df_long <- tidyr::pivot_longer(df_noCount, cols = c(mirPercent, rbpPercent), names_to = "interGroup", values_to = "value")
df_long$group <- factor(df_long$group, levels = c("non", pancancer_driver_list))
df_long$value <- round(as.numeric(df_long$value), 1)
p <- ggplot(df_long, aes(x = group, y = value, fill = interGroup)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    x = "drivers group", y = "Prop. of lncRNAs with interactions", fill = "Type of interactions",
    tittle = "Interactions of drivers"
  ) +
  scale_fill_manual(values = c("#f6b47b", "#249fca"), labels = c("miRNA-lncRNA", "RNA binding protein")) +
  geom_text(aes(label = value), vjust = -0.3, position = position_dodge(width = 0.7), size = 13) +
  theme_classic() +
  scale_x_discrete(labels = c(
    "Non-drivers", "All drivers", ">= 2", ">= 5",
    ">= 7"
  )) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 50),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 38, face = "bold"),
    # axis.text.x = element_blank(),
    axis.text.y = element_text(size = 40),
    axis.title.y = element_text(size = 43, vjust = 0.5),
    legend.text = element_text(size = 35),
    legend.title = element_text(size = 33),
    title = element_text(size = 50)
  )
p <- p + scale_y_break(breaks = c(15, 70), scales = 0.5, space = 0.3, ticklabels = c(75, 90)) +
  scale_y_continuous(limits = c(0, 100))
#  save
ggsave("./output/interactions.tiff",
       height = 9, width = 21, dpi = 600
)
