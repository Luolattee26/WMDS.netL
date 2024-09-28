# This code is used to viz the NMF components of cancers
# if u want to reproduce the result, please run the code in the following order
# and remember to change the file path to your own path (and cancer type)



# init code
rm(list = ls())
source("./code/utils_TCGAmodel.R")
out_path <- "./output/NMF/UCEC/"


# read data
load("./output/NMF/UCEC/result_FPKM_3group.Rdata")
# get group info
group_info_reordered <- survival_LUAD$group[match(rownames(survival_LUAD), colnames(nmf.input))]

# components heatmap
W <- basis(result2)
# normalization
W <- t(apply(W, 1, scale))
colnames(W) <- paste("sample group", 1:ncol(W))
group_info_reordered <- survival_LUAD$group[match(rownames(survival_LUAD), colnames(W))]
tiff(
  filename = "./output/NMF/UCEC/heatmap_components.tiff",
  width = 20, height = 20, units = "cm",
  res = 600
)
ht_opt$TITLE_PADDING <- unit(c(4, 4), "points")
hr <- Heatmap(W,
              column_split = ncol(W),
              show_row_names = F, show_column_names = T,
              cluster_rows = T,
              heatmap_legend_param = list(
                title = "normalized basis",
                title_gp = gpar(fontsize = 17),
                labels_gp = gpar(fontsize = 15)
              ),
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize = 15),
              column_title_gp = gpar(fill = c("red", "blue"), font = 1:2)
)
draw(hr, column_title = "NMF basis of drivers in UCEC",
     column_title_gp = gpar(fontsize = 25))
dev.off()
