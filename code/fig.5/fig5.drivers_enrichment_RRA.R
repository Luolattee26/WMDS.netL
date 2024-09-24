# This code is used to get the enriched GO and KEGG pathways
# for the important drivers selected by RRA algorithm.



# init code
rm(list = ls())
source("./code/utils.R")


# input lnc and mRNA
lnc_mRNA <- read_excel("./data/WMDS_latest/final_network_nodes_edges.xlsx",
                       col_names = c("node1", "node2", "pvalue")
)
pan_7 <- read.csv("./data/pancancer_driver/pan_dri_7.csv")
pan_7_list <- pan_7$lnc

# whether to use full list of drivers
load("./data/important_lnc_in_cox.Rdata")
useFull <- F
if (!useFull) {
  k <- pan_7_list %in% lst_ID
  table(k)
  pan_7_list <- pan_7_list[k]
}

# get related mRNA for each lnc in pan_7_list
res_list <- vector("list")
for (lnc in pan_7_list) {
  tmp_lnc_mRNA <- lnc_mRNA[lnc_mRNA$node1 == lnc, ]
  res_list[[lnc]] <- tmp_lnc_mRNA[order(tmp_lnc_mRNA$pvalue), ]$node2
}
# get freq
freq <- as.data.frame(table(unlist(res_list)))
# RRA algorithms and add freq
ag <- aggregateRanks(res_list)
ag$Freq <- freq[match(ag$Name, freq$Var1), 2]
# select the important mRNA
useImportant <- TRUE
if (useImportant) {
  important_mRNA <- subset(ag, Freq >= 2 | Score < 1)
  dim(important_mRNA)
} else {
  important_mRNA <- ag
}

# ID translation
ID <- bitr(important_mRNA$Name,
           fromType = "ENSEMBL",
           toType = "ENTREZID",
           OrgDb = org.Hs.eg.db
)
non_duplicates_idx <- which(!(duplicated(ID$ENSEMBL)))
ID <- ID[non_duplicates_idx, ]

# enrichment
# GO
go_experiment <- enrichGO(
  gene = ID$ENTREZID, OrgDb = "org.Hs.eg.db",
  ont = "all", readable = T, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1
)
go_result_experiment <- as.data.frame(go_experiment@result)
goBP <- subset(go_result_experiment, subset = (ONTOLOGY == "BP"))[1:100, ]
goCC <- subset(go_result_experiment, subset = (ONTOLOGY == "CC"))[1:100, ]
goMF <- subset(go_result_experiment, subset = (ONTOLOGY == "MF"))[1:100, ]
# KEGG
KEGG_experiment <- enrichKEGG(
  gene = ID$ENTREZID, organism = "hsa",
  keyType = "kegg", pAdjustMethod = "BH", pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)
KEGG_result_experiment <- setReadable(KEGG_experiment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
KEGG_result_experiment <- data.frame(KEGG_result_experiment@result)
kegg <- KEGG_result_experiment[1:100, ]


# save enrichment result
save(
  list = c("goBP", "goCC", "goMF", "kegg"),
  file = "./output/enrichment/important_drivers_RRAselected_enrichment.Rdata"
)
