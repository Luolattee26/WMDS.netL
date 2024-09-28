# This code is used to get the NMF subtypes of cancers through drivers
# if u want to reproduce the result, please run the code in the following order
# and remember to change the file path to your own path (and cancer type)



rm(list = ls())
source("./code/utils_TCGAmodel.R")


# input TCGA exp data
raw_data <- read.table("./data/TCGA_14cancer/fpkm/TCGA-UCEC.htseq_fpkm.tsv",
                       header = T, sep = "\t", check.names = F)
fixed_samples <- seprate_samples(raw_data,
                                 data_preprocess = FALSE
)
normal_df <- fixed_samples[[1]]
tumor_df <- fixed_samples[[2]]
rownames(normal_df) <- raw_data[, 1]
rownames(tumor_df) <- raw_data[, 1]
# FPKM to TPM
normal_df <- as.matrix(normal_df)
tumor_df <- as.matrix(tumor_df)

# survival data
# all cancer
file_path <- "./data/survival_data/"
fileName <- dir(file_path)
cancer_list <- fileName2cancerList(fileName)
filePath <- sapply(fileName, function(x) {
  paste(file_path, x, sep = "")
})
survival_data <- lapply(filePath, function(x) {
  read.table(x, header = T, sep = "\t", row.names = 1)
})
names(survival_data) <- cancer_list
survival_LUAD <- survival_data[["UCEC"]]
# select the sample with survival data
key <- rownames(survival_LUAD) %in% colnames(tumor_df)
table(key)
survival_LUAD <- survival_LUAD[key, ]
# os time >= 30
k1 <- survival_LUAD$OS.time >= 30
table(k1)
# no na
k2 <- !(is.na(survival_LUAD$OS.time) | is.na(survival_LUAD$OS))
table(k2)
survival_LUAD <- survival_LUAD[k1 & k2, ]
colnames(survival_LUAD)[c(1, 3)] <- c("event", "time")

# input driver lnc data
pan_7_driver <- read.csv("./data/pancancer_driver/pan_dri_7.csv")
lncRNA_list <- unlist(pan_7_driver$lnc)
all_lnc <- read_excel("./data/WMDS_latest/LNC_name.xlsx", col_names = F)
# get driver lnc list with version
ID_noV <- lapply(rownames(tumor_df), function(x) {
  str_split(x, "\\.")[[1]][1]
})
k <- ID_noV %in% lncRNA_list
lncRNA_list <- rownames(tumor_df)[k]

# prepare NMF input data
k <- intersect(colnames(tumor_df), rownames(survival_LUAD))
nmf.input <- tumor_df[lncRNA_list, k]
# for loop to get best cluster numbers
ranks <- 2:6
path <- "./output/NMF/UCEC/"
f <- paste(path, "/", "rank_for_loop.rdata", sep = "")
random_seed <- 2024
if (!file.exists(f)) {
  result <- nmf(nmf.input, ranks, seed = random_seed)
  save(result, file = f)
}
# visualization the for loop result
load(f)
plot(result)

# the best number is 3
rank <- 3
result_best <- nmf(nmf.input,
                   rank = rank,
                   seed = random_seed
)
# get features
index <- extractFeatures(result_best, "max")
# use the features to cluster again
# and remove the NA feature of a group
nmf.input2 <- nmf.input[na.omit(unlist(index)), ]
result2 <- nmf(nmf.input2,
               rank = rank,
               seed = random_seed
)
# or just use the result of best k
result2 <- result_best

# km clustering
W <- basis(result2)
H <- coef(result2)
kmeans_result <- kmeans(t(H), centers = rank)
table(kmeans_result$cluster)
group <- kmeans_result$cluster

# visualization
# km plot
# match the order
group <- group[match(rownames(survival_LUAD), names(group))]
identical(names(group), rownames(survival_LUAD))
survival_LUAD$group <- group
sfit <- survfit(Surv(time, event) ~ group,
                data = survival_LUAD
)
ggsurvplot(sfit, pval = T, palette = "jco")
# save result
# here just use the result of best k
save(
  list = c("survival_LUAD", "result2", "nmf.input"),
  file = "./output/NMF/UCEC/result_FPKM_3group.Rdata"
)
