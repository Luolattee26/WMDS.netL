# This script is used to generate the figure 3b in the main text.
# u can get the length of the transcript from the gtf file.



# init code
rm(list = ls())
source("./code/utils.R")
output_path <- "./output/"
filePath <- "./data/pancancer_driver/"


# input gtf file
gtf_file <- import("./data/gencode.v22.annotation.gtf")
gtf_df <- as.data.frame(gtf_file)

# import lncDATA
fileName1 <- paste(filePath, "pan_dri_1.csv", sep = "")
fileName2 <- paste(filePath, "pan_dri_2.csv", sep = "")
fileName5 <- paste(filePath, "pan_dri_5.csv", sep = "")
fileName7 <- paste(filePath, "pan_dri_7.csv", sep = "")
# lnc lists
# >= 7
pan_7_driver <- read.csv(fileName7)
pan_7_driver_list <- unlist(pan_7_driver$lnc)
# all lncs
all_lnc <- read_excel("./data/WMDS_latest/LNC_name.xlsx", col_names = F)
all_lnc_list <- unlist(all_lnc$...1)
# >= 2
pan_2_driver <- read.csv(fileName2)
pan_2_driver_list <- pan_2_driver$lnc
# >= 1
all_driver <- read.csv(fileName1)
all_driver_list <- all_driver$lnc
# >= 5
pan_5_driver <- read.csv(fileName5)
pan_5_driver_list <- unlist(pan_5_driver$lnc)

# transcript length
# delete the transcript which is not lncRNA
gtf_all_trans <- gtf_df[gtf_df$type == "transcript", ]
Retained_lst <- c(
  "non_coding", "3prime_overlapping_ncrna",
  "antisense", "lincRNA", "retained_intron",
  "sense_intronic",
  "sense_overlapping", "macro_lncRNA"
)
k <- gtf_all_trans$transcript_type %in% Retained_lst
table(k)
gtf_all_trans <- gtf_all_trans[k, ]
gtf_all_exons <- gtf_df[gtf_df$type == "exon", ]
Retained_lst <- c(
  "non_coding", "3prime_overlapping_ncrna",
  "antisense", "lincRNA", "retained_intron",
  "sense_intronic",
  "sense_overlapping", "macro_lncRNA"
)
k <- gtf_all_exons$transcript_type %in% Retained_lst
table(k)
gtf_all_exons <- gtf_all_exons[k, ]
# ID translation
ID_list <- unlist(lapply(gtf_all_trans$gene_id, function(x) {
  str_split_1(x, "\\.")[[1]]
}))

# get driver gtf data
# >= 7
k <- ID_list %in% pan_7_driver_list
gtf_7_driver <- gtf_all_trans[k, ]
# >= 5
k <- ID_list %in% pan_5_driver_list
gtf_5_driver <- gtf_all_trans[k, ]
# >= 2
k <- ID_list %in% pan_2_driver_list
gtf_2_driver <- gtf_all_trans[k, ]
# >= 1
k <- ID_list %in% all_driver_list
gtf_1_driver <- gtf_all_trans[k, ]
# non-drivers
k <- all_lnc_list %in% all_driver_list
table(k)
all_non_list <- all_lnc_list[!k]
k <- ID_list %in% all_non_list
gtf_non_driver <- gtf_all_trans[k, ]

# length calculation
# first get a list
cal_list <- list(gtf_1_driver, gtf_2_driver, gtf_5_driver, gtf_7_driver, gtf_non_driver)
cal_res <- list()
index <- 0
for (single_type in cal_list) {
  index <- index + 1
  trans_length <- lapply(single_type$transcript_id, function(x) {
    len <- sum(gtf_all_exons[gtf_all_exons$transcript_id == x, ]$width)
    return(len)
  })
  cal_res[index] <- list(trans_length)
}
# data prepare(to data_frame)
df_length <- data.frame(
  width = c(
    unlist(cal_res[[1]]), unlist(cal_res[[2]]), unlist(cal_res[[3]]),
    unlist(cal_res[[4]]), unlist(cal_res[[5]])
  ),
  group = c(
    rep("All drivers", length(cal_res[[1]])),
    rep(">= 2", length(cal_res[[2]])),
    rep(">= 5", length(cal_res[[3]])),
    rep(">= 7", length(cal_res[[4]])),
    rep("Non-drivers", length(cal_res[[5]]))
  )
)
df_length$group <- factor(df_length$group,
                          levels = c("Non-drivers", "All drivers",
                                     ">= 2", ">= 5", ">= 7"))

# plot
# statistic test
stat.test.width <- df_length %>%
  wilcox_test(
    width ~ group,
    p.adjust.method = "BH"
  )
# plot of width
# set y
y_min <- 0
y_q1 <- quantile(df_length$width, 0.25)
y_q3 <- quantile(df_length$width, 0.75)
y_max <- max(df_length$width)
y_diff <- (y_max - y_min) / 3
colors <- c(
  "#eb4b3a", "#48bad0", "#1a9781",
  "#355783", "#ef9a80"
)
trans <- scales::trans_new(
  "custom",
  transform = function(x) ifelse(x <= 2, x / 4, 0.5 + (x - 2) / 180),
  inverse = function(x) ifelse(x <= 0.5, x * 4, 2 + (x - 0.5) * 180)
)
p <- ggplot(df_length) +
  geom_boxplot(aes(group, (as.numeric(width) / 1000), color = group),
               outlier.shape = 21
  ) +
  scale_color_manual(values = c(
    "#eb4b3a", "#48bad0", "#1a9781",
    "#355783", "#ef9a80"
  )) +
  # scale_x_discrete(labels = c('Non-drivers', 'All drivers', '>= 2', '>= 5', '>= 7')) +
  xlab("") +
  ylab("Transcript length") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 30, vjust = 0.8, face = "bold",
      color = colors, size = 13, hjust = 0.65
    ),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18)
  )
# scale_y_continuous(trans = trans, breaks = c(0, 2, 92), labels = c(0, 2, 92))
# scale_y_continuous(breaks = c(0, 25, 50, 75, 92),
#                    limits = c(0, 94))
# add stat.test result
x_value <- c(2, 3, 4, 5, 3, 4, 5, 4, 5, 5)
y_value_data <- unlist(lapply(levels(df_length$group), function(x) {
  return((max(df_length[df_length$group == x, 1])) / 1000)
}))[2:5]
names(y_value_data) <- levels(df_length$group)[2:5]
y_value <- rep(y_value_data, 4:1) + 0.5
y_value <- y_value + c(0.015 * 1:4, 0.015 * 1:3, 0.015 * 1:2, 0.015)
color_value <- c(
  colors[1], colors[1], colors[1], colors[1],
  colors[2], colors[2], colors[2],
  colors[3], colors[3],
  colors[4]
)
for (i in 1:nrow(stat.test.width)) {
  if (stat.test.width$p.adj.signif[i] != "ns") {
    y_tmp <- y_value[i]
    p <- p + annotate(
      geom = "text",
      label = stat.test.width$p.adj.signif[i],
      x = x_value[i],
      y = y_tmp,
      color = color_value[i],
      size = 5
    )
  }
}
p <- p + scale_y_break(breaks = c(2, 90), scales = 0.3, space = 0.3,
                       ticklabels = c(90, 94)) +
  scale_y_continuous(
    breaks = c(0, 2, 4, 25, 50, 75),
    limits = c(0, 94)
  )
# save
ggsave(paste(output_path, "/", "fancy_width_breaked.tiff", sep = ""),
       height = 7.4, width = 8, dpi = 600, units = "cm"
)
