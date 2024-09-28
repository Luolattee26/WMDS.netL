# This code is used to visualize the survival curve of cancer patients
# Based on the NMF subtypes of cancer patients



rm(list = ls())
source("./code/utils_TCGAmodel.R")
out_path <- "./output/NMF/UCEC/"


# read data
load("./output/NMF/UCEC/result_FPKM_3group.Rdata")
# plot
sfit <- survfit(Surv(time / 30, event) ~ group,
                data = survival_LUAD
)
sur_p <- ggsurvplot(sfit,
                    pval = T, palette = "jco", conf.int = TRUE,
                    risk.table = TRUE
)
# use ggplot2 to customize
# plot body
sur_p$plot <- sur_p$plot +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15)
  ) +
  labs(title = "UCEC", x = "Time", y = "Survival probability")
# risk table
sur_p$table <- sur_p$table +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 17)
  )
# save
# plot
ggsave(
  plot = sur_p$plot, "./output/NMF/UCEC/survplot_FPKM_3group.tiff",
  width = 15, height = 12, units = "cm", dpi = 600
)
ggsave(
  plot = sur_p$table, "./output/NMF/UCEC/survtable_FPKM_3group.tiff",
  width = 15, height = 3, units = "cm", dpi = 600
)
