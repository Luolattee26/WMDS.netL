library(dplyr)
library(tidyverse)
if (!require("rtracklayer")) {
  BiocManager::install("rtracklayer")
} else {
  library("rtracklayer")
}
library(readxl)
library(survival)
library(survminer)
library(Rtsne)
library(paletteer)
library(stringr)
library(stringi)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(edgeR)
library(NMF)
library(FactoMineR)
library(factoextra)
library(sva)
if (!require("ConsensusClusterPlus")) {
  BiocManager::install("ConsensusClusterPlus")
} else {
  library("ConsensusClusterPlus")
}
library(ggplot2)
library(stats)
library(patchwork)
library(pROC)
library(boot)
library(caret)
library(latex2exp)
library(rstatix)
library(ggbreak)


seprate_samples <- function(data, data_preprocess = TRUE) {
  if (is.character(data)) {
    samples_all <- read.table(data, header = T, sep = "\t", check.names = F)
  } else {
    samples_all <- data
  }
  samples_all <- samples_all[, -1]
  tumor_list <- vector("list")
  normal_index <- vector("list")
  n <- 0
  list_index <- 0
  if (data_preprocess) {
    print("the input data is log2(FPKM + 1)")
    print("convert it into log2(TPM + 1)")
    samples_all <- log2(as.data.frame(apply((2^samples_all - 1), 2, FPKM2TPM)) + 1)
  } else {
    print("the input data is log2(FPKM + 1)")
    print("stay unchanged in log2(FPKM + 1)")
    samples_all <- samples_all
  }
  for (i in colnames(samples_all)) {
    n <- n + 1
    i_split <- str_split(i, "-")
    tissue_barcode <- substr(i_split[[1]][4], 1, 2)
    if (as.numeric(tissue_barcode) == 11) {
      list_index <- list_index + 1
      normal_index[list_index] <- n
    }
  }
  normal_index <- as.vector(unlist(normal_index))
  normal_data <- samples_all[, normal_index]
  n <- 0
  list_index <- 0
  for (i in colnames(samples_all)) {
    n <- n + 1
    i_split <- str_split(i, "-")
    tissue_barcode <- substr(i_split[[1]][4], 1, 2)
    if (as.numeric(tissue_barcode) == 01) {
      list_index <- list_index + 1
      tumor_list[list_index] <- n
    }
  }
  tumor_list <- as.vector(unlist(tumor_list))
  tumor_data <- samples_all[, tumor_list]
  data <- list(normal_data, tumor_data)
  return(data)
}

fileName2cancerList <- function(fileName_list) {
  tmp <- lapply(fileName_list, function(x) {
    tmp1 <- str_split_1(x, "\\.")[1]
    tmp2 <- str_split_1(tmp1, "-")[2]
    return(tmp2)
  })
  return(unlist(tmp))
}

# get symbol
get_symbol_list <- function(ID_list, all_lnc_file) {
  list_index <- 0
  lst <- vector("list")
  for (i in ID_list) {
    symbol_index <- 0
    for (j in all_lnc_file$...1) {
      symbol_index <- symbol_index + 1
      i_split <- str_split(i, "\\.")
      if (j == i_split[[1]][1]) {
        list_index <- list_index + 1
        lst[list_index] <- all_lnc_file[symbol_index, 2]
      }
    }
  }
  return(lst)
}

# add_exp_feature
add_exp_feature <- function(tumor_exp_file, target_list, percent) {
  exp_group <- vector("list", length = length(target_list))
  n <- 0
  for (i in 1:length(target_list)) {
    n <- n + 1
    lnc <- target_list[[i]]
    exp_list <- unlist(tumor_exp_file[lnc, ])
    lst <- t(as.data.frame(tumor_exp_file[lnc, ]))
    thresholds_1 <- quantile(exp_list, probs = percent)
    thresholds_2 <- quantile(exp_list, probs = 1 - percent)
    high <- max(thresholds_1, thresholds_2)
    min <- min(thresholds_1, thresholds_2)
    # print(high)
    group <- vector("list")
    for (i in 1:length(lst)) {
      if (high == min) {
        if (exp_list[i] == 0) {
          group[i] <- NA
        } else if (exp_list[i] > high) {
          group[i] <- 1
        } else if (exp_list[i] <= high) {
          group[i] <- 0
        } else {
          group[i] <- NA
        }
      } else {
        if (exp_list[i] == 0) {
          group[i] <- NA
        } else if (exp_list[i] >= high) {
          group[i] <- 1
        } else if (exp_list[i] <= min) {
          group[i] <- 0
        } else {
          group[i] <- NA
        }
      }
    }
    lst <- rbind(lst, group)
    exp_group[[n]] <- lst
  }
  names(exp_group) <- lncRNA_list
  return(exp_group)
}

# mix survival and group
mix <- function(survival_file, group_data) {
  replacement <- 0
  survival_mixed <- survival_file
  col_index <- ncol(survival_mixed)
  for (gene in names(group_data)) {
    survival_mixed[gene] <- survival_mixed$OS.time
    data <- (group_data[[gene]])
    col_index <- col_index + 1
    for (j in 1:nrow(survival_mixed)) {
      rowname <- survival_mixed[j, 1]
      n <- 0
      for (name in colnames(data)) {
        n <- n + 1
        if (name == rowname) {
          replacement <- data[2, n]
        }
      }
      if (is.null(replacement)) {
        survival_mixed[j, col_index] <- NA
      } else {
        survival_mixed[j, col_index] <- replacement
      }
    }
  }
  return(survival_mixed)
}

# mix survival and expression
mix_exp <- function(survival_file, group_data) {
  replacement_exp <- 0
  survival_mixed <- survival_file
  col_index <- ncol(survival_mixed)
  for (gene in names(group_data)) {
    survival_mixed[gene] <- survival_mixed$OS.time
    data <- (group_data[[gene]])
    col_index <- col_index + 1
    for (j in 1:nrow(survival_mixed)) {
      rowname <- survival_mixed[j, 1]
      n <- 0
      for (name in colnames(data)) {
        n <- n + 1
        if (name == rowname) {
          replacement_exp <- data[1, n]
          # print(name)
          # print(rowname)
          # print(replacement_exp)
        }
      }
      if (is.null(replacement_exp)) {
        survival_mixed[j, col_index] <- NA
      } else {
        survival_mixed[j, col_index] <- replacement_exp
      }
    }
  }
  return(survival_mixed)
}

mix_group <- function(group, survival_data) {
  for (i in names(group)) {

  }
}

# get ID from symbol
get_ID_list <- function(symbol_list, all_lnc_file) {
  list_index <- 0
  lst <- vector("list", length = length(symbol_list))
  for (i in symbol_list) {
    symbol_index <- 0
    for (j in all_lnc_file$...2) {
      symbol_index <- symbol_index + 1
      if (j == i) {
        list_index <- list_index + 1
        lst[list_index] <- all_lnc_file[symbol_index, 1]
      }
    }
  }
  return(lst)
}

# get TCGA exp data for ML through the pan-cancer cox result
get_samples <- function(exp_data, cancer_type, threshold, cox_data, path, method = "best", padj = FALSE,
                        padj.method = "BH", data_preprocess = TRUE, useTop = FALSE) {
  all_lnc <- read_excel("./data/WMDS_latest/LNC_name.xlsx", col_names = F)
  rownames(exp_data) <- unlist(exp_data[1])
  exp_data <- exp_data[-1]
  if (data_preprocess) {
    print("the input data is log2(FPKM + 1)")
    print("convert it into log2(TPM + 1)")
    exp_data <- log2(as.data.frame(apply((2^exp_data - 1), 2, FPKM2TPM)) + 1)
  } else {
    print("the input data is log2(FPKM + 1)")
    print("stay unchanged in log2(FPKM + 1)")
    exp_data <- exp_data
  }
  if (method == "best") {
    k <- lapply(colnames(exp_data), function(x) {
      if (grepl("GTE", x)) {
        return(TRUE)
      } else {
        barcode <- substr(str_split(x, "-")[[1]][4], 1, 3)
        if (barcode == "01A" | barcode == "11A") {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }
    })
    k <- unlist(k)
    tmp <- exp_data[k]
    status <- lapply(colnames(tmp), function(x) {
      if (grepl("GTE", x)) {
        return(1)
      } else {
        barcode <- substr(str_split(x, "-")[[1]][4], 1, 2)
        if (barcode == "01") {
          return(0)
        } else {
          return(1)
        }
      }
    })
    # tmp <- log2(cpm(2 ^ tmp - 1) + 1)
    # return(tmp)
    if (padj) {
      cox_data[cancer_type] <- p.adjust(unlist(cox_data[cancer_type]), method = padj.method)
    } else {
      print("no p-value adjust")
    }
    if (!useTop) {
      print("use all features")
      sig_cox <- subset(cox_data, cox_data[cancer_type] < threshold)
      symbol_list <- rownames(sig_cox)
      ID_list <- get_ID_list(symbol_list, all_lnc)
    } else {
      print(paste("use top ", useTop, " features", sep = ""))
      # sig_cox <- subset(cox_data, cox_data[cancer_type] < threshold)
      sig_cox <- cox_data[cox_data[cancer_type] < threshold, ]
      symbol_list <- rownames(sig_cox[order(sig_cox[, cancer_type]), ][1:useTop, ])
      ID_list <- get_ID_list(symbol_list, all_lnc)
    }

    ID_no_Version <- lapply(rownames(tmp), function(x) {
      str_split(x, "\\.")[[1]][1]
    })
    k <- ID_no_Version %in% ID_list
    tmp <- tmp[k, ]
    tmp <- rbind(tmp, status)
    rownames(tmp)[nrow(tmp)] <- "status"
    write.csv(tmp, file = paste(path, "/", cancer_type, "_MLdata.csv", sep = ""), quote = F)
    return(tmp)
  }
}

cox_sig_count <- function(cox_data, threshold, padj = FALSE, padj.method = "BH") {
  if (padj) {
    sig_count <- lapply(colnames(cox_data), function(x) {
      tmp <- p.adjust(unlist(cox_data[x]), method = padj.method)
      return(sum(tmp < threshold))
    })
  } else {
    sig_count <- lapply(colnames(cox_data), function(x) {
      tmp <- cox_data[x]
      return(sum(tmp < threshold))
    })
  }
  names(sig_count) <- colnames(cox_data)
  return(sig_count)
}

runPCA <- function(merged_data, path = NULL) {
  label <- lapply(rownames(merged_data), function(x) {
    if (grepl("GTE", x)) {
      return("GTEx")
    } else {
      return("TCGA")
    }
  })
  colname <- colnames(merged_data)[-length(colnames(merged_data))]
  merged_data$batch <- as.factor(unlist(label))
  data <- merged_data[, 1:(dim(merged_data)[2] - 2)]
  data <- data.frame(data)
  data <- t(apply(data, 1, as.numeric))
  colnames(data) <- colname
  pca <- PCA(data, graph = T)

  # PCA plot
  fviz_pca_ind(pca,
    geom.ind = "point",
    addEllipses = TRUE,
    legend.title = "Batch",
    addlabels = TRUE, habillage = merged_data$batch
  )
}

runComBat <- function(merged_data, cancer_type, path = NULL) {
  label <- lapply(rownames(merged_data), function(x) {
    if (grepl("GTE", x)) {
      return("GTEx")
    } else {
      return("TCGA")
    }
  })
  colname <- colnames(merged_data)[-length(colnames(merged_data))]
  merged_data$batch <- as.factor(unlist(label))
  data <- merged_data[, 1:(dim(merged_data)[2] - 2)]
  data <- data.frame(data)
  data <- t(apply(data, 1, as.numeric))
  colnames(data) <- colname
  batch <- unlist(merged_data$batch)
  status <- unlist(merged_data$status)
  design.mat <- model.matrix(~ factor(status))
  data <- ComBat(dat = t(data), batch = batch, mod = design.mat)
  data <- as.data.frame(data)
  data <- rbind(data, status)
  rownames(data)[nrow(data)] <- "status"
  if (!is.null(path)) {
    write.csv(data, file = paste(path, "/", cancer_type, "_MLdata.csv", sep = ""), quote = F)
  }
  return(data)
}

FPKM2TPM <- function(fpkm) {
  fpkm / sum(fpkm) * 1e6
}

get_tissue_sample <- function(gtex_phenotype, tissues) {
  print("there are all the tissue type in this data:")
  print((table(unlist(gtex_phenotype[3]))))
  tissues <- tissues
  sample_list <- vector("list", length = length(tissues))
  index <- 0
  for (tissue in tissues) {
    index <- index + 1
    samples_fit <- list(gtex_phenotype[gtex_phenotype[3] == tissue, 1])
    sample_list[index] <- samples_fit
  }
  names(sample_list) <- tissues
  return(sample_list)
}

get_exp_data_single <- function(driver_list, tumor_data, normal_data, method = "mean", all_non_file, with_non_drivers = TRUE) {
  if (with_non_drivers) {
    exp_list <- vector("list", 4)
    names(exp_list) <- c("tumor", "normal", "tumor_non", "normal_non")
  } else {
    exp_list <- vector("list", 2)
    names(exp_list) <- c("tumor", "normal")
  }
  load(all_non_file)
  ID_noV <- unlist(lapply(rownames(tumor_data), function(x) {
    return(str_split_1(x, "\\.")[1])
  }))
  k <- ID_noV %in% driver_list
  tmp_tmuor <- tumor_data[k, ]
  tmp_normal <- normal_data[k, ]
  k <- ID_noV %in% all_non
  tmp_tmuor_non <- tumor_data[k, ]
  tmp_normal_non <- normal_data[k, ]

  if (with_non_drivers) {
    if (method == "mean") {
      exp_tumor <- apply(tmp_tmuor, 1, mean)
      exp_normal <- apply(tmp_normal, 1, mean)
      exp_tumor_non <- apply(tmp_tmuor_non, 1, mean)
      exp_normal_non <- apply(tmp_normal_non, 1, mean)
    } else if (method == "median") {
      exp_tumor <- apply(tmp_tmuor, 1, median)
      exp_normal <- apply(tmp_normal, 1, median)
      exp_tumor_non <- apply(tmp_tmuor_non, 1, median)
      exp_normal_non <- apply(tmp_normal_non, 1, median)
    } else {
      print("the method argument only supports mean or median")
    }
    exp_list[1] <- list(exp_tumor)
    exp_list[2] <- list(exp_normal)
    exp_list[3] <- list(exp_tumor_non)
    exp_list[4] <- list(exp_normal_non)
  } else {
    if (method == "mean") {
      exp_tumor <- apply(tmp_tmuor, 1, mean)
      exp_normal <- apply(tmp_normal, 1, mean)
    } else if (method == "median") {
      exp_tumor <- apply(tmp_tmuor, 1, median)
      exp_normal <- apply(tmp_normal, 1, median)
    } else {
      print("the method argument only supports mean or median")
    }
    exp_list[1] <- list(exp_tumor)
    exp_list[2] <- list(exp_normal)
  }
  return(exp_list)
}

exp_plot <- function(exp_list, name) {
  group1 <- exp_list[[1]]
  group2 <- exp_list[[2]]
  group3 <- exp_list[[3]]
  group4 <- exp_list[[4]]
  cdf1 <- ecdf(group1)
  cdf2 <- ecdf(group2)
  cdf3 <- ecdf(group3)
  cdf4 <- ecdf(group4)
  x <- seq(min(group1, group2, group3, group4), max(group1, group2, group3, group4), length.out = 10000)
  y1 <- cdf1(x)
  y2 <- cdf2(x)
  y3 <- cdf3(x)
  y4 <- cdf4(x)

  data <- data.frame(x = c(x, x, x, x), y = c(y1, y2, y3, y4), group = c(
    rep("drivers(tumor)", length(x)),
    rep("drivers(normal)", length(x)),
    rep("non-drivers(tumor)", length(x)),
    rep("non-drivers(normal)", length(x))
  ))

  ks_test12 <- ks.test(group1, group2)
  ks_test13 <- ks.test(group1, group3)
  ks_test34 <- ks.test(group3, group4)

  p <- ggplot(data, aes(x = x, y = y, color = group)) +
    geom_line() +
    geom_text(aes(x = max(x), y = 0.2, label = paste("drivers(tumor) vs drivers(normal):", format.pval(ks_test12$p.value, digits = 3))),
      color = "black", hjust = 1, size = 3.5
    ) +
    geom_text(aes(x = max(x), y = 0.1, label = paste("drivers(tumor) vs non-drivers(tumor):", format.pval(ks_test13$p.value, digits = 3))), color = "black", hjust = 1, size = 3.5) +
    geom_text(aes(x = max(x), y = 0, label = paste("non-drivers(tumor) vs non-drivers(normal):", format.pval(ks_test34$p.value, digits = 3))), color = "black", hjust = 1, size = 3.5) +
    labs(title = name, x = "log2(TPM + 1)", y = "Cumulative Probability") +
    theme_classic() +
    theme(
      title = element_text(size = 15),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 14),
      axis.text.x = element_text(size = 13),
      axis.title.x = element_text(size = 14),
      axis.text.y = element_text(size = 13),
      axis.title.y = element_text(size = 14),
    )
  return(p)
}

pancancer_exp_plot <- function(exp_data_path, ouput_path, driver_type = 7) {
  gc()
  files <- list.files(exp_data_path)
  plot_list <- vector("list", length = length(files))
  name_list <- vector("list", length = length(files))
  index <- 0
  for (file in files) {
    gc()
    index <- index + 1
    cancer <- sub("^TCGA-(\\w+).*", "\\1", file)
    name_list[index] <- cancer
    print(paste("running on dataset of: ", cancer, sep = ""))
    # input TCGA exp data
    path_tcga <- paste("./data/TCGA_14cancer/fpkm/TCGA-", cancer, ".htseq_fpkm.tsv", sep = "")
    raw_data <- read.table(path_tcga, header = T, sep = "\t", check.names = F)
    fixed_samples <- seprate_samples(raw_data)
    normal_df <- fixed_samples[[1]]
    tumor_df <- fixed_samples[[2]]
    rownames(normal_df) <- raw_data[, 1]
    rownames(tumor_df) <- raw_data[, 1]
    # FPKM to TPM
    normal_df <- as.matrix(normal_df)
    tumor_df <- as.matrix(tumor_df)

    # input driver lnc data
    driver <- paste("./data/pancancer_driver/pan_dri_", driver_type, ".csv", sep = "")
    pan_7_driver <- read.csv(driver)
    lncRNA_list <- unlist(pan_7_driver$lnc)

    # get exp list
    exp_list <- get_exp_data_single(
      driver_list = lncRNA_list,
      tumor_data = tumor_df,
      normal_data = normal_df,
      method = "mean",
      all_non_file = "./data/pancancer_driver/all_non_drivers.Rdata",
      with_non_drivers = TRUE
    )
    p <- exp_plot(exp_list, cancer)
    output_file_name <- paste(ouput_path, "/", cancer, ".tiff", sep = "")
    ggsave(filename = output_file_name, width = 6, height = 6, dpi = 600, plot = p)
  }
}

# violin data prepare
violin_data <- function(tcga_path, driver_list, pattern = 1) {
  files <- list.files(path = tcga_path, full.names = TRUE, pattern = "*.tsv")
  name_list <- lapply(files, function(file) {
    name <- sub(".*/TCGA-(\\w+).*", "\\1", file)
    return(name)
  })
  print("reading files list: ")
  print(files)

  data_list <- lapply(files, function(file) {
    print(paste("is reading: ", file, sep = ""))
    raw_data <- read.table(file, header = T, sep = "\t", check.names = F)
    fixed_samples <- seprate_samples(raw_data)
    normal_df <- fixed_samples[[1]]
    tumor_df <- fixed_samples[[2]]
    rownames(normal_df) <- raw_data[, 1]
    rownames(tumor_df) <- raw_data[, 1]
    # FPKM to TPM
    normal_df <- as.matrix(normal_df)
    tumor_df <- as.matrix(tumor_df)
    tmp_res <- list(tumor_df, normal_df)
    names(tmp_res) <- c("tumor", "normal")
    return(tmp_res)
  })
  names(data_list) <- name_list

  print("the pattern supported:")
  print("1--drivers in tumor vs. non-drivers in tumor")
  print("2--drivers in normal vs. non-drivers in normal")
  print("3--drivers in tumor vs. drivers in normal")
  print("4--non-drivers in tumor vs. non-drivers in normal")

  if (pattern == 1) {
    print("now is pattern 1")
    index <- 0
    result <- data.frame()
    for (cancer in data_list) {
      gc()
      index <- index + 1
      # cancer[[1]] is the tumor data
      tmp_data <- cancer[[1]]
      ID_noV <- unlist(lapply(rownames(tmp_data), function(x) {
        return(str_split_1(x, "\\.")[1])
      }))
      k <- ID_noV %in% driver_list
      tmp_driver <- tmp_data[k, ]
      tmp_non <- tmp_data[!k, ]
      exp_driver <- unlist(apply(tmp_driver, 1, mean))
      exp_non <- unlist(apply(tmp_non, 1, mean))
      tmp_res <- data.frame(
        value = c(exp_driver, exp_non),
        group = c(rep("drivers", length(exp_driver)), rep("non-drivers", length(exp_non))),
        cancer = rep(name_list[[index]], length(c(exp_driver, exp_non)))
      )
      # return(tmp_res)
      result <- rbind(result, tmp_res)
    }
    return(result)
  } else if (pattern == 2) {
    print("now is pattern 2")
    index <- 0
    result <- data.frame()
    for (cancer in data_list) {
      gc()
      index <- index + 1
      # cancer[[2]] is the normal data
      tmp_data <- cancer[[2]]
      ID_noV <- unlist(lapply(rownames(tmp_data), function(x) {
        return(str_split_1(x, "\\.")[1])
      }))
      k <- ID_noV %in% driver_list
      tmp_driver <- tmp_data[k, ]
      tmp_non <- tmp_data[!k, ]
      exp_driver <- unlist(apply(tmp_driver, 1, mean))
      exp_non <- unlist(apply(tmp_non, 1, mean))
      tmp_res <- data.frame(
        value = c(exp_driver, exp_non),
        group = c(rep("drivers", length(exp_driver)), rep("non-drivers", length(exp_non))),
        cancer = rep(name_list[[index]], length(c(exp_driver, exp_non)))
      )
      # return(tmp_res)
      result <- rbind(result, tmp_res)
    }
    return(result)
  } else if (pattern == 3) {
    print("now is pattern 3")
    index <- 0
    result <- data.frame()
    for (cancer in data_list) {
      gc()
      index <- index + 1
      # cancer[[1]] is the tumor data
      # cancer[[2]] is the normal data
      tmp_data_1 <- cancer[[1]]
      tmp_data_2 <- cancer[[2]]
      ID_noV <- unlist(lapply(rownames(tmp_data_1), function(x) {
        return(str_split_1(x, "\\.")[1])
      }))
      k <- ID_noV %in% driver_list
      tmp_driver <- tmp_data_1[k, ]
      tmp_non <- tmp_data_2[k, ]
      exp_driver <- unlist(apply(tmp_driver, 1, mean))
      exp_non <- unlist(apply(tmp_non, 1, mean))
      tmp_res <- data.frame(
        value = c(exp_driver, exp_non),
        group = c(rep("tumor", length(exp_driver)), rep("normal", length(exp_non))),
        cancer = rep(name_list[[index]], length(c(exp_driver, exp_non)))
      )
      # return(tmp_res)
      result <- rbind(result, tmp_res)
    }
    return(result)
  } else if (pattern == 4) {
    print("now is pattern 4")
    index <- 0
    result <- data.frame()
    for (cancer in data_list) {
      gc()
      index <- index + 1
      # cancer[[1]] is the tumor data
      # cancer[[2]] is the normal data
      tmp_data_1 <- cancer[[1]]
      tmp_data_2 <- cancer[[2]]
      ID_noV <- unlist(lapply(rownames(tmp_data_1), function(x) {
        return(str_split_1(x, "\\.")[1])
      }))
      k <- ID_noV %in% driver_list
      tmp_driver <- tmp_data_1[!k, ]
      tmp_non <- tmp_data_2[!k, ]
      exp_driver <- unlist(apply(tmp_driver, 1, mean))
      exp_non <- unlist(apply(tmp_non, 1, mean))
      tmp_res <- data.frame(
        value = c(exp_driver, exp_non),
        group = c(rep("tumor", length(exp_driver)), rep("normal", length(exp_non))),
        cancer = rep(name_list[[index]], length(c(exp_driver, exp_non)))
      )
      # return(tmp_res)
      result <- rbind(result, tmp_res)
    }
    return(result)
  }
}

cox_sig_count_forLNC <- function(cox_data, threshold, padj = FALSE, padj.method = "BH") {
  if (padj) {
    sig_count <- lapply(rownames(cox_data), function(x) {
      tmp <- p.adjust(unlist(cox_data[x, ]), method = padj.method)
      return(sum(tmp < threshold))
    })
  } else {
    sig_count <- lapply(rownames(cox_data), function(x) {
      tmp <- cox_data[x, ]
      return(sum(tmp < threshold))
    })
  }
  names(sig_count) <- rownames(cox_data)
  return(sig_count)
}

get_heatmap_data_list <- function(input_path) {
  files <- list.files(input_path)
  file_list <- lapply(files, function(x) {
    tmp_data <- read.csv(paste(input_path, "/", x, sep = ""),
      sep = "\t"
    )
    return(tmp_data)
  })
  file_name <- lapply(files, function(x) {
    str_split_1(x, "\\.")[[1]]
  })
  names(file_list) <- file_name
  return(file_list)
}

Get_SUBCELLULAR_LOCALIZATION_exp <- function(sub_location_data, target_gtf_df) {
  result <- vector("list")
  index <- 0
  for (data in sub_location_data) {
    index <- index + 1
    if (dim(data)[2] == 19) {
      data$id <- lapply(data$gene_id, function(x) {
        return(str_split_1(x, "\\.")[[1]])
      })
      k <- data$id %in% target_gtf_df$ID
      tmp_data <- data[k, ]
      result[[index]] <- mean(tmp_data[, 6])
    } else {
      split_res <- lapply(data$target_id, function(x) {
        return(str_split_1(x, "\\|")[2])
      })
      split_res <- na.omit(split_res)
      data$id <- lapply(split_res, function(x) {
        if (!is.na(x)) {
          return(str_split_1(x, "\\.")[[1]])
        }
      })
      k <- data$id %in% target_gtf_df$ID
      tmp_data <- data[k, ]
      result[[index]] <- mean(tmp_data[, 5])
    }
  }
  names(result) <- names(sub_location_data)
  return(result)
}
