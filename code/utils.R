library(stringr)
library(readxl)
library(edgeR)
library(dplyr)
library(tidyverse)
library(rtracklayer)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggsignif)
library(rstatix)
library(ggprism)
library(ggbreak)
library(RobustRankAggreg)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)


conservation_plot <- function(
    path, output_path, standard = "all non-driver", stat_method = "wilcox.test",
    width = 6, height = 8, dpi = 300) {
  files <- list.files(path = path, pattern = "*.csv", full.names = TRUE)
  name_list <- lapply(files, function(file) {
    name <- gsub("\\.csv$", "", basename(file))
    return(name)
  })
  print("reading files: ")
  print(files)
  # data import
  data_list <- lapply(files, function(file) {
    if (grepl("best", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phastCon"))
      tmp <- tmp[, -1]
      tmp <- tmp %>% filter(!(tmp[, 2] == "No data"))
      return(tmp)
    } else if (grepl("fraction", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phyloP"))
      tmp <- tmp[, -1]
      return(tmp)
    } else if (grepl("mean", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phyloP"))
      tmp <- tmp[, -1]
      tmp <- na.omit(tmp)
      return(tmp)
    }
  })
  names(data_list) <- name_list

  # all box plot phastCon
  k <- !(as.data.frame(data_list["best_more_all"])$ID %in% as.data.frame(data_list["best_more1"])$ID)
  table(k)
  non_driver_all <- as.data.frame(data_list["best_more_all"])[k, ]
  all_phastCon_mixed <- data.frame("ID" = NA, "phastCon" = NA)
  group_list <- c()
  for (name in names(data_list)) {
    if (grepl("best", name)) {
      tmp <- as.data.frame(data_list[name])
      colnames(tmp) <- c("ID", "phastCon")
      all_phastCon_mixed <- rbind(all_phastCon_mixed, tmp)
      if (grepl("1", name)) {
        if (grepl("0", name)) {
          group <- rep("more than 10", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        } else {
          group <- rep("all drivers", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        }
      } else if (grepl("all", name)) {
        group <- rep("all non-driver", dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      } else {
        group <- rep(paste("more than", gsub("\\D", "", name)), dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      }
    }
  }
  all_phastCon_mixed <- all_phastCon_mixed[-1, ]
  all_phastCon_mixed$group <- group_list
  comparisons <- lapply(unique(group_list), function(x) {
    if (x != standard) {
      return(c(x, standard))
    }
  })
  comparisons <- Filter(Negate(is.null), comparisons)
  p_phastCon_mixed <- ggplot(all_phastCon_mixed, aes(x = reorder(group, -as.numeric(phastCon)), y = as.numeric(phastCon), fill = group)) +
    geom_boxplot() +
    stat_compare_means(
      comparisons = comparisons,
      label = "p.format", method = stat_method
    ) +
    labs(x = "Driver group", y = "best 200bp phastCon score") +
    theme_minimal()
  output_file_name <- paste(output_path, "/best_phastCon_score.jpg", sep = "")
  ggsave(filename = output_file_name, width = width, height = height, dpi = dpi, plot = p_phastCon_mixed)

  # all box plot phyloP
  k <- !(as.data.frame(data_list["fraction_more_all"])$ID %in% as.data.frame(data_list["fraction_more1"])$ID)
  table(k)
  non_driver_all <- as.data.frame(data_list["fraction_more_all"])[k, ]
  all_phyloP_mixed <- data.frame("ID" = NA, "phyloP" = NA)
  group_list <- c()
  for (name in names(data_list)) {
    if (grepl("fraction", name)) {
      tmp <- as.data.frame(data_list[name])
      colnames(tmp) <- c("ID", "phyloP")
      all_phyloP_mixed <- rbind(all_phyloP_mixed, tmp)
      if (grepl("1", name)) {
        if (grepl("0", name)) {
          group <- rep("more than 10", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        } else {
          group <- rep("all drivers", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        }
      } else if (grepl("all", name)) {
        group <- rep("all non-driver", dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      } else {
        group <- rep(paste("more than", gsub("\\D", "", name)), dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      }
    }
  }
  all_phyloP_mixed <- all_phyloP_mixed[-1, ]
  all_phyloP_mixed$group <- group_list
  comparisons <- lapply(unique(group_list), function(x) {
    if (x != standard) {
      return(c(x, standard))
    }
  })
  comparisons <- Filter(Negate(is.null), comparisons)
  p_phyloP_mixed <- ggplot(all_phyloP_mixed, aes(x = reorder(group, -as.numeric(phyloP)), y = as.numeric(phyloP), fill = group)) +
    geom_boxplot() +
    stat_compare_means(
      comparisons = comparisons,
      label = "p.format", method = stat_method
    ) +
    labs(x = "Driver group", y = "fraction of phyloP score < 0.01") +
    theme_minimal()
  output_file_name <- paste(output_path, "/fraction_phyloP_scores.jpg", sep = "")
  ggsave(filename = output_file_name, width = width, height = height, dpi = dpi, plot = p_phyloP_mixed)

  # mean
  k <- !(as.data.frame(data_list["meanP_more_all"])$ID %in% as.data.frame(data_list["meanP_more1"])$ID)
  table(k)
  non_driver_all <- as.data.frame(data_list["meanP_more_all"])[k, ]
  normal_phyloP_mixed <- data.frame("ID" = NA, "phyloP" = NA)
  group_list <- c()
  for (name in names(data_list)) {
    if (grepl("meanP", name)) {
      tmp <- as.data.frame(data_list[name])
      colnames(tmp) <- c("ID", "phyloP")
      normal_phyloP_mixed <- rbind(normal_phyloP_mixed, tmp)
      if (grepl("1", name)) {
        if (grepl("0", name)) {
          group <- rep("more than 10", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        } else {
          group <- rep("all drivers", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        }
      } else if (grepl("all", name)) {
        group <- rep("all non-driver", dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      } else {
        group <- rep(paste("more than", gsub("\\D", "", name)), dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      }
    }
  }
  normal_phyloP_mixed <- normal_phyloP_mixed[-1, ]
  normal_phyloP_mixed$group <- group_list
  comparisons <- lapply(unique(group_list), function(x) {
    if (x != standard) {
      return(c(x, standard))
    }
  })
  print(unique(group_list))
  comparisons <- Filter(Negate(is.null), comparisons)
  mean_phyloP_mixed <- ggplot(normal_phyloP_mixed, aes(x = reorder(group, -as.numeric(phyloP)), y = as.numeric(phyloP), fill = group)) +
    geom_boxplot() +
    stat_compare_means(
      comparisons = comparisons,
      label = "p.format", method = stat_method
    ) +
    labs(x = "Driver group", y = "mean phyloP") +
    theme_minimal()
  output_file_name <- paste(output_path, "/mean_phyloP_scores.jpg", sep = "")
  ggsave(filename = output_file_name, width = width, height = height, dpi = dpi, plot = mean_phyloP_mixed)
}

seprate_samples <- function(data) {
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

JS_calc <- function(E, Et) {
  left <- (E + Et) / 2
  H_left <- sapply(left, function(x) {
    -x * log(x)
  })
  H_right1 <- -1 / 2 * sapply(E, function(x) {
    -x * log(x)
  })
  H_right2 <- -1 / 2 * sapply(Et, function(x) {
    -x * log(x)
  })
  H_right2 <- na.omit(H_right2)
  return(sum(H_left) + sum(H_right1) + sum(H_right2))
}

lnc_cancer_TSscores <- function(input_file_path, cancer_list, pan_cancer_lnc, output_path, csv_method = F) {
  input_file_path <- paste(input_file_path, "/", sep = "")
  fileName <- dir(input_file_path)
  fileName_list <- sapply(fileName, function(x) {
    str_split(x, "\\.")[[1]][1]
  })
  filePath <- sapply(fileName, function(x) {
    paste(input_file_path, x, sep = "")
  })
  print("reading exp file start")
  start_time <- Sys.time()
  exp_data <- lapply(filePath, function(x) {
    print(paste("reading: ", x, sep = ""))
    read.table(x, header = T, sep = "\t", check.names = F)
  })
  names(exp_data) <- fileName_list
  end_time <- Sys.time()
  print("reading exp file done")
  run_time <- end_time - start_time
  print(run_time)

  print("reading pan-cancer lnc data")
  if (csv_method) {
    pan_lnc <- read.csv(pan_cancer_lnc)
    pan_lnc <- pan_lnc[-1]
    colnames(pan_lnc) <- cancer_list
  } else {
    pan_lnc <- read_xlsx(pan_cancer_lnc, col_names = F)
    colnames(pan_lnc) <- cancer_list
  }

  # if (type(pan_cancer_lnc) == 'str') {
  #     pan_lnc <- read_xlsx(pan_cancer_lnc, col_names = F)
  #     colnames(pan_lnc) <- cancer_list
  # } else {
  #     if (tmp_method == T) {
  #         pan_lnc <- pan_cancer_lnc
  #         print(colnames(pan_lnc))
  #     } else {
  #         pan_lnc <- pan_cancer_lnc
  #         colnames(pan_lnc) <- cancer_list
  #     }
  # }

  calc_progress <- 0
  for (singleCancer in cancer_list) {
    calc_progress <- calc_progress + 1
    print(paste("calc lncs of ", singleCancer, ", the progress: ", calc_progress / length(cancer_list) * 100, sep = ""))
    lnc_list <- pan_lnc[, singleCancer]
    lnc_list <- na.omit(lnc_list)

    # get id with version
    ID <- exp_data[[1]]$Ensembl_ID
    ID_noV <- sapply(ID, function(x) {
      str_split(x, "\\.")[[1]][1]
    })
    k <- ID_noV %in% unlist(lnc_list)
    ID_V <- exp_data[[1]]$Ensembl_ID[k]

    # calculate E vector
    E_df <- data.frame(a = rep("test", length(lnc_list)))
    for (cancer in exp_data) {
      fixed <- seprate_samples(cancer)
      tumor <- fixed[[2]]
      rownames(tumor) <- cancer[, 1]
      tumor <- 2^tumor - 1
      tumor_target <- tumor[ID_V, ]
      E_df <- cbind(E_df, rowMeans(tumor_target))
    }
    E_df <- E_df[, -1]
    colnames(E_df) <- names(exp_data)

    n <- ncol(E_df)
    Et_df <- diag(1, nrow = n, ncol = n)
    colnames(Et_df) <- colnames(E_df)

    TS_result <- data.frame()
    for (lnc in rownames(E_df)) {
      single_lncTS <- list()
      E <- E_df[lnc, ]
      E <- E / sum(E)
      index <- 0
      for (t in 1:nrow(Et_df)) {
        index <- index + 1
        Et <- Et_df[t, ]
        single_lncTS[index] <- 1 - as.numeric(JS_calc(E, Et))**0.5
      }
      TS_result <- rbind(TS_result, single_lncTS)
    }
    colnames(TS_result) <- colnames(E_df)
    rownames(TS_result) <- rownames(E_df)

    max_list <- list()
    specific_list <- list()
    index <- 0
    for (i in rownames(TS_result)) {
      index <- index + 1
      max_list[[index]] <- max(TS_result[i, ])
      specific_list[[index]] <- colnames(TS_result)[which.max(TS_result[i, ])]
    }
    maxTS <- unlist(max_list)
    specific <- unlist(specific_list)
    TS_result <- cbind(TS_result, maxTS)
    TS_result <- cbind(na.omit(TS_result), specific)

    fileName <- paste(output_path, "/", singleCancer, ".csv", sep = "")
    write.csv(TS_result, fileName)
  }
}

single_cancer_TSplot <- function(TS_result_path, output_path) {
  file_path <- paste(TS_result_path, "/", sep = "")
  fileName <- dir(file_path)
  fileName_list <- sapply(fileName, function(x) {
    str_split(x, "\\.")[[1]][1]
  })
  filePath <- sapply(fileName, function(x) {
    paste(file_path, x, sep = "")
  })
  TS_data <- lapply(filePath, function(x) {
    read.table(x, header = T, sep = ",", check.names = F)
  })
  names(TS_data) <- fileName_list

  for (i in fileName_list) {
    tmp <- TS_data[[i]]
    tissue <- tmp$specific
    tissue_count <- round(table(tissue) / length(tissue) * 100, 1)
    most_frequent_cancer <- names(tissue_count)[which.max(tissue_count)]
    colors <- rep("skyblue", length(tissue_count))
    colors[which(names(tissue_count) == most_frequent_cancer)] <- "red"

    data <- data.frame(cancer_type = names(tissue_count), specific_count = tissue_count)

    ggplot(data, aes(x = cancer_type, y = specific_count.Freq, fill = cancer_type)) +
      geom_bar(stat = "identity") +
      ggtitle(i) +
      xlab("cancer type") +
      ylab("percent of cancer") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(values = colors) +
      geom_text(aes(label = specific_count.Freq), vjust = -0.5) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(), # 去除背景
        legend.position = "none"
      ) # 去除图例
    ggsave(paste(output_path, "/", i, ".pdf", sep = ""))
  }
}

group_threshold_select <- function(max = 14, threshold, groups) {
  element_out <- c()
  for (count in c(threshold:max)) {
    element_out <- c(element_out, unlist(groups[count]))
  }
  return(element_out)
}

pandri_TSscore_file <- function(pan_cancer_lnc, all_lnc, method = c(1, 2, 5, 7, 10), output_path) {
  pan_cancer_lnc_data <- read_xlsx(pan_cancer_lnc, col_names = F)
  all_lnc <- read_xlsx(all_lnc, col_names = F)
  pan_lnc_list <- unlist(pan_cancer_lnc_data)
  # pan_lnc_list <- na.omit(unlist(pan_cancer_lnc))
  counts <- table(pan_lnc_list)
  groups <- split(names(counts), counts)
  tmp_list <- list()
  index <- 0
  for (pandri_count in method) {
    index <- index + 1
    tmp_list[index] <- list(group_threshold_select(threshold = pandri_count, groups = groups))
  }
  max_length <- max(unlist(lapply(tmp_list, length)))

  pandri_df <- data.frame(tmp = rep(NA, max_length))
  for (pandri in tmp_list) {
    length(pandri) <- max_length
    pandri_df <- cbind(pandri_df, pandri)
  }
  pandri_df <- pandri_df[-1]

  tmp <- c()
  for (pandri_count in method) {
    tmp <- c(tmp, paste("pandri", pandri_count, sep = ""))
  }
  colnames(pandri_df) <- tmp

  k <- unlist(all_lnc[, 1]) %in% unlist(pandri_df[, 1])
  all_no_lnc <- all_lnc[!k, 1]
  colnames(all_no_lnc) <- c("non-drivers")

  # pandri file
  write.csv(pandri_df, paste(output_path, "/", "pandri_", paste(method, collapse = ""), ".csv", sep = ""))
  # non_drivers file
  write.csv(all_no_lnc, paste(output_path, "/", "nondri_all", ".csv", sep = ""))
}

TS_panplot <- function(
    TS_result_path, standard = "all non-driver", width = 6, height = 8, dpi = 300,
    output_path) {
  files <- list.files(path = TS_result_path, pattern = "*.csv", full.names = TRUE)
  name_list <- lapply(files, function(file) {
    name <- gsub("\\.csv$", "", basename(file))
    return(name)
  })
  TS_data <- lapply(files, function(x) {
    read.csv(x, header = T, sep = ",", check.names = F)
  })
  names(TS_data) <- name_list

  index <- 0
  df_long <- data.frame()
  for (data in TS_data) {
    index <- index + 1
    tmp <- data$maxTS
    label <- rep(names(TS_data)[index], dim(data)[1])
    tmp <- cbind(tmp, label)
    df_long <- rbind(df_long, tmp)
  }

  group_list <- c()
  for (name in names(TS_data)) {
    if (grepl("1", name)) {
      if (grepl("0", name)) {
        group <- rep("more than 10", dim(as.data.frame(TS_data[name]))[1])
        group_list <- c(group_list, group)
      } else {
        group <- rep("all drivers", dim(as.data.frame(TS_data[name]))[1])
        group_list <- c(group_list, group)
      }
    } else if (grepl("non", name)) {
      group <- rep("all non-driver", dim(as.data.frame(TS_data[name]))[1])
      group_list <- c(group_list, group)
    } else {
      group <- rep(paste("more than", gsub("\\D", "", name)), dim(as.data.frame(TS_data[name]))[1])
      group_list <- c(group_list, group)
    }
  }

  comparisons <- lapply(unique(group_list), function(x) {
    if (x != standard) {
      return(c(x, standard))
    }
  })

  df_long$label <- group_list
  comparisons <- Filter(Negate(is.null), comparisons)


  p <- ggplot(df_long, aes(x = reorder(label, -as.numeric(tmp)), y = as.numeric(tmp), fill = label)) +
    geom_boxplot() +
    stat_compare_means(
      comparisons = comparisons,
      label = "p.format", method = "wilcox.test"
    ) +
    labs(x = "Driver group", y = "TS") +
    theme_minimal()

  output_file_name <- paste(output_path, "/pandri_TS_plot.jpg", sep = "")
  ggsave(filename = output_file_name, width = width, height = height, dpi = dpi, plot = p)
}

phastCon_data_df <- function(input_path) {
  files <- list.files(path = input_path, pattern = "*.csv", full.names = TRUE)
  name_list <- lapply(files, function(file) {
    name <- gsub("\\.csv$", "", basename(file))
    return(name)
  })
  print("reading files: ")
  print(files)
  # data import
  data_list <- lapply(files, function(file) {
    if (grepl("best", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phastCon"))
      tmp <- tmp[, -1]
      tmp <- tmp %>% filter(!(tmp[, 2] == "No data"))
      return(tmp)
    } else if (grepl("fraction", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phyloP"))
      tmp <- tmp[, -1]
      return(tmp)
    } else if (grepl("mean", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phyloP"))
      tmp <- tmp[, -1]
      tmp <- na.omit(tmp)
      return(tmp)
    }
  })
  names(data_list) <- name_list

  # all box plot phastCon
  k <- !(as.data.frame(data_list["best_more_all"])$ID %in% as.data.frame(data_list["best_more1"])$ID)
  table(k)
  non_driver_all <- as.data.frame(data_list["best_more_all"])[k, ]
  all_phastCon_mixed <- data.frame("ID" = NA, "phastCon" = NA)
  group_list <- c()
  for (name in names(data_list)) {
    if (grepl("best", name)) {
      tmp <- as.data.frame(data_list[name])
      colnames(tmp) <- c("ID", "phastCon")
      all_phastCon_mixed <- rbind(all_phastCon_mixed, tmp)
      if (grepl("1", name)) {
        if (grepl("0", name)) {
          group <- rep("more than 10", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        } else {
          group <- rep("all drivers", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        }
      } else if (grepl("all", name)) {
        group <- rep("all non-driver", dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      } else {
        group <- rep(paste("more than", gsub("\\D", "", name)), dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      }
    }
  }
  all_phastCon_mixed <- all_phastCon_mixed[-1, ]
  all_phastCon_mixed$group <- group_list

  return(all_phastCon_mixed[-1])
}

fraction_data_df <- function(input_path) {
  files <- list.files(path = input_path, pattern = "*.csv", full.names = TRUE)
  name_list <- lapply(files, function(file) {
    name <- gsub("\\.csv$", "", basename(file))
    return(name)
  })
  print("reading files: ")
  print(files)
  # data import
  data_list <- lapply(files, function(file) {
    if (grepl("best", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phastCon"))
      tmp <- tmp[, -1]
      tmp <- tmp %>% filter(!(tmp[, 2] == "No data"))
      return(tmp)
    } else if (grepl("fraction", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phyloP"))
      tmp <- tmp[, -1]
      return(tmp)
    } else if (grepl("mean", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phyloP"))
      tmp <- tmp[, -1]
      tmp <- na.omit(tmp)
      return(tmp)
    }
  })
  names(data_list) <- name_list

  # all box plot phastCon
  k <- !(as.data.frame(data_list["fraction_more_all"])$ID %in% as.data.frame(data_list["fraction_more1"])$ID)
  table(k)
  non_driver_all <- as.data.frame(data_list["fraction_more_all"])[k, ]
  all_phastCon_mixed <- data.frame("ID" = NA, "fraction" = NA)
  group_list <- c()
  for (name in names(data_list)) {
    if (grepl("fraction", name)) {
      tmp <- as.data.frame(data_list[name])
      colnames(tmp) <- c("ID", "fraction")
      all_phastCon_mixed <- rbind(all_phastCon_mixed, tmp)
      if (grepl("1", name)) {
        if (grepl("0", name)) {
          group <- rep("more than 10", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        } else {
          group <- rep("all drivers", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        }
      } else if (grepl("all", name)) {
        group <- rep("all non-driver", dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      } else {
        group <- rep(paste("more than", gsub("\\D", "", name)), dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      }
    }
  }
  all_phastCon_mixed <- all_phastCon_mixed[-1, ]
  all_phastCon_mixed$group <- group_list

  return(all_phastCon_mixed[-1])
}

meanP_data_df <- function(input_path) {
  files <- list.files(path = input_path, pattern = "*.csv", full.names = TRUE)
  name_list <- lapply(files, function(file) {
    name <- gsub("\\.csv$", "", basename(file))
    return(name)
  })
  print("reading files: ")
  print(files)
  # data import
  data_list <- lapply(files, function(file) {
    if (grepl("best", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phastCon"))
      tmp <- tmp[, -1]
      tmp <- tmp %>% filter(!(tmp[, 2] == "No data"))
      return(tmp)
    } else if (grepl("fraction", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phyloP"))
      tmp <- tmp[, -1]
      return(tmp)
    } else if (grepl("mean", file)) {
      tmp <- read.csv(file, header = T, col.names = c("x", "ID", "phyloP"))
      tmp <- tmp[, -1]
      tmp <- na.omit(tmp)
      return(tmp)
    }
  })
  names(data_list) <- name_list

  # all box plot phastCon
  k <- !(as.data.frame(data_list["meanP_more_all"])$ID %in% as.data.frame(data_list["meanP_more1"])$ID)
  table(k)
  non_driver_all <- as.data.frame(data_list["meanP_more_all"])[k, ]
  all_phastCon_mixed <- data.frame("ID" = NA, "meanP" = NA)
  group_list <- c()
  for (name in names(data_list)) {
    if (grepl("meanP", name)) {
      tmp <- as.data.frame(data_list[name])
      colnames(tmp) <- c("ID", "meanP")
      all_phastCon_mixed <- rbind(all_phastCon_mixed, tmp)
      if (grepl("1", name)) {
        if (grepl("0", name)) {
          group <- rep("more than 10", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        } else {
          group <- rep("all drivers", dim(as.data.frame(data_list[name]))[1])
          group_list <- c(group_list, group)
        }
      } else if (grepl("all", name)) {
        group <- rep("all non-driver", dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      } else {
        group <- rep(paste("more than", gsub("\\D", "", name)), dim(as.data.frame(data_list[name]))[1])
        group_list <- c(group_list, group)
      }
    }
  }
  all_phastCon_mixed <- all_phastCon_mixed[-1, ]
  all_phastCon_mixed$group <- group_list

  return(all_phastCon_mixed[-1])
}

pandriTS_data_df <- function(TS_result_path) {
  files <- list.files(path = TS_result_path, pattern = "*.csv", full.names = TRUE)
  name_list <- lapply(files, function(file) {
    name <- gsub("\\.csv$", "", basename(file))
    return(name)
  })
  TS_data <- lapply(files, function(x) {
    read.csv(x, header = T, sep = ",", check.names = F)
  })
  names(TS_data) <- name_list

  index <- 0
  df_long <- data.frame()
  for (data in TS_data) {
    index <- index + 1
    tmp <- data$maxTS
    label <- rep(names(TS_data)[index], dim(data)[1])
    tmp <- cbind(tmp, label)
    df_long <- rbind(df_long, tmp)
  }

  group_list <- c()
  for (name in names(TS_data)) {
    if (grepl("1", name)) {
      if (grepl("0", name)) {
        group <- rep("more than 10", dim(as.data.frame(TS_data[name]))[1])
        group_list <- c(group_list, group)
      } else {
        group <- rep("all drivers", dim(as.data.frame(TS_data[name]))[1])
        group_list <- c(group_list, group)
      }
    } else if (grepl("non", name)) {
      group <- rep("all non-driver", dim(as.data.frame(TS_data[name]))[1])
      group_list <- c(group_list, group)
    } else {
      group <- rep(paste("more than", gsub("\\D", "", name)), dim(as.data.frame(TS_data[name]))[1])
      group_list <- c(group_list, group)
    }
  }
  df_long$label <- group_list

  return(df_long)
}

group_creator <- function(df) {
  group <- c()
  combinations <- expand.grid(levels(factor(df$group2)), levels(factor(df$cancer)))
  index <- 0
  standard_mask <- c()
  for (i in 1:nrow(combinations)) {
    index <- index + 1
    group <- append(group, paste("r", index, sep = ""))
    if (grepl("non", combinations[i, 1])) {
      standard_mask <- append(standard_mask, paste("r", index, sep = ""))
    }
  }
  print("the standard mask is:")
  print(standard_mask)
  combinations$group <- group
  group_ <- c()
  for (i in 1:nrow(df)) {
    tmp_1 <- df[i, "group2"]
    tmp_2 <- df[i, "cancer"]
    group_ <- append(group_, subset(combinations, Var1 == tmp_1 & Var2 == tmp_2)$group)
  }
  print(combinations)
  return(group_)
}

multi_group_boxplot <- function(data, output_path, height = 7, width = 14, dpi = 300) {
  data$values <- as.numeric(data$values)
  data_long <- data

  x_value <- c(0.7, 1.3)
  # 12group, 11 times
  for (i in 1:9) {
    x_value <- append(x_value, x_value[((i - 1) * 2 + 1):(i * 2)] + 1)
  }

  med_value <- data %>%
    group_by(group) %>%
    summarise(median_value = median(values))


  tmp_data <- data.frame(
    x_value = x_value,
    med_value = rep(med_value$median_value, each = 2),
    group = rep(med_value$group, each = 2)
  )


  data_long$group <- factor(data_long$group,
    levels = c("r5", "r4", "r3", "r1", "r2", "r10", "r9", "r8", "r6", "r7")
  )
  tmp_data <- arrange(tmp_data, match(group, c("r5", "r4", "r3", "r1", "r2", "r10", "r9", "r8", "r6", "r7")))
  x_value <- c(0.7, 1.3)
  # 12group, 11 times
  for (i in 1:9) {
    x_value <- append(x_value, x_value[((i - 1) * 2 + 1):(i * 2)] + 1)
  }
  tmp_data$x_value <- x_value
  p <- ggplot(data_long) +
    geom_boxplot(aes(group, values, fill = group2, color = group2),
      width = 0.5,
      outlier.shape = NA
    ) +
    geom_line(data = tmp_data, aes(x_value, med_value, group = group), color = "#ffffff") +
    geom_rect(aes(xmin = 0, xmax = 11, ymin = 1.9, ymax = 2.0), fill = "#eaeae0") +
    geom_vline(xintercept = c(5.5), color = "#bcbdbf", alpha = 0.8) +
    scale_x_discrete(labels = rep(c(
      ">= 7", ">= 5",
      ">= 2", "All drivers", "Non-drivers"
    ), 2)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
    scale_fill_manual(values = c("#696d8c", "#bd8874", "#e6e3b7", "#ac9e8a", "#8794a6")) +
    scale_color_manual(values = c("#696d8c", "#bd8874", "#e6e3b7", "#ac9e8a", "#8794a6")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 20)
    ) +
    ggtitle("Sequence conservation & Tissue specific") +
    geom_signif(aes(group, values),
      comparisons = list(
        c("r1", "r2"),
        c("r3", "r2"),
        c("r4", "r2"),
        c("r5", "r2"),
        c("r8", "r7"),
        c("r9", "r7"),
        c("r10", "r7"),
        c("r6", "r7")
      ),
      map_signif_level = T,
      vjust = 0.68,
      tip_length = rep(0, 8),
      textsize = 8,
      y_position = c(1.05, 1.1, 1.15, 1.2, 0.55, 0.6, 0.65, 0.7)
    ) +
    annotate("text", x = c(3, 8), y = 1.95, label = unique(data_long$cancer), size = 7)

  ggsave(paste(output_path, "/", "multi_group_boxplot.jpg", sep = ""),
    height = height, width = width, dpi = dpi
  )
}

# calculate CNRCI fraction
calc_RCI <- function(driver_list, all_lnc, all_driver_list, RCI_data) {
  all_lnc_list <- all_lnc$...1
  k <- all_lnc_list %in% all_driver_list
  all_non <- all_lnc_list[!k]
  k_non <- RCI_data$ENSEMBL.ID %in% all_non
  k_driver <- RCI_data$ENSEMBL.ID %in% driver_list
  df_driver <- RCI_data[k_driver, ]
  df_non <- RCI_data[k_non, ]
  driver_fraction <- list()
  index <- 0
  for (single_lnc in driver_list) {
    index <- index + 1
    single_lnc_df <- df_driver[df_driver$ENSEMBL.ID == single_lnc, ]
    value_list <- single_lnc_df$Value
    neg_count <- sum(unlist(value_list) < 0)
    driver_fraction[index] <- neg_count / length(value_list)
    driver_fraction <- as.numeric(driver_fraction)
  }
  driver_fraction <- driver_fraction[!is.na(driver_fraction)]
  non_fraction <- list()
  index <- 0
  for (single_non in all_non) {
    index <- index + 1
    single_non_df <- df_non[df_non$ENSEMBL.ID == single_non, ]
    value_list <- single_non_df$Value
    neg_count <- sum(unlist(value_list) < 0)
    non_fraction[index] <- neg_count / length(value_list)
    non_fraction <- as.numeric(non_fraction)
  }
  non_fraction <- non_fraction[!is.na(non_fraction)]
  wilcox.test(driver_fraction, non_fraction)
}
