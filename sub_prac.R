#######################
#### 파일 불러오기 ####
#######################

file_list <- list.files(
  path = "Raw_data",
  pattern = "\\.csv$",
  full.names = FALSE
)

file_list <- sub("\\.csv$", "", file_list)

data_list <- lapply(file_list, function(fname) {
  read.csv(
    file.path("Raw_data", paste0(fname, ".csv")),
    header = TRUE
  )
})

names(data_list) <- file_list

all_genes <- unique(unlist(lapply(data_list, function(df) df$genes)))

######################################################
#### FC cutoff별 PCIT input 자동 생성 (FDR < 0.01) ####
######################################################

fdr_cutoff <- 0.01
fdr_label <- "FDR001"

fc_cutoffs <- c(
  FC2  = 1,
  FC4  = 2,
  FC8  = 3,
  FC16 = 4,
  FC32 = 5
)

dir.create("PCIT_Input", showWarnings = FALSE, recursive = TRUE)
dir.create("DEG_Filter_Summary", showWarnings = FALSE, recursive = TRUE)

fc_overall_summary <- data.frame()

for (fc_name in names(fc_cutoffs)) {
  
  logfc_cutoff <- fc_cutoffs[[fc_name]]
  
  cat("\n============================\n")
  cat("Processing:", fc_name, "\n")
  cat("logFC cutoff:", logfc_cutoff, "\n")
  cat("FDR cutoff:", fdr_cutoff, "\n")
  cat("============================\n")
  
  deg_genes <- unique(unlist(lapply(data_list, function(df) {
    df$genes[
      abs(df$logFC) >= logfc_cutoff &
        df$FDR < fdr_cutoff
    ]
  })))
  
  merged_deg <- data.frame(genes = deg_genes)
  
  removed_genes <- setdiff(all_genes, deg_genes)
  removed_percent <- round(length(removed_genes) / length(all_genes) * 100, 2)
  
  cat(fc_name, "DEG union gene count:", length(deg_genes), "\n")
  cat("전체 unique genes:", length(all_genes), "\n")
  cat("제외된 unique genes:", length(removed_genes), "\n")
  cat("제외율:", removed_percent, "%\n")
  
  deg_summary <- data.frame(
    Sample = file_list,
    Total_genes = sapply(data_list, nrow),
    Passed = sapply(data_list, function(df) {
      sum(
        abs(df$logFC) >= logfc_cutoff &
          df$FDR < fdr_cutoff,
        na.rm = TRUE
      )
    })
  )
  
  deg_summary$Removed_genes <- deg_summary$Total_genes - deg_summary$Passed
  deg_summary$Removed_percent <- round(
    deg_summary$Removed_genes / deg_summary$Total_genes * 100,
    2
  )
  deg_summary$Passed_percent <- round(
    deg_summary$Passed / deg_summary$Total_genes * 100,
    2
  )
  
  write.csv(
    deg_summary,
    file.path(
      "DEG_Filter_Summary",
      paste0("DEG_filter_summary_", fc_name, "_", fdr_label, ".csv")
    ),
    row.names = FALSE,
    quote = FALSE
  )
  
  pcit_input <- merged_deg
  
  for (i in seq_along(file_list)) {
    df <- data_list[[file_list[i]]][, c("genes", "logFC")]
    colnames(df)[2] <- file_list[i]
    
    pcit_input <- merge(
      pcit_input,
      df,
      by = "genes",
      all.x = TRUE
    )
  }
  
  pcit_input[is.na(pcit_input)] <- 0
  
  cat("PCIT input gene count:", nrow(pcit_input), "\n")
  cat("PCIT input condition count:", ncol(pcit_input) - 1, "\n")
  
  write.table(
    pcit_input,
    file.path(
      "PCIT_Input",
      paste0("PCIT_input_", fc_name, "_", fdr_label, ".txt")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  fc_overall_summary <- rbind(
    fc_overall_summary,
    data.frame(
      Cutoff = fc_name,
      logFC_cutoff = logfc_cutoff,
      FDR_cutoff = fdr_cutoff,
      Total_unique_genes = length(all_genes),
      Passed_unique_genes = length(deg_genes),
      Removed_unique_genes = length(removed_genes),
      Removed_percent = removed_percent
    )
  )
}

write.csv(
  fc_overall_summary,
  paste0("PCIT_Input_FC_cutoff_overall_summary_", fdr_label, ".csv"),
  row.names = FALSE,
  quote = FALSE
)

print(fc_overall_summary)

###############################
##### Attribute file 생성 #####
###############################

attribute <- data_list[[file_list[1]]][, c("genes", "logFC", "FDR")]

colnames(attribute) <- c(
  "genes",
  paste0(file_list[1], "_FC"),
  paste0(file_list[1], "_FDR")
)

for (i in 2:length(file_list)) {
  
  temp_df <- data_list[[file_list[i]]][, c("genes", "logFC", "FDR")]
  
  colnames(temp_df) <- c(
    "genes",
    paste0(file_list[i], "_FC"),
    paste0(file_list[i], "_FDR")
  )
  
  attribute <- merge(
    attribute,
    temp_df,
    by = "genes",
    all = TRUE
  )
}

attribute[is.na(attribute)] <- 0

fc_cols <- grep("_FC$", colnames(attribute))
fdr_cols <- grep("_FDR$", colnames(attribute))

attribute$absMAX_FC <- apply(
  abs(attribute[, fc_cols, drop = FALSE]),
  1,
  max,
  na.rm = TRUE
)

max_fc_idx <- apply(
  abs(attribute[, fc_cols, drop = FALSE]),
  1,
  which.max
)

attribute$FC_sample <- sub(
  "_FC$",
  "",
  colnames(attribute)[fc_cols][max_fc_idx]
)

attribute$Min_FDR <- apply(
  attribute[, fdr_cols, drop = FALSE],
  1,
  min,
  na.rm = TRUE
)

min_fdr_idx <- apply(
  attribute[, fdr_cols, drop = FALSE],
  1,
  which.min
)

attribute$FDR_sample <- sub(
  "_FDR$",
  "",
  colnames(attribute)[fdr_cols][min_fdr_idx]
)

write.table(
  attribute,
  paste0("Attribute_FC_FDR_summary_", fdr_label, ".txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nDone.\n")
cat("Generated PCIT input files in: PCIT_Input/\n")
cat("Generated DEG summary files in: DEG_Filter_Summary/\n")
cat("Generated attribute file:", paste0("Attribute_FC_FDR_summary_", fdr_label, ".txt"), "\n")

