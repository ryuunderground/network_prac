# CRAN 패키지
cran_packages <- c("clipr", "stringr","calibrate", "ggplot2", "gridExtra", "ggrepel", "dplyr")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Bioconductor 패키지
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_packages <- c("edgeR", "limma", "locfit")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

#루트 디렉토리 설정
root <- "~/Documents/GitHub/network_prac"
setwd(root)

endo_Samples <- c("endo03", "endo06", "endo09", "endo12", "endo15", "endo18")

for (date in endo_Samples){
  # 메타데이터 불러오기
  meta_folder <- "./inputs"
  dir.create(meta_folder, showWarnings = FALSE, recursive=TRUE)
  date_csv <- paste0(date,".csv")
  meta_path <- file.path(meta_folder, date_csv)
  meta_data <- read.csv(meta_path, header = TRUE)
  
  # 샘플 이름 리스트
  sample_names <- as.character(meta_data[, 1])
  
  # Expression list 초기화
  expression_list <- list()
  
  # 반복적으로 샘플 읽기
  for (i in seq_along(sample_names)) {
    file_path <- file.path(root, "Samples_endo/", sample_names[i])
    temp <- read.csv(file_path, sep = "", header = FALSE)
    
    if (i == 1) {
      gene_symbols <- as.character(temp[, 1])  # 첫 샘플에서 gene 이름 추출
    }
    
    expression_list[[i]] <- temp[, 2]  # 발현량만 저장
  }
  
  # 발현량을 데이터프레임으로 결합
  Expression <- as.data.frame(do.call(cbind, expression_list))
  rownames(Expression) <- gene_symbols
  colnames(Expression) <- sample_names
  
  # 확인
  head(Expression)
  dim(Expression)
  head(gene_symbols)
  
  ###########################################
  ### 1. 0값인 유전자 제거 (모든 샘플에서)
  ###########################################
  
  nonzero_gene_idx <- rowSums(Expression) != 0
  Expression <- Expression[nonzero_gene_idx, ]
  gene_symbols <- gene_symbols[nonzero_gene_idx]  # Gene_Symbol 변수명을 맞춰줬다면
  
  cat("남은 유전자 수:", nrow(Expression), "\n")
  
  ###########################################
  ### 2. 그룹 정보 설정
  ###########################################
  
  # 그룹 벡터 생성 (예: A, B, C 또는 Control, Treatment 등)
  group <- as.factor(meta_data[, 2])
  
  # 그룹이 2개일 경우 Control 기준 재정렬 (선택적 처리)
  if ("Control" %in% levels(group)) {
    group <- relevel(group, ref = "Control")
  }
  
  print(table(group))  # 그룹 분포 확인
  
  ###########################################
  ### 3. DGEList 생성 및 필터링
  ###########################################
  
  # DGEList 생성
  dge <- DGEList(counts = Expression, group = group)
  
  # 저발현 유전자 필터링
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # 필터링 후 DGEList 재생성 (선택 사항)
  dge <- DGEList(counts = dge$counts, group = group)
  
  # 확인
  cat("필터링 후 유전자 수:", nrow(dge$counts), "\n")
  
  # 1. 정규화 계수 계산 (한 번만!)
  dge <- calcNormFactors(dge)
  
  # 2. 정규화 계수 확인
  cat("Normalization factors:\n")
  print(dge$samples$norm.factors)
  
  # 3. 디자인 행렬 생성 (Treat는 기존 그룹 변수명과 일치해야 함)
  design <- model.matrix(~ group)
  head(design)
  summary(design)
  
  # 4. TMM-normalized CPM 계산
  # - 로그 변환 안 한 값
  cpm_tmm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  # - 로그 변환한 값
  log_cpm_tmm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
  
  # 5. 저장 (디렉토리 먼저 확인 또는 생성 필요)
  tmm_folder <- paste0("./TMM_CSV/", date)
  dir.create(tmm_folder, showWarnings = FALSE, recursive=TRUE)
  tmm_path <- file.path(tmm_folder, paste0(date,"_TMM_endo.csv"))
  log_tmm_path <- file.path(tmm_folder, paste0(date, "_logTMM_endo.csv"))
  write.csv(cpm_tmm, tmm_path, row.names = TRUE, quote = FALSE)
  write.csv(log_cpm_tmm, log_tmm_path, row.names = TRUE, quote = FALSE)
  
  #################################################
  ###### Estimate dispersion & fit GLM model ######
  #################################################
  
  # 1. Dispersion 추정
  bcv <- 0.4
  fit <- glmFit(dge, design, dispersion = bcv^2)
  
  # LRT 실행: 디자인 행렬의 두 번째 열(coef=2)을 기준으로 테스트
  lrt_result <- glmLRT(fit, coef = 2)
  
  # 유전자 전체 결과 테이블 추출 (PValue 기준 정렬)
  result_table <- topTags(lrt_result, n = nrow(dge), sort.by = "PValue")$table
  head(result_table)
  # 필요한 열만 추출 (logFC, PValue, FDR)
  deg_table <- result_table[, c("logFC", "PValue", "FDR")]
  
  # 확인
  head(deg_table)
  
  # 저장 (선택 사항)
  deg_folder <- paste0("./DEG_Results/", date)
  dir.create(deg_folder, showWarnings = FALSE, recursive=TRUE)
  deg_path <- file.path(deg_folder, paste0(date, "_DEG_table_endo.csv"))
  write.csv(deg_table, deg_path, row.names = TRUE, quote = FALSE)
  
  # 1. Annotation 파일 경로 및 불러오기
  annotation_path <- file.path(root, "annotation/pig.txt")
  gene_annot <- read.csv(annotation_path, header = TRUE, stringsAsFactors = FALSE)
  
  # 2. Gene symbol vector
  gene_symbols <- rownames(result_table)  # result_table은 glmLRT 결과에서 만든 테이블
  
  # 3. Symbol 기준으로 Gene name 매칭
  gene_names <- gene_annot[match(gene_symbols, gene_annot[, 1]), 2]
  
  # 4. 결과 테이블 결합
  result_annotated <- cbind(Gene = gene_symbols, GeneName = gene_names, result_table)
  dim(result_annotated)
  head(result_annotated)
  
  # 5. DEG 필터링: FDR < 0.05 & |logFC| ≥ 1 & 이름 있는 것만
  is_up <- result_annotated$logFC >= 1 & result_annotated$FDR < 0.05 & result_annotated$GeneName != ""
  is_down <- result_annotated$logFC <= -1 & result_annotated$FDR < 0.05 & result_annotated$GeneName != ""
  
  DEG_up   <- result_annotated[is_up, ]
  DEG_down <- result_annotated[is_down, ]
  DEG      <- rbind(DEG_up, DEG_down)

  # 6. DEG 이름 추출 및 클립보드 복사
  DEG_names <- DEG$GeneName
  clipr::write_clip(DEG_names)
  
  # 7. 결과 저장
  result_annotated <- result_annotated[c("GeneName", "Gene", "logFC", "PValue", "FDR")]
  colnames(result_annotated)[colnames(result_annotated) == "GeneName"] <- "Gene_name"
  colnames(result_annotated)[colnames(result_annotated) == "Gene"] <- "genes"
  dir.create("./DEG_Results", showWarnings = FALSE, recursive=TRUE)
  write.csv(result_annotated, file.path(root, "./DEG_Results/", paste0(date, ".csv")),
            row.names = FALSE, quote = FALSE)
  write.csv(DEG_names, file.path(root, "./DEG_Results/", paste0(date, "_endo_gene_list.csv")), row.names = FALSE, quote = FALSE)
  
  # 8. 결과 확인
  cat("⬆️  Up-regulated genes:", nrow(DEG_up), "\n")
  cat("⬇️  Down-regulated genes:", nrow(DEG_down), "\n")
  cat("🧬 Total DEG (with names):", length(DEG_names), "\n")
  
  # 라벨 겹침 방지
  library(ggplot2)
  library(ggrepel)
  
  ControlColor = "#FF7ED4"
  endoColor = "#ff8c00"
  
  Result_df <- result_annotated
  
  # 1. 유의한 유전자만 색 지정
  Result_df$group <- "NS"
  Result_df$group[Result_df$logFC >= 1 & Result_df$FDR < 0.05] <- "UP"
  Result_df$group[Result_df$logFC <= -1 & Result_df$FDR < 0.05] <- "DOWN"
  
  # 2. factor 순서 지정 (NS를 넣지 않음)
  Result_df$group <- factor(Result_df$group, levels = c("UP", "DOWN", "NS"))
  y_cutoff <- -log10(max(Result_df$PValue[Result_df$FDR < 0.05], na.rm = TRUE))
  
  # 3. Volcano plot
  p <- ggplot(Result_df, aes(x = logFC, y = -log10(PValue))) +
    geom_point(
      data = subset(Result_df, group == "NS"),
      aes(x = logFC, y = -log10(PValue)),
      color = "gray80", size = 1.5, alpha = 0.5,
      show.legend = FALSE
    ) +
    geom_point(
      data = subset(Result_df, group == "UP"),
      aes(color = group), size = 2
    ) +
    geom_point(
      data = subset(Result_df, group == "DOWN"),
      aes(color = group), size = 2
    ) +
    scale_color_manual(
      values = c("UP" = ControlColor, "DOWN" = endoColor),
      name = "Group"
    ) +
    theme_bw() +
    theme(legend.position = "right") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = y_cutoff, linetype = "dashed", color = "black") +
    labs(
      title = "Volcano Plot",
      x = "log2 Fold Change",
      y = "-log10(P-value)"
    ) +
    coord_cartesian(xlim = c(-17, 17), ylim = c(0, 30))
  
  plot(p)
  
  # Save
  ggsave(filename = file.path(root, "./DEG_Results/", date, paste0(date, "_endo_volcano.pdf")), plot = p, width = 8, height = 8)
}
