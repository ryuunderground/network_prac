################################################
### 1. 패키지 설치 및 로딩 (최소화 & 조건부) ###
################################################
# R 버전 확인 및 업데이트

# Window
#if (!requireNamespace("installr", quietly = TRUE)) {
#  install.packages("installr")
#}
#library(installr)
#if (!check.for.updates.R()) {
#  install.R()
#}

# Mac
# brew update && brew upgrade r && Rscript -e 'cat("✅ Updated R version: ", getRversion(), "\n")'

# CRAN 패키지
cran_packages <- c("clipr", "stringr","calibrate", "ggplot2", "gridExtra", "ggrepel")

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
bioc_packages <- c("edgeR", "limma", "locfit", "plotly")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

#루트 디렉토리 설정
root <- "~/Documents/github/network_prac"
setwd(root)

####################################
### Build Gene Expression Matrix ###
####################################
# 메타데이터 불러오기
meta_data <- read.csv("./inputs/network_prac_inputs_group.csv", header = TRUE)

# 샘플 이름 리스트
sample_names <- as.character(meta_data[, 1])


# Expression list 초기화
expression_list <- list()

# 반복적으로 샘플 읽기
for (i in seq_along(sample_names)) {
  file_path <- file.path(root, "Samples/", sample_names[i])
  temp <- read.table(file_path, sep = "", header = FALSE)
  
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

###################################
###### Remove all zero gene #######
######### for mRNA data ###########
###################################

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
if ("endo" %in% levels(group)) {
  group <- relevel(group, ref = "endo")
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

##################################
##### TMM Normalization #########
##################################
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
dir.create("./TMM_CSV", showWarnings = FALSE)
write.csv(cpm_tmm, "./TMM_CSV/TMM.csv", row.names = TRUE, quote = FALSE)
write.csv(log_cpm_tmm, "./TMM_CSV/logTMM.csv", row.names = TRUE, quote = FALSE)

###########################################
### Estimate dispersion & fit GLM model ###
###########################################
# 1. Dispersion 추정 (공통 + 트렌드 + 유전자별 모두 포함)
dge <- estimateDisp(dge, design)

# 결과 확인
plotBCV(dge)  # Biological Coefficient of Variation plot

# 2. GLM 적합
fit <- glmFit(dge, design)

# 3. 적합 결과 요약
summary(fit)

#########################################
### MDS plot data & axis label 계산 #####
#########################################
# Treat 정의
Treat <- dge$samples$group
levels(Treat)

# MDS 계산 (top 100000개 유전자 기준)
MDS_data <- plotMDS(dge, labels=sample_names, col = as.numeric(Treat), top = 100000)
title("MDS Plot")

# 3차원 MDS
# 필요한 패키지 로드
library(edgeR)
library(plotly)
# sample_names 편집
library(stringr)
sample_names <- str_extract(sample_names,"^[^-]+")
sample_names <- ifelse(
  str_detect(sample_names, ".txt"),
  str_extract(sample_names, "^[^.txt]+"),  # ".txt"가 나오기 전까지 추출
  sample_names  # ".txt"가 없으면 원래 문자열 유지
)

# 1. 거리 계산 (top = 100000은 대부분의 유전자 사용)
mds <- plotMDS(dge, top=100000, plot=FALSE)

# 2. MDS 좌표 추출 (최초 3개의 차원)
mds_coords <- data.frame(
  Sample = sample_names,
  Group = Treat,
  Dim1 = mds$x,
  Dim2 = mds$y,
  Dim3 = mds$distance.matrix[,3]  # 3번째 차원은 distance matrix에서 추출 필요
)

# 3. 3차원 산점도 그리기
fig <- plot_ly(
  mds_coords,
  x = ~Dim1, y = ~Dim2, z = ~Dim3,
  type = "scatter3d",
  mode = "markers+text",
  text = ~Sample,
  color = ~Group,
  colors = c("#ff8c00", "#6499E9", "#228b22"),  # 필요에 맞게 수정
  marker = list(size = 5)
)

fig <- fig %>% layout(
  title = "3D MDS Plot",
  scene = list(
    xaxis = list(title = "Dim1"),
    yaxis = list(title = "Dim2"),
    zaxis = list(title = "Dim3")
  )
)

fig

# 고유값 합
three_groups_eigen <- sum(MDS_data$eigen.values)

# PC1, PC2 백분율 계산 및 축 레이블 생성
pc1_percent <- round(100 * MDS_data$eigen.values[1] / three_groups_eigen)
pc2_percent <- round(100 * MDS_data$eigen.values[2] / three_groups_eigen)
pc3_percent <- round(100 * MDS_data$eigen.values[3] / three_groups_eigen)

PC1 <- paste0("PC1 (", pc1_percent, "%)")
PC2 <- paste0("PC2 (", pc2_percent, "%)")
PC3 <- paste0("PC3 (", pc3_percent, "%)")

#######################
#### MDS Plotting #####
#######################

library(ggplot2)
library(patchwork)

# 1. 그룹 정보 및 MDS 좌표 정리
Treat <- dge$samples$group

plot_data <- data.frame(
  Group = Treat,
  X = MDS_data$x,
  Y = MDS_data$y,
  Sample = sample_names
)

# 2. 그룹별 색상 지정
manual_colors <- c(
  "endo" = "#ff8c00",
  "ovary" = "#6499E9",
  "ovid" = "#228b22"
)

# 3. 축 범위 및 축 라벨
x_lim <- c(-2.5, 2.5)
y_lim <- c(-2.5, 2.5)

x_lab <- PC1
y_lab <- PC2

# 4. 메인 MDS scatter plot
p1 <- ggplot(plot_data, aes(x = X, y = Y, color = Group)) +
  geom_point(size = 3) +
  # 점마다 sample 이름 표시
  geom_text(
    aes(label = Sample),
    vjust = -0.8,
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = manual_colors) +
  scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  xlab(x_lab) +
  ylab(y_lab) +
  theme_bw() +
  theme(
    legend.position = "left",
    panel.border = element_rect(colour = "gray87"),
    plot.margin = margin(0, 0, 0, 0)
  )

# 5. X축 marginal density
p2 <- ggplot(plot_data, aes(x = X, fill = Group, color = Group)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = manual_colors) +
  scale_color_manual(values = manual_colors) +
  scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,10), expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

# 6. Y축 marginal density
p3 <- ggplot(plot_data, aes(x = Y, fill = Group, color = Group)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = manual_colors) +
  scale_color_manual(values = manual_colors) +
  scale_x_continuous(limits = y_lim, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,40), expand = c(0, 0)) +
  coord_flip() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

# 7. 빈 패널
empty <- ggplot() + theme_void()

# 8. 패널 조합
# widths 두 번째 값이 p3 폭
# heights 첫 번째 값이 p2 높이
g <- (p2 + empty) /
  (p1 + p3) +
  plot_layout(
    widths = c(5, 0.9),
    heights = c(0.9, 5)
  )

# 9. 플롯 확인
plot(g)

# 10. PDF 저장
# pdf_path <- file.path(root, "DEG_Results/three_groups_MDS_no_text.pdf")
pdf_path <- file.path(root, "DEG_Results/three_groups_MDS.pdf")
ggsave(filename = pdf_path, plot = g, width = 16, height = 8)

