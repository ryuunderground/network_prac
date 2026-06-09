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
root <- "~/Documents/network_prac"
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
library(gridExtra)

# 1. Treat 벡터: factor로 그룹 지정 (A/B/C 등 자동 인식)
Treat <- dge$samples$group
plot_data <- data.frame(
  Group = Treat, 
  X = MDS_data$x, 
  Y = MDS_data$y
)

endoColor = "#ff8c00"
ovaryColor = "#6499E9"
ovidColor = "#228b22"

# 2. 축 레이블: PC1, PC2
x_lab <- PC1  # 예: "PC1 (42%)"
y_lab <- PC2

# 3. 테마 함수
theme0 <- function(...) theme(
  legend.position = "none",
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.spacing = unit(0, "null"),
  axis.ticks = element_blank(),
  axis.text.x = element_text(margin = margin()),
  axis.text.y = element_text(margin = margin()),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.length = unit(0, "null"),
  panel.border = element_rect(color = NA),
  ...
)


# 4. MDS Scatter Plot
# plot_data에 sample_names 추가

plot_data$Sample <- sample_names

p1 <- ggplot(plot_data, aes(x = X, y = Y, color = Group)) +
  geom_point(size = 3) +
  # 점마다 일자 여부
  geom_text(aes(label = Sample), vjust = -0.8, size = 3.5, show.legend = FALSE) +
  xlab(x_lab) + ylab(y_lab) +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  scale_color_manual(
    values = c(
      "endo" = endoColor,
      "ovary" = ovaryColor,
      "ovid" = ovidColor
    )
  ) +
  theme_bw() +
  theme(legend.position = "left") +
  theme(panel.border = element_rect(colour = "gray87")) +
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5))

# 5. X축 밀도
p2 <- ggplot(plot_data, aes(x = X, fill = Group, color = Group)) +
  geom_density(alpha = 0.7) +
  scale_x_continuous(breaks = NULL, expand = c(0.02, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(0.02, 0)) +
  theme_bw() + theme0(plot.margin = unit(c(1, -0.5, -0.5, 8), "lines"))

# 6. Y축 밀도
p3 <- ggplot(plot_data, aes(x = Y, fill = Group, color = Group)) +
  geom_density(alpha = 0.7) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, expand = c(0.02, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(0.02, 0)) +
  theme_bw() + theme0(plot.margin = unit(c(-0.55, 3, 1.3, -0.2), "lines"))

# 7. 플롯 화면에 표시
g <- grid.arrange(
  arrangeGrob(p2, ncol = 2, widths = c(3, 1)),
  arrangeGrob(p1, p3, ncol = 2, widths = c(3, 1)),
  heights = c(1, 3))
plot(g)


# 8. PDF 저장
# pdf_path <- file.path(root, "DEG_Results/three_groups_MDS_no_text.pdf")
pdf_path <- file.path(root, "DEG_Results/three_groups_MDS.pdf")
ggsave(filename = pdf_path, plot = g, width = 8, height = 8)
