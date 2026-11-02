# 用固定MF参考生成一个小型演示，不修改CSV或包内RDA。
library(PlasmixProtQC)
data("reference_dataset", package = "PlasmixProtQC")
data("reference_dataset_quali", package = "PlasmixProtQC")
ref <- subset(reference_dataset, platform == "DIA" & sample_pair == "M/F" & UseForPCC)
features <- unique(subset(reference_dataset_quali, platform == "DIA")$feature)
ratio <- ref$assigned_value[match(features, ref$feature)]
ratio[is.na(ratio)] <- 1
meta <- expand.grid(sample = c("M", "F", "P"), replicate = 1:3, stringsAsFactors = FALSE)
meta$library <- paste(meta$sample, meta$replicate, sep = "_")
meta$platform <- "DIA"
set.seed(123)
expr <- sapply(meta$sample, function(s) {
  factor <- switch(s, M = ratio, F = rep(1, length(ratio)), P = (ratio + 1) / 2)
  20 + log2(factor) + rnorm(length(factor), sd = .02)
})
colnames(expr) <- meta$library
expr <- data.frame(Feature = features, expr, check.names = FALSE)
out <- tempfile("plasmix-demo-"); dir.create(out)
write.csv(expr, file.path(out, "expression.csv"), row.names = FALSE)
write.csv(meta[, c("library", "sample", "platform")], file.path(out, "metadata.csv"), row.names = FALSE)
qc <- qc_conclusion(file.path(out, "expression.csv"), file.path(out, "metadata.csv"), output_dir = out, plot = FALSE)
print(qc$qc_metrics_table)
message("演示输出：", out, "。蛋白数不足1000，因此即使MF PCC高也不能判整体合格。")
