#' Statistics for basic information
#' @param expr_dt A expression profile
#' @param meta_dt A metadata file
#' @import stats
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @importFrom ggplot2 geom_abline
#' @export
qc_info <- function(expr_dt, meta_dt) {
  # Load data
  m <- meta_dt
  d <- expr_dt
  s <- meta_dt$sample

  # Replace zero by NA
  d[d == 0] <- NA

  # 1. Feature Count
  uniq_pro <- unique(d[, 1])
  stat_num <- length(uniq_pro)

  # 2. Missing %
  d_vals <- d[, 2:ncol(d)]
  d_all_num <- nrow(d_vals) * ncol(d_vals)
  d_missing_num <- length(which(is.na(d_vals)))
  prop_missing <- 0
  if (d_all_num > 0) prop_missing <- d_missing_num * 100 / d_all_num
  stat_missing <- round(prop_missing, 3)

  # 3. Correlation & CV
  samples <- table(s)
  rep_samples <- samples[samples > 1]

  stat_acor <- NA
  stat_cv <- NA

  if (length(rep_samples) > 0) {
    # Correlation (Log data correlation is fine)
    d_mtx <- d[, 2:ncol(d)]
    common_libs <- intersect(colnames(d_mtx), m$library)
    if (length(common_libs) >= 2) {
      d_mtx <- d_mtx[, common_libs]
      cor_res <- cor(d_mtx, use = "pairwise.complete.obs", method = "pearson")
      stat_acor <- round(median(cor_res[upper.tri(cor_res)], na.rm = TRUE), 3)
    }

    # CV Calculation
    d_long <- melt(d, id.vars = colnames(d)[1])
    d_long <- na.omit(d_long)
    d_long <- merge(d_long, m, by.x = "variable", by.y = "library")

    # 强制假设输入数据是 Log2 的，先转回 Linear 再算 CV
    d_long$value_linear <- 2^(d_long$value)

    d_cv <- aggregate(
      value_linear ~ Feature + sample,
      data = d_long,
      FUN = function(x) sd(x) / mean(x)
    )
    stat_cv <- round(median(d_cv$value_linear, na.rm = T) * 100, 3)
  }

  return(c(stat_num, stat_missing, stat_acor, stat_cv))
}

#' Get color mapping
#' @export
get_sample_colors <- function(samples) {
  cols <- c("M" = "#3171b8", "Y" = "#53a949", "P" = "#6d3390", "X" = "#ffc65d", "F" = "#ae231c")
  res <- cols[samples]
  res[is.na(res)] <- "gray50"
  return(res)
}

#' Calculating SNR
#' @import stats
#' @import utils
#' @importFrom data.table data.table setkey
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme labs guides guide_legend ggsave scale_x_continuous scale_y_continuous
#' @importFrom ggthemes theme_few
#' @export
qc_snr <- function(expr_dt, meta_dt, output_dir = NULL, plot = TRUE) {
  # 数据准备
  expr_df <- data.frame(expr_dt[, 2:ncol(expr_dt)], row.names = expr_dt[, 1])
  expr_df[is.na(expr_df)] <- 0

  # 对齐 Metadata
  libs <- colnames(expr_df)
  meta_sub <- meta_dt[match(libs, meta_dt$library), ]
  group <- meta_sub$sample

  # --- 修改 1: PCA 缩放逻辑调整 ---
  # 检查平台信息
  current_platform <- unique(meta_sub$platform)

  # 默认不 Scale (SomaScan/Olink 都是相对定量单位)
  do_scale <- FALSE

  # 只有 DIA 平台才开启 Scale
  if (length(current_platform) > 0 && any(grepl("DIA", current_platform, ignore.case = TRUE))) {
    do_scale <- TRUE
  }

  # PCA
  expr_t <- t(expr_df)
  vars <- apply(expr_t, 2, var)
  expr_t <- expr_t[, vars > 0 & !is.na(vars), drop = FALSE]

  if (ncol(expr_t) < 2) stop("Too few features for PCA.")

  # 应用 Scale 参数
  pca <- prcomp(expr_t, scale. = do_scale)

  pcs <- as.data.frame(pca$x)
  pcs$sample <- group

  # 计算 SNR
  prop <- summary(pca)$importance[2, 1:2]

  n <- nrow(pcs)
  dist_val <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      d2 <- prop[1] * (pcs[i, 1] - pcs[j, 1])^2 + prop[2] * (pcs[i, 2] - pcs[j, 2])^2
      dist_val[i, j] <- d2
    }
  }

  intra_dists <- c()
  inter_dists <- c()

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (group[i] == group[j]) {
        intra_dists <- c(intra_dists, dist_val[i, j])
      } else {
        inter_dists <- c(inter_dists, dist_val[i, j])
      }
    }
  }

  snr <- mean(inter_dists) / mean(intra_dists)
  # 保留完整精度用于质量判定；只在图表展示时四舍五入。
  snr_db <- 10 * log10(snr)

  # Plot
  p <- NULL
  if (plot) {
    cols <- get_sample_colors(unique(group))
    p <- ggplot(pcs, aes(x = PC1, y = PC2, color = sample)) +
      geom_point(size = 6) +
      scale_color_manual(values = cols) +
      theme_few() +
      labs(
        title = paste0("SNR = ", round(snr_db, 2), " dB"),
        x = sprintf("PC1 (%.1f%%)", prop[1] * 100),
        y = sprintf("PC2 (%.1f%%)", prop[2] * 100)
      )
  }

  if (!is.null(output_dir)) {
    if (plot) ggsave(file.path(output_dir, "pca_plot.png"), p, width = 6, height = 5)
    write.table(pcs, file.path(output_dir, "pca_table.tsv"), sep = "\t", quote = F)
  }

  return(list(SNR = snr_db, snr_plot = p))
}

#' Calculating ratio PCC against a fixed reference panel
#' @param expr_dt Log2 expression table with feature identifiers in the first column.
#' @param meta_dt Sample metadata containing library, sample and platform.
#' @param output_dir Optional output directory, retained for API compatibility.
#' @param plot Whether to draw the ratio scatter plot.
#' @param show_sample_pairs Retained for API compatibility.
#' @param sample_pair Ratio to evaluate; M/F is used by the standard report.
#' @importFrom utils data
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
#' @export
qc_cor <- function(expr_dt, meta_dt, output_dir = NULL, plot = FALSE,
                   show_sample_pairs = TRUE, sample_pair = "M/F") {
  expr_dt <- as.data.frame(expr_dt); meta_dt <- as.data.frame(meta_dt)
  plat <- unique(meta_dt$platform)
  if (length(plat) != 1L) stop("One platform is required per report.")
  parts <- strsplit(sample_pair, "/", fixed = TRUE)[[1]]
  if (length(parts) != 2L || parts[1] == parts[2]) stop("Invalid sample_pair.")
  utils::data("reference_dataset", package = "PlasmixProtQC", envir = environment())
  ref <- reference_dataset[reference_dataset$platform == plat & reference_dataset$sample_pair == sample_pair, , drop = FALSE]
  if (!all(c("UseForPCC", "ReferenceVersion") %in% names(ref))) stop("The bundled reference needs updating.")
  version <- unique(as.character(ref$ReferenceVersion))
  ref <- ref[ref$UseForPCC %in% TRUE, , drop = FALSE]
  if (anyDuplicated(ref$feature)) stop("Duplicate fixed reference features.")
  md <- meta_dt[match(names(expr_dt)[-1], meta_dt$library), , drop = FALSE]
  if (anyNA(md$library) || anyDuplicated(md$library) || anyDuplicated(expr_dt[[1]])) stop("Input feature/library alignment is invalid.")
  mat <- as.matrix(expr_dt[, -1, drop = FALSE])
  rownames(mat) <- expr_dt[[1]]
  groups <- lapply(parts, function(x) which(md$sample == x))
  common <- intersect(rownames(mat), ref$feature)
  merged <- data.frame(Feature = character(), Sample.Pair = character(), logFC.Test = numeric(), assigned_value = numeric())
  status <- "sample_group_missing"
  if (all(lengths(groups) > 0L)) {
    linear <- 2^mat[common, , drop = FALSE]
    ratio <- rowMeans(linear[, groups[[1]], drop = FALSE], na.rm = TRUE) /
      rowMeans(linear[, groups[[2]], drop = FALSE], na.rm = TRUE)
    assigned <- ref$assigned_value[match(common, ref$feature)]
    keep <- is.finite(ratio) & ratio > 0 & is.finite(assigned) & assigned > 0
    merged <- data.frame(Feature = common[keep], Sample.Pair = rep(sample_pair, sum(keep)),
                         logFC.Test = ratio[keep], assigned_value = assigned[keep])
    status <- if (nrow(merged) < 10L) "insufficient_features" else "computed"
  }
  n <- nrow(merged); value <- NA_real_
  if (status == "computed") {
    if (sd(merged$logFC.Test) > 0 && sd(merged$assigned_value) > 0) value <- cor(merged$logFC.Test, merged$assigned_value)
    else status <- "zero_variance"
  }
  p <- NULL
  if (plot && n > 0L) p <- ggplot(merged, aes(x = assigned_value, y = logFC.Test)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_point(color = "#3171b8", alpha = .7, size = 2) + ggthemes::theme_few() +
    labs(title = sprintf("%s PCC = %.3f (n=%d)", sample_pair, value, n),
         x = paste(sample_pair, "reference ratio"), y = paste(sample_pair, "test ratio"))
  list(COR = value, N = n, Pair = sample_pair, ReferenceVersion = paste(version, collapse = ";"),
       Status = status, cor_plot = p, logfc = merged)
}

#' Recall
#' @export
qc_recall <- function(expr_dt, meta_dt) {
  plat <- unique(meta_dt$platform)
  utils::data("reference_dataset_quali", package = "PlasmixProtQC", envir = environment())
  ref <- reference_dataset_quali[reference_dataset_quali$platform == plat, ]
  if (nrow(ref) == 0) {
    return(NA)
  }

  samps <- unique(meta_dt$sample)
  common <- intersect(samps, gsub("Quartet ", "", unique(ref$sample)))

  if (length(common) == 0) {
    return(NA)
  }

  vals <- c()
  for (s in common) {
    feats_ref <- ref$feature[gsub("Quartet ", "", ref$sample) == s]
    libs <- meta_dt$library[meta_dt$sample == s]

    for (l in libs) {
      values <- expr_dt[[l]]
      detected <- if (plat == "Olink") is.finite(values) else !is.na(values) & values > 0
      det <- expr_dt[[1]][detected]
      recall <- length(intersect(det, feats_ref)) / length(feats_ref)
      vals <- c(vals, recall)
    }
  }
  return(mean(vals))
}
