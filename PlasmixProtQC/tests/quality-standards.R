library(PlasmixProtQC)

# Exercise the public report API and inspect its actual Word table, including
# failed/missing metrics, strict boundaries, formatting and absent plots.
metric_names <- c(
  "Number of features", "Recall", "Signal-to-Noise Ratio (SNR)",
  "Relative Correlation with Reference Datasets (RC)"
)
passing_values <- c(1001, 0.91, 5.1, 0.81)
thresholds <- c(1000, 0.90, 5, 0.80)
report_dir <- tempfile("qc-standards-")
dir.create(report_dir)

check_report <- function(values, expected, failed_column = NULL,
                         invalid = FALSE, omit = NULL, duplicate = NULL) {
  metrics <- data.frame(Quality_Metrics = metric_names, Value = values)
  if (!is.null(omit)) metrics <- metrics[-omit, , drop = FALSE]
  if (!is.null(duplicate)) metrics <- rbind(metrics, metrics[duplicate, ])
  result <- list(qc_metrics_table = metrics, input_meta = data.frame(),
                 snr_plot = NULL, cor_plot = NULL)
  path <- file.path(report_dir, "report.docx")
  suppressMessages(generate_protein_report(
    result, report_dir = report_dir, report_name = "report.docx",
    batch_name = "Threshold test"
  ))
  extracted <- file.path(report_dir, "unpacked")
  unzip(path, files = "word/document.xml", exdir = extracted)
  document <- xml2::read_xml(file.path(extracted, "word/document.xml"))
  ns <- xml2::xml_ns(document)
  rows <- xml2::xml_find_all(document, ".//w:tbl[1]/w:tr", ns)
  standards <- xml2::xml_find_all(rows[[2]], "./w:tc", ns)
  cells <- xml2::xml_find_all(rows[[3]], "./w:tc", ns)
  text <- function(node) paste(xml2::xml_text(xml2::xml_find_all(node, ".//w:t", ns)), collapse = "")
  stopifnot(identical(vapply(standards[2:5], text, character(1)),
                      c(">1,000", ">0.90", ">5", ">0.80")))
  stopifnot(identical(text(cells[[6]]), expected))
  if (!is.null(failed_column)) {
    cell <- cells[[failed_column + 1L]]
    if (invalid) stopifnot(identical(text(cell), "-"))
    else stopifnot(grepl("↓", text(cell), fixed = TRUE))
    red <- xml2::xml_find_all(cell, ".//w:color[@w:val='FF0000']", ns)
    stopifnot(length(red) > 0L)
  }
}

check_report(passing_values, "Yes")

# Report tables keep numeric PCC only; n remains available in plot details.
report_input <- list(qc_metrics_table = data.frame(Quality_Metrics = metric_names, Value = passing_values),
  input_meta = data.frame(), snr_plot = NULL, cor_plot = NULL,
  detail_results = list(cor_results = list(N = 22L, Pair = "M/F")))
suppressMessages(generate_protein_report(report_input, report_dir = report_dir, report_name = "template.docx"))
unzip(file.path(report_dir, "template.docx"), files = "word/document.xml", exdir = report_dir)
report_xml <- xml2::read_xml(file.path(report_dir, "word/document.xml"))
report_text <- paste(xml2::xml_text(xml2::xml_find_all(report_xml, ".//w:tbl//w:t", xml2::xml_ns(report_xml))), collapse = "")
stopifnot(!grepl("n=", report_text, fixed = TRUE),
  !grepl("PlasmixBody|PlasmixHeading", as.character(report_xml)))
check_report(thresholds, "No", failed_column = 1L)
for (i in seq_along(thresholds)) {
  values <- passing_values
  values[i] <- thresholds[i]
  check_report(values, "No", failed_column = i)
  values[i] <- thresholds[i] - c(1, 0.01, 0.01, 0.01)[i]
  check_report(values, "No", failed_column = i)
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    values[i] <- bad
    check_report(values, "No", failed_column = i, invalid = TRUE)
  }
  check_report(passing_values, "No", failed_column = i,
               invalid = TRUE, omit = i)
  check_report(passing_values, "No", failed_column = i,
               invalid = TRUE, duplicate = i)
}

# Values just above a threshold must pass even if the display rounds down to it.
check_report(c(1001, 0.90001, 5.0001, 0.8001), "Yes")
check_report(c(1001, 0.91, 4.9999, 0.81), "No", failed_column = 3L)
check_report(c(1001, 0.91, 5.1, 0.7999), "No", failed_column = 4L)

# Verify that the metric calculators themselves preserve precision before
# qc_allmetrics() passes values to the report.
expr <- data.frame(Feature = c("A", "B", "C"),
                   P1 = c(1, 2, 4), P2 = c(1.1, 2.05, 3.9),
                   M1 = c(2, 3, 5), M2 = c(2.15, 2.9, 5.2))
meta <- data.frame(library = names(expr)[-1], sample = c("P", "P", "M", "M"),
                   platform = "Olink")
actual_snr <- qc_snr(expr, meta, plot = FALSE)$SNR
pca <- prcomp(t(as.matrix(expr[-1])), scale. = FALSE)
weights <- summary(pca)$importance[2, 1:2]
distances <- as.matrix(dist(pca$x[, 1:2] %*% diag(sqrt(weights))))^2
intra <- c(distances[1, 2], distances[3, 4])
inter <- distances[1:2, 3:4]
expected_snr <- 10 * log10(mean(inter) / mean(intra))
stopifnot(isTRUE(all.equal(actual_snr, expected_snr, tolerance = 1e-12)))
stopifnot(actual_snr != round(actual_snr, 2))

reference <- new.env()
utils::data("reference_dataset", package = "PlasmixProtQC", envir = reference)
ref <- subset(reference$reference_dataset, platform == "Olink" & sample_pair == "M/F" & UseForPCC)
ref <- ref[is.finite(ref$assigned_value) & ref$assigned_value > 0, ]
ref <- ref[!duplicated(ref$feature), ]
ref <- head(ref, 30)
stopifnot(nrow(ref) >= 10)
test_ratio <- ref$assigned_value + 0.02 * seq_len(nrow(ref))
cor_expr <- data.frame(Feature = ref$feature, F1 = 2, M1 = log2(4 * test_ratio), P1 = 15)
cor_meta <- data.frame(library = c("F1", "M1", "P1"), sample = c("F", "M", "P"), platform = "Olink")
cor_result <- qc_cor(cor_expr, cor_meta, plot = FALSE)
expected_cor <- with(cor_result$logfc, cor(logFC.Test, assigned_value, use = "complete.obs"))
stopifnot(isTRUE(all.equal(cor_result$COR, expected_cor, tolerance = 1e-12)))
stopifnot(cor_result$COR != round(cor_result$COR, 3))
stopifnot(cor_result$Pair == "M/F", cor_result$N == 30, cor_result$ReferenceVersion == "FDU_v2.0_20261101",
          isTRUE(all.equal(cor_result$COR, cor(test_ratio, ref$assigned_value), tolerance = 1e-12)))
missing_f <- qc_cor(cor_expr[, c("Feature", "M1", "P1")], cor_meta[cor_meta$sample != "F", ], plot = FALSE)
stopifnot(is.na(missing_f$COR), missing_f$N == 0, missing_f$Status == "sample_group_missing")
small <- qc_cor(cor_expr[1:9, ], cor_meta, plot = FALSE)
stopifnot(is.na(small$COR), small$N == 9, small$Status == "insufficient_features")
invalid <- cor_expr; invalid$M1[1] <- NA_real_; invalid$F1[2] <- Inf
filtered <- qc_cor(invalid, cor_meta, plot = FALSE)
stopifnot(filtered$N == 28, all(is.finite(filtered$logfc$logFC.Test)))
unused <- subset(reference$reference_dataset, platform == "Olink" & sample_pair == "M/F" & !UseForPCC)
extra <- data.frame(Feature = head(unused$feature, 5), F1 = 2, M1 = 20, P1 = 0)
unchanged <- qc_cor(rbind(cor_expr, extra), cor_meta, plot = FALSE)
stopifnot(unchanged$N == 30, identical(unchanged$COR, cor_result$COR))
# Negative NPX is valid reported coverage; it must not lower Olink Recall.
nominal <- new.env(); utils::data("reference_dataset_quali", package = "PlasmixProtQC", envir = nominal)
nom <- subset(nominal$reference_dataset_quali, platform == "Olink" & sample == "M")
coverage <- data.frame(Feature = nom$feature, M1 = -1)
coverage_meta <- data.frame(library = "M1", sample = "M", platform = "Olink")
stopifnot(qc_recall(coverage, coverage_meta) == 1)

unlink(report_dir, recursive = TRUE)
cat("2026 quality standards: report boundaries, invalid values, formatting and metric precision passed.\n")
