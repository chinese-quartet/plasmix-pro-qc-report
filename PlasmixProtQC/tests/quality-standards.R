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
ref <- subset(reference$reference_dataset, platform == "Olink" & sample_pair == "M/P")
ref <- ref[is.finite(ref$assigned_value) & ref$assigned_value > 0, ]
ref <- ref[!duplicated(ref$feature), ]
ref <- head(ref, 30)
stopifnot(nrow(ref) >= 3)
test_ratio <- ref$assigned_value + 0.02 * seq_len(nrow(ref))
cor_expr <- data.frame(Feature = ref$feature, P1 = 2, M1 = log2(4 * test_ratio))
cor_meta <- data.frame(library = c("P1", "M1"), sample = c("P", "M"), platform = "Olink")
cor_result <- qc_cor(cor_expr, cor_meta, plot = FALSE)
expected_cor <- with(cor_result$logfc, cor(logFC.Test, assigned_value, use = "complete.obs"))
stopifnot(isTRUE(all.equal(cor_result$COR, expected_cor, tolerance = 1e-12)))
stopifnot(cor_result$COR != round(cor_result$COR, 3))

unlink(report_dir, recursive = TRUE)
cat("2026 quality standards: report boundaries, invalid values, formatting and metric precision passed.\n")
