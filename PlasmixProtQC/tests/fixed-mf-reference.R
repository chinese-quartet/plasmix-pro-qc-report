library(PlasmixProtQC)
e <- new.env()
utils::data("reference_dataset", package = "PlasmixProtQC", envir = e)
utils::data("reference_dataset_quali", package = "PlasmixProtQC", envir = e)
ref <- e$reference_dataset
stopifnot(!anyDuplicated(ref[c("platform", "sample_pair", "feature")]),
          all(is.finite(ref$assigned_value) & ref$assigned_value > 0),
          identical(unique(ref$ReferenceVersion), "FDU_v2.0_20261101"))
expected <- c(DIA = 25L, Olink = 134L, SomaScan = 687L)
for (platform_name in names(expected)) {
  selected <- ref[ref$platform == platform_name & ref$sample_pair == "M/F" & ref$UseForPCC, ]
  stopifnot(nrow(selected) == expected[[platform_name]])
  expr <- data.frame(Feature = selected$feature, F1 = 10,
                     M1 = 10 + log2(selected$assigned_value), P1 = 20)
  meta <- data.frame(library = c("F1", "M1", "P1"), sample = c("F", "M", "P"), platform = platform_name)
  result <- qc_cor(expr, meta)
  stopifnot(result$N == nrow(selected), abs(result$COR - 1) < 1e-12, result$Pair == "M/F")
  # Simultaneously reorder columns and leave metadata order unchanged.
  permuted <- qc_cor(expr[, c("Feature", "P1", "F1", "M1")], meta)
  stopifnot(identical(result$N, permuted$N), abs(result$COR - permuted$COR) < 1e-12)
}
cat("Fixed MF panels, bundled reference and input alignment passed.\n")
