# 从两张发布CSV重建两份主要内置RDA；对象名保持不变。
script <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script)) dirname(normalizePath(sub("^--file=", "", script[1]))) else getwd()
if (!file.exists(file.path(root, "标称特性.csv"))) root <- dirname(root)
reference_dataset_quali <- read.csv(file.path(root, "标称特性.csv"), check.names = FALSE)
reference_dataset <- read.csv(file.path(root, "特性量值.csv"), check.names = FALSE)
quali_keys <- c("platform", "sample", "feature")
quant_keys <- c("platform", "sample_pair", "feature")
stopifnot(all(c(quali_keys, "ReferenceVersion") %in% names(reference_dataset_quali)),
  all(c(quant_keys, "assigned_value", "UseForPCC", "ReferenceVersion") %in% names(reference_dataset)),
  !anyDuplicated(reference_dataset_quali[quali_keys]), !anyDuplicated(reference_dataset[quant_keys]),
  all(is.finite(reference_dataset$assigned_value) & reference_dataset$assigned_value > 0),
  is.logical(reference_dataset$UseForPCC), !anyNA(reference_dataset$UseForPCC))
out <- file.path(root, "PlasmixProtQC/data")
save(reference_dataset_quali, file = file.path(out, "reference_dataset_quali.rda"), version = 2)
save(reference_dataset, file = file.path(out, "reference_dataset.rda"), version = 2)
message("已更新两份主要内置RDA；M/F报告名单由UseForPCC固定。")
