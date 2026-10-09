#' @title Generate Protein Quality Report (Chinese Version)
#' @description Generates a report following the ShenKang Plasmix template structure.
#' @param qc_result list, the output from qc_conclusion()
#' @param report_template character, path to a .docx file to use as style template (optional)
#' @param report_dir character, directory to save the report
#' @param report_name character, filename of the report
#' @param batch_name character, optional, manual batch name. If NULL, tries to infer from metadata.
#' @importFrom dplyr %>%
#' @importFrom flextable flextable theme_box align width bold bg color body_add_flextable set_header_labels
#' @importFrom officer read_docx body_add_par body_add_gg body_add_break body_add_fpar fpar fp_text prop_section
#' @export
generate_protein_report <- function(qc_result,
                                    report_template = NULL,
                                    report_dir = NULL,
                                    report_name = NULL,
                                    batch_name = NULL) {
  # --- 1. 参数检查与初始化 ---
  if (is.null(qc_result)) stop("qc_result is required.")

  if (is.null(report_dir)) {
    report_dir <- file.path(getwd(), "output")
    dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)
  }

  if (is.null(report_name)) {
    report_name <- paste0("Plasmix_QC_Report_", format(Sys.Date(), "%Y%m%d"), ".docx")
  }
  output_file <- file.path(report_dir, report_name)

  # --- 2. 提取并处理数据 ---
  raw_table <- qc_result$qc_metrics_table
  meta <- qc_result$input_meta

  # 尝试获取批次名 (优先用参数，其次查 metadata 的 study_id，最后用默认值)
  if (is.null(batch_name)) {
    if ("study_id" %in% colnames(meta) && length(unique(meta$study_id)) == 1) {
      batch_name <- unique(meta$study_id)[1]
    } else {
      batch_name <- paste0(format(Sys.Date(), "%Y%m%d"), "_data")
    }
  }

  # 提取单个指标；缺失、重复或无法转换的值不能用于判定通过。
  metric_value <- function(metric_name) {
    values <- raw_table$Value[raw_table$Quality_Metrics == metric_name]
    if (length(values) != 1L) return(NA_real_)
    suppressWarnings(as.numeric(as.character(values)))
  }
  val_prot <- metric_value("Number of features")
  val_recall <- metric_value("Recall")
  val_snr <- metric_value("Signal-to-Noise Ratio (SNR)")
  val_rc <- metric_value("Relative Correlation with Reference Datasets (RC)")

  # --- 3. 2026 年质量标准：四项均须严格大于阈值 ---

  # 3.1 Protein (Target > 1000)
  if (!is.finite(val_prot)) {
    txt_prot <- "-"
    pass_prot <- FALSE
  } else {
    txt_prot <- format(val_prot, big.mark = ",")
    if (val_prot <= 1000) {
      txt_prot <- paste0(txt_prot, " ↓")
      pass_prot <- FALSE
    } else {
      pass_prot <- TRUE
    }
  }

  # 3.2 Recall (Target > 0.90)
  if (!is.finite(val_recall)) {
    txt_recall <- "-"
    pass_recall <- FALSE
  } else {
    txt_recall <- sprintf("%.3f", val_recall)
    if (val_recall <= 0.90) {
      txt_recall <- paste0(txt_recall, " ↓")
      pass_recall <- FALSE
    } else {
      pass_recall <- TRUE
    }
  }

  # 3.3 SNR (Target > 5 dB；使用未四舍五入的值判定)
  if (!is.finite(val_snr)) {
    txt_snr <- "-"
    pass_snr <- FALSE
  } else {
    txt_snr <- sprintf("%.2f", val_snr)
    if (val_snr <= 5) {
      txt_snr <- paste0(txt_snr, " ↓")
      pass_snr <- FALSE
    } else {
      pass_snr <- TRUE
    }
  }

  # 3.4 PCC (内部指标名为 RC，Target > 0.80)
  if (!is.finite(val_rc)) {
    txt_rc <- "-"
    pass_rc <- FALSE
  } else {
    txt_rc <- sprintf("%.3f", val_rc)
    if (val_rc <= 0.80) {
      txt_rc <- paste0(txt_rc, " ↓")
      pass_rc <- FALSE
    } else {
      pass_rc <- TRUE
    }
  }

  # 3.5 整体判定
  # 四项全部通过才可判定为 Yes；任一项不达标或缺失均为 No。
  is_all_pass <- pass_prot && pass_recall && pass_snr && pass_rc
  overall_quality <- ifelse(is_all_pass, "Yes", "No")
  
  # --- 4. 构建表格数据 ---
  # 定义列名
  cols <- c("样本组", "蛋白鉴定数", "标称特性召回率", "信噪比 (dB)", "Pearson相关系数", "是否通过")

  # 第一行：本报告采用的 2026 年质量标准
  row_std <- c("2026 年质量标准", ">1,000", ">0.90", ">5", ">0.80", "-")

  # 第二行：实际数据
  row_dat <- c(batch_name, txt_prot, txt_recall, txt_snr, txt_rc, overall_quality)

  df_rep <- data.frame(rbind(row_std, row_dat), stringsAsFactors = FALSE)
  colnames(df_rep) <- cols

  # 创建 Flextable
  ft <- flextable(df_rep) %>%
    theme_box() %>%
    flextable::font(part = "all", fontname = "Times New Roman") %>%
    align(align = "center", part = "all") %>%
    width(width = 1.2) %>%
    # 表头加粗背景灰
    bold(part = "header") %>%
    bg(part = "header", bg = "#EFEFEF") %>%
    # 第一行(标准)字体颜色变灰
    color(i = 1, color = "gray40") %>%
    # 第二行(数据)加粗
    bold(i = 2) %>%
    # 未通过的指标均标红，包括缺失值。
    color(i = 2, j = 2, color = ifelse(pass_prot, "black", "red")) %>%
    color(i = 2, j = 3, color = ifelse(pass_recall, "black", "red")) %>%
    color(i = 2, j = 4, color = ifelse(pass_snr, "black", "red")) %>%
    color(i = 2, j = 5, color = ifelse(pass_rc, "black", "red")) %>%
    # ----------------------------------------------------
    # 最后一列整体质量颜色 (Yes=黑, No=红)
    # 根据结果给"整体质量"上色
    color(i = 2, j = 6, color = ifelse(overall_quality == "Yes", "black", "red"))

  # --- 5. 文档文本内容 (来自提供的模板) ---
  txt_summary <- "本报告基于多项组学关键质量控制指标，总结了 Plasmix血浆参考物质所生成数据的质量情况。质量控制流程从用户输入血浆蛋白组表达矩阵开始，分别计算蛋白鉴定数、标称特性召回率 (Recall)、信噪比 (SNR) 和与参考数据集的 Pearson 相关系数 (PCC，内部指标名为 RC)。按照本报告采用的 2026 年质量要求，仅当蛋白鉴定数 >1000、Recall >0.90、SNR >5 dB 且 PCC >0.80 时，整体判定为通过；任一项等于或低于阈值、缺失或非有限值均判定为不通过。判定采用未四舍五入的计算值，表格数值仅为格式化展示。"

  txt_def_title <- "质量控制指标"
  txt_def_1 <- "蛋白鉴定数（Number of identified proteins）：按输入表达矩阵第一列唯一 Feature 的数量统计，沿用蛋白鉴定数名称。该计数口径不额外排除所有样本信号均缺失的 Feature；请在输入数据准备阶段完成相应的检测和鉴定筛选。"
  txt_def_2 <- "标称特性召回率（Recall of nominal characteristics）：定义为标称特性蛋白在测试数据集中被成功检测到的比例。"
  txt_def_3 <- "信噪比（Signal-to-Noise Ratio, SNR）：用于刻画某一检测平台、实验室或批次区分不同生物样本组之间内在生物学差异（“信号”）与同一样本组技术重复变异（“噪声”）的能力。SNR 越高，表明区分多组样本差异的能力越强。"
  txt_def_4 <- "与参考数据集的 Pearson 相关系数（Pearson correlation coefficient, PCC）：将输入的 Log2 信号还原到线性空间，对每个样本组的技术重复取算术均值，再计算测试样本组相对于参考样本组的比值。与参考数据集按特征和样本对匹配后，汇总计算比值表达谱的 Pearson 相关系数，用于表征相对定量的一致性。内部指标名沿用 RC。"

  txt_refs <- c(
    "1. Zheng Y, et al. Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials. Nature Biotechnology, 2024.",
    "2. Liu Y, et al. Harmonizing plasma proteomics data with the sample-to-reference ratio approach. BioRxiv (2026)",
    "3. 上海临床队列组学检测工作指引（征求意见稿）, 2025/11/26"
  )

  txt_disclaimer <- "本数据质量报告仅针对所评估的特定数据集提供分析结果，仅供信息参考之用。尽管已尽最大努力确保分析结果的准确性和可靠性，但本报告按“现状（AS IS）”提供，不附带任何形式的明示或暗示担保。报告作者及发布方不对基于本报告内容所采取的任何行动承担责任。本报告中的结论不应被视为对任何产品或流程质量的最终判定，也不应用于关键应用场景、商业决策或法规合规用途，除非经过专业核查和独立验证。对于分析结果的正确性、准确性、可靠性或适用性，不作任何明示或暗示的保证。"

  # --- 6. 生成 Word 文档 ---

  # 如果提供了模板，就加载模板；否则新建空白文档
  if (!is.null(report_template) && file.exists(report_template)) {
    doc <- read_docx(path = report_template)
  } else {
    doc <- read_docx()
  }

  # 文档构建顺序：
  # Title -> Summary -> Table -> Definitions -> Refs -> Disclaimer -> Plots

  # 标题
  doc <- doc %>%
    body_add_par("Plasmix血浆蛋白组质量报告", style = "heading 1") %>%
    body_add_par("", style = "Normal") # 空行

  # 摘要
  doc <- doc %>%
    body_add_par("摘要", style = "heading 2") %>%
    body_add_par(txt_summary, style = "Normal") %>%
    body_add_par("", style = "Normal")

  # 表格
  doc <- doc %>%
    body_add_flextable(ft) %>%
    body_add_par("", style = "Normal")

  # 定义
  doc <- doc %>%
    body_add_par(txt_def_title, style = "heading 2") %>%
    body_add_par(txt_def_1, style = "Normal") %>%
    body_add_par(txt_def_2, style = "Normal") %>%
    body_add_par(txt_def_3, style = "Normal") %>%
    body_add_par(txt_def_4, style = "Normal") %>%
    body_add_par("", style = "Normal")

  # 参考文献
  doc <- doc %>%
    body_add_par("参考文献", style = "heading 2")
  for (ref in txt_refs) {
    doc <- doc %>% body_add_par(ref, style = "Normal")
  }
  doc <- doc %>% body_add_par("", style = "Normal")

  # 免责声明
  doc <- doc %>%
    body_add_par("免责声明", style = "heading 2") %>%
    body_add_par(txt_disclaimer, style = "Normal") %>%
    body_add_break() # 分页

  # 图表部分
  # Signal-to-Noise Ratio
  doc <- doc %>% body_add_par("Signal-to-Noise Ratio", style = "heading 2")
  if (!is.null(qc_result$snr_plot)) {
    doc <- doc %>% body_add_gg(value = qc_result$snr_plot, style = "centered", width = 6, height = 5)
  } else {
    doc <- doc %>% body_add_par("信噪比图未生成，请查看质量指标表。", style = "Normal")
  }
  doc <- doc %>% body_add_par("", style = "Normal")

  # Correlation
  doc <- doc %>% body_add_par("Pearson Correlation Coefficient", style = "heading 2")
  if (!is.null(qc_result$cor_plot)) {
    doc <- doc %>% body_add_gg(value = qc_result$cor_plot, style = "centered", width = 6, height = 5)
  } else {
    doc <- doc %>% body_add_par("相关性图未生成，请查看质量指标表。", style = "Normal")
  }
  doc <- doc %>% body_add_par("", style = "Normal")

  # --- 7. 保存 ---
  print(doc, target = output_file)
  message(paste("Report generated successfully at:", output_file))
}
