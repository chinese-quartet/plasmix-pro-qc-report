# PlasmixProtQC: Plasmix 血浆蛋白组学质量控制工具包

**PlasmixProtQC** 是一个用于评估血浆蛋白组学数据质量的 R 包。它专为 Plasmix 参考物质设计，支持 **SomaScan**、**Olink** 和 **DIA** 等多种主流检测平台。通过计算关键质量指标（如 SNR、Recall、RC），并生成可视化的 Word 报告，帮助用户快速判断批次数据的可靠性。


## 功能特性

* **多平台支持**：自动适配 SomaScan、Olink 和 DIA 数据的特性（如 PCA 缩放逻辑）。
* **自动化质控流程**：一键计算所有关键指标。
* **交互式与报告输出**：既可在 R 环境中查看结果，也能生成标准化的中文 Word 质量报告。
* **核心指标**：
* **Recall (标称特性召回率)**：评估检测的覆盖深度。
* **SNR (信噪比)**：评估区分不同生物样本组的能力。
* **PCC (Pearson 相关系数，内部指标名为 RC)**：基于比值（Ratio-based）评估相对定量一致性，与内置参考数据集进行比对。

## 2026 年质量判定要求

Word 报告的“是否通过”列要求以下四项**同时达标**：

| 指标 | 通过条件 |
| --- | --- |
| 标称特性召回率 Recall | **>0.90** |
| 信噪比 SNR | **>5 dB** |
| 相对定量一致性 PCC（内部指标名 RC） | **>0.80** |
| 蛋白鉴定数 | **>1000** |

所有条件均为严格大于；等于阈值不通过。任一指标缺失、非有限（`NA`、`NaN`、`Inf`）或不达标时，整体结果为 `No`，只有四项全部达标时才为 `Yes`。不达标的数值显示红色及向下箭头，缺失或非有限值显示红色 `-`。

阈值比较使用未四舍五入的计算值；表格及图标题中的小数位数仅用于展示。SNR 使用 `10 * log10(组间平均距离 / 组内平均距离)` 的 dB 值。蛋白鉴定数沿用输入表达矩阵第一列唯一 Feature 的计数口径。


## 参考数据集与维护

当前参考集为 **FDU_v2.0_20261101**，适用于FDU材料的DIA、Olink HT和SomaScan 11K。参考集版本与R包软件版本分别管理。

继续使用两份主要内置RDA：`reference_dataset_quali.rda`用于Recall；`reference_dataset.rda`保留完整S/P和直接定值的MF，以及`UseForPCC`和`ReferenceVersion`。主报告依据参考集的UseForPCC标记纳入条目，DIA25、Olink134、SomaScan687个蛋白，不随提交批次重新筛选。PCC后括号n为实际有效交集。

两张图采用相同字体与主题。SNR散点缩小并增加细白边，保持样本组颜色；PCC散点单色显示PCC及实际n，横纵轴使用同一范围和等比例坐标。

需要更新参考数据时，替换仓库根目录两张CSV，再运行`Rscript data_genaration.R`重建这两份RDA；不要编辑二进制文件，也不要从S/P共识值相除生成MF。数据打包脚本与标准集构建代码分开，报告运行不会重新构建参考集。

Olink负NPX有效，Recall使用有限值报告覆盖；DIA/SomaScan沿用提交值>0。输入需为log2尺度。无提交LOD时不能把报告覆盖解释为高于LOD。秩转换数据不能按常规信号解释其比例PCC。

## 安装

软件包位于仓库的PlasmixProtQC子目录：

```r
install.packages("remotes")  # 已安装时跳过
remotes::install_github("chinese-quartet/plasmix-pro-qc-report", subdir = "PlasmixProtQC")
```

开发本地副本时，在仓库根目录运行`devtools::load_all("PlasmixProtQC")`。

## 数据准备

在使用前，请确保您的输入数据符合以下格式：

### 1. 表达矩阵文件 (Expression Matrix)

* 格式支持 `.txt` (制表符分隔) 或 `.csv`。
* **第一列**：必须是特征 ID (如 Protein Name, UniProt ID, SeqId)。
* **后续列**：样本的 Library ID。
* **数值**：建议输入仪器导出的 **Log2 转换后的信号值**（程序内部会自动处理还原计算）。

### 2. 元数据文件 (Metadata)

* 必须包含以下三列（列名不区分大小写）：
* `library`: 对应表达矩阵中的列名。
* `sample`: 样本类型，正式报告必须包含 **M和F**。P、X、Y可继续用于其他指标；P不作为主PCC分母。
* `platform`: 检测平台名称，必须准确填写为 **`SomaScan`**, **`Olink`** 或 **`DIA`** (这将决定参考数据集的选择和计算逻辑)。


## 快速开始

以下示例展示了如何对 SomaScan、Olink 或 DIA 数据进行质控分析并生成报告。

### 第一步：加载包与定义路径

```r
library(PlasmixProtQC)

# 定义输出目录
output_dir <- './output/'

# --- 场景 A: SomaScan 数据 ---
expr_file <- system.file("extdata", "plasmix_somascan_test_expr.csv", package = "PlasmixProtQC")
meta_file <- system.file("extdata", "plasmix_somascan_test_meta.csv", package = "PlasmixProtQC")

# --- 场景 B: Olink 数据 ---
# expr_file <- system.file("extdata", "plasmix_olink_test_expr.csv", package = "PlasmixProtQC")
# meta_file <- system.file("extdata", "plasmix_olink_test_meta.csv", package = "PlasmixProtQC")

# --- 场景 C: DIA 数据 ---
# expr_file <- system.file("extdata", "plasmix_dia_test_expr.csv", package = "PlasmixProtQC")
# meta_file <- system.file("extdata", "plasmix_dia_test_meta.csv", package = "PlasmixProtQC")
```

### 第二步：运行质控分析

使用 `qc_conclusion` 函数进行核心计算。该函数会自动执行数据清洗、标准化、统计计算和绘图。

```r
# 运行分析
qc_res <- qc_conclusion(
  exp_path = expr_file,
  meta_path = meta_file,
  output_dir = output_dir
)

# 查看分析结果概览
message("\n=== 分析完成 ===")
print("结论表:")
print(qc_res$qc_metrics_table)

# 在 RStudio 中预览关键图表 (可选)
# print(qc_res$snr_plot)  # 查看 SNR PCA 图
# print(qc_res$cor_plot)  # 查看 RC 相关性散点图
```

### 第三步：生成 Word 报告

使用 `generate_protein_report` 函数生成包含图表和详细解读的 `.docx` 报告。

```r
# 获取包内自带的报告模板路径
template_path <- system.file("extdata", "Plasmix_template.docx", package = "PlasmixProtQC")

# 生成报告
generate_protein_report(
  qc_result = qc_res,              # 上一步的分析结果
  report_template = template_path, # Word 模板路径
  report_dir = output_dir,         # 输出目录
  report_name = "Plasmix_QC_Report.docx", # 输出文件名
  batch_name = "20260108测试批次"   # (可选) 指定报告中显示的批次名称
)
```


## 输出结果说明

运行结束后，`output_dir` 目录下将生成以下文件：

1. **conclusion_table.tsv**: 包含各项指标得分的表格。
2. **pca_plot.png / pca_table.tsv**: SNR 分析的 PCA 图与坐标数据。
3. **Plasmix_QC_Report.docx**: 最终的质控报告，包含：
* **摘要与结论**：Pass / Fail 判定。
* **指标详情**：蛋白鉴定数、Recall、SNR、RC 的具体数值。
* **可视化图表**：PCA 散点图和 RC 回归图。


## 常见问题
* **关于有效交集**：PCC使用参考数据与提交数据的有效交集。缺少M或F、有效交集少于10个或零方差时，PCC缺失，整体判定不通过。

---

**维护者**: [Yanming Xie/Fudan University]
**许可证**: MIT
