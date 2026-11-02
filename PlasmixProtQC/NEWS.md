# PlasmixProtQC 1.1.1

- 统一两张图的字体与主题；SNR点缩小并增加细白边。
- PCC图采用相同横纵轴范围与等比例坐标，显示PCC及实际n。
- 撤回此前额外加入的正文说明，恢复原报告文字；参考数据和计算不变。

# PlasmixProtQC 1.1.0

- 内置参考集更新为FDU_v2.0_20261101，继续使用reference_dataset和reference_dataset_quali两个主要RDA对象。
- 正式报告采用固定M/F名单的线性比例PCC并显示实际n；保留S/P定值及按样本对读取的接口。
- Olink Recall按有限NPX报告覆盖计算，负NPX有效；SNR算法保持不变。
- 检出数>1000、Recall>0.90、SNR>5、MF PCC>0.80四项均满足才合格，阈值判断使用完整精度。
- 数据打包脚本保留固定名单及参考集版本；演示脚本不再覆盖真实RDA。
