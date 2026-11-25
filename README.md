# LIDC-IDRI 肺结节 CT 影像组学：良恶性评分无增量与医师分歧（灰区）有增量的双终点研究

**Radiomics on LIDC-IDRI pulmonary nodules: no size-adjusted increment over diameter for malignancy prediction, but a robust increment for predicting radiologist disagreement.**

## 摘要 / Abstract

在 LIDC-IDRI 的 302 例 / 602 个带医师标注结节的公开队列上，用 PyRadiomics 107 维特征 + 嵌套交叉验证（按病人分组）做了两个终点的对照实验：

1. **终点一（负结果）**：预测医师恶性评分（≥4 vs ≤2）时，影像组学 AUC 0.975 ± 0.003，但 log(直径) 单变量基线即达 0.980——该评分与结节直径近共线（Spearman ρ = 0.768），组学在此标签上无尺寸增量。在这个标签上宣称"组学 AUC > 0.9"属于把体积估计包装成恶性预测。
2. **终点二（正结果）**：预测医师分歧（评分 = 3 分的"灰区" vs 1/2/4/5 分"确信"）时，纹理特征 AUC 0.732 ± 0.007，直径基线仅 0.401（倒挂），size-adjusted 增量 ΔAUC = +0.334（95% CI +0.259 ~ +0.410）。该增量在两个独立子集中存活：≤15 mm 窄尺寸带 ΔAUC +0.149（CI +0.071 ~ +0.226）；≥4 位医师标注子集 ΔAUC +0.266（CI +0.124 ~ +0.396）。

**结论**：LIDC malignancy 评分携带的可预测信号大部分是尺寸；而"医师为什么犹豫"与尺寸近似正交，结节形态/纹理本身携带可复现的不确定性信号。提示以医师分歧为终点的组学研究比以主观恶性评分为终点更能分离形态学信息。

## 关键结果

| 分析队列 | n（结节/病人） | 纹理 AUC | 直径基线 | ΔAUC（95% CI） |
|---|---|---|---|---|
| 良恶性终点（负结果） | 322 / 206 | 0.975 ± 0.003 | 0.980 | ≈ 0（无增量） |
| 灰区 · 全样本 | 602 / 301 | 0.732 ± 0.007 | 0.401 | +0.334（+0.259 ~ +0.410） |
| 灰区 · ≤15 mm 窄带 | 471 / 246 | 0.671 ± 0.006 | 0.522 | +0.149（+0.071 ~ +0.226） |
| 灰区 · ≥4 位医师标注 | 184 / 34 | 0.706 ± 0.024 | 0.443 | +0.266（+0.124 ~ +0.396） |

方法：LASSO-L1 逻辑回归 + 嵌套 CV（外层 5 折 × 5 次重复，GroupKFold 按病人分组，折内 LASSO 选特征），病人级聚类 bootstrap 2000 次出 CI，置换检验 100 次验证显著性（p = 0.010）。完整结果文件见 [`results/`](results/)（各队列的 metrics.json / report.txt / roc.png / selected_features.csv）。

## 数据声明 / Data statement

- 数据来源：[LIDC-IDRI](https://www.cancerimagingarchive.net/collections/lidc-idri/)（TCIA 公开数据集）及其标注 XML；本仓库**不分发任何图像或衍生特征表**，使用者需自行从 TCIA 获取并遵守其数据使用协议。
- 仓库内 `results/` 仅含聚合统计（metrics.json）、图表、以及不含病人信息的特征入选频率表（selected_features.csv）。

## 仓库结构

```
├── download_lidc.py     # TCIA 原始通道下载 CT+SEG（支持代理、指定病例）
├── download_hf.py       # HuggingFace 镜像下载（分层随机抽样，推荐 --sample 100 --seed 42）
├── extract_case.py      # 逐例特征提取：每位医师标注一行，一阶 + 形状 + GLCM（<PID>_features.csv）
├── diagnose_match.py    # 逐例医师标注与结节元数据的匹配核验（ Dice 矩阵 / 匹配诊断）
├── build_model.py       # 核心建模：双终点 + 嵌套 CV + size-adjusted 增量检验 + 置换检验 + ROC
└── results/             # 各分析队列的完整结果产物（metrics.json / report.txt / roc.png / selected_features.csv）
```

## 复现步骤

```bash
# 0) 环境（开发环境：Python 3.9.23）
pip install -r requirements.txt

# 1) 数据获取（二选一）
python download_lidc.py --limit 100 --proxy http://127.0.0.1:7890   # TCIA 原始通道
python download_hf.py --sample 100 --seed 42 --verify               # HF 镜像，分层抽样

# 2) 特征提取（下载脚本对每例自动调用；也可单独执行）
python extract_case.py "E:/lidc_data/lidc_idri/LIDC-IDRI-0001"

# 3) 建模（--endpoint mal 良恶性终点 / grey 灰区分歧终点；--add-size 强制直径入模做增量检验）
python build_model.py --repeats 5 --n-perm 100 --out out_model
python build_model.py --endpoint grey --repeats 5 --n-perm 100 --out out_model_grey
python build_model.py --endpoint grey --add-size --repeats 5 --n-perm 0 --out out_model_grey_size
python build_model.py --endpoint grey --max-diameter 15 --add-size --repeats 3 --n-perm 0 --out out_model_grey_sub15
python build_model.py --csv features/nodules_features_gt4.csv --endpoint grey --add-size --repeats 5 --n-perm 0 --out out_model_grey_gt4_size
```

**已知复现缺口（如实声明）**：逐例 `<PID>_features.csv` 到逐结节汇总表 `nodules_features.csv` 的合并（拼医师评分、标注数、共识直径）当时以一次性 pandas 脚本完成，未随仓库发布。`build_model.py` 的输入列结构可从 `results/*/selected_features.csv` 与 `report.txt` 反推；补发布合并脚本在 TODO 列表中。

## 主要局限

1. LIDC malignancy 评分为医师主观 1–5 分，本研究报告"对评分的预测"而非"对真实恶性的预测"。
2. 标注数 ≤2 的结节占 50.5%，其"3 分"来自单人犹豫；≥4 标注敏感性分析部分缓解（但该子集仅 34 位病人，CI 宽）。
3. class_weight=balanced 改变概率量级，校准曲线未检验——概率用于临床决策前需补校准。
4. 单中心公开数据集，无外部验证队列。

## 致谢与 AI 使用声明 / Acknowledgements & AI Use Disclosure

- 中文：本研究的数据处理流水线搭建、代码调试与文稿整理过程使用大语言模型（Claude）作为辅助工具。研究问题与双终点设计由作者提出，全部数据获取与模型训练由作者本人执行，所有结果经作者逐项核对。作者对内容的正确性承担全部责任。依据 ICMJE 及相关期刊政策，AI 工具不列入作者署名。
- English: Large language models (Claude) were used as auxiliary tools for pipeline construction, code debugging, and manuscript preparation. The research question and dual-endpoint design were conceived by the author; all data acquisition and model training were performed by the author; all results were verified by the author, who takes full responsibility for the content. In accordance with ICMJE guidance, AI tools are not listed as authors.

## 引用本仓库 / How to cite

引用格式见 [`CITATION.cff`](CITATION.cff)。

**DOI: [10.5281/zenodo.24531302](https://doi.org/10.5281/zenodo.24531302)**（Concept DOI，代表所有版本；v1.0.0 版本 DOI：[10.5281/zenodo.24531303](https://doi.org/10.5281/zenodo.24531303)）

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24531302.svg)](https://doi.org/10.5281/zenodo.24531302)

If you use this repository, please cite it via `CITATION.cff` or the Zenodo Concept DOI above.

## Reproducibility note (2026-11-25)

A full rerun of `build_model.py --endpoint grey` reproduced the results to within AUC ±0.0016 (well inside the bootstrap CI width of ±0.043). Root cause of the residual wobble was `solver="liblinear"` being invoked without a `random_state`; both L1 calls are now seeded (`random_state=0`), and two consecutive smoke runs are bit-exact. The committed `results/` remain the authoritative v1.0.0 outputs.

## License

代码 MIT；数据遵循 LIDC-IDRI / TCIA 使用条款。
