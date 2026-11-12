# 公共数据 bulk 外推验证 —— 执行报告

> 日期：2026-11-08
> 方案来源：用户提交的「零湿实验表型外推方案」（H1 免疫抑制 / H2 纤维化 / H3 耐药）
> 结论一句话：**stroma 模块在 bulk 中的可迁移性 = desmoplasia 的可迁移性（ρ≈0.97，6 队列一致）；它对免疫排斥、CD8 浸润与 pCR 的全部关联在剥离纤维化后归零；neural-intrinsic 模块的「正交性」在原稿之外不可复现。**

---

## 1. 数据可用性核查（逐项判定）

方案列出 9 个公共数据集 + 3 项工具。逐项核对 GEO 元数据、平台注释与样本特性后：

| # | 资源 | 方案中的描述 | 实际核查结果 | 判定 |
|---|---|---|---|---|
| A1 | GSE213688 | 12 donor / 15 capture Visium，11 cancer capture 可用 | 15 张 Visium（7 CL-like、4 genomic-unstable TNBC、4 metaplastic），h5 + positions + scalefactors 齐全 | ✅ 可用（但非纯 TNBC，n 效能有限） |
| A2 | GSE293199 | TNBC PDX + Xenium 280-gene panel，做「伪空间」 | 单样本 TNBC **PDX**（OmiCLIP 基准）。PDX 间质为**小鼠来源**，人源探针不捕获 → **该数据集中不存在人 myCAF 变量** | ❌ 不可用（不是降级问题，是变量不存在） |
| A3 | ARTEMIS / MDA Xenium | 101 例 TNBC 穿刺，部分公开 | 仅 scRNA-seq 公开（PRJNA1041570），44 例空间数据未公开 | ❌ 不可用 |
| B1 | GSE176078 | ~100k cells，算 myCAF/iCAF/LRRC15⁺ | 26 例、2.7 GB 已在本地，含 `celltype_minor`（CAFs myCAF-like / MSC iCAF-like） | ✅ 已用完 |
| B2 | GSE243526 | TNBC scRNA，14 例 | 40 样本 **scATAC + scRNA 混合**，全乳腺癌类型，非 TNBC 专属 | ⚠️ 不符描述，未采用 |
| B3 | GSE279219 | TNBC scRNA，~14 例 | 24 例 TNBC **FMD（模拟禁食）临床试验 bulk RNA-seq**（GPL24676 NovaSeq），非单细胞 | ❌ 不可用 |
| C1 | TCGA-BRCA | TNBC ~160–200 例 | 1100 例有 RNA-seq；ER/PR/HER2 三阴 **n=122** | ✅ 可用（n 比预期小） |
| C2 | METABRIC | 已用，扩展交互项 | 1976 例，原稿已跑完 | ✅ 已有 |
| C3 | GSE163882 | n=222 TNBC NAC，pCR 标签 | **222 例全乳腺癌**，其中 **TNBC 仅 90 例（pCR 38 / RD 52）**；RRASeq（AmpliSeq 20k panel），标签齐全 | ✅ 主队列（但 TNBC 子集 n=90，非 222） |
| C4 | GSE76275 | TNBC 新辅助 bulk，响应标签 | 265 例 TNBC **未治疗**肿瘤（Burstein 2015 分型队列），**无 NAC、无 pCR 标签** | ❌ 不可用（无响应标签） |
| C5 | GSE25065 | TNBC NAC bulk，响应标签 | 198 例 HER2- NAC（紫杉+蒽环），pCR 42 / RD 140 / NA 16，含 PAM50、grade、stage | ✅ 复制队列 |
| C6 | GSE91061 | ICB 队列，跨适应症 | 109 例**黑色素瘤** anti-CTLA4/PD-1，非乳腺 | ⚠️ 跨适应症，未纳入主分析 |
| C7 | IMvigor210 | ICB 队列 | 仅有 R 包 `IMvigor210CoreBiologies`，本机**无 R** | ❌ 不可执行 |

**方案未列但更合适的资源（已发现）**

| 资源 | 内容 | 状态 |
|---|---|---|
| **GSE319641（NeoTRIP）** | 280 例高危 TNBC 随机接受 CT（n=142）或 CT+atezolizumab（n=138），401 个 baseline/D1C2 FFPE 转录组（TPM）已在本地 | ✅ 数据在本地；但 GEO 样本特性**只有 tissue/timepoint/patient/batch，无 pCR**，pCR 需从论文附表补 |
| GSE350571 | Xenium Prime 5K（5000 基因，含 UCHL1/TUBB3），6 例乳腺肿瘤 + 配对 Visium/Visium HD | 未下载（数 GB），可作为后续单细胞分辨率验证 |

**工具可行性**

- **CIBERSORTx**：需注册并申请 token，无公开 API → 未使用，改为本地 mean-z marker 打分。
- **TIDE**：仅网页提交，无法批处理 → 改本地 proxy（见 §3）。
- **MSigDB / GSVA**：本项目已有 MSigDB 缓存与打分模块，可用。

---

## 2. 对方案设计的五处修正

### 2.1 Stage 1 的 Go/No-go 阈值没有筛选力

方案设「stroma–CAF/ECM ρ ≥ 0.5 则通过」。原稿 METABRIC 已测得 **ρ = 0.975**，必然通过。本次 6 个队列全部 **0.964–0.973**（见 §3.1）。真正有信息量的判定应设在**偏相关残差**上：控制 CAF/ECM 后 stroma 模块剩余方差仅 1 − 0.97² ≈ **5%**，这决定了 Stage 2 的所有"矫正后"检验都几乎没有统计空间。

### 2.2 Stage 2 的预期结论不成立

方案预期「部分免疫抑制信号由 desmoplasia 共享，但 LRRC15/SPARC 子集在矫正后仍相关」。由于 ρ(stroma, CAF/ECM) = 0.97，**"矫正后仍相关"在数学上不可能**。实测结果方向完全相反：矫正后 stroma 残余与 CD8 **正相关**（+0.10 ~ +0.39），与免疫抑制假设相反。唯一可行的替代是把 LRRC15/SPARC 操作化为**对 CAF/ECM 回归的残差**再检验，而不是模块层面。

### 2.3 Stage 3 的组间比较设计不可执行

方案写「bulk stroma-module 高低组 ↔ 外部 scRNA 里 myCAF 比例」。bulk 与 scRNA 是**两套不同的样本集**，无法逐例对应，所谓"高低组"不存在。已改为在 scRNA 内部做 pseudobulk 相关：检验「模块评分 ↔ myCAF 细胞占比」（n=26）。

### 2.4 GSE293199 应完全剔除

PDX 的间质是小鼠细胞，人源 Xenium 探针不捕获；且该样本无神经束。所谓"把 Xenium 细胞按 55µm 网格聚成伪空间、反卷积复算 stroma 模块"在此数据上得到的是人肿瘤细胞信号，不含 myCAF。该分析不能补"蛋白级神经锚定"的短板。

### 2.5 外推的主角应该是 neural-intrinsic，而非 stroma

方案把边界划在「神经指令端（NGF/CGRP）不在此方案内」，但原稿的分解恰恰指出：stroma 模块与 desmoplasia 重合（ρ=0.975），neural-intrinsic 模块才与 desmoplasia 正交（METABRIC ρ=0.040）。若外推只用 stroma 模块，无论结果阳性阴性都只能得到"desmoplasia 的故事"。**本次分析因此双模块并行**——而结果（§3.1）显示 neural 模块的"正交性"本身也不可复现，这反而是更值得写进讨论的发现。

---

## 3. 结果

### 3.1 Stage 1 —— 模块迁移与正交性（6 队列）

Spearman ρ，括号内为 1000× bootstrap 95% CI。

| 队列 | n | stroma ~ CAF/ECM | neural ~ CAF/ECM | NPP58 ~ CAF/ECM |
|---|---|---|---|---|
| TCGA-BRCA 全队列 | 1100 | **0.973** (0.967–0.977) | 0.198 (0.140–0.257) | 0.370 (0.318–0.420) |
| TCGA-TNBC | 122 | **0.964** (0.946–0.974) | 0.252 (0.073–0.417) | 0.172 (−0.022–0.345) |
| GSE163882 全队列 | 222 | **0.973** (0.962–0.979) | 0.406 (0.283–0.512) | 0.394 (0.270–0.503) |
| GSE163882-TNBC | 90 | **0.970** (0.946–0.981) | **0.569** (0.408–0.690) | 0.434 (0.232–0.614) |
| GSE25065 全队列 | 198 | **0.972** (0.960–0.979) | 0.124 (−0.013–0.256) | 0.538 (0.427–0.634) |
| GSE25065-TNBC | 66 | **0.970** (0.942–0.979) | **−0.170** (−0.396–0.067) | 0.597 (0.390–0.739) |
| *METABRIC（原稿）* | 1976 | *0.975* | *0.040* | *0.54* |

**读法**

1. **stroma 模块 = desmoplasia 的重述，这一点跨平台稳如磐石**：Illumina HT-12（METABRIC）、Agilent（TCGA 芯片版对照）、Illumina RNA-seq（TCGA）、Affymetrix U133A（GSE25065）、AmpliSeq 靶向测序（GSE163882）五个平台一致落在 0.96–0.98。
2. **neural 模块的"正交性"不能复现**：原稿的 0.040（P=0.075，不显著）在另外 5 个队列变成 0.12 / 0.20 / 0.25 / 0.41 / **0.57**，且多数显著。这不是低检出基因造成的技术噪声——按检出率分层后，高检出子集的 ρ 反而更高（GSE163882：全 39 基因 0.406，检出 ≥50% 的 21 基因 **0.488**，检出 <50% 的 18 基因 0.117；TCGA：0.252 / 0.210 / 0.213）。
3. **对主稿的直接影响**：3.5 节与 Note 10 中「neural-intrinsic 模块正交于 CAF/ECM（ρ = 0.040）」应改写为「在 METABRIC（Illumina HT-12）中近乎正交，但在 RNA-seq/AmpliSeq 队列中 ρ 为 0.20–0.57，提示该正交性为平台或队列构成依赖，不能作为普适性质陈述」。这一条如果被审稿人独立发现会很被动，主动改写更稳。

### 3.2 Stage 2 —— 免疫抑制（5 队列）

偏相关控制 CAF/ECM。完整 53×5 行见 `stage2_immune_association.csv`，核心四行：

| 队列 | stroma ~ Excl_proxy（raw → 偏） | stroma ~ CD8（raw → 偏） | stroma ~ TGF-β（raw → 偏） | neural ~ CD8（偏） |
|---|---|---|---|---|
| TCGA-BRCA (n=1100) | 0.846 → **0.071** | 0.086 → **+0.205** | 0.715 → 0.366 | +0.115 |
| TCGA-TNBC (n=122) | 0.794 → **0.020** | −0.078 → **+0.215** | 0.670 → 0.331 | −0.161 |
| GSE163882 (n=222) | 0.846 → **0.014** | 0.102 → **+0.363** | 0.725 → 0.569 | +0.164 |
| GSE163882-TNBC (n=90) | 0.826 → **0.047** | −0.025 → **+0.394** | 0.704 → 0.562 | −0.033 |
| GSE25065 (n=198) | 0.740 → **0.103** | −0.045 → **+0.103** | 0.594 → 0.291 | −0.092 |

**结论**

- **H1（免疫抑制）被否**。免疫排斥 proxy 与 stroma 模块的 raw 相关高达 0.74–0.87，但**100% 由 desmoplasia 承担**：剥离 CAF/ECM 后塌到 0.01–0.10。
- **方向反转**：矫正后 stroma 残余与 CD8 正相关（+0.10 ~ +0.39，5/5 队列一致）。即"超出纤维化"的那部分基质程序伴随更多 CD8，而不是更少。
- **TGF-β 是唯一矫正后仍稳健的读数**（0.29–0.57），但 **TGFB1 同时属于 perineural-stroma 模块与 TGFb_core 打分**，存在基因共享所致的循环，不能作为独立证据；若需保留必须剔除共享基因后重算。
- **M2 / SPP1 巨噬细胞**：raw 0.37–0.63，矫正后仍 0.23–0.46（部分由 TGFB1 共享贡献）。
- neural 模块与 CD8 仅在 TCGA-TNBC 呈弱负相关（−0.16），GSE163882-TNBC 为 −0.03，方向不一致。

### 3.3 Stage 3 —— 单细胞 myCAF 锚定（GSE176078，n=26）

| 模块 | ~ myCAF 占全部细胞 | ~ myCAF 占 CAF | ~ CAF 占全部细胞 |
|---|---|---|---|
| Perineural_stroma | **+0.719 (P=3.5×10⁻⁵)** | +0.290 (P=0.15) | +0.733 |
| CAF_ECM | **+0.762 (P=6.1×10⁻⁶)** | +0.354 (P=0.076) | +0.750 |
| LRRC15_myCAF | **+0.748 (P=1.1×10⁻⁵)** | +0.307 (P=0.13) | +0.749 |
| Neural_intrinsic | +0.385 (P=0.052) | −0.012 (P=0.95) | +0.482 (P=0.013) |
| NPP58 | +0.393 (P=0.047) | +0.121 (P=0.56) | +0.465 (P=0.017) |

TNBC 亚组（n=10）：stroma +0.600（P=0.067）、CAF/ECM +0.648（P=0.043）、neural +0.261（P=0.47）。

**结论**

- stroma 模块确实读得到 myCAF **总量**（0.72），但与 CAF_ECM 的 0.76 几乎无差别 → 它读的是"CAF 有多少"，不是"myCAF 占多少"。
- **所有模块与 myCAF 占 CAF 比例均不相关**（0.29–0.35，全不显著）→ 模块不携带 CAF 亚型构成信息。原稿「模块归因于 myCAF 输出」的表述需要限定为「myCAF 数量」而非「myCAF 状态」。
- neural 模块与 myCAF 总量正相关（0.385，P=0.052），与 §3.1 的非正交性方向一致。

### 3.4 Stage 4 —— NAC 病理完全缓解（pCR）

logistic：logit(pCR) ~ 模块 z 分数（per SD）；AUC 为 5 折 × 5 次重复交叉验证。

| 模块 | GSE163882 全 (n=222, 80 pCR) | GSE163882 TNBC (n=90, 38) | GSE25065 (n=182, 42) |
|---|---|---|---|
| NPP58 | 0.81 (0.62–1.07), AUC 0.61 | 0.86 (0.56–1.31), AUC 0.55 | 0.92 (0.65–1.30), AUC 0.42 |
| Perineural_stroma | 1.04 (0.79–1.37), AUC 0.57 | 1.04 (0.68–1.58), AUC 0.51 | 0.82 (0.57–1.17), AUC 0.53 |
| Neural_intrinsic | 0.96 (0.72–1.28), AUC 0.58 | 0.98 (0.65–1.49), AUC 0.51 | 1.07 (0.76–1.51), AUC 0.42 |
| CAF/ECM | 0.97 (0.74–1.28), AUC 0.58 | 0.95 (0.62–1.44), AUC 0.52 | 0.78 (0.54–1.12), AUC 0.54 |
| LRRC15 myCAF | 0.92 (0.69–1.22), AUC 0.59 | 0.94 (0.62–1.44), AUC 0.53 | 0.74 (0.51–1.07), AUC 0.56 |
| **CD8_cytotox** | **1.69 (1.21–2.36), P=0.002, AUC 0.68** | **1.70 (1.03–2.80), P=0.037, AUC 0.67** | 1.18 (0.84–1.66), AUC 0.53 |
| TLS | 1.30 (0.99–1.72), P=0.057, AUC 0.63 | 1.27 (0.82–1.95), AUC 0.59 | **1.62 (1.15–2.30), P=0.006, AUC 0.62** |
| TGF-β core | 1.41 (1.06–1.86), P=0.017, AUC 0.62 | 1.35 (0.88–2.08), AUC 0.55 | 0.87 (0.61–1.22), AUC 0.48 |

**结论**：**H3（耐药）被否**。没有任何空间先验模块预测 pCR；三个队列中唯一复现方向的是 CD8 细胞毒（GSE163882 两个子集阳性，GSE25065 同向但未达显著）。这与原稿"唯一稳定的独立信号是细胞毒 T 细胞含量"一致。

### 3.5 Stage 6 —— 随机基线校准（510 组 size-matched 基因集）

| 队列 | 模块 | 观测 AUC | 随机集中位数 | 随机集 P95 | 经验 P |
|---|---|---|---|---|---|
| GSE163882-TNBC | Perineural_stroma | 0.493 | 0.564 | 0.619 | 0.971 |
| GSE163882-TNBC | Neural_intrinsic | 0.485 | 0.558 | 0.632 | 0.875 |
| GSE163882-TNBC | CAF/ECM | 0.462 | 0.546 | 0.624 | 0.633 |
| GSE163882-TNBC | NPP58 | 0.484 | 0.567 | 0.629 | 0.900 |
| GSE25065 | Perineural_stroma | 0.454 | 0.524 | 0.639 | 0.529 |

**这是最硬的一条阴性证据**：模块的观测 AUC（0.45–0.51）**低于**同规模随机基因集的中位数（0.51–0.57）。即"用这些基因算分"并不比随便挑同样多的基因更好——签名在 pCR 上无信息量。

---

## 4. 对主稿的修改建议

| 原稿位置 | 现表述 | 建议 |
|---|---|---|
| 3.5 节 / Note 10.4 | 「neural-intrinsic 模块与 CAF/ECM 正交（ρ = 0.040，不显著）」 | 限定为「在 METABRIC（Illumina HT-12）中近乎正交；在 TCGA RNA-seq 与 AmpliSeq 队列中 ρ = 0.20–0.57，提示平台依赖」 |
| 3.7 节 / Note 12 | 「模块归因于 myCAF 输出」 | 限定为「myCAF 数量而非亚型构成」；补一句"模块与 myCAF 占 CAF 比例不相关（ρ = 0.29，n=26）" |
| 讨论 | 「bulk 代理被 desmoplasia 主导」 | 可强化：新增 6 队列 ρ = 0.96–0.98 的外部证据，以及免疫排斥 proxy 矫正后归零（0.85 → 0.07） |
| 新增 | — | 建议以「空间先验模块在 bulk 中不预测 NAC 缓解」作为一条独立的阴性结果写入讨论或补充材料（三队列 + 510 随机基线） |

---

## 5. 产出文件

```
01_data/raw/tcga_brca/tcga_brca_target_genes.csv     TCGA 178 基因 × 1100 样本（cBioPortal 逐基因）
01_data/raw/tcga_brca/tcga_brca_clin.csv             ER/PR/HER2/PAM50/OS
01_data/raw/gse163882/                               222 样本 TPM + series matrix 标签
01_data/raw/gse25065/                                198 样本表达 + pCR/RCB/PAM50 标签

02_scripts/bulk_lib.py                               评分/偏相关/队列装载
02_scripts/run_bulk_stage1.py                        Stage 1
02_scripts/run_bulk_stage2.py                        Stage 2
02_scripts/run_bulk_stage3.py                        Stage 3
02_scripts/run_bulk_stage4.py                        Stage 4 + 6
02_scripts/_fetch_tcga_cbioportal.py                 TCGA 抓取

03_results/tables/bulk_extrapolation/
    stage1_module_correlations.csv                   6 队列 × 9 组模块相关 + bootstrap CI
    stage1_partial_correlations.csv                  偏相关（控制 CAF/ECM）
    stage1_neural_detection_sensitivity.csv          neural 正交性的检出率分层
    stage2_immune_association.csv                    5 队列 × 53 关联（raw + partial）
    stage2_tertile_effects.csv                       三分位 Hedges g
    stage3_scrna_anchor.csv / stage3_scrna_per_sample.csv
    stage4_pcr_association.csv                       3 队列 × 8 模块
    stage6_random_baseline.csv                       510 组随机基线
```

## 6. 未完成 / 需决策

1. **GSE213688（Visium 外部复现）**：650 MB RAW.tar 已部分下载，未解包。作为"同平台复现"仍有价值（可验证神经邻近 spot 比例与衰减长度），但不影响本报告的结论。
2. **GSE350571（Xenium Prime 5K）**：唯一能真正回答"LRRC15⁺ myCAF 是否在神经束 150µm 内"的公开数据（5000 基因含 UCHL1/TUBB3 + 配对 Visium HD）。体积数 GB–数十 GB，需确认是否拉取。
3. **GSE319641（NeoTRIP）**：本地已有 401 样本 baseline/D1C2 TPM。可做"治疗诱导的模块变化"（D1C2 vs baseline），但无 pCR 标签，需从论文附表补。
4. **蛋白级验证（mIF/CyCIF）**：公共资源中无任何含神经标记的乳腺 mIF/CyCIF 队列（OHSU 102 例 42-marker 面板不含 PGP9.5/S100/TUJ1；复旦 465 例无公开数据）。此路线只剩湿实验，或降级为"补充材料给 12 色 panel 设计 + 分析代码 + limitation"。


---

# 追加（同日）：下游链执行结果（Stage 7–9）

用户指示"按你的执行"，随后完成了正交化、纤维化定量、生存与稳健性分析。

## 1. 正交化（Stage 7，`run_bulk_stage7_orthogonal.py`）

- 模块被 CAF/ECM 解释的秩方差比例：Perineural_stroma **93.0–94.6%**（5 队列）、LRRC15_myCAF 93.5–96.1%、SPARC_POSTN 81.8–86.4%、COLLAGEN_MATURATION 51–68%、Neural_intrinsic 1.5–32.4%、NPP58 2.9–28.9%。
- 三级衰减（raw → +CAF/ECM → +纯度代理）：**纯度代理（EPCAM/KRT 上皮评分）不是混杂**；TGF-β 轴保留 50–80%（残差 ρ 0.29–0.57，5/5 队列 BH q<1.5e-3）；M2/SPP1 次之；**CD8 无负相关**（raw −0.08~+0.10）→ 空间 CD8 减少是局部几何现象，bulk 不排斥。
- TCGA-BRCA 生存（n=761，96 事件）：NPP58 残差 HR 0.80 (p=0.043，校正后不显著)。
- 修复 statsmodels PHReg 在 HR≈1 时 bse 退化的 bug（`bulk_lib.cox_breslow` 自实现 Newton-Raphson，与 PHReg 在正常输入下 HR 一致到 4 位小数）。

## 2. 纤维化"量 vs 质"（Stage 8，`run_bulk_stage8_fibrosis.py`）

- 8 个基质程序全部对 OS 做 M0→M3 递增 Cox：**只有 SPARC/POSTN 核心残差**显著（HR 1.297 [1.073,1.567]，p=0.0072，校正 age/node/Basal/CAF总量）；CAF 总量、胶原核心/成熟、机械转导、EMT、TGF-β 程序全部 null（p 0.28–0.98）。
- TCGA-TNBC（19 事件）同向（HR 1.51–1.53，p≈0.052）。
- 单基因分解：无单基因显著 → 程序级信号。

## 3. 稳健性（Stage 9，`run_bulk_stage9_robustness.py`）

- +MKI67 校正：HR 1.281 (p=0.016) —— 非增殖混杂。
- Bootstrap 1000×：HR 1.259 [1.061, 1.511]。
- **510 组随机 8 基因集基线：经验 P = 0.080（TCGA-BRCA）、0.38（TNBC）→ 未过双校准门槛，稿件定为 exploratory**。
- KM 上三分位 vs 其余：log-rank p=0.37；5 年 OS 77.0% vs 81.3%。

## 4. 耐药（pCR）链（Stage 7 残差版）

- 原始模块全部不预测 pCR（Stage 4，随机基线 P 0.42–0.98）。
- 残差版：PNI-stroma 残差 3 队列方向一致（OR 1.32/1.54/1.16）但均 p>0.05；SPARC/POSTN 残差 GSE163882 显著（1.37*/1.68*）但 GSE25065 反向（0.80）→ **按复制失败处理**。

## 5. 稿件更新

- `manuscript_main.md`：3.5 / 摘要 / Discussion 4.2 / Key Points 的"neural-intrinsic 与 CAF/ECM 正交（ρ=0.04）"全部限定为 METABRIC 特异（外部 0.12–0.57）；4.3 节新增下游链一段（TGF-β 残留伴联、CD8 无负相关、基质总量预后 null、SPARC/POSTN 残差 exploratory HR 1.30）；核心正文压回 **5,999 词**（≤6,000 硬线），摘要 315。
- 新增 `supplementary/Supplementary_Note_13_bulk_extrapolation.md`（N13.1–N13.8）。
- 新增图 `SuppFigS13_bulk_extrapolation.png/pdf`（6 panel：正交化 / TGF-β 存活 / pCR 森林 / OS 森林 / KM / 随机基线）。
- `Supplementary_Tables.xlsx` 新增 sheet S22–S27（stage1/4/6/7/8/9 全部结果表）。
- 主稿 docx 重生成为 `manuscript_main_BIB_v2.docx`（原 BIB.docx 被编辑器占用，关闭后可将 v2 改名覆盖）。

## 6. 三条链的最终结论

| 链 | 结论 |
|---|---|
| PNI ↔ 纤维化 | bulk 里 94% 不可分；唯一独立的"质"信号 = SPARC/POSTN 核心（探索性） |
| PNI ↔ 免疫抑制 | 非 CD8 排斥；是 TGF-β 程序 + M2/SPP1 极化（5 队列 q<0.05 复现） |
| PNI ↔ 耐药 | 模块与残差均不能预测 pCR；CD8 细胞毒是唯一可复制阳性（AUC 0.67） |


---

# 追加（同日 23:15）：Stage 10–12 —— 免疫区室 / 中介 / METABRIC / 跨癌种

针对用户第二轮清单，逐项对账后补跑真正缺失的 4 块。

## 0. 清单纠错（用户"未做"清单 vs 实际）
- pCR（GSE163882）**已做**（Stage 4：全 null，随机基线 P 0.42–0.98；TNBC n=90 而非 222）
- bulk 随机基因集基线**已做**（Stage 6/9，各 510 组）
- LOX/YAP-TAZ/EMT/TGF-β 程序定量**已做**（Stage 8 的 PROGRAMS；LOX 残差 HR 1.20 p=0.053）
- myCAF 单细胞锚定**已做**（Stage 3）
- CIBERSORTx/xCell/TIDE：**均不可离线运行**（token/网页/R），以本地等价实现替代并如实标注

## 1. Stage 10 免疫区室（`run_bulk_stage10_deconv.py`）
- **反卷积不可识别**：256 基因靶向面板 + GSE176078 自建 10 区室签名，三种标度化 NNLS 全部退化（免疫权重→0 或均匀分布，R²<0.02）。改用区室特异标记富集分（117 特异基因，xCell 式）。
- **效度验证**：富集分 vs 配对单细胞真实占比 ρ=0.79–0.92（26 例 pseudobulk，全 P<2e-6）。
- 组成：stroma 模块与所有白细胞区室**正相关**（TCGA-BRCA：内皮 +0.58、巨噬 +0.40、CD4 +0.20、CD8 残差后 +0.28）→ 场是"白细胞丰富+基质丰富"，不是排斥。
- **中介**（2000 bootstrap）：残差→CD8 总效应为正（c=+0.31/+0.36）；巨噬中介 64–72%、Treg 61–84%、CD8T 75–88%；TGF-β 仅 26%（TNBC 9%，CI 含 0）。无抑制可中介。

## 2. Stage 11 METABRIC（`run_bulk_stage11_metabric.py`）
- 交互项全 null（HR_int 0.97–1.05，p 0.12–0.88）。
- 亚型分层：Basal 0.95 [0.79,1.15]、LumA 0.91 p=0.06，无一亚型显著。
- 构建 bootstrap：stroma ρ=0.976 [0.955,0.982] —— "签名≡纤维化"不是选基因的偶然；neural 0.036 [-0.050,0.119]，而**随机基因集中位 0.352** → 神经模块比随机集更正交。

## 3. Stage 12 随机基线正交性 + 跨癌种（`run_bulk_stage12_pancancer.py`）
- **核心方法学发现**：ρ 不能和 0 比，要和同队列 size-matched 随机集比。7 队列中 neural 模块在 5 个显著低于随机基线（z=-2.0~-11.1），GSE163882/GSE25065 无法区分；stroma 模块在全部 7 个高于基线（z=+6.5~+10.0）。主稿 3.5/摘要/Discussion 已改为"相对正交"表述。
- 跨癌种主链复现：CAF/ECM 解释 94.6–96.2%（PAAD 96.0%、HNSC 96.2%）；残差~TGF-β 0.25–0.58；残差~CD8 全为正（0.12–0.39）；残差对 OS 全 null（HR 0.70–1.02）。

## 4. 产物
- 新脚本 4 个；新数据 `01_data/raw/{brca,paad,hnsc}_tcga/`（cBioPortal 逐基因面板 + OS 临床）
- 新表 stage10_*/stage11_*/stage12_*（9 个 CSV）；xlsx 新 sheet S28–S31
- 新图 **SuppFigS14_immune_pancancer.png/pdf**（A 随机基线 z / B 区室效度 / C 中介森林 / D 跨癌种热图）
- Note 13 扩至 N13.12（新增 N13.9–N13.13，改写 N13.2 为"相对正交"）
- 主稿：摘要/3.5/4.2/4.3 更新为 null 校准表述 + 中介 + 跨癌种一句；core 5,998 词
- docx：`manuscript_main_BIB.docx` 与 `_v3.docx`（23:10 同步）；`_v2.docx` 为旧版（被预览面板占用未能覆盖）

---

# 追加（2026-11-09）：Stage 13 —— GSE91061 ICB 探索 + GSE350571 Xenium 启动

## 0. 用户方案对账（第三轮）
用户提交的"下一步方案"6 个方向中 5 个已完成（pCR / 反卷积等效实现 / METABRIC / 跨癌种 / 随机集基线；GSE76275 与 GSE243526 经核查数据不符描述，已剔除并记录）。真正遗留 3 件：GSE91061 ICB 探索（本轮完成）、OSF 分析计划草稿（已完成 `04_manuscript/osf_analysis_plan.md`）、GSE350571 Xenium 空间验证（下载中）。

## 1. Stage 13 —— GSE91061 ICB 探索（`run_bulk_stage13_icb.py`）
- 队列：Riaz et al. Cell 2017（PMID 29033130），BMS038 nivolumab±ipilimumab 黑色素瘤，109 样本/65 患者；GEO 无臂别/无生存 → 终点仅 RECIST（PD/SD/PRCR）。Pre-treatment 可评估 n=49（10 响应）。
- NCBI FTP 当日两度全站故障（500/503），以单连接断点续传完成 16.4MB FPKM；12 并行分块触发限流，教训：对 NCBI 并行度 ≤4。
- Entrez→symbol 经 mygene.info（symbol→entrez 反向查询，192 个模块基因命中 191，缓存 `symbol2entrez.csv`）。
- **结果**：全部 PNI/stroma 模块 OR 跨 1（Perineural_stroma OR 1.06 [0.37–3.03], AUC 0.52, P=0.91）；510 随机集基线经验 P ≥0.22；配对 pre→on delta 全 null。唯一超随机基线的读数 = CD8 细胞毒（AUC 0.62, 经验 P=0.016）——与化疗 pCR 链结论完全同构。Treg 名义 P=0.047 不过 Fisher/BH，按噪声处理。
- 方案中的"TNBC-like basal 筛选"不可行：黑色素瘤 basal 角蛋白 IQR 宽度仅 0.10，已在 Note 13.14 声明。
- 稿件整合：Note 13 新增 N13.14/N13.15；主稿 3.8 ICB 段改写为"no TNBC cohort; null cross-indication probe"（4.3 与 4.4 同步微调，词数守恒）；xlsx 新增 S32。

## 2. GSE350571 —— Xenium Prime 5K 神经邻近验证（进行中）
- 面板核查（go/no-go 通过）：神经标记 15/16 在面板（UCHL1/TUBB3/NCAM1/MPZ/NEFH/PRNP/SOX10/GAP43/RET/NTRK1/NGF/NGFR/CALCA/TAC1/NRP1，仅缺 S100B）；模块覆盖 Perineural_stroma 47/60、Neural 26/39、CAF_ECM 10/14、LRRC15 myCAF 6/8；免疫读数 Treg/M2/SPP1/TLS 全齐。COL1A1/2、COL3A1、SPARC、CD3D、KRT18 确实不在面板（Prime 5K 设计取舍）。
- 样本：6 例 FFPE 乳腺（P17/P22/P24/P26/P29/P32），TLS-rich 选择性取样（偏倚已声明）。Zenodo 细胞级数据被墙，改走 GEO 每样本单文件下载（跳过形态学图像）。
- 分析脚本已就绪：`run_bulk_stage14_xenium.py`——神经束 = 神经纤维转录本密度聚类核（UCHL1/TUBB3/NEFH/PRNP/GAP43/MPZ；排除 SOX10——乳腺中标记 basal 肿瘤；排除 NCAM1——NK 假阳性）；构成比 prox(≤150µm)/dist(≥300µm) + 环形平移 999 次置换；per-cell 模块评分（GeneGetter 惰性取列防 dense 化）。
- 当前状态：P17 cells.parquet 完成（7.3MB）；h5/transcripts 下载中（NCBI 限速，预计数小时）。

## 3. OSF 分析计划
`04_manuscript/osf_analysis_plan.md`：Stage 1–14 全链结构化（假设/数据/统计/多重性规则/占位符），明确标注 post-hoc 透明记录而非前瞻预注册。

---

# 追加（2026-11-09 晚）：Stage 15/16 —— 治疗语境生存 + 基因组桥

用户问"还能不能做转移、免疫、复发等表型"，对账后选定两项（用户确认 ①+②）。

## Stage 15 — GSE25065 NAC 治疗语境 DRFS（`run_bulk_stage15_gse25065_drfs.py`）
- 队列：198 HER2- NAC（taxane+anthracycline），DRFS 45 事件/中位 3.2 年；四层：ALL/TNBC-IHC(n=64)/RD-only(n=140)/pCR-only(n=42,3 事件)
- **结果**：全部 PNI/stroma 模块在所有层 null（HR 1.09–1.25，随机基线经验 P 0.33–0.65）；TNBC 层 HR 0.84–0.99。唯一信号 = CD8 细胞毒**反向**（ALL HR 1.33 P=0.048；RD-only HR 1.52 校正 P=0.016）——与 GSE12276 CD8 脑转移 HR 1.43 同向，与"残余病中 TIL 不良"的临床文献一致；与 CD8 预测 pCR 获益并存 = "CD8 预测响应，非响应者中标记残存侵袭病"的一致图景
- 修 bug 两处：AJCC 罗马数字前缀匹配把 IIIA-C 吞成 II（近恒定协变量→完全分离 HR=1.9e11）；RD 层恒定协变量剔除
- 输出：stage15_drfs_cox.csv、stage15_random_baseline.csv

## Stage 16 — TCGA-BRCA 基因型桥（`run_bulk_stage16_genotype.py` + `_fetch_tcga_mutations.py`）
- 数据：cBioPortal DETAILED 抓 60,864 条突变（分页），MUTATION_COUNT 与逐条非同义计数 ρ=1.000 交叉验证；982 样本 × TP53(301)/PIK3CA(317)/BRCA1(13)/BRCA2(15)
- **TMB**：stroma 模块原始负相关（LRRC15 −0.11、TGFb −0.14）在纯度+MKI67 校正后全部消失（|ρ|≤0.045）= 肿瘤含量伪影；唯一存活 = CD8(+0.10)/IFNγ(+0.10) 偏相关（P≈0.002，复制 TMB→免疫浸润经典轴，阳性对照成立）
- **驱动基因**：TP53-mut → LRRC15 myCAF 降低(q=0.028)+免疫读数低（已知）；PIK3CA-mut → stroma 模块升高（q=1.7e-5）但面板无 ESR1/ERBB2/PGR 无法排亚型混杂，仅描述性；BRCA1/2 合并(n=27) → CD8 升高（Δ=+0.32 P=0.005，HRD 生物学复制）
- **交互**：stroma × TP53 → OS HR 1.02 P=0.87
- 输出：stage16_tmb_correlation.csv、stage16_mutation_status.csv、stage16_tp53_os.csv

## 稿件整合
- Note 13 新增 N13.16/N13.17/N13.18；主稿 4.3 扩一句（DRFS null + 基因组轴 null）；xlsx 新增 S33
- 词数：发现上轮遗留超线（6,117），系统性压缩 4.2/3.6/3.3/3.5/3.8 冗余后回到 **6,000（硬线）**
- docx 重新生成：manuscript_main_BIB.docx = manuscript_main_BIB_v4.docx（v3 仍被预览占用）

---

# 追加（2026-11-09 深夜）：pCR 结论升入主稿正文

用户指出"就差一个 TNBC 专属的治疗响应验证"——对账确认 GSE163882 pCR logistic（Stage 4）早已完成（Note 13.4、S22-S27、随机基线经验 P 0.42-0.98），但**主稿正文从未给出 GSE163882 的数字**（只散在补充材料）。本轮修正：4.3 段落显式写入 TNBC n=90 的 stroma 模块 pCR null（OR 1.04, 0.68–1.58; AUC 0.51 = 随机集中位）与 CD8 唯一阳性（OR 1.70; AUC 0.67），与 DRFS/基因组轴 null 并列成"化疗耐药全链 null"一句。3.2/4.4 等量压缩腾挪词数，core 回到 6,000。docx 重新生成为 manuscript_main_BIB_v5.docx。

---

# 追加（2026-11-09 深夜二）：第四轮方案 8 条对账 + Stage 17-19 + 校准勘误

## 0. 用户第四轮清单对账
8 条中：④ GSE163882 pCR 已完成（Stage 4/N13.4/主稿 4.3）；⑤ GSE194040 探测后发现为 I-SPY2-990 资源（987 例/10 臂/含 pembro n=69+对照 210/pCR 标签）→ 升级为 Stage 20 在途；①②③⑥⑦⑧ 本轮完成（Stage 17/18/19）。

## 1. 关键勘误：axis 层随机基线 bug（N13.22）
- **bug**：`run_axis_gradient.py`（及新 stage17 初版）中 `agg(size=("rho","size"))` 把聚合函数 "size"（=样本计数）当成了基因集大小 → bucket 按样本数分组、size-matching 失效、所有经验 P 撞分辨率地板。**原稿主链（511 组 Note 8 主校准）解析基因集大小的代码正确，不受影响。**
- **受影响**：Note 12 表格、TableS20（且 S20 原为 30 例中期版本）、主稿 3.7 Line 1 与结论段"模块超过 size-matched 基线"的声明。
- **修正后**：单基因标记（SPARC/TNC/ACTA2/POSTN/LRRC15, FDR 0.030）、胶原核心（0.042）、apCAF 轴（0.030）仍超基线；iCAF 边缘（0.070）；**60 基因模块（FDR 0.145）与 20 基因 myCAF 轴（0.238）不再超过**——集合级平均本身即足以产生 −0.06 平滑梯度（基线 |ρ| 0.012@1 基因 → 0.059@60 基因）。CD8_T/pan-T 描述性正方向（远离神经）超基线（0.012, 16/17 例）。
- **稿件动作**：Note 12 表格+叙述重写（含勘误声明）；TableS20 以 43 例修正版重生成；主稿 3.7 Line 1 与结论改为"perineural 特异超额由单个 matrisome 基因携带，集合级处于基线"——与 bulk 层模块≡desmoplasia（ρ=0.97）互为空间镜像，反而强化中心论点。

## 2. Stage 17（空间功能表型，`run_bulk_stage17_tlscaf.py` + 修正校准）
- **TLS-CAF 共编程：否定**。CAF-TLS 输出（CXCL12/FAP/POSTN）原始梯度 −0.031 不超随机基线（P=0.098）；TLS 标记场可测 16/43，原始距离上略远离神经（+0.037）；耦合检验：TLS+ spots 内 CAFout~dist = +0.06（远离方向 P=0.005），远端 TLS+ 的 CAF 输出高于近端（+0.071 vs −0.006），near-vs-far Fisher-z null（P=0.83）。主稿"TLS 近神经富集"是 CAF 含量校准后（相对其基质背景）的陈述，与原始距离不矛盾，已在 N13.19 澄清。
- **胶原交联成熟度：否定**。Coll_mat（9 基因）−0.026 P_dir=0.49；LOX-only −0.010 P=0.90；EMILIN1/FBN1 均低于 size-1 基线尾部。"成熟胶原更陡衰减"不成立。
- **缺氧：唯一方向性发现（反向）**。10 基因缺氧评分随距离**升高**（+0.043；核心 4 基因 +0.041；CA9/SLC2A1/VEGFA +0.027~0.031，各 P_dir=1/51 地板，族 FDR 0.051）→ 近神经 150µm 域**不是**缺氧纤维化核心，是血管化 perineural cuff 样结构。

## 3. Stage 18（bulk CD8 效应/趋化解耦 + EMILIN1，`run_bulk_stage18_decouple.py`）
- **解耦假设：否定**。CYT/CD8/IFNγ 与趋化轴（CXCL9/10/11、CCL5、CCL19/21、CXCL13）在全部三队列同向正相关（偏 ρ +0.10~+0.39）。GSE25065 中趋化轴超随机基线（+0.256, empP=0.020；CXCL9/10/11 +0.237, P=0.046）。
- **GSE163882 平台警示**：随机集中位本身 +0.30（AmpliSeq TPM 的文库组成共享因子）——任何未校准偏相关在该平台都不可解读，本身是方法学发现。
- **EMILIN1 不复制**（GSE163882 −0.16 名义不超基线；GSE25065 +0.18 empP=0.14）。METABRIC 随机池为 210 靶向基因（偏倚），仅报告偏相关。

## 4. Stage 19（Schwann 评分投影，`run_bulk_stage19_schwann.py`）
- SC 评分（SOX10/MPZ/MBP/PMP22/PLP1/NCAM1/PRNP/NEFL/NEFH/S100B）与 neural-intrinsic 模块 ρ=0.52–0.77（三队列），远高于与 stroma 模块（0.25–0.53，仅 TCGA 超基线 z=+2.6~2.9）。
- SOX10 单基因与 basal 角蛋白 ρ=0.62–0.74 → bulk 单基因"神经标记"是亚型混杂。
- 结论：bulk "神经信号" = Schwann/胶质成分（neural 模块）+ 基质成分（stroma 模块），均非神经元固有信号。

## 5. 稿件动作
- Note 13 新增 N13.19–N13.23（含 N13.22 勘误）；Note 12 表格/叙述重写；TableS20 重生成（43 例修正版）；xlsx S34。
- 主稿：3.7 Line 1 与结论段修正（模块不再声称超基线）；4.3 新增三句（无解耦/缺氧反向/Schwann 分解）；2.3/2.4/2.8/3.1/3.3/3.8/4.2 系统压缩；core = 6,000（硬线）。
- docx：待重新生成为 v6。

## 6. 在途
- GSE194040 I-SPY2（Stage 20）：基因级矩阵后台下载中（~12.5MB+），临床元数据抓取中。
- GSE350571 Xenium（Stage 14）：P17 h5 57M/113.6M，持续限速。

---

# 追加（2026-11-10 凌晨）：Stage 20 I-SPY2 完成 + 全部整合

## Stage 20 — I-SPY2 chemo vs IO pCR（`run_bulk_stage20_ispy2.py`）
- 数据：GSE194040 family.soft（988 样本/987 患者/14 臂）+ 基因级矩阵 19,134×988（49.3MB，FTP range-append 下载）；患者号对齐 99.8%。
- **TNBC 全臂（n=362, 142 pCR）**：stroma 模块 null（OR 0.90, AUC 0.47=随机中位 0.435，经验 P=0.998）；neural-intrinsic 名义保护方向（OR 0.75, P=0.009，未做随机校准，仅记录）；**CD8 又是唯一阳性（OR 1.51, 1.21–1.89, AUC 0.61）**。
- **分臂（TNBC）**：紫杉醇单药（n=85, 17 pCR）stroma OR 1.48 (0.86–2.56)——"纤维化耐药"方向但不显著；pembro（n=29, 19 pCR）stroma OR 0.42 (0.17–1.07)、CAF/ECM OR 0.35 (0.12–0.98)、CD8 OR 3.33 (AUC 0.74)。
- **stroma × pembro 交互（n=114）：OR 0.28 (0.10–0.84), P=0.023**——名义显著：化疗臂的 stroma 劣势在 IO 臂衰减并反转。与用户假设⑤方向一致。三重限制如实声明：pembro 臂 TNBC 仅 29 例；CAF/ECM 在该臂方向反直觉；臂级未做随机校准。全段标注 hypothesis-generating、cross-regimen（I-SPY2 为平行随机、非头对头）。

## 整合完成
- Note 13 → N13.25（N13.24 I-SPY2）；主稿 3.8 ICB 段改写为含 I-SPY2 split；4.4 第 5 条限制更新；4.1 与结论段结构性重复消除（-17 词）；core = 5,983（硬线内余量 17）。
- xlsx：S34（stage17/18/19 五表）+ S35（I-SPY2）；补充材料清单同步。
- docx：manuscript_main_BIB_v6.docx。
- OSF 分析计划补 Stage 15–20。

## 遗留（唯一）
- GSE350571 Xenium（Stage 14）：P17 h5 69M/113.6M 下载中（NCBI 限速），transcripts 与其余 5 样本未开始。分析脚本已就绪（run_bulk_stage14_xenium.py）。
