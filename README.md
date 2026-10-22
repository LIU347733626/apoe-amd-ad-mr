# FINAL_RELEASE_20260910 — 终稿复现包

**建立**：2026-09-10 · **对应手稿**：`MANUSCRIPT/00_manuscript_NC_v67_20260910.md`（NC 工作稿 v67）
**目的**：把所有**进入终稿**的干净代码集中一处，配数据清单与环境指纹，使关键结果可被第三方重跑验证。

> 本包是**终稿审核版**。历史上被取代/已撤回/仅调试用的脚本**不在包内**（它们仍在项目 `revision/` 下，但不属于终稿证据链）。

## 目录结构

```
FINAL_RELEASE_20260910/
├── README.md                 ← 本文件（复现指南）
├── DATA_MANIFEST.md          ← 数据来源、样本量、结果文件与生成脚本
├── ENV.md                    ← 环境指纹（R 与包版本）
├── SCRIPTS/                  ← 66 个干净脚本，按分析层分组 + MANIFEST.tsv（含 sha1）
│   ├── 01_primary_mr/        ← 主分析链（含黄金锚点生产脚本）
│   ├── 02_coloc_finemap/     ← 共定位与精细定位
│   ├── 03_mechanism/         ← 机制层（单细胞/受体轴/CellChat/扰动）
│   ├── 04_figures/           ← 出图脚本
│   ├── 05_gates_and_audit/   ← 门禁与审计（可复用的验证工具）
│   └── 06_manuscript_patches/← 手稿版本的手术脚本（每步可复现）
├── MANUSCRIPT/               ← 手稿快照（NC v67 / NA v6 / 投稿构建版）
└── FIGURES/                  ← 14 张定稿图（pdf+png，均 ≤180 mm）
```

## 复现步骤

**前置**：原始数据不随包分发（体积大且有唯一真源），路径见 `DATA_MANIFEST.md`。环境见 `ENV.md`（R 4.5.3 + gsmr **1.0.6**）。

```bash
cd /home/huan/projects/08_APOE_AMD_AD_MR
bash check_env.sh                     # 1) 环境门禁：缺包/缺数据当场拦截
bash revision/anchor_check.sh         # 2) 黄金锚点：跑对没跑对，先看这一条
```

**黄金锚点（必须逐位复现，±1e-6）**：
```
APOE→AMD  IVW β = +0.1818 , p = 7.89e-11     (FinnGen R13, 11 clumped UKB-PPP cis-pQTL)
APOE→AD   IVW β = −0.9418 , p = 4.08e-18
```

**分步复现**（脚本一律走 `bash run_r.sh < 脚本.R`，它自动设线程并跑 figgate）：

```bash
# 主分析
bash run_r.sh < SCRIPTS/01_primary_mr/20_r13.R
bash run_r.sh < SCRIPTS/01_primary_mr/24_replication_iamdgc.R      # IAMDGC 第一验证
bash run_r.sh < SCRIPTS/01_primary_mr/74_mvp_ancestry_mr.R          # MVP 第二验证
bash run_r.sh < SCRIPTS/01_primary_mr/98_second_validation_gsmr.R   # 第二验证层 GSMR
bash run_r.sh < SCRIPTS/01_primary_mr/96_leaveoneout_official11.R   # 工具杠杆检验
# 出图
bash run_r.sh < SCRIPTS/04_figures/Figure1_ggplot.R
# 门禁（手稿数字溯源 + 派生值重算 + 跨层一致性）
python3 SCRIPTS/05_gates_and_audit/verify_claims.py  <手稿路径>
python3 SCRIPTS/05_gates_and_audit/regenerate_claims.py <手稿路径>
python3 SCRIPTS/05_gates_and_audit/audit_crosslayer.py
```

## 本轮（2026-09-10）审计修复的缺陷（已全部落在包内脚本里）

| 编号 | 缺陷 | 影响 |
|---|---|---|
| A1 | `24_replication_iamdgc.R` 的 `grep` 输出无表头却传 `header=TRUE` → 吃掉首个工具 | IAMDGC 复制由 10 工具误算；修复后 11 工具：IVW β=0.2132, P=8.5e-8 |
| 横向同族 | `60c` / `59` / `60` / `fig3_prep` 同类 `header` 问题 | 各丢 1 个变异（PP 值不变） |
| A5 | 3 处 `mr_presso()` 缺 `data=`（无默认值）→ 报错被 tryCatch 吞掉 | MR-PRESSO **从未运行**；另修列名与 `[["Main MR results"]]` |
| A4 | `40a_*_worker.R` Egger 行 Beta/SE 取斜率、P 取**截距** | 面板 Egger 的 P 不同源 |
| A6 | `35_r13.R` 仍写已作废的 `rs7442`（正确为 `rs7412`） | 重跑会覆盖已修正的 ε2 note |
| A7 | `53b_clumped_mr.R` 用 `_amd_snp`（R 非法标识符）+ 未定义变量 | 脚本**语法级不可运行** |
| B2 | 7 处 `mr_presso()` 缺 `seed=` | 同数据多次结果不同 |
| 口径 | Second validation 两层原为 IVW，与其余行的 GSMR 不一致 | 已补跑 GSMR（`98_...R`），全图统一 |

**共定位层的一处撤销**：AD 侧曾报告 "multi-signal (suSiE) PP3 = 1.0"，经查是 `p > 0` 过滤器删除了 APOE 区 31 个 p 被钳制为 `0.0` 的最强变异（含 rs429358）造成的伪影，且产出脚本实为 `coloc.abf` 而非 `coloc.suSiE`。已整层撤销（`06_manuscript_patches/patch_v61.py`）。

## 已知限制（写作与引用时须留意）

1. **p = 0 是双精度下溢**（|z| ≳ 37.5），不是格式问题。任何 `p > 0` / `p != 0` 的"清洗"都在**精确删除最强变异**。正确做法：按有限 β 与正 SE 保留。
2. **`grep -f <rsid列表>` 的输出恒无表头**（表头行不匹配 rsID）——提取时一律 `header=FALSE` 并断言工具覆盖数。
3. **暴露文件（UKB-PPP APOE）有 49 个变异 p 字面为 0** → 自行 clumping 时索引 SNP 由并列打破规则决定，工具集**不可唯一确定**。终稿主分析用的是**官方给定**的 clumped 文件，故可复现；自建 clumping 的结果不可与之直接比较。
4. **`AD_r13.ma` 的 p 列在 APOE 区被钳制为字面 `0.0`**（31 个），而 `IAMDGC_late.ma` 的字面 `0` 是另一批（160 个）——判断数值 0 必须用数值比较或正则 `^-?0([.]0*)?$`，字符串 `=="0"` 会全中一个、全漏另一个。
5. **Bellenguez 2022 的 APOE 区受限**：缺 rs429358/rs7412，且 12 个经典 tag SNP 全部不存在；发布文件 per-variant N 上限 487,511 ≠ 论文的 788,989 → 仅作方向性验证。
6. **Jansen 2019 与 FinnGen/IAMDGC 不在同一效应尺度**（|β| 小 3.8–6.1 倍）→ 只作方向/显著性验证。
7. **AD 侧共定位在方法上不可分离**：ε4 既是最强的水平调节因子又是编码变异；多信号（suSiE）在本数据上退化（全部 credible set 塌缩为单变异）。

## 验证状态（建立本包时）

- `anchor_check.sh` 黄金锚点：**通过**（+0.1818 / −0.9418）
- `verify_claims`（手稿数字溯源，计数类已收紧为 ±1 精确匹配）：**UNMATCHED = 0**
- `regenerate_claims`（派生值重算逐字比对）：**28/28 consistent**
- `audit_crosslayer`（手稿↔结果↔图↔表）：**117/117 PASS**
- 14 张定稿图：**figgate 双裁判 PASS**，宽度 **全部 ≤180 mm**
- 全部入包脚本：`Rscript -e 'parse(...)'` 通过

---
*本包由 `scripts/build_final_release.py` 生成（脚本清单与 sha1 见 `SCRIPTS/MANIFEST.tsv`），可重跑重建。*
