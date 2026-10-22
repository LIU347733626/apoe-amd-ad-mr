# 环境指纹（FINAL_RELEASE_20260910）

建立于 2026-09-10，实测输出：

```
R: R version 4.5.3 (2026-03-11) 
  data.table               1.18.4
  ggplot2                  4.0.3
  MendelianRandomization   0.10.0
  TwoSampleMR              0.7.9
  MRPRESSO                 1.0
  coloc                    5.2.3
  susieR                   0.14.2
  gsmr                     1.0.6
  GenomicSEM               0.0.5
  fgsea                    1.36.2
  patchwork                1.3.2
  Cairo                    1.7.0
gsmr API nargs: 14
```

## 版本锁定（重要）

- **gsmr 必须 1.0.6**：脚本按 1.0.6 的参数顺序写；1.1.1 会多一个位置参数导致全错。
- **R 4.5.3**（biomni_r conda 环境）。
- 分析一律用 `bash run_r.sh < 脚本.R>`（自动设线程 + 跑 figgate 图质检）。
- Python 侧（单细胞/机制）：`bioinfo` 环境，numpy 锁 1.26.4（<2）。
