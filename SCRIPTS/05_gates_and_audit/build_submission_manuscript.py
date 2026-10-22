#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_submission_manuscript.py — 生成"投稿干净版"手稿（2026-09-09）

工作稿里保留 `> **Version note ...` 版本注记（内部审计用，含中文）；投稿版必须去掉它们。
本脚本：
  1) 读工作稿（默认 v23），删掉所有 Version note 行；
  2) 校验：无 CJK 字符、无 Version note 残留、参考文献编号连续、正文引用编号都在范围内、
     图/补充图编号都有引用、占位符清单打印（供人工确认）；
  3) 写 `final_analysis/00_manuscript_NC_SUBMISSION_20260909.md` 与
     `NC_package/manuscript/manuscript_SUBMISSION.md`。
用法: python3 scripts/build_submission_manuscript.py [工作稿路径]
"""
import io, os, re, sys, glob

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v64_20260910.md")
OUT1 = os.path.join(PROJ, "final_analysis", "00_manuscript_NC_SUBMISSION_20260909.md")
OUT2 = os.path.join(PROJ, "NC_package", "manuscript", "manuscript_SUBMISSION.md")

txt = io.open(SRC, encoding="utf-8").read()
lines = txt.split("\n")

# 1) 删版本注记（含紧随其后的空行）
keep = []
skip_next_blank = False
for ln in lines:
    if ln.startswith("> **Version note"):
        skip_next_blank = True
        continue
    if skip_next_blank and ln.strip() == "":
        skip_next_blank = False
        continue
    skip_next_blank = False
    keep.append(ln)
clean = "\n".join(keep)
# 清理可能出现的连续空行
clean = re.sub(r"\n{4,}", "\n\n\n", clean)

# 2) 校验
problems = []
cjk = re.findall(r"[\u4e00-\u9fff]+", clean)
if cjk:
    problems.append(f"CJK 残留 {len(cjk)} 处: {cjk[:5]}")
if "Version note" in clean:
    problems.append("Version note 残留")

refs = re.findall(r"^(\d+)\. ", clean, flags=re.M)
ref_nums = [int(x) for x in refs]
if ref_nums:
    if ref_nums != list(range(1, len(ref_nums) + 1)):
        problems.append(f"参考文献编号不连续: {ref_nums[:5]}...{ref_nums[-5:]} (共 {len(ref_nums)})")
    maxref = max(ref_nums)
    cited = {int(x) for x in re.findall(r"\[(\d+)(?:,\d+)*\]", clean)}
    cited |= {int(y) for x in re.findall(r"\[(\d+(?:,\d+)+)\]", clean) for y in x.split(",")}
    out_of_range = sorted(c for c in cited if c > maxref)
    if out_of_range:
        problems.append(f"正文引用编号超出文献表: {out_of_range}")
    orphans = sorted(set(range(1, maxref + 1)) - cited)
    if orphans:
        problems.append(f"文献表里未被正文引用的条目: {orphans}")

# A20 (2026-09-09): 主图 1-6（Figure 6 于 v50 入稿）+ 补充图 1-7
# （Supplementary Fig. 8 / ScPagwas 已于 v51 撤回，不再要求引用）
for n in [f"Figure {i}" for i in range(1, 7)] + [f"Supplementary Fig. {i}" for i in range(1, 9)]:
    if n not in clean:
        problems.append(f"未引用 {n}")

# A20/P9：正文引用的 Supplementary Table Sxx 必须在随包 supplementary/ 下有文件
supp_files = {os.path.basename(p) for p in glob.glob(os.path.join(PROJ, "NC_package/supplementary/*"))}
cited_tables = sorted(set(re.findall(r"Supplementary Table (S\d+[a-d]?)", clean)))
missing_tables = [t for t in cited_tables if not any(x.startswith("Table" + t) for x in supp_files)]
if missing_tables:
    problems.append(f"引用的补充表在包内缺文件: {missing_tables}")

placeholders = [l.strip()[:160] for l in clean.split("\n") if re.search(r"\[(To be|TBD|待)", l)]

# 3) 写出
io.open(OUT1, "w", encoding="utf-8").write(clean)
os.makedirs(os.path.dirname(OUT2), exist_ok=True)
io.open(OUT2, "w", encoding="utf-8").write(clean)

print(f"source: {SRC}")
print(f"written: {OUT1}")
print(f"written: {OUT2}")
print(f"lines: {len(clean.split(chr(10)))} | chars: {len(clean)}")
print(f"references: {len(ref_nums)}")
if problems:
    print("\nPROBLEMS:")
    for p in problems:
        print("  -", p)
else:
    print("\nALL CHECKS PASS")
if placeholders:
    print("\nPLACEHOLDERS (需人工填):")
    for p in placeholders:
        print("  -", p)
