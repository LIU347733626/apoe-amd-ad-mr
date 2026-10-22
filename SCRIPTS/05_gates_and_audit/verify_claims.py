#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_claims.py — 手稿数字溯源门禁（监督机制 ①，2026-09-10）

规则：手稿里出现的**每一个**数字，必须满足下列之一，否则报 UNMATCHED：
  (1) 能在本地结果文件（CSV/TSV/JSON）里找到容差内的匹配  → MATCHED
  (2) 出现在"外部来源"语境里（同行含 [引用编号] / PMID / NCT / doi / 队列名+年）→ EXTERNAL
  (3) 在显式白名单里（定义性常数：P 阈值 5e-8、染色体坐标、版本号、作者年等）→ WHITELISTED

用法: python3 scripts/verify_claims.py [手稿路径]
输出: docs/audit_20260909/CLAIM_PROVENANCE_<date>.md + 控制台摘要
不通过的判定标准：UNMATCHED 数 > 0 即 exit 1（门禁）。
"""
import io, os, re, sys, json, csv, glob, math, datetime

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
MS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v64_20260910.md")
OUTDIR = os.path.join(PROJ, "docs", "audit_20260909")
WHITELIST = os.path.join(PROJ, "scripts", "claim_whitelist.json")

# v2（2026-09-10）：语料改为**手稿实际引用的结果文件清单**（scripts/claim_sources.txt），
# 否则 4800 万 token 的宽语料会让任何数字都"匹配上"——通过一切的门禁等于没有门禁。
CORPUS_MANIFEST = os.path.join(PROJ, "scripts", "claim_sources.txt")
MAX_MB = 60

# ---------------------------------------------------------------- 数值提取
SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-", "⁺": "+"}
def norm_sup(s):
    return "".join(SUP.get(ch, ch) for ch in s)

NUM_RE = re.compile(
    r"(?<![A-Za-z0-9_.])"
    r"(-?\d{1,3}(?:,\d{3})+(?:\.\d+)?"
    r"|-?\d+\.\d+"
    r"|-?\d+"
    r")(?:\s*[×x]\s*10\s*([⁻⁺\-+]?\d+))?"
    r"\s*(%|kcal/mol)?")

def to_float(tok, exp, pct):
    v = float(tok.replace(",", ""))
    if exp:
        v = v * (10 ** int(norm_sup(exp)))
    if pct:
        v = v / 100.0                      # 百分比统一存为比例，匹配时双向尝试
    return v

def numbers_in(text, keep_context=True):
    """返回 [(value, 原字符串, 上下文)]"""
    out = []
    for m in NUM_RE.finditer(text):
        raw = m.group(0)
        tok, exp, pct = m.group(1), m.group(2), m.group(3)
        try:
            v = to_float(tok, exp, pct)
        except Exception:
            continue
        ctx = text[max(0, m.start() - 160): m.end() + 60].replace("\n", " ") if keep_context else ""
        out.append((v, raw.strip(), ctx, bool(pct), bool(exp)))
    return out

# ---------------------------------------------------------------- 语料索引
def corpus_numbers():
    idx = []          # (value, "file:col")
    files = []
    if os.path.exists(CORPUS_MANIFEST):
        for line in io.open(CORPUS_MANIFEST, encoding="utf-8"):
            line = line.strip()
            if line and os.path.isfile(os.path.join(PROJ, line)):
                files.append(os.path.join(PROJ, line))
    for p in files:
            if os.path.getsize(p) > MAX_MB * 1e6:
                continue
            ext = os.path.splitext(p)[1].lower()
            rel = os.path.relpath(p, PROJ)
            try:
                if ext in (".csv", ".tsv", ".txt"):
                    with io.open(p, encoding="utf-8", errors="ignore") as f:
                        head = f.readline()
                        if not head:
                            continue
                        delim = "\t" if ext == ".tsv" or head.count("\t") > head.count(",") else ","
                        cols = next(csv.reader([head], delimiter=delim))
                        for r in csv.DictReader(f, fieldnames=cols, delimiter=delim):
                            for c, val in r.items():
                                if val is None:
                                    continue
                                for mm in re.finditer(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", val):
                                    try:
                                        idx.append((float(mm.group(0)), f"{rel}:{c}"))
                                    except Exception:
                                        pass
                elif ext == ".json":
                    txt = io.open(p, encoding="utf-8", errors="ignore").read()
                    for mm in re.finditer(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", txt):
                        try:
                            idx.append((float(mm.group(0)), f"{rel}"))
                        except Exception:
                            pass
                elif ext in (".md", ".log"):
                    txt = io.open(p, encoding="utf-8", errors="ignore").read()
                    for mm in re.finditer(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", txt):
                        try:
                            idx.append((float(mm.group(0)), f"{rel}"))
                        except Exception:
                            pass
            except Exception:
                continue
    return idx

def match(v, idx, pct, is_count=False):
    """容差：相对 0.6%（兼容四舍五入到 1–3 位小数）或绝对 1e-9（极小值）。

    FIX(2026-09-10, M-1): 计数类数字（无小数点的整数 token，如样本量/病例数/SNP 数）
    是精确事实，不能用 0.6% 相对容差——该容差在 472,750 上宽达 ±2,836、在 13,947 上
    达 ±84，会让"写成 13,900"也判 MATCHED。计数改用 ±1 的精确匹配。
    """
    cands = [v, v * 100.0] if pct else [v]
    hits = []
    for cv in cands:
        tol = 1.0 if is_count else max(abs(cv) * 0.006, 1e-9)
        for c, src in idx:
            if abs(c - cv) <= tol:
                hits.append(src)
                if len(hits) >= 4:
                    return hits
        if hits:
            return hits
    return hits

# ---------------------------------------------------------------- 主流程
# FIX(2026-09-10): 原正则只认 FinnGen/MVP/PMID/doi，导致"引用外部数据集/文献的计数"
# （IAMDGC N、deCODE N、FAERS 报告数、ADGC 分层 n、Bellenguez 论文计数等）在语料收窄后
# 被误报为 UNMATCHED；此前它们是靠 0.6% 宽相对容差匹配到语料里的【无关数字】而假通过。
EXTERNAL_MARK = re.compile(
    r"\[\d+(?:,\d+)*\]|PMID|NCT\d|doi|clinicaltrials|GSE\d+|GCST\d*|dbGaP|"
    r"FinnGen|MVP|IAMDGC|deCODE|Olink|SomaScan|UK Biobank|UKB|FAERS|CVARD|ADGC|ADSP|"
    r"Jansen|Bellenguez|Kunkle|Fritsche|Rajabli|Strunz|EyeGEx|GTEx|HapMap|Million Veteran|"
    r"Alzheimers Dement|JAMA Ophthalmol|Brain|Nat Genet|Lancet")

def main():
    ms = io.open(MS, encoding="utf-8").read()
    body = ms[:ms.index("## Figure Legends")] if "## Figure Legends" in ms else ms
    body = "\n".join(l for l in body.split("\n") if not l.startswith("> **Version note"))
    # 参考文献表与方法里的"外部数据来源"整段豁免（由 EXTERNAL_MARK 判定）
    idx = corpus_numbers()
    print(f"[index] corpus numeric tokens: {len(idx)}")

    wl = {}
    if os.path.exists(WHITELIST):
        wl = json.load(io.open(WHITELIST, encoding="utf-8")).get("whitelist", {})

    rows, n_ok, n_ext, n_wl, n_bad, n_triv = [], 0, 0, 0, 0, 0
    for v, raw, ctx, pct, exp in numbers_in(body):
        key = raw
        if (not pct) and (not exp) and ("." not in raw.replace(",", "")) and abs(v) <= 10:
            rows.append((raw, v, "TRIVIAL", "个位数整数，非主张（不计入门禁）")); n_triv += 1; continue
        if key in wl:
            rows.append((raw, v, "WHITELISTED", wl[key])); n_wl += 1; continue
        is_count = (not pct) and (not exp) and ("." not in raw.replace(",", ""))
        hits = match(v, idx, pct, is_count=is_count)
        if hits:
            rows.append((raw, v, "MATCHED", "; ".join(sorted(set(hits))[:3]))); n_ok += 1; continue
        if EXTERNAL_MARK.search(ctx):
            rows.append((raw, v, "EXTERNAL", ctx[:150])); n_ext += 1; continue
        rows.append((raw, v, "UNMATCHED", ctx[:200])); n_bad += 1

    date = datetime.date.today().isoformat()
    out = os.path.join(OUTDIR, f"CLAIM_PROVENANCE_{date.replace('-','')}.md")
    with io.open(out, "w", encoding="utf-8") as f:
        f.write(f"# 手稿数字溯源门禁（{date}）\n\n")
        f.write(f"手稿: `{os.path.relpath(MS, PROJ)}`\n\n")
        f.write(f"| 判定 | 数量 |\n|---|---|\n| MATCHED（在引用的结果文件中找到） | {n_ok} |\n"
                f"| EXTERNAL（外部来源语境） | {n_ext} |\n| WHITELISTED（定义性常数） | {n_wl} |\n"
                f"| TRIVIAL（个位数整数，非主张） | {n_triv} |\n| **UNMATCHED（未溯源）** | **{n_bad}** |\n\n")
        f.write("## UNMATCHED 明细（需人工处理）\n\n")
        for raw, v, st, note in rows:
            if st == "UNMATCHED":
                f.write(f"- `{raw}`（值 {v:g}）… {note}\n")
        f.write("\n## WHITELISTED 明细\n\n")
        for raw, v, st, note in rows:
            if st == "WHITELISTED":
                f.write(f"- `{raw}` — {note}\n")
    print(f"MATCHED={n_ok} EXTERNAL={n_ext} WHITELISTED={n_wl} TRIVIAL={n_triv} UNMATCHED={n_bad}")
    print("report:", out)
    return 1 if n_bad else 0

if __name__ == "__main__":
    sys.exit(main())
