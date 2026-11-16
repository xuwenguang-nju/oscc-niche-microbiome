#!/usr/bin/env python3
"""抓取 ORMCKB（口腔微生物组-疾病知识库）并提取与 OSCC 相关的部分。

来源：http://sysbio.org.cn/ormckb/（West China Hospital / Sichuan University）
接口：无官方批量下载；页面为 PHP + DataTables（导出按钮是客户端 Excel/CSV/PDF）。
      本脚本因此直接抓浏览页表格与逐类群详情页。

原始 HTML 全部保存到 data/raw/ormckb/（只读留档），tidy 结果写入 results/tables/。

用法：python analysis/scripts/fetch_ormckb.py
"""
from __future__ import annotations

import re
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

BASE = "http://sysbio.org.cn/ormckb/"
RAW = Path("data/raw/ormckb")
TAB = Path("results/tables")
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

# 浏览页的表格端点（取自 Browse.html 的链接）
BROWSE_TABLES = {
    "phylum": "php/Phylum.php",
    "class": "php/Class.php",
    "order": "php/Order.php",
    "family": "php/Family.php",
    "genus": "php/Genus.php",
    "species": "php/Species.php",
    "disease_oral": "php/oreal.php",
    "disease_systemic": "php/Systemic.php",
    "disease_other_1": "php/other2.php",
    "disease_other_2": "php/other3.php",
    "disease_other_3": "php/other4.php",
    "intervention_treatment": "php/1treatment.php",
    "intervention_drug": "php/2drug.php",
    "intervention_food": "php/3FD.php",
    "intervention_smoking": "php/4smoking.php",
    "intervention_age": "php/5age.php",
    "bacteriostat_extracts": "php/extracts.php",
    "bacteriostat_compounds": "php/compounds.php",
}

OSCC_PATTERN = (
    r"Oral Squamous|Oral Cancer|Tongue Cancer|Tongue Squamous|Head And Neck Squamous"
    r"|Malignant Head and Neck"
)


def fetch(path: str, data: dict | None = None, tries: int = 3) -> str:
    for attempt in range(tries):
        try:
            if data is None:
                req = urllib.request.Request(BASE + path, headers=UA)
            else:
                req = urllib.request.Request(
                    BASE + path,
                    data=urllib.parse.urlencode(data).encode(),
                    headers={**UA, "Content-Type": "application/x-www-form-urlencoded",
                             "Referer": BASE + "Browse.html"},
                )
            with urllib.request.urlopen(req, timeout=120) as rh:
                return rh.read().decode("utf-8", errors="replace")
        except Exception:
            if attempt == tries - 1:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


def strip_html(cell: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", cell)).strip()


def parse_table(html: str) -> pd.DataFrame:
    m = re.search(r"<table[^>]*>(.*?)</table>", html, flags=re.S | re.I)
    if not m:
        return pd.DataFrame()
    body = m.group(1)
    heads = [strip_html(c) for c in re.findall(r"<th[^>]*>(.*?)</th>", body, flags=re.S | re.I)]
    rows = []
    for tr in re.findall(r"<tr>(.*?)</tr>", body, flags=re.S | re.I):
        tds = re.findall(r"<td[^>]*>(.*?)</td>", tr, flags=re.S | re.I)
        if not tds:
            continue
        cells = [strip_html(c) for c in tds]
        if len(cells) != len(heads):
            continue
        rec = dict(zip(heads, cells))
        rid = re.search(r"id=['\"]([^'\"]+)['\"]", tr)
        rec["entry_id"] = rid.group(1) if rid else None
        rows.append(rec)
    return pd.DataFrame(rows)


def parse_detail_tables(html: str) -> list[tuple[list[str], list[list[str]], str]]:
    """返回 (表头, 仅 td 的行, 该表原始 HTML)，第三项用于解析 th/td 键值型样本表。"""
    out = []
    for t in re.findall(r"<table[^>]*>(.*?)</table>", html, flags=re.S | re.I):
        heads = [strip_html(c) for c in re.findall(r"<th[^>]*>(.*?)</th>", t, flags=re.S | re.I)]
        rows = []
        for tr in re.findall(r"<tr>(.*?)</tr>", t, flags=re.S | re.I):
            tds = re.findall(r"<td[^>]*>(.*?)</td>", tr, flags=re.S | re.I)
            if tds:
                rows.append([strip_html(c) for c in tds])
        out.append((heads, rows, t))
    return out


def parse_th_td_pairs(table_html: str) -> list[tuple[str, str]]:
    """样本元数据表是 th(键)/td(值) 交替布局；按文档顺序把相邻的 th→td 抽成对。"""
    pairs: list[tuple[str, str]] = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", table_html, flags=re.S | re.I):
        cells = re.findall(r"<(th|td)[^>]*>(.*?)</\1>", tr, flags=re.S | re.I)
        pending_key = None
        for tag, content in cells:
            value = strip_html(content)
            if tag.lower() == "th":
                pending_key = value
            elif pending_key is not None:
                pairs.append((pending_key, value))
                pending_key = None
    return pairs


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    (RAW / "_retrieved_utc.txt").write_text(stamp + "\n" + BASE + "\n", encoding="utf-8")

    # 1) 浏览页全表
    frames = {}
    for name, path in BROWSE_TABLES.items():
        html = fetch(path)
        (RAW / f"{name}.html").write_text(html, encoding="utf-8")
        df = parse_table(html)
        df["table"] = name
        frames[name] = df
        print(f"{name:<18} rows={len(df):<6} -> data/raw/ormckb/{name}.html")

    all_rows = pd.concat([f for f in frames.values() if len(f)], ignore_index=True)
    all_rows.to_csv(TAB / "ormckb_taxonomy_associations.csv", index=False, encoding="utf-8-sig")
    print(f"合并 tidy 表 -> results/tables/ormckb_taxonomy_associations.csv ({len(all_rows)} 行)")

    # 2) OSCC/HNSCC 子集（属水平为主）
    genus = frames["genus"]
    oscc = genus[genus["Disease"].str.contains(OSCC_PATTERN, case=False, na=False)].copy()
    oscc["taxid"] = oscc["entry_id"].str.split("-").str[1]
    oscc.to_csv(TAB / "ormckb_oscc_genus_entries.csv", index=False, encoding="utf-8-sig")
    print(f"属水平 OSCC/HNSCC 条目 {len(oscc)} 条 / {oscc['Microbe'].nunique()} 类群 / "
          f"{oscc['PMID'].nunique()} 篇")

    # 3) 方向投票表
    vote = (oscc.groupby(["Microbe", "Alteration"])["PMID"].nunique().unstack(fill_value=0)
            .rename(columns={"Up": "n_papers_up", "Down": "n_papers_down"}))
    for c in ("n_papers_up", "n_papers_down"):
        if c not in vote:
            vote[c] = 0
    vote["n_papers"] = vote["n_papers_up"] + vote["n_papers_down"]
    vote["direction"] = [
        "conflict" if (u > 0 and d > 0) else ("Up" if u > 0 else "Down")
        for u, d in zip(vote["n_papers_up"], vote["n_papers_down"])
    ]
    vote = vote.sort_values("n_papers", ascending=False)
    vote.to_csv(TAB / "ormckb_oscc_genus_votes.csv", encoding="utf-8-sig")
    print(f"方向投票表 -> results/tables/ormckb_oscc_genus_votes.csv "
          f"(conflict={(vote['direction']=='conflict').sum()})")

    # 4) 逐类群详情页：拿 OSCC 行的 p 值 + 测序方法
    detail_rows = []
    taxids = sorted(t for t in oscc["taxid"].dropna().unique())
    for tid in taxids:
        entry = oscc[oscc["taxid"] == tid].iloc[0]["entry_id"]
        html = fetch(f"php/info.php?id={entry}")
        (RAW / f"detail_{entry}.html").write_text(html, encoding="utf-8")
        tables = parse_detail_tables(html)
        for heads, rows, _raw in tables:
            if "TaID" in heads and "Pvalue" in heads:
                for r in rows:
                    rec = dict(zip(heads, r))
                    if re.search(OSCC_PATTERN, rec.get("Dame", "") or "", re.I):
                        detail_rows.append(rec)
        time.sleep(0.4)
    det = pd.DataFrame(detail_rows).drop_duplicates()
    if len(det):
        det.to_csv(TAB / "ormckb_oscc_entry_pvalues.csv", index=False, encoding="utf-8-sig")
        n_p = det["Pvalue"].astype(str).str.contains(r"\d", na=False).sum()
        print(f"OSCC 条目 p 值表 {len(det)} 行，其中报告了数值 p 的 {n_p} 行")
    else:
        print("明细表为空（可能表格结构变动）")

    # 5) 每篇 OSCC 论文的样本类型（详情页的 Sample type 字段）
    pub_rows = []
    for pmid, grp in oscc.groupby("PMID"):
        entry = grp.iloc[0]["entry_id"]
        html = fetch(f"php/info.php?id={entry}")
        tables = parse_detail_tables(html)
        pub = {}
        for heads, _rows, raw in tables:
            if "PMID" in heads[:2] and "Title" in heads:
                pub = dict(parse_th_td_pairs(raw))  # 与样本表同为 th/td 键值布局
                break
        sample_types, sources, descriptions = [], [], []
        for _heads, _rows, raw in tables:
            for key, value in parse_th_td_pairs(raw):
                if key == "Sample type" and value and value != "Null":
                    sample_types.append(value)
                elif key == "Sample source" and value and value != "Null":
                    sources.append(value)
                elif key == "Description" and value and value != "Null":
                    descriptions.append(value)
        pub_rows.append({
            "pmid": pmid,
            "title": pub.get("Title", ""),
            "journal": pub.get("Journal", ""),
            "year": pub.get("Year", ""),
            "article_type": pub.get("Article type", ""),
            "sample_types": ";".join(sorted(set(sample_types))),
            "sample_sources": ";".join(sorted(set(sources))),
            "sample_groups": ";".join(sorted(set(descriptions)))[:300],
            "n_genera_reported": grp["Microbe"].nunique(),
        })
        time.sleep(0.4)
    pubs = pd.DataFrame(pub_rows).sort_values("year")
    pubs.to_csv(TAB / "ormckb_oscc_publications.csv", index=False, encoding="utf-8-sig")
    print(f"论文表 -> results/tables/ormckb_oscc_publications.csv ({len(pubs)} 篇)")
    print(pubs["sample_types"].value_counts().to_string())

    # 6) 检索工具返回的 OSCC 条目（含"干预"轴，浏览表不含干预信息）
    html = fetch("php/search.php", {"micro": "", "disease": "Oral Squamous Cell Carcinoma",
                                    "change": "", "submit": ""})
    (RAW / "search_oscc.html").write_text(html, encoding="utf-8")
    sdf = parse_table(html)
    sdf.to_csv(TAB / "ormckb_oscc_search_entries.csv", index=False, encoding="utf-8-sig")
    print(f"检索接口 OSCC 条目 {len(sdf)} 行 -> results/tables/ormckb_oscc_search_entries.csv")

    # 7) 29 个 OSCC 候选属在 ORMCKB 中的泛病种覆盖（说明"属级关联并不 OSCC 特异"）
    cand = ["Fusobacterium", "Porphyromonas", "Prevotella", "Streptococcus", "Veillonella",
            "Parvimonas", "Peptostreptococcus", "Capnocytophaga", "Treponema", "Actinomyces",
            "Rothia", "Neisseria", "Haemophilus", "Gemella", "Eikenella", "Aggregatibacter",
            "Tannerella", "AlloPrevotella", "Selenomonas", "Leptotrichia", "Solobacterium",
            "Filifactor", "Olsenella", "Dialister", "Megasphaera", "Campylobacter",
            "Stomatobaculum", "Lactobacillus", "Bifidobacterium"]
    genus_df = genus.copy()
    genus_df["Microbe"] = genus_df["Microbe"].str.replace(r"\s+", " ", regex=True).str.strip()
    rows = []
    for g in cand:
        sub = genus_df[genus_df["Microbe"].str.lower() == g.lower()]
        rows.append({
            "genus": g,
            "n_diseases_genus_level": sub["Disease"].nunique() if len(sub) else 0,
            "n_oscc_rows_in_genus_table": int(sub["Disease"].str.contains(OSCC_PATTERN, case=False, na=False).sum()),
        })
    cov = pd.DataFrame(rows).sort_values("n_diseases_genus_level", ascending=False)
    cov.to_csv(TAB / "ormckb_candidate_genus_coverage.csv", index=False, encoding="utf-8-sig")
    print(f"候选属覆盖 -> results/tables/ormckb_candidate_genus_coverage.csv "
          f"(疾病数 {cov['n_diseases_genus_level'].min()}–{cov['n_diseases_genus_level'].max()})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
