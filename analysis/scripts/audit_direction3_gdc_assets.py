#!/usr/bin/env python3
"""审核方向三（宿主胚系变异 -> 瘤内菌丰度 -> 肿瘤分子表型）所需的 TCGA-HNSC 资产。

只查询 GDC API 的元数据（文件计数 / 文件体积 / 访问级别 / 样本类型），不下载任何受控数据。
输出 results/tables/direction3_asset_audit.csv，供设计报告引用与复核。

用法：
    python analysis/scripts/audit_direction3_gdc_assets.py
"""
from __future__ import annotations

import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://api.gdc.cancer.gov"
PROJECT = "TCGA-HNSC"
OUT = Path("results/tables/direction3_asset_audit.csv")


def api(path: str, params: dict) -> dict:
    url = f"{BASE}/{path}?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return json.loads(fh.read())
        except Exception as exc:  # 网络抖动重试
            if attempt == 3:
                raise
            print(f"  retry {attempt + 1}/3: {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


def in_filter(field: str, values: list[str]) -> dict:
    return {"op": "in", "content": {"field": field, "value": values}}


def hnsc(*conditions: dict) -> dict:
    return {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT]}},
            *conditions,
        ],
    }


def facet(filters: dict, facet_fields: str) -> dict[str, list[tuple[str, int]]]:
    body = api("files", {"size": 0, "facets": facet_fields, "filters": json.dumps(filters)})
    out = {}
    for name, agg in body["data"]["aggregations"].items():
        out[name] = sorted(((b["key"], b["doc_count"]) for b in agg["buckets"]), key=lambda x: -x[1])
    return out


def summarise(filters: dict) -> tuple[int, float]:
    body = api(
        "files",
        {
            "size": 3000,
            "fields": "file_id,file_size",
            "filters": json.dumps(filters),
        },
    )
    hits = body["data"]["hits"]
    return len(hits), sum(h.get("file_size") or 0 for h in hits) / 1024**3


def main() -> int:
    release = api("status", {})["data_release"]
    project = api(
        f"projects/{PROJECT}",
        {"expand": "summary,summary.experimental_strategies", "pretty": "false"},
    )["data"]
    n_cases = project["summary"]["case_count"]
    strategies = {
        s["experimental_strategy"]: s["case_count"] for s in project["summary"]["experimental_strategies"]
    }

    rows: list[dict] = []

    def add(asset: str, desc: str, filters: dict, link: str, role: str, note: str = "") -> None:
        n, gb = summarise(filters)
        acc = facet(filters, "access").get("access", [])
        access = ", ".join(f"{k}={v}" for k, v in acc) if acc else "n/a"
        st = facet(filters, "cases.samples.sample_type").get("cases.samples.sample_type", [])
        breakdown = "; ".join(f"{k}={v}" for k, v in st) if st else "n/a"
        rows.append(
            {
                "asset": asset,
                "description": desc,
                "n_files": n,
                "total_gb": round(gb, 1),
                "access": access,
                "case_sample_type_breakdown": breakdown,
                "direction3_link": link,
                "role_in_pipeline": role,
                "note": note,
            }
        )

    # 链路①：宿主胚系变异
    add(
        "germline_snp6_calls",
        "Genotyping Array 的 Simple Germline Variation（SNP6 Birdseed 基因型）",
        hnsc(in_filter("data_type", ["Simple Germline Variation"])),
        "1-germline",
        "mQTL 基因型（首选：全基因组标记，插补后覆盖调控区）",
        "全部受控访问",
    )
    add(
        "snp6_raw_intensities",
        "SNP6 原始荧光强度（CEL），含血源正常样本",
        hnsc(in_filter("data_type", ["Raw Intensities"])),
        "1-germline",
        "自行调用基因型（CEL 体积远小于 BAM，本地可承受）",
        "受控访问；配对的胚系来源为 blood derived normal",
    )
    add(
        "wxs_normals",
        "WXS 比对读段（Blood Derived Normal）",
        hnsc(
            in_filter("data_type", ["Aligned Reads"]),
            in_filter("experimental_strategy", ["WXS"]),
            in_filter("cases.samples.sample_type", ["Blood Derived Normal"]),
        ),
        "1-germline",
        "自行调用胚系变异（体积大，需云端处理）",
        "受控访问",
    )
    add(
        "rna_aligned",
        "RNA-seq 比对读段（BAM）",
        hnsc(in_filter("data_type", ["Aligned Reads"]), in_filter("experimental_strategy", ["RNA-Seq"])),
        "1-germline",
        "RNA-seq 推断基因型 / HLA 分型的预备路线",
        "受控访问，不构成开放替代方案",
    )

    # 链路②③：微生物丰度 + 宿主分子表型（开放）
    add(
        "rna_expression_open",
        "RNA-seq 基因表达定量",
        hnsc(in_filter("data_type", ["Gene Expression Quantification"]), in_filter("access", ["open"])),
        "2-3-abundance-to-phenotype",
        "宿主分子表型（表达 / 免疫浸润 / 通路）",
        "开放访问",
    )
    add(
        "methylation_open",
        "甲基化 beta 值",
        hnsc(in_filter("data_type", ["Methylation Beta Value"]), in_filter("access", ["open"])),
        "2-3-abundance-to-phenotype",
        "表观层分子表型",
        "开放访问",
    )
    add(
        "somatic_maf_open",
        "体细胞突变（Masked Somatic Mutation）",
        hnsc(in_filter("data_type", ["Masked Somatic Mutation"]), in_filter("access", ["open"])),
        "3-phenotype",
        "突变特征与驱动基因状态",
        "开放访问；Annotated Somatic Mutation 为受控",
    )
    add(
        "clinical_open",
        "临床与随访",
        hnsc(in_filter("data_type", ["Clinical Supplement"]), in_filter("access", ["open"])),
        "all-covariates",
        "分期、部位、吸烟饮酒、HPV/p16、生存",
        "开放访问",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"GDC data release: {release}")
    print(f"{PROJECT}: {n_cases} cases; strategies={strategies}")
    print(f"written -> {OUT}")
    for r in rows:
        print(
            f"  {r['asset']:<22} n={r['n_files']:<5} {r['total_gb']:>8.1f} GB  "
            f"{r['access']:<18} {r['case_sample_type_breakdown']}"
        )
    germline = [r for r in rows if r["direction3_link"].startswith("1-")]
    controlled = [r for r in germline if "controlled" in r["access"]]
    print(f"link-1 assets: {len(germline)}, all controlled: {len(controlled)} -> dbGaP authorisation required")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
