#!/usr/bin/env python3
"""Detail probe for the one candidate that fits S6's requirement.

S6 needs host expression and microbial signal from the SAME samples. Standard 16S
cohorts (ours included) have no host transcriptome, and TCGA microbial matrices
(TCMA) were unreachable from this host, so the only structure found so far is a
total-RNA-seq study: GSE227919 (oral pre-malignant lesions and OSCC). Total RNA-seq
captures both poly-A host transcripts and prokaryotic reads, so one library yields
both sides of the coupling.

This script checks what is actually deposited: platform, sample count, group
labels, and whether the supplementary files are raw counts or only a matrix.

Run: python3 probe_gse227919.py
"""
from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
UA = {"User-Agent": "wisp-science/1.0 (S6 probe)"}


def get(url: str, timeout: int = 40) -> tuple[int, bytes]:
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, b""
    except Exception as e:                                    # noqa: BLE001
        return -1, str(e).encode()


def main() -> None:
    out: dict = {"gse": "GSE227919"}

    # 1) series metadata via GEO SOFT (text, small)
    code, body = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                     "?acc=GSE227919&targ=self&form=text&view=brief")
    if code == 200:
        txt = body.decode("utf-8", "replace")
        for key in ["!Series_title", "!Series_summary", "!Series_overall_design",
                    "!Series_platform_id", "!Series_sample_id"]:
            vals = [l.split("=", 1)[1].strip() for l in txt.splitlines()
                    if l.startswith(key + " ")]
            if key == "!Series_sample_id":
                out["n_samples"] = len(vals)
            else:
                out[key.replace("!Series_", "")] = vals[:1] or ""
    out["soft_code"] = code

    # 2) supplementary file listing
    code, body = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                     "?acc=GSE227919&targ=self&form=text&view=full")
    if code == 200:
        txt = body.decode("utf-8", "replace")
        sup = re.findall(r"!Series_supplementary_file\s*=\s*(\S+)", txt)
        out["supplementary_files"] = [s.rsplit("/", 1)[-1] for s in sup]
        # data processing hints
        out["data_processing"] = [l.split("=", 1)[1].strip()[:160]
                                 for l in txt.splitlines()
                                 if l.startswith("!Sample_data_processing")][:8]
        chars = [l.split("=", 1)[1].strip() for l in txt.splitlines()
                 if l.startswith("!Sample_characteristics_ch1")]
        out["characteristics_examples"] = sorted(set(chars))[:20]
        src = [l.split("=", 1)[1].strip() for l in txt.splitlines()
               if l.startswith("!Sample_source_name_ch1")]
        out["source_names"] = sorted(set(src))[:20]
        lib = [l.split("=", 1)[1].strip() for l in txt.splitlines()
               if l.startswith("!Sample_library_strategy")]
        out["library_strategy"] = sorted(set(lib))
    out["full_code"] = code

    dst = ROOT / "results" / "tables" / "s6_gse227919_probe.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")

    for k, v in out.items():
        if isinstance(v, list) and len(v) > 8:
            print(f"{k} ({len(v)}):")
            for x in v[:8]:
                print("   ", x)
        else:
            print(f"{k}: {v}")
    print(f"\nwrote {dst}")


if __name__ == "__main__":
    main()
