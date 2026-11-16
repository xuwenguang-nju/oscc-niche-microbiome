"""Download the DADA2 taxonomy reference sets.

Two references are fetched, because they answer different questions:

  SILVA 138.2 (toGenus)   the cross-cohort default. Needed for a consistent
                          genus-level table across all six cohorts.
  eHOMD V15.22            the Human Oral Microbiome Database, formatted for
                          DADA2. Oral-specific, and the reference the original
                          PRJNA412445 study itself used, so it is the right
                          check for whether SILVA loses oral taxa that eHOMD
                          resolves.

Both are large (~139 MB and ~139 KB gzipped). Downloads are idempotent.

Outputs
-------
data/reference/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz
data/reference/ehomd/eHOMD_RefSeq_dada2_V15.22.fasta.gz
data/reference/REFERENCE_SOURCES.md   provenance record (URL, size, sha256)

Usage
-----
python analysis/modules/01-preprocess/scripts/fetch_references.py
"""

from __future__ import annotations

import hashlib
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
REF = ROOT / "data" / "reference"
UA = {"User-Agent": "wisp-science/1.0 (reference fetch)"}

TARGETS = [
    {
        "name": "SILVA 138.2 toGenus (DADA2 formatted)",
        "dest": REF / "silva" / "silva_nr99_v138.2_toGenus_trainset.fa.gz",
        "url": "https://zenodo.org/records/20955974/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1",
        "record": "https://zenodo.org/records/20955974",
        "licence": "SILVA (see upstream record)",
        "note": "used by assignTaxonomy(); genus level is the analysis unit",
    },
    {
        "name": "eHOMD RefSeq V15.22 (DADA2 formatted)",
        "dest": REF / "ehomd" / "eHOMD_RefSeq_dada2_V15.22.fasta.gz",
        "url": "https://zenodo.org/records/7380292/files/eHOMD_RefSeq_dada2_V15.22.fasta.gz?download=1",
        "record": "https://zenodo.org/records/7380292",
        "licence": "see upstream record",
        "note": "oral-specific reference; cross-check for SILVA genus assignments",
    },
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> str:
    if dest.exists() and dest.stat().st_size > 0:
        return "present"
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=1800) as r, \
            tmp.open("wb") as fh:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    tmp.replace(dest)
    return "downloaded"


def main() -> None:
    rows = []
    for t in TARGETS:
        status = fetch(t["url"], t["dest"])
        size = t["dest"].stat().st_size
        digest = sha256(t["dest"])
        rows.append({**t, "status": status, "bytes": size, "sha256": digest})
        print(f"{status:11} {t['dest'].relative_to(ROOT)}  {size / 1e6:.1f} MB")

    md = [
        "# Reference databases",
        "",
        "Fetched by `analysis/modules/01-preprocess/scripts/fetch_references.py`.",
        "Checksums are of the downloaded files as used in this project.",
        "",
    ]
    for r in rows:
        md += [
            f"## {r['name']}",
            f"- file: `{r['dest'].relative_to(ROOT)}`",
            f"- source: {r['record']}",
            f"- download URL: {r['url']}",
            f"- size: {r['bytes']} bytes",
            f"- sha256: `{r['sha256']}`",
            f"- licence: {r['licence']}",
            f"- role: {r['note']}",
            "",
        ]
    (REF / "REFERENCE_SOURCES.md").write_text("\n".join(md), encoding="utf-8")
    print("\nwrote", (REF / "REFERENCE_SOURCES.md").relative_to(ROOT))


if __name__ == "__main__":
    main()
