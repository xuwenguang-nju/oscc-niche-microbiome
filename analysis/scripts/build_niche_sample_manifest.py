"""Build the analysis-ready niche x patient sample manifest for the OSCC
intratumoral-microbiome source-attribution study.

Every rule below was decoded from the submitter's own metadata and cross-checked
against their own isolate / isolation_source labels (see
analysis/out/niche_patient_matrix_v2.txt and niche_final_metadata.txt).
Nothing is inferred beyond what the metadata states; ambiguous assignments are
flagged with niche_confidence != "high".

Output: results/tables/oscc_niche_sample_manifest.csv
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "wisp-science/1.0 (manifest builder)"}

FIELDS = ("run_accession,sample_accession,sample_title,sample_alias,isolate,"
          "isolation_source,read_count,library_strategy,instrument_model,fastq_ftp")


def get(url: str) -> bytes:
    last = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=180) as r:
                return r.read()
        except Exception as e:  # transient EBI 5xx
            last = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"failed after retries: {last}")


def ena_runs(acc: str) -> list[dict]:
    u = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + acc +
         "&result=read_run&fields=" + FIELDS + "&format=json&limit=0")
    return json.loads(get(u))


def sra_source_material(acc: str) -> dict[str, str]:
    """run/sample accession -> source_material_id, from NCBI SRA XML."""
    q = urllib.parse.urlencode({"db": "sra", "term": f"{acc}[BioProject]",
                                "retmax": 400, "retmode": "json"})
    ids = json.loads(get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + q))[
        "esearchresult"]["idlist"]
    out: dict[str, str] = {}
    for i in range(0, len(ids), 60):
        d = urllib.parse.urlencode({"db": "sra", "id": ",".join(ids[i:i + 60]),
                                    "retmode": "xml"}).encode()
        req = urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                                     data=d, headers=UA)
        with urllib.request.urlopen(req, timeout=180) as r:
            root = ET.fromstring(r.read())
        for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
            run = pkg.find(".//RUN")
            s = pkg.find(".//SAMPLE")
            if run is None or s is None:
                continue
            smid = ""
            for x in s.findall(".//SAMPLE_ATTRIBUTE"):
                if (x.findtext("TAG") or "") == "source_material_id":
                    smid = x.findtext("VALUE") or ""
            out[run.get("accession") or ""] = smid
        time.sleep(0.5)
    return out


ROWS: list[dict] = []


def add(study, run, niche, patient, conf, note, r):
    ROWS.append({
        "study": study,
        "run": run,
        "sample": r.get("sample_accession", ""),
        "niche": niche,
        "patient_id": patient,
        "niche_confidence": conf,
        "note": note,
        "reads": r.get("read_count", ""),
        "platform": r.get("instrument_model", ""),
        "strategy": r.get("library_strategy", ""),
        "fastq_ftp": r.get("fastq_ftp", ""),
    })


# ---------------- 1. PRJNA866676: 6 niches, alias = <prefix><patient> ----------------
NICHE_866676 = {
    "C": ("tumor", "high"),
    "CR": ("paracancerous", "high"),
    "CS": ("cancer_surface", "high"),
    "N": ("contralateral_normal_mucosa", "high"),
    "S": ("saliva", "high"),
    "T": ("tongue_coat", "high"),
}
for r in ena_runs("PRJNA866676"):
    m = re.match(r"^([A-Z]{1,2})(\d{1,3})$", (r.get("sample_alias") or "").strip())
    if not m:
        continue
    pref, pid = m.group(1), m.group(2)
    niche, conf = NICHE_866676.get(pref, ("unassigned", "low"))
    add("PRJNA866676", r["run_accession"], niche, f"JN{pid}", conf,
        f"alias={r['sample_alias']};submitter_label={r.get('isolate','')}", r)
    time.sleep(0)

# ---------------- 2. PRJNA1415359: P<n>-{T,N,Pre,Post} + AP<n>-pre-s ----------------
smid = sra_source_material("PRJNA1415359")
NICHE_1415359 = {"T": "tumor", "N": "adjacent_normal_tissue",
                 "Pre": "preop_saliva", "Post": "postop_saliva"}
for r in ena_runs("PRJNA1415359"):
    run = r["run_accession"]
    sm = smid.get(run, "") or (r.get("sample_alias") or "")
    m = re.match(r"^P(\d{1,3})-(T|N|Pre|Post)$", sm)
    if m:
        add("PRJNA1415359", run, NICHE_1415359[m.group(2)], f"KR{m.group(1)}", "high",
            f"source_material_id={sm}", r)
        continue
    m = re.match(r"^AP(\d{1,3})-pre-s$", sm)
    if m:
        add("PRJNA1415359", run, "saliva", f"AP{m.group(1)}", "high",
            "AP cohort: 仅术前唾液,无组织", r)
        continue
    add("PRJNA1415359", run, "unassigned", "", "low", f"source_material_id={sm}", r)

# ---------------- 3. PRJNA412445: <patient><a-e>, b/c ambiguous ----------------
NICHE_412445 = {
    "a": ("subgingival_plaque", "high"),
    "b": ("lesion_side_mucosa", "medium"),
    "c": ("contralateral_normal_mucosa", "medium"),
    "d": ("saliva", "high"),
    "e": ("tumor", "high"),
}
for r in ena_runs("PRJNA412445"):
    m = re.match(r"^(\d{1,3})([a-e])$", (r.get("sample_alias") or "").strip())
    if not m:
        continue
    pid, letter = m.group(1), m.group(2)
    niche, conf = NICHE_412445[letter]
    note = f"alias={r['sample_alias']};submitter_label={r.get('isolation_source','')}"
    if letter in ("b", "c"):
        note += ";b/c 在元数据中同为 mucous,此处按原文描述(病灶侧 vs 对侧)推断,需核对原文"
    add("PRJNA412445", r["run_accession"], niche, f"SC{pid}", conf, note, r)

# ---------------- 4. PRJNA666746: title 'OSCC <n>' / 'Normal <n>' ----------------
for r in ena_runs("PRJNA666746"):
    m = re.search(r"(OSCC|Normal)\s*(\d{1,3})", r.get("sample_title") or "")
    if not m:
        continue
    niche = "tumor" if m.group(1) == "OSCC" else "matched_normal"
    add("PRJNA666746", r["run_accession"], niche, f"IN{m.group(2)}", "high",
        f"title={r['sample_title']};isolate={r.get('isolate','')}", r)

# ---------------- 5. PRJNA1415268: T<n>/N<n> (Seoul) ----------------
smid2 = sra_source_material("PRJNA1415268")
for r in ena_runs("PRJNA1415268"):
    sm = smid2.get(r["run_accession"], "")
    m = re.match(r"^([TN])(\d{1,3})$", sm)
    if not m:
        add("PRJNA1415268", r["run_accession"], "unassigned", "", "low",
            f"source_material_id={sm}", r)
        continue
    niche = "tumor" if m.group(1) == "T" else "adjacent_normal_tissue"
    add("PRJNA1415268", r["run_accession"], niche, f"KS{m.group(2)}", "high",
        f"source_material_id={sm}", r)

# ---------------- 6. PRJNA597251: CT<n>/PT<n> ----------------
for r in ena_runs("PRJNA597251"):
    m = re.match(r"^(CT|PT)(\d{1,3})$", (r.get("sample_alias") or "").strip())
    if not m:
        continue
    niche = "tumor" if m.group(1) == "CT" else "paracancerous"
    add("PRJNA597251", r["run_accession"], niche, f"QD{m.group(2)}", "high",
        f"alias={r['sample_alias']};label={r.get('isolation_source','')}", r)

# ---------------- write ----------------
COLS = ["study", "run", "sample", "niche", "patient_id", "niche_confidence",
        "note", "reads", "platform", "strategy", "fastq_ftp"]
out = TAB / "oscc_niche_sample_manifest.csv"
with out.open("w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS)
    w.writeheader()
    w.writerows(ROWS)

print(f"manifest rows: {len(ROWS)} -> {out}")
print("\nby study and niche:")
tab: Counter = Counter((r["study"], r["niche"]) for r in ROWS)
for (s, n), c in sorted(tab.items()):
    print(f"  {s:<14} {n:<32} {c:>4}")
print("\npatients per study with a tumour sample:")
for s in sorted({r["study"] for r in ROWS}):
    pts = {r["patient_id"] for r in ROWS if r["study"] == s and r["niche"] == "tumor"}
    print(f"  {s:<14} tumor patients = {len(pts)}")
