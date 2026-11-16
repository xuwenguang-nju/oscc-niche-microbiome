"""Render every page of the tables PDF to PNG so the layout can be inspected."""

import pathlib

import pymupdf

# parents: [0]=scripts [1]=15-tables [2]=modules [3]=analysis [4]=project root
ROOT = pathlib.Path(__file__).resolve().parents[4]
PDF = ROOT / "analysis" / "out" / "tables-review.pdf"
OUT = ROOT / "analysis" / "out" / "pages"

OUT.mkdir(parents=True, exist_ok=True)
for old in OUT.glob("*.png"):
    old.unlink()

doc = pymupdf.open(PDF)
print(f"pdf: {PDF.name}  pages={doc.page_count}")

for i, page in enumerate(doc, 1):
    pix = page.get_pixmap(dpi=110)
    p = OUT / f"page{i:02d}.png"
    pix.save(p)
    w_mm = page.rect.width / 72 * 25.4
    h_mm = page.rect.height / 72 * 25.4
    print(f"  page {i}: {pix.width}x{pix.height} px  ({w_mm:.0f} x {h_mm:.0f} mm)  {p.stat().st_size//1024} KB")

doc.close()
