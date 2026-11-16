"""Render an existing plotting script at an exact physical width.

Why this exists
---------------
Every plotting script in this project sets a nominal `figsize` and then saves with
`bbox_inches="tight"`. `tight` recomputes the canvas from the content, so the delivered PNG
routinely ends up 15-90 mm wider than the nominal 180 mm figure width -- which is why 13 of 16
figures were over the journal limit even though "figsize=(7.09, ...)" was set everywhere.

Rewriting eight scripts to thread a scale factor through every `add_gridspec` call is fragile.
Instead this helper re-runs a script with `plt.figure`/`plt.subplots` monkey-patched so that the
requested figsize is multiplied by `scale`, captures the figure that the script would have saved,
measures its true tight bounding box, and iterates until the delivered width matches the target.

Font sizes are NOT scaled: point sizes stay as authored, so shrinking the canvas makes text
relatively larger, which is the desired behaviour at column width. That also means the width is
not a linear function of `scale`, hence the iteration.

Usage
-----
    from figlib import render_at_width
    render_at_width(script, out_png, mm=180)
"""
from __future__ import annotations

import runpy
from pathlib import Path


def project_root(start: Path | str | None = None) -> Path:
    """Walk up until a directory contains both `analysis/` and `results/`.

    Scripts sit at analysis/modules/<module>/scripts/<name>.py, i.e. four levels below the
    project root, but `parents[N]` silently resolves to analysis/ if a script moves. Searching
    for a marker is stable against that, and the fallback keeps a visible failure.
    """
    here = Path(start).resolve() if start else Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "analysis").is_dir() and (cand / "results").is_dir():
            return cand
    return here.parents[4]

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.figure import Figure                # noqa: E402


def _run_capturing(script: Path, scale: float):
    """Execute `script` with figure sizes scaled; return every Figure it tried to save.

    Several scripts emit more than one figure (a main panel and its sibling). Capturing a list
    keeps both so a multi-figure script can be laid out figure by figure.
    """
    held: dict[str, list] = {"figs": []}
    orig_figure = plt.figure
    orig_savefig = Figure.savefig

    def scaled_figure(*args, **kwargs):
        fs = kwargs.get("figsize")
        if fs:
            kwargs["figsize"] = (fs[0] * scale, fs[1] * scale)
        return orig_figure(*args, **kwargs)

    def capture_save(self, fname, *args, **kwargs):
        # Titles must be centred HERE: scripts often call panel_crops() after savefig,
        # which strips axes titles, so a post-script hook would be a silent no-op.
        centre_panel_titles(self)
        # Let the script's own savefig write normally: several scripts call
        # out.stat() right afterwards to print the file size, and swallowing the write
        # makes them crash. These are intermediates; the caller cleans them up.
        try:
            orig_savefig(self, fname, *args, **kwargs)
        except Exception:                                        # noqa: BLE001
            pass
        held["figs"].append(self)

    plt.figure = scaled_figure           # plt.subplots routes through plt.figure
    Figure.savefig = capture_save
    try:
        runpy.run_path(str(script), run_name="__main__")
    finally:
        plt.figure = orig_figure
        Figure.savefig = orig_savefig
    return held["figs"]


def centre_panel_titles(fig) -> None:
    """Centre MULTI-LINE panel titles.

    matplotlib's Text default is multialignment="left": with ha="center" the title BLOCK is
    centred but its lines are left-aligned inside the block, so a two-line title looks
    left-aligned and the standalone panel letter ends up hugging its first character.
    Applied to every axes (including titles created inside loops) right before rendering.
    """
    for ax in fig.axes:
        if ax.get_title():
            ax.title.set_multialignment("center")


def render_at_width(script: Path | str, out_png: Path | str, mm: float = 180.0,
                    dpi: int = 330, pad_in: float = 0.085, max_iter: int = 8,
                    tol_mm: float = 0.8, verbose: bool = True, which: int = 0) -> dict:
    # pad_in = 0.045 in ~ 1.14 mm of white margin on every side. The panel letter sits
    # ABOVE the title, so at the previous 0.01 in it touched the top edge of the canvas.
    """Iterate the figsize scale until the tight-bbox width equals `mm`.

    `which` selects the figure index for scripts that save more than one.
    """
    script, out_png = Path(script), Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    target_in = mm / 25.4 - 2 * pad_in
    scale = 1.0
    best = None
    history: list[float] = []

    for it in range(1, max_iter + 1):
        figs = _run_capturing(script, scale)
        if not figs:
            raise RuntimeError(f"{script.name}: script did not save a figure")
        if which >= len(figs):
            raise RuntimeError(f"{script.name}: only {len(figs)} figure(s) saved, "
                               f"requested index {which}")
        for extra in figs:
            if extra is not figs[which]:
                plt.close(extra)
        fig = figs[which]
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        if not isinstance(fig.canvas, FigureCanvasAgg):
            fig.set_canvas(FigureCanvasAgg(fig))
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        w_in = bb.width + 2 * pad_in          # the delivered canvas includes the margin
        err_mm = (w_in - target_in) * 25.4
        if verbose:
            print(f"    iter {it}: scale={scale:.3f} -> {w_in*25.4:6.1f} mm "
                  f"(err {err_mm:+.1f})")
        best = (fig, bb, w_in)
        if abs(err_mm) <= tol_mm:
            break
        plt.close(fig)
        prev = scale
        scale *= target_in / w_in
        if it >= 3 and scale < prev and w_in >= history[-1]:
            raise RuntimeError(
                f"width does not converge: shrinking the canvas makes the tight bbox "
                f"wider ({history[-1]*25.4:.0f} -> {w_in*25.4:.0f} mm). This figure's "
                f"density needs a layout change, not a scale change.")
        history.append(w_in)
    else:
        fig, bb, w_in = best

    # pad_inches only applies when bbox_inches='tight'; with an explicit Bbox matplotlib
    # ignores it, so the margin has to be added to the Bbox itself. This is what gives the
    # panel letter (drawn ABOVE the title) room instead of touching the canvas edge.
    from matplotlib.transforms import Bbox
    bb_padded = Bbox.from_extents(bb.x0 - pad_in,
                                  bb.y0 - pad_in,
                                  bb.x1 + pad_in,
                                  bb.y1 + pad_in)
    fig.savefig(out_png, dpi=dpi, bbox_inches=bb_padded,
                facecolor="white")
    plt.close(fig)
    from PIL import Image
    im = Image.open(out_png)
    real_mm = im.width / dpi * 25.4
    if verbose:
        print(f"    -> {out_png.name}: {im.width}x{im.height} px = {real_mm:.0f} mm")
    return {"out": str(out_png), "px": f"{im.width}x{im.height}",
            "mm": round(real_mm, 1), "kb": round(out_png.stat().st_size / 1024)}


def panel_letter_fixed(ax, letter, pad_in=0.30, dy=1.02, case="upper", fontsize=None):
    """Bold panel letter at a FIXED distance left of the axes' left edge.

    panel_letter() offsets by a fraction of the axes' width, so a narrow column puts the
    letter almost on the axis and a wide one pushes it far out -- the first column of each
    figure then sits at a different x. This puts the letter at a constant inches offset in
    FIGURE coordinates instead, so the letters of the left-hand column line up across
    figures that share a page width.

    pad_in: inches from the axes' left spine to the letter's left edge.
    """
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    fig = ax.figure
    box = ax.get_position()
    W, H = fig.get_size_inches()
    x = box.x0 - pad_in / W
    # y follows THIS panel's top edge (dy is in axes-fraction units, as in panel_letter),
    # not the figure top. Passing a figure-fraction dy put every letter in one row at the
    # top of the figure -- which is what happened to Fig. 5.
    y = box.y1 + (dy - 1.0) * box.height
    ax.text(x, y, s, transform=fig.transFigure,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")
    return x
