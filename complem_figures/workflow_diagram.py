"""
workflow_diagram.py
============================
Publication workflow figure of the catalog construction.

Draws the seven processing stages top to bottom, data products alternating
with processes, the five source agencies with their time spans on top, and the
side inputs (station inventory, magnitude conversion models, 1-D velocity
model, external picks) feeding the stages that use them. Static content: no
pipeline file is read. Drawn at print size (double column by default) and
written as a vector PDF and a 300 dpi PNG.

``--horizontal`` draws the slide version instead: 16:9, stages left to right,
data products as arrow labels, no citations.

Usage
-----
    python complem_figures/workflow_diagram.py --width-cm 17.4
    python complem_figures/workflow_diagram.py --horizontal

    # Minimal / all defaults
    python complem_figures/workflow_diagram.py
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path

_MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)
_DEFAULT_OUTPUT = os.path.join(_MODULE_DIR, "workflow", "workflow.pdf")
_DEFAULT_OUTPUT_H = os.path.join(_MODULE_DIR, "workflow",
                                 "workflow_horizontal.pdf")
SLIDE_W_CM = 33.867                # PowerPoint 16:9 slide width

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

# One style per node role; fills from Okabe-Ito tints, spread in luminance so
# the roles stay apart in greyscale print.
STYLES = {
    "input":   dict(fc="#D4E9F7", ec="#2F77AD", lw=1.1, round=0.15),
    "process": dict(fc="#F3CF8E", ec="#9A6812", lw=0.8, round=0.0),
    "sub":     dict(fc="#FBEACB", ec="#9A6812", lw=0.7, round=0.0),
    "data":    dict(fc="#FFFFFF", ec="#707070", lw=0.8, round=0.15),
    "final":   dict(fc="#9ED6C0", ec="#16694F", lw=1.5, round=0.15),
}
BADGE_COLOR = "#7A4F0C"
ARROW_COLOR = "#333333"
BUS_COLOR = "#8C8C8C"
CITE_COLOR = "#555555"
FS_PROCESS = 8.0
FS_DATA = 7.5
FS_SMALL = 7.0

# ---------------------------------------------------------------------------
# Layout (cm, y pointing down)
# ---------------------------------------------------------------------------

MAIN_X, MAIN_W = 9.6, 7.0          # main chain
LEFT_X, LEFT_W = 2.35, 4.1         # station inventory branch
RIGHT_X, RIGHT_W = 15.3, 3.5      # magnitude models, velocity model, picks
GAP = 0.45                         # vertical gap between consecutive nodes
TOP = 0.15

AGENCIES = [("OMP", "1978–2019"), ("RESIF", "2020–2025"),
            ("ICGC", "2020–2025"), ("IGN", "2020–2025"),
            ("LDG", "2020–2025")]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

class _Box:
    """Axis-aligned box given by centre and size, with anchor helpers."""

    def __init__(self, x, y, w, h):
        self.x, self.y, self.w, self.h = x, y, w, h

    @property
    def top(self):
        return (self.x, self.y - self.h / 2)

    @property
    def bottom(self):
        return (self.x, self.y + self.h / 2)

    @property
    def left(self):
        return (self.x - self.w / 2, self.y)

    @property
    def right(self):
        return (self.x + self.w / 2, self.y)


def _rect(ax, box, role, stack=0, z=2):
    """Draw a box in the style of ``role``; ``stack`` adds offset copies."""
    st = STYLES[role]
    style = (f"round,pad=0,rounding_size={st['round']}" if st["round"]
             else "square,pad=0")
    for k in range(stack, -1, -1):
        d = 0.09 * k
        ax.add_patch(FancyBboxPatch(
            (box.x - box.w / 2 + d, box.y - box.h / 2 - d), box.w, box.h,
            boxstyle=style, fc=st["fc"], ec=st["ec"], lw=st["lw"],
            zorder=z - 0.1 * k))


def _node(ax, box, role, label, fs=None, weight="normal", stack=0):
    """Box with a centred label."""
    _rect(ax, box, role, stack=stack)
    ax.text(box.x, box.y, label, ha="center", va="center",
            fontsize=fs or FS_DATA, fontweight=weight, zorder=5,
            linespacing=1.15)
    return box


def _badge(ax, x, y, n):
    """Numbered stage badge."""
    ax.add_patch(Circle((x, y), 0.19, fc=BADGE_COLOR, ec="none", zorder=6))
    ax.text(x, y, str(n), ha="center", va="center", fontsize=FS_SMALL,
            fontweight="bold", color="white", zorder=7)


def _process(ax, box, n, label):
    """Process box: number badge on the left, name centred right of it."""
    _rect(ax, box, "process")
    ax.text(box.x + 0.25, box.y, label, ha="center", va="center",
            fontsize=FS_PROCESS, fontweight="bold", zorder=5)
    _badge(ax, box.x - box.w / 2 + 0.32, box.y, n)
    return box


def _container(ax, box, n, label, cite, header_h=0.55):
    """Stage drawn as a stack of 6 zone boxes holding sub-steps."""
    _rect(ax, box, "process", stack=2)
    ax.text(box.x + box.w / 2 - 0.15, box.y + box.h / 2 - 0.12, cite,
            ha="right", va="bottom", fontsize=FS_SMALL, color=CITE_COLOR,
            zorder=5)
    yh = box.y - box.h / 2 + header_h / 2 + 0.05
    ax.text(box.x, yh, label, ha="center", va="center", fontsize=FS_PROCESS,
            fontweight="bold", zorder=5)
    _badge(ax, box.x - box.w / 2 + 0.32, yh, n)
    # Outside, just above the top-right corner of the back copy (offset 2 x 0.09).
    ax.text(box.x + box.w / 2 + 0.18, box.y - box.h / 2 - 0.18 - 0.06,
            "× 6 zones", ha="right", va="bottom", fontsize=FS_DATA,
            fontweight="bold", color=BADGE_COLOR, zorder=5)


def _arrow(ax, *points, color=ARROW_COLOR, lw=0.8, both=False):
    """Polyline through ``points`` ending in an arrowhead."""
    path = Path(list(points), [Path.MOVETO] + [Path.LINETO] * (len(points) - 1))
    ax.add_patch(FancyArrowPatch(
        path=path, arrowstyle="<|-|>" if both else "-|>", mutation_scale=7,
        color=color, lw=lw, shrinkA=0, shrinkB=0, joinstyle="miter",
        capstyle="butt", zorder=3))


def _line(ax, *points, color=ARROW_COLOR, lw=0.8):
    xs, ys = zip(*points)
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="butt",
            solid_joinstyle="miter", zorder=3)


def _dot(ax, x, y, color):
    ax.add_patch(Circle((x, y), 0.045, fc=color, ec="none", zorder=4))


def _rescale(ax, s):
    """Scale text, line widths and arrowheads drawn on ``ax`` by ``s``."""
    for t in ax.texts:
        t.set_fontsize(t.get_fontsize() * s)
    for p in ax.patches:
        p.set_linewidth(p.get_linewidth() * s)
        if isinstance(p, FancyArrowPatch):
            p.set_mutation_scale(p.get_mutation_scale() * s)
    for line in ax.lines:
        line.set_linewidth(line.get_linewidth() * s)


def _save(fig, output, dpi):
    """Write ``output`` (PDF) and a PNG beside it, then close the figure."""
    png = os.path.splitext(output)[0] + ".png"
    fig.savefig(output)
    fig.savefig(png, dpi=dpi)
    plt.close(fig)
    print(f"Wrote {output}\nWrote {png}")
    return {"output": output, "png": png}


def _stack_main(heights):
    """Centre y of each main-chain node, stacked from TOP with GAP."""
    ys, y = [], TOP
    for h in heights:
        ys.append(y + h / 2)
        y += h + GAP
    return ys, y - GAP


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot_workflow(output=None, width_cm=17.4, dpi=300):
    """
    Draw the catalog-construction workflow figure.

    Parameters
    ----------
    output : str, optional
        Output PDF path; a PNG is written next to it. Default
        ``complem_figures/workflow/workflow.pdf``.
    width_cm : float
        Printed width. The layout is designed for 17.4 cm and is scaled
        uniformly for any other width.
    dpi : int
        PNG resolution.

    Returns
    -------
    dict
        ``{'output': pdf_path, 'png': png_path}``.
    """
    output = output or _DEFAULT_OUTPUT
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)

    # Main-chain heights: sources, P2, D, P3, D, P4, D, P5, D, P6, D, P7, final
    heights = [0.8, 0.65, 0.55, 0.65, 0.55, 2.15, 0.55, 0.65, 0.55, 2.35,
               0.55, 0.65, 0.95]
    ys, bottom = _stack_main(heights)
    (y_src, y_p2, y_dsrc, y_p3, y_duni, y_p4, y_drel, y_p5, y_daug, y_p6,
     y_dssst, y_p7, y_final) = ys
    W, H = 17.4, bottom + 0.15

    scale = width_cm / W
    fig = plt.figure(figsize=(W * scale / 2.54, H * scale / 2.54))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.set_axis_off()

    data_w = 3.9

    # --- sources and stage 2 ------------------------------------------------
    n, aw, ag = len(AGENCIES), 1.75, 0.2
    x0 = MAIN_X - (n * aw + (n - 1) * ag) / 2 + aw / 2
    y_collect = y_src + heights[0] / 2 + GAP / 2
    src_boxes = []
    for i, (name, span) in enumerate(AGENCIES):
        b = _Box(x0 + i * (aw + ag), y_src, aw, heights[0])
        _rect(ax, b, "input")
        ax.text(b.x, b.y - 0.15, name, ha="center", va="center",
                fontsize=FS_DATA, fontweight="bold", zorder=5)
        ax.text(b.x, b.y + 0.17, span, ha="center", va="center",
                fontsize=FS_SMALL, zorder=5)
        src_boxes.append(b)
        _line(ax, b.bottom, (b.x, y_collect))
    _line(ax, (src_boxes[0].x, y_collect), (src_boxes[-1].x, y_collect))

    p2 = _process(ax, _Box(MAIN_X, y_p2, MAIN_W, heights[1]), 2,
                  "Bulletin retrieval & conversion")
    _dot(ax, MAIN_X, y_collect, ARROW_COLOR)
    _arrow(ax, (MAIN_X, y_collect), p2.top)

    d_src = _node(ax, _Box(MAIN_X, y_dsrc, data_w, heights[2]), "data",
                  "Source bulletins")
    _arrow(ax, p2.bottom, d_src.top)

    # --- stage 1: station inventory ------------------------------------------
    meta = _Box(LEFT_X, y_src, LEFT_W, heights[0])
    _rect(ax, meta, "input")
    ax.text(meta.x, meta.y, "All station inventories", ha="center", va="center",
            fontsize=FS_DATA, fontweight="bold", zorder=5)
    p1 = _process(ax, _Box(LEFT_X, y_p2, LEFT_W, heights[1]), 1,
                  "Station inventory fusion")
    _arrow(ax, meta.bottom, p1.top)
    inv = _node(ax, _Box(LEFT_X, y_dsrc, LEFT_W, heights[2]), "data",
                "Unified station inventory")
    _arrow(ax, p1.bottom, inv.top)

    # --- stage 3 -------------------------------------------------------------
    p3 = _process(ax, _Box(MAIN_X, y_p3, MAIN_W, heights[3]), 3,
                  "Bulletin harmonization")
    _arrow(ax, d_src.bottom, p3.top)
    mag = _node(ax, _Box(RIGHT_X, y_p3, RIGHT_W, 0.75), "data",
                "Magnitude\nconversion models")
    _arrow(ax, p3.right, mag.left, both=True)

    d_uni = _node(ax, _Box(MAIN_X, y_duni, data_w, heights[4]), "data",
                  "Unified bulletin")
    _arrow(ax, p3.bottom, d_uni.top)

    # --- stage 4: NonLinLoc, two passes, 6 zones ------------------------------
    p4 = _Box(MAIN_X, y_p4, MAIN_W, heights[5])
    _container(ax, p4, 4, "Absolute relocation (NonLinLoc)",
               "Lomax et al. (2000)")
    _arrow(ax, d_uni.bottom, p4.top)
    y_sub = p4.top[1] + 1.175
    pass1 = _node(ax, _Box(MAIN_X - 2.35, y_sub, 1.5, 0.6), "sub", "Pass 1")
    delays = _node(ax, _Box(MAIN_X, y_sub, 2.1, 0.75), "data",
                   "Station delays\n(from residuals)")
    pass2 = _node(ax, _Box(MAIN_X + 2.35, y_sub, 1.5, 0.6), "sub", "Pass 2")
    _arrow(ax, pass1.right, delays.left)
    _arrow(ax, delays.right, pass2.left)

    vel = _Box(RIGHT_X, y_p4, RIGHT_W, 1.1)
    _rect(ax, vel, "input")
    ax.text(vel.x, vel.y - 0.27, "1-D velocity model", ha="center",
            va="center", fontsize=FS_DATA, fontweight="bold", zorder=5)
    ax.text(vel.x, vel.y + 0.15,
            "Souriau & Pauchet (1998)\nVp/Vs: Derode et al. (2023)",
            ha="center", va="center", fontsize=FS_SMALL, color=CITE_COLOR,
            linespacing=1.15, zorder=5)
    _arrow(ax, vel.left, p4.right)

    d_rel = _node(ax, _Box(MAIN_X, y_drel, data_w, heights[6]), "data",
                  "Relocated bulletin")
    _arrow(ax, p4.bottom, d_rel.top)

    # --- stage 5: external picks ----------------------------------------------
    p5 = _process(ax, _Box(MAIN_X, y_p5, MAIN_W, heights[7]), 5,
                  "Pick augmentation")
    _arrow(ax, d_rel.bottom, p5.top)
    ext = _node(ax, _Box(RIGHT_X, y_drel, RIGHT_W, heights[6]), "input",
                "External picks", weight="bold")
    x_ext = MAIN_X + MAIN_W / 2 - 0.6
    _arrow(ax, ext.left, (x_ext, ext.y), (x_ext, p5.top[1]))

    d_aug = _node(ax, _Box(MAIN_X, y_daug, data_w, heights[8]), "data",
                  "Augmented bulletin")
    _arrow(ax, p5.bottom, d_aug.top)

    # --- stage 6: SSST, iterative, 6 zones ------------------------------------
    p6 = _Box(MAIN_X, y_p6, MAIN_W, heights[9])
    _container(ax, p6, 6, "SSST relocation (NonLinLoc)",
               "Lomax & Savvaidis (2022)")
    _arrow(ax, d_aug.bottom, p6.top)
    y_sub = p6.top[1] + 1.075
    loc = _node(ax, _Box(MAIN_X - 1.45, y_sub, 1.7, 0.6), "sub", "Location")
    ssst = _node(ax, _Box(MAIN_X + 1.45, y_sub, 2.3, 0.75), "data",
                 "Source-specific\nstation terms")
    _arrow(ax, loc.right, ssst.left)
    y_loop = y_sub + 0.62
    _arrow(ax, ssst.bottom, (ssst.x, y_loop), (loc.x, y_loop), loc.bottom)
    ax.text(MAIN_X, y_loop + 0.05, "iterate", ha="center", va="top",
            fontsize=FS_SMALL, style="italic", zorder=5)

    d_ssst = _node(ax, _Box(MAIN_X, y_dssst, data_w, heights[10]), "data",
                   "Hypocenters + location PDFs")
    _arrow(ax, p6.bottom, d_ssst.top)

    # --- stage 7: export --------------------------------------------------------
    p7 = _process(ax, _Box(MAIN_X, y_p7, MAIN_W, heights[11]), 7,
                  "Catalog export")
    _arrow(ax, d_ssst.bottom, p7.top)
    fin = _Box(MAIN_X, y_final, 5.6, heights[12])
    _rect(ax, fin, "final")
    ax.text(fin.x, fin.y - 0.17, "Final catalog (QuakeML)", ha="center",
            va="center", fontsize=FS_PROCESS, fontweight="bold", zorder=5)
    ax.text(fin.x, fin.y + 0.2,
            "PDF-expectation hypocenters · quality flag", ha="center",
            va="center", fontsize=FS_SMALL, zorder=5)
    _arrow(ax, p7.bottom, fin.top)

    # --- buses: station inventory (left), velocity model (right) --------------
    targets = [p3.left, p4.left, p5.left, p6.left, p7.left]
    _line(ax, inv.bottom, (LEFT_X, targets[-1][1]), color=BUS_COLOR, lw=0.7)
    for i, (tx, ty) in enumerate(targets):
        if i < len(targets) - 1:
            _dot(ax, LEFT_X, ty, BUS_COLOR)
        _arrow(ax, (LEFT_X, ty), (tx, ty), color=BUS_COLOR, lw=0.7)

    x_bus = W - 0.15
    _line(ax, vel.right, (x_bus, vel.y), (x_bus, p6.y), color=BUS_COLOR,
          lw=0.7)
    _dot(ax, x_bus, p5.y, BUS_COLOR)
    _arrow(ax, (x_bus, p5.y), p5.right, color=BUS_COLOR, lw=0.7)
    _arrow(ax, (x_bus, p6.y), p6.right, color=BUS_COLOR, lw=0.7)

    # Sizes are set for 17.4 cm; keep them in proportion at other widths.
    _rescale(ax, scale)
    return _save(fig, output, dpi)


def plot_workflow_horizontal(output=None, width_cm=SLIDE_W_CM, dpi=300):
    """
    Draw the slide (16:9, left-to-right) version of the workflow figure.

    Stages 2-7 run in one row, intermediate data products are labels on the
    arrows, sources and the station-inventory branch sit above, the side
    inputs and the final catalog below. No citations. The top band is left
    empty for a slide title.

    Parameters
    ----------
    output : str, optional
        Output PDF path; a PNG is written next to it. Default
        ``complem_figures/workflow/workflow_horizontal.pdf``.
    width_cm : float
        Slide width (default 33.867 cm, the PowerPoint 16:9 slide). The
        height follows the 16:9 ratio.
    dpi : int
        PNG resolution.

    Returns
    -------
    dict
        ``{'output': pdf_path, 'png': png_path}``.
    """
    output = output or _DEFAULT_OUTPUT_H
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)

    # Designed with the paper's sizes on a 22.58 cm canvas, then scaled up
    # (x1.5 on a 33.9 cm slide) so text, lines and arrowheads grow together.
    W = 22.58
    H = W * 9 / 16
    scale = width_cm / W
    fig = plt.figure(figsize=(width_cm / 2.54, H * scale / 2.54))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.set_axis_off()

    y0 = 2.0                            # top band left free for a title
    y_src, y_collect, y_bus = y0 + 1.0, y0 + 1.75, y0 + 2.5
    y_row, y_side, y_vbus = y0 + 5.4, y0 + 8.6, y0 + 9.7
    h_proc, h_cont = 1.15, 4.1

    # --- main row: stages 2..7 ------------------------------------------------
    widths = [2.2, 2.2, 2.5, 2.2, 2.5, 2.2]
    margin = 0.2
    gap = (W - 2 * margin - sum(widths)) / (len(widths) - 1)
    xs, x = [], margin
    for w in widths:
        xs.append(x + w / 2)
        x += w + gap

    def corner_process(box, n, label):
        _rect(ax, box, "process")
        ax.text(box.x, box.y, label, ha="center", va="center",
                fontsize=FS_PROCESS, fontweight="bold", linespacing=1.1,
                zorder=5)
        _badge(ax, box.x - box.w / 2, box.top[1], n)
        return box

    def corner_container(box, n, label):
        _rect(ax, box, "process", stack=2)
        ax.text(box.x, box.top[1] + 0.62, label, ha="center", va="center",
                fontsize=FS_PROCESS, fontweight="bold", linespacing=1.1,
                zorder=5)
        _badge(ax, box.x - box.w / 2, box.top[1], n)
        ax.text(box.x + box.w / 2 + 0.18, box.top[1] - 0.24, "× 6 zones",
                ha="right", va="bottom", fontsize=FS_DATA, fontweight="bold",
                color=BADGE_COLOR, zorder=5)
        return box

    p2 = corner_process(_Box(xs[0], y_row, widths[0], h_proc), 2,
                        "Bulletin\nretrieval &\nconversion")
    p3 = corner_process(_Box(xs[1], y_row, widths[1], h_proc), 3,
                        "Bulletin\nharmonization")
    p4 = corner_container(_Box(xs[2], y_row, widths[2], h_cont), 4,
                          "Absolute\nrelocation\n(NonLinLoc)")
    p5 = corner_process(_Box(xs[3], y_row, widths[3], h_proc), 5,
                        "Pick\naugmentation")
    p6 = corner_container(_Box(xs[4], y_row, widths[4], h_cont), 6,
                          "SSST\nrelocation\n(NonLinLoc)")
    p7 = corner_process(_Box(xs[5], y_row, widths[5], h_proc), 7,
                        "Catalog\nexport")
    row = [p2, p3, p4, p5, p6, p7]

    # Stage 4: two passes with station delays, stacked vertically.
    top = p4.top[1]
    pass1 = _node(ax, _Box(p4.x, top + 1.6, 1.6, 0.5), "sub", "Pass 1")
    delays = _node(ax, _Box(p4.x, top + 2.55, 1.8, 0.75), "data",
                   "Station\ndelays")
    pass2 = _node(ax, _Box(p4.x, top + 3.5, 1.6, 0.5), "sub", "Pass 2")
    _arrow(ax, pass1.bottom, delays.top)
    _arrow(ax, delays.bottom, pass2.top)

    # Stage 6: location <-> station terms, iterated.
    top = p6.top[1]
    loc = _node(ax, _Box(p6.x, top + 1.65, 1.7, 0.5), "sub", "Location")
    terms = _node(ax, _Box(p6.x, top + 2.9, 2.2, 0.75), "data",
                  "Source-specific\nstation terms")
    _arrow(ax, (p6.x - 0.35, loc.bottom[1]), (p6.x - 0.35, terms.top[1]))
    _arrow(ax, (p6.x + 0.35, terms.top[1]), (p6.x + 0.35, loc.bottom[1]))
    ax.text(p6.x + 0.45, (loc.bottom[1] + terms.top[1]) / 2, "iterate",
            ha="left", va="center", fontsize=FS_SMALL, style="italic",
            zorder=5)

    # Arrows between stages, labelled underneath (the corner badges sit above)
    # with the data product they carry.
    products = ["Source\nbulletins", "Unified\nbulletin", "Relocated\nbulletin",
                "Augmented\nbulletin", "Hypocenters\n+ PDFs"]
    for a, b, label in zip(row[:-1], row[1:], products):
        _arrow(ax, (a.right[0], y_row), (b.left[0], y_row))
        x_mid = (a.right[0] + 0.09 + b.left[0]) / 2
        ax.text(x_mid, y_row + 0.12, label, ha="center", va="top",
                fontsize=FS_SMALL, linespacing=1.1, zorder=5)

    # --- top band: sources -> stage 2, inventory branch -> bus ----------------
    aw, ag = 1.7, 0.18
    src = [_Box(margin + aw / 2 + i * (aw + ag), y_src, aw, 0.8)
           for i in range(len(AGENCIES))]
    for b, (name, span) in zip(src, AGENCIES):
        _rect(ax, b, "input")
        ax.text(b.x, b.y - 0.15, name, ha="center", va="center",
                fontsize=FS_DATA, fontweight="bold", zorder=5)
        ax.text(b.x, b.y + 0.17, span, ha="center", va="center",
                fontsize=FS_SMALL, zorder=5)
        _line(ax, b.bottom, (b.x, y_collect))
    _line(ax, (src[0].x, y_collect), (src[-1].x, y_collect))
    _dot(ax, p2.x, y_collect, ARROW_COLOR)
    _arrow(ax, (p2.x, y_collect), p2.top)

    inv = _node(ax, _Box(19.2, y_src, 2.2, 0.8), "data",
                "Unified station\ninventory")
    p1 = corner_process(_Box(inv.x - 1.1 - 0.45 - 1.45, y_src, 2.9, 0.8), 1,
                        "Station inventory\nfusion")
    allinv = _Box(p1.x - 1.45 - 0.45 - 1.1, y_src, 2.2, 0.8)
    _rect(ax, allinv, "input")
    ax.text(allinv.x, allinv.y, "All station\ninventories", ha="center",
            va="center", fontsize=FS_DATA, fontweight="bold", linespacing=1.1,
            zorder=5)
    _arrow(ax, allinv.right, p1.left)
    _arrow(ax, p1.right, inv.left)

    _line(ax, inv.bottom, (inv.x, y_bus), color=BUS_COLOR, lw=0.7)
    _line(ax, (p3.x, y_bus), (p7.x, y_bus), color=BUS_COLOR, lw=0.7)
    _dot(ax, inv.x, y_bus, BUS_COLOR)
    for i, p in enumerate([p3, p4, p5, p6, p7]):
        if 0 < i < 4:
            _dot(ax, p.x, y_bus, BUS_COLOR)
        _arrow(ax, (p.x, y_bus), p.top, color=BUS_COLOR, lw=0.7)

    # --- bottom band: side inputs, velocity bus, final catalog ----------------
    mag = _node(ax, _Box(p3.x, y_side, 2.4, 0.75), "data",
                "Magnitude\nconversion models")
    _arrow(ax, p3.bottom, mag.top, both=True)

    vel = _Box(p4.x, y_side, 2.2, 0.75)
    _rect(ax, vel, "input")
    ax.text(vel.x, vel.y, "1-D velocity\nmodel", ha="center", va="center",
            fontsize=FS_DATA, fontweight="bold", linespacing=1.1, zorder=5)
    _arrow(ax, vel.top, p4.bottom)

    ext = _Box(p5.x + 0.55, y_side, 2.0, 0.6)
    _rect(ax, ext, "input")
    ax.text(ext.x, ext.y, "External picks", ha="center", va="center",
            fontsize=FS_DATA, fontweight="bold", zorder=5)
    _arrow(ax, ext.top, (ext.x, p5.bottom[1]))

    x_p5 = p5.x - 0.75
    _line(ax, vel.bottom, (vel.x, y_vbus), (p6.x, y_vbus), color=BUS_COLOR,
          lw=0.7)
    _dot(ax, x_p5, y_vbus, BUS_COLOR)
    _arrow(ax, (x_p5, y_vbus), (x_p5, p5.bottom[1]), color=BUS_COLOR, lw=0.7)
    _arrow(ax, (p6.x, y_vbus), p6.bottom, color=BUS_COLOR, lw=0.7)

    fin = _Box(min(p7.x, W - 0.15 - 1.5), y_side + 0.35, 3.0, 1.5)
    _rect(ax, fin, "final")
    ax.text(fin.x, fin.y - 0.3, "Final catalog\n(QuakeML)", ha="center",
            va="center", fontsize=FS_PROCESS, fontweight="bold",
            linespacing=1.1, zorder=5)
    ax.text(fin.x, fin.y + 0.4, "expectation hypocenters\n+ quality flag",
            ha="center", va="center", fontsize=FS_SMALL, linespacing=1.1,
            zorder=5)
    _arrow(ax, (fin.x, p7.bottom[1]), fin.top)

    _rescale(ax, scale)
    return _save(fig, output, dpi)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Draw the catalog-construction workflow figure.")
    parser.add_argument("--output", default=None,
                        help="output PDF (PNG written alongside); default "
                             "complem_figures/workflow/workflow.pdf, or "
                             "workflow_horizontal.pdf with --horizontal")
    parser.add_argument("--width-cm", type=float, default=None,
                        help="width in cm (default 17.4, or 33.867 = a 16:9 "
                             "slide with --horizontal)")
    parser.add_argument("--horizontal", action="store_true",
                        help="slide version: 16:9, left to right, no citations")
    args = parser.parse_args()
    if args.horizontal:
        plot_workflow_horizontal(args.output, args.width_cm or SLIDE_W_CM)
    else:
        plot_workflow(args.output, args.width_cm or 17.4)


if __name__ == "__main__":
    main()
