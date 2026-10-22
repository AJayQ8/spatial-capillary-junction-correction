#!/usr/bin/env python3
"""Build a compact, editable journal schematic for the spatial correction."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from matplotlib.patches import Circle, Polygon, Rectangle

WIDTH, HEIGHT = 6.70, 2.70
INK = "#252D32"
MUTED = "#626D75"
INTERFACE = "#747F86"
BORDER = "#A2ABB0"
WINDOW = "#277C82"
SOURCE = "#765A83"
CAPILLARY = "#2B6F91"
BULK = "#4D585F"
PHASE = {
    0: "#E4EDF3",
    1: "#E8EEE9",
    2: "#F1E9E0",
}


def text(ax, x, y, s, *, size, color=INK, ha="left", va="center", weight="normal"):
    return ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va,
                   fontweight=weight, zorder=8)


def panel_title(ax, x, letter, title):
    text(ax, x, 2.47, f"({letter})  {title}", size=9.5, weight="semibold")


def draw_y_field(ax, left: float, *, window: bool = False):
    """Same rectangular, straight-Y field in panels a and b."""
    xmin = left + 0.08
    xmax = xmin + 1.90
    ymin, ymax = 0.83, 2.08
    ymid = (ymin + ymax) / 2
    node = (xmin + 0.92, ymid)
    upper = (node[0] + 0.60, ymax)
    lower = (node[0] + 0.60, ymin)

    # Three adjoining regions, separated by a simple Y-shaped interface.
    regions = [
        (0, [(xmin, ymid), (xmin, ymax), upper, node]),
        (1, [upper, (xmax, ymax), (xmax, ymin), lower, node]),
        (2, [(xmin, ymin), (xmin, ymid), node, lower]),
    ]
    for phase, vertices in regions:
        ax.add_patch(Polygon(vertices, closed=True, facecolor=PHASE[phase],
                             edgecolor="none", zorder=1))

    for start, end in [((xmin, ymid), node), (node, upper), (node, lower)]:
        ax.plot([start[0], end[0]], [start[1], end[1]], color=INTERFACE,
                linewidth=0.62, solid_capstyle="round", zorder=3)
    ax.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,
                           facecolor="none", edgecolor=BORDER,
                           linewidth=0.55, zorder=4))

    # Labels are centered inside their respective phase regions.
    text(ax, xmin + 0.39, ymid + 0.34, r"$\phi_0$", size=8.0,
         color=INK, ha="center")
    text(ax, xmax - 0.31, ymid + 0.01, r"$\phi_1$", size=8.0,
         color=INK, ha="center")
    text(ax, xmin + 0.39, ymid - 0.34, r"$\phi_2$", size=8.0,
         color=INK, ha="center")

    # Compact ternary source marker at the junction.
    ax.add_patch(Circle(node, 0.042, facecolor=SOURCE, edgecolor="white",
                        linewidth=0.45, zorder=6))

    if window:
        # Receiver sits on the left binary arm; the square contains the source.
        receiver = (xmin + 0.64, ymid)
        half = 0.47
        ax.add_patch(Rectangle((receiver[0] - half, receiver[1] - half),
                               2 * half, 2 * half,
                               facecolor="none", edgecolor=WINDOW,
                               linewidth=0.95, zorder=5))
        ax.add_patch(Circle(receiver, 0.030, facecolor="white",
                            edgecolor=WINDOW, linewidth=0.9, zorder=7))
        text(ax, receiver[0] - 0.045, receiver[1] + 0.080, r"$x$",
             size=8.1, color=INK, ha="right")


def equation_row(ax, x, y, expressions, font=10.5):
    """Naturally space math terms as a single horizontally packed equation."""
    children = [
        TextArea(expr, textprops={"fontsize": font, "color": color,
                                  "fontfamily": "DejaVu Sans"})
        for expr, color in expressions
    ]
    row = HPacker(children=children, align="baseline", pad=0, sep=4.0)
    ax.add_artist(AnnotationBbox(row, (x, y), xycoords="data",
                                 box_alignment=(0.5, 0.5), frameon=False,
                                 pad=0, annotation_clip=False, zorder=8))


def draw_force_rule(ax, center: float):
    # A vertical before/after sequence leaves room for naturally spaced math.
    text(ax, center, 2.02, "native", size=8.0, color=MUTED, ha="center")
    equation_row(ax, center, 1.75,
                [(r"$C$", CAPILLARY), (r"$+$", INK), (r"$B$", BULK)],
                font=11.0)
    ax.annotate("", xy=(center, 1.13), xytext=(center, 1.45),
                arrowprops={"arrowstyle": "-|>", "color": MUTED,
                            "lw": 0.75, "mutation_scale": 7.0}, zorder=7)
    text(ax, center, 0.98, "corrected", size=8.0, color=MUTED, ha="center")
    equation_row(ax, center, 0.70,
                [(r"$(1+\chi H_\rho)C$", CAPILLARY),
                 (r"$+$", INK), (r"$B$", BULK)],
                font=10.5)
    text(ax, center, 0.28, "Bulk term unchanged", size=8.0,
         color=MUTED, ha="center")


def make_figure():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(WIDTH, HEIGHT))
    fig.patch.set_facecolor("white")
    ax.set_position((0, 0, 1, 1))
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(0, HEIGHT)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    gap = 0.18
    margin = 0.09
    panel_w = (WIDTH - 2 * margin - 2 * gap) / 3
    lefts = [margin + i * (panel_w + gap) for i in range(3)]
    for left, letter, title in zip(lefts, "abc",
                                  ["Ternary source", "Square neighborhood", "Force correction"]):
        panel_title(ax, left, letter, title)

    draw_y_field(ax, lefts[0])
    draw_y_field(ax, lefts[1], window=True)
    draw_force_rule(ax, lefts[2] + panel_w / 2)

    text(ax, lefts[0] + panel_w / 2, 0.48,
         r"$P=27\,\phi_0\phi_1\phi_2$", size=9.8, ha="center")
    text(ax, lefts[1] + panel_w / 2, 0.52,
         r"$H_\rho(x)=\max_{\|y-x\|_\infty\leq\rho}P(y)$",
         size=8.8, ha="center")
    text(ax, lefts[1] + panel_w / 2, 0.28,
         r"$P(x)=0,\quad H_\rho(x)>0$", size=8.7, ha="center", color=MUTED)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig = make_figure()
    for extension in ("pdf", "svg"):
        fig.savefig(args.output_dir / f"method_schematic.{extension}",
                    facecolor="white", bbox_inches="tight", pad_inches=0.015)
    fig.savefig(args.output_dir / "method_schematic.png",
                facecolor="white", bbox_inches="tight", pad_inches=0.015,
                dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
