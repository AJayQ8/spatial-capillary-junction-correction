#!/usr/bin/env python3
"""Rebuild manuscript result figures from the ten frozen local CSV exports.

No simulation, fitting or private laboratory files are used.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "figure_data"
FIGURES = ROOT / "build" / "figures"

COLORS = {
    "native": "#0072B2",
    "whole": "#D55E00",
    "local": "#009E73",
    "grose": "#CC79A7",
    "dark": "#263238",
    "gray": "#6B7280",
    "light": "#E8EDF2",
    "red": "#B23A48",
    "gold": "#E6AB02",
}


ROUTE_LABELS = {
    "native": "Native B2 (matched D9 dt)",
    "whole": "Whole halo",
    "local": "Local memory",
    "grose": "Grose A1 (source-based reconstruction)",
}


ROUTE_MARKERS = {"native": "o", "whole": "s", "local": "^", "grose": "D"}


def read_export_csv(name: str) -> list[dict[str, Any]]:
    """Read a compact package-local export for a portable figure-only rebuild."""
    path = DATA / name
    with path.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    converted: list[dict[str, Any]] = []
    for row in rows:
        item: dict[str, Any] = {}
        for key, value in row.items():
            if value == "":
                item[key] = None
            elif value == "True":
                item[key] = True
            elif value == "False":
                item[key] = False
            else:
                try:
                    item[key] = float(value)
                except (TypeError, ValueError):
                    item[key] = value
        converted.append(item)
    return converted



def save_figure(fig: plt.Figure, path: Path, dpi: int = 320) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".tmp" + path.suffix)
    fmt = path.suffix.lstrip(".")
    metadata = {"CreationDate": None, "ModDate": None} if fmt == "pdf" else None
    fig.savefig(
        tmp,
        format=fmt,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        metadata=metadata,
    )
    with tmp.open("rb") as handle:
        os.fsync(handle.fileno())
    tmp.replace(path)



def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.4,
            "axes.titlesize": 9.3,
            "axes.labelsize": 8.5,
            "legend.fontsize": 7.2,
            "xtick.labelsize": 7.7,
            "ytick.labelsize": 7.7,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.45,
            "lines.markersize": 4.1,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": False,
        }
    )



def make_n20_figure(n20: dict[str, Any], profile_rows: list[dict[str, Any]], shape_rows: list[dict[str, Any]], metric_rows: list[dict[str, Any]]) -> plt.Figure:
    fig = plt.figure(figsize=(6.55, 5.35), constrained_layout=False)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.0))
    axes = {
        "profile": fig.add_subplot(grid[0, :]),
        "shape": fig.add_subplot(grid[1, 0]),
        "speed": fig.add_subplot(grid[1, 1]),
    }
    ax = axes["profile"]
    x = [r["x_W"] for r in profile_rows]
    y0 = [r["original_B2_terminal_aligned_y_W"] for r in profile_rows]
    y1 = [r["corrected_terminal_aligned_y_W"] for r in profile_rows]
    y_sharp = [r["sharp_reference_y_W"] for r in profile_rows]
    ax.plot(x, y0, color=COLORS["native"], label="Original B2")
    ax.plot(x, y1, color=COLORS["whole"], label="Corrected B2")
    ax.plot(x, y_sharp, color=COLORS["dark"], ls="--", lw=1.15, label="Sharp reference")
    ax.set_title("(a) Terminal fixed-92 aligned profile at $t_{rel}=0.25$", loc="center", fontweight="bold")
    ax.set_xlabel("$x/W$")
    ax.set_ylabel("Aligned interface $y/W$")
    ax.grid(alpha=0.22)
    # The profile has a deep central trough; put the key in the clear upper
    # center rather than covering the aligned minima.
    ax.legend(frameon=False, ncol=3, loc="upper center", fontsize=6.9)

    ax = axes["shape"]
    for key, label, color in (("original_B2", "Original B2", COLORS["native"]), ("corrected_V008", "Corrected B2", COLORS["whole"])):
        rows = [r for r in shape_rows if r["route_key"] == key]
        ax.plot([r["relative_time_from_t1"] for r in rows], [r["fixed92_sharp_shape_rms_W"] for r in rows], marker="o", ms=2.9, color=color, label=label)
    ax.axhline(0.05, color=COLORS["red"], lw=1, ls="--", label="0.05 W study target")
    ax.set_title("(b) Fixed-92 RMS shape error", loc="center", fontweight="bold")
    ax.set_xlabel("Relative time from warm source, $t_{rel}$")
    ax.set_ylabel("RMS shape error (W)")
    ax.set_xlim(0, 0.25)
    ax.set_ylim(0, 0.125)
    ax.grid(alpha=0.22)
    # Keep the small study-screen label away from the two trajectories.
    ax.text(0.012, 0.052, "0.05 W target", ha="left", va="bottom", fontsize=6.1,
            color=COLORS["red"], bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 0.8})
    # Direct labels avoid covering the early relaxation curves in this compact panel.
    rows = [r for r in shape_rows if r["route_key"] == "original_B2"]
    if rows:
        last = rows[-1]
        ax.annotate("Original B2", (last["relative_time_from_t1"], last["fixed92_sharp_shape_rms_W"]),
                    xytext=(-5, 8), textcoords="offset points", ha="right", fontsize=6.5, color=COLORS["native"])
    rows = [r for r in shape_rows if r["route_key"] == "corrected_V008"]
    if rows:
        last = rows[-1]
        ax.annotate("Corrected B2", (last["relative_time_from_t1"], last["fixed92_sharp_shape_rms_W"]),
                    xytext=(-5, 8), textcoords="offset points", ha="right", fontsize=6.5, color=COLORS["whole"])

    ax = axes["speed"]
    xpos = [0, 1]
    vals = [r["speed_error_percent"] for r in metric_rows]
    labels = [r["route"] for r in metric_rows]
    colors = [COLORS["native"], COLORS["whole"]]
    ax.bar(xpos, vals, color=colors, width=0.56, edgecolor="white", linewidth=0.8)
    ax.axhline(5.0, color=COLORS["red"], lw=1, ls="--")
    for x0, value in zip(xpos, vals):
        ax.text(x0, value + 0.45, f"{value:.3f}%", ha="center", va="bottom", fontsize=7.4)
    ax.set_xticks(xpos, labels)
    ax.set_title("(c) Late area-speed error", loc="center", fontweight="bold")
    ax.set_ylabel("Relative error vs sharp speed (%)")
    ax.set_ylim(0, max(vals) * 1.25)
    ax.grid(axis="y", alpha=0.22)
    ax.text(0.98, 0.94, "5% study target", ha="right", va="top", transform=ax.transAxes, fontsize=6.8, color=COLORS["red"])
    fig.subplots_adjust(left=0.12, right=0.985, bottom=0.12, top=0.965, wspace=0.36, hspace=0.62)
    return fig



def make_h24_figure(shape_rows: list[dict[str, Any]], metric_rows: list[dict[str, Any]]) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.65), constrained_layout=False)
    routes = ["native", "whole", "local", "grose"]
    for col, pair in enumerate(("short", "long")):
        ax = axes[0, col]
        for route in routes:
            rows = sorted(
                (r for r in shape_rows if r["route_key"] == route and r["pair"] == pair and r["hausdorff_sampled_W"] is not None),
                key=lambda r: r["time_W2"],
            )
            if not rows:
                continue
            times = [r["time_W2"] for r in rows]
            values = [r["hausdorff_sampled_W"] for r in rows]
            lo = [max(0.0, v - r["hausdorff_lower_W"]) for v, r in zip(values, rows)]
            hi = [max(0.0, r["hausdorff_upper_W"] - v) for v, r in zip(values, rows)]
            ax.errorbar(times, values, yerr=[lo, hi], color=COLORS[route], marker=ROUTE_MARKERS[route], ms=3.0, capsize=1.3, lw=1.05, label=ROUTE_LABELS[route])
        ax.axhline(0.05, color=COLORS["red"], lw=0.9, ls="--")
        ax.set_title(f"({'a' if col == 0 else 'b'}) {pair.capitalize()} pair: sampled Hausdorff error", loc="center", fontweight="bold")
        ax.set_xlabel("Physical time ($W^2$)")
        ax.set_ylabel("Sampled Hausdorff error (W)")
        ax.grid(alpha=0.2)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

    pair_centers = {"short": 0, "long": 1}
    offsets = {route: (i - 1.5) * 0.14 for i, route in enumerate(routes)}
    for row_idx, metric_key, title, ylabel, threshold, panel in (
        (1, "speed_error_percent", "Relative area-speed error", "Error vs sharp speed (%)", 5.0, "(c)"),
        (1, "T99_upper_error_percent", "T99 upper-bracket error", "Upper relative error (%)", 5.0, "(d)"),
    ):
        # Two bottom axes, selected by metric for each pass.
        ax = axes[1, 0] if metric_key == "speed_error_percent" else axes[1, 1]
        for route in routes:
            xs, ys = [], []
            for pair in ("short", "long"):
                rec = next(r for r in metric_rows if r["route_key"] == route and r["pair"] == pair)
                xs.append(pair_centers[pair] + offsets[route])
                ys.append(rec[metric_key])
            ax.plot(xs, ys, marker=ROUTE_MARKERS[route], color=COLORS[route], lw=0, ms=4.7, label=ROUTE_LABELS[route])
        ax.axhline(threshold, color=COLORS["red"], lw=0.9, ls="--")
        ax.set_xticks([0, 1], ["Short", "Long"])
        ax.set_title(f"{panel} {title}", loc="center", fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.2)
        ax.set_xlim(-0.38, 1.38)
        ax.set_ylim(bottom=0)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.11, top=0.79, wspace=0.34, hspace=0.62)
    fig.text(0.5, 0.985, "Driven junction-pair comparison",
             ha="center", va="top", fontsize=9.0, fontweight="bold")
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.925), ncol=2, frameon=False, fontsize=7.0)
    return fig



def make_robustness_figure(rows: list[dict[str, Any]]) -> plt.Figure:
    fig, axes = plt.subplots(1, 2, figsize=(6.65, 3.45), constrained_layout=False)
    x = list(range(len(rows)))
    labels = ["(27,10)", "(20,10)", "(34,10)", "(27,8)", "(27,12)"]
    colors = [COLORS["dark"], COLORS["native"], COLORS["native"], COLORS["whole"], COLORS["whole"]]
    markers = ["D", "o", "o", "^", "^"]
    for ax, key, ylabel, threshold, panel in (
        (axes[0], "speed_error_percent", "Late area-speed error (%)", 5.0, "(a)"),
        (axes[1], "terminal_fixed92_rms_W", "Terminal fixed-92 RMS (W)", 0.05, "(b)"),
    ):
        values = [r[key] for r in rows]
        for i, (v, c, m) in enumerate(zip(values, colors, markers)):
            ax.scatter(i, v, s=43, marker=m, facecolor=c, edgecolor="white", linewidth=0.7, zorder=3)
            fmt = f"{v:.2f}%" if key == "speed_error_percent" else f"{v:.4f}"
            ax.annotate(fmt, (i, v), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=6.8, color=COLORS["dark"])
        ax.set_xticks(x, labels)
        # Leave enough space for the full value labels over the first/last points.
        ax.set_xlim(-0.52, 4.52)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{panel} {ylabel}", loc="center", fontweight="bold")
        ax.grid(axis="y", alpha=0.2)
        if key == "speed_error_percent":
            ax.set_ylim(0.9, 1.5)
        else:
            ax.set_ylim(0.0032, 0.0052)
    fig.subplots_adjust(left=0.12, right=0.985, bottom=0.30, top=0.76, wspace=0.38)
    fig.text(0.5, 0.97, "Sensitivity to correction strength and radius", ha="center", va="top", fontsize=9.4, fontweight="bold")
    fig.text(0.5, 0.20, r"Settings shown as $(\chi, r/h)$; each varies one parameter from the common developed-front state.",
             ha="center", fontsize=7.0, color=COLORS["gray"])
    fig.text(0.5, 0.09, "All five settings stayed below the 5% speed and 0.05 W shape study targets; no joint or continuous-range claim.",
             ha="center", fontsize=6.8, color=COLORS["gray"])
    return fig



def make_openphase_figure(shape_rows: list[dict[str, Any]], metric_rows: list[dict[str, Any]]) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.3), constrained_layout=False)
    ax = axes[0, 0]
    for pair, color, marker in (("short", COLORS["native"], "o"), ("long", COLORS["whole"], "s")):
        rows = sorted((r for r in shape_rows if r["pair"] == pair and r["hausdorff_sampled_W"] is not None), key=lambda r: r["time_W2"])
        values = [r["hausdorff_sampled_W"] for r in rows]
        lo = [max(0.0, v - r["hausdorff_lower_W"]) for v, r in zip(values, rows)]
        hi = [max(0.0, r["hausdorff_upper_W"] - v) for v, r in zip(values, rows)]
        ax.errorbar([r["time_W2"] for r in rows], values, yerr=[lo, hi], marker=marker, color=color, capsize=1.5, lw=1.1, ms=3.4, label=f"Halo {pair}")
    ax.axhline(0.05, color=COLORS["red"], lw=0.9, ls="--", label="0.05 W study target")
    long_rows = [r for r in shape_rows if r["pair"] == "long" and r["hausdorff_sampled_W"] is not None]
    late = max(long_rows, key=lambda r: r["hausdorff_sampled_W"])
    ax.annotate(
        f"late long sample\n{late['hausdorff_sampled_W']:.5f} W",
        (late["time_W2"], late["hausdorff_sampled_W"]), xytext=(0.075, 0.071), textcoords="data",
        ha="left", va="top", arrowprops={"arrowstyle": "-", "lw": 0.7}, fontsize=6.7,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.7},
    )
    ax.set_title("(a) Halo shape error at recorded samples", loc="center", fontweight="bold")
    ax.set_xlabel("Physical time ($W^2$)")
    ax.set_ylabel("Sampled Hausdorff error (W)")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=6.6, loc="upper left")
    ax.set_ylim(0, 0.076)

    pairs = ["short", "long"]
    pair_x = {"short": 0, "long": 1}
    route_colors = {"Native OpenPhase": COLORS["native"], "OpenPhase halo": COLORS["whole"]}
    for ax, key, title, ylabel, threshold, panel in (
        (axes[0, 1], "worst_nonempty_shape_W", "Maximum nonempty sampled shape error", "Hausdorff error (W)", 0.05, "(b)"),
        (axes[1, 0], "speed_error_percent", "Area-speed error", "Relative error (%)", 5.0, "(c)"),
        (axes[1, 1], "T99_upper_error_percent", "T99 upper-bracket error", "Relative error (%)", 5.0, "(d)"),
    ):
        ax.set_title(f"{panel} {title}", loc="center", fontweight="bold")
        for route, color, offset, marker in (("Native OpenPhase", COLORS["native"], -0.11, "o"), ("OpenPhase halo", COLORS["whole"], 0.11, "s")):
            vals = [next(r[key] for r in metric_rows if r["route"] == route and r["pair"] == pair) for pair in pairs]
            ax.plot([pair_x[p] + offset for p in pairs], vals, ls="", marker=marker, ms=5.2, color=color, label=route)
        ax.axhline(threshold, color=COLORS["red"], lw=0.9, ls="--")
        ax.set_xticks([0, 1], ["Short", "Long"])
        ax.set_ylabel(ylabel)
        ax.set_xlim(-0.38, 1.38)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=0.2)
        legend_loc = "upper left" if key == "worst_nonempty_shape_W" else "center left"
        ax.legend(frameon=False, fontsize=6.7, loc=legend_loc)
    fig.suptitle("OpenPhase transfer: one cold finite-grid recipe at unchanged base step", y=0.985, fontsize=9.5, fontweight="bold")
    fig.subplots_adjust(left=0.12, right=0.985, bottom=0.105, top=0.89, wspace=0.38, hspace=0.62)
    return fig



def make_graphical_abstract(n20_metrics: list[dict[str, Any]], open_metrics: list[dict[str, Any]]) -> plt.Figure:
    n20 = {r["route"]: r for r in n20_metrics}
    n_native = next(r for r in open_metrics if r["route"] == "Native OpenPhase" and r["pair"] == "short")
    n_halo = next(r for r in open_metrics if r["route"] == "OpenPhase halo" and r["pair"] == "short")
    short_halo_shape = next(r for r in open_metrics if r["route"] == "OpenPhase halo" and r["pair"] == "short")["worst_nonempty_shape_W"]
    long_halo_shape = next(r for r in open_metrics if r["route"] == "OpenPhase halo" and r["pair"] == "long")["worst_nonempty_shape_W"]

    fig = plt.figure(figsize=(11.5, 4.3), constrained_layout=False)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.78, 1.15, 1.55])

    # Small schematic: a ternary junction with a square support halo.
    ax = fig.add_subplot(gs[0, 0])
    ax.set_aspect("equal")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.axis("off")
    ax.add_patch(Rectangle((-0.64, -0.64), 1.28, 1.28, fill=False, ec=COLORS["whole"], lw=1.5, ls="--"))
    ax.plot([-0.9, 0], [0, 0], color=COLORS["native"], lw=3)
    ax.plot([0, 0.62], [0, 0.9], color=COLORS["whole"], lw=3)
    ax.plot([0, 0.62], [0, -0.9], color=COLORS["local"], lw=3)
    ax.add_patch(Circle((0, 0), 0.07, fc=COLORS["dark"], ec="white", lw=0.6, zorder=4))
    ax.text(-0.48, 0.23, "phase 0", ha="center", va="center", fontsize=7.2, color=COLORS["native"])
    ax.text(0.41, 0.59, "phase 1", ha="center", va="center", fontsize=7.2, color=COLORS["whole"])
    ax.text(0.43, -0.59, "phase 2", ha="center", va="center", fontsize=7.2, color=COLORS["local"])
    ax.text(0, 0.87, r"$m=1+27H$", ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0, -0.88, "max-square ternary halo", ha="center", va="center", fontsize=7.8)

    ax = fig.add_subplot(gs[0, 1])
    ax.set_title("Matched high-drive B2 benchmark", loc="left", fontweight="bold", fontsize=10)
    ax.axis("off")
    ax.text(0.02, 0.78, r"$q=\sqrt{3},\ p=-50$; same grid, $dt$, warm source", transform=ax.transAxes, fontsize=8.0, color=COLORS["gray"])
    ax.text(0.02, 0.53, "Speed error", transform=ax.transAxes, fontsize=9.2, fontweight="bold")
    a = n20["Original B2"]["speed_error_percent"]
    b = n20["Corrected B2"]["speed_error_percent"]
    ax.text(0.02, 0.37, f"{a:.2f}%", transform=ax.transAxes, fontsize=17, color=COLORS["native"], fontweight="bold")
    ax.text(0.37, 0.37, "→", transform=ax.transAxes, fontsize=18, color=COLORS["gray"], ha="center")
    ax.text(0.51, 0.37, f"{b:.3f}%", transform=ax.transAxes, fontsize=17, color=COLORS["whole"], fontweight="bold")
    ax.text(0.02, 0.20, "Terminal aligned RMS", transform=ax.transAxes, fontsize=8.6, fontweight="bold")
    y0 = n20["Original B2"]["terminal_fixed92_rms_W"]
    y1 = n20["Corrected B2"]["terminal_fixed92_rms_W"]
    ax.text(0.02, 0.04, f"{y0:.3f} W  →  {y1:.3f} W", transform=ax.transAxes, fontsize=12.1, color=COLORS["dark"])

    ax = fig.add_subplot(gs[0, 2])
    ax.set_title("OpenPhase: positive transfer, visible limit", loc="left", fontweight="bold", fontsize=10)
    ax.axis("off")
    ax.text(0.02, 0.82, "Speed error", transform=ax.transAxes, fontsize=8.8, fontweight="bold")
    ax.text(0.02, 0.67, f"{n_native['speed_error_percent']:.2f}%", transform=ax.transAxes, fontsize=13, color=COLORS["native"], fontweight="bold")
    ax.text(0.25, 0.67, "→", transform=ax.transAxes, fontsize=15, color=COLORS["gray"])
    ax.text(0.37, 0.67, f"{n_halo['speed_error_percent']:.2f}%", transform=ax.transAxes, fontsize=13, color=COLORS["whole"], fontweight="bold")
    ax.text(0.02, 0.49, "Halo sampled shape errors (W)", transform=ax.transAxes, fontsize=8.8, fontweight="bold")
    ax.text(0.02, 0.35, f"Short: {short_halo_shape:.3f} W (5/5 ≤ 0.05 W)", transform=ax.transAxes, fontsize=8.9, color=COLORS["dark"])
    ax.text(0.02, 0.21, f"Long:   {long_halo_shape:.3f} W (4/5 ≤ 0.05 W)", transform=ax.transAxes, fontsize=8.9, color=COLORS["dark"])
    ax.plot([0.02, 0.95], [0.13, 0.13], transform=ax.transAxes, color=COLORS["red"], ls="--", lw=1.0)
    ax.text(0.02, 0.03, "One late OpenPhase shape sample exceeds 0.05 W", transform=ax.transAxes, fontsize=7.4, color=COLORS["red"])

    fig.subplots_adjust(left=0.025, right=0.98, bottom=0.08, top=0.82, wspace=0.24)
    fig.text(0.5, 0.98, "Spatial capillary correction for driven junctions", ha="center", va="top", fontsize=12, fontweight="bold")
    return fig



def save_all_figures(
    profiles: list[dict[str, Any]],
    n20_shape: list[dict[str, Any]],
    n20_metrics: list[dict[str, Any]],
    h24_shape: list[dict[str, Any]],
    h24_metrics: list[dict[str, Any]],
    robustness: list[dict[str, Any]],
    open_shape: list[dict[str, Any]],
    open_metrics: list[dict[str, Any]],
) -> list[Path]:
    figs = [
        ("fig01_n20_matched.pdf", make_n20_figure({}, profiles, n20_shape, n20_metrics)),
        ("fig02_h24_routes_grose.pdf", make_h24_figure(h24_shape, h24_metrics)),
        ("fig03_parameter_robustness.pdf", make_robustness_figure(robustness)),
        ("fig04_openphase_transfer.pdf", make_openphase_figure(open_shape, open_metrics)),
        ("graphical_abstract.pdf", make_graphical_abstract(n20_metrics, open_metrics)),
    ]
    produced: list[Path] = []
    for pdf_name, fig in figs:
        pdf = FIGURES / pdf_name
        png = FIGURES / pdf_name.replace(".pdf", ".png")
        dpi = 190 if pdf_name == "graphical_abstract.pdf" else 320
        save_figure(fig, pdf)
        save_figure(fig, png, dpi=dpi)
        produced.extend((pdf, png))
        plt.close(fig)
    return produced



def main():
    global FIGURES
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=FIGURES)
    args=parser.parse_args(); FIGURES=args.output_dir.resolve(); FIGURES.mkdir(parents=True,exist_ok=True)
    setup_style()
    profiles=read_export_csv('n20_terminal_fixed92_profiles.csv')
    shapes=read_export_csv('n20_shape_history.csv')
    metrics=read_export_csv('n20_route_metrics.csv')
    h24shape=read_export_csv('h24_shape_samples.csv')
    h24metrics=read_export_csv('h24_speed_t99_metrics.csv')
    robustness=read_export_csv('halo_parameter_robustness.csv')
    opshape=read_export_csv('openphase_halo_shape_samples.csv')
    opmetrics=read_export_csv('openphase_route_metrics.csv')
    produced=save_all_figures(profiles,shapes,metrics,h24shape,h24metrics,robustness,opshape,opmetrics)
    manifest={'scope':'Plot-only regeneration from frozen CSVs; no PDE updates or refitting.', 'tables':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(DATA.glob('*.csv'))},'outputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in produced}}
    (FIGURES/'figure_build_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'status':'PASS','figures':len(produced)//2,'output_dir':str(FIGURES)},indent=2))

if __name__=='__main__': main()
