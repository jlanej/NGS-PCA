#!/usr/bin/env python3
"""Generate a publication-quality NGS-PCA pipeline workflow diagram.

This script produces a vector-quality (PDF + PNG) figure suitable for
supplementary material or methods sections.  It is designed to be run
inside the companion Docker image (docs/Dockerfile.diagram) so that
font rendering and library versions are fully reproducible.

Usage:
    python docs/generate_workflow_diagram.py [output_directory]

If *output_directory* is omitted the images are written to ``docs/``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---------------------------------------------------------------------------
# Color palette – muted, color-blind-friendly (adapted from Paul Tol)
# ---------------------------------------------------------------------------
CLR_INPUT = "#4477AA"       # blue  – input data
CLR_TOOL = "#228833"        # green – external tool (mosdepth)
CLR_CORE = "#EE6677"        # red   – NGS-PCA core steps
CLR_OUTPUT = "#CCBB44"      # gold  – output files
CLR_NEXT = "#AA3377"        # purple – downstream / future
CLR_ARROW = "#555555"       # grey  – arrows
CLR_BG = "#FFFFFF"          # white background
CLR_BORDER = "#333333"      # dark border
CLR_TEXT = "#222222"         # text color
CLR_SUBTITLE = "#555555"    # subtitle / annotation color
CLR_ANNOT_BG = "#F0F0F0"    # annotation box background


def _rounded_box(ax, x, y, w, h, text, subtext, facecolor, *,
                 fontsize=10, sub_fontsize=7.5, text_color="white",
                 corner_radius=0.15, linewidth=1.0, edgecolor=CLR_BORDER):
    """Draw a rounded-rectangle box with a title and optional subtitle."""
    box = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle=f"round,pad=0,rounding_size={corner_radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        zorder=3,
    )
    ax.add_patch(box)
    if subtext:
        ax.text(x, y + 0.12, text, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", color=text_color,
                zorder=4)
        ax.text(x, y - 0.15, subtext, ha="center", va="center",
                fontsize=sub_fontsize, color=text_color, zorder=4,
                style="italic")
    else:
        ax.text(x, y, text, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", color=text_color,
                zorder=4)
    return box


def _arrow(ax, x1, y1, x2, y2, *, color=CLR_ARROW, linewidth=1.5,
           connectionstyle="arc3,rad=0", shrinkA=8, shrinkB=8):
    """Draw an arrow between two points."""
    arrow = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="->,head_length=6,head_width=4",
        color=color,
        linewidth=linewidth,
        connectionstyle=connectionstyle,
        shrinkA=shrinkA, shrinkB=shrinkB,
        zorder=2,
    )
    ax.add_patch(arrow)
    return arrow


def _stage_bracket(ax, x_margin, y_top, y_bot, label, *, color=CLR_SUBTITLE):
    """Draw a thin bracket + label on the left margin for stage grouping."""
    tick = 0.06
    ax.plot([x_margin, x_margin], [y_top, y_bot], color=color,
            linewidth=1.2, solid_capstyle="round", zorder=1)
    ax.plot([x_margin, x_margin + tick], [y_top, y_top], color=color,
            linewidth=1.2, solid_capstyle="round", zorder=1)
    ax.plot([x_margin, x_margin + tick], [y_bot, y_bot], color=color,
            linewidth=1.2, solid_capstyle="round", zorder=1)
    ax.text(x_margin - 0.06, (y_top + y_bot) / 2, label,
            ha="right", va="center", fontsize=7.5, color=color,
            rotation=90, fontweight="bold", zorder=1)


def generate_diagram(output_dir: str | Path) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9.5, 13))
    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-0.5, 11.5)
    ax.axis("off")
    fig.patch.set_facecolor(CLR_BG)
    ax.set_facecolor(CLR_BG)

    # ── Title ──────────────────────────────────────────────────────────────
    ax.text(0, 11.15, "NGS-PCA Pipeline",
            ha="center", va="center", fontsize=17, fontweight="bold",
            color=CLR_TEXT, zorder=5)
    ax.text(0, 10.85,
            "Principal component analysis of sequencing coverage\n"
            "via randomized singular value decomposition",
            ha="center", va="center", fontsize=9, color=CLR_SUBTITLE,
            zorder=5, linespacing=1.4)

    # ── Y coordinates (top to bottom) ─────────────────────────────────────
    y_input = 10.1
    y_mosdepth = 9.0
    y_bed = 7.95
    y_region = 7.05
    y_norm = 6.0
    y_svd = 4.9
    y_out = 3.6
    y_next = 1.9

    box_w = 2.8
    box_h = 0.62

    # bracket x position – well to the left of the boxes
    bracket_x = -(box_w / 2) - 0.35

    # ── 1. Input ───────────────────────────────────────────────────────────
    _rounded_box(ax, 0, y_input, box_w, box_h,
                 "BAM / CRAM Files",
                 "aligned sequencing reads",
                 CLR_INPUT)

    # ── 2. mosdepth ───────────────────────────────────────────────────────
    _arrow(ax, 0, y_input - box_h / 2, 0, y_mosdepth + box_h / 2)
    _rounded_box(ax, 0, y_mosdepth, box_w, box_h,
                 "mosdepth",
                 "coverage in fixed-width bins (default 1 kb)",
                 CLR_TOOL)
    _stage_bracket(ax, bracket_x,
                   y_input + box_h / 2 - 0.05,
                   y_mosdepth - box_h / 2 + 0.05,
                   "Coverage\nComputation")

    # Reference genome annotation (left side of mosdepth)
    ax.annotate(
        "Reference genome\n(FASTA)",
        xy=(-(box_w / 2) - 0.05, y_mosdepth),
        xytext=(-(box_w / 2) - 0.85, y_mosdepth + 0.45),
        fontsize=6.5, color=CLR_SUBTITLE, ha="center", va="center",
        arrowprops=dict(arrowstyle="->", color=CLR_SUBTITLE,
                        linewidth=0.8, connectionstyle="arc3,rad=0.15"),
        zorder=4,
        bbox=dict(boxstyle="round,pad=0.3", facecolor=CLR_ANNOT_BG,
                  edgecolor=CLR_SUBTITLE, linewidth=0.6),
    )

    # ── 3. BED output ─────────────────────────────────────────────────────
    _arrow(ax, 0, y_mosdepth - box_h / 2, 0, y_bed + box_h / 2)
    _rounded_box(ax, 0, y_bed, box_w, box_h,
                 "*.regions.bed.gz",
                 "per-sample depth \u00d7 genomic bin",
                 CLR_INPUT, text_color="white")

    # ── 4. Region selection ───────────────────────────────────────────────
    _arrow(ax, 0, y_bed - box_h / 2, 0, y_region + box_h / 2)
    _rounded_box(ax, 0, y_region, box_w, box_h,
                 "Region Selection",
                 "autosomal bins; exclude SV/mappability blacklist",
                 CLR_CORE)

    # Exclusion BED annotation (right side)
    ax.annotate(
        "Exclusion BED\n(SV blacklist,\nlow mappability,\nseg. dups)",
        xy=(box_w / 2 + 0.05, y_region),
        xytext=(box_w / 2 + 0.85, y_region + 0.40),
        fontsize=6.5, color=CLR_SUBTITLE, ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color=CLR_SUBTITLE,
                        linewidth=0.8, connectionstyle="arc3,rad=-0.15"),
        zorder=4,
        bbox=dict(boxstyle="round,pad=0.3", facecolor=CLR_ANNOT_BG,
                  edgecolor=CLR_SUBTITLE, linewidth=0.6),
    )

    # ── 5. Normalization ──────────────────────────────────────────────────
    _arrow(ax, 0, y_region - box_h / 2, 0, y_norm + box_h / 2)
    _rounded_box(ax, 0, y_norm, box_w, box_h,
                 "Normalization",
                 "log\u2082 fold-change \u2192 median-centre across samples",
                 CLR_CORE)

    # Normalization detail annotation (left side)
    ax.annotate(
        "Per-sample:\n  log\u2082(depth / median)\n"
        "Per-bin:\n  subtract row median",
        xy=(-(box_w / 2) - 0.05, y_norm),
        xytext=(-(box_w / 2) - 0.85, y_norm - 0.40),
        fontsize=6.5, color=CLR_SUBTITLE, ha="right", va="center",
        arrowprops=dict(arrowstyle="->", color=CLR_SUBTITLE,
                        linewidth=0.8, connectionstyle="arc3,rad=-0.15"),
        zorder=4,
        bbox=dict(boxstyle="round,pad=0.3", facecolor=CLR_ANNOT_BG,
                  edgecolor=CLR_SUBTITLE, linewidth=0.6),
    )

    # ── 6. Randomized SVD ─────────────────────────────────────────────────
    _arrow(ax, 0, y_norm - box_h / 2, 0, y_svd + box_h / 2)
    _rounded_box(ax, 0, y_svd, box_w, box_h,
                 "Randomized SVD",
                 "Halko et al. 2011  \u2022  power iterations",
                 CLR_CORE)

    _stage_bracket(ax, bracket_x,
                   y_region + box_h / 2 - 0.05,
                   y_svd - box_h / 2 + 0.05,
                   "NGS-PCA")

    # SVD detail annotation (right side)
    ax.annotate(
        "Random projection \u2192\n"
        "QR power iterations \u2192\n"
        "Truncated SVD on\n"
        "reduced matrix  (A \u2248 U\u03A3V\u1d40)",
        xy=(box_w / 2 + 0.05, y_svd),
        xytext=(box_w / 2 + 0.85, y_svd - 0.40),
        fontsize=6.5, color=CLR_SUBTITLE, ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color=CLR_SUBTITLE,
                        linewidth=0.8, connectionstyle="arc3,rad=0.15"),
        zorder=4,
        bbox=dict(boxstyle="round,pad=0.3", facecolor=CLR_ANNOT_BG,
                  edgecolor=CLR_SUBTITLE, linewidth=0.6),
    )

    # ── 7. Output ─────────────────────────────────────────────────────────
    _arrow(ax, 0, y_svd - box_h / 2, 0, y_out + box_h / 2 + 0.3)

    out_box_h = 1.1
    out_box = FancyBboxPatch(
        (-box_w / 2, y_out - out_box_h / 2), box_w, out_box_h,
        boxstyle="round,pad=0,rounding_size=0.15",
        facecolor=CLR_OUTPUT,
        edgecolor=CLR_BORDER,
        linewidth=1.0,
        zorder=3,
    )
    ax.add_patch(out_box)

    ax.text(0, y_out + 0.35, "Output Files",
            ha="center", va="center", fontsize=10, fontweight="bold",
            color=CLR_TEXT, zorder=4)

    outputs = [
        ("svd.pcs.txt", "sample \u00d7 PC scores"),
        ("svd.loadings.txt", "bin \u00d7 PC loadings"),
        ("svd.singularvalues.txt", "variance per PC"),
        ("svd.bins.txt", "retained genomic bins"),
        ("svd.samples.txt", "sample identifiers"),
    ]
    col_file = -0.55
    col_desc = 0.65
    for i, (fname, desc) in enumerate(outputs):
        row_y = y_out + 0.12 - i * 0.18
        ax.text(col_file, row_y, fname, ha="right", va="center",
                fontsize=6.5, fontfamily="monospace", fontweight="bold",
                color=CLR_TEXT, zorder=4)
        ax.text(col_desc, row_y, desc, ha="left", va="center",
                fontsize=6.5, color=CLR_TEXT, zorder=4)

    # ── 8. Downstream / Next Steps ────────────────────────────────────────
    _arrow(ax, 0, y_out - out_box_h / 2, 0, y_next + 0.7)

    next_box_h = 1.35
    next_box = FancyBboxPatch(
        (-box_w / 2, y_next - next_box_h / 2), box_w, next_box_h,
        boxstyle="round,pad=0,rounding_size=0.15",
        facecolor=CLR_NEXT,
        edgecolor=CLR_BORDER,
        linewidth=1.0,
        linestyle="--",
        zorder=3,
    )
    ax.add_patch(next_box)

    ax.text(0, y_next + 0.45, "Potential Downstream Analyses",
            ha="center", va="center", fontsize=9.5, fontweight="bold",
            color="white", zorder=4)

    next_items = [
        "\u2022 Sequencing batch-effect detection & mitigation",
        "\u2022 Sample outlier identification",
    ]
    for i, item in enumerate(next_items):
        ax.text(0, y_next + 0.16 - i * 0.18, item,
                ha="center", va="center", fontsize=7, color="white",
                zorder=4)

    # ── Legend ─────────────────────────────────────────────────────────────
    legend_handles = [
        mpatches.Patch(facecolor=CLR_INPUT, edgecolor=CLR_BORDER,
                       label="Input data"),
        mpatches.Patch(facecolor=CLR_TOOL, edgecolor=CLR_BORDER,
                       label="External tool"),
        mpatches.Patch(facecolor=CLR_CORE, edgecolor=CLR_BORDER,
                       label="NGS-PCA core"),
        mpatches.Patch(facecolor=CLR_OUTPUT, edgecolor=CLR_BORDER,
                       label="Output files"),
        mpatches.Patch(facecolor=CLR_NEXT, edgecolor=CLR_BORDER,
                       label="Downstream"),
    ]
    ax.legend(handles=legend_handles, loc="lower center",
              ncol=5, fontsize=7, frameon=True, framealpha=0.9,
              edgecolor="#CCCCCC", fancybox=True,
              bbox_to_anchor=(0.5, -0.01))

    # ── Citation footer ───────────────────────────────────────────────────
    ax.text(0, -0.3,
            "Halko, Martinsson & Tropp (2011) SIAM Rev 53:217\u2013288  \u2022  "
            "Pedersen & Quinlan (2018) Bioinformatics 34:867\u2013868",
            ha="center", va="center", fontsize=5.5, color="#999999",
            zorder=5)

    # ── Save ──────────────────────────────────────────────────────────────
    for fmt in ("png", "pdf"):
        out_path = output_dir / f"pipeline_workflow.{fmt}"
        fig.savefig(
            out_path,
            dpi=300,
            bbox_inches="tight",
            facecolor=CLR_BG,
            pad_inches=0.25,
        )
        print(f"Saved {out_path}")

    plt.close(fig)


if __name__ == "__main__":
    dest = sys.argv[1] if len(sys.argv) > 1 else str(
        Path(__file__).resolve().parent)
    generate_diagram(dest)
