#!/usr/bin/env python3

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"

GENES = [
    "TYK2",
    "STAT4",
    "FCGR2A",
    "TNFSF4",
    "IL12A",
    "PTPN22",
    "ITGAM",
    "CLEC16A",
    "IRF5",
    "HLA-DR",
]

TARGETS = [
    "JAK-STAT",
    "Complement",
    "OX40 / IL-12",
    "FcR modulation",
    "mTOR pathway",
    "Phosphatase",
]

MATRIX = np.array(
    [
        [90, 0, 0, 0, 0, 0],
        [85, 0, 0, 0, 0, 0],
        [0, 0, 0, 85, 0, 0],
        [0, 65, 80, 0, 0, 0],
        [35, 0, 80, 0, 0, 0],
        [0, 0, 0, 0, 0, 80],
        [0, 75, 0, 0, 0, 0],
        [0, 0, 0, 0, 75, 0],
        [45, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
    ],
    dtype=float,
)


def draw_stage_legend(ax):
    ax.text(0.0, 0.98, "Clinical stage", ha="left", va="top", fontsize=6.4, fontweight="bold", color="#111111")
    items = [
        ("Approved", "4", "#1f9d7a"),
        ("Phase III", "2", "#cf741c"),
        ("Phase II", "3", "#e04b3f"),
        ("Investigational", "1", "#6b7280"),
    ]
    y = 0.82
    for label, count, color in items:
        ax.add_patch(Rectangle((0.02, y - 0.035), 0.026, 0.07, facecolor=color, edgecolor="#111111", linewidth=0.35))
        ax.text(0.075, y + 0.012, label, ha="left", va="center", fontsize=5.8, color="#111111", fontweight="bold")
        ax.text(0.075, y - 0.035, f"({count})", ha="left", va="center", fontsize=5.2, color="#4b5563")
        y -= 0.18


def draw_approach_legend(ax):
    ax.text(0.52, 0.98, "Therapeutic\napproaches", ha="left", va="top", fontsize=6.4, fontweight="bold", color="#111111")
    items = [
        ("JAK-STAT", "JAK inhibitors", "#b91c1c"),
        ("Complement", "C5a inhibition", "#1f9d7a"),
        ("OX40/IL-12", "Co-stimulation & cytokine", "#2f7fb8"),
        ("FcR", "Fc receptor", "#cf741c"),
        ("mTOR", "mTOR pathway", "#4f46e5"),
        ("Phosphatase", "Tyrosine", "#0f766e"),
    ]
    y = 0.82
    for label, detail, color in items:
        ax.scatter([0.55], [y], s=24, color=color, edgecolor="#111111", linewidth=0.35, zorder=3)
        ax.text(0.59, y + 0.018, label, ha="left", va="center", fontsize=5.7, color="#111111", fontweight="bold")
        ax.text(0.59, y - 0.034, detail, ha="left", va="center", fontsize=4.9, color="#4b5563")
        y -= 0.14


def make_plot():
    cmap = LinearSegmentedColormap.from_list(
        "premium_reds",
        ["#fff7f3", "#fdd5c8", "#fca082", "#ef5a49", "#b91222"],
    )
    norm = Normalize(vmin=0, vmax=100)

    fig = plt.figure(figsize=(7.8, 4.55), dpi=600, facecolor="white")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.026, 0.60], wspace=0.13)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    lax = fig.add_subplot(gs[0, 2])

    im = ax.imshow(MATRIX, cmap=cmap, norm=norm, aspect="auto")

    ax.set_title("Therapeutic Target Landscape for SLE Susceptibility Loci", loc="left", fontsize=8.5, fontweight="bold", color="#111111", pad=8)
    ax.set_xlabel("Therapeutic pathway/target", fontsize=7.3, fontweight="bold", color="#111111", labelpad=8)
    ax.set_ylabel("SLE susceptibility gene", fontsize=7.3, fontweight="bold", color="#111111", labelpad=8)

    ax.set_xticks(np.arange(len(TARGETS)))
    ax.set_xticklabels(TARGETS, rotation=45, ha="right", fontsize=6.1, color="#334155")
    ax.set_yticks(np.arange(len(GENES)))
    ax.set_yticklabels(GENES, fontsize=6.4, color="#111111")

    ax.set_xticks(np.arange(-0.5, len(TARGETS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(GENES), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.3)
    ax.tick_params(which="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for y in range(MATRIX.shape[0]):
        for x in range(MATRIX.shape[1]):
            value = int(MATRIX[y, x])
            color = "white" if value >= 60 else "#334155"
            ax.text(x, y, str(value), ha="center", va="center", fontsize=5.8, color=color, fontweight="bold" if value else "normal")

    cb = fig.colorbar(im, cax=cax)
    cb.ax.set_title("Evidence\nscore", fontsize=5.9, color="#111111", fontweight="bold", pad=7)
    cb.ax.tick_params(labelsize=5.2, colors="#334155", length=0, pad=2)
    cb.outline.set_visible(False)

    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.axis("off")
    draw_stage_legend(lax)
    draw_approach_legend(lax)

    fig.subplots_adjust(left=0.095, right=0.985, top=0.88, bottom=0.24)
    return fig


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = make_plot()
    fig.savefig(OUT_DIR / "Figure_7.png", dpi=600, facecolor="white")
    fig.savefig(OUT_DIR / "Figure_7.tiff", dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("Saved therapeutic target landscape at 600 dpi.")


if __name__ == "__main__":
    main()
