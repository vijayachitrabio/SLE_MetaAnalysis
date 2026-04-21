from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "SLE_Publication_Package" / "Supplementary" / "Figures"


def add_box(ax, x, y, w, h, title, body, fill, edge, title_color):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        linewidth=2.0,
        edgecolor=edge,
        facecolor=fill,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h * 0.68,
        title,
        ha="center",
        va="center",
        fontsize=17,
        weight="bold",
        color=title_color,
        family="DejaVu Serif",
    )
    ax.text(
        x + w / 2,
        y + h * 0.34,
        body,
        ha="center",
        va="center",
        fontsize=11.5,
        color="#2d3648",
        linespacing=1.45,
        family="DejaVu Sans",
    )


def add_arrow(ax, x1, y1, x2, y2):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=18,
            linewidth=2.2,
            color="#5c6b82",
            shrinkA=8,
            shrinkB=8,
        )
    )


def build_figure():
    fig = plt.figure(figsize=(16, 9), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    bg = "#fcfcfd"
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)

    ax.text(
        0.5,
        0.94,
        "Supplementary Workflow of the European SLE GWAS Meta-analysis",
        ha="center",
        va="center",
        fontsize=24,
        weight="bold",
        color="#172235",
        family="DejaVu Serif",
    )
    ax.text(
        0.5,
        0.90,
        "From cohort harmonization to replication, causal prioritization, and therapeutic interpretation",
        ha="center",
        va="center",
        fontsize=12.5,
        color="#64748b",
        family="DejaVu Sans",
    )

    # top input box
    add_box(
        ax,
        0.18,
        0.75,
        0.64,
        0.11,
        "Input Cohorts and Study Design",
        "Bentham 2015 and FinnGen R12 summary statistics for discovery\nSpanish cohort used for independent North-to-South replication",
        fill="#edf4fb",
        edge="#9cb9d9",
        title_color="#2f5d87",
    )

    # five-step workflow
    xs = [0.05, 0.24, 0.43, 0.62, 0.81]
    colors = [
        ("#dff2ed", "#78b8a8", "#2f7f73"),
        ("#e6f0fb", "#8fb0d7", "#355f93"),
        ("#f7eed8", "#d1b067", "#8b6822"),
        ("#f6e5ea", "#cb8ea1", "#9d4f6a"),
        ("#ebe8f8", "#a59bda", "#6557b8"),
    ]
    titles = [
        "1. Discovery",
        "2. Refinement",
        "3. Functional Mapping",
        "4. Validation",
        "5. Final Synthesis",
    ]
    bodies = [
        "IVW meta-analysis\nGWAS-significant locus detection\nManhattan and QQ visualization",
        "QC summary tables\nannotation and gene mapping\nLD pruning and novelty checks",
        "Pathway enrichment\neQTL integration\nLAVA and COLOC prioritization",
        "Spanish replication\nheterogeneity and sensitivity analyses\nregional locus review",
        "PheWAS and pleiotropy profiling\ntherapeutic target mapping\nhigh-confidence locus reporting",
    ]

    y = 0.42
    w = 0.14
    h = 0.20
    for x, (fill, edge, title_color), title, body in zip(xs, colors, titles, bodies):
        add_box(ax, x, y, w, h, title, body, fill, edge, title_color)

    for i in range(len(xs) - 1):
        add_arrow(ax, xs[i] + w, y + h / 2, xs[i + 1], y + h / 2)

    # output band
    add_box(
        ax,
        0.13,
        0.15,
        0.74,
        0.12,
        "Key Outputs",
        "47 independent loci | 25 putatively novel signals | 15 high-confidence targets | "
        "replication evidence | pathway, eQTL, pleiotropy, and therapeutic interpretation",
        fill="#f3f5f8",
        edge="#c7d0dc",
        title_color="#41516a",
    )

    add_arrow(ax, 0.50, 0.75, 0.50, 0.62)
    add_arrow(ax, 0.50, 0.42, 0.50, 0.27)

    return fig


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = build_figure()
    fig.savefig(
        OUT_DIR / "Supplementary_Pipeline_Workflow_Figure.png",
        dpi=300,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
    )
    fig.savefig(
        OUT_DIR / "Supplementary_Pipeline_Workflow_Figure.pdf",
        dpi=300,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
