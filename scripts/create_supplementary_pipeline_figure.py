from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "SLE_Publication_Package" / "Supplementary" / "Figures"


def add_box(ax, x, y, w, h, title, body, fill, edge, title_color, body_size=10.4):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.022",
        linewidth=1.8,
        edgecolor=edge,
        facecolor=fill,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h * 0.70,
        title,
        ha="center",
        va="center",
        fontsize=13.8,
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
        fontsize=body_size,
        color="#334155",
        linespacing=1.35,
        family="DejaVu Sans",
    )


def add_arrow(ax, x1, y1, x2, y2, color="#66768f"):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=16,
            linewidth=1.9,
            color=color,
            shrinkA=6,
            shrinkB=6,
        )
    )


def build_figure():
    fig = plt.figure(figsize=(18, 7.6), dpi=300)
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
        fontsize=22,
        weight="bold",
        color="#1b263b",
        family="DejaVu Serif",
    )
    ax.text(
        0.5,
        0.905,
        "From cohort harmonization to replication, causal prioritization, and therapeutic interpretation",
        ha="center",
        va="center",
        fontsize=11.5,
        color="#64748b",
        family="DejaVu Sans",
    )

    add_box(
        ax,
        0.17,
        0.75,
        0.66,
        0.105,
        "Input Cohorts and Study Design",
        "Bentham 2015 and FinnGen R12 used for discovery\nSpanish cohort used for independent North-to-South replication",
        fill="#eef4fb",
        edge="#a9bfdc",
        title_color="#305a86",
        body_size=10.6,
    )

    titles = [
        "1. Discovery",
        "2. Refinement",
        "3. Functional Mapping",
        "4. Validation",
        "5. Final Synthesis",
    ]
    bodies = [
        "IVW meta-analysis\nGenome-wide signal discovery\nManhattan and QQ plots",
        "QC summary tables\nAnnotation and gene mapping\nLD pruning and novelty review",
        "Pathway enrichment\neQTL integration\nLAVA and COLOC support",
        "Spanish replication\nHeterogeneity checks\nRegional locus review",
        "PheWAS profiling\nTherapeutic mapping\nHigh-confidence locus prioritization",
    ]
    colors = [
        ("#e6f5ef", "#84bea9", "#2f7c6d"),
        ("#ebf3fb", "#98b7dc", "#335f93"),
        ("#faf0da", "#d5b672", "#8c6a22"),
        ("#f7e8ec", "#d29bad", "#9b526b"),
        ("#efeaf9", "#b1a5de", "#6657b5"),
    ]

    x0 = 0.055
    gap = 0.02
    w = 0.156
    h = 0.17
    y = 0.43
    xs = [x0 + i * (w + gap) for i in range(5)]

    for x, title, body, (fill, edge, title_color) in zip(xs, titles, bodies, colors):
        add_box(ax, x, y, w, h, title, body, fill, edge, title_color)

    for i in range(4):
        add_arrow(ax, xs[i] + w, y + h / 2, xs[i + 1], y + h / 2)

    add_arrow(ax, 0.5, 0.75, 0.5, 0.60)
    add_arrow(ax, 0.5, 0.43, 0.5, 0.28)

    add_box(
        ax,
        0.15,
        0.14,
        0.70,
        0.10,
        "Key Outputs",
        "47 independent loci | 25 putatively novel signals | 15 high-confidence targets | "
        "replication, pathway, eQTL, pleiotropy, and therapeutic evidence",
        fill="#f4f6f9",
        edge="#c9d2de",
        title_color="#415268",
        body_size=10.2,
    )

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
