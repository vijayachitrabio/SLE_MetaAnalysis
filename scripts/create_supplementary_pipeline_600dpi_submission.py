from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "SLE_submission"


def add_box(
    ax,
    x,
    y,
    w,
    h,
    title,
    body,
    body_size=8.6,
    title_size=10.8,
    title_y=0.69,
    body_y=0.33,
):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.006,rounding_size=0.010",
        linewidth=1.25,
        edgecolor="#111111",
        facecolor="#ffffff",
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h * title_y,
        title,
        ha="center",
        va="center",
        fontsize=title_size,
        weight="bold",
        color="#111111",
        family="DejaVu Sans",
        linespacing=1.05,
    )
    ax.text(
        x + w / 2,
        y + h * body_y,
        body,
        ha="center",
        va="center",
        fontsize=body_size,
        color="#222222",
        linespacing=1.25,
        family="DejaVu Sans",
    )


def add_arrow(ax, x1, y1, x2, y2, scale=12):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=scale,
            linewidth=1.20,
            color="#111111",
            shrinkA=3.5,
            shrinkB=3.5,
        )
    )


def build_figure():
    fig = plt.figure(figsize=(9.2, 2.55), dpi=600)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("#ffffff")
    ax.set_facecolor("#ffffff")

    add_box(
        ax,
        0.16,
        0.78,
        0.68,
        0.145,
        "Data integration",
        "Discovery: Bentham 2015 + FinnGen M13_SLE  |  Replication: Julià Spanish cohort",
        body_size=8.3,
    )

    titles = [
        "1. Discovery",
        "2. Locus\nrefinement",
        "3. Functional\nprioritization",
        "4. Validation",
        "5. Cross-trait /\ntranslation",
    ]
    bodies = [
        "IVW meta-analysis\nSignal discovery",
        "LD pruning\nIndependent loci\nLAVA local architecture",
        "COLOC with GTEx v10 eQTLs\nPathway enrichment\nImmune chromatin-state overlap",
        "Spanish replication\nProxy-SNP evaluation\nSensitivity analysis",
        "Bivariate LAVA\nTherapeutic annotation",
    ]

    x0 = 0.045
    gap = 0.020
    w = 0.172
    h = 0.285
    y = 0.350
    xs = [x0 + i * (w + gap) for i in range(5)]

    for x, title, body in zip(xs, titles, bodies):
        add_box(
            ax,
            x,
            y,
            w,
            h,
            title,
            body,
            body_size=7.6,
            title_size=9.4,
            title_y=0.79,
            body_y=0.36,
        )

    for i in range(4):
        add_arrow(ax, xs[i] + w, y + h / 2, xs[i + 1], y + h / 2)

    add_arrow(ax, 0.5, 0.78, 0.5, y + h, scale=12)
    add_arrow(ax, 0.5, y, 0.5, 0.225, scale=12)

    add_box(
        ax,
        0.16,
        0.075,
        0.68,
        0.145,
        "Prioritized output",
        "47 genome-wide significant loci  |  9 prioritized candidate regions  |  shared autoimmune architecture  |  immune regulatory pathways",
        body_size=7.9,
    )

    return fig


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = build_figure()
    for suffix, kwargs in {
        ".png": {"dpi": 600},
        ".tiff": {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}},
        ".pdf": {},
    }.items():
        fig.savefig(
            OUT_DIR / f"Supplementary_Figure_1{suffix}",
            bbox_inches="tight",
            pad_inches=0.055,
            facecolor="white",
            **kwargs,
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
