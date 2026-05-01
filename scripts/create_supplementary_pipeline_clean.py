# scripts/create_supplementary_pipeline_clean.py
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "SLE_Publication_Package" / "Supplementary" / "Figures"

def add_box(ax, x, y, w, h, title, body, body_size=9.2):
    # Pure white box with clean black outline
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.005,rounding_size=0.01",
                           linewidth=1.5, edgecolor='#000000', facecolor='#ffffff')
    ax.add_patch(patch)
    # Bold black title
    ax.text(x + w/2, y + h * 0.72, title, ha='center', va='center', fontsize=12, 
            weight='bold', color='#000000', family='DejaVu Sans')
    # Plain black body text
    ax.text(x + w/2, y + h * 0.32, body, ha='center', va='center', fontsize=body_size, 
            color='#000000', family='DejaVu Sans', linespacing=1.4)

def add_arrow(ax, x1, y1, x2, y2):
    # Simple black arrow
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=14,
                                linewidth=1.5, color='#000000', shrinkA=4, shrinkB=4))

def build_minimalist_figure():
    fig = plt.figure(figsize=(15, 4.2), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    fig.patch.set_facecolor('#ffffff')

    # 1. Input Section
    add_box(ax, 0.15, 0.78, 0.70, 0.15, "Data Integration", 
            "Discovery: Bentham (2015) & FinnGen R12  |  Replication: Spanish Cohort")

    # 2. Main Workflow Steps (Minimalist)
    titles = ["1. Discovery", "2. Annotation", "3. Functional", "4. Validation", "5. Translation"]
    bodies = [
        "IVW Meta-analysis\nSignal Discovery",
        "Novelty Audit\nLD-Pruning",
        "LAVA / COLOC\neQTL Mapping",
        "Replication\nHeterogeneity",
        "PheWAS / Drug\nPrioritization"
    ]

    x0, gap, w, h, y = 0.04, 0.02, 0.172, 0.30, 0.36
    xs = [x0 + i * (w + gap) for i in range(5)]

    for x, title, body in zip(xs, titles, bodies):
        add_box(ax, x, y, w, h, title, body)

    for i in range(4):
        add_arrow(ax, xs[i] + w, y + h/2, xs[i+1], y + h/2)

    add_arrow(ax, 0.5, 0.78, 0.5, 0.66)
    add_arrow(ax, 0.5, 0.36, 0.5, 0.21)

    # 3. Final Summary Section
    add_box(ax, 0.15, 0.06, 0.70, 0.15, "Validated Output", 
            "47 Independent Loci  |  25 Putatively Novel Signals  |  15 High-Confidence Therapeutic Targets")

    return fig

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = build_minimalist_figure()
    fig.savefig(OUT_DIR / "Supplementary_Pipeline_Workflow_Clean.png", bbox_inches='tight', pad_inches=0.1)
    fig.savefig(OUT_DIR / "Supplementary_Pipeline_Workflow_Clean.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

if __name__ == "__main__":
    main()
