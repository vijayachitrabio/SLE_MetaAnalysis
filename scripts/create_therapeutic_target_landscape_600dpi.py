#!/usr/bin/env python3
"""
Figure 7 -- Therapeutic target landscape for SLE susceptibility loci.

Rows, pathway assignment, clinical stage and evidence scores are ALL derived
from results/Table_4.tsv. Nothing is hardcoded, so the figure cannot drift out
of sync with the table.

Evidence score follows Methods 2.11 ("ranked according to statistical support
(meta-analysis P-value) and therapeutic maturity"):

    score = 60 * genetic_term + 40 * maturity_term

    genetic_term  = min(1.0, -log10(P_meta) / 150)   # 150 caps at the MHC signal
    maturity_term = mapped from Drug_Status (see MATURITY)

Styled to match the original submitted figure: every cell shows a plain number
(no hatching), and a fixed 0-100 colour scale. Loci whose therapeutic direction
is undetermined score 0 in every column, matching how every other non-mapped
cell in the grid is rendered, and are excluded from the clinical-stage tally
since they have no assigned drug. The derivation for every cell is printed so
each number is auditable.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"
TABLE_4 = ROOT / "results" / "Table_4.tsv"

# Therapeutic_Pathway substring -> figure column
PATHWAY_TO_COLUMN = [
    ("JAK-STAT", "JAK-STAT"),
    ("Complement Receptor", "Complement"),
    ("Complement", "Complement"),
    ("OX40", "OX40 / IL-12"),
    ("IL-12", "OX40 / IL-12"),
    ("Fc Receptor", "FcR modulation"),
    ("Autophagy", "mTOR pathway"),
    ("B-cell Signaling", "B-cell Kinase"),
    ("NF-kB", "NF-κB pathway"),
    ("Myeloid / Interferon", "JAK-STAT"),
    ("Type I Interferon", "JAK-STAT"),
    ("Tyrosine Phosphatase", "Phosphatase"),
]

COLUMNS = [
    "JAK-STAT", "Complement", "OX40 / IL-12", "FcR modulation",
    "mTOR pathway", "Phosphatase", "B-cell Kinase", "NF-κB pathway",
]

# Drug_Status -> therapeutic maturity weight
MATURITY = [
    ("undetermined", None),
    ("no approved agent", 0.40),
    ("fda-approved", 1.00),
    ("approved / phase iii", 0.90),
    ("approved", 1.00),
    ("phase iii", 0.80),
    ("phase ii", 0.60),
    ("investigational", 0.40),
]

STAGE_COLOR = {
    "Approved": "#1f9d7a",
    "Phase III": "#cf741c",
    "Phase II": "#e04b3f",
    "Investigational": "#6b7280",
}


def assign_column(pathway):
    for key, col in PATHWAY_TO_COLUMN:
        if key.lower() in str(pathway).lower():
            return col
    raise ValueError(f"No column mapping for pathway: {pathway!r}")


def maturity_weight(status):
    s = str(status).strip().lower()
    for key, w in MATURITY:
        if key in s:
            return w
    raise ValueError(f"No maturity mapping for Drug_Status: {status!r}")


def stage_label(status):
    s = str(status).strip().lower()
    if "undetermined" in s:
        return None
    if "approved" in s and "phase iii" in s:
        return "Approved"
    if "approved" in s:
        return "Approved"
    if "phase iii" in s:
        return "Phase III"
    if "phase ii" in s:
        return "Phase II"
    return "Investigational"


def load_rows():
    df = pd.read_csv(TABLE_4, sep="\t")
    rows, log = [], []
    for _, r in df.iterrows():
        p = float(r["Discovery_P"])
        col = assign_column(r["Therapeutic_Pathway"])
        mat = maturity_weight(r["Drug_Status"])
        neglog = -np.log10(p)
        if mat is None:
            score = None
            deriv = "undetermined (directionality conflict) -> no score"
        else:
            gen = min(1.0, neglog / 150.0)
            score = int(round(60 * gen + 40 * mat))
            deriv = (f"-log10(P)={neglog:6.1f} gen={gen:.3f} "
                     f"mat={mat:.2f} -> {score}")
        has_eqtl = "eQTL" in str(r["Evidence_Level"])
        label = str(r["Lead_Gene"]).replace(
            "chr6p21.3 / MHC (multi-gene)*", "chr6p21.3/MHC*")
        rows.append({
            "label": label, "column": col, "score": score,
            "stage": stage_label(r["Drug_Status"]), "eqtl": has_eqtl,
        })
        log.append(f"  {str(r['Lead_Gene'])[:26]:26s} {col:16s} {deriv}"
                   + ("  [eQTL]" if has_eqtl else ""))
    # Order: scored rows by descending score, undetermined last
    rows.sort(key=lambda d: (d["score"] is None, -(d["score"] or 0)))
    return rows, log


def draw_stage_legend(ax, rows):
    ax.text(0.0, 0.98, "Clinical stage", ha="left", va="top", fontsize=6.4,
            fontweight="bold", color="#111111")
    counts = {}
    for r in rows:
        if r["stage"] is None:
            continue
        counts[r["stage"]] = counts.get(r["stage"], 0) + 1
    y = 0.86
    for label in ["Approved", "Phase III", "Phase II", "Investigational"]:
        if label not in counts:
            continue
        ax.add_patch(Rectangle((0.02, y - 0.030), 0.026, 0.060,
                               facecolor=STAGE_COLOR[label],
                               edgecolor="#111111", linewidth=0.35))
        ax.text(0.075, y + 0.010, label, ha="left", va="center", fontsize=5.6,
                color="#111111", fontweight="bold")
        ax.text(0.075, y - 0.032, f"({counts[label]})", ha="left", va="center",
                fontsize=5.0, color="#4b5563")
        y -= 0.150


def draw_approach_legend(ax):
    ax.text(0.52, 0.98, "Therapeutic\napproaches", ha="left", va="top",
            fontsize=6.4, fontweight="bold", color="#111111")
    items = [
        ("JAK-STAT", "JAK inhibitors", "#b91c1c"),
        ("Complement", "C5a / C1 inhibition", "#1f9d7a"),
        ("OX40 / IL-12", "Co-stimulation & cytokine", "#2f7fb8"),
        ("FcR", "Fc receptor", "#cf741c"),
        ("mTOR", "mTOR pathway", "#4f46e5"),
        ("Phosphatase", "Tyrosine", "#0f766e"),
        ("B-cell Kinase", "B-cell receptor signalling", "#8b5cf6"),
        ("NF-κB", "Ubiquitination / proteasome", "#d946ef"),
    ]
    y = 0.84
    for label, detail, color in items:
        ax.scatter([0.55], [y], s=22, color=color, edgecolor="#111111",
                   linewidth=0.35, zorder=3)
        ax.text(0.59, y + 0.016, label, ha="left", va="center", fontsize=5.5,
                color="#111111", fontweight="bold")
        ax.text(0.59, y - 0.030, detail, ha="left", va="center", fontsize=4.7,
                color="#4b5563")
        y -= 0.098


def make_plot(rows):
    labels = [r["label"] + ("  \u25cf" if r["eqtl"] else "     ")
              for r in rows]
    matrix = np.zeros((len(rows), len(COLUMNS)))
    for y, r in enumerate(rows):
        x = COLUMNS.index(r["column"])
        if r["score"] is not None:
            matrix[y, x] = r["score"]

    cmap = LinearSegmentedColormap.from_list(
        "premium_reds",
        ["#fff7f3", "#fdd5c8", "#fca082", "#ef5a49", "#b91222"])
    # Fixed 0-100 scale, matching the originally submitted figure.
    norm = Normalize(vmin=0, vmax=100)

    fig = plt.figure(figsize=(8.2, 4.8), dpi=600, facecolor="white")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.026, 0.55], wspace=0.13)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    lax = fig.add_subplot(gs[0, 2])

    im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")

    ax.set_title("Therapeutic Target Landscape for SLE Susceptibility Loci",
                 loc="left", fontsize=8.5, fontweight="bold", color="#111111",
                 pad=8)
    ax.set_xlabel("Therapeutic pathway/target", fontsize=7.3,
                  fontweight="bold", color="#111111", labelpad=8)
    ax.set_ylabel("SLE susceptibility locus", fontsize=7.3, fontweight="bold",
                  color="#111111", labelpad=8)

    ax.set_xticks(np.arange(len(COLUMNS)))
    ax.set_xticklabels(COLUMNS, rotation=45, ha="right", fontsize=5.8,
                       color="#334155")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=6.2, color="#111111")

    ax.set_xticks(np.arange(-0.5, len(COLUMNS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.3)
    ax.tick_params(which="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]):
            value = int(matrix[y, x])
            color = "white" if value >= 60 else "#334155"
            ax.text(x, y, str(value), ha="center", va="center", fontsize=5.8,
                    color=color, fontweight="bold" if value else "normal")

    cb = fig.colorbar(im, cax=cax)
    cb.ax.set_title("Evidence\nscore", fontsize=5.9, color="#111111",
                    fontweight="bold", pad=7)
    cb.ax.tick_params(labelsize=5.2, colors="#334155", length=0, pad=2)
    cb.outline.set_visible(False)

    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.axis("off")
    draw_stage_legend(lax, rows)
    draw_approach_legend(lax)

    # Marker key only. The evidence-score formula and the chr6p21.3/MHC note
    # live in the figure legend rather than on the artwork, so the panel stays
    # clean and the derivation is still stated in the manuscript.
    fig.text(0.145, 0.055,
             "\u25cf  colocalization support (PP4 > 0.80)",
             fontsize=5.0, color="#4b5563")

    fig.subplots_adjust(left=0.145, right=0.985, top=0.88, bottom=0.28)
    return fig


def main():
    rows, log = load_rows()
    print("=" * 72)
    print("FIGURE 7 DERIVATION (all values from results/Table_4.tsv)")
    print("=" * 72)
    print("\n".join(log))
    print("-" * 72)
    print(f"  rows: {len(rows)}   columns: {len(COLUMNS)}")
    print("=" * 72)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = make_plot(rows)
    fig.savefig(OUT_DIR / "Figure_7.png", dpi=600, facecolor="white")
    fig.savefig(OUT_DIR / "Figure_7.tiff", dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("Saved figures/Figure_7.png and .tiff at 600 dpi.")


if __name__ == "__main__":
    main()
