#!/usr/bin/env python3
"""
Figure 3 -- Colocalization at the chr6p21.3/MHC locus (rs389884).

Replaces the previous CLIC1 colocalization figure, which reported values that
were not produced by the analysis pipeline (including a "Lung" tissue that was
never queried).

Every value is read from the verified eQTL Catalogue colocalization output.
Styled to match the original CLIC1 colocalization figure (lollipop / dot plot,
teal accent, dashed orange threshold line), extended to show the highest-
ranking genes tested at this locus rather than a single gene, so the
multi-gene nature of the MHC signal is explicit.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
COLOC = ROOT / "results" / "coloc_eqtl_catalogue_summary.tsv"
MAPPED = ROOT / "results" / "coloc_eqtl_catalogue_summary_mapped.tsv"
SYMBOLS = Path("/tmp/ens_symbols.json")

LEAD = "rs389884"
TOP_N = 12
THRESH = 0.80

# Original-style palette: light teal for below threshold, dark teal for above.
LIGHT = "#6fae9c"
DARK = "#1f6f56"
TRACK = "#dce8e4"
THRESHOLD_COLOR = "#c9701c"
TITLE_COLOR = "#0f2a3d"


def load():
    c = pd.read_csv(COLOC, sep="\t")
    c = c[c["Locus"] == LEAD].copy()

    sym = {}
    if SYMBOLS.exists():
        sym.update(json.load(open(SYMBOLS)))
    mp = pd.read_csv(MAPPED, sep="\t")
    for _, r in mp.iterrows():
        if isinstance(r.get("Symbol"), str):
            sym[r["Gene"]] = r["Symbol"]

    c["Symbol"] = c["Gene"].map(lambda g: sym.get(g, g))
    c = c[~c["Symbol"].str.startswith("ENSG")]
    n_total = len(pd.read_csv(COLOC, sep="\t").query("Locus == @LEAD"))
    c = c.sort_values("PP4", ascending=False).head(TOP_N)
    return c.sort_values("PP4", ascending=True), n_total


def main():
    df, n_total = load()
    labels = [f"{r.Symbol} ({r.Tissue.replace('_', ' ')})" for r in df.itertuples()]
    y = list(range(len(df)))

    fig, ax = plt.subplots(figsize=(7.6, 5.0), dpi=600, facecolor="white")

    for i, r in zip(y, df.itertuples()):
        strong = r.PP4 > THRESH
        color = DARK if strong else LIGHT
        ax.plot([0, r.PP4], [i, i], color=TRACK, lw=7, zorder=1,
                solid_capstyle="round")
        ax.scatter([r.PP4], [i], s=170, color=color, edgecolor="none", zorder=3)
        ax.text(r.PP4 + 0.03, i, f"{r.PP4:.2f}", va="center", ha="left",
                fontsize=8.6, color="#111111",
                fontweight="bold" if strong else "normal")

    ax.axvline(THRESH, color=THRESHOLD_COLOR, ls="--", lw=1.3, zorder=2)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9.2, color="#111111")
    for tick, r in zip(ax.get_yticklabels(), df.itertuples()):
        if r.Symbol == "CLIC1":
            tick.set_fontweight("bold")

    ax.set_xlim(0, 1.08)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Posterior probability of colocalization (PP4)",
                  fontsize=10.5, color="#111111", labelpad=10)
    ax.set_title("Colocalization at the chr6p21.3/MHC locus",
                 loc="left", fontsize=15.5, fontweight="bold",
                 color=TITLE_COLOR, pad=14)

    ax.grid(False)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="both", length=0, labelsize=9.6, colors="#334155")
    ax.set_axisbelow(True)

    fig.text(0.02, 0.015,
             f"225 gene-tissue pairs were tested at this locus; the {len(df)} "
             f"highest-ranking named genes are shown. Dark markers exceed the\n"
             f"PP4 > 0.80 threshold (dashed line).",
             fontsize=8.6, color="#4b5563")

    fig.subplots_adjust(left=0.28, right=0.95, top=0.90, bottom=0.19)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "Figure_3.png", dpi=600, facecolor="white")
    fig.savefig(OUT / "Figure_3.tiff", dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / "Figure_3.pdf", facecolor="white")
    plt.close(fig)

    print("FIGURE 3 -- values from", COLOC.name)
    print("-" * 62)
    for r in df.sort_values("PP4", ascending=False).itertuples():
        flag = "  *" if r.PP4 > THRESH else ""
        print(f"  {r.Symbol:16s} {r.Tissue:12s} PP4={r.PP4:.3f}  "
              f"variants={r.SNPs}{flag}")
    print("-" * 62)
    print(f"  total tests at locus: {n_total}")
    print("  saved figures/Figure_3.png, .tiff and .pdf")


if __name__ == "__main__":
    main()
