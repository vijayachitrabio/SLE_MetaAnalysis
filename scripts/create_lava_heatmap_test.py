from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "final_lava_consolidated_loci.tsv"
OUT_PNG = ROOT / "figures" / "lava_heatmap_test.png"
OUT_PDF = ROOT / "figures" / "lava_heatmap_test.pdf"


def load_data():
    df = pd.read_csv(INPUT, sep="\t")
    df = df[df["Final_Assessment"] == "HIGH CONFIDENCE"].copy()
    df["display"] = df["Gene"].fillna("").replace("", np.nan)
    df["display"] = df["display"].fillna(df["RSID"])
    df["minus_log10_p"] = -np.log10(pd.to_numeric(df["P_meta"], errors="coerce"))
    df["LAVA_rg"] = pd.to_numeric(df["LAVA_rg"], errors="coerce")
    df["LAVA_P"] = pd.to_numeric(df["LAVA_P"], errors="coerce")
    df = df.sort_values(["minus_log10_p", "CHR", "BP"], ascending=[False, True, True]).reset_index(drop=True)
    return df


def draw_cell(ax, x, y, color, text=None, text_color="#111827", lw=0.8):
    rect = patches.Rectangle((x, y), 1, 1, facecolor=color, edgecolor="white", linewidth=lw)
    ax.add_patch(rect)
    if text is not None:
        ax.text(x + 0.5, y + 0.5, text, ha="center", va="center", fontsize=10, color=text_color, weight="bold")


def make_plot(df):
    n = len(df)
    fig = plt.figure(figsize=(max(12, n * 0.62), 4.8), dpi=300, facecolor="#fcfcfd")
    ax = fig.add_axes([0.08, 0.23, 0.84, 0.50])
    ax.set_xlim(0, n)
    ax.set_ylim(0, 2)
    ax.invert_yaxis()
    ax.set_facecolor("#fcfcfd")

    sig_norm = Normalize(vmin=df["minus_log10_p"].min(), vmax=df["minus_log10_p"].max())
    rg_norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    sig_cmap = plt.cm.YlOrBr
    rg_cmap = plt.cm.PuOr

    for i, row in enumerate(df.itertuples(index=False)):
        draw_cell(ax, i, 0, sig_cmap(sig_norm(row.minus_log10_p)))

        if np.isnan(row.LAVA_rg):
            draw_cell(ax, i, 1, "#eef2f7", text="NA", text_color="#64748b")
        else:
            star = ""
            if not np.isnan(row.LAVA_P):
                if row.LAVA_P < 0.001:
                    star = "***"
                elif row.LAVA_P < 0.05:
                    star = "*"
            draw_cell(ax, i, 1, rg_cmap(rg_norm(row.LAVA_rg)), text=star if star else None)

    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_xticklabels(df["display"], rotation=45, ha="right", fontsize=10, color="#334155")
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(["Discovery\n-log10(P)", "LAVA local\nrg"], fontsize=11, color="#334155", weight="bold")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.text(
        0.08,
        0.92,
        "Prototype Heatmap of High-Confidence SLE Loci",
        fontsize=18,
        fontweight="bold",
        color="#18263d",
        family="DejaVu Serif",
    )
    fig.text(
        0.08,
        0.875,
        "Adapted from the LDSC/LAVA-style idea using available data in this project: discovery significance on top and local discovery-replication LAVA correlation on bottom",
        fontsize=10.8,
        color="#64748b",
    )

    # Colorbars
    cax1 = fig.add_axes([0.08, 0.12, 0.24, 0.035])
    cb1 = plt.colorbar(ScalarMappable(norm=sig_norm, cmap=sig_cmap), cax=cax1, orientation="horizontal")
    cb1.outline.set_visible(False)
    cb1.ax.tick_params(labelsize=9, colors="#475569")
    cb1.set_label("Discovery significance", fontsize=10, color="#334155")

    cax2 = fig.add_axes([0.40, 0.12, 0.24, 0.035])
    cb2 = plt.colorbar(ScalarMappable(norm=rg_norm, cmap=rg_cmap), cax=cax2, orientation="horizontal")
    cb2.outline.set_visible(False)
    cb2.ax.tick_params(labelsize=9, colors="#475569")
    cb2.set_label("LAVA local correlation (rg)", fontsize=10, color="#334155")

    fig.text(0.72, 0.12, "* nominal LAVA P < 0.05   *** Bonferroni-like strong support", fontsize=9.7, color="#334155")

    return fig


def main():
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    df = load_data()
    fig = make_plot(df)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    main()
