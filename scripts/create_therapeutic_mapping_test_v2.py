from pathlib import Path
from textwrap import fill

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "therapeutic_mapping_summary.tsv"
OUT_PNG = ROOT / "figures" / "therapeutic_mapping_test_v2.png"
OUT_PDF = ROOT / "figures" / "therapeutic_mapping_test_v2.pdf"

STATUS_PRIORITY = {
    "FDA-approved SLE 2022": 1,
    "Approved": 2,
    "Approved / Phase III": 3,
    "Approved / SLE studies": 4,
    "FDA-approved (RA/SLE)": 5,
    "Phase II SLE": 6,
    "Investigational SLE": 7,
}

STATUS_COLORS = {
    "FDA-approved SLE 2022": "#215ca8",
    "Approved": "#2f8f5b",
    "Approved / Phase III": "#5f84c9",
    "Approved / SLE studies": "#8c6bb1",
    "FDA-approved (RA/SLE)": "#2e9d87",
    "Phase II SLE": "#d98a29",
    "Investigational SLE": "#8b8f97",
}


def clean_drugs(values):
    seen = set()
    items = []
    for value in values.dropna():
        for part in str(value).split(","):
            item = part.strip().replace("(Olumiant)", "").strip()
            if item and item not in seen:
                seen.add(item)
                items.append(item)
    return ", ".join(items)


def choose_row(group):
    ranked = group.copy()
    ranked["status_rank"] = ranked["Drug_Status"].map(STATUS_PRIORITY).fillna(99)
    ranked = ranked.sort_values(["status_rank", "Discovery_P"], ascending=[True, True])
    best = ranked.iloc[0].copy()
    best["Lead_Gene"] = group.name
    best["Specific_Drug"] = clean_drugs(group["Specific_Drug"])
    best["Novelty"] = "Novel" if (group["Novelty"] == "Novel").any() else group["Novelty"].iloc[0]
    return best


def load_data():
    df = pd.read_csv(INPUT, sep="\t")
    df["Discovery_P"] = pd.to_numeric(df["Discovery_P"], errors="coerce")
    out = df.groupby("Lead_Gene", group_keys=False).apply(choose_row).reset_index(drop=True)
    out["minus_log10_p"] = -np.log10(out["Discovery_P"])
    out = out.sort_values("minus_log10_p", ascending=False).reset_index(drop=True)
    return out


def make_plot(df):
    n = len(df)
    fig = plt.figure(figsize=(15.5, 6.8), dpi=300, facecolor="#fcfcfd")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.8, 2.1], wspace=0.05)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharey=ax1)

    yvals = np.arange(n)
    for i, y in enumerate(yvals):
        bg = "#f5f7fa" if i % 2 == 0 else "#ffffff"
        for ax in (ax1, ax2, ax3):
            ax.axhspan(y - 0.5, y + 0.5, color=bg, zorder=0)

    max_x = df["minus_log10_p"].max() * 1.08
    for y, row in zip(yvals, df.itertuples(index=False)):
        color = STATUS_COLORS.get(row.Drug_Status, "#6b7280")
        edge = "#d4a62a" if row.Novelty == "Novel" else "#243046"
        ax1.hlines(y, 0, row.minus_log10_p, color="#c5cfdb", linewidth=2.0, zorder=1)
        ax1.scatter(row.minus_log10_p, y, s=180, color=color, edgecolor=edge, linewidth=2.1, zorder=3)

        ax2.text(0.03, y, fill(str(row.Therapeutic_Pathway), width=26),
                 va="center", ha="left", fontsize=11.2, color="#334155", linespacing=1.35)

        ax3.text(0.02, y - 0.12, fill(str(row.Specific_Drug), width=31),
                 va="center", ha="left", fontsize=11.3, color="#111827", fontweight="bold", linespacing=1.25)
        ax3.text(0.02, y + 0.22, fill(str(row.Evidence_Level), width=34),
                 va="center", ha="left", fontsize=9.7, color="#64748b")
        ax3.text(0.98, y, row.Drug_Status,
                 va="center", ha="right", fontsize=10.1, color="white",
                 bbox=dict(boxstyle="round,pad=0.28,rounding_size=0.12", fc=color, ec="none"))

    ax1.set_xlim(0, max_x)
    ax1.set_ylim(-0.6, n - 0.4)
    ax1.invert_yaxis()
    ax1.set_yticks(yvals)
    ax1.set_yticklabels(df["Lead_Gene"], fontsize=12, fontweight="bold", color="#1f2937")
    ax1.set_xlabel(r"$-log_{10}(P_{discovery})$", fontsize=12, color="#334155")
    ax1.set_title("Statistical Support", fontsize=14, fontweight="bold", color="#334155", pad=10)
    ax1.grid(axis="x", color="#e2e8f0", linewidth=0.8)
    ax1.tick_params(axis="x", labelsize=10, colors="#475569")
    ax1.tick_params(axis="y", length=0)

    ax2.set_xlim(0, 1)
    ax2.set_xticks([])
    ax2.tick_params(axis="y", left=False, labelleft=False)
    ax2.set_title("Therapeutic Biology", fontsize=14, fontweight="bold", color="#334155", pad=10)

    ax3.set_xlim(0, 1)
    ax3.set_xticks([])
    ax3.tick_params(axis="y", left=False, labelleft=False)
    ax3.set_title("Drug Context and Status", fontsize=14, fontweight="bold", color="#334155", pad=10)

    for ax in (ax1, ax2, ax3):
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color("#cbd5e1")

    fig.text(0.05, 0.965, "Therapeutic Mapping of Replicated SLE Risk Loci",
             fontsize=20, fontweight="bold", color="#18263d", family="DejaVu Serif")
    fig.text(0.05, 0.925,
             "Test redesign: gene-centered summary linking discovery strength, targetable biology, and clinical positioning",
             fontsize=11.5, color="#64748b")

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#d1d5db",
               markeredgecolor="#d4a62a", markeredgewidth=2, markersize=9, label="Putatively novel locus"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#d1d5db",
               markeredgecolor="#243046", markeredgewidth=2, markersize=9, label="Previously reported locus"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.05, 0.02),
               frameon=False, fontsize=10, ncol=2)
    fig.subplots_adjust(top=0.83, bottom=0.12, left=0.12, right=0.985)
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
