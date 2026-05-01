from pathlib import Path
from textwrap import fill

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "therapeutic_mapping_summary.tsv"
OUT_DIR = ROOT / "figures"

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

CLASS_SHORT = {
    "JAK Inhibitors": "JAK",
    "JAK1/2 Inhibitors": "JAK1/2",
    "Complement Inhibitors": "Complement",
    "Anti-IFN / JAK Pathway": "IFN/JAK",
    "Anti-OX40L": "OX40L",
    "IVIG / FcγR modulation": "Fc/IVIG",
    "Anti-IL-12/23": "IL-12/23",
    "mTOR Inhibitors": "mTOR",
}

CLASS_ORDER = ["JAK", "JAK1/2", "Complement", "IFN/JAK", "OX40L", "Fc/IVIG", "IL-12/23", "mTOR"]


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
    out["Class_Short"] = out["Drug_Class"].map(CLASS_SHORT).fillna(out["Drug_Class"])
    return out.sort_values("minus_log10_p", ascending=False).reset_index(drop=True)


def save(fig, stem):
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(OUT_DIR / f"{stem}.pdf", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def lane_plot(df):
    n = len(df)
    fig, ax = plt.subplots(figsize=(15, 6.6), dpi=300, facecolor="#fcfcfd")
    ax.set_facecolor("#fcfcfd")
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.7, n - 0.3)
    ax.invert_yaxis()
    ax.axis("off")

    for i, row in enumerate(df.itertuples(index=False)):
        bg = "#f5f7fa" if i % 2 == 0 else "#ffffff"
        ax.axhspan(i - 0.5, i + 0.5, color=bg, zorder=0)

        edge = "#d4a62a" if row.Novelty == "Novel" else "#243046"
        color = STATUS_COLORS.get(row.Drug_Status, "#6b7280")
        size = 70 + row.minus_log10_p * 2.0

        ax.text(2, i, row.Lead_Gene, va="center", ha="left", fontsize=13, fontweight="bold", color="#1f2937")
        ax.scatter(18, i, s=size, color=color, edgecolor=edge, linewidth=2, zorder=3)
        ax.text(22, i, f"-log10(P) {row.minus_log10_p:.1f}", va="center", ha="left", fontsize=10.5, color="#475569")

        ax.plot([38, 62], [i, i], color="#cbd5e1", linewidth=4, solid_capstyle="round")
        ax.text(50, i - 0.17, row.Class_Short, va="center", ha="center", fontsize=11.3, color="#334155", fontweight="bold")
        ax.text(50, i + 0.18, fill(str(row.Therapeutic_Pathway), width=26), va="center", ha="center", fontsize=9.5, color="#64748b")

        ax.text(66, i - 0.12, fill(str(row.Specific_Drug), width=34), va="center", ha="left", fontsize=11.2, color="#111827", fontweight="bold")
        ax.text(66, i + 0.23, fill(str(row.Evidence_Level), width=34), va="center", ha="left", fontsize=9.4, color="#64748b")
        ax.text(98.5, i, row.Drug_Status, va="center", ha="right", fontsize=10.2, color="white",
                bbox=dict(boxstyle="round,pad=0.28,rounding_size=0.12", fc=color, ec="none"))

    ax.text(2, -1.05, "Therapeutic Mapping of Replicated SLE Risk Loci", fontsize=21, fontweight="bold",
            color="#18263d", family="DejaVu Serif")
    ax.text(2, -0.72, "Test v3A: lane-style summary balancing figure aesthetics with manuscript readability",
            fontsize=11.4, color="#64748b")

    ax.text(2, -0.15, "Gene", fontsize=12, fontweight="bold", color="#334155")
    ax.text(18, -0.15, "Strength", fontsize=12, fontweight="bold", color="#334155", ha="center")
    ax.text(50, -0.15, "Target Class and Biology", fontsize=12, fontweight="bold", color="#334155", ha="center")
    ax.text(66, -0.15, "Drug Context", fontsize=12, fontweight="bold", color="#334155")

    return fig


def matrix_plot(df):
    n = len(df)
    xlabels = CLASS_ORDER
    xmap = {label: i for i, label in enumerate(xlabels)}
    fig, ax = plt.subplots(figsize=(13.5, 6.8), dpi=300, facecolor="#fcfcfd")
    ax.set_facecolor("#fcfcfd")

    for i in range(n):
        bg = "#f5f7fa" if i % 2 == 0 else "#ffffff"
        ax.axhspan(i - 0.5, i + 0.5, color=bg, zorder=0)

    for row_i, row in enumerate(df.itertuples(index=False)):
        color = STATUS_COLORS.get(row.Drug_Status, "#6b7280")
        edge = "#d4a62a" if row.Novelty == "Novel" else "#243046"
        x = xmap.get(row.Class_Short, 0)
        ax.scatter(x, row_i, s=90 + row.minus_log10_p * 3, color=color, edgecolor=edge, linewidth=2.2, zorder=3)
        ax.text(len(xlabels) - 0.05, row_i - 0.14, fill(str(row.Specific_Drug), width=28),
                ha="left", va="center", fontsize=10.6, color="#111827", fontweight="bold")
        ax.text(len(xlabels) - 0.05, row_i + 0.19, row.Drug_Status,
                ha="left", va="center", fontsize=9.6, color=color)

    ax.set_xlim(-0.6, len(xlabels) + 2.7)
    ax.set_ylim(-0.7, n - 0.3)
    ax.invert_yaxis()
    ax.set_yticks(range(n))
    ax.set_yticklabels(df["Lead_Gene"], fontsize=12, fontweight="bold", color="#1f2937")
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels, fontsize=10.5, color="#334155")
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#cbd5e1")

    ax.set_title("Test v3B: simplified therapy-class matrix", loc="left",
                 fontsize=11.4, color="#64748b", pad=36)
    fig.text(0.06, 0.955, "Therapeutic Mapping of Replicated SLE Risk Loci",
             fontsize=21, fontweight="bold", color="#18263d", family="DejaVu Serif")
    fig.text(0.06, 0.915, "Genes mapped to therapy classes with cleaner status encoding and fewer competing elements",
             fontsize=11.4, color="#64748b")

    legend_y = -0.95
    ax.scatter(-0.4, legend_y, s=140, color="#d1d5db", edgecolor="#d4a62a", linewidth=2)
    ax.text(-0.15, legend_y, "Putatively novel locus", va="center", fontsize=10.2, color="#334155")
    ax.scatter(1.9, legend_y, s=140, color="#d1d5db", edgecolor="#243046", linewidth=2)
    ax.text(2.15, legend_y, "Previously reported locus", va="center", fontsize=10.2, color="#334155")

    return fig


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_data()
    save(lane_plot(df), "therapeutic_mapping_test_v3_lane")
    save(matrix_plot(df), "therapeutic_mapping_test_v3_matrix")


if __name__ == "__main__":
    main()
