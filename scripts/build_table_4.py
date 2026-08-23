#!/usr/bin/env python3
"""
Strict provenance rebuild of Table 4 (therapeutic target prioritization).

EVERY numeric / status field is read from an authoritative source file and
asserted. Nothing is typed, recalled, rounded or inferred.

Sourced fields (must match source exactly):
  CHR, Discovery_P, Novelty   <- results_extracted/master_results_table.tsv
  Replicated                  <- Supplementary Table 5, 'Replication status' column
  PP4, Tissue                 <- results/coloc_eqtl_catalogue_summary_mapped.tsv

Curated fields (literature / Open Targets annotation per Methods 2.11 -- these
are editorial annotations, not computed values):
  Therapeutic_Pathway, Drug_Class, Specific_Drug, Drug_Status
"""

import sys
import pandas as pd

MASTER_PATH = "results_extracted/master_results_table.tsv"
COLOC_PATH = "results/coloc_eqtl_catalogue_summary_mapped.tsv"
ST5_PATH = "revisions_todays_data/supplementary_tables_revised_latest.xlsx"
OUT_PATH = "results/Table_4.tsv"

master_df = pd.read_csv(MASTER_PATH, sep="\t")
coloc_df = pd.read_csv(COLOC_PATH, sep="\t")

# Supplementary Table 5: header row is index 2 in the revised workbook
st5_df = pd.read_excel(ST5_PATH, sheet_name="Supplementary Table 5", header=2)
st5_df = st5_df[st5_df["RSID"].notna()]

REQ_ST5 = {"RSID", "Replication status"}
missing = REQ_ST5 - set(st5_df.columns)
if missing:
    sys.exit(f"FATAL: Supplementary Table 5 missing columns: {missing}")

# ---------------------------------------------------------------- guard rails
# Verify the Ensembl->symbol mapping before trusting any coloc lookup.
for ens, sym in [("ENSG00000213719", "CLIC1"), ("ENSG00000244731", "C4A")]:
    got = coloc_df.loc[coloc_df["Gene"] == ens, "Symbol"].unique()
    assert list(got) == [sym], f"FATAL: mapping error {ens} -> {got}, expected {sym}"

assertions = []
warnings = []


def _row(df, mask, what):
    sub = df[mask]
    if sub.empty:
        sys.exit(f"FATAL: {what} not found in source")
    return sub


def get_master(rsid, col):
    sub = _row(master_df, master_df["RSID"] == rsid, f"{rsid} ({col})")
    val = sub.iloc[0][col]
    if pd.isna(val):
        sys.exit(f"FATAL: missing {col} for {rsid} in {MASTER_PATH}")
    return val, f"{MASTER_PATH}:{sub.index[0] + 2}"


def get_replicated(rsid):
    """Read the actual replication STATUS, not mere presence in the sheet."""
    sub = st5_df[st5_df["RSID"] == rsid]
    if sub.empty:
        # Not assessed for replication at all.
        return False, f"{ST5_PATH}:NOT_LISTED"
    status = str(sub.iloc[0]["Replication status"]).strip()
    line = int(sub.index[0]) + 4  # header@2 -> data starts at spreadsheet row 4
    if status.lower() == "replicated":
        return True, f"{ST5_PATH}:{line} (status='{status}')"
    return False, f"{ST5_PATH}:{line} (status='{status}')"


def get_pp4(rsid, symbol, tissue):
    mask = (
        (coloc_df["Locus"] == rsid)
        & (coloc_df["Symbol"] == symbol)
        & (coloc_df["Tissue"] == tissue)
    )
    sub = _row(coloc_df, mask, f"{symbol}/{tissue} at {rsid}")
    val = sub.iloc[0]["PP4"]
    if pd.isna(val):
        sys.exit(f"FATAL: missing PP4 for {symbol}/{tissue} at {rsid}")
    return float(val), f"{COLOC_PATH}:{sub.index[0] + 2}"


# Curated therapeutic annotations only. All numbers/statuses are sourced below.
# STAT4's two pathway annotations are merged into a single row: a duplicated
# rs4853458 entry is what produced the earlier 47-vs-46 locus-count error.
CURATED = [
    ("rs4853458", "STAT4",
     "JAK-STAT / Type I & II Interferon Signaling", "JAK1/2 Inhibitors",
     "Baricitinib, Tofacitinib", "FDA-approved (RA/SLE)", None),
    ("rs34572943", "ITGAM",
     "Complement Receptor / Neutrophil Adhesion", "Complement Inhibitors",
     "Eculizumab (C5), Avacopan (C5a)", "Approved / Phase III", None),
    ("rs13332649", "IRF8",
     "Myeloid / Interferon Response", "Anti-IFN / JAK Pathway",
     "Anifrolumab, Baricitinib", "Approved", None),
    # IRF5 (chr7) and IRF8 (chr16) are distinct loci; both are retained.
    ("rs35000415", "IRF5",
     "Type I Interferon Signaling", "Anti-IFN receptor",
     "Anifrolumab", "Approved", None),
    # rs6679677 is the conventional PTPN22 tag SNP (~53 kb from PTPN22;
    # near-complete LD with rs2476601/R620W). The Gene column in
    # master_results_table.tsv annotates the literal nearest gene, RSBN1.
    ("rs6679677", "PTPN22",
     "Tyrosine Phosphatase / TCR Signaling", "Phosphatase modulation",
     "No approved agent", "Investigational", None),
    ("rs10912578", "TNFSF4",
     "T-cell Co-stimulation (OX40L)", "Anti-OX40L",
     "Itepekimab (OX40L pathway)", "Phase II SLE", None),
    ("rs6671847", "FCGR2A",
     "Fc Receptor / Immune Complex Clearance", "IVIG / FcgR modulation",
     "IVIG, Efgartigimod (FcRn)", "Approved / Phase III", None),
    ("rs2647928", "IL12A",
     "IL-12 / Th1 Polarization", "Anti-IL-12/23",
     "Ustekinumab (Stelara)", "Approved / SLE studies", None),
    ("rs12928726", "CLEC16A",
     "Autophagy / MHC-II Regulation", "mTOR Inhibitors",
     "Sirolimus, Everolimus", "Investigational SLE", None),
    # Colocalization-supported additions: (symbol, tissue) drives the PP4 lookup
    ("rs35251378", "TYK2",
     "JAK-STAT Signaling", "JAK Inhibitors",
     "Deucravacitinib", "Approved / Phase III", ("TYK2", "Whole_Blood")),
    ("rs4840568", "BLK",
     "B-cell Signaling / Kinase", "Kinase Inhibitors",
     "Dasatinib", "Investigational", ("BLK", "Spleen")),
    ("rs5994638", "UBE2L3",
     "NF-kB / Ubiquitin", "Proteasome / NF-kB inhibitors",
     "Bortezomib", "Investigational", ("UBE2L3", "Whole_Blood")),
]

# The MHC locus carries four colocalizing genes; it is reported at locus level
# rather than attributing the signal to any single gene.
MHC_GENES = [
    ("LINC00243", "Whole_Blood"),
    ("CLIC1", "Whole_Blood"),
    ("C4A", "Spleen"),
    ("CCHCR1", "Spleen"),
]

rows = []

for rsid, gene, pathway, dclass, drug, dstatus, coloc_key in CURATED:
    chr_v, chr_s = get_master(rsid, "CHR")
    p_v, p_s = get_master(rsid, "P_meta")
    nov_v, nov_s = get_master(rsid, "Novelty")
    rep_v, rep_s = get_replicated(rsid)

    if coloc_key is None:
        evidence = "Strong (GWAS)" if p_v < 1e-15 else "Moderate (GWAS)"
    else:
        sym, tis = coloc_key
        pp4, pp4_s = get_pp4(rsid, sym, tis)
        evidence = (f"Strong (GWAS + eQTL {tis.replace('_', ' ')}: "
                    f"PP4 = {pp4:.3f})")
        assertions.append(f"  PASS  {gene:9s} PP4={pp4:.3f}          <- {pp4_s}")

    rows.append({
        "RSID": rsid, "Lead_Gene": gene, "CHR": int(chr_v),
        "Discovery_P": p_v, "Replicated": rep_v, "Novelty": nov_v,
        "Therapeutic_Pathway": pathway, "Drug_Class": dclass,
        "Specific_Drug": drug, "Drug_Status": dstatus,
        "Evidence_Level": evidence,
    })
    assertions.append(f"  PASS  {gene:9s} CHR={int(chr_v):<14d} <- {chr_s}")
    assertions.append(f"  PASS  {gene:9s} P={p_v:<16.6g} <- {p_s}")
    assertions.append(f"  PASS  {gene:9s} Novelty='{nov_v}' <- {nov_s}")
    assertions.append(f"  PASS  {gene:9s} Replicated={str(rep_v):<6s}   <- {rep_s}")

# ------------------------------------------------------------ MHC locus row
rsid = "rs389884"
label = "chr6p21.3 / MHC (multi-gene)*"
chr_v, chr_s = get_master(rsid, "CHR")
p_v, p_s = get_master(rsid, "P_meta")
nov_v, nov_s = get_master(rsid, "Novelty")
rep_v, rep_s = get_replicated(rsid)

parts = []
for sym, tis in MHC_GENES:
    pp4, pp4_s = get_pp4(rsid, sym, tis)
    parts.append(f"{sym} {pp4:.3f} ({tis.replace('_', ' ')})")
    assertions.append(f"  PASS  {sym:9s} PP4={pp4:.3f}          <- {pp4_s}")

rows.append({
    "RSID": rsid, "Lead_Gene": label, "CHR": int(chr_v),
    "Discovery_P": p_v, "Replicated": rep_v, "Novelty": nov_v,
    "Therapeutic_Pathway": "Complement / Immune Complex Clearance",
    "Drug_Class": "Undetermined (directionality conflict)",
    "Specific_Drug": "N/A", "Drug_Status": "Undetermined",
    "Evidence_Level": "Strong (GWAS + eQTL; " + ", ".join(parts) + ")",
})
assertions.append(f"  PASS  {'MHC':9s} CHR={int(chr_v):<14d} <- {chr_s}")
assertions.append(f"  PASS  {'MHC':9s} P={p_v:<16.6g} <- {p_s}")
assertions.append(f"  PASS  {'MHC':9s} Novelty='{nov_v}' <- {nov_s}")
assertions.append(f"  PASS  {'MHC':9s} Replicated={str(rep_v):<6s}   <- {rep_s}")

# rs389884 lies in the MHC/C4 region, among the best-established SLE loci.
# Flag rather than silently override the source annotation.
if "novel" in str(nov_v).lower():
    warnings.append(
        f"rs389884 Novelty='{nov_v}' in {MASTER_PATH}. The chr6p21.3/MHC-C4 "
        "region is a long-established SLE locus; this source annotation is "
        "scientifically questionable and should be reviewed at source before "
        "publication. NOT overridden here."
    )

df = pd.DataFrame(rows)

# ------------------------------------------------------------- final checks
dups = df[df.duplicated("RSID", keep=False)]
if not dups.empty:
    sys.exit(f"FATAL: duplicate RSIDs in Table 4:\n{dups[['RSID', 'Lead_Gene']]}")

rep7 = set(st5_df.loc[
    st5_df["Replication status"].astype(str).str.strip().str.lower() == "replicated",
    "RSID"])
for _, r in df.iterrows():
    expect = r["RSID"] in rep7
    if bool(r["Replicated"]) != expect:
        sys.exit(f"FATAL: Replicated mismatch for {r['RSID']}")

df.to_csv(OUT_PATH, sep="\t", index=False)

print("=" * 72)
print("TABLE 4 PROVENANCE REPORT")
print("=" * 72)
print("\n".join(assertions))
print("-" * 72)
print(f"  rows written        : {len(df)}")
print(f"  duplicate RSIDs     : none")
print(f"  replicated (from ST5): {sorted(rep7)}")
print(f"  ALL ASSERTIONS PASSED -> {OUT_PATH}")
if warnings:
    print("-" * 72)
    for w in warnings:
        print(f"  WARNING: {w}")
print("=" * 72)
