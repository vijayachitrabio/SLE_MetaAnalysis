#!/usr/bin/env python3
import os
import pandas as pd
import pysam
import math
from pathlib import Path

# Paths
ROOT_DIR = Path("/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/GWAs_meta_analysis/AMH_MEnopause/SLE_MetaAnalysis")
LOCI_FILE = ROOT_DIR / "results_extracted/top_loci_summary_table.tsv"
GWAS_FILE = ROOT_DIR / "results/discovery_meta_results.tsv"
OUTPUT_DIR = ROOT_DIR / "results_extracted/coloc_inputs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# eQTL Catalogue GTEx v8 HTTPS Tabix URLs
URLS = {
    "Whole_Blood": "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000356/QTD000356.all.tsv.gz",
    "Spleen": "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000326/QTD000326.all.tsv.gz"
}

# eQTL columns
EQTL_COLS = [
    "molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", 
    "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2", 
    "molecular_trait_object_id", "gene_id", "median_tpm", "rsid"
]

print("Loading GWAS summary stats (this may take a moment)...")
gwas_df = pd.read_csv(GWAS_FILE, sep="\t", usecols=["CHR", "BP", "RSID", "OA", "EA", "BETA_meta", "SE_meta", "P_meta"])
gwas_df = gwas_df.dropna(subset=["BETA_meta", "SE_meta", "RSID"])
gwas_df["CHR"] = gwas_df["CHR"].astype(str)

print("Loading Top Loci...")
loci_df = pd.read_csv(LOCI_FILE, sep="\t")
print(f"Loaded {len(loci_df)} loci.")

for tissue, url in URLS.items():
    print(f"\nOpening Tabix for {tissue}...")
    try:
        tb = pysam.TabixFile(url)
    except Exception as e:
        print(f"Failed to open tabix for {tissue}: {e}")
        continue
    
    for idx, row in loci_df.iterrows():
        rsid = row["RSID"]
        chr_str = str(row["CHR"])
        pos = int(row["BP"])
        
        start = max(1, pos - 250000)
        end = pos + 250000
        
        print(f"  Fetching {tissue} for {rsid} ({chr_str}:{start}-{end})...")
        
        try:
            records = tb.fetch(chr_str, start, end)
        except ValueError:
            print(f"    Region {chr_str}:{start}-{end} not in tabix.")
            continue
            
        eqtl_data = []
        for rec in records:
            vals = rec.split("\t")
            if len(vals) < len(EQTL_COLS): continue
            eqtl_data.append(vals)
            
        if not eqtl_data:
            print(f"    No eQTLs found in region.")
            continue
            
        eqtl_df = pd.DataFrame(eqtl_data, columns=EQTL_COLS)
        eqtl_df["position"] = eqtl_df["position"].astype(int)
        
        gwas_region = gwas_df[(gwas_df["CHR"] == chr_str) & 
                              (gwas_df["BP"] >= start) & 
                              (gwas_df["BP"] <= end)].copy()
                              
        if len(gwas_region) == 0:
            print(f"    No GWAS variants in region.")
            continue
            
        merged = pd.merge(gwas_region, eqtl_df, left_on="BP", right_on="position", how="inner")
        
        agree_mask = (merged["EA"] == merged["alt"]) & (merged["OA"] == merged["ref"])
        flip_mask = (merged["EA"] == merged["ref"]) & (merged["OA"] == merged["alt"])
        
        retained = merged[agree_mask | flip_mask].copy()
        
        if len(retained) == 0:
            print(f"    No variants matched alleles.")
            continue
            
        retained["beta"] = pd.to_numeric(retained["beta"])
        retained["se"] = pd.to_numeric(retained["se"])
        retained["pvalue"] = pd.to_numeric(retained["pvalue"])
        retained["maf"] = pd.to_numeric(retained["maf"])
        
        is_flipped = (retained["EA"] == retained["ref"]) & (retained["OA"] == retained["alt"])
        retained.loc[is_flipped, "beta"] = -retained.loc[is_flipped, "beta"]
        
        final_df = retained[[
            "RSID", "BP", "molecular_trait_id", "gene_id", 
            "BETA_meta", "SE_meta", "P_meta",
            "beta", "se", "pvalue", "maf"
        ]].copy()
        
        final_df.columns = [
            "snp", "position", "gene", "ensembl_gene",
            "gwas_beta", "gwas_se", "gwas_pval",
            "eqtl_beta", "eqtl_se", "eqtl_pval", "eqtl_maf"
        ]
        
        final_df = final_df.dropna()
        final_df["gwas_N"] = 388655
        final_df["eqtl_N"] = 670 if tissue == "Whole_Blood" else 227
        
        out_file = OUTPUT_DIR / f"{rsid}_{tissue}.tsv"
        final_df.to_csv(out_file, sep="\t", index=False)
        
        genes_tested = final_df["gene"].nunique()
        print(f"    Saved {len(final_df)} variants across {genes_tested} genes for {rsid} in {tissue}.")

print("\nDone harmonization!")
