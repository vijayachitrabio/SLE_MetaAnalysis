#!/usr/bin/env python3
"""
Step 3 -- report distances, ranked.

Reads the FULL output of step 2. Distances were computed there from the complete
set of in-window SLE rsIDs and their positions, so no truncation can affect them.
(An earlier version of this script parsed a 60-character preview column and was
therefore unsafe -- that column no longer exists.)
"""
import os, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "our_outputs", "novelty_positional_FULL.tsv")
if not os.path.exists(SRC):
    sys.exit(f"missing input: {SRC}  (run 02_positional_novelty.py first)")

df = pd.read_csv(SRC, sep="\t")
hits = df[df.distance_kb.notna()].sort_values("distance_kb")
cols = ["RSID", "CHR", "Gene", "nearest_sle_snp", "distance_kb", "sle_snps_in_1Mb"]
print("=== loci with a reported SLE association within 500 kb, ranked by proximity ===")
print(hits[cols].to_string(index=False))
print(f"\n  <7 kb    (treated as same signal) : {(hits.distance_kb < 7).sum()}")
print(f"  7-100 kb (LD unresolved)          : {((hits.distance_kb >= 7) & (hits.distance_kb < 100)).sum()}")
print(f"  >=100 kb                          : {(hits.distance_kb >= 100).sum()}")
print(f"  no SLE association within 500 kb  : {df.distance_kb.isna().sum()}")
print("\nNOTE: distance is a proxy. The rigorous test is LD (r2) between each lead")
print("and the nearest reported SLE variant in 1000G EUR. That was NOT run.")
