#!/usr/bin/env python3
"""
Step 2 of the novelty check -- positional test.

For each locus labelled Putative Novel, fetch every GWAS Catalog SNP within
+/-500 kb and intersect with the SLE-associated rsID set from step 1.

Captures the FULL list of matching rsIDs together with their positions, so the
nearest-neighbour distance in step 3 cannot be affected by list truncation.

All paths are relative to this script's folder -- runnable from anywhere.
"""
import urllib.request, json, time, os, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "our_outputs"); os.makedirs(OUT, exist_ok=True)
MASTER = os.path.join(HERE, "master_results_table.tsv")
SLEJSON = os.path.join(OUT, "sle_associated_rsids_MONDO_0007915.json")
for f in (MASTER, SLEJSON):
    if not os.path.exists(f):
        sys.exit(f"missing input: {f}  (run 01_collect_sle_rsids.py first)")

WINDOW = 500_000
sle = set(json.load(open(SLEJSON)))
d = pd.read_csv(MASTER, sep="\t")
# Select every locus in the novelty-assessment set by LABEL, never by a hardcoded
# rsID list -- otherwise a locus reclassified in a previous round drops out of reruns.
LABELS = {
    "Putative Novel",
    "Putative Novel (near SLE signal; LD unresolved)",
    "Known (near reported SLE signal)",
}
p = d[d.Novelty.astype(str).str.strip().isin(LABELS)]
if len(p) != 21:
    print(f"WARNING: expected 21 loci in the novelty set, selected {len(p)}. "
          f"Labels present: {sorted(set(d.Novelty.astype(str)) & LABELS)}")
print(f"loci to check: {len(p)}  (SLE rsIDs in set: {len(sle)})")

# SSL note: some machines lack a usable CA bundle for urllib and fail with
# CERTIFICATE_VERIFY_FAILED. We prefer `requests` (which ships certifi) when it is
# installed and fall back to urllib otherwise.
#   fix if needed:  python3 -m pip install requests certifi
try:
    import requests

    def _fetch(url, timeout):
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        return r.json()
except ImportError:
    def _fetch(url, timeout):
        with urllib.request.urlopen(url, timeout=timeout) as f:
            return json.load(f)


def get(url, tries=4):
    for i in range(tries):
        try:
                return _fetch(url, 90)
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(3)

rows = []
for _, r in p.iterrows():
    rsid, chrom, bp = r["RSID"], int(r["CHR"]), int(r["BP"])
    hits, n, page = {}, 0, 0
    url = ("https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/"
           f"search/findByChromBpLocationRange?chrom={chrom}"
           f"&bpStart={max(1, bp - WINDOW)}&bpEnd={bp + WINDOW}&size=1000")
    while url and page < 10:
        j = get(url); page += 1
        for s in j.get("_embedded", {}).get("singleNucleotidePolymorphisms", []):
            n += 1
            name = s.get("rsId")
            if name in sle:
                for loc in s.get("locations", []):
                    if str(loc.get("chromosomeName")) == str(chrom):
                        hits[name] = int(loc["chromosomePosition"])
        nxt = j.get("_links", {}).get("next", {}).get("href")
        url = nxt if nxt and nxt != url else None
    if hits:
        near = min(hits.items(), key=lambda kv: abs(kv[1] - bp))
        dist = round(abs(near[1] - bp) / 1000, 1)
    else:
        near, dist = (None, None), None
    rows.append(dict(RSID=rsid, CHR=chrom, BP=bp,
                     Gene=(r["Gene"] if pd.notna(r["Gene"]) else ""),
                     catalog_snps_in_1Mb=n, sle_snps_in_1Mb=len(hits),
                     nearest_sle_snp=near[0], distance_kb=dist,
                     all_sle_snps_in_window=";".join(sorted(hits))))   # FULL, untruncated
    print(f"{rsid:12s} n_sle={len(hits):<3d} nearest={str(near[0]):14s} "
          f"{(str(dist)+' kb') if dist is not None else '-':>10s}", flush=True)

res = pd.DataFrame(rows)
out = os.path.join(OUT, "novelty_positional_FULL.tsv")
res.to_csv(out, sep="\t", index=False)
print("\nwrote", out)
ok = res[res.distance_kb.notna()]
print(f"loci with an SLE association within 500 kb: {len(ok)} of {len(res)}")
print(f"  <7 kb   : {(ok.distance_kb < 7).sum()}  -> {sorted(ok[ok.distance_kb < 7].RSID)}")
print(f"  7-100 kb: {((ok.distance_kb >= 7) & (ok.distance_kb < 100)).sum()}")
print(f"  >=100 kb: {(ok.distance_kb >= 100).sum()}")
