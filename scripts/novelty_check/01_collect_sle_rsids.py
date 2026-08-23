#!/usr/bin/env python3
"""
Step 1 of the novelty check.

Collect every rsID with a reported systemic lupus erythematosus association
from the GWAS Catalog REST API.

IMPORTANT ENDPOINT NOTES (both verified 2026-08-23):
  * efoTraits/MONDO_0007915/associations  -> HTTP 200  (use this)
  * efoTraits/EFO_0002690/associations    -> HTTP 404  (does NOT work)
Writes results/novelty_check/sle_associated_rsids_MONDO_0007915.json
"""
import urllib.request, json, re, time, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "our_outputs")
os.makedirs(OUTDIR, exist_ok=True)

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


def get(url, tries=3):
    for i in range(tries):
        try:
                return _fetch(url, 90)
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(3)

sle = set()
url = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0007915/associations?size=1000"
page = 0
while url and page < 12:
    d = get(url); page += 1
    for a in d.get("_embedded", {}).get("associations", []):
        for loc in a.get("loci", []):
            for ra in loc.get("strongestRiskAlleles", []):
                m = re.match(r"(rs\d+)", ra.get("riskAlleleName", ""))
                if m:
                    sle.add(m.group(1))
    nxt = d.get("_links", {}).get("next", {}).get("href")
    url = nxt if nxt and nxt != url else None

print(f"SLE-associated rsIDs collected: {len(sle)}")
# positive controls -- textbook SLE variants that MUST be present
for probe in ["rs2476601", "rs7574865", "rs2004640"]:
    print(f"  control {probe} present: {probe in sle}")
out = os.path.join(OUTDIR, "sle_associated_rsids_MONDO_0007915.json")
json.dump(sorted(sle), open(out, "w"))
print("wrote", out)
