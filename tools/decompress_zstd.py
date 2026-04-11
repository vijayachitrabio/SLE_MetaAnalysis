import zstandard as zstd
import os

src = "data/ld_ref/1000G_EUR.tgz"
out = "data/ld_ref/1000G_EUR.tar"

print(f"Decompressing {src} to {out}...")
dctx = zstd.ZstdDecompressor()
with open(src, "rb") as ifh:
    with open(out, "wb") as ofh:
        dctx.copy_stream(ifh, ofh)
print("Done. Checking file type of the result...")
os.system(f"file {out}")
