#!/usr/bin/env python3
"""Decompress a zstd-compressed tar archive."""
import sys, tarfile
import zstandard as zstd

src = "data/ld_ref/1000G_EUR.tgz"
dst = "data/ld_ref/"

print(f"Decompressing {src} → {dst}")
with open(src, "rb") as f:
    dctx = zstd.ZstdDecompressor()
    with dctx.stream_reader(f) as sr:
        with tarfile.open(fileobj=sr, mode="r|*") as tar:
            print("Extracting files...")
            tar.extractall(dst)
print("Done.")
