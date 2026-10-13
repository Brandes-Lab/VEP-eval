#!/usr/bin/env python3
"""Download the benchmark data used by the figure notebooks.

The two input files are hosted on Hugging Face (https://huggingface.co/datasets/brandeslab/VEP-eval).
This script downloads them from a pinned dataset revision, checks their SHA-256 checksums, and saves
them under the paths the notebooks expect:

    data/processed/clinvar_benchmark_updated.csv   <- clinvar_benchmark.csv  (Figures 1-3, Supp. Figs 2-4, 6)
    cosmic/cosmic_all_scores.csv                   <- cosmic_benchmark.csv   (Supplementary Figure 5)

Run from the repository root:

    python scripts/download_data.py
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

REPO = "brandeslab/VEP-eval"
REVISION = "3bdb5477b65362293cf4adc41ca0ce7810b90f74"  # dataset revision used for the published figures

FILES = [
    ("clinvar_benchmark.csv", Path("data/processed/clinvar_benchmark_updated.csv"),
     "3ddd9e3f54522ff460a717ebff087119fca37eb2a76a843bcd4531f6f300bd36"),
    ("cosmic_benchmark.csv", Path("cosmic/cosmic_all_scores.csv"),
     "1c9dc1eb475a7ef2c0379bac7847fd263122964e4de1aceaa30993c3c1883e29"),
]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def download(name, target, expected):
    if target.exists() and sha256(target) == expected:
        print(f"{target}: already present, checksum OK")
        return
    url = f"https://huggingface.co/datasets/{REPO}/resolve/{REVISION}/{name}"
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + ".part")
    print(f"Downloading {url}")
    with urllib.request.urlopen(url) as response, open(partial, "wb") as handle:
        total = int(response.headers.get("Content-Length", 0))
        done = 0
        for block in iter(lambda: response.read(1 << 20), b""):
            handle.write(block)
            done += len(block)
            if total:
                print(f"\r  {done / 1e6:7.1f} / {total / 1e6:.1f} MB", end="", flush=True)
    print()
    actual = sha256(partial)
    if actual != expected:
        partial.unlink()
        sys.exit(f"Checksum mismatch for {name}: expected {expected}, got {actual}")
    partial.replace(target)
    print(f"{target}: saved, checksum OK")


if __name__ == "__main__":
    for name, target, expected in FILES:
        download(name, target, expected)
