"""
make_html_report.py — Microbesomics_BS HTML report

TODO: implement full report. Placeholder creates an empty HTML file so
      the Snakemake pipeline can complete end-to-end while this script
      is developed separately.

Arguments:
    --out   output HTML path
"""

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--out", required=True)
args = parser.parse_args()

os.makedirs(os.path.dirname(args.out), exist_ok=True)

with open(args.out, "w") as f:
    f.write("<html><body><h1>Microbesomics_BS — Report (placeholder)</h1></body></html>\n")

print(f"Written: {args.out}")
