#!/usr/bin/env python3
"""Prefix cell barcode (CB) tags with sample identifiers for multi-sample runs.

This prevents barcode collisions when IsoQuant (or SCOTCH) processes multiple
samples that may share identical cell barcodes.

Usage:
    python3 00b_prefix_barcodes.py <data_dir>

Reads:   multi_s1_5M.bam, multi_s2_5M.bam, multi_s3_5M.bam
Writes:  multi_s1_5M.prefixed.bam, multi_s2_5M.prefixed.bam, multi_s3_5M.prefixed.bam

The CB tag is modified: CB:Z:ACGTACGT-1 -> CB:Z:S1_ACGTACGT-1
"""

import os
import sys
import pysam


def prefix_barcodes(input_bam, output_bam, prefix):
    """Rewrite BAM with prefixed CB tags."""
    infile = pysam.AlignmentFile(input_bam, "rb")
    outfile = pysam.AlignmentFile(output_bam, "wb", header=infile.header)

    n_total = 0
    n_prefixed = 0

    for read in infile:
        n_total += 1
        if read.has_tag("CB"):
            cb = read.get_tag("CB")
            read.set_tag("CB", f"{prefix}_{cb}", value_type="Z")
            n_prefixed += 1
        outfile.write(read)

    infile.close()
    outfile.close()

    # Index output
    pysam.index(output_bam)

    print(f"  {input_bam} -> {output_bam}")
    print(f"  {n_prefixed}/{n_total} reads had CB tag prefixed with '{prefix}_'")


def main():
    data_dir = sys.argv[1]

    samples = [
        ("multi_s1_5M.bam", "S1"),
        ("multi_s2_5M.bam", "S2"),
        ("multi_s3_5M.bam", "S3"),
    ]

    for bam_name, prefix in samples:
        input_bam = os.path.join(data_dir, bam_name)
        output_bam = os.path.join(data_dir, bam_name.replace(".bam", ".prefixed.bam"))

        if not os.path.exists(input_bam):
            print(f"WARNING: {input_bam} not found, skipping")
            continue

        if os.path.exists(output_bam):
            print(f"Skipping {output_bam} (already exists)")
            continue

        print(f"Prefixing {bam_name} with '{prefix}'...")
        prefix_barcodes(input_bam, output_bam, prefix)

    print("Done.")


if __name__ == "__main__":
    main()
