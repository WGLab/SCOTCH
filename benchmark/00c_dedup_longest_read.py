#!/usr/bin/env python3
"""Deduplicate BAM by keeping only the longest read per CB+UMI combination.

SCOTCH performs this deduplication internally during its annotation step
(see annotation.py: bam_info_to_dict). For fair benchmarking, IsoQuant
inputs must be pre-filtered the same way, and the filtering time is counted
toward IsoQuant's total runtime.

For 10x-ONT, the default tags are:
  - Cell barcode: CB
  - UMI: UB

Logic:
  1. First pass: scan all reads, group by CB+UB, record the longest read name
  2. Second pass: write only the longest read per CB+UB group

Usage:
    python3 00c_dedup_longest_read.py <input.bam> <output.bam> [--cb CB --umi UB]
"""

import argparse
import pysam


def dedup_longest(input_bam, output_bam, cb_tag="CB", umi_tag="UB"):
    """Keep only the longest read per CB+UMI combination.

    Only primary alignments (not supplementary/secondary) are considered when
    choosing the longest read. In pass 2, all alignment records (including
    supplementary) for the winning read are kept, while all records for
    non-winning reads are discarded. This matches SCOTCH's behavior in
    annotation.py:bam_info_to_dict().
    """

    # --- Pass 1: find the longest PRIMARY read for each CB+UMI ---
    print(f"Pass 1: scanning {input_bam} for longest reads per {cb_tag}+{umi_tag}...")
    infile = pysam.AlignmentFile(input_bam, "rb")

    # cbumi -> (query_name, aligned_length)
    best = {}
    n_total = 0
    n_primary = 0
    n_with_tags = 0

    for read in infile:
        n_total += 1
        # Skip supplementary (0x800) and secondary (0x100) alignments for selection
        if read.is_supplementary or read.is_secondary:
            continue
        n_primary += 1
        if not read.has_tag(cb_tag) or not read.has_tag(umi_tag):
            continue
        n_with_tags += 1

        cb = read.get_tag(cb_tag)
        umi = read.get_tag(umi_tag)
        cbumi = f"{cb}_{umi}"
        length = read.query_alignment_length or 0

        if cbumi not in best or length > best[cbumi][1]:
            best[cbumi] = (read.query_name, length)

    infile.close()

    # Build set of winning query names
    keep_qnames = set(qname for qname, _ in best.values())

    print(f"  Total alignments: {n_total}")
    print(f"  Primary alignments: {n_primary}")
    print(f"  Primary with {cb_tag}+{umi_tag}: {n_with_tags}")
    print(f"  Unique {cb_tag}+{umi_tag} groups: {len(best)}")

    # --- Pass 2: write only records belonging to winning reads ---
    print(f"Pass 2: writing deduplicated reads to {output_bam}...")
    infile = pysam.AlignmentFile(input_bam, "rb")
    outfile = pysam.AlignmentFile(output_bam, "wb", header=infile.header)

    n_kept = 0
    n_skipped = 0

    for read in infile:
        if not read.has_tag(cb_tag) or not read.has_tag(umi_tag):
            # Keep reads without CB/UMI tags (e.g., unmapped) as-is
            outfile.write(read)
            n_kept += 1
            continue

        # Keep all alignment records (primary + supplementary) for the winning read
        if read.query_name in keep_qnames:
            outfile.write(read)
            n_kept += 1
        else:
            n_skipped += 1

    infile.close()
    outfile.close()

    # Index output
    pysam.index(output_bam)

    print(f"  Kept: {n_kept}, Skipped: {n_skipped}")
    print(f"  Dedup ratio: {n_skipped / max(n_with_tags, 1) * 100:.1f}% reads removed")
    return n_kept, n_skipped


def main():
    parser = argparse.ArgumentParser(
        description="Keep only the longest read per CB+UMI (same as SCOTCH annotation step)")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output deduplicated BAM file")
    parser.add_argument("--cb", default="CB", help="Cell barcode BAM tag (default: CB)")
    parser.add_argument("--umi", default="UB", help="UMI BAM tag (default: UB)")
    args = parser.parse_args()

    dedup_longest(args.input_bam, args.output_bam, args.cb, args.umi)
    print("Done.")


if __name__ == "__main__":
    main()
