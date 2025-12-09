#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
simulate_reads_with_art.py

Find core_clusters.fasta in the provided input directory, keep sequences longer
than min length (default 450bp), simulate paired-end reads with art_illumina,
and write outputs to the input directory with prefix "ART-simulated-reads_".

Example:
  python simulate_reads_with_art.py --input-dir /path/to/sample
"""

import argparse
import os
import subprocess
import sys
import tempfile
import shutil

def parse_args():
    parser = argparse.ArgumentParser(description="Filter core_clusters.fasta by length and simulate reads with art_illumina.")
    parser.add_argument("--input-dir", "-i", required=True, help="Directory that contains core_clusters.fasta and where outputs will be written.")
    parser.add_argument(
        "--art-path", 
        "-a", 
        default="art_illumina",
        help="Path to art_illumina executable (default: use 'art_illumina' from system PATH)."
    )
    parser.add_argument("--min-len", type=int, default=450, help="Minimum sequence length to keep (default: 450).")
    parser.add_argument("--coverage", "-f", default="5", help="Fold coverage for ART -f (default 5).")
    parser.add_argument("--art-params", default="", help="Optional extra parameters to pass to art_illumina.")
    return parser.parse_args()

def filter_fasta_by_length(input_fasta, output_fasta, min_len):
    written = 0
    with open(input_fasta, "r") as inf, open(output_fasta, "w") as outf:
        header = None
        seq_lines = []
        for line in inf:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_lines)
                    if len(seq) > min_len:
                        outf.write(header + "\n")
                        for i in range(0, len(seq), 80):
                            outf.write(seq[i:i+80] + "\n")
                        written += 1
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # final record
        if header is not None:
            seq = "".join(seq_lines)
            if len(seq) > min_len:
                outf.write(header + "\n")
                for i in range(0, len(seq), 80):
                    outf.write(seq[i:i+80] + "\n")
                written += 1
    return written

def run_art(art_path, input_fasta, output_prefix, coverage="5", extra_params=""):
    cmd = [
        art_path,
        "-ss", "HS25",
        "-i", input_fasta,
        "-l", "150",
        "-f", str(coverage),
        "-o", output_prefix,
        "-qL", "40",
        "-p",
        "-m", "300",
        "-s", "10",
        "-k", "0"
    ]

    if extra_params:
        cmd.extend(extra_params.split())

    print("Running art_illumina:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)

def main():
    args = parse_args()

    input_dir = os.path.abspath(args.input_dir)
    art_path = args.art_path
    min_len = args.min_len
    coverage = args.coverage
    extra = args.art_params.strip()

    if not os.path.isdir(input_dir):
        print(f"ERROR: Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(2)

    # check art_illumina availability
    if shutil.which(art_path) is None:
        print(f"ERROR: art_illumina not found: {art_path}", file=sys.stderr)
        print("Make sure art_illumina is installed or provide --art-path /full/path/to/art_illumina")
        sys.exit(2)

    core_fasta = os.path.join(input_dir, "core_clusters.fasta")
    if not os.path.exists(core_fasta):
        print(f"ERROR: core_clusters.fasta not found in {input_dir}", file=sys.stderr)
        sys.exit(2)

    # temporary filtered fasta
    fd, filtered_fasta = tempfile.mkstemp(prefix="filtered_core_", suffix=".fasta", dir=input_dir)
    os.close(fd)

    try:
        print(f"Filtering sequences longer than {min_len} bp ...")
        kept = filter_fasta_by_length(core_fasta, filtered_fasta, min_len)
        print(f"Kept {kept} sequences.")

        if kept == 0:
            print("No sequences longer than threshold. Exiting.")
            os.remove(filtered_fasta)
            sys.exit(0)

        output_prefix = os.path.join(input_dir, "ART-simulated-reads_")
        run_art(art_path, filtered_fasta, output_prefix, coverage=coverage, extra_params=extra)

        print("Simulation completed.")
        print(f"Output files written to directory: {input_dir}")

    finally:
        if os.path.exists(filtered_fasta):
            os.remove(filtered_fasta)

if __name__ == "__main__":
    main()
