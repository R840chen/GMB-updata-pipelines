#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
select_and_simulate.py

1) Select core clusters from annotation_dir following coreness rules.
2) Write cat-seqs.fasta and core_clusters.fasta into outdir.
3) Run ART illumina to simulate reads from core_clusters.fasta and write outputs into outdir.

Usage example:
  python select_and_simulate.py \
    -i /path/to/annotation_dir \
    -o /path/to/outdir \
    --art /full/path/to/art_illumina \
    --threads 4

If ART is in PATH, you can omit --art and it will try "art_illumina".
"""
import argparse
import sys
import os
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
import subprocess

def parse_args():
    p = argparse.ArgumentParser(description="Select core clusters then simulate reads with ART.")
    p.add_argument("--annotation_dir", "-i", required=True,
                   help="Directory containing *-UniCluster-updated.ffn files (recursive).")
    p.add_argument("--outdir", "-o", required=True,
                   help="Output directory where cat-seqs.fasta, core_clusters.fasta and simulated reads will be written.")
    p.add_argument("--min_len", type=int, default=450,
                   help="Minimum length for selected core gene (default 450).")
    p.add_argument("--art", type=str, default="art_illumina",
                   help="Path to art_illumina executable (default: art_illumina in PATH).")
    p.add_argument("--art_profile", type=str, default="HS25",
                   help="ART profile (e.g. HS25).")
    p.add_argument("--read_len", type=int, default=150, help="Simulated read length (default 150).")
    p.add_argument("--fold_coverage", type=float, default=5.0, help="Fold coverage to simulate (default 5).")
    p.add_argument("--frag_mean", type=int, default=300, help="Fragment mean size for paired-end (default 300).")
    p.add_argument("--frag_std", type=int, default=10, help="Fragment std dev (default 10).")
    p.add_argument("--qual_lower", type=int, default=40, help="Base quality lower bound (default 40).")
    p.add_argument("--paired", action="store_true", help="Simulate paired-end reads (default: paired if flag set).")
    p.add_argument("--single", dest="paired", action="store_false", help="Simulate single-end (opposite).")
    p.add_argument("--coverage_mode", type=str, default="fold", choices=["fold"],
                   help="Coverage mode (currently only 'fold' supported meaning --f in ART).")
    return p.parse_args()

def decide_threshold(genomeN):
    if genomeN <= 3:
        return 2.0
    if genomeN == 4 or genomeN == 5:
        return 3.0
    if genomeN > 5 and genomeN <= 100:
        return float(genomeN) * 0.6
    return float(genomeN) * 0.5

def retrieve_corelist_and_cat(annotation_dir: Path, cat_path: Path):
    AllC = []
    GenomeN = 0
    annotation_dir = Path(annotation_dir)

    with open(cat_path, "w") as catfh:
        for f in sorted(annotation_dir.rglob("*-UniCluster-updated.ffn")):
            GenomeN += 1
            clusters = []
            try:
                for rec in SeqIO.parse(str(f), "fasta"):
                    parts = rec.description.split("|")
                    if len(parts) >= 2:
                        cl = parts[1]
                    else:
                        cl = "NA"
                    clusters.append(cl)
                    catfh.write(">" + str(rec.description).strip() + "\n" + str(rec.seq).strip() + "\n")
            except Exception as e:
                print(f"[WARN] Failed parsing {f}: {e}", file=sys.stderr)
                continue

            Uclusters = set(clusters)
            for item in Uclusters:
                if clusters.count(item) == 1:
                    AllC.append(item)

    return GenomeN, AllC

def build_cluster_and_seq_libs(cat_fasta_path: Path):
    SeqsLib = {}
    ClusterLib = defaultdict(list)
    for rec in SeqIO.parse(str(cat_fasta_path), "fasta"):
        fullname = str(rec.description).strip()
        parts = fullname.split("|")
        if len(parts) >= 2:
            cluster = parts[1]
        else:
            cluster = "NA"
        SeqsLib[fullname] = str(rec.seq).strip()
        ClusterLib[cluster].append(fullname)
    return SeqsLib, ClusterLib

def choose_longest_per_cluster(cluster_list, ClusterLib, SeqsLib, min_len):
    results = []
    for cl in cluster_list:
        members = ClusterLib.get(cl, [])
        if not members:
            continue
        best = max(members, key=lambda gid: len(SeqsLib.get(gid, "")))
        best_seq = SeqsLib.get(best, "")
        if best_seq and len(best_seq) >= min_len:
            results.append((best, best_seq))
    return results

def run_art_simulation(art_path, fasta_input: Path, outdir: Path, prefix_base: str,
                       profile="HS25", read_len=150, fold_coverage=5.0, paired=True,
                       frag_mean=300, frag_std=10, qual_lower=40):
    # ensure art exists
    art_exec = shutil_which(art_path)
    if art_exec is None:
        raise FileNotFoundError(f"ART executable not found: {art_path}")

    out_prefix = str(outdir / (prefix_base + "_"))
    # build command
    cmd = [
        art_exec,
        "-ss", profile,
        "-i", str(fasta_input),
        "-l", str(read_len),
        "-f", str(fold_coverage),
        "-o", out_prefix,
        "-qL", str(qual_lower)
    ]
    if paired:
        cmd.append("-p")
        cmd += ["-m", str(frag_mean), "-s", str(frag_std)]
    # run
    print("[INFO] Running ART:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("[INFO] ART finished; outputs with prefix:", out_prefix)

def shutil_which(name):
    # minimal replacement for shutil.which to avoid importing shutil if missing
    # but we can use shutil.which
    import shutil
    return shutil.which(name)

def main():
    args = parse_args()
    ann_dir = Path(args.annotation_dir)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    cat_file = out / "cat-seqs.fasta"
    core_fa = out / "core_clusters.fasta"

    print(f"[INFO] Scanning annotations in: {ann_dir}")
    GenomeN, AllC = retrieve_corelist_and_cat(ann_dir, cat_file)
    if GenomeN == 0:
        print("[ERROR] No *-UniCluster-updated.ffn files found in the annotation directory.", file=sys.stderr)
        sys.exit(2)
    print(f"[INFO] Found {GenomeN} genome files; collected {len(AllC)} per-genome single-occurrence cluster entries.")

    threshold = decide_threshold(GenomeN)
    print(f"[INFO] Using threshold = {threshold} (GenomeN = {GenomeN})")
    unique_clusters = set(AllC)
    passed = [cl for cl in unique_clusters if AllC.count(cl) >= threshold]
    print(f"[INFO] {len(passed)} clusters passed coreness threshold.")

    SeqsLib, ClusterLib = build_cluster_and_seq_libs(cat_file)
    chosen = choose_longest_per_cluster(passed, ClusterLib, SeqsLib, args.min_len)
    print(f"[INFO] {len(chosen)} sequences will be written to {core_fa}")

    with open(core_fa, "w") as ofh:
        for gid, seq in chosen:
            ofh.write(f">{gid}\n{seq}\n")

    # run ART on core_fa (if any sequences)
    if len(chosen) == 0:
        print("[WARN] No sequences passed selection; skipping ART simulation.")
        return

    # prefix base for output files: use outdir basename
    prefix_base = "ART-simulated-reads"
    try:
        run_art_simulation(
            args.art,
            core_fa,
            out,
            prefix_base,
            profile=args.art_profile,
            read_len=args.read_len,
            fold_coverage=args.fold_coverage,
            paired=args.paired,
            frag_mean=args.frag_mean,
            frag_std=args.frag_std,
            qual_lower=args.qual_lower
        )
    except FileNotFoundError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(3)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] ART failed: {e}", file=sys.stderr)
        sys.exit(4)

    print("[INFO] Done: core selection and read simulation complete.")

if __name__ == "__main__":
    main()
