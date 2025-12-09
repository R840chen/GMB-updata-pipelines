#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bowtie2_pipeline.py (modified to use bam convert)

Features / changes:
 - SAM to BAM conversion now uses 'bam convert --in ... --out ...'
 - Other steps (sorting, indexing, extracting mapped reads) remain unchanged
 - Full English script without special characters
"""
import os
import sys
import glob
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

# ---------- utilities ----------
def safe_mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def find_paired_fastqs(fq_dir):
    files = sorted(glob.glob(os.path.join(fq_dir, "*")))
    pairs = {}
    for f in files:
        b = os.path.basename(f)
        if b.endswith("_1.fq") or b.endswith("_1.fq.gz") or b.endswith("_R1.fastq") or b.endswith("_R1.fastq.gz") or b.endswith("_1.fastq") or b.endswith("_1.fastq.gz"):
            prefix = b.rsplit("_1", 1)[0]
            pairs.setdefault(prefix, {})["r1"] = f
        elif b.endswith("_2.fq") or b.endswith("_2.fq.gz") or b.endswith("_R2.fastq") or b.endswith("_R2.fastq.gz") or b.endswith("_2.fastq") or b.endswith("_2.fastq.gz"):
            prefix = b.rsplit("_2", 1)[0]
            pairs.setdefault(prefix, {})["r2"] = f
    complete = {p: io for p, io in pairs.items() if "r1" in io and "r2" in io}
    return complete

def detect_index_base(index_dir, index_name):
    if os.path.exists(index_name):
        folder = index_name
    else:
        folder = os.path.join(index_dir, index_name)

    candidates = [
        os.path.join(folder, "bowtie2_" + os.path.basename(index_name)),
        os.path.join(folder, os.path.basename(index_name)),
        os.path.join(folder, os.path.basename(index_name) + ".bt2"),
        os.path.join(folder, "index"),
        folder
    ]

    for c in candidates:
        if os.path.exists(c + ".1.bt2") or os.path.exists(c + ".rev.1.bt2"):
            return c, False
        if os.path.exists(c + ".1.bt2l") or os.path.exists(c + ".rev.1.bt2l"):
            return c, True

    for c in candidates:
        parent = os.path.dirname(c) or "."
        prefix = os.path.basename(c)
        if os.path.isdir(parent):
            for entry in os.listdir(parent):
                if entry.startswith(prefix + ".") and (".bt2l" in entry or ".bt2" in entry):
                    is_large = ".bt2l" in entry
                    return os.path.join(parent, prefix), is_large

    raise FileNotFoundError(
        f"No Bowtie2 index pieces found for index '{index_name}' under '{index_dir}'. Checked candidates: {candidates}"
    )

def run_command(cmd_list):
    try:
        p = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, text=True)
        return p.returncode, p.stdout, p.stderr
    except Exception as e:
        return 1, "", str(e)

# ---------- core pipeline functions ----------
def generate_bowtie_commands(pairs, index_names, index_dir, txt_output_dir, sam_output_base, threads_bowtie, bowtie_extra_opts):
    planned = []
    safe_mkdir(txt_output_dir)
    safe_mkdir(sam_output_base)

    for prefix, io in sorted(pairs.items()):
        fq1 = io["r1"]
        fq2 = io["r2"]
        for idx in index_names:
            try:
                idx_base, is_large = detect_index_base(index_dir, idx)
            except FileNotFoundError as e:
                print(f"[ERROR] {e}", file=sys.stderr)
                sys.exit(2)

            sam_folder = os.path.join(sam_output_base, idx + "_bowtie_result")
            safe_mkdir(sam_folder)
            sam_file = prefix + "_" + idx + ".sam"
            sam_path = os.path.join(sam_folder, sam_file)

            cmd = ["bowtie2", "-p", str(threads_bowtie)]
            if is_large:
                cmd.append("--large-index")
            cmd += ["-x", idx_base, "-1", fq1, "-2", fq2, "-S", sam_path]
            if bowtie_extra_opts:
                cmd.extend(bowtie_extra_opts.split())

            txtfile = os.path.join(txt_output_dir, idx + "_bowtie2_command.txt")
            with open(txtfile, "a") as of:
                of.write(" ".join(cmd) + "\n")

            planned.append({
                "prefix": prefix,
                "index_name": idx,
                "index_base": idx_base,
                "is_large_index": is_large,
                "fq1": fq1,
                "fq2": fq2,
                "sam_path": sam_path,
                "cmd": cmd
            })
    return planned

def sam_to_mapped_txt(sam_path, threads_samtools):
    folder = os.path.dirname(sam_path)
    base = os.path.basename(sam_path).rsplit(".sam", 1)[0]
    bam_path = os.path.join(folder, base + ".bam")
    sorted_bam = os.path.join(folder, "sorted_" + base + ".bam")
    mapped_txt = os.path.join(folder, base + "_mapped_reads.txt")

    # SAM -> BAM using bam convert
    rc1, out1, err1 = run_command(["bam", "convert", "--in", sam_path, "--out", bam_path])
    if rc1 != 0:
        return mapped_txt, (rc1, err1)

    # sort with samtools
    rc2, out2, err2 = run_command(["samtools", "sort", bam_path, "-o", sorted_bam, "-@", str(threads_samtools)])
    if rc2 != 0:
        return mapped_txt, (rc2, err2)

    # index with samtools
    rc3, out3, err3 = run_command(["samtools", "index", sorted_bam, "-@", str(threads_samtools)])
    if rc3 != 0:
        return mapped_txt, (rc3, err3)

    # extract mapped reads
    proc = subprocess.Popen(["samtools", "view", "-F", "4", sorted_bam, "-@", str(threads_samtools)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    with open(mapped_txt, "w") as outf:
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5:
                qname = fields[0]
                rname = fields[2]
                pos = fields[3]
                mapq = fields[4]
                outf.write("\t".join([qname, rname, pos, mapq]) + "\n")
    proc.stdout.close()
    _, err4 = proc.communicate()
    rc4 = proc.returncode if proc.returncode is not None else 0
    if rc4 != 0:
        return mapped_txt, (rc4, err4)

    return mapped_txt, (0, "")

def process_sam_files_for_index(sam_folder, threads_samtools, max_workers):
    sam_files = sorted(glob.glob(os.path.join(sam_folder, "*.sam")))
    results = []
    if not sam_files:
        return results

    with ThreadPoolExecutor(max_workers=max_workers) as exe:
        future_map = {exe.submit(sam_to_mapped_txt, s, threads_samtools): s for s in sam_files}
        for fut in as_completed(future_map):
            sam = future_map[fut]
            mapped_txt, (rc, err) = fut.result()
            results.append((sam, mapped_txt, rc, err))
    return results

# ---------- main ----------
def main():
    parser = argparse.ArgumentParser(description="Bowtie2 + bam convert integrated pipeline")
    parser.add_argument("--fq_dir", required=True, help="Directory with paired fq files")
    parser.add_argument("--index_dir", required=True, help="Directory containing index subfolders")
    parser.add_argument("--index", action="append", help="Specific index name(s) to use (repeatable)")
    parser.add_argument("--txt_out_dir", required=True, help="Directory to write bowtie2 command txt files")
    parser.add_argument("--sam_out_dir", required=True, help="Base directory to write sam outputs (per index subfolder)")
    parser.add_argument("--threads_bowtie", type=int, default=10, help="Threads per bowtie2 process")
    parser.add_argument("--threads_samtools", type=int, default=1, help="Threads for samtools steps")
    parser.add_argument("--bowtie_parallel", type=int, default=8, help="Number of concurrent bowtie2 jobs")
    parser.add_argument("--samtool_parallel", type=int, default=4, help="Number of concurrent sam processing jobs")
    parser.add_argument("--mode", choices=["generate", "run", "all"], default="all",
                        help="generate only, run only, or run all")
    parser.add_argument("--clean_intermediate", action="store_true", help="Remove sam/bam/bai after processing")
    parser.add_argument("--bowtie_extra_opts", default="--sensitive", help="Extra options for bowtie2")
    args = parser.parse_args()

    pairs = find_paired_fastqs(args.fq_dir)
    if not pairs:
        print("No fq pairs found in", args.fq_dir, file=sys.stderr)
        sys.exit(1)
    print("Found", len(pairs), "paired fq files")

    if args.index:
        index_folders = args.index
        for idx in index_folders:
            candidate_folder = os.path.join(args.index_dir, idx)
            if not (os.path.isdir(candidate_folder) or os.path.exists(idx)):
                print(f"[ERROR] Specified index '{idx}' not found", file=sys.stderr)
                sys.exit(2)
    else:
        index_folders = sorted([d for d in os.listdir(args.index_dir) if os.path.isdir(os.path.join(args.index_dir, d))])
        if not index_folders:
            print("No index folders found in", args.index_dir, file=sys.stderr)
            sys.exit(1)

    print("Using index names:", index_folders)

    planned = generate_bowtie_commands(
        pairs, index_folders, args.index_dir,
        args.txt_out_dir, args.sam_out_dir,
        args.threads_bowtie, args.bowtie_extra_opts
    )
    print("Generated", len(planned), "bowtie2 command entries and wrote per-index txt files to", args.txt_out_dir)

    if args.mode == "generate":
        print("Mode generate: done")
        sys.exit(0)

    if args.mode in ("run", "all"):
        print("Running bowtie2 commands with parallel workers:", args.bowtie_parallel)
        results = []
        with ThreadPoolExecutor(max_workers=args.bowtie_parallel) as exe:
            future_map = {exe.submit(run_command, p["cmd"]): p for p in planned}
            for fut in as_completed(future_map):
                p = future_map[fut]
                rc, out, err = fut.result()
                results.append((p, rc, out, err))
                if rc != 0:
                    print("Bowtie2 failed for:", p["sam_path"], "rc:", rc, "err:", err, file=sys.stderr)

        print("Bowtie2 completed for all planned tasks")
        print("Processing sam files to extract mapped reads using bam convert + samtools")

        all_sam_folders = sorted([os.path.join(args.sam_out_dir, d) for d in os.listdir(args.sam_out_dir) if os.path.isdir(os.path.join(args.sam_out_dir, d))])
        for sf in all_sam_folders:
            print("Processing folder:", sf)
            res = process_sam_files_for_index(sf, args.threads_samtools, args.samtool_parallel)
            for sam, mapped_txt, rc, err in res:
                if rc != 0:
                    print("Error processing", sam, "rc:", rc, "err:", err, file=sys.stderr)
                else:
                    print("Mapped reads written to", mapped_txt)

        if args.clean_intermediate:
            print("Cleaning intermediate sam/bam/bai files")
            for sf in all_sam_folders:
                for ext in ("*.sam", "*.bam", "*.bai"):
                    for f in glob.glob(os.path.join(sf, ext)):
                        try:
                            os.remove(f)
                        except Exception:
                            pass

    print("Pipeline finished")

if __name__ == "__main__":
    main()
