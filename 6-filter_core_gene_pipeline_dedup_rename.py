#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
6-filter_core_gene_pipeline_dedup_rename.py

Combined pipeline:
 1) Dedup merged mapped summaries (keep highest MAPQ per gene id), collect mapped gene ids,
    and filter original core-gene fasta files to remove those mapped gene ids.
 2) On the filtered core-gene fasta files, run length filtering (<= max_len) and
    prioritize/dedupe (bac120), keep between min_keep and max_keep genes per file,
    then write outputs named by species_name (normalized) with .ffn extension into outdir.

Simplified CLI (single output path):
  --merged_folder   : folder with merged mapped summary txt files (input)
  --core_fasta_dir  : folder containing original core-gene fasta/ffn files (input)
  --tax_infofile    : tax TSV used to map filenames -> species_name
  --outdir          : single output folder (used for dedup sorted files, filtered fastas and final species .ffn)
  --bac120_tsv      : optional bac120 TSV (default 'bac120-annotated-single-copy-gene.tsv' in cwd)
  --threads         : parallelism for dedup/filter step (default 8)
  --max_len, --max_keep, --min_keep, --dryrun as before.

Notes:
 - This script always executes both steps (no flags). If you want a dry run, use --dryrun.
 - Requires Biopython for robust fasta parsing. If missing the script still tries a fallback parser.
"""

from __future__ import annotations
import argparse
import os
import sys
import re
import csv
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import List, Tuple, Dict

# Try import SeqIO (preferred); provide fallback if absent
try:
    from Bio import SeqIO
except Exception:
    SeqIO = None

# ---------------- utilities ----------------
def safe_mkdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def normalize(s: str) -> str:
    if s is None:
        return ""
    return re.sub(r'\s+', '_', s.strip().lower())

# ---------------- Step1: dedup / collect mapped gene ids / filter fasta ----------------
def parse_line_for_read_and_mapq(line: str):
    s = line.rstrip("\n")
    if not s:
        return None, 0, s
    parts = s.split("\t")
    read_name = None
    mapq = 0
    try:
        if len(parts) >= 5:
            read_name = parts[0]
            mapq = int(parts[4]) if parts[4].isdigit() else (int(parts[-1]) if parts[-1].isdigit() else 0)
        elif len(parts) >= 4:
            read_name = parts[0]
            mapq = int(parts[3]) if parts[3].isdigit() else 0
        else:
            read_name = parts[0]
            mapq = 0
    except Exception:
        if parts:
            read_name = parts[0]
        mapq = 0
    return read_name, mapq, s

def extract_gene_id_from_readname(read_name: str):
    if not read_name:
        return None
    if "-" in read_name:
        return read_name.rsplit("-", 1)[0]
    return read_name

def dedup_merged_file(input_summary_path: str, out_sorted_folder: str):
    safe_mkdir(Path(out_sorted_folder))
    base = os.path.basename(input_summary_path)
    out_name = base.replace(".txt", "_sorted.txt")
    out_path = os.path.join(out_sorted_folder, out_name)
    best = {}
    with open(input_summary_path, "r") as inf:
        for line in inf:
            read_name, mapq, raw_line = parse_line_for_read_and_mapq(line)
            if not read_name:
                continue
            gene_id = extract_gene_id_from_readname(read_name)
            if not gene_id:
                continue
            prev = best.get(gene_id)
            if prev is None or mapq > prev[0]:
                best[gene_id] = (mapq, raw_line)
    with open(out_path, "w") as outf:
        for mapq, line in best.values():
            if not line.endswith("\n"):
                outf.write(line + "\n")
            else:
                outf.write(line)
    return out_path

def dedup_all_merged_parallel(merged_folder: str, dedup_out_folder: str, threads: int):
    safe_mkdir(Path(dedup_out_folder))
    merged_files = sorted([os.path.join(merged_folder, f) for f in os.listdir(merged_folder) if f.endswith(".txt") and not f.endswith("_sorted.txt")])
    results = []
    if not merged_files:
        return results
    with ProcessPoolExecutor(max_workers=threads) as exe:
        future_map = {exe.submit(dedup_merged_file, mf, dedup_out_folder): mf for mf in merged_files}
        for fut in as_completed(future_map):
            mf = future_map[fut]
            try:
                outp = fut.result()
                results.append(outp)
            except Exception:
                continue
    return results

def collect_gene_ids_from_sorted(sorted_folder: str):
    gene_ids = set()
    for fname in sorted(os.listdir(sorted_folder)):
        if not fname.endswith("_sorted.txt"):
            continue
        p = os.path.join(sorted_folder, fname)
        with open(p, "r") as inf:
            for line in inf:
                read_name, mapq, raw_line = parse_line_for_read_and_mapq(line)
                if not read_name:
                    continue
                gene_id = extract_gene_id_from_readname(read_name)
                if gene_id:
                    gene_ids.add(gene_id)
    return gene_ids

def filter_fasta_by_geneids(fasta_in_path: str, fasta_out_path: str, geneids_to_remove: set):
    kept = 0
    removed = 0
    with open(fasta_in_path, "r") as inf, open(fasta_out_path, "w") as outf:
        current_header = None
        seq_lines = []
        for line in inf:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    header_id = current_header[1:].split()[0]
                    if header_id not in geneids_to_remove:
                        outf.write(current_header + "\n")
                        if seq_lines:
                            outf.write("\n".join(seq_lines) + "\n")
                        else:
                            outf.write("\n")
                        kept += 1
                    else:
                        removed += 1
                current_header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_header is not None:
            header_id = current_header[1:].split()[0]
            if header_id not in geneids_to_remove:
                outf.write(current_header + "\n")
                if seq_lines:
                    outf.write("\n".join(seq_lines) + "\n")
                else:
                    outf.write("\n")
                kept += 1
            else:
                removed += 1
    return kept, removed

def filter_all_fastas_parallel(fasta_folder: str, out_folder: str, geneids_to_remove: set, threads: int):
    safe_mkdir(Path(out_folder))
    exts = (".fasta", ".fa", ".fna", ".ffn")
    fasta_files = sorted([f for f in os.listdir(fasta_folder) if any(f.endswith(e) for e in exts)])
    results = []
    if not fasta_files:
        return results
    with ThreadPoolExecutor(max_workers=threads) as exe:
        future_map = {}
        for f in fasta_files:
            in_path = os.path.join(fasta_folder, f)
            out_path = os.path.join(out_folder, f)
            future_map[exe.submit(filter_fasta_by_geneids, in_path, out_path, geneids_to_remove)] = (f, out_path)
        for fut in as_completed(future_map):
            f, outp = future_map[fut]
            try:
                kept, removed = fut.result()
                results.append((f, kept, removed, outp))
            except Exception:
                results.append((f, 0, 0, None))
    return results

# ---------------- Step2: rename/filter/prioritize ----------------
def read_tax_tsv(tsv_path: Path) -> List[Dict[str,str]]:
    rows = []
    with tsv_path.open('r', encoding='utf-8') as fh:
        first = fh.readline()
        if not first:
            return rows
        header_tokens = [t.strip().lower() for t in first.strip().split('\t')]
        known_cols = {'species_name','strain_id','taxonomy','taxid_string','genome_len','merged_alts','filename','file'}
        has_header = any(h in known_cols for h in header_tokens)
        fh.seek(0)
        reader = csv.reader(fh, delimiter='\t')
        if has_header:
            header = next(reader)
            header = [h.strip() for h in header]
            for rec in reader:
                if not rec or all(not x.strip() for x in rec):
                    continue
                if len(rec) < len(header):
                    rec += [''] * (len(header) - len(rec))
                d = {header[i]: rec[i].strip() for i in range(len(header))}
                row = {
                    'maybe_filename': d.get('filename') or d.get('file') or '',
                    'species_name': d.get('species_name') or d.get('species') or '',
                    'strain_id': d.get('strain_id') or d.get('strain') or '',
                    'taxonomy': d.get('taxonomy') or '',
                    'taxid_string': d.get('taxid_string') or d.get('taxid') or '',
                    'genome_len': d.get('genome_len') or '',
                    'merged_alts_raw': d.get('merged_alts') or '',
                    'raw_line': '\t'.join([d.get(h,'') for h in header])
                }
                rows.append(row)
        else:
            fh.seek(0)
            for rec in reader:
                if not rec or all(not x.strip() for x in rec):
                    continue
                if len(rec) >= 6:
                    row = {
                        'maybe_filename': rec[0].strip(),
                        'species_name': rec[1].strip(),
                        'strain_id': rec[2].strip(),
                        'taxonomy': rec[3].strip(),
                        'taxid_string': rec[4].strip(),
                        'genome_len': rec[5].strip(),
                        'merged_alts_raw': '\t'.join(rec[6:]).strip() if len(rec) > 6 else '',
                        'raw_line': '\t'.join([c.strip() for c in rec])
                    }
                elif len(rec) >= 5:
                    row = {
                        'maybe_filename': '',
                        'species_name': rec[0].strip(),
                        'strain_id': rec[1].strip(),
                        'taxonomy': rec[2].strip(),
                        'taxid_string': rec[3].strip(),
                        'genome_len': rec[4].strip(),
                        'merged_alts_raw': '\t'.join(rec[5:]).strip() if len(rec) > 5 else '',
                        'raw_line': '\t'.join([c.strip() for c in rec])
                    }
                else:
                    padded = rec + [''] * (5 - len(rec))
                    row = {
                        'maybe_filename': '',
                        'species_name': padded[0].strip(),
                        'strain_id': padded[1].strip(),
                        'taxonomy': padded[2].strip(),
                        'taxid_string': padded[3].strip(),
                        'genome_len': padded[4].strip(),
                        'merged_alts_raw': '',
                        'raw_line': '\t'.join([c.strip() for c in padded])
                    }
                rows.append(row)
    return rows

def find_candidate_files(fasta_dir: Path) -> List[Path]:
    exts = ('.ffn', '.fasta', '.fa', '.fna', '.fas')
    return [p for p in sorted(fasta_dir.iterdir()) if p.is_file() and p.suffix.lower() in exts]

def file_matches_token(p: Path, token: str) -> bool:
    if not token:
        return False
    return normalize(token) in normalize(p.name) or normalize(token) in normalize(p.stem)

def try_match_for_row(row: dict, files: List[Path], id_col: str|None) -> List[Path]:
    candidates: List[Path] = []
    if id_col and id_col in row and row[id_col].strip():
        tok = row[id_col].strip()
        candidates = [p for p in files if file_matches_token(p, tok)]
        if candidates:
            return candidates
    if row.get('strain_id','').strip():
        candidates = [p for p in files if file_matches_token(p, row['strain_id'])]
        if candidates:
            return candidates
    if row.get('taxid_string','').strip():
        tok = row['taxid_string'].split('|')[0]
        candidates = [p for p in files if file_matches_token(p, tok)]
        if candidates:
            return candidates
    if row.get('species_name','').strip():
        candidates = [p for p in files if file_matches_token(p, row['species_name'])]
        if candidates:
            return candidates
    core_matches = [p for p in files if re.search(r'core_clusters', p.name, flags=re.I)]
    if core_matches:
        return core_matches
    return []

def parse_fasta_to_records(path: Path) -> List[Tuple[str,str]]:
    recs = []
    if SeqIO is not None:
        for rec in SeqIO.parse(str(path), "fasta"):
            header = f">{rec.id}" + ((" " + rec.description[len(rec.id):].strip()) if rec.description and len(rec.description) > len(rec.id) else "")
            recs.append((header, str(rec.seq)))
    else:
        header = None
        seq_lines = []
        with path.open('r', encoding='utf-8') as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if header:
                        recs.append((header, "".join(seq_lines)))
                    header = line
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if header:
                recs.append((header, "".join(seq_lines)))
    return recs

def extract_function_from_header(header: str) -> str|None:
    h = header.lstrip('>')
    if ' ' not in h:
        return None
    after = h.split(' ',1)[1]
    cut_tokens = ["|UniRef90","|UniClust90","|UniRef50","|UniRef100"]
    func = after
    for t in cut_tokens:
        if t in func:
            func = func.split(t)[0]
            break
    func = func.strip()
    return func if func else None

def is_hypothetical(func: str|None, hypo_keyword: str="hypothetical protein") -> bool:
    if not func:
        return False
    return hypo_keyword.lower() in func.lower()

def load_bac120_set(tsv_path: Path) -> set:
    s = set()
    if not tsv_path or not tsv_path.exists():
        return s
    try:
        with tsv_path.open('r', encoding='utf-8') as fh:
            for line in fh:
                line=line.rstrip("\n")
                if not line:
                    continue
                parts=line.split("\t")
                if len(parts) >= 2:
                    name = parts[1].strip()
                    if name:
                        s.add(name)
    except Exception:
        pass
    return s

def filter_by_length(records: List[Tuple[str,str]], max_len: int) -> List[Tuple[str,str]]:
    return [(h,s) for (h,s) in records if len(s) <= max_len]

def prioritize_and_dedupe(records: List[Tuple[str,str]], bac120_set:set, max_keep:int, min_keep:int) -> List[Tuple[str,str]]:
    if not records:
        return []
    bac_map = {}
    other_map = {}
    for header, seq in records:
        func = extract_function_from_header(header) or ""
        if is_hypothetical(func):
            key = header
        else:
            key = func
        target = bac_map if (func in bac120_set) else other_map
        if key not in target:
            target[key] = (header, seq)
        else:
            if len(seq) > len(target[key][1]):
                target[key] = (header, seq)
    bac_list = list(bac_map.values())
    other_list = list(other_map.values())
    bac_list.sort(key=lambda x: len(x[1]), reverse=True)
    other_list.sort(key=lambda x: len(x[1]), reverse=True)
    combined = bac_list + other_list
    top = combined[:max_keep]
    if len(top) < min_keep:
        existing_headers = set(h for h,_ in top)
        all_sorted = sorted(records, key=lambda x: len(x[1]), reverse=True)
        for header, seq in all_sorted:
            if header in existing_headers:
                continue
            top.append((header, seq))
            existing_headers.add(header)
            if len(top) >= min_keep:
                break
    return top

def write_fasta_records(out_path: Path, records: List[Tuple[str,str]], width:int=60):
    with out_path.open('w', encoding='utf-8') as fh:
        for header, seq in records:
            fh.write(header + "\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i+width] + "\n")

# ---------------- main ----------------
def main():
    p = argparse.ArgumentParser(description="Simplified pipeline: dedup/filter mapped genes then rename+filter core clusters. Single outdir used for intermediate and final outputs.")
    p.add_argument("--merged_folder", required=True, help="Folder containing merged summary txt files (input)")
    p.add_argument("--core_fasta_dir", required=True, help="Folder containing original core-gene fasta files (input)")
    # alias for compatibility
    p.add_argument("--fasta_dir", help="alias for core_fasta_dir", default=None)
    p.add_argument("--mapped_fasta_dir", help="alias for core_fasta_dir", default=None)
    p.add_argument("--tax_infofile", required=True, help="Tax TSV file")
    p.add_argument("--outdir", required=True, help="Single output folder to hold dedup outputs, filtered fastas, and final species .ffn files")
    p.add_argument("--bac120_tsv", default="bac120-annotated-single-copy-gene.tsv", help="Bac120 TSV (default bac120-annotated-single-copy-gene.tsv in cwd)")
    p.add_argument("--threads", type=int, default=8, help="Threads for parallel operations (default 8)")
    p.add_argument("--max_len", type=int, default=6000, help="Max sequence length to keep (default 6000)")
    p.add_argument("--max_keep", type=int, default=200, help="Max genes to keep per file (default 200)")
    p.add_argument("--min_keep", type=int, default=10, help="Min genes to keep per file (default 10)")
    p.add_argument("--id-col", default=None, help="Optional column name in tax_infofile to match filenames")
    p.add_argument("--dryrun", action="store_true", help="Dry run (no files written)")
    args = p.parse_args()

    # unify core_fasta_dir from possible aliases
    core_fasta_dir = args.core_fasta_dir
    if args.fasta_dir:
        core_fasta_dir = args.fasta_dir
    if args.mapped_fasta_dir:
        core_fasta_dir = args.mapped_fasta_dir

    merged_folder = args.merged_folder
    outdir = Path(args.outdir).resolve()
    safe_mkdir(outdir)

    dedup_out = str(outdir)  # dedup outputs go to outdir
    mapped_fasta_out = str(outdir)  # filtered fastas also written into outdir (overwritten names preserved)
    mapped_fasta_dir = str(Path(core_fasta_dir).resolve())

    threads = max(1, args.threads)

    # Step1
    print("STEP1: dedup merged summaries")
    deduped = dedup_all_merged_parallel(merged_folder, dedup_out, threads)
    print("  Dedup produced", len(deduped), "sorted files (in outdir)")

    print("STEP1: collect mapped gene ids")
    geneids = collect_gene_ids_from_sorted(dedup_out)
    print("  Collected", len(geneids), "unique gene ids to remove from fasta files")

    print("STEP1: filter core-gene fasta files (remove mapped gene ids) --- input:", mapped_fasta_dir, " -> output:", mapped_fasta_out)
    if args.dryrun:
        print("[DRYRUN] Would filter fasta files to remove mapped gene ids (not performing in dryrun).")
    else:
        results = filter_all_fastas_parallel(mapped_fasta_dir, mapped_fasta_out, geneids, threads)
        for fname, kept, removed, outp in results:
            print("   ", fname, "kept:", kept, "removed:", removed)

    # Step2
    print("\nSTEP2: rename + length-filter + prioritize/dedupe (use tax_infofile and bac120)")
    tax_file = Path(args.tax_infofile).resolve()
    if not tax_file.exists():
        print("ERROR: tax_infofile not found:", tax_file, file=sys.stderr)
        sys.exit(2)
    rows = read_tax_tsv(tax_file)
    if not rows:
        print("ERROR: no rows parsed from tax_infofile", file=sys.stderr)
        sys.exit(2)

    # load filtered fastas from mapped_fasta_out (which equals outdir)
    fasta_input_dir = Path(mapped_fasta_out)
    files = find_candidate_files(fasta_input_dir)
    if not files:
        print("ERROR: no fasta/ffn files found in filtered folder:", fasta_input_dir, file=sys.stderr)
        sys.exit(2)

    # map rows to files
    mapped = []
    unmapped = []
    for i, row in enumerate(rows, start=1):
        cands = try_match_for_row(row, files, args.id_col)
        if not cands:
            unmapped.append((i, row))
            continue
        mapped.append((row, cands[0]))

    print(f"  files found: {len(files)}; tax rows: {len(rows)}; mapped rows: {len(mapped)}; unmapped rows: {len(unmapped)}")
    if unmapped:
        print("[WARN] example unmapped rows:")
        for u in unmapped[:10]:
            idx, r = u
            print("   ", idx, r.get('species_name'), r.get('strain_id'), r.get('maybe_filename'))

    bac120_path = Path(args.bac120_tsv).resolve()
    bac120_set = set()
    if bac120_path.exists():
        bac120_set = load_bac120_set(bac120_path)
        print("  Loaded", len(bac120_set), "bac120 names from", bac120_path)
    else:
        print("  Bac120 file not found at", bac120_path, "- continuing without bac120 prioritization")

    source_use_counts = {}
    for idx, (row, src_path) in enumerate(mapped, start=1):
        species = row.get('species_name') or f"unnamed_row_{idx}"
        safe_basename = normalize(species)
        count = source_use_counts.get(src_path, 0) + 1
        source_use_counts[src_path] = count
        dest_name = f"{safe_basename}.ffn" if count == 1 else f"{safe_basename}_{count}.ffn"
        dest_path = outdir / dest_name
        print(f"[STEP2] {idx}: {src_path.name} -> {dest_name}")

        try:
            recs = parse_fasta_to_records(src_path)
        except Exception as e:
            print(f"[ERROR] parse failed for {src_path}: {e}")
            continue

        before_n = len(recs)
        recs_lenflt = filter_by_length(recs, args.max_len)
        after_len_n = len(recs_lenflt)
        print(f"   length filter: {before_n} -> {after_len_n} (max_len={args.max_len})")
        top = prioritize_and_dedupe(recs_lenflt, bac120_set, args.max_keep, args.min_keep)
        print(f"   after prioritize/dedupe: {len(top)} (min_keep={args.min_keep}, max_keep={args.max_keep})")

        if args.dryrun:
            print(f"   [DRYRUN] Would write {len(top)} records to {dest_path}")
        else:
            try:
                write_fasta_records(dest_path, top, width=60)
                print(f"   Wrote {dest_path} (kept {len(top)})")
            except Exception as e:
                print(f"[ERROR] write failed for {dest_path}: {e}")

    print("\nAll done. Outputs are in:", outdir)
    if args.dryrun:
        print("[DRYRUN] No files were written.")

if __name__ == "__main__":
    main()
