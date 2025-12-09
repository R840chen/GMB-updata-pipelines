#!/usr/bin/env python3
"""
merged_drep_pipeline.py

This script performs the following steps:
 1. Copy new genome files from --new-dir into --input-dir (convert extensions to .fna), avoiding name collisions.
 2. Build a genome list file for dRep (one absolute path per line) at the location specified by --genome-list or default to <input-dir>/genome_list.txt.
 3. Run dRep compare using the provided --output-dir and thread count.
 4. (Optional) After dRep completes, parse dRep output data_tables/Cdb.csv (or CdbF.csv) to find "unique" primary_clusters that are represented only by genomes from the new set, or clusters with a single member. Copy those genome files to --unique-dest and optionally write their paths to --unique-list.

Usage example:
  python merged_drep_pipeline.py \
    --input-dir /path/to/main/genomes \
    --new-dir /path/to/new/genomes \
    --output-dir /path/to/drep-output \
    --threads 64 \
    --genome-list /path/to/genome_list.txt \
    --extract-unique --unique-dest /path/to/unique-output \
    --unique-list /path/to/unique_genomes.txt

Notes:
 - If --genome-list is not given, genome_list.txt will be written to <input-dir>/genome_list.txt
 - If --extract-unique is used, the script will look for data_tables/Cdb.csv (preferred) or CdbF.csv in the dRep output directory after dRep finishes.
 - The script treats files in --new-dir as the "new set" even if the files were copied into --input-dir earlier.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
import pandas as pd

def find_genome_files(root, ext):
    ext = ext.lower().lstrip('.')
    matches = []
    for p in Path(root).rglob(f'*.{ext}'):
        if p.is_file():
            matches.append(p)
    return matches

def copy_new_files(new_dir, input_dir, suffix):
    new_dir = Path(new_dir)
    input_dir = Path(input_dir)
    candidates = find_genome_files(new_dir, suffix)
    if not candidates:
        raise SystemExit(f"No files with extension .{suffix} found under {new_dir}")
    copied = []
    for src in candidates:
        dest_name = src.stem + '.fna'
        dest_path = input_dir.joinpath(dest_name)
        i = 1
        while dest_path.exists():
            dest_path = input_dir.joinpath(f'{src.stem}_copy{i}.fna')
            i += 1
        shutil.copy2(src, dest_path)
        copied.append(dest_path)
    return copied

def write_genome_list(input_dir, genome_list_path):
    all_fna = find_genome_files(input_dir, 'fna')
    if not all_fna:
        raise SystemExit(f'Error: no .fna files found under {input_dir}')
    with open(genome_list_path, 'w') as fh:
        for p in sorted(all_fna):
            fh.write(str(p) + '\n')
    return genome_list_path, len(all_fna)

def run_drep(output_dir, genome_list_path, threads):
    drep_cmd = [
        'dRep', 'compare', str(output_dir),
        '-p', str(threads),
        '-g', str(genome_list_path),
        '-pa', '0.95',
        '--multiround_primary_clustering',
        '--primary_chunksize', '3000',
        '-d',
        '-nc', '0.6'
    ]
    print('Constructed dRep command:')
    print(' '.join(drep_cmd))
    try:
        proc = subprocess.run(drep_cmd, check=False)
        return proc.returncode
    except FileNotFoundError:
        print('Error: dRep executable not found in PATH. Please ensure dRep is installed and available.', file=sys.stderr)
        return 127

def find_cdb_file(drep_out):
    dt_dir = Path(drep_out) / 'data_tables'
    candidates = [dt_dir / 'Cdb.csv', dt_dir / 'CdbF.csv']
    for c in candidates:
        if c.exists():
            return c
    return None

def read_cdb(cdb_path):
    df = pd.read_csv(cdb_path, dtype=str, low_memory=False)
    df.columns = [c.strip() for c in df.columns]
    if 'primary_cluster' not in df.columns or 'genome' not in df.columns:
        raise ValueError(f"Cdb file {cdb_path} does not contain required columns 'primary_cluster' and 'genome'")
    return df[['primary_cluster', 'genome']].copy()

def build_cluster_map(df):
    cluster_map = {}
    for _, row in df.iterrows():
        pc = str(row['primary_cluster']).strip()
        g = str(row['genome']).strip()
        cluster_map.setdefault(pc, []).append(g)
    return cluster_map

def match_genome_entry(genome_entry: str, file_path: Path):
    g = genome_entry.strip()
    if g == file_path.name:
        return True
    if g.endswith(file_path.name):
        return True
    if g == file_path.stem or g.endswith(file_path.stem):
        return True
    for ext in ('fna','fasta'):
        if g.endswith(f"{file_path.stem}.{ext}"):
            return True
    if file_path.stem in g:
        return True
    return False

def extract_unique_newgenomes(drep_out, new_dir, input_dir, unique_dest, dry_run=False, verbose=False):
    cdb = find_cdb_file(drep_out)
    if cdb is None:
        raise SystemExit(f"Error: could not find Cdb.csv or CdbF.csv in {drep_out}/data_tables")
    if verbose:
        print(f"Using Cdb file: {cdb}")
    df = read_cdb(cdb)
    cluster_map = build_cluster_map(df)
    genome_to_cluster = {}
    for pc, members in cluster_map.items():
        for g in members:
            genome_to_cluster.setdefault(g, []).append(pc)
    new_files = sorted([p for p in Path(new_dir).rglob('*') if p.is_file() and p.suffix.lower() in ('.fna','.fasta')])
    if not new_files:
        print(f"Warning: no .fna/.fasta in new-dir {new_dir}")
        return []
    cdb_genomes = list(genome_to_cluster.keys())
    unique_copied = []
    for newf in new_files:
        matched_entries = [g for g in cdb_genomes if match_genome_entry(g, newf)]
        if not matched_entries:
            if verbose:
                print(f"[NOT FOUND IN CDB] {newf}")
            continue
        unique_for_this_file = False
        for g_entry in matched_entries:
            clusters = genome_to_cluster.get(g_entry, [])
            for pc in clusters:
                members = cluster_map.get(pc, [])
                members_count = len(members)
                member_is_from_newdir = []
                for m in members:
                    matched_in_new = any(match_genome_entry(m, nf) for nf in new_files)
                    member_is_from_newdir.append(matched_in_new)
                all_members_from_newdir = all(member_is_from_newdir)
                if members_count == 1 or all_members_from_newdir:
                    unique_for_this_file = True
                    if verbose:
                        print(f"[UNIQUE] {newf.name} -> cluster {pc} (members={members_count}, all_from_newdir={all_members_from_newdir})")
                    break
            if unique_for_this_file:
                break
        if not unique_for_this_file:
            if verbose:
                print(f"[NOT UNIQUE] {newf.name}")
            continue
        candidate_src = None
        p_in_input = input_dir / newf.name
        if p_in_input.exists():
            candidate_src = p_in_input
        else:
            matches_in_input = list(Path(input_dir).rglob(f"{newf.stem}*"))
            if matches_in_input:
                preferred = None
                for m in matches_in_input:
                    if m.suffix.lower() == '.fna':
                        preferred = m
                        break
                candidate_src = preferred or matches_in_input[0]
        if candidate_src is None:
            candidate_src = newf
        dest = Path(unique_dest).resolve()
        dest.mkdir(parents=True, exist_ok=True)
        dest_path = dest / candidate_src.name
        if dry_run:
            print(f"[DRY-RUN] Would copy: {candidate_src} -> {dest_path}")
            unique_copied.append((candidate_src, dest_path))
        else:
            shutil.copy2(candidate_src, dest_path)
            unique_copied.append((candidate_src, dest_path))
            if verbose:
                print(f"Copied: {candidate_src} -> {dest_path}")
    return unique_copied

def main():
    parser = argparse.ArgumentParser(description='Prepare genomes, run dRep, and optionally extract unique new genomes')
    parser.add_argument('--input-dir', required=True, help='Main genomes directory (existing genomes, where new ones will be copied)')
    parser.add_argument('--new-dir', required=True, help='Directory containing new genome files to copy into the main directory')
    parser.add_argument('--suffix', default='fna', choices=['fna','fasta'], help='Suffix of the new genome files in new-dir (default: fna)')
    parser.add_argument('--output-dir', required=True, help='Output directory for dRep results (work directory). Must be empty or non-existent prior to running dRep')
    parser.add_argument('--threads', type=int, default=64, help='Number of threads for dRep')
    parser.add_argument('--genome-list', default=None, help='Path to write the genome list file for dRep. Default: <input-dir>/genome_list.txt')
    parser.add_argument('--dry-run', action='store_true', help='If set, do not run dRep; only copy files and write the genome list')
    parser.add_argument('--extract-unique', action='store_true', help='After dRep, extract new genomes that occupy unique primary clusters')
    parser.add_argument('--unique-dest', default=None, help='Where to copy unique new genomes. Default: <input-dir>/unique-new-genomes')
    parser.add_argument('--unique-list', default=None, help='Path to write the list of unique genome files. Default: <unique-dest>/unique_genomes.txt')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    args = parser.parse_args()

    input_dir = Path(args.input_dir).expanduser().resolve()
    new_dir = Path(args.new_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    suffix = args.suffix.lower().lstrip('.')
    threads = args.threads

    genome_list_path = Path(args.genome_list).expanduser().resolve() if args.genome_list else input_dir.joinpath('genome_list.txt')

    if not input_dir.is_dir():
        print(f'Error: input-dir does not exist or is not a directory: {input_dir}', file=sys.stderr)
        sys.exit(1)
    if not new_dir.is_dir():
        print(f'Error: new-dir does not exist or is not a directory: {new_dir}', file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    copied = copy_new_files(new_dir, input_dir, suffix)
    print(f'Copied {len(copied)} files into {input_dir} (converted to .fna).')

    genome_list_path, total = write_genome_list(input_dir, genome_list_path)
    print(f'Wrote genome list to {genome_list_path} with {total} entries.')

    if args.dry_run:
        print('Dry run requested; exiting before running dRep.')
        sys.exit(0)

    rc = run_drep(output_dir, genome_list_path, threads)
    if rc != 0:
        print(f'dRep exited with return code {rc}', file=sys.stderr)
        sys.exit(rc)

    if args.extract_unique:
        unique_dest = Path(args.unique_dest) if args.unique_dest else input_dir.joinpath('unique-new-genomes')
        unique = extract_unique_newgenomes(output_dir, new_dir, input_dir, unique_dest, dry_run=False, verbose=args.verbose)
        print('Extraction finished. Unique copied count:', len(unique))
        if unique:
            for s,d in unique:
                print(f'  {s} -> {d}')

        # Write unique genome list
        unique_list_path = Path(args.unique_list).expanduser().resolve() if args.unique_list else unique_dest.joinpath('unique_genomes.txt')
        with open(unique_list_path, 'w') as fh:
            for _, dest_path in unique:
                fh.write(str(dest_path) + '\n')
        print(f'Unique genome list written to: {unique_list_path}')

    print('All done.')

if __name__ == '__main__':
    main()
