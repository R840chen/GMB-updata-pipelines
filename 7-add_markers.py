#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gmb_add_coregenes.py

Add new core-gene ffn files into an existing GMB (MetaPhlAn-like) database.

Main steps:
 1) Read new tax TSV (user provided). Detect whether first column is filename or species_name.
 2) Validate that for each TSV row a corresponding .ffn exists (match by filename or by strain_id/species_name/taxid).
 3) Merge new tax TSV rows into existing GMB tax TSV (skip duplicates by strain_id if present).
 4) Load existing GMB pkl (bz2 pickle). If not present, start a new db skeleton.
 5) For each matched .ffn file, add its marker records into db['markers'], add db['taxonomy'] entries, and merged_taxon if present.
 6) Append .ffn contents into the GMB fasta file (Python-safe copy), unless --no-append.
 7) Run bowtie2-build to rebuild index using --bt2-threads (unless --no-build).
 8) Write out new pkl (orig_pkl -> orig_pkl_NEW.pkl) and merged tax TSV (orig_tax -> orig_tax_merged.tsv) plus a log file.

Usage example:
  python3 gmb_add_coregenes.py \
    --ffn_dir /path/to/ffn_dir \
    --new_tax_tsv /path/to/new_tax.tsv \
    --gmb_tax_tsv /path/to/gmb_tax.tsv \
    --gmb_fasta /path/to/gmb_db.fasta \
    --gmb_pkl /path/to/gmb_db.pkl \
    --out_dir /path/to/output_dir \
    --bt2_threads 16

Notes:
 - The expected canonical tax TSV columns used by GMB are:
     species_name, strain_id, taxonomy, taxid_string, genome_len [, merged_alts]
   But the new_tax_tsv may optionally include a first column that is the ffn filename.
 - The script is conservative on marker id conflicts: if a marker id already exists in db
   but with a different taxon, the script will log a conflict and will NOT overwrite.
 - Biopython is required (SeqIO). The script checks import and will exit with a helpful message if missing.
"""

from __future__ import annotations
import argparse
import os
import sys
import bz2
import pickle
import shutil
import subprocess
from pathlib import Path
import csv
import re
from datetime import datetime

# Try import SeqIO and fail fast with hint
try:
    from Bio import SeqIO
except Exception as e:
    print("ERROR: Biopython is required (pip install biopython or conda install -c bioconda biopython).")
    raise

# ---------- helpers ----------
def safe_mkdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def normalize_token(s: str) -> str:
    return re.sub(r'\s+', '_', s.strip().lower()) if s else ''

def read_tax_tsv_flexible(path: Path):
    """
    Read a tax TSV flexibly. Return list of rows (ordered) where each row is dict with keys:
      maybe_filename, species_name, strain_id, taxonomy, taxid_string, genome_len, merged_alts_raw, raw_line
    The function supports files with or without header. If header exists and column names are present,
    columns are mapped by name. Otherwise the code assumes at least 5 useful columns in order:
      species_name, strain_id, taxonomy, taxid_string, genome_len [, merged_alts]
    If a first column looks like a filename (endswith .ffn/.fa/.fasta/.fna) or if its values match ffn files later,
    consumer code can use the 'maybe_filename' field.
    """
    rows = []
    with path.open('r', encoding='utf-8') as fh:
        sample = fh.readline()
        if not sample:
            return rows
        # detect header by checking tokens against known column names
        fh.seek(0)
        reader = csv.reader(fh, delimiter='\t')
        first = next(reader)
        header_names = [c.strip().lower() for c in first]
        known = {'species_name','strain_id','taxonomy','taxid_string','genome_len','merged_alts','filename'}
        has_header = any(h in known for h in header_names)
        fh.seek(0)
        reader = csv.reader(fh, delimiter='\t')
        if has_header:
            hdr = next(reader)
            hdr = [h.strip() for h in hdr]
            for rec in reader:
                if not rec or all(not x.strip() for x in rec):
                    continue
                # pad
                if len(rec) < len(hdr):
                    rec += [''] * (len(hdr) - len(rec))
                d = {hdr[i]: rec[i].strip() for i in range(len(hdr))}
                # map into canonical keys
                row = {
                    'maybe_filename': d.get('filename') or d.get('file') or d.get(hdr[0]) if ('filename' in hdr or 'file' in hdr) else None,
                    'species_name': d.get('species_name') or d.get('species') or '',
                    'strain_id': d.get('strain_id') or d.get('strain') or '',
                    'taxonomy': d.get('taxonomy') or '',
                    'taxid_string': d.get('taxid_string') or d.get('taxid') or '',
                    'genome_len': d.get('genome_len') or '',
                    'merged_alts_raw': d.get('merged_alts') or ''
                }
                row['raw_line'] = '\t'.join([d.get(k,'') for k in hdr])
                rows.append(row)
        else:
            # no header, assume canonical order
            fh.seek(0)
            for rec in reader:
                if not rec or all(not x.strip() for x in rec):
                    continue
                # if len >=6 maybe first col is filename: heuristics below (we still map)
                if len(rec) >= 6:
                    # assume: maybe_filename, species_name, strain_id, taxonomy, taxid_string, genome_len [, merged_alts]
                    row = {
                        'maybe_filename': rec[0].strip(),
                        'species_name': rec[1].strip(),
                        'strain_id': rec[2].strip(),
                        'taxonomy': rec[3].strip(),
                        'taxid_string': rec[4].strip(),
                        'genome_len': rec[5].strip(),
                        'merged_alts_raw': '\t'.join(rec[6:]).strip() if len(rec) > 6 else ''
                    }
                elif len(rec) >= 5:
                    # assume: species_name, strain_id, taxonomy, taxid_string, genome_len [, merged_alts]
                    row = {
                        'maybe_filename': '',
                        'species_name': rec[0].strip(),
                        'strain_id': rec[1].strip(),
                        'taxonomy': rec[2].strip(),
                        'taxid_string': rec[3].strip(),
                        'genome_len': rec[4].strip(),
                        'merged_alts_raw': '\t'.join(rec[5:]).strip() if len(rec) > 5 else ''
                    }
                else:
                    # too few columns, create best-effort mapping
                    padded = rec + [''] * (5 - len(rec))
                    row = {
                        'maybe_filename': '',
                        'species_name': padded[0].strip(),
                        'strain_id': padded[1].strip(),
                        'taxonomy': padded[2].strip(),
                        'taxid_string': padded[3].strip(),
                        'genome_len': padded[4].strip(),
                        'merged_alts_raw': ''
                    }
                row['raw_line'] = '\t'.join(rec)
                rows.append(row)
    return rows

def parse_merged_alts(field: str):
    """Parse merged_alts_raw into list of triples (tax, taxids, count_int or None)"""
    out = []
    if not field:
        return out
    for item in field.split('|||'):
        item = item.strip()
        if not item:
            continue
        parts = item.split(':::')
        if len(parts) != 3:
            # malformed, ignore
            continue
        a = parts[0].strip()
        b = parts[1].strip()
        try:
            c = int(parts[2].strip())
        except Exception:
            c = None
        out.append((a, b, c))
    return out

def load_existing_db(pkl_path: Path):
    """Load existing bz2 pkl or return skeleton db"""
    if not pkl_path.exists() or pkl_path.stat().st_size == 0:
        return {"markers": {}, "taxonomy": {}, "merged_taxon": {}}
    try:
        with bz2.open(str(pkl_path), "rb") as fh:
            db = pickle.load(fh)
            if not isinstance(db, dict) or 'markers' not in db or 'taxonomy' not in db:
                return {"markers": {}, "taxonomy": {}, "merged_taxon": {}}
            if 'merged_taxon' not in db:
                db['merged_taxon'] = {}
            return db
    except Exception as e:
        print(f"[WARN] Failed to load existing pkl {pkl_path}: {e}. Starting from empty db.")
        return {"markers": {}, "taxonomy": {}, "merged_taxon": {}}

def append_fastas_python(ffn_paths: list[Path], db_fasta: Path):
    """
    Append multiple ffn files into db_fasta (safe with spaces).
    Returns tuple (copied_count, total_bytes)
    """
    copied = 0
    total_bytes = 0
    db_fasta.parent.mkdir(parents=True, exist_ok=True)
    # ensure file exists
    if not db_fasta.exists():
        db_fasta.touch()
    with db_fasta.open('ab') as outfh:
        for p in ffn_paths:
            # skip if same file as db_fasta
            if p.resolve() == db_fasta.resolve():
                continue
            try:
                with p.open('rb') as inf:
                    shutil.copyfileobj(inf, outfh)
                    total_bytes += p.stat().st_size
                copied += 1
            except Exception as e:
                print(f"[WARN] Failed to append {p}: {e}")
    return copied, total_bytes

def write_new_pkl(db: dict, orig_pkl: Path) -> Path:
    out_pkl = orig_pkl.with_name(orig_pkl.stem + "_NEW.pkl")
    with bz2.BZ2File(str(out_pkl), "w") as of:
        pickle.dump(db, of, pickle.HIGHEST_PROTOCOL)
    return out_pkl

def write_merged_tax_tsv(existing_tax: Path, merged_lines: list[str], out_path: Path):
    """
    existing_tax: original gmb tax tsv (may be used to preserve header)
    merged_lines: list of new rows (raw lines) that were appended
    out_path: write merged tsv file
    Behavior: simple append of merged_lines to existing_tax content (does not dedupe by text).
    """
    # if existing_tax exists, copy to out_path, then append new lines
    if existing_tax.exists():
        shutil.copy2(str(existing_tax), str(out_path))
    else:
        # create out_path with a simple header if none exists
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open('w', encoding='utf-8') as fh:
            fh.write("")  # empty
    if merged_lines:
        with out_path.open('a', encoding='utf-8') as fh:
            for line in merged_lines:
                fh.write(line.rstrip("\n") + "\n")
    return out_path

def run_bowtie2_build(fasta_path: Path, out_prefix: str, threads: int, large_index: bool=True):
    cmd = ["bowtie2-build"]
    if threads and isinstance(threads, int) and threads > 0:
        cmd += ["--threads", str(threads)]
    if large_index:
        cmd += ["--large-index"]
    cmd += [str(fasta_path), out_prefix]
    print("[INFO] Running:", " ".join(cmd))
    try:
        res = subprocess.run(cmd, check=False)
        return res.returncode
    except Exception as e:
        print(f"[ERROR] bowtie2-build failed to start: {e}")
        return 255

# ---------- main flow ----------
def main():
    p = argparse.ArgumentParser(description="Add core gene ffn files into GMB (fasta+pkl+tax tsv merge) and rebuild bowtie2 index.")
    p.add_argument("--ffn_dir", required=True, help="Directory containing new core gene .ffn files")
    p.add_argument("--new_tax_tsv", required=True, help="User-supplied tax TSV for the new species (may include filename as first column)")
    p.add_argument("--gmb_tax_tsv", required=True, help="Existing GMB tax TSV used to build the current DB")
    p.add_argument("--gmb_fasta", required=True, help="Existing GMB database fasta to append to")
    p.add_argument("--gmb_pkl", required=True, help="Existing GMB pkl file (bz2 pickle) to update")
    p.add_argument("--out_dir", required=True, help="Output directory to place new pkl, merged tax tsv and logs")
    p.add_argument("--bt2_threads", type=int, default=8, help="Threads for bowtie2-build (default 8)")
    p.add_argument("--no-append", action="store_true", help="Do not append .ffn to gmb_fasta (for testing)")
    p.add_argument("--no-build", action="store_true", help="Do not run bowtie2-build (for testing)")
    p.add_argument("--dryrun", action="store_true", help="Dry run: validate and print actions but do not modify files")
    args = p.parse_args()

    ffn_dir = Path(args.ffn_dir).resolve()
    new_tax_tsv = Path(args.new_tax_tsv).resolve()
    gmb_tax_tsv = Path(args.gmb_tax_tsv).resolve()
    gmb_fasta = Path(args.gmb_fasta).resolve()
    gmb_pkl = Path(args.gmb_pkl).resolve()
    out_dir = Path(args.out_dir).resolve()
    safe_mkdir(out_dir)

    log = {
        "start_time": datetime.now().isoformat(),
        "warnings": [],
        "errors": [],
        "mapped_rows": [],
        "unmapped_rows": [],
        "marker_conflicts": [],
        "added_markers_total": 0,
        "appended_files": [],
        "appended_bytes": 0,
        "bt2_returncode": None
    }

    # validate inputs
    if not ffn_dir.exists() or not ffn_dir.is_dir():
        print("ERROR: ffn_dir not found:", ffn_dir); sys.exit(2)
    if not new_tax_tsv.exists() or not new_tax_tsv.is_file():
        print("ERROR: new_tax_tsv not found:", new_tax_tsv); sys.exit(2)
    if not gmb_tax_tsv.exists() or not gmb_tax_tsv.is_file():
        print("ERROR: gmb_tax_tsv not found:", gmb_tax_tsv); sys.exit(2)
    # gmb_fasta and gmb_pkl may or may not exist; we handle both cases

    # read new tax tsv flexibly
    new_rows = read_tax_tsv_flexible(new_tax_tsv)
    if not new_rows:
        print("ERROR: no rows parsed from new_tax_tsv"); sys.exit(2)

    # list ffn files
    ffn_files = sorted([p for p in ffn_dir.iterdir() if p.is_file() and p.suffix.lower() in ('.ffn','.fa','.fasta','.fna')])
    ffn_stems = {p.name: p for p in ffn_files}
    ffn_stem_noext = {p.stem: p for p in ffn_files}

    # Determine mapping between rows and ffn files
    mapped = []      # list of tuples (rowdict, ffn_path)
    unmapped = []    # list of rowdicts
    for idx, row in enumerate(new_rows, start=1):
        # build canonical keys
        maybe_fname = row.get('maybe_filename','') or ''
        species = row.get('species_name','').strip()
        strain = row.get('strain_id','').strip()
        taxid = row.get('taxid_string','').strip()
        # first check maybe_filename exact match
        chosen = None
        if maybe_fname:
            # try raw match (with extension)
            if maybe_fname in ffn_stems:
                chosen = ffn_stems[maybe_fname]
            else:
                # try if user provided filename without extension
                if maybe_fname in ffn_stem_noext:
                    chosen = ffn_stem_noext[maybe_fname]
        # then try strain id in filename
        if chosen is None and strain:
            for key,p in ffn_stems.items():
                if strain in key:
                    chosen = p
                    break
        # then try species normalized in filename
        if chosen is None and species:
            ns = normalize_token(species)
            for key,p in ffn_stems.items():
                if ns in normalize_token(key) or ns in normalize_token(p.stem):
                    chosen = p
                    break
        # then try taxid token
        if chosen is None and taxid:
            t0 = taxid.split('|')[0]
            for key,p in ffn_stems.items():
                if t0 in key or t0 in p.stem:
                    chosen = p
                    break
        # fallback: if only one ffn file and one row, assign it
        if chosen is None and len(ffn_files)==1 and len(new_rows)==1:
            chosen = ffn_files[0]
        if chosen:
            mapped.append((row, chosen))
            log['mapped_rows'].append((idx, chosen.name))
        else:
            unmapped.append((idx, row))
            log['unmapped_rows'].append((idx, species, strain, maybe_fname))

    # report mapping summary
    print(f"Found {len(ffn_files)} ffn files in {ffn_dir}")
    print(f"Rows in new_tax_tsv: {len(new_rows)}; mapped: {len(mapped)}; unmapped: {len(unmapped)}")
    if unmapped:
        print("[WARN] Some rows did not map to any ffn file. Examples:")
        for u in unmapped[:10]:
            print("  row", u[0], "species:", u[1].get('species_name'), "strain:", u[1].get('strain_id'), "maybe_filename:", u[1].get('maybe_filename'))

    # Merge tax TSVs: skip duplicates by strain_id if present
    existing_strain_ids = set()
    # read existing gmb_tax_tsv strains
    with gmb_tax_tsv.open('r', encoding='utf-8') as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                existing_strain_ids.add(cols[1].strip())
    # prepare merged lines to append (raw_line from new_rows)
    merged_lines_to_append = []
    appended_count = 0
    skipped_duplicates = 0
    for row in new_rows:
        sid = row.get('strain_id','').strip()
        raw = row.get('raw_line') or '\t'.join([
            row.get('maybe_filename',''), row.get('species_name',''), row.get('strain_id',''),
            row.get('taxonomy',''), row.get('taxid_string',''), row.get('genome_len','')
        ])
        if sid and sid in existing_strain_ids:
            skipped_duplicates += 1
            continue
        merged_lines_to_append.append(raw)
        appended_count += 1

    merged_tax_out = out_dir / (gmb_tax_tsv.stem + "_merged.tsv")
    if args.dryrun:
        print("[DRYRUN] Would append", appended_count, "rows to", gmb_tax_tsv, "(skipping", skipped_duplicates, "duplicates by strain_id)")
    else:
        write_merged_tax_tsv(gmb_tax_tsv, merged_lines_to_append, merged_tax_out)
        print("Wrote merged tax tsv to:", merged_tax_out)

    # Load db
    db = load_existing_db(gmb_pkl)
    added_markers_total = 0

    # Process mapped ffn files: add taxonomy entries, merged_taxon, and markers
    for row, fpath in mapped:
        species = row.get('species_name','').strip()
        strain = row.get('strain_id','').strip() or fpath.stem  # fallback to filename if no strain_id
        taxonomy = row.get('taxonomy','').strip()
        taxids = row.get('taxid_string','').strip()
        genome_len = row.get('genome_len','').strip()
        merged_alts_raw = row.get('merged_alts_raw','').strip()
        merged_alts = parse_merged_alts(merged_alts_raw)

        # tmp tax key
        tmp_tax_key = taxonomy + "|t__" + strain

        # add taxonomy if absent
        if tmp_tax_key not in db.get("taxonomy", {}):
            try:
                glen = int(genome_len) if genome_len and genome_len.isdigit() else None
            except Exception:
                glen = None
            db.setdefault("taxonomy", {})[tmp_tax_key] = (taxids, glen)
        else:
            # overlap
            log['warnings'].append(f"taxonomy key already present: {tmp_tax_key}")

        # add merged_alts into db['merged_taxon']
        if merged_alts:
            key = (tmp_tax_key, taxids)
            db.setdefault('merged_taxon', {})
            if key not in db['merged_taxon']:
                db['merged_taxon'][key] = []
            existing = set((a[0], a[1]) for a in db['merged_taxon'][key])
            for alt_tax, alt_taxids, alt_count in merged_alts:
                if (alt_tax, alt_taxids) in existing:
                    log.setdefault('merged_dups', []).append((str(fpath.name), alt_tax, alt_taxids))
                    continue
                db['merged_taxon'][key].append((alt_tax, alt_taxids, alt_count))
                log.setdefault('merged_added', []).append((str(fpath.name), alt_tax, alt_taxids, alt_count))

        # add markers from ffn into db
        added_local = 0
        conflict_local = 0
        try:
            for rec in SeqIO.parse(str(fpath), "fasta"):
                gid = rec.id.split()[0]
                if gid in db.get('markers', {}):
                    existing_tax = db['markers'][gid].get('taxon')
                    if existing_tax == tmp_tax_key:
                        continue
                    # conflict: log and skip
                    log.setdefault('marker_conflicts', []).append((gid, str(fpath.name), existing_tax, tmp_tax_key))
                    conflict_local += 1
                    continue
                # add marker
                seq_len = len(str(rec.seq).strip())
                db.setdefault('markers', {})[gid] = {
                    'clade': 't__' + strain,
                    'ext': [],
                    'len': int(seq_len),
                    'taxon': tmp_tax_key
                }
                added_local += 1
        except Exception as e:
            log['errors'].append(f"Failed to parse ffn {fpath}: {e}")
            continue

        added_markers_total += added_local
        log.setdefault('per_file', []).append((str(fpath.name), strain, added_local, conflict_local))
        print(f"Processed {fpath.name}: added {added_local} markers, conflicts {conflict_local}")

    log['added_markers_total'] = added_markers_total

    # Append ffn files into gmb_fasta (python-safe)
    ffn_paths_to_append = [p for _, p in mapped]
    if args.no_append:
        print("[INFO] --no-append set: skip appending ffn to gmb_fasta")
    else:
        if args.dryrun:
            print("[DRYRUN] Would append", len(ffn_paths_to_append), "files to", gmb_fasta)
        else:
            copied, total_bytes = append_fastas_python(ffn_paths_to_append, gmb_fasta)
            log['appended_files'] = [p.name for p in ffn_paths_to_append]
            log['appended_bytes'] = total_bytes
            print(f"Appended {copied} files ({total_bytes} bytes) to {gmb_fasta}")

    # Rebuild bowtie2 index
    if args.no_build:
        print("[INFO] --no-build set: skip bowtie2-build")
    else:
        # out prefix: gmb_fasta path without extension with _NEW suffix
        out_prefix = str(gmb_fasta).replace('.fasta', '_NEW')
        if args.dryrun:
            print("[DRYRUN] Would run bowtie2-build on", gmb_fasta, "-> prefix", out_prefix, "threads", args.bt2_threads)
        else:
            rc = run_bowtie2_build(gmb_fasta, out_prefix, args.bt2_threads, large_index=True)
            log['bt2_returncode'] = rc
            if rc != 0:
                print(f"[WARN] bowtie2-build returned nonzero code {rc}")

    # Write new pkl
    if args.dryrun:
        print("[DRYRUN] Would write new pkl replacing", gmb_pkl, "->", gmb_pkl.with_name(gmb_pkl.stem + "_NEW.pkl"))
    else:
        out_pkl = write_new_pkl(db, gmb_pkl)
        print("Wrote new pkl to:", out_pkl)

    # Write merged tax tsv already done above if not dryrun; if dryrun, report what would be written
    if args.dryrun:
        print("[DRYRUN] Merged tax tsv would be written to", out_dir / (gmb_tax_tsv.stem + "_merged.tsv"))
    else:
        print("Merged tax tsv written to:", merged_tax_out)

    # write log file
    logfile = out_dir / (gmb_pkl.stem + "_addmarkers.log")
    with logfile.open('w', encoding='utf-8') as lf:
        lf.write(f"Run at: {datetime.now().isoformat()}\n")
        lf.write(f"Input new_tax_tsv: {new_tax_tsv}\n")
        lf.write(f"ffn_dir: {ffn_dir}\n")
        lf.write("\nSummary:\n")
        lf.write(f"  total ffn files found: {len(ffn_files)}\n")
        lf.write(f"  rows in new_tax_tsv: {len(new_rows)}\n")
        lf.write(f"  mapped rows: {len(mapped)}\n")
        lf.write(f"  unmapped rows: {len(unmapped)}\n")
        lf.write(f"  markers added total: {added_markers_total}\n")
        lf.write(f"  appended bytes: {log.get('appended_bytes',0)}\n")
        lf.write(f"  bowtie2 return code: {log.get('bt2_returncode')}\n\n")
        if log.get('per_file'):
            lf.write("Per-file details (name, strain, added, conflicts):\n")
            for rec in log['per_file']:
                lf.write("  " + "\t".join(map(str,rec)) + "\n")
        if log.get('marker_conflicts'):
            lf.write("\nMarker conflicts (first 50):\n")
            for rec in log['marker_conflicts'][:50]:
                lf.write("  " + "\t".join(map(str,rec)) + "\n")
        if log.get('merged_added'):
            lf.write("\nMerged alternatives added (sample):\n")
            for rec in log['merged_added'][:50]:
                lf.write("  " + "\t".join(map(str,rec)) + "\n")
        if log.get('merged_dups'):
            lf.write("\nMerged alternatives duplicates skipped (sample):\n")
            for rec in log['merged_dups'][:50]:
                lf.write("  " + "\t".join(map(str,rec)) + "\n")
        if log.get('errors'):
            lf.write("\nErrors:\n")
            for e in log['errors']:
                lf.write("  " + str(e) + "\n")
        if log.get('warnings'):
            lf.write("\nWarnings:\n")
            for w in log['warnings']:
                lf.write("  " + str(w) + "\n")
    print("Wrote log to:", logfile)
    print("Done.")

if __name__ == "__main__":
    main()
