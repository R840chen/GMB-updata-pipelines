#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Add new marker info to MetaPhlAn pkl database (Scheme A: prefer combined ids).
- Prefer combined token species_X_cluster_Y extracted from filenames.
- Parse optional 6th column `merged_alts` in tax_infofile:
    alt_taxonomy:::alt_taxids:::alt_count|||alt2_taxonomy:::alt2_taxids:::alt2_count
- Add entries into db['merged_taxon'] using key (taxonomy_with_t__, taxids).
- Provide --no-append and --no-build flags to avoid touching fasta/index during testing.
"""
import os
import re
import bz2
import pickle
import argparse
from Bio import SeqIO
from datetime import datetime

def get_args():
    parser = argparse.ArgumentParser(description="Add new marker info to MetaPhlAn pkl database (prefer combined ids).")
    parser.add_argument("-t", "--tax_infofile", required=True, type=str,
                        help="Taxonomy info file (tab-separated). Cols: species_name, strain_id, taxonomy, taxid_string, genome_len [, merged_alts]")
    parser.add_argument("-f", "--fold_filtered_markers", required=True, type=str,
                        help="Folder containing filtered marker seqs (*.ffn). Filenames include species_ and/or cluster_.")
    parser.add_argument("-d", "--database_fasta_path", required=True, type=str,
                        help="Absolute path to mpa database fasta file (will be appended unless --no-append).")
    parser.add_argument("-p", "--pkl_filename", required=True, type=str,
                        help="MetaPhlAn database pkl file path (existing or new). Output will be saved as *_NEW.pkl")
    parser.add_argument("--no-append", action="store_true", help="Do not append .ffn to database fasta (for testing merged_taxon only)")
    parser.add_argument("--no-build", action="store_true", help="Do not run bowtie2-build index (for testing only)")
    return parser.parse_args()

def normalize_str(s):
    if s is None:
        return s
    # strip and remove invisible whitespace
    return re.sub(r'\s+', '', s.strip())

def load_tax_info(path):
    """
    Load tax_infofile and create optional aliases when strain_id is species_X_cluster_Y.
    Returns (tax_dict, dup_info)
      tax_dict[strain_id] = {species_name, taxonomy, taxids, genome_len, merged_alts}
    merged_alts: list of (alt_taxonomy, alt_taxids, alt_count_int)
    """
    tax = {}
    dup_info = []
    with open(path, 'r', encoding='utf-8') as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                print(f"[WARN] Skipping invalid line {lineno} (not enough columns).")
                continue
            species_name = cols[0].strip()
            raw_id = cols[1].strip()
            strain_id = normalize_str(raw_id)
            taxonomy_str = cols[2].strip()
            taxid_str = cols[3].strip()
            genome_len = cols[4].strip()
            merged_alts_raw = cols[5].strip() if len(cols) > 5 else ""

            # parse merged_alts_raw into list of triples
            merged_alts = []
            if merged_alts_raw:
                for item in merged_alts_raw.split('|||'):
                    item = item.strip()
                    if not item:
                        continue
                    parts = item.split(':::')
                    if len(parts) != 3:
                        print(f"[WARN] malformed merged_alts at line {lineno}: '{item}'")
                        continue
                    alt_taxonomy = parts[0].strip()
                    alt_taxids = parts[1].strip()
                    alt_count_raw = parts[2].strip()
                    try:
                        alt_count_int = int(alt_count_raw)
                    except Exception:
                        alt_count_int = None
                    merged_alts.append((alt_taxonomy, alt_taxids, alt_count_int))

            if strain_id in tax:
                dup_info.append((strain_id, lineno))
                continue

            entry = {
                "species_name": species_name,
                "taxonomy": taxonomy_str,
                "taxids": taxid_str,
                "genome_len": genome_len,
                "merged_alts": merged_alts
            }
            tax[strain_id] = entry

            # If combined form species_X_cluster_Y, add aliases species_X and cluster_Y if absent
            m = re.match(r'^(species_\d+)_cluster_(\d+)$', strain_id)
            if m:
                species_alias = m.group(1)
                cluster_alias = f"cluster_{m.group(2)}"
                if species_alias not in tax:
                    tax[species_alias] = entry
                if cluster_alias not in tax:
                    tax[cluster_alias] = entry

    return tax, dup_info

def extract_candidate_token_from_filename(fname):
    """
    Prefer combined token species_X_cluster_Y; else species_X; else cluster_Y.
    Returns candidate token string or None.
    """
    m = re.search(r'(species_\d+_cluster_\d+)', fname)
    if m:
        return m.group(1)
    m = re.search(r'(species_\d+)', fname)
    if m:
        return m.group(1)
    m = re.search(r'(cluster_\d+)', fname)
    if m:
        return m.group(1)
    return None

def load_existing_db(pkl_path):
    """
    Load existing pkl bz2 database if present. Otherwise return empty structure.
    Ensure db contains 'merged_taxon'.
    """
    if not os.path.exists(pkl_path) or os.path.getsize(pkl_path) == 0:
        return {"markers": {}, "taxonomy": {}, "merged_taxon": {}}
    try:
        with bz2.open(pkl_path, "rb") as fh:
            db = pickle.load(fh)
            if not isinstance(db, dict) or 'markers' not in db or 'taxonomy' not in db:
                return {"markers": {}, "taxonomy": {}, "merged_taxon": {}}
            if 'merged_taxon' not in db:
                db['merged_taxon'] = {}
            return db
    except Exception as e:
        print(f"[WARN] Failed to load existing pkl ({pkl_path}): {e}. Starting from empty DB.")
        return {"markers": {}, "taxonomy": {}, "merged_taxon": {}}

def process_files(db, maker_dir, tax_dict, log):
    """
    Process files: species phase then cluster phase. Add taxonomy, markers, merged_taxon entries.
    """
    all_files = sorted([f for f in os.listdir(maker_dir) if f.endswith(".ffn")])
    species_files = [f for f in all_files if "species_" in f]
    cluster_files = [f for f in all_files if "cluster_" in f and f not in species_files]

    processed_files = []

    def add_taxonomy_if_absent(tmp_taxonomy_key, taxids, genome_len):
        if tmp_taxonomy_key in db["taxonomy"]:
            return False
        try:
            glen = int(genome_len) if genome_len and str(genome_len).isdigit() else None
        except Exception:
            glen = None
        db["taxonomy"][tmp_taxonomy_key] = (taxids, glen)
        return True

    def add_markers_from_file(filepath, strain_id, tmp_taxon_key):
        conflicts = 0
        added = 0
        for rec in SeqIO.parse(filepath, "fasta"):
            gid = rec.id.split()[0]
            if gid in db["markers"]:
                existing_taxon = db["markers"][gid].get("taxon")
                if existing_taxon == tmp_taxon_key:
                    continue
                else:
                    log.setdefault("marker_conflicts", []).append((gid, filepath, existing_taxon, tmp_taxon_key))
                    conflicts += 1
                    continue
            else:
                seq_len = len(str(rec.seq).strip())
                db["markers"][gid] = {
                    "clade": "t__" + strain_id,
                    "ext": [],
                    "len": int(seq_len),
                    "taxon": tmp_taxon_key
                }
                added += 1
        return added, conflicts

    def add_merged_taxon_entries(tmp_tax_key, taxinfo):
        """
        tmp_tax_key: taxonomy + "|t__" + sid (string)
        Add alternatives to db['merged_taxon'] keyed by (tmp_tax_key, taxids)
        """
        main_key = tmp_tax_key
        main_taxids = taxinfo["taxids"]
        key = (main_key, main_taxids)
        if key not in db['merged_taxon']:
            db['merged_taxon'][key] = []
        existing_alts = set((alt[0], alt[1]) for alt in db['merged_taxon'][key])
        for alt_tax, alt_taxids, alt_count in taxinfo.get("merged_alts", []):
            if (alt_tax, alt_taxids) in existing_alts:
                log.setdefault("merged_alts_duplicates", []).append((tmp_tax_key, alt_tax, alt_taxids))
                continue
            db['merged_taxon'][key].append((alt_tax, alt_taxids, alt_count))
            log.setdefault("merged_alts_added", []).append((tmp_tax_key, alt_tax, alt_taxids, alt_count))

    def process_list(file_list, phase_name):
        for fname in file_list:
            filepath = os.path.join(maker_dir, fname)
            sid_candidate = extract_candidate_token_from_filename(fname)
            if not sid_candidate:
                log.setdefault("no_id_files", []).append(fname)
                continue

            # Attempt exact match first; if fails and candidate is combined, try aliases
            sid = None
            tried = []
            if sid_candidate in tax_dict:
                sid = sid_candidate
            else:
                tried.append(sid_candidate)
                if '_cluster_' in sid_candidate:
                    sp = sid_candidate.split('_cluster_')[0]
                    cl = 'cluster_' + sid_candidate.split('_cluster_')[1]
                    # try species alias then cluster alias
                    for c in (sp, cl):
                        if c in tax_dict:
                            sid = c
                            tried.append(c)
                            break

            if sid is None:
                log.setdefault("id_not_in_taxinfo", []).append((fname, tried))
                continue

            taxinfo = tax_dict[sid]
            tmp_tax_key = taxinfo["taxonomy"] + "|t__" + sid

            # Add merged alternatives first (non-destructive)
            if taxinfo.get("merged_alts"):
                add_merged_taxon_entries(tmp_tax_key, taxinfo)

            # Add taxonomy if absent
            added_tax = add_taxonomy_if_absent(tmp_tax_key, taxinfo["taxids"], taxinfo["genome_len"])
            if not added_tax:
                log.setdefault("tax_overlaps", []).append((fname, sid, tmp_tax_key))

            # Add markers
            added_count, conflict_count = add_markers_from_file(filepath, sid, tmp_tax_key)
            processed_files.append((fname, sid, added_count, conflict_count))

    process_list(species_files, "species")
    process_list(cluster_files, "cluster")
    return processed_files

def append_fastas(maker_dir, db_fasta_path):
    pattern = os.path.join(maker_dir, "*.ffn")
    cmd = f"cat {pattern} >> {db_fasta_path}"
    return os.system(cmd)

def write_new_pkl(db, orig_pkl_path):
    out_pkl = orig_pkl_path.replace(".pkl", "_NEW.pkl")
    with bz2.BZ2File(out_pkl, "w") as of:
        pickle.dump(db, of, pickle.HIGHEST_PROTOCOL)
    return out_pkl

def write_logs(log, base_prefix):
    now = datetime.now().strftime("%Y%m%d_%H%M%S")
    logfile = f"{base_prefix}.log"
    with open(logfile, "w") as fh:
        fh.write(f"Log generated: {now}\n\n")
        fh.write("Files without ID (no species_/cluster_ token found):\n")
        for x in log.get("no_id_files", []):
            fh.write(f"  {x}\n")
        fh.write("\nFiles whose ID not found in tax_infofile (filename and attempted candidates):\n")
        for fname, tried in log.get("id_not_in_taxinfo", []):
            fh.write(f"  {fname}\ttried: {tried}\n")
        fh.write("\nTaxonomy overlaps (taxonomy key already existed):\n")
        for entry in log.get("tax_overlaps", []):
            fh.write(f"  {entry}\n")
        fh.write("\nMarker ID conflicts (same marker id existed with different taxon):\n")
        for gid, fname, old_tax, new_tax in log.get("marker_conflicts", []):
            fh.write(f"  marker_id: {gid}\tfile: {fname}\n    existing_taxon: {old_tax}\n    attempted_taxon: {new_tax}\n")
        fh.write("\nMerged_taxon entries added:\n")
        for entry in log.get("merged_alts_added", []):
            fh.write(f"  main_tmp_tax: {entry[0]}\n    alt_tax: {entry[1]}\n    alt_taxids: {entry[2]}\n    alt_count: {entry[3]}\n")
        fh.write("\nMerged_taxon duplicates skipped:\n")
        for entry in log.get("merged_alts_duplicates", []):
            fh.write(f"  main_tmp_tax: {entry[0]}\n    alt_tax: {entry[1]}\n    alt_taxids: {entry[2]}\n")
    return logfile

def main():
    args = get_args()
    tax_dict, dup_ids = load_tax_info(args.tax_infofile)
    if dup_ids:
        print("[WARN] Duplicate strain_id entries found in tax_infofile (showing first occurrences):")
        for sid, lineno in dup_ids[:20]:
            print(" ", sid, "at line", lineno)

    print("[DEBUG] tax_infofile sample keys (first 40):", list(tax_dict.keys())[:40])

    db = load_existing_db(args.pkl_filename)

    log = {
        "no_id_files": [],
        "id_not_in_taxinfo": [],
        "tax_overlaps": [],
        "marker_conflicts": [],
        # merged_alts entries added/duplicates logged dynamically
    }

    processed = process_files(db, args.fold_filtered_markers, tax_dict, log)

    if not args.no_append:
        rc = append_fastas(args.fold_filtered_markers, args.database_fasta_path)
        if rc != 0:
            print(f"[WARN] 'cat' command returned non-zero exit code: {rc}")
    else:
        print("[INFO] --no-append specified: skipping fasta append")

    if not args.no_build:
        bt2_cmd = f"bowtie2-build --threads 64 --large-index {args.database_fasta_path} {args.database_fasta_path.replace('.fasta','_NEW')}"
        os.system(bt2_cmd)
    else:
        print("[INFO] --no-build specified: skipping bowtie2-build")

    out_pkl = write_new_pkl(db, args.pkl_filename)
    base_prefix = os.path.splitext(os.path.basename(args.pkl_filename))[0] + "_addmarkers"
    logfile = write_logs(log, base_prefix)

    print("===== Summary =====")
    print(f"Original pkl: {args.pkl_filename}")
    print(f"New pkl: {out_pkl}")
    print(f"DB fasta appended to: {args.database_fasta_path}")
    print(f"Log file: {logfile}")
    print(f"Processed files: {len(processed)}")
    for fname, sid, added, conflicts in processed[:50]:
        print(f"  {fname}\t{sid}\tadded_markers:{added}\tmarker_conflicts:{conflicts}")
    if len(processed) > 50:
        print(f"  ... (and {len(processed)-50} more files)")

if __name__ == "__main__":
    main()
