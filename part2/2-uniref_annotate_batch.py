#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch annotate FFN files with UniRef databases.

Behavior:
 - For each ffn file in input_dir create a per-sample output folder under outdir.
 - Use absolute paths for input and annotator so relative path resolution is stable.
 - Run annotator with cwd set to the per-sample output folder so outputs are written there.
 - After annotator finishes move translated and annotated files into sample_out.
 - Remove common hits and temporary files in likely locations.
 - Diamond options are passed as a single argument string to the annotator via --diamond-options.
"""
import argparse
import subprocess
from pathlib import Path
import shutil
import sys
import os
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Batch annotate FFN files with UniRef databases.")
    parser.add_argument("--input_dir", "-i", required=True, help="Directory containing FFN files")
    parser.add_argument("--outdir", "-o", required=True, help="Base output directory")
    parser.add_argument("--uniref90", required=True, help="UniRef90 Diamond DB")
    parser.add_argument("--uniref50", required=True, help="UniRef50 Diamond DB")
    parser.add_argument("--annotator", default="uniref_annotator.py", help="Annotator executable or script (path or name)")
    parser.add_argument("--threads", type=int, default=16, help="Threads to pass to diamond via annotator")
    return parser.parse_args()

def resolve_annotator(path_str):
    """
    Resolve annotator to an absolute path if possible.
    Resolution order:
      1) if absolute path and file exists -> use it
      2) if relative and exists relative to current working directory -> use that
      3) if exists relative to the script directory -> use that
      4) if found in PATH via shutil.which -> use that (returns absolute path)
      else raise FileNotFoundError
    """
    p = Path(path_str)
    # 1) absolute or explicit relative as given
    if p.is_absolute() and p.exists():
        return p.resolve()
    # 2) as given relative to cwd
    cand = Path.cwd() / path_str
    if cand.exists():
        return cand.resolve()
    # 3) relative to the directory where this script file lives
    try:
        script_dir = Path(__file__).resolve().parent
        cand2 = script_dir / path_str
        if cand2.exists():
            return cand2.resolve()
    except Exception:
        pass
    # 4) found in PATH
    which = shutil.which(path_str)
    if which:
        return Path(which).resolve()
    # not found
    raise FileNotFoundError(f"Annotator executable not found: {path_str}")

def safe_move(src_path, dest_dir):
    """
    Move src_path to dest_dir preserving basename.
    Overwrite existing file at destination.
    Return True on success, False otherwise.
    """
    try:
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / src_path.name
        if dest.exists():
            dest.unlink()
        shutil.move(str(src_path), str(dest))
        return True
    except Exception:
        return False

def find_and_move_outputs(input_path, sample_out):
    """
    Look for translated and annotated outputs in several likely locations and move them into sample_out.
    Locations searched: sample_out, current working directory, input parent directory.
    Name patterns checked include:
      - <input_path>.translated
      - <input_path>.annotated
      - <input_name>.translated
      - <input_name>.annotated
      - <input_stem>.translated
      - <input_stem>.annotated
    """
    input_path = Path(input_path)
    input_parent = input_path.parent
    name = input_path.name
    stem = input_path.stem

    suffixes = [".translated", ".annotated"]

    search_dirs = [Path.cwd(), input_parent, sample_out]

    for suf in suffixes:
        variants = [
            input_parent / (name + suf),
            input_parent / (stem + suf),
            Path.cwd() / (name + suf),
            Path.cwd() / (stem + suf),
            sample_out / (name + suf),
            sample_out / (stem + suf)
        ]
        for p in variants:
            if p.exists():
                try:
                    # move if not already in sample_out
                    if p.parent.resolve() != sample_out.resolve():
                        safe_move(p, sample_out)
                except Exception:
                    pass

def remove_hits_and_temps(input_path, sample_out):
    """
    Remove common hits and temporary files that annotator or diamond may create.
    Search locations: cwd, input parent, sample_out
    Use glob patterns to be robust.
    """
    input_path = Path(input_path)
    input_parent = input_path.parent
    stem = input_path.stem
    name = input_path.name

    search_dirs = [Path.cwd(), input_parent, sample_out]

    # common patterns to remove
    patterns = [
        f"{name}.hits", f"{stem}.hits",
        f"{name}.m8", f"{stem}.m8",
        f"{name}.daa", f"{stem}.daa",
        f"{name}.blast", f"{stem}.blast",
        f"{name}.blastout", f"{stem}.blastout",
        f"{stem}.diamond.*", f"{stem}.uniref90.*", f"{stem}.uniref50.*",
        f"{stem}*.hits", f"{stem}*.m8", f"{stem}*.daa", f"{stem}*.blast", f"{stem}*.diamond*",
        f"{name}.tmp", f"{stem}.tmp", f"{stem}*.tmp"
    ]

    for d in search_dirs:
        for pat in patterns:
            try:
                for p in d.glob(pat):
                    try:
                        if p.is_file():
                            p.unlink()
                        elif p.is_dir():
                            shutil.rmtree(p)
                    except Exception:
                        pass
            except Exception:
                pass

def main():
    args = parse_args()
    inp = Path(args.input_dir)
    if not inp.exists() or not inp.is_dir():
        print("Input directory not found:", args.input_dir, file=sys.stderr)
        sys.exit(2)

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # resolve annotator to an absolute path or fail early with helpful message
    try:
        annotator_path = resolve_annotator(args.annotator)
    except FileNotFoundError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        print("Please provide annotator as an absolute path or place it in the same directory as this script or in your PATH.", file=sys.stderr)
        sys.exit(2)

    # collect ffn files (case-insensitive suffix check) and ensure they are files
    ffn_files = [p for p in sorted(inp.iterdir()) if p.is_file() and p.suffix.lower() in [".ffn", ".fa", ".fasta"]]
    print(len(ffn_files), "FFN files found.")

    for f in ffn_files:
        sample_out = out / f.stem
        sample_out.mkdir(parents=True, exist_ok=True)

        # diamond options passed as a single argument string
        diamond_opts = f"--threads {args.threads}"

        # use absolute path for input to avoid cwd-relative resolution issues
        input_abs = f.resolve()

        # Build command to run annotator via the python interpreter to avoid execute-permission issues
        cmd = [
            sys.executable,
            str(annotator_path),
            str(input_abs),
            "--seqtype", "cds",
            "--uniref90db", args.uniref90,
            "--uniref50db", args.uniref50,
            "--diamond-options", diamond_opts
        ]

        # log the command and working directory
        print("Running:", " ".join([str(x) for x in cmd]), "in", str(sample_out))

        try:
            # Run annotator with working directory set to sample_out
            subprocess.run(cmd, check=True, cwd=str(sample_out))
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] annotator failed for {f.name}: returncode={e.returncode}", file=sys.stderr)
            # attempt to collect outputs and clean up even after failure
            find_and_move_outputs(input_abs, sample_out)
            remove_hits_and_temps(input_abs, sample_out)
            continue
        except FileNotFoundError:
            print(f"[ERROR] annotator executable not found after resolution: {annotator_path}", file=sys.stderr)
            sys.exit(2)

        # Move translated and annotated files into sample_out if found elsewhere
        find_and_move_outputs(input_abs, sample_out)

        # Remove hits and temporary files related to this input
        remove_hits_and_temps(input_abs, sample_out)

    print("All done.")

if __name__ == "__main__":
    main()
