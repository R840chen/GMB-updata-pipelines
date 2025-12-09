#!/usr/bin/env python3
import argparse
import subprocess
import os
from pathlib import Path
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_args():
    parser = argparse.ArgumentParser(description="Batch run Prokka on genome files and collect FFN outputs.")
    parser.add_argument("--input_dir", "-i", required=True, help="Input directory of genome files.")
    parser.add_argument("--outdir", "-o", required=True, help="Output directory.")
    parser.add_argument("--extension", default="fna", help="File extension to search for (default: fna).")
    parser.add_argument("--cpus", type=int, default=8, help="CPUs per Prokka job.")
    parser.add_argument("--prokka_path", default="prokka", help="Path to Prokka executable.")
    parser.add_argument("--threads", type=int, default=4, help="Parallel Prokka jobs.")
    return parser.parse_args()

def run_prokka(prokka_path, infile, sample_dir, prefix, cpus):
    try:
        outdir = sample_dir / "prokka"
        if outdir.exists():
            shutil.rmtree(outdir)
        cmd = [prokka_path, str(infile), "--outdir", str(outdir), "--prefix", prefix, "--cpus", str(cpus)]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)
        ffn = outdir / f"{prefix}.ffn"
        if ffn.exists():
            return ffn
        print("Warning: FFN output missing for", prefix)
        return None
    except subprocess.CalledProcessError:
        print("Prokka failed for", infile)
        return None

def main():
    args = parse_args()
    inp = Path(args.input_dir)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    ffn_out = out / "ffn_outputs"
    ffn_out.mkdir(exist_ok=True)
    ann_dir = out / "prokka_annotations"
    ann_dir.mkdir(exist_ok=True)

    files = [p for p in inp.iterdir() if p.is_file() and p.suffix == f".{args.extension}"]
    print(len(files), "genome files detected.")

    futures = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        for f in files:
            prefix = f.stem
            sample_dir = ann_dir / prefix
            sample_dir.mkdir(exist_ok=True)
            futures.append(executor.submit(run_prokka, args.prokka_path, f, sample_dir, prefix, args.cpus))

        for fut in as_completed(futures):
            result = fut.result()
            if result:
                shutil.copy2(result, ffn_out / result.name)

if __name__ == "__main__":
    main()
