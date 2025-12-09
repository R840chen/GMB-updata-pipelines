# GMB-update-pipeline

**A pipeline to update the Human Gut Microbial Gene Biobank (GMB)**

---

## Overview

`GMB-update-pipeline` is a collection of scripts and workflows designed to update the Human Gut Microbial Gene Biobank (GMB). The pipeline has two main parts:

1. **Presence detection (Part I)** — determine whether a provided genome/species is already represented in the GMB dataset, using genome distance comparisons (dRep).
2. **Marker-gene selection & update (Part II)** — predict genes from input genomes, annotate them (UniRef90/50), identify core genes and select marker genes to update the GMB dataset.

Each part contains one or more scripts. Part I contains a script that computes genome distances against the GMB representative cluster genomes using `dRep`. Part II is a higher-level workflow assembled from multiple scripts (11 small components) handling gene prediction, annotation, core gene identification and marker-gene selection.

---

## Table of Contents

* [Overview](#overview)
* [Repository Structure](#repository-structure)
* [Datasets Required](#datasets-required)
* [Software Dependencies](#software-dependencies)
* [Installation](#installation)
* [Usage](#usage)

  * [Part I — Presence detection](#part-i--presence-detection)
  * [Part II — Marker-gene selection & update](#part-ii--marker-gene-selection--update)
* [Notes and Recommendations](#notes-and-recommendations)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)

---

## Repository Structure

```
GMB-update-pipeline/
├── part1_presence_detection/        # scripts and helpers for Part I (dRep-based)
├── part2_marker_selection/          # main integrated workflow + component scripts (11 scripts)
├── docs/                            # documentation, usage examples
├── envs/                            # optional conda environment YAML files
├── examples/                        # example inputs and expected outputs
└── README.md                        # this file
```

---

## Datasets Required

You will need the following datasets (place paths or download instructions in `docs/` or an `examples/` folder):

* **GMB representative clusters (11,167 genomes)** — used by Part I to decide whether a query genome/species is already included. (Provide download path or link in `docs/`.)
* **UniRef90 and UniRef50 databases** — used to annotate predicted proteins in Part II.
* **GMB dataset (marker profiles / reference)** — used to profile species composition and relative abundance of metagenomes and to update the GMB marker set.

> **Note:** Replace the placeholder `xxx` with the real download URLs or local paths in `docs/DATASETS.md`.

---

## Software Dependencies

### Part I (Presence detection)

* `dRep` (compare mode)
* Optional / recommended: `mash`, `mummer` (nucmer), `fastANI`, `centrifuge`

Suggested install commands:

```bash
pip install drep
conda install -c bioconda mash mummer fastani centrifuge -y
```

### Part II (Marker-gene selection & update)

* `prokka`
* `diamond`
* `mmseqs2`
* `art_illumina` (or `art`)
* `bowtie2`
* `samtools`
* `bamutil`
* `biopython`

Suggested install command:

```bash
conda create -n gmb-update -c bioconda prokka diamond mmseqs2 art bowtie2 samtools bamutil biopython -y
conda activate gmb-update
```

> **Important:** The integrated script for Part II assumes all required tools are available within a single conda environment. If you cannot install all tools in one environment, run the component scripts separately in suitable environments as documented in `docs/COMPONENT_USAGE.md`.

---

## Installation

1. Clone this repository:

```bash
git clone https://github.com/<your-org>/GMB-update-pipeline.git
cd GMB-update-pipeline
```

2. Prepare datasets (see `docs/DATASETS.md`) and ensure required software is installed and on your `PATH`.

3. (Optional) Create the recommended conda environment for Part II:

```bash
conda env create -f envs/gmb-update.yml   # if provided
conda activate gmb-update
```

---

## Usage

### Part I — Presence detection

This part checks whether provided genomes or species are already represented in GMB.

**Example:**

```bash
# run the presence-detection script (example)
python part1_presence_detection/check_presence.py \
  --input-genomes /path/to/query_genomes/ \
  --gmb-rep-genomes /path/to/11167_representatives.txt \
  --out results/presence_report.tsv
```

**What it does:**

* Computes genome distances against the GMB representative genomes using `dRep`.
* Outputs a tabular report listing whether each query genome is represented, with distance/ANI values.

### Part II — Marker-gene selection & update

This part performs gene prediction, annotation, core-gene detection and marker selection to update the GMB marker set.

**High-level example (integrated workflow):**

```bash
# run the integrated Part II workflow (example)
python part2_marker_selection/run_pipeline.py \
  --input-genomes /path/to/genomes/ \
  --uniref90 /path/to/uniref90.fasta \
  --uniref50 /path/to/uniref50.fasta \
  --gmb-reference /path/to/gmb_reference/ \
  --outdir results/part2_output
```

**If tools are installed in different environments:**

* Use the separated component scripts in `part2_marker_selection/components/` and follow `docs/COMPONENT_USAGE.md` to run each step sequentially.

---

## Notes and Recommendations

* Provide absolute paths for input/output arguments to avoid ambiguity.
* Use the recommended conda environment to ensure reproducible behavior.
* Large databases (UniRef, representative genomes) may be large; make sure you have sufficient disk space.
* Where possible, run compute-intensive steps (e.g., annotation, all-vs-all comparisons) on a cluster or machine with multiple cores.

---

## Contributing

Contributions are welcome. Please open an issue or submit a pull request. If you are adding features or changing behavior, update documentation in `docs/` and include tests/examples in `examples/`.

---

## License

Specify your license here (e.g., MIT, GPL-3.0). Add a `LICENSE` file at the repository root.

---

## Contact

If you have questions, open an issue or contact the maintainer: `chenc` (or add an email/GitHub handle).

---

*Last updated: replace with date when you publish this README.*
