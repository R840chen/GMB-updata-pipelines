# GMB-update-pipeline
A pipelines used to update the human gut microbial gene biobank from genome to a new dataset that used to profiling the metagenomes information

This pipeline includes two parts: identify whether the provided species included in the GMB dataset or not; select maker genes from given genomes to update the GMB dataset.

Each part includes different steps. For the part one, it includes one script used to calculate the genome distance with the representative cluster genomes in GMB with dRep; as to the part two, it included one domain script which composed with 11 scripts, and the scripts takes different roles from gene prediction, uniref database annotation, core gene identification and maker gene selection to update the GMB dataset.

Installation

Dataset: 11167 representative cluster genomes dataset; uniref90 database; uniref50 database and GMB dataset.

11167 representative cluster genome dataset: xxx. This dataset were used to identify whether the provided species included in the GMB dataset or not, part 1 need. 

uniref90 database and uniref50 database: xxx. This database were used to annotate the predicted genes from the given genomes.

GMB dataset: xxx. This dataset were used to profile the species composition and relative abundance of human gut metagenomes, part 2 need.

Software required in Part 1: dRep(compare model) 
Suggestion download way:

$ pip install drep

$ conda install -c bioconda mash mummer fastani centrifuge -y

Software required in Part 2: prokka diamond mmseqs2 art bowtie2 samtools bamutil biopython
Suggestion download way:

$ conda install -c bioconda prokka diamond mmseqs2 art bowtie2 samtools bamutil biopython -y

The part 2 integrated script required the software in one conda environment, if you need to update the GMB-modified dataset, but the software not in one environment, the separated scripts were suggested to you and the usage of each script in xxx and you can take a usage of them step by step.


