# GMB-updata-pipelines
A pipelines used to update the human gut microbial gene biobank from genome to a new dataset that used to profiling the metagenomes information

This pipeline includes two parts: identify whether the provided species included in the GMB dataset or not; select maker genes from given genomes to update the GMB dataset.

Each part includes different steps. For the part one, it includes one script used to calculate the genome distance with the representative cluster genomes in GMB with dRep; as to the part two, it included one domain script which composed with 11 scripts, and the scripts takes different roles from gene prediction, uniref database annotation, core gene identification and maker gene selection to update the GMB dataset.
