<img src="shamrock-logo.png" align="right" width="200" style="margin-left: 15px;">

# SHAMROCK: Splitting Homeologue Ancestors by Mapping Repeats with Overlapping Common Kmers

kmer-based separation of allotetraploid parental subgenomes. This readme was written for `v0.3.0`.

## Introduction

SHAMROCK is a kmer-based tool for disentangling homeologous chromosomes in allotetraploid genomes by subtracting shared sequence content, 
identifying "allokmers" unique to each homoelogue, and clustering chromosomes into subgenomes based on shared allokmer content.


## Installation and setup

To run SHAMROCK, you will need to have [KMC](https://github.com/refresh-bio/KMC) and R installed.
R will need the `tidyverse` and `r-markdown` libraries. For default settings, [compleasm](https://github.com/huangnengCSU/compleasm)
should also be installed on the system.

These can be set up easily in a conda environment.

```
conda create -n shamrock
conda activate shamrock
mamba install -c conda-forge -c bioconda kmc r-base r-tidyverse r-rmarkdown compleasm
```

(If you don't use `mamba`, then replace the last command with `conda`.)

To run the examples, or to use [Telociraptor](https://github.com/slimsuite/telociraptor) for formatting input, you will also want to have
python3 on your system, and clone the [Telociraptor](https://github.com/slimsuite/telociraptor) repo next to the 
[SHAMROCK](https://github.com/slimsuite/shamrock) repo.

```
├── shamrock/
└── telociraptor/
```

For example, to set up and run in the directory `~/code`:

```
cd ~/code/
git clone https://github.com/slimsuite/shamrock
git clone https://github.com/slimsuite/telociraptor

# Test run
~/code/shamrock.sh
python3 ~/code/telociraptor/code/telociraptor.sh --help
```

## Running SHAMROCK

SHAMROCK is designed to run as a simple shell script on an input assembly in fasta format:

```
./shamrock/shamrock.sh <FASTA>
```

### Input format

The assembly should be in fasta format with a **single line per sequence** and the names of chromosomes starting `chrXX`,
where `XX` is a number from `01` to `N`. [Telociraptor](https://github.com/slimsuite/telociraptor) can be used to generate this format 
(see `get_examples.sh` and `chromformat.sh`).

```
>chr01 OX638062 ENA|OX638062|OX638062.1 Trifolium dubium genome assembly, chromosome: 1
CCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCT...
>chr02 OX638063 ENA|OX638063|OX638063.1 Trifolium dubium genome assembly, chromosome: 2
GTGAGTGTGTTTTATAGAATTAGTTTGGATAAACTTTTGAGCAAACACACCAAAAGAAAA...
>chr03 OX638064 ENA|OX638064|OX638064.1 Trifolium dubium genome assembly, chromosome: 3
TAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAA...
>chr04 OX638065 ENA|OX638065|OX638065.1 Trifolium dubium genome assembly, chromosome: 4
ACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCCTAAACCCTAAACCCTAAACC...
...
```

**NOTE:** This needs a bit of work, and is currently limited to assemblies with 10-99 chromosomes, as it assumes a two digit counter.
For assemblies with fewer than 10 chromosomes, rename sequences using `sed -i s/>chr/>chr0/ $SEQIN`. 
For 100 or more chromosomes, edit the `shamrock.sh` file to have `chr%03d` in place of `chr%02d`.

**NOTE:** Telociraptor identifies chromosomes based on a minimum length criterion - 10 Mbp by default.
This can be altered with the `minchrom=INT` setting.

### Primary outputs

The main outputs are a PDF of the clustering, along with two fasta files representing the parental partitioning.

### Known issues

The initial pairing of the chromosomes based on maximum shared kmers can go wrong, 
particularly when there has been a lot of rearrangements - or very long chromosomes - where
a same chromosome can be the "best" hit for many other chromosomes. Using HiC-guided manual 
assignment might be more effective. 

## Example

Partitioning of the [Shamrock genome](https://pmc.ncbi.nlm.nih.gov/articles/PMC11384199/), _Trifolium dubium_:

![Shamrock partitioning](example/drTriDubi3.shamrock.png)

Parent 1: chr01, chr02, chr03, chr04, chr07, chr08, chr10.

Parent 2: chr05, chr06, chr09, chr11, chr12, chr13, chr14, chr15.


## Citation

SHAMROCK has not yet been published. In the meantime, please cite this GitHub if you find it useful.

## Future Plans

SHAMROCK is a work in progress. Please post questions/requests as GitHub issues. One likely issue is where part of one chromosome has translocated to another 
chromosome since the ploidy event. A future

### A note on the name

SHAMROCK is a rather contrived acronym, as the inspirational species for its development was the 
[genome of the Irish Shamrock](https://pmc.ncbi.nlm.nih.gov/articles/PMC11384199/), _Trifolium dubium_.
Suggestions for a better expansion of the acronym are welcome!


