# ProkFunFind: A pipeline to detect genes of functional interest

## Overview
  ProkFunFind is a pipeline that can be used to predict genes related to function of interest.
  Criteria for function detection include presence of gene, and genomic co-localization.

## Prerequisites
+ [Biopython](https://biopython.org/)
+ [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
+ [HMMER](http://eddylab.org/software/hmmer/hmmer.tar.gz)
+ [InterProScan](https://github.com/ebi-pf-team/interproscan)

## Installation

```
conda create -n ProkFun python=3
conda activate ProkFun
git clone git@github.com:nlm-irp-jianglab/ProkFunFind.git
pip install ./ProkFunFind
```

ProkFunFind is also available as a docker image at https://hub.docker.com/repository/docker/keithdt/prokfunfind/general

## Quick Start
```
prokfunfind.py -f config.yaml -gtab genome.table -o output.folder/prefix
```

## Example Datasets
We have included multiple search configuration files from our publications 
that have utilized ProkFunFind to search for functions in microbial genome 
datasets. These search configurations can be found in the data/ directory
of this ProkFunFind git repository. 

To run these searches on these example set of genomes you can run the
commands listed in the following table from within the respective functions
directory. 

|Function Directory | Description | Publication (PMID) | Search command |
|-------------------|-------------|--------------------|----------------|
|Bilirubin_reductase | Bilirubin reductase enzyme BilR | 36798240 | prokfunfind -f ./config.yaml -o ./test -g ../genome-table.tsv |
|Cysteine_degradation | Cysteine degredation pathway | 34745023 | prokfunfind -f ./config.yaml -o ./test -g ../genome-table.tsv |
|Histidine_decarboxylase | Histidine decarboxylase | 34563136 | prokfunfind -f ./config.yaml -o ./test -g ../genome-table.tsv |
|Equol | Equol biosynthesis operon | 35247986 | prokfunfind -f ./config.yaml -o ./test -g ../genome-table.tsv |
|Flagella | Flagellar genes | TBD | prokfunfind -f ./config-{bac or pseudo}.yaml -o ./test -g ../genome-table.tsv |


## ProkFunAnnotate
The ProkFunAnnotate pipelile has been developed to assist in the production
of the required annotation files for ProkFunFind searches. This annotation
pipeline can be found in the following git repository:
- https://github.com/nlm-irp-jianglab/ProkFunAnnotate

## Documentation and Tutorial
The documentation for ProkFunFind can be found here:
- https://prokfunfind.readthedocs.io/en/latest/

A tutorial describing how to run ProkFunFind searches can be found here:
- https://prokfunfind.readthedocs.io/en/docs-and-tests/tutorial/1-intro.html

## Scripts
A script, gen-function-itol.py, for generated iTol tree annotation files from 
the TSV output of a ProkFunFind search is provided in the scripts/ directory of 
the ProkFunFind repository: 
- https://github.com/nlm-irp-jianglab/ProkFunFind/tree/master/scripts
- This script takes the path to the ProkFunFind output directory, the character used to separate the gene from genome name, and a single column file with the feature IDs to plot as input. 

## Updates

* v0.1.1
    * add support for InterProScan
    * fix bug in cluster label

* v0.1.2
    * add support for InterProScan Tsv formatted output(will use tsv-formatted file in absence of xml-formatted file)

* v0.2.0
    * add support for pangenome

* v0.2.1
    * fix system.json process error and modify pangenome for speed improvement

* v0.2.2
    * add support for kofam
