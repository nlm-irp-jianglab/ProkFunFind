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

## Quick Start
```
prokfunfind.py -f config.yaml -gtab genome.table -o output.folder/prefix
```

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
