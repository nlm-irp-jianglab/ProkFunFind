# GutFunFind: A pipeline to detect genes of functional interest

## Overview
  GutFunFind is a pipeline that can be used to predict genes related to function of interest.
  Criteria for function detection include presence of gene, and genomic co-localization.

## Prerequisites
+ [Biopython](https://biopython.org/)
+ [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
+ [HMMER](http://eddylab.org/software/hmmer/hmmer.tar.gz)
+ [InterProScan](https://github.com/ebi-pf-team/interproscan)

## Installation

```
conda create -n GutFun python=3
conda activate GutFun
git clone git@github.com:XiaofangJ/GutFunFind.git
cd GutFunFind && python setup.py install
```

## Quick Start
```
run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g input.folder/MGYG-HGUT-02439 -o output.folder/prefix
```
