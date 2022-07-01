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
conda create -n GutFun python=3
conda activate GutFun
git clone git@github.com:XiaofangJ/ProkFunFind.git
cd ProkFunFind && python setup.py install
```

## Quick Start
```
prokfunfind.py -b data.folder -f Mucic_and_Saccharic_Acid -g input.folder/MGYG-HGUT-02439 -o output.folder/prefix
```

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
