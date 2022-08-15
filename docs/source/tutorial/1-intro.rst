********************
ProkFunFind Tutorial
********************

This tutorial will walk you through how to use the ProkFunFind tool to search
for putative functions within microbial genomes. This tutorial will cover
the query format used by ProkFunFind, how to generate the annotation files,
and how to perform searches using ProkFunFind. The tutorial has been designed
so that each section can be run independently of the others. Each section will
also include pre-generated input and output files for you to look at while
going through the tutorial.

The tutorial materials are available on Github here: `PFF Tutorial <https://github.com/nlm-irp-jianglab/prokfunfind_tutorial>`_

These materials can be downloaded using the following command:

.. code-block::

   git clone https://github.com/nlm-irp-jianglab/ProkFunFind.git

All tutorial steps can be run from within the root directory of the repository.
From this root directory you should see the following files and directories:

================   ===================================================
File/Directory     Content
================   ===================================================
genomes/           Directory containing genomes to search during the
                   tutorial
----------------   ---------------------------------------------------
queries/           Directory containing query input files organized
                   in subdirectories
----------------   ---------------------------------------------------
precomputed-out/   Directory containing precomputed output from each
                   tutorial step
----------------   ---------------------------------------------------
out/               Directory where tutorial output will be generated
----------------   ---------------------------------------------------
ebi-list.tsv       Search genome table for EBI genome tutorial
----------------   ---------------------------------------------------
genome-list.tsv    Search genome table for main tutorial materials
================   ===================================================

A complete set of precomputed output from each step in the tutorial is
included in the `precomputed-out/` directory. This directory mirrors what
would be generated during the tutorial in the `out/` directory.

For instructions on how to install ProkFunFind please see: :doc:`../installation`

.. toctree::
  :maxdepth: 1
  :caption: Tutorial Sections

  3-queries
  4-seqsearch
  5-annotsearch
  6-mixedsearch
  7-pfa
  8-ebi
