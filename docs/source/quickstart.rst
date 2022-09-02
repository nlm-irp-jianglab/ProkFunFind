.. ProkFunFind -

.. _quickstart:

************
Quick Start
************
In order to run ProkFunFind on your genome as soon as you have installed it, you can simply follow the next steps:


* Type:
  ``prokfunfind -h``
  to see all options available. All command-line options are described below:

.. literalinclude:: help.txt

* On a completely assembled genome (for example the representative EBI genome) run a command like the following:
  ``prokfunfind -f config.yaml -g genome.inputs -o output.prefix``

This command will will detect the function ``Function`` in a genome or set of genomes designated
in a genome table (``genome.inputs``)  based on the function defined in the
``config.yaml`` file and outputs the result in ``output.prefix``

To test your installation of ProkFunFind you can run the following command on
some sample data included with the repository:

.. code-block::

  prokfunfind -f {path to PFF repository}/ProkFunFind/data/HDC/config.yaml --gtab {path to PFF repository}/ProkFunFind/data/genome-table.tsv -o tmp

This command should print an output to the screen that looks like:

.. code-block::

  INFO:root:Checking configuration files
  INFO:root:Searching for function
  INFO:root:Identifying gene clusters
  INFO:root:Summarizing function presence and genes
  Failed to detect function: HDC in genome {path to ProkFunFind}/ProkFunFind/data/./test-genome//GTDB18040
  1 out of 2 essential components present
  0 out of 0 nonessential components present

This command will also produce multiple output files with the names tmp.GTDB18040*

See the other sections of the documents for descriptions of each of these files. 
