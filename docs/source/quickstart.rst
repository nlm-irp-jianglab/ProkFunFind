.. ProkFunFind -

.. _quickstart:

************
Quick Start
************
In order to run ProkFunFind on your genome as soon as you have installed it, you can simply follow the next steps:


* Type:
  ``prokfunfind.py -h``
  to see all options available. All command-line options are described below:

.. literalinclude:: help.txt

* On a completely assembled genome (for example the representative EBI genome) run a command like the following:
  ``prokfunfind.py rep -b data.folder -f Function -g genome.inputs -o output.prefix``
This command will will detect the function ``Function`` in a genome or set of genomes designated
in a genome table (``genome.inputs``)  based on the functional configuration files placed
in the folder ``data.folder`` and outputs the result in ``output.prefix``
