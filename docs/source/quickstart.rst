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
