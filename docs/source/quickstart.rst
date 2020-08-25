.. GutFunFind - 

.. _quickstart:

.. Attention::

   We add pan-genome support and the command line is different from previous version.  

     1. To run the original command, please add "`rep`" after `run_GutFunFind.py`
     2. Use relative path instead of absolute path in all the configuration files (`base.dir` is removed)

************
Quick Start
************
  In order to run GutFunFind on your genome as soon as you have installed it, you can simply follow the next steps:

  
  * Type:
    ``run_GutFunFind.py -h``
    to see all options available. All command-line options are described below:
  
  .. literalinclude:: help.txt
  
  * On a completely assembled genome (for example the representative UHGG genome)
  
    ``run_GutFunFind.py rep -b data.folder -f Mucic_and_Saccharic_Acid -g input.folder/MGYG-HGUT-02439 -o output.prefix`` 
  
    will detect the function ``Mucic_and_Saccharic_Acid`` in a complete genome(``input.folder/MGYG-HGUT-02439``)  based on functional configuration files placed in the folder ``data.folder`` and outputs the result in ``output.prefix``

  * On preconfigure pangenome set 
  
    ``run_GutFunFind.py pan  -b GutFunFind/data  -f Flagellar  -g UHGG/all_genomes/MGYG-HGUT-014/MGYG-HGUT-01463 -p UHGG/uhgg_catalogue/MGYG-HGUT-014/MGYG-HGUT-01463/pan-genome -d output.dir``
  
    will detect the function ``Flagellar`` in all the genomes in folder (``UHGG/all_genomes/MGYG-HGUT-014/MGYG-HGUT-01463``) based on functional configuration files placed in the folder ``data.folder`` and outputs the result in ``output.prefix``
    These genomes are belong to the same species and the pangenome info is stored in ``UHGG/uhgg_catalogue/MGYG-HGUT-014/MGYG-HGUT-01463/pan-genome``

