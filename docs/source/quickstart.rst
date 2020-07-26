.. GutFunFind - 

.. _quickstart:


************
Quick Start
************

In order to run GutFunFind on your favorite dataset as soon as you have installed it, you can simply follow the next steps:


* Type:
  ``run_GutFunFind.py -h``
  to see all options available. All command-line options are described below:

.. literalinclude:: help.txt

* On a completely assembled genome 

  ``run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g input.folder/MGYG-HGUT-02439 -o output.prefix`` 

  will detect the function ``Mucic_and_Saccharic_Acid`` in a complete genome(``input.folder/MGYG-HGUT-02439``)  based on functional configuration files placed in the folder ``data.folder`` and outputs the result in ``output.prefix``

