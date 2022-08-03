.. ProkFunFind - Detection of genes of functional interest in genomes

.. _installation:

************
Installation
************

ProkFunFind dependencies
========================

ProkFunFind has several dependencies:

 - The Blast suite of program (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
 - The program *Hmmer* version 3.1 (http://hmmer.janelia.org/).

All these programs should be installed (e.g., in the ``$PATH``) in order to use ProkFunFind.

Additional python packages will be installed during the ProkFunFind installation
including:

  - biopython_

  - scikit-learn_

  - configparser_

  - typing_

  - importlib_

  - argparse_

  - numpy_

  - six_

  .. _biopython: http://biopython.org/DIST/docs/tutorial/Tutorial.html
  .. _scikit-learn: https://scikit-learn.org/stable/
  .. _configparser: https://docs.python.org/3/library/configparser.html
  .. _typing: https://docs.python.org/3/library/typing.html
  .. _importlib: https://docs.python.org/3/library/importlib.html
  .. _argparse: https://docs.python.org/3/library/argparse.html
  .. _numpy: https://numpy.org/
  .. _six: https://github.com/benjaminp/six


Archive overview
=================

  * Folder structure in ProkFunFind package:

  * bin => contain the executable (prokfunfind.py)

  * data => the preconfigured data for each function

  * doc => the documentation in html and pdf

  * test => all needed for unit tests

  * ProkFunFind => the ProkFunFind python library

  * setup.py => the installation script



Installation Instructions
=========================

Download ProkFunFind
""""""""""""""""""""

  .. code-block:: bash

     git clone git@github.com:nlm-irp-jianglab/ProkFunFind.git


Install ProkFunFind
""""""""""""""""""""
We recommend installing ProkFunFind within a python virtual environment.

  .. code-block:: bash

     conda create -n ProkFun python=3
     conda activate ProkFun
     cd ProkFunFind
     python setup.py install
