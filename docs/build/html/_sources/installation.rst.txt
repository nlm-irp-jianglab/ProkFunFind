.. GutFunFind - Detection of genes of functional interest in genomes

.. _installation:

************
Installation
************

========================
GutFunFind dependencies
========================

  GutFunFind has several dependencies:
  
   - The Blast suite of program (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
   - The program *Hmmer* version 3.1 (http://hmmer.janelia.org/).
  
  All these programs should be installed (e.g., in the ``PATH``) in order to use GutFunFind.
  
  .. Attention::
  
     More program will be needed if more detect and cluster methods were added.
  
======================
Installation procedure
======================


Archive overview
=================
  
  Fold structure in GutFunFind package:
  
  * bin => contain the executable (run_GutFunFind.py)
  * data => the preconfigured data for each function
  * doc => the documentation in html and pdf
  * test => all needed for unit tests
  * GutFunFind => the GutFunFind python library
  * setup.py => the installation script


Install required dependence
============================

Install programs 
"""""""""""""""""
  Each different detection method rely on different programs. Make sure the program is installed.

Install python package
"""""""""""""""""""""""
  * bcbio-gff_
  * biopython_
  * scikit-learn_
  * configparser_
  * typing_
  * importlib_
  * argparse_
  * numpy_
  
  .. _bcbio-gff: https://github.com/chapmanb/bcbb/tree/master/gff
  .. _biopython: http://biopython.org/DIST/docs/tutorial/Tutorial.html
  .. _scikit-learn: https://scikit-learn.org/stable/
  .. _configparser: https://docs.python.org/3/library/configparser.html
  .. _typing: https://docs.python.org/3/library/typing.html
  .. _importlib: https://docs.python.org/3/library/importlib.html
  .. _argparse: https://docs.python.org/3/library/argparse.html
  .. _numpy: https://numpy.org/


Perform the installation
=========================

Download GutFunFind 
""""""""""""""""""""
  
  .. code-block:: bash
  
     git clone git@github.com:XiaofangJ/GutFunFind.git


Install GutFunFind
""""""""""""""""""""

  .. code-block:: bash

     conda create -n GutFun python=3
     conda activate GutFun
     cd GutFunFind 
     python setup.py install

Test GutFunFind
""""""""""""""""

  .. Attention::
  
     Still working on how to test the package. 
     It is **highly recommended** to run tests before performing the full installation.
     To test the libraries you just build: 
  
