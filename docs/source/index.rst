.. GutFunFind documentation master file, created by
   sphinx-quickstart on Fri Jul 24 19:50:07 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GutFunFind's documentation!
======================================
  GutFunFind is a pipeline that can be used to predict genes related to function of interest.
  Criteria for function detection include presence of detected genes, and genomic co-localization.

  **In order to identify the correct genes, the user needs:**

   - Select the correct method and filter parameters to detect genes

     - Blast
     - Interproscan(Run in advance)
     - Hmmer
     - ...

   - Choose the correct bait sequences/domain signatures/HMM protein profiles for the components of interest
   - Specify the definition of gene cluster and cluster method 

     - DBSCAN
     - ...

   - Define the gene components for the function system in JSON format


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   inputs
   outputs
   run_parallel
   todo


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
