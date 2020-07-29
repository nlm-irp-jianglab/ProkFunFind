.. GutFunFind - Detection of genes of functional interest in genomes

.. _run_parallel:

*****************
Batch processing
*****************

Once the configuration files in the data/function_name are made you can run GutFunFind on every species in UHGG, a comprehensive database of gut microbiome genomes. 

At the NIH, UHGG is located here 
``Path to UHGG``

At UMD, UHGG is located here 
``Path to UHGG``




for example, if you are using the UHGG genomes, you will generate a file like below (let's name it ``cmd.list``):

::

  run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g path_to_MGYG-HGUT-00001 -o output_for_MGYG-HGUT-00001
  run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g path_to_MGYG-HGUT-00002 -o output_for_MGYG-HGUT-00002
  run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g path_to_MGYG-HGUT-00003 -o output_for_MGYG-HGUT-00003
  run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g path_to_MGYG-HGUT-00004 -o output_for_MGYG-HGUT-00004
  run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g path_to_MGYG-HGUT-00005 -o output_for_MGYG-HGUT-00005
  run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g path_to_MGYG-HGUT-00006 -o output_for_MGYG-HGUT-00006
  run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g path_to_MGYG-HGUT-00007 -o output_for_MGYG-HGUT-00007
  run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g path_to_MGYG-HGUT-00008 -o output_for_MGYG-HGUT-00008
  run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g path_to_MGYG-HGUT-00009 -o output_for_MGYG-HGUT-00009
  run_GutFunFind.py -b data.folder -f Mucic_and_Saccharic_Acid -g path_to_MGYG-HGUT-00010 -o output_for_MGYG-HGUT-00010
  ...


Then depends on the computer or server you used, you can submit the job in batch. 

If you are on biowulf, you can submit the job with 
``swarm -f cmd.list -t 8``

If you are on single computer, you can run the job with parallel_ 
``cat -f cmd.list|parallel -j 8`` (if there are 64 cores).

.. _parallel: http://ftp.gnu.org/gnu/parallel/
