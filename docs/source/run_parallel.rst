.. ProkFunFind - Detection of genes of functional interest in genomes

.. _run_parallel:

*****************
Batch processing
*****************

Once the configuration files in the data/function_name are made you can run ProkFunFind on every species in UHGG, a comprehensive database of gut microbiome genomes. 

At the NIH, UHGG is located here 
``/data/Irp-jiang/share/DB_Share/UHGG/``

At UMD, UHGG is located here 
``/fs/cbcb-lab/hall/data/UHGG/``




for example, if you are using the UHGG genomes, you will generate a file like below (let's name it ``cmd.list``):

::

  prokfunfind.py rep -b ProkFunFind/data -f Mucic_and_Saccharic_Acid -g /data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00001/genome/MGYG-HGUT-00001 -o output_for_MGYG-HGUT-00001
  prokfunfind.py rep -b ProkFunFind/data -f Mucic_and_Saccharic_Acid -g /data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00002/genome/MGYG-HGUT-00002 -o output_for_MGYG-HGUT-00002
  prokfunfind.py rep -b ProkFunFind/data -f Mucic_and_Saccharic_Acid -g /data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00003/genome/MGYG-HGUT-00003 -o output_for_MGYG-HGUT-00003
  prokfunfind.py rep -b ProkFunFind/data -f Mucic_and_Saccharic_Acid -g /data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00004/genome/MGYG-HGUT-00004 -o output_for_MGYG-HGUT-00004
  prokfunfind.py rep -b ProkFunFind/data -f Mucic_and_Saccharic_Acid -g /data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00005/genome/MGYG-HGUT-00005 -o output_for_MGYG-HGUT-00005
  prokfunfind.py rep -b ProkFunFind/data -f Mucic_and_Saccharic_Acid -g /data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00006/genome/MGYG-HGUT-00006 -o output_for_MGYG-HGUT-00006
  prokfunfind.py rep -b ProkFunFind/data -f Mucic_and_Saccharic_Acid -g /data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00007/genome/MGYG-HGUT-00007 -o output_for_MGYG-HGUT-00007
  prokfunfind.py rep -b ProkFunFind/data -f Mucic_and_Saccharic_Acid -g /data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00008/genome/MGYG-HGUT-00008 -o output_for_MGYG-HGUT-00008
  prokfunfind.py rep -b ProkFunFind/data -f Mucic_and_Saccharic_Acid -g /data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00009/genome/MGYG-HGUT-00009 -o output_for_MGYG-HGUT-00009
  prokfunfind.py rep -b ProkFunFind/data -f Mucic_and_Saccharic_Acid -g /data/Irp-jiang/share/DB_Share/UHGG/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00010/genome/MGYG-HGUT-00010 -o output_for_MGYG-HGUT-00010
  ...


Then depends on the computer or server you used, you can submit the job in batch. 

If you are on biowulf, you can submit the job with 
``swarm -f cmd.list -t 8``

If you are on single computer, you can run the job with parallel_ 
``cat -f cmd.list|parallel -j 8`` (if there are 64 cores).

.. _parallel: http://ftp.gnu.org/gnu/parallel/












