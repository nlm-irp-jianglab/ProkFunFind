.. GutFunFind

.. _inputs:


*************
Input dataset
*************

====================
Command-line options
====================

  .. literalinclude:: help.txt


Single genome mode(``rep``)
===========================

  ::

    $ run_GutFunFind.py rep -h
    > usage: run_GutFunFind.py rep [-h] -b  -f  -g  -o
    > 
    > optional arguments:
    >   -h, --help            show this help message and exit
    >   -b , --databasedir    The base dir of function
    >   -f , --function       Name of the function
    >   -g , --genomeprefix   The prefix of genome
    >   -o , --outputprefix   The output prefix


Pangenome mode(``pan``)
=======================

  ::

    $ run_GutFunFind.py pan  -h
    > usage: run_GutFunFind.py pan [-h] -b  -f  -p  -g  -d
    > 
    > optional arguments:
    >   -h, --help           show this help message and exit
    >   -b , --databasedir   The base dir of function
    >   -f , --function      Name of the function
    >   -p , --pangenome     The path to the pan-genome info dir
    >   -g , --genomeset     the dir of gzipped gff file
    >   -d , --outputdir     The output directory

=============================
Input genomic data sets
=============================

Single genome mode(``rep``)
===========================

  ``-g`` should be followed by the prefix of genome ``${prefix}``.
  
  The folder should include: genome sequences(``${prefix}.fna``), protein sequences(``${prefix}.faa``), annotation(``${prefix}.gff``) and Interproscan xml(``${prefix}.xml``: optional; or Interproscan tsv ``${prefix}.tsv``).
  
  In the example below:
  
  ::
  
    $ ls input.folder/
    > MGYG-HGUT-03390.faa   # protein sequences
    > MGYG-HGUT-03390.fna   # genome sequences
    > MGYG-HGUT-03390.gff   # annotations
    > MGYG-HGUT-03390.xml   # Interproscan xml output (only need for interproscan)
    > MGYG-HGUT-03390_InterProScan.tsv   # Interproscan tsv output (only need if Interproscan xml-formatted file is absent)
  
  The prefix should be ``input.folder/MGYG-HGUT-03390``
  
  If we need to analyze other genome(not included in the UHGG database),  we can use Prokka_ to annotate the genome.
  
  .. _Prokka: https://github.com/tseemann/prokka


Pangenome mode(``pan``)
=======================

  .. Attention::
  
     Not every representive species of UHGG has `pan-genome/` folder. Check before you run the command line. 


  ``-p`` should be followed by the path to the dir of pangenome information

  In the example below, ``-p`` should be followed by ``UHGG/uhgg_catalogue/MGYG-HGUT-014/MGYG-HGUT-01463/pan-genome/``.

  ::

    $ ls UHGG/uhgg_catalogue/MGYG-HGUT-014/MGYG-HGUT-01463/pan-genome/
    > genes_presence-absence_locus.csv  # csv file that specify the pangenome information
    > pan-genome.faa                    # protein sequences in the pangenome (required for blastp)
    > pan-genome_InterProScan.tsv       # Interproscan tsv file (required for interproscan)

  ``-g`` should be followed by the path to a dir of all gzipped gff files

  In the example below, ``-g`` should be followed by ``UHGG/all_genomes/MGYG-HGUT-014/MGYG-HGUT-01463``

  ::

    $ tree UHGG/all_genomes/MGYG-HGUT-014/MGYG-HGUT-01463
    > UHGG/all_genomes/MGYG-HGUT-014/MGYG-HGUT-01463
    > └── genomes1
    >     ├── GUT_GENOME096417.gff.gz
    >     ├── GUT_GENOME179203.gff.gz
    >     └── GUT_GENOME179220.gff.gz
    > 
    > 1 directory, 3 files

=============================
Input function configuration 
=============================

  ``-b`` should be followed by the data folder(``${data}``) that contains the configuration files for all functions.
  
  ::
  
    $ ls data/
    > Flagellar                   # Where the configuration file for Flagellar stored
    > Mucic_and_Saccharic_Acid    # Where the configuration file for Mucic_and_Saccharic_Acid stored
  
  
  ``-f`` should be followed by the name of the function(``${function}``). 
  
  ::
  
    $ ls data/Mucic_and_Saccharic_Acid/
    > bait.fa
    > bait.fa.phr
    > bait.fa.pin
    > bait.fa.psq
    > cluster.ini
    > config.ini
    > detect.ini
    > filter.ini
    > ortho_query_pair.tsv
    > system.json


.. Attention::

   Please remember to make bait.fa file blastable by running command line `makeblastdb -in bait.fa -dbtype prot`


=================================
Configuration file specification
=================================

config.ini
==========
  
  ::
  
    [main]
    detect.tool    = blast
    detect.config  = detect.ini
    cluster.tool   = DBSCAN
    cluster.config = cluster.ini
    system.file    = system.json
  
  
  
  ===============  ==============================================================================
  Name              Description
  ===============  ==============================================================================
  detect.tool       The method used to detect the genes
                    option:
                   
                    * blast
                    * hmmer
                    * interproscan
  ---------------  ------------------------------------------------------------------------------
  detect.config     The name of the configuration file that store the detect method information
  ---------------  ------------------------------------------------------------------------------
  cluster.tool      The method used to cluster the genes
                    option:
                   
                    * DBSCAN
  ---------------  ------------------------------------------------------------------------------
  system.file       The name of the file that describe the structure of the function system
  ===============  ==============================================================================


detect.ini
==========
  
Blast Configuration
--------------------

  ::
  
     [blast]
     blast.query    = bait.fa
     blast.exec     = blastp
     blast.evalue   = 1e-4
     blast.threads  = 8
     filter.config  = filter.ini
     map.ortho_pair = ortho_query_pair.tsv
  
  
  ===============  ================================================================================================================================
  Name              Description
  ===============  ================================================================================================================================
  ``[blast]``       The header of the detect configuration. Should be consistent with ``detect.tool`` in the ``config.ini`` file.
  ---------------  --------------------------------------------------------------------------------------------------------------------------------
  blast.exec        The executable tool will be passed to the cmd to run blast
  ---------------  --------------------------------------------------------------------------------------------------------------------------------
  blast.evalue      The evalue will be passed to the cmd to run blast
  ---------------  --------------------------------------------------------------------------------------------------------------------------------
  blast.threads     The number of threads will be passed to the cmd to run blast ([TODO]_: optional)
  ---------------  --------------------------------------------------------------------------------------------------------------------------------
  filter.config     The name of the configuration file that store the filter configuration 
  ---------------  --------------------------------------------------------------------------------------------------------------------------------
  map.ortho_pair    The name of the file that specify how the name(unique) of sequence in ``blast.query`` corrspond to  *orthoID*

                    An example of the map.ortho_query_pair files(separated by tab):
                   
                    ::
                   
                      $ cat ortho_query_pair.tsv
                      > gudD	ecoli_gudD
                      > gudP	ecoli_gudP
                      > garK	ecoli_garK
                      > garD	ecoli_garD
                      > garL	ecoli_garL
                      > garP	ecoli_garP
                      > garR	ecoli_garR
                      > gudD	cclostridioforme_GudD1
                      > garD	cclostridioforme_GarD
                      > gudA	cclostridioforme_gudA
                      > gudB	cclostridioforme_gudB
                      > gudC	cclostridioforme_gudC
                      > gudD	cclostridioforme_GudD2
                      > garL	cclostridioforme_GarL
                      > garR	cclostridioforme_GarR
  ===============  ================================================================================================================================
    

**filter.ini**
  
    ::
    
       [filter.global]
       evalue = 1e-6
       ident_pct = 30
  
       [filter.local]
       filter_file = hit_filter.tab
    
    ====================  =================================================================================================================
    Name                  Description
    ====================  =================================================================================================================
    ``[filter.global]``    Use to specify filter criteria that will apply to all hits
    --------------------  -----------------------------------------------------------------------------------------------------------------
     evalue                Use to specify filter evalue(maximal) criteria that will apply to all hits
    --------------------  -----------------------------------------------------------------------------------------------------------------
     ident_pct             Use to specify filter identity(minimal) criteria that will apply to all hits
    --------------------  -----------------------------------------------------------------------------------------------------------------
    ``[filter.local]``     Use to specify filter criteria for individual hit
    --------------------  -----------------------------------------------------------------------------------------------------------------
     filter_file           The relative path the the file containing filter information for individual hit
  
  
                           All the four columns:
  
                           1. hit_name(should be the same as access name of bait.fa) 
                           2. Attributes that can be used as criteria:
                              ``evalue/ident_pct/hit_start/hit_end/bitscore``
                           3. operator:">", "<", ">=", "<=", "==","!="
                           4. value that will beused as cutoff
  
                           An example of the filter_file file(separated by tab):
  
                           :: 
                            
                              $ cat hit_filter.tab
                              > cclostridioforme_GarR	evalue	<=	1e-110
                              > cclostridioforme_GarR	ident_pct	>=	50
  
    ====================  =================================================================================================================

.. Note::

   The parameters in ``detect.inc`` and ``filter.ini`` is detection method specific.

Interproscan Configuration
---------------------------

  ::
  
     [interproscan]
     orthoID_domain_precision = domain_precision.txt


  ==========================  =================================================================================================================
  Name                        Description
  ==========================  =================================================================================================================
  ``[Interproscan]``          The header of the detect configuration. Should be consistent with ``detect.tool`` in the ``config.ini`` file.
  --------------------------  -----------------------------------------------------------------------------------------------------------------
  orthoID_domain_precision    The name of the file that specify the precision of the domain corrspond to  *orthoID*

                              An example(separated by tab):

                              ::

                                $ cat domain_precision.txt
                                > K00575	G3DSA:1.10.155.10	0.908991
                                > K00575	PF01739	0.705724
                                > K00575	PF03705	0.704411
                                > K00575	PIRSF000410	0.99
                                > K00575	PR00996	0.708515
                                > K00575	PS50123	0.706645
                                > K00575	PTHR24422	0.634774
                                > ...
  ==========================  =================================================================================================================
    

cluster.ini
===========

  ::
  
     [DBSCAN]
     # Parameter pass to sklearn.cluster.DBSCAN
     cluster.eps         = 4
     # Parameter pass to sklearn.cluster.DBSCAN; The number of function-related-genes (or total weight) in a neighborhood for a point to be considered as a core point.
     cluster.min_samples = 1
  
  ====================  =================================================================================================================
  Name                  Description
  ====================  =================================================================================================================
  ``[DBSCAN]``          The header of the cluster configuration. Should be consistent with ``cluster.tool`` in the ``config.ini`` file.
  --------------------  -----------------------------------------------------------------------------------------------------------------
  cluster.eps           Parameters required for DBSCAN to run
  cluster.min_samples  
  ====================  =================================================================================================================
  
.. Note::

   The parameters in ``cluster.inc`` is cluster method specific. Currently DBSCAN is the only detection method supported.
  
system.json
===========
  
  Json formatted file that specify how the components are organized to perform a function.
  
  +-----------------------------------+------------------------------------+
  |  Example Structure                |     JSON formatted file            |
  +===================================+====================================+
  | .. image:: images/GutFunFind.jpg  |  .. literalinclude:: example.json  |
  |    :width: 550px                  |     :language: JSON                |
  |    :align: left                   |                                    |
  |    :alt: alternate text           |                                    |
  +-----------------------------------+------------------------------------+
  
  
  ======================  ========================================================
  Name                    Description
  ======================  ========================================================
  name/orthoID:(*str*)    The name of the components/ The orthoID 
  ----------------------  --------------------------------------------------------
  components:(*list*)      The list of subcomponents
  ----------------------  --------------------------------------------------------
  presence:(*option*)     "essential", "nonessential" or ([TODO]_) "forbidden"
  ----------------------  --------------------------------------------------------
  analogs:(*dict*)        Followed an equivalent component
  ======================  ========================================================


.. [TODO] To implementation later.
