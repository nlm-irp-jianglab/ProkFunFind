.. ProkFunFind

.. _inputs:


*************
Input dataset
*************


Command-line options
####################

  .. literalinclude:: help.txt


A typical ProkFunFind command looks like the following::

   prokfunfind -d queries -f fun -g ./genome-list.tsv -o ./out/search-out.

The options provide the following information:

====================  =================================================================================================================
Option                Description
====================  =================================================================================================================
-d, --databasedir     This option is used to provide the base directory where the query directories are stored.
--------------------  -----------------------------------------------------------------------------------------------------------------
-f, --function        This option is used to provide the name of the query directory.
--------------------  -----------------------------------------------------------------------------------------------------------------
-o, --outputprefix    This option is used to give the output file name prefix.
--------------------  -----------------------------------------------------------------------------------------------------------------
-g, --gtab            This option is used to provide the path to an input search genome table.
====================  =================================================================================================================



Input genomic data sets
########################

The genome input data should be organized so that all of the files associated
with one genome are in the same directory. Multiple genomes can be stored in the
same directory.

The directoy should include the genome sequences(``${prefix}.fna``),
protein sequences(``${prefix}.faa``), annotation(``${prefix}.gff``) and
associated annotation files.
The annotation files can include InterProScan tsv annotations (``${prefix}.tsv``),
KOfamScan tsv output (``${prefix}.kofam.tsv``), or EGGNog-mapper tsv output
(``${prefix}.emapper.annotations``). Note not all annotation files are required
to run ProkFunFind, for example if a search is being performed only using protein
sequences and profile HMMs, then no additional annotation files are required (see
below for more information).

The annotation files can be generated using the EGGNog-mapper, KOFamScan, and
InterProScan programs. The output form these programs does not require any
additional formatting to be used by ProkFunFind. Example commands to generate
the anntoation files are:

.. code-block::

  interproscan -t p -iprlookup --goterms --pathways -f tsv {fasta}.faa IPR

  exec_annotation --tmp-dir {path to tmp directory} -f detail-tsv  \
    -p {path to KofamScan profiles} -o {path to output file} {fasta}.faa

  emapper.py -i $infile -o {prefix} --temp_dir {path to tmp directory} \
    --data_dir {path to EGGNog-mapper data} --tax_scope {taxonomic scope} \
    --override --cpu {CPUS} --output_dir {path to output}


See the example below for the general format of the genome input files:

.. code-block::

  $ ls input.folder/
  MGYG-HGUT-03390.faa   # protein sequences
  MGYG-HGUT-03390.fna   # genome sequences
  MGYG-HGUT-03390.gff   # annotations
  MGYG-HGUT-03390_InterProScan.tsv   # Interproscan tsv output
  MGYG-HGUT-03390.kofam.tsv   # KOfamScan annotations
  MGYG-HGUT-03390.emapper.annotations   # EGGNog-mapper annotations

The genomes input is provided to the `ProkFunFind` program through a tab separated
two column table that includes the genome file prefixes and the paths to the
genome directories. This table should be provided through the `--gtab` argument:

.. code-block::

  MGYG-HGUT-03390	./input.directory/
  MGYG-HGUT-03391	./input.directory/

Any number of genomes inputs can be provided in this file and ProkFunFind will
run all searches on the provided genomes sequentially.


Input function configuration
############################

``-b`` should be followed by the data folder(``${data}``) that contains the configuration files for all functions.

.. code-block::

  $ ls data/
  Flagellar                   # Where the search files for a search
  Mucic_and_Saccharic_Acid    # Where the search files for Mucic_and_Saccharic_Acid stored


``-f`` should be followed by the name of the function directoy (``${function}``).

.. code-block::

  $ ls data/Mucic_and_Saccharic_Acid/
  config.ini
  search-term.tsv
  system.json
  query.fa
  query.fa.phr
  query.fa.pin
  query.fa.psq
  query.hmm


.. NOTE::

 Please remember to make fasta file blastable by running the command
 `makeblastdb -in {query.fasta} -dbtype prot`

 Similarly if profile HMMs are being used in the search they need to
 prepared for the search using the command `hmmpress {query.hmm}`



Configuration File
##################

config.ini
**********
The configuration files ``config.ini`` is where the settings for the ProkFunFind
search are specified. This file is made up of a main section and multiple other
sections related to specfic search approachces and filtering.

.. code-block::

    [main]
    cluster.tool   = DBSCAN
    system.file    = system.json
    faa_suffix     = .faa
    gff_suffix     = .gff3
    fna_suffix     = .fna
    search_terms = search_terms.tsv

    [DBSCAN]
    cluster.eps         = 4
    cluster.min_samples = 1

    [blast]
    blast.query    = bait.fa
    blast.exec     = blastp
    blast.evalue   = 1e-4
    blast.threads  = 1
    evalue = 1e-3
    ident_pct = 30
    bitscore = 30
    filter_file = hit_filter.tab

    [kofamscan]
    annot_suffix = .kofam.tsv



main
****
The main section of the configuration file contains general information about
the annotation file suffixes and points to the feature model file and search
terms table.

.. code-block::

  [main]
  cluster.tool   = DBSCAN
  system.file    = system.json
  faa_suffix     = .faa
  gff_suffix     = .gff3
  fna_suffix     = .fna
  search_terms = search_terms.tsv

===============  ==============================================================================
Name              Description
===============  ==============================================================================
search_terms      The name of the file that relates search term IDs and query IDs (see below)
---------------  ------------------------------------------------------------------------------
cluster.tool      The method used to cluster the genes
                  options:

                  * DBSCAN
---------------  ------------------------------------------------------------------------------
system.file       The name of the file that describe the structure of the function system
---------------  ------------------------------------------------------------------------------
faa_suffix        The suffix of the fasta file that contains the predicted amino acid
                  gene sequences
---------------  ------------------------------------------------------------------------------
fna_suffix        The suffix of the fasta file that contains the genome sequence(s)
---------------  ------------------------------------------------------------------------------
gff_suffix        The suffix of the file that contains the GFF gene annotations for the genome
===============  ==============================================================================


DBSCAN
******
If multiple hits are found in the genomes during the ProkFunFind searches, the
hits will be checked to see if they are in the same genomic region. This is done
using Density-Based Spatial Clustering of Applications with Noise (DBSCAN). For
more information on the scikit-learn DBSCAN implementation see `DBSCAN`_.

.. code-block::

     [DBSCAN]
     cluster.eps         = 4
     cluster.min_samples = 1

====================  =================================================================================================================
Name                  Description
====================  =================================================================================================================
cluster.eps           How close two genes should be in order for them to be considered to be in the same cluster. Distance is in
                      number of genes.
--------------------  -----------------------------------------------------------------------------------------------------------------
cluster.min_samples   Minimum number of genes of interest within range set by cluster.eps required for a given gene to be considered
                      a core member of a cluster.
====================  =================================================================================================================

.. _DBSCAN: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html


Search Approach Settings
************************
The remaining sections of the configuration file are used to defined search
approach specific settings. The settings allowed in each section are detailed
below.

'blast'
^^^^^^^
.. code-block::

   [blast]
   blast.query    = bait.fa
   blast.exec     = blastp
   blast.evalue   = 1e-4
   blast.threads  = 1
   evalue = 1e-3
   ident_pct = 30
   filter_file = hit_filter.tab


===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
blast.query       The name of the protein fasta file containing the query sequences. This fasta file needs to be indexed using the 'makeblastdb'
                  command.
---------------  --------------------------------------------------------------------------------------------------------------------------------
blast.exec        The executable tool will be passed to the cmd to run blast. Currently blastp is the only supported blast method.
---------------  --------------------------------------------------------------------------------------------------------------------------------
blast.evalue      The evalue will be passed to the cmd to run blast. Only hits below this will be returned from the blast program. Default is 10.
---------------  --------------------------------------------------------------------------------------------------------------------------------
blast.threads     The number of threads will be passed to the cmd to run blast. Default is 1.
---------------  --------------------------------------------------------------------------------------------------------------------------------
evalue            The evalue threshold used to filter the blast results after they are generated. This does not affect the raw BLAST output, but
                  is instead used to filter the results after they are generated. Default is 0.01
---------------  --------------------------------------------------------------------------------------------------------------------------------
ident_pct         The identity threshold used to filter blast hits. The default value is 30 (30% identity).
---------------  --------------------------------------------------------------------------------------------------------------------------------
filter_file       The file name of additional filtering settings for specific search terms (see filter file section below). Optional
===============  ================================================================================================================================

'hmmer'
^^^^^^^

.. code-block::

    ['hmmer']
    hmmer.query    = Hdc.hmm
    hmmer.exec     = hmmscan
    hmmer.evalue   = 1e-4
    hmmer.threads  = 1
    evalue = 1e-3
    bitscore = 0
    filter_file = hit_filter.tab

===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
hmmer.query       The name of the profile HMM file file.
---------------  --------------------------------------------------------------------------------------------------------------------------------
hmmer.exec        The executable tool will be passed to the cmd to run blast. Currently hmmscan is the only supported HMMER method.
---------------  --------------------------------------------------------------------------------------------------------------------------------
hmmer.evalue      The evalue will be passed to the cmd to run hmmscan. Only hits below this will be returned from the hmmscan program.
                  Default is 10.
---------------  --------------------------------------------------------------------------------------------------------------------------------
hmmer.threads     The number of threads will be passed to the cmd to run hmmscan. Default is the number of cpu cores detected on your machine.
---------------  --------------------------------------------------------------------------------------------------------------------------------
evalue            The evalue threshold used to filter the hmmscan results after they are generated. This does not affect the raw hmmscan
                  output, but is instead used to filter the results after they are generated. Default is 0.01
---------------  --------------------------------------------------------------------------------------------------------------------------------
bitscore         The bitscore threshold used to filter blast hits. The default value is 0.
---------------  --------------------------------------------------------------------------------------------------------------------------------
filter_file       The file name of additional filtering settings for specific search terms (see filter file section below). Optional
===============  ================================================================================================================================


'kofamscan'
^^^^^^^^^^^

.. code-block::

    [kofamscan]
    annot_suffix = .kofam.tsv
    evalue = 1e-3
    threshold = 1
    filter_file = hit_filter.tab

===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
annot_suffix      The file extension for the kofamscan prediction output.
---------------  --------------------------------------------------------------------------------------------------------------------------------
evalue            The evalue threshold used to filter the kofamscan results. Default is 0.01
---------------  --------------------------------------------------------------------------------------------------------------------------------
threshold         The threshold value is used to adjust the score thresholds which are used to determine if a kofamscan prediction is
                  significant or not. Kofamscan assigns a prediction score to each protein query for each KO number. If the score is above a
                  predetermined value for that KO, then the protein is putatively assigned to that KO. This score can be adjusted using this
                  threshold setting, which will be used to multiply the score needed to make it more or less strict.
                  Example:
                  .. code-block::

                    K00001  gene1  score: 10    KO_value: 12
                    - if the threshold is set to 1, then this gene would not be assigned to K00001
                    - if the threshold is set to 0.5, then the KO_value needed would be adjusted to 6 (12*0.5), resulting in the gene being
                      assigned to K00001
---------------  --------------------------------------------------------------------------------------------------------------------------------
filter_file       The file name of additional filtering settings for specific search terms (see filter file section below). Optional
===============  ================================================================================================================================

'interproscan'
^^^^^^^^^^^^^^

.. code-block::

  [interproscan]
  annot_suffix = _InterProScan.tsv

===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
annot_suffix      The name of the profile HMM file file.
===============  ================================================================================================================================


'emapper'
^^^^^^^^^^^

.. code-block::

    [emapper]
    annot_suffix = .emapper.annotations
    evalue = 1e-3
    filter_file = hit_filter.tab

===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
annot_suffix      The file extension for the EGGNog-mapper prediction output.
---------------  --------------------------------------------------------------------------------------------------------------------------------
evalue            The evalue threshold used to filter the EGGNog-mapper results. Default is 0.01
---------------  --------------------------------------------------------------------------------------------------------------------------------
filter_file       The file name of additional filtering settings for specific search terms (see filter file section below). Optional
===============  ================================================================================================================================


Filter file
###########
Separate search term specific filtering files can be provided as tab separated
tables that specify specific filtering parameters for any query. These
settings will be applied instead of the global filtering parameters that are set
in the configuration file. Any of the filtering values that are allowed in the
configuration file can be used in the filtering file. Filtering files can
be provided for each search approach being used through the filter_file
properties in the configuration sections.

An example of a filtering file can be seen here:

.. code-block::

    seq1  ident_pct  >=  50
    PF0001  evalue  <=  1e-100

The file consists of a four column, tab separated table. The first column
contains the IDs of the query (e.g., sequence ID, PFAM ID, Profile ID).
The second column contains the property that you want to filter by. The
fitlering properties allowed for each feature are listed in the configuration
file section of these docs. The fourth column contains the filtering logic
(>, <, >=, <=). The last column contains the value that will be used for the
filtering.

search terms
#############
The search terms file specifies the relationship between individual queries and
the broader search term IDs. This file is a three column table consiting of the
search terms IDs, query IDs, and search methods.

.. code-block::

    gene1  seq1   blast
    gene1  PFAM1  interproscan
    gene2  COG1   emapper

system
#######

Json formatted file that specify how the components are organized to perform a function.

  .. literalinclude:: example.json



======================  ========================================================
Name                    Description
======================  ========================================================
name/queryID:(*str*)    The name of the components/ The orthoID
----------------------  --------------------------------------------------------
components:(*list*)      The list of subcomponents
----------------------  --------------------------------------------------------
presence:(*option*)     "essential", "nonessential"
----------------------  --------------------------------------------------------
analogs:(*dict*)        Followed an equivalent component
======================  ========================================================
