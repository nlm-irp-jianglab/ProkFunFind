.. ProkFunFind

.. _inputs:


*************
Input dataset
*************


Command-line options
####################

  .. literalinclude:: help.txt


A typical ProkFunFind command looks like the following::

   prokfunfind -f queries/fun/config.yaml -g ./genome-list.tsv -o ./out/search-out.

The options provide the following information:

====================  =================================================================================================================
Option                Description
====================  =================================================================================================================
-f, --function        This option is used to provide the path to the configuration yaml file.
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
genome directories. The path can be given as an absolute path, the full path to 
the genome file directory, or as a relative path from the gtab table file, specified
as a relative path starting with either './' or '../'.

This table should be provided through the `--gtab` argument:

.. code-block::

  MGYG-HGUT-03390	./input.directory/
  MGYG-HGUT-03391	./input.directory/

Any number of genomes inputs can be provided in this file and ProkFunFind will
run all searches on the provided genomes sequentially.

These input data can be generated using the `ProkFunAnnotate <https://github.com/nlm-irp-jianglab/ProkFunAnnotate>`_
snakemake pipeline. This pipeline can be run to generate the EggNOG-mapper and KOFamScan
data used by snakemake from a set of input genomes. The input file needed to run the
snakemake pipeline is the same one used by the ProkFunFind --gtab argument. This annotation pipeline can be used to
generate the annotation files for genomes starging from just a genome sequence fasta file, or a genome sequence
fasta file with a GenBank file, like can be downloaded from NCBI's Genbank database. 


Input function configuration
############################
``-f`` should be followed by the path to the configuration yaml file (``${function}``).

.. code-block::

  $ ls data/Mucic_and_Saccharic_Acid/
  config.yaml
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
The configuration file (``config.yaml``) is the main input for each query. This file is used to
provide general settings like file extensions and filtering thresholds, and to
provide the definition of the function being searched for. This file is split
into two sections separated by '---' with the first section containing
the

configuration section
**********************
The configuration section of the ``config.yaml`` file is where the settings for the ProkFunFind
search are specified. This file is made up of a main section and multiple other
sections related to the specific search approaches and filtering.

.. code-block::

    ---
    main:
      cluster_tool: DBSCAN
      faa_suffix: .faa
      gff_suffix: .gff
      fna_suffix: .fna
    DBSCAN:
      cluster_eps: 4
      cluster_min_samples: 2
    hmmer:
      hmmer_query: ./query.hmm
      hmmer_exec: hmmscan
      hmmer_threads: 1
      evalue: 1e-3
    blast:
      blast_query: ./query.fa
      blast_exec: blastp
      blast_threads: 1
      evalue: 1e-3
    kofamscan:
      annot_suffix: .kofam.tsv
      threshold: 0.5
    emapper:
      annot_suffix: .emapper.annotations
    interproscan:
      annot_suffix: _InterProScan.tsv




main
****
The main section of the configuration file contains general information about
the annotation file suffixes and points to the feature model file and search
terms table.

.. code-block::

  main:
    cluster_tool: DBSCAN
    faa_suffix: .faa
    gff_suffix: .gff
    fna_suffix: .fna

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

  DBSCAN:
    cluster_eps: 4
    cluster_min_samples: 2

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


    blast:
      blast_query: ./bait.fa
      blast_exec: blastp
      blast_evalue: 1e-4
      blast_threads: 1
      evalue: 1e-3
      ident_pct: 30



===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
blast_query       The name of the protein fasta file containing the query sequences. This fasta file needs to be indexed using the 'makeblastdb'
                  command. This can be given as an absolute path to the query sequence file, or as a relative path 
                  starting with './' or '../'
---------------  --------------------------------------------------------------------------------------------------------------------------------
blast_exec        The executable tool will be passed to the cmd to run blast. Currently blastp is the only supported blast method.
---------------  --------------------------------------------------------------------------------------------------------------------------------
blast_evalue      The evalue will be passed to the cmd to run blast. Only hits below this will be returned from the blast program. Default is 10.
---------------  --------------------------------------------------------------------------------------------------------------------------------
blast_threads     The number of threads will be passed to the cmd to run blast. Default is 1.
---------------  --------------------------------------------------------------------------------------------------------------------------------
evalue            The evalue threshold used to filter the blast results after they are generated. This does not affect the raw BLAST output, but
                  is instead used to filter the results after they are generated. Default is 0.01
---------------  --------------------------------------------------------------------------------------------------------------------------------
ident_pct         The identity threshold used to filter blast hits. The default value is 30 (30% identity).
===============  ================================================================================================================================

'hmmer'
^^^^^^^

.. code-block::

    hmmer:
      hmmer_query: ./Hdc.hmm
      hmmer_exec: hmmscan
      hmmer_evalue: 1e-4
      hmmer_threads: 1
      evalue: 1e-3
      bitscore: 0

===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
hmmer.query       The name of the profile HMM file. The HMM file should be indexed with hmmpress. The path can be given as an absolute path
                  or a relative path starting with './' or '../'. 
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
===============  ================================================================================================================================


'kofamscan'
^^^^^^^^^^^

.. code-block::

    kofamscan:
      annot_suffix: .kofam.tsv
      evalue: 1e-3
      threshold: 1

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
===============  ================================================================================================================================

KofamScan Annotation File
"""""""""""""""""""""""""""
The KofamScan annotation file is a tabular file generated using teh kofam_scan tool provided here (https://github.com/takaram/kofam_scan).
The tabular oputput generated by KofamScan consists of seven columns where each line represents a putative annotation for a gene
with certain lines being marked with a '*' if they pass predetermined score threshold for that KO. An example KofamScan annotation file can be seen 
below. 
In this example multiple possible annotations are reproted for GCA_001563995.1_00002, with only one, K07333, being considered signfiicant based on 
the score being higher than the threshold for that KO. Adjusting the threshold parameter in the ProkFunFind configuration file will multiple the 
value for a signifcant hit by the value given in the configuration file. For example a threshold setting of 0.5, would make it so that the required
threshold for K12511 in the output below change to 64.735 (i.e., 129.47 * 0.5). This can potentially affect what annotations are considered signficant
in the annotation file. 

.. code-block:: 

    #       gene name               KO      thrshld score   E-value         "KO definition"
    #       ---------               ------  ------- ------  ---------       -------------
    *       GCA_001563995.1_00002   K07333  99.57   216.9   2.1e-65         "archaeal flagellar protein FlaJ"
            GCA_001563995.1_00002   K12511  129.47  79.5    1e-23           "tight adherence protein C"
            GCA_001563995.1_00002   K12510  83.50   50.3    1e-14           "tight adherence protein B"
    ...


'interproscan'
^^^^^^^^^^^^^^

.. code-block::

  interproscan:
    annot_suffix: _InterProScan.tsv
    evalue: 1e-3

===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
annot_suffix      the file extension for the InterProScan annotation file. 
---------------  --------------------------------------------------------------------------------------------------------------------------------
evalue            The evalue threshold used to filter the InterProScan results. Default is 0.01
===============  ================================================================================================================================

InterProScan Annotation File
"""""""""""""""""""""""""""""
The InterProScan annotation file is a tab separated table generated by the InterProScan command line tool. Information on installing and running
InterProScan is available here: https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html
The annotation file generated by InterProScan consists of 13 columns, with the key ones for ProkFunFind being the gene ID (1), type of annotation (4),
annotation (5), and evalue (9). An example of an InterProScan annotaiton file can be seen below: 

.. code-block:: 

  GCA_001563995.1_00003   1fde1e3b11f2b914cce637e9d04ad07a        140     Pfam    PF07790 Archaeal Type IV pilin, N-terminal      4       48      1.6E-5  T       02-05-2023      IPR012859       Archaeal Type IV pilin, N-terminal
  GCA_001563995.1_00003   1fde1e3b11f2b914cce637e9d04ad07a        140     TIGRFAM TIGR02537       arch_flag_Nterm: archaeal flagellin N-terminal-like domain      3       25      2.3E-6  T       02-05-2023      IPR013373       Flagellin
  ...

'prokka'
^^^^^^^^^^^

.. code-block::

  prokka:
    annot_suffix: .prokka.tsv

===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
annot_suffix      The file extension for the Prokka tabular annotation file output.
===============  ================================================================================================================================

Prokka Annotation File
"""""""""""""""""""""""""""""
Prokka can provide preliminary annotations including COG assignments for a genome (see the Prokka tool and documentation here: https://github.com/tseemann/prokka).
The Prokka annotation file is a genreated as part of the Prokka pipeline with a '.tsv' extension. If you are using multiple annotaiton files
with ProkFunFind, be careful about possible overlapping generic file extensions like '.tsv'. It would be recommended to rename the file to something
more descriptive like '.prokka.tsv'. The Prokka annotaiton file consists of 7 columns including the gene name (1), and COG (6). No e-values or bit scores 
are provided for these annotaitons, so no filtering of these results is possible in ProkFunFind. An example Prokka annotaiton file can be seen below: 

.. code-block::

  IMMGDCFF_00012  CDS     876     yisK            COG0179 putative protein YisK
  IMMGDCFF_00013  CDS     1308    asnS_1  6.1.1.22        COG0017 Asparagine--tRNA ligase

'bakta'
^^^^^^^^^^^

.. code-block::

  bakta:
    annot_suffix: .bakta.tsv

===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
annot_suffix      The file extension for the Bakta tabular annotation file output. 
===============  ================================================================================================================================

Bakta Annotation File
"""""""""""""""""""""""""""""
Bakta can provide preliminary anntoation files including KEGG KOs and COGs. The Bakta tool is available here: https://github.com/oschwengers/bakta
The '.tsv' output file from bakta provides information on the called genes and their putative annotations. The gene IDs are provided in column 5 while
the annotations can be found in a comma separated list in column 7. An example bakta annotaiton file can be seen below: 

.. code-block:: 

  contig_11       cds     18802   22275   -       KBDKFM_00100    bcsC    cellulose synthase complex outer membrane protein BcsC  COG:COG5010, COG:UW, GO:0019867, GO:0030244, KEGG:K20543, RefSeq:WP_167582918.1, SO:0001217, UniParc:UPI001430EF41, UniRef:UniRef100_A0A7D7DNT3, UniRef:UniRef50_P37650, UniRef:UniRef90_P37650
  contig_11       cds     22257   23363   -       KBDKFM_00105    bcsZ    cellulose synthase complex periplasmic endoglucanase BcsZ       COG:COG3405, COG:G, EC:3.2.1.4, GO:0005576, GO:0008810, GO:0030245, KEGG:K20542, RefSeq:WP_167582921.1, SO:0001217, UniParc:UPI00142F9A65, UniRef:UniRef100_A0A7H9QRD3, UniRef:UniRef50_P37651, UniRef:UniRef90_P37651

'emapper'
^^^^^^^^^^^

.. code-block::

  emapper:
    annot_suffix: .emapper.annotations
    evalue: 1e-3

===============  ================================================================================================================================
Name              Description
===============  ================================================================================================================================
annot_suffix      The file extension for the EGGNog-mapper prediction output.
---------------  --------------------------------------------------------------------------------------------------------------------------------
evalue            The evalue threshold used to filter the EGGNog-mapper results. Default is 0.01
===============  ================================================================================================================================

EggNOG-mapper Annotation File
"""""""""""""""""""""""""""""""
The EggNOG-mapper anontation tool produces a tabular output that provides annotaiton information including COG assignments for a given genome. Details
on how to install and run EggNOG-mapper can be found here: https://github.com/eggnogdb/eggnog-mapper
The EggNog mapper annotation file includes the gene IDs in column 1, the e-values in column 3, and the COG annotaitons in column 5. The COG annotations
include the base level COGs as well as COGs defined at different taxonomic levels. ProkFunFind searches can be performed with any of these COGs at any
taxonomic level. An example EggNOG-mapper annotation file can be seen below: 

.. code-block::

  ##
  #query  seed_ortholog   evalue  score   eggNOG_OGs      max_annot_lvl   COG_category    Description      Preferred_name  GOs     EC      KEGG_ko KEGG_Pathway    KEGG_Module     KEGG_Reaction    KEGG_rclass     BRITE   KEGG_TC CAZy    BiGG_Reaction   PFAMs
  GCA_001563995.1_00002   694440.JOMF01000005_gene106     3.8e-27 129.4   COG2064@1|root,arCOG01808@2157|Archaea,2XT4Q@28890|Euryarchaeota,2NAP1@224756|Methanomicrobia    2157|Archaea     N       Type II secretion system (T2SS), protein F      -       -       -       ko:K07333        -       -       -       -       ko00000,ko02035,ko02044 -       -       -T2SSF
  GCA_001563995.1_00012   1343739.PAP_05820       4e-147  528.5   COG1855@1|root,arCOG04116@2157|Archaea,2XSZY@28890|Euryarchaeota,242R1@183968|Thermococci        2157|Archaea    VK homology RNA-binding domain   -       GO:0005575,GO:0005618,GO:0005623,GO:0030312,GO:0044464,GO:0071944        -       ko:K06865       -       -       -       -       ko00000 --       -       Intein_splicing,KH_1,KH_2,LAGLIDADG_3,PIN,PIN_4,T2SSE

Function definition
####################
The second part of the configuration file contains the definition of the
function of interest. Functions are defined in the YAML format in a hierarchical
structure. An example of a function definition can be seen below:

.. code-block::

    ---
    name: Equol Gene Cluster
    components:
    - name: Equol Production Pathway
      presence: essential
      components:
      - geneID: DHDR
        description: Dihydrodaidzein reductase
        presence: essential
        terms:
        - id: GCF_000422625.1_00043
          method: blast
          ident_pct: 90
          evalue: 0.00001
      - geneID: THDR
        description: Tetrahydrodaidzein reductase
        presence: essential
        terms:
        - id: COG1053
          method: emapper




Functions are defined in a nested structure. Each component is
defined with a name property, an optional description property, and
a presence property which defines if that component is essential or
nonessential for the overall function.

======================  ========================================================
Name                    Description
======================  ========================================================
name/geneID:(*str*)      The name of the components/ The gene ID
----------------------  --------------------------------------------------------
components:(*list*)      The list of subcomponents
----------------------  --------------------------------------------------------
presence:(*option*)     "essential", "nonessential"
----------------------  --------------------------------------------------------
analogs:(*dict*)        Followed an equivalent component
======================  ========================================================


Components are ultimately associated with geneIDs, which have the same
set of properties as higher level components, but also have search terms
associated with them. In the example below the geneID 'DHDR' is associated
with a sequence as a search term:

.. code-block::

  - geneID: DHDR
    description: Dihydrodaidzein reductase
    presence: essential
    terms:
    - id: GCF_000422625.1_00043
      method: blast
      ident_pct: 90
      evalue: 0.00001

Search terms consist of a search term ID, the method associated with searching
for this term, and additional filtering parameters. Any of the filtering parameters
applicable to a given search term can be set for individual search terms in this
way. See the configuration settings in the above sections for info on what
filtering parameters are applicable for each approach.
