.. ProkFunFind

.. _ebi:


*************************************
Using Precomputed EBI Annotation Data
*************************************

The European Molecular Biology Laboratory-European Bioinformatics Institute (EMBL-EBI)
provides additional annotation information with all of its sepcies genomes in the MGnify
database. these can be accessed through their FTP site here: 
`MGnify FTP <http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/>`
Each genome in the MGnify Genomes species catalogues contains all of the 
information needed to perform searches by protein sequences, profile HMMs,
protein domains, and COG IDs. 
An example can be seen with the `MGYG000000101 <http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/species_catalogue/MGYG0000001/MGYG000000101/genome/>` genome in the human gut
species collection. The files included with this genome are: 

.. code-block:: 

   MGYG000000101.faa *Protein Fasta File
   MGYG000000101.fna *Genome Fasta File
   MGYG000000101.fna.fai
   MGYG000000101.gff *Genome GFF File
   MGYG000000101_InterProScan.tsv *InterProScan Domain Predictions
   MGYG000000101_annotation_coverage.tsv
   MGYG000000101_cazy_summary.tsv
   MGYG000000101_cog_summary.tsv
   MGYG000000101_eggNOG.tsv *EGGNog-mapper Predictions
   MGYG000000101_kegg_classes.tsv
   MGYG000000101_kegg_modules.tsv
   MGYG000000101_rRNAs.fasta

The input files needed ProkFunFind searches are marked with '*' in the
above file list. 

# EBI Search Configuration Files
The EBI genomes are all formatted in the same way, making it easy to
set up a standard configuration file to use with any search being
done on these genomes. 

The configuration (config.ini) file would look something like this: 

.. code-block:: 

   [main]
   system.file = system.json
   cluster.tool = DBSCAN
   search_terms = search_terms.tsv
   faa_suffix = .faa
   fna_suffix = .fna
   gff_suffix = .gff3

   [DBSCAN]
   cluster.eps = 3
   cluster.min_samples = 1

   [blast]
   blast.threads = 1
   blast.query = query.fa

   [hmmer]
   hmmer.threads = 1   
   hmmer.query = query.hmm

   [emapper]
   annot_suffix = _eggNOG.tsv

   [interproscan]
   annot_suffix = _InterProScan.tsv

With this configuration file setup you can perform searches on any EMBL-EBI
genome data that you download. 

For a working example of a search using an EMBL-EBI genome please see the
'Working with Genome Databases' section of the tutorial here: ` <>`
   
   
