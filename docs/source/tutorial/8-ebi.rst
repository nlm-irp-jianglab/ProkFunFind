*************************
Using Existing Databases
*************************

Various databases have started to include precomputed and standardized
annotation files alongside their genome collections. These can be used with
ProkFunFind as well, providing a convenient way to survey for functions among
these large genome collections.

One example of this can be seen in the EMBL-EBI MGnify Genomes collections. This
collection can be found here:
`MGnify Genomes FTP <http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/>`_
Within these genomes collections each of the genomes included in the species
catalogues (for example this one: `MGYG000000001 <http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/species_catalogue/MGYG0000000/MGYG000000001/genome/>`_).

The EBI species catalogue genomes come with EGGNog-mapper output in the `*_eggNOG.tsv`
file and InterProScan predictions in the `*_InterProScan.tsv` file. This information
along with the genome fasta, protein fasta, and gff files included for each genome
let us perform searches with the EBI genomes based on sequences, profile HMMs,
protein domains, and COGs.

Setting up configuration files for EBI Genomes
###############################################
The configuration files that would be used for EBI genome searches are largely
the same as the ones that have been used for the other searches in this tutorial.
The only major changes needed compared to the other examples are altering the
file extensions to fit the conventions used by EBI. These changes are going
to be made in the config.ini file.

An set of query files used to search for the equol gene cluster have been set up
in the `./queries/ebi-search/` directory. The `./queries/ebi-search/config.ini`
file has the updated extension for the EGGNog-mapper table
for the EBI file formats:

.. code-block::

  [main]
  cluster.tool = DBSCAN
  system.file = systems.json
  search_terms = search-terms.tsv
  faa_suffix = .faa
  gff_suffix = .gff
  fna_suffix = .fna

  [DBSCAN]
  cluster.eps = 4
  cluster.min_samples = 2

  [hmmer]
  hmmer.query = query.hmm
  hmmer.exec = hmmscan
  hmmer.threads = 1
  evalue = 1e-3

  [blast]
  blast.query = query.fa
  blast.exec = blastp
  blast.threads = 1
  evalue = 1e-3

  [emapper]
  annot_suffix = _eggNOG.tsv

  [interproscan]
  annot_suffix = _InterProScan.tsv

This configuration file sets up a search for the equol gene cluster using a
mixture of COGs, profile HMMs, sequences, and protein domains designated in
the `./queries/ebi-search/search-terms.tsv` file:

.. code-block::

  DZNR	DZNR	hmmer
  DEVR	DEVR	hmmer
  DDRC	DDRC	hmmer
  DHDR	GCF_000422625.1_00043	blast
  FIXX	GCF_000422625.1_00039	blast
  FIXC	GCF_000422625.1_00040	blast
  HYDF	PF01926	interproscan
  HNDD	PF02256	interproscan
  HYDE	PF04055	interproscan
  HYDG	PTHR43583	interproscan
  NADO	COG1894	emapper
  THDR	COG1053	emapper
  FIXB	COG2025	emapper
  FIXA	COG2086	emapper

In this file the search term ID is specified in the first column,
the specific query item (sequence ID, profile ID, domain, or COG) is
specified in the second column, and the search method is specified in
the last column.

Search
######
An example search can be performed from the root directory of the tutorial
repository using the following command:

.. code-block::

  prokfunfind-tutorial % prokfunfind -d queries -f ebi-search --gtab ./ebi-list.tsv --outputprefix ./out/ebi-search/ebi

The results of this search show that only 2 of the non-essential components of
the functions were found in the search:

.. code-block::

  INFO:root:Checking configuration files
  INFO:root:Searching for function
  INFO:root:Identifying gene clusters
  INFO:root:Summarizing function presence and genes
  Failed to detect function: Equol Gene Cluster in genome ./genomes//MGYG000000001
  0 out of 1 essential components present
  2 out of 3 nonessential components present

The remaining output is the same as has been detailed in the previous tutorial
sections.
