*****************************
Mixed Searches
*****************************

This tutorial section will cover how to perform searches using multiple types
of search terms. In the previous tutorial examples you performed searches using
only one type of supported search term in each search. While these types of
searches can be useful, in multiple of the tutorial examples you saw how using
just one type of search term can result in false positive or false negative
search results. A benefit of using the ProkFunFind search approach comes from
the ability to utilize multiple different types of queries during a search. For
example one gene may be best represented by a profile HMM while another may
be part of a well defined KO or COG group.

# Mixed Search Terms
In ProkFunFind mixed searches you can incorporate any of the support search
term types in any combination. The only limitation is based on what
files you have available for each genome. For example in order to perform searches
based on protein sequences and COG identifiers you will need to have the protein
fasta files and EGGNog-mapper results available for each genome being searched.

## Queries

The two most important files for configuring a mixed search in ProkFunFind are
the configuration file and search-terms file. The configuration file is going to
be formatted like it has been in previous examples, but it needs to include
sections for each of the search term types that are being used. An example can
be seen in the `./queries/mixed-search/config.ini` file:

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

  [kofamscan]
  annot_suffix = .kofam.tsv
  threshold = 0.5

  [emapper]
  annot_suffix = .emapper.annotations

  [interproscan]
  annot_suffix = _InterProScan.tsv


The other important file is the search-terms file. In this file each query ID
is associated with the individual search terms and the search approaches. In
this file multiple search terms can be associated with the same query, allowing
for queries to be identified through multiple approaches simultaneously. The
search terms file for this tutorial search can be seen in the `./queries/mixed-search/search-terms.tsv`:

.. code-block::

  DZNR	DZNR	hmmer
  DEVR	DEVR	hmmer
  DDRC	DDRC	hmmer
  DHDR	GCF_000422625.1_00043	blast
  FIXX	K03855	kofamscan
  FIXC	K00313	kofamscan
  HYDF	PF01926	interproscan
  HNDD	PF02256	interproscan
  HYDE	PF04055	interproscan
  HYDG	PTHR43583	interproscan
  NADO	COG1894	emapper
  THDR	COG1053	emapper
  FIXB	COG2025	emapper
  FIXA	COG2086	emapper

This search is going to use a mix of 3 profile HMMs, 1 protein sequence, 2 KOs,
4 domain signatures, and 4 COGs.

## Search
To perform the search from the root directory of the tutorial repository you can
run the following command:

.. code-block::

  prokfunfind -d queries -f mixed-search --gtab ./genome-list.tsv --outputprefix ./out/mixed-search/mixed

This command will return an initial summary of the component presence and
absence:

.. code-block::

  INFO:root:Checking configuration files
  INFO:root:Searching for function
  INFO:root:Identifying gene clusters
  INFO:root:Summarizing function presence and genes
  Detected function: Equol Gene Cluster in genome ./genomes//GTDB18040
  1 out of 1 essential components present
  3 out of 3 nonessential components present
  INFO:root:Searching for function
  INFO:root:Identifying gene clusters
  INFO:root:Summarizing function presence and genes
  Failed to detect function: Equol Gene Cluster in genome ./genomes//GTDB26128
  0 out of 1 essential components present
  3 out of 3 nonessential components present

## Output
The output is the same as what is produced by other searches. Because the
search is done using multiple search terms it can also be useful to check the
output to see what search terms are producing hits to certain genes in the
results. This information can be found in the gff output of the search. For
this search the output can be seen in the `./out/mixed-search/mixed.GTDB18040.annot.gff`
file:

.. code-block::

  GCF_000478885.1_1	GuFunFind	CDS	7382	8305	.	-	.	ID=GCF_000478885.1_00007;Name=HYDE;ClusterID=Cl_NA;Target=PF04055;evalue=3.2e-19
  GCF_000478885.1_1	GuFunFind	CDS	28201	29646	.	+	.	ID=GCF_000478885.1_00024;Name=HYDE;ClusterID=Cl_NA;Target=PF04055;evalue=2e-15
  GCF_000478885.1_1	GuFunFind	CDS	261233	262642	.	+	.	ID=GCF_000478885.1_00150;Name=HYDE;ClusterID=Cl_NA;Target=PF04055;evalue=7.7e-22
  GCF_000478885.1_1	GuFunFind	CDS	288712	290358	.	-	.	ID=GCF_000478885.1_00174;Name=DEVR;ClusterID=Cl_NA;Target=DEVR;evalue=1.2e-07

In this output the Target property in column 9 provides what specific search term
ID produced the hit to that gene. 
