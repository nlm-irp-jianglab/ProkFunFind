*****************************
Sequence and Profile Searches
*****************************

This tutorial section will walk you through how to set up a search using BLAST or
HMMER within ProkFunFind. The materials for this section of the tutorial are
going to be in the `queries/blast-search/`, `queries/hmmer-search/`, and
the `genomes/` directories of the prokfunfind-tutorial repository. Pre-generated
output is provided for all commands in the `precomputed-out/` directory

These two search approaches on their own are no different from running BLAST or
hmmscan and then applying different filtering to the results. In this case
ProkFunFind is just providing an easy way to perform the search and facilitate
the filtering within on platform. A major benefit of having access to these
search approaches within ProkFunFind comes from the ability to combine these
searches with other annotation-based searches. This process will be detailed
in the later :doc:`Mixed Searches <./6-mixedsearch>` section of the tutorial.

BLAST-based Searches
######################

Input
*****
To perform a BLAST-based search with ProkFunFind you will need a protein
sequences for your search terms and need the genome fasta, protein fasta, and
GFF files for your search genomes.

For this tutorial the genome information is already prepared in the `genomes/`
directory. In this directory there are two genomes from the GTDB database for
*Adlercreutzia equolifaciens* and *Adlercreutzie celatus_A* named GTDB18040
and GTDB26128 respectively. Both of these genomes have pre-generated
annotation information. For this blast search the essential files are the
'.gff', '.fna', and '.faa' files.

The query for the search consists of the standard configuration file and
function defintion along with a protein fasta file containing the query sequences.
The protein fasta file needs to be indexed using the `makeblastdb` tool before
performing a search (this has already been done for these sequences).
These input files are in the `queries/blast-search/` directory.

Search
******
To perform a BLAST-based search for the equol gene cluster you can use the following
ProkFunFind command from the base directory of the tutorial repository.

.. code-block::

   prokfunfind -d queries -f blast-search -g ./genome-list.tsv --outputprefix ./out/blast-blast/blast

This command will first print an overall summary of how many function
components were detected in each genome. This output will summarize how many
essential function components were detected and how many nonessential components
were detected:

.. code-block::

  INFO:root:Checking configuration files
  INFO:root:Searching for function
  INFO:root:Identifying gene clusters
  INFO:root:Summarizing function presence and genes
  Detected function: Equol Gene Cluster in genome ./genomes/GTDB18040
  1 out of 1 essential components present
  3 out of 3 nonessential components present
  INFO:root:Searching for function
  INFO:root:Identifying gene clusters
  INFO:root:Summarizing function presence and genes
  Failed to detect function: Equol Gene Cluster in genome ./genomes/GTDB26128
  0 out of 1 essential components present
  1 out of 3 nonessential components present

The printed output does not provide the specific information
about what components were detected or missing though, this information is
instead provided in the multiple output files that are generated for each
genome. From this information you can see that in the first genome, GTDB18040,
the only essential component of the function and the 3 non-essential components
were all detected in the genome. In contrast the summary for the second genome,
GTDB26128, only had 1 of the non-essential components detected.

This summary can provide a good overview of your search and can be a good first
thing to look at when analyzing collections of genomes, but it should not be
relied on as a final assessment of if a set of genes is present or not in a
genome. Instead you will want to examine the additional output files to get
a better idea of what genes were found, if those genes were found in gene
clusters, and what function components they belonged to.

Output
******
The search command will produce three output files for each genome in the
./out/blast-search/ directory.
* A GFF file
* A TSV file
* A blast output tabular file
* A genome pkl file
* A json summary file of the component presence and absence.

For more information on each of the output formats please see the `Outputs <outputs.rst>`
documentation section. Each tutorial section will provide additional information
about some of the output files.

The first place to look for more specific information on your search is going to
be in the .tsv file. This file provides a more detailed summary of what genes
were detected in the search and what

If you look at the `./out/blast-search/blast.GTDB18040.tsv` file you will see
the following table summary of the genes that were detected:

.. code-block::

  Gene_Name	Cluster_ID	Functions
  GCF_000478885.1_02267	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/hydrogenase maturase/HYDF
  GCF_000478885.1_02268	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/hydrogenase maturase/HYDG
  GCF_000478885.1_02269	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/hydrogenase maturase/HYDE
  GCF_000478885.1_02270	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/hydrogenase maturase/HNDD
  GCF_000478885.1_02271	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/hydrogenase maturase/NADO
  GCF_000478885.1_02272	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/fix electron transport/FIXX
  GCF_000478885.1_02273	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/fix electron transport/FIXC
  GCF_000478885.1_02274	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/Equol Production Pathway/DZNR
  GCF_000478885.1_02275	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/other genes/HYPO
  GCF_000478885.1_02276	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/Equol Production Pathway/DHDR
  GCF_000478885.1_02277	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/Equol Production Pathway/THDR
  GCF_000478885.1_02278	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/fix electron transport/FIXB
  GCF_000478885.1_02279	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/fix electron transport/FIXA
  GCF_000478885.1_02280	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/Equol Production Pathway/DDRC
  GCF_000478885.1_02281	GCF_000478885.1_1:Cl_0	Equol Gene Cluster/other genes/DEVR

For each gene that was detected in the search this table provides information
about the genomic contig and gene cluster they are from and the function component
they were found to be a part of. In this case all of the genes were found to
be in the same cluster, given the ID Cl_0, this indicates that they are found
in the same genomic region.

In contrast if you look at the './out/blast-search/blast.GTDB26128.tsv' you
can see that only one set of genes, the hydrogenase maturase component, was
detected in that genome.


Other things to try
*******************
If you want to explore this search approach more consider doing some of the
following and seeing how it affects the search results for each genome.
* Change the blast filtering parameters in the config.ini file to make them
  more or less stringent.
* Modify some of the essentiality requirements in the systems.json file and
  observe how that changes the search results.


HMMER-based searches
####################
A HMMER-based search using ProkFunFind is performed in a similar way to a BLAST-based
search, but instead of using protein sequences as your search terms, you use
profile HMMs.

Profile HMMs are probabilistic models of the conservation of a set of sequences.
They can be used with tools like `HMMER's` `hmmscan` to perform searches against
protein databases to find sequences that are similar to the profile. For more
information on how to generate and use profile HMMs please see the HMMER
documentation here: `HMMER Docs <eddylab.org/software/hmmer/Userguide.pdf>`_

For this search the query files can be found in the `queries/hmmer-search/`
directory. The queries consist of the standard configuration and systems files,
and the profile HMMs are contained in the `profiles.hmm` file.

Search
*******
To run the profile HMM search within ProkFunFind you can use the following command
from the base directory of the tutorial repository:

.. code-block::

   prokfunfind -d queries -f hmmer-search -g ./genome-list.tsv --outputprefix ./out/hmmer-search/hmmer

The same summary output is printed to the screen as in the BLAST tutorial. This
output provides a simple summary of the component presence and absence in the
genomes being searched.

In this search we see a slightly different result compared to the BLAST search:

.. code-block::

  INFO:root:Checking configuration files
  INFO:root:Searching for function
  INFO:root:Identifying gene clusters
  INFO:root:Summarizing function presence and genes
  Detected function:Equol Gene Cluster in genome ./genomes//GTDB18040
  1 out of 1 essential components present
  3 out of 3 nonessential components present
  INFO:root:Searching for function
  INFO:root:Identifying gene clusters
  INFO:root:Summarizing function presence and genes
  Failed to detect function:Equol Gene Cluster in genome ./genomes//GTDB26128
  0 out of 1 essential components present
  3 out of 3 nonessential components present

In this search you can still detect all four components in the GTDB18040 genome,
but in the second, GTDB26128, genome we detect 3 non-essential components, where
in the blast search we only detected one. This highlights one of the motivations
behind ProkFunFind, demonstrating that single search approaches, for example
just a BLAST search, may not be sufficient to get a full picture of the presence
or absence of functions.

Output
******

For this section of the tutorial we are going to focus on the GFF output files
from the search. These can be seen in the `./out/hmmer-search/*.gff` files. The
output from the first genome can be seen in the `./out/hmmer-search/hmmer.GTDB18040.annot.gff`

.. code-block::

  ...
  GCF_000478885.1_1	GuFunFind	CDS	2774610	2776022	.	-	.	ID=GCF_000478885.1_02267;Name=HYDF;ClusterID=Cl_36;Target=HYDF;evalue=5.9e-242
  GCF_000478885.1_1	GuFunFind	CDS	2776166	2777611	.	-	.	ID=GCF_000478885.1_02268;Name=HYDG;ClusterID=Cl_36;Target=HYDG;evalue=4.2e-306
  GCF_000478885.1_1	GuFunFind	CDS	2777598	2778668	.	-	.	ID=GCF_000478885.1_02269;Name=HYDE;ClusterID=Cl_36;Target=HYDE;evalue=2.5e-222
  GCF_000478885.1_1	GuFunFind	CDS	2778770	2780563	.	-	.	ID=GCF_000478885.1_02270;Name=HNDD;ClusterID=Cl_36;Target=HNDD;evalue=0.0
  GCF_000478885.1_1	GuFunFind	CDS	2780557	2782395	.	-	.	ID=GCF_000478885.1_02271;Name=NADO;ClusterID=Cl_36;Target=NADO;evalue=0.0
  GCF_000478885.1_1	GuFunFind	CDS	2782612	2782923	.	-	.	ID=GCF_000478885.1_02272;Name=FIXX;ClusterID=Cl_36;Target=FIXX;evalue=2.6e-72
  GCF_000478885.1_1	GuFunFind	CDS	2782920	2784233	.	-	.	ID=GCF_000478885.1_02273;Name=FIXC;ClusterID=Cl_36;Target=FIXC;evalue=3.7e-302
  GCF_000478885.1_1	GuFunFind	CDS	2784304	2786232	.	-	.	ID=GCF_000478885.1_02274;Name=DZNR;ClusterID=Cl_36;Target=DZNR;evalue=0.0
  GCF_000478885.1_1	GuFunFind	CDS	2786295	2786774	.	-	.	ID=GCF_000478885.1_02275;Name=HYPO;ClusterID=Cl_36;Target=HYPO;evalue=3.8e-85
  GCF_000478885.1_1	GuFunFind	CDS	2786868	2787716	.	-	.	ID=GCF_000478885.1_02276;Name=DHDR;ClusterID=Cl_36;Target=DHDR;evalue=5.3e-195
  GCF_000478885.1_1	GuFunFind	CDS	2787796	2789259	.	-	.	ID=GCF_000478885.1_02277;Name=THDR;ClusterID=Cl_36;Target=THDR;evalue=0.0
  GCF_000478885.1_1	GuFunFind	CDS	2789323	2790237	.	-	.	ID=GCF_000478885.1_02278;Name=FIXB;ClusterID=Cl_36;Target=FIXB;evalue=1.5e-167
  GCF_000478885.1_1	GuFunFind	CDS	2790267	2790986	.	-	.	ID=GCF_000478885.1_02279;Name=FIXA;ClusterID=Cl_36;Target=FIXA;evalue=2.3e-139
  GCF_000478885.1_1	GuFunFind	CDS	2791008	2791460	.	-	.	ID=GCF_000478885.1_02280;Name=DDRC;ClusterID=Cl_36;Target=DDRC;evalue=2.2e-92
  GCF_000478885.1_1	GuFunFind	CDS	2791670	2792440	.	-	.	ID=GCF_000478885.1_02281;Name=DEVR;ClusterID=Cl_36;Target=DEVR;evalue=1.3e-150
  ...

This output is a standard GFF format table that provides information about the
genes and their locations on the genome, along with the annotation information
related to what components and genomic clusters they are a part of.

The first thing to note is that, in contrast to the BLAST-based search which
only returned one hit per search term, the HMMER search identifies 160 hits.
Depending on the type of function beign searched for you may expect alot of
just a few hits, and the total number of hits may be a good inital way to assess
if your filtering parameters are too strict or too lenient.

Other things to try
*******************
Try adjusting the e-value threshold in the config.ini file to get fewer
hits returned by the search, but still return hits to the actual equol gene
cluster (genes 02268-02281).
