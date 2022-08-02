Annotation-based Searches
*****************************

This tutotorial section focuses on performing searches based on precomputed
annotation data for a genome. These methods allow you to take an annotation
file for a genome of interest and perform searches on that using common
annotation features including KEGG Orthology (KO) identifiers, Clusters of
Orthologous genes (COGs) identifiers, and protein domain IDs.

This section of the tutorial will focus on how to run searches for
each kind of annotation feature individually, but they can be mixed
together as well to allow for more complex searches. These type of
mixed query searches are covered in the Mixed Search tutorial section.

# Protein Domain-Based Searches
Protein domains are structural or functional regions of a larger protein. Knowing
what domains are present on a predicted protein can provide valuable clues about
what that proteins functions might be. Similarly, searching for protein domains
within the proteins predicted for a genome can help with the idenfication of
specific genes and pathways.

In ProkFunFind we currently support annotations of protein domains using the
`InterProScan <interproscan-docs.readthedocs.io/en/lastest/index.html>`.
ProkFunFind can incoperate these annotations in the tsv output format from
InterProScan. An example of InterProScan output for the tutorial genomes can
be seen in the `./genomes/GTDB18040_InterProScan.tsv`:

.. code-block::

   GCF_000478885.1_00001	64a2914933b3ef78b443cd67acb6aabf	542	Pfam	PF00308	Bacterial dnaA  protein	196	417	3.8E-65	T	14-10-2020	IPR013317	Chromosomal replication initiator protein DnaA
   GCF_000478885.1_00001	64a2914933b3ef78b443cd67acb6aabf	542	Pfam	PF08299	Bacterial dnaA protein helix-turn-helix	450	516	6.1E-24	T	14-10-2020	IPR013159	Chromosomal replication initiator, DnaA C-terminal	GO:0005524|GO:0006270|GO:0006275|GO:0043565
   GCF_000478885.1_00001	64a2914933b3ef78b443cd67acb6aabf	542	Pfam	PF11638	DnaA N-terminal domain	7	71	5.2E-11	T	14-10-2020	IPR024633	DnaA N-terminal domain


This format consists of multiple columns, but the essential ones
for ProkFunFind are going to be column 1, the protein acession, and column 5,
the domain signature ID. ProkFunFind supports searching by any of the domain
signature types reported by InterProScan.

## Queries
The query files for this search are included in the `./queries/ipr-search/`
directory. Unlike the sequence and profile based searches, an additional query
file is not required for domain-based searches. Instead the domains of
interest are defined directly in the search terms file. In this example the
search terms are defined in the `./queries/ipr-search/search-terms.tsv` file.

.. code-block::

  DZNR	PTHR42917	interproscan
  DEVR	PF00196	interproscan
  DDRC	PF13669	interproscan
  DHDR	PF13561	interproscan
  FIXC	PF03486	interproscan
  FIXX	PTHR43082	interproscan
  HYDF	PF01926	interproscan
  HNDD	PF02256	interproscan
  HYDE	PF04055	interproscan
  HYDG	PTHR43583	interproscan
  NADO	PF07992	interproscan
  THDR	PF00890	interproscan
  FIXA	PF01012	interproscan
  FIXB	PF00766	interproscan


In this file the query ID is provided in column 1, the specific search term
ID (in this case a protein domain accession) is provided in column 2, and the
search method (interproscan) is provided in column 3. As described in the Queries
tutorial section, multiple domains can be associated with one query ID by
providing them on multiple lines.

## Search
The domain-based search can be run with the following command from the base
directory of the tutorial repository.

.. code-block::

  prokfunfind -d queries -f ipr-search -g ./genome-list.tsv --outputprefix ./out/ipr-search/ipr

The output printed to the screen summarizes the overall presence and abseance of
different components in the function definition:

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
  Detected function: Equol Gene Cluster in genome ./genomes//GTDB26128
  1 out of 1 essential components present
  3 out of 3 nonessential components present

In this search we find that both genomes have all components of the function
detected. If you have done the previous tutorial section related to
`Sequence-based Searches <4-seqsearch.rst>` then you will have seen that these
genes actually do not appear to be present in the second, GTDB16128, genome.
The reason why these results differ when basing the search largely off of PFAM
domains is that a domain-based search can be less specific than a sequence or
profile based search. Multiple different proteins can have the same domain,
especially when dealing with domains that are related to common functions like
cofactor binding or dehydrogenase functions. In the next section we will take
look at the search output and see if there are any clues about that we could
use to interpret this result.

## Output

In this example we searched for a set of genes related to a the production
of the metabolite equol. It is often common for genes related to the same
function to be found in the same region of the genome, and based on what is
known about the equol production genes this is true for this set of genes.
We can use this information to help parse through the search results and
get a better idea of if we actually find the gene cluster in both genomes.

The best output to look at for this is going to be the `./out/ipr-search/*.tsv`
outputs. Specifically you want to look at the clusters that the genes are put
into in this output. ProkFunFind uses the DBSCAN algorithm to cluster putative
gene hits into groups on the genome. This grouping is based on the distance
(in number of genes) between two putative hits, meaning that genes that the
clusters that are determined are made of up genes that are in close proximity
to each other.

If you look at the `./out/ipr-search/ipr.GTDB18040.tsv` you will see the cluster
information in the second column. This cluster information gives the
contig ID of the genome assembly and the cluster ID. Genes with the
ID Cl_NA were not found to be part of a cluster. The other clusters are assinged
numerical IDs based on the order they are present in the genome.

The clusters identified in the GTDB18040 tend to be small, consisting of only
3 to four genes in many cases. But There is one cluster that consists of
all 15 genes from the search:

.. code-block::

  GCF_000478885.1_02267	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/hydrogenase maturase/HYDF
  GCF_000478885.1_02268	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/hydrogenase maturase/HYDG
  GCF_000478885.1_02269	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/hydrogenase maturase/HYDE
  GCF_000478885.1_02270	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/hydrogenase maturase/HNDD
  GCF_000478885.1_02271	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/hydrogenase maturase/NADO
  GCF_000478885.1_02272	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/fix electron transport/FIXX
  GCF_000478885.1_02273	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/fix electron transport/FIXC
  GCF_000478885.1_02274	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/Equol Production Pathway/DZNR
  GCF_000478885.1_02276	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/Equol Production Pathway/DHDR
  GCF_000478885.1_02277	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/Equol Production Pathway/THDR
  GCF_000478885.1_02278	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/fix electron transport/FIXB
  GCF_000478885.1_02279	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/fix electron transport/FIXA
  GCF_000478885.1_02280	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/Equol Production Pathway/DDRC
  GCF_000478885.1_02281	GCF_000478885.1_1:Cl_35	Equol Gene Cluster/other genes/DEVR

In contrast when looking at the output for the GTDB26128 genome,
`./out/ipr-search/ipr.GTDB26128.tsv`, you can see that a majority of the
clusters are small and even the largest ones like Cl_28 consist of multiple hits
to the same genes. This provides an indication that despite putative hits to
all of the genes being identified, there do not seem to be any 'real' looking
clusters.

.. code-block::

  GCF_011405655.1_01933	GCF_011405655.1_1:Cl_28	Equol Gene Cluster/Equol Production Pathway/THDR
  GCF_011405655.1_01934	GCF_011405655.1_1:Cl_28	Equol Gene Cluster/other genes/DEVR
  GCF_011405655.1_01935	GCF_011405655.1_1:Cl_28	Equol Gene Cluster/Equol Production Pathway/THDR
  GCF_011405655.1_01936	GCF_011405655.1_1:Cl_28	Equol Gene Cluster/other genes/DEVR
  GCF_011405655.1_01939	GCF_011405655.1_1:Cl_28	Equol Gene Cluster/other genes/DEVR
  GCF_011405655.1_01943	GCF_011405655.1_1:Cl_28	Equol Gene Cluster/other genes/DEVR
  GCF_011405655.1_01944	GCF_011405655.1_1:Cl_28	Equol Gene Cluster/Equol Production Pathway/THDR


The use of this clustering information to identify high quality putative hits
is highly dependent on the features being searched. While genes being in the
same gene cluster can be an indication of related function, this is not always
true. Many metabolic pathways consist of genes that are not found in the
same gene cluster, so your interpretation of these results may vary based on
your scientific question.


# KEGG Orthology-Based Searches
The KEGG database groups genes into manually defined functional ortholog groups.
The KO database has become a popular resource to link genes to their functions
within larger metabolic pathways and subsystems. For more information on the
KO database see `KEGG Ortholog <genome.jp/kegg/ko.html>`.

In ProkFunFind the KO assignments are parsed from KofamScan tabular output. An
example of this output for the tutorial genomes can be seen in `./genomes/GTDB18040.kofam.tsv`:

.. code-block::

   *	GCF_000478885.1_00001	K02313	130.33	443.5	1.7e-133	"chromosomal replication initiator protein"
  	  GCF_000478885.1_00001	K10763	171.70	89.6	2.1e-26	"DnaA-homolog protein"
  	  GCF_000478885.1_00001	K02315	138.67	64.8	8.1e-19	"DNA replication protein DnaC"


## Queries
KO based searches are done using KO identifiers as search terms. More information
on how KO identifiers are assigned and full references of all KO identifiers please
see the KEGG database here: `KEGG <https://www.genome.jp/kegg/ko.html>`.

For this query KO identifiers for each of the components of the equol gene clusters
were assigned KO identifiers. This can be seen in the
`./queries/kofam-search/search-terms.tsv` file:

.. code-block::

  HYDF	K03977	kofamscan
  HYDG	K03150	kofamscan
  HYDE	K01012	kofamscan
  HNDD	K18332	kofamscan
  FIXX	K03855	kofamscan
  FIXC	K00313	kofamscan
  DZNR	K00219	kofamscan
  DEVR	K07695	kofamscan
  DDRC	K05606	kofamscan
  DHDR	K18009	kofamscan
  FIXB	K03522	kofamscan
  FIXA	K03521	kofamscan
  HYPO	K02004	kofamscan
  NADO	K15022	kofamscan
  THDR	K00244	kofamscan

Not all of the genes being used in the query for this tutorial are have great
matches to the current KO groups defined by KEGG. Because of this you also have
to make the search a little more lenient by adjusting the threshold filtering
property in the `./queries/kofam-search/config.ini` `[kofamscan]` section:

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

  [kofamscan]
  annot_suffix = .kofam.tsv
  threshold = 0.5


For the KO assignment in kofamscan, a match score is calculated for each gene
to KO pair. This score is then compared to an predetermined score for
each KO. The threshold parameter allows you to adjust that score requirement.
The score will be multiplied by the value provided in the threshold argument,
requiring either a higher or lower score for a KO assignment. In this case
setting the threhsold parameter to 0.5 would make the score half as strict.
This score threshold and the evalue parameter may need to be adjusted in
different searches to fine tune your search, especially when there are not
great KO matches for your genes of interest.

## Search

The KO based search can be done from the root directory of the tutorial
repository using the following command.

.. code-block::

  prokfunfind -d queries -f kofam-search -g ./genome-list.tsv --outputprefix ./out/kofam-search/kofam


Based on this search we can detect all four components in the first genome,
but only the three non-essential components in the second genome:

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

The output for this type of search is the same as the other approaches providing
information about the putative gene hits and clusters of genes found during the
search.

# COG-Based Searches

ProkFunFind also supports searching by Clusters of Orthologous Genes (COGs). COGs
are widely used ortholog groupings. For ProkFunFind searches we use EGGNog-mapper
as the annotation tool to assign COGs. The pregenerated output for this tutorial
can be seen in the `./genomes/*.emapper.annotations` files:

.. code-block::

   GCF_011405655.1_00003	1384484.AEQU_2159	2.78e-73	220.0	COG2198@1|root,COG2198@2|Bacteria,2HVH2@201174|Actinobacteria,4CWUG@84998|Coriobacteriia	2|Bacteria	T	Hpt domain	-	-	-	-	-	-	-	-	-	-	-	-	Hpt
   GCF_011405655.1_00004	1384484.AEQU_2160	0.0	1390.0	COG2199@1|root,COG3437@1|root,COG2199@2|Bacteria,COG3437@2|Bacteria,2I49F@201174|Actinobacteria,4CUE6@84998|Coriobacteriia	2|Bacteria	T	HD domain	-	-	-	ko:K07814	-	-	-	-	ko00000,ko02022	-	-	-	GGDEF,HD,HD_5,Response_reg
   GCF_011405655.1_00005	1384484.AEQU_2161	9.13e-303	825.0	COG1541@1|root,COG1541@2|Bacteria,2GJC7@201174|Actinobacteria,4CUT2@84998|Coriobacteriia	2|Bacteria	H	AMP-binding enzyme C-terminal domain	paaK-3	-	6.2.1.30	ko:K01912	ko00360,ko01120,ko05111,map00360,map01120,map05111	-	R02539	RC00004,RC00014	ko00000,ko00001,ko01000	-	-	-	AMP-binding,AMP-binding_C_2

The orthology assignments in this output can be seen in the fifth column of the output.
This column gives ortholog assignments at different taxonomic levels in this output
and any of these IDs can be used to search through ProkFunFind.

## Query

Similar to the KO-based search, the COG based searches define the queries based on the ortholog IDs, in
this case COG IDs. The search term input can be found in the `./queries/emap-search/search-terms.tsv` file:

.. code-block::

  HYDF	COG1160	emapper
  HYDG	2IKFZ	emapper
  HYDE	2HSHP	emapper
  HNDD	COG3383	emapper
  FIXX	COG2440	emapper
  FIXC	COG0644	emapper
  DZNR	COG1902	emapper
  DEVR	COG2197	emapper
  DDRC	COG0346	emapper
  DHDR	COG1028	emapper
  FIXB	COG2025	emapper
  FIXA	COG2086	emapper
  HYPO	29K0U	emapper
  NADO	COG1894	emapper
  THDR	COG1053	emapper

Similarly to the KO-based search, many of the queries in this example search do not have great COG
matches, so a mix of COGs and ortholog groups at higher levels are used in this search.

Additionally, because ortholog groups can have varying levels of specificity and our search terms are
not perfect matches to each COG group this search will be performed using an additional search term
specific filtering file. This kind of input file can be used to add individual filtering parameters
to the search, for example setting different evalue thresholds for different COGs.

The filtering file can be found in the `./queries/emap-search/filter.tsv` file:

.. code-block::

  COG1160	evalue	<=	1e-100
  2IKFZ	evalue	<=	1e-100
  2HSHP	evalue	<=	1e-200
  COG3383	evalue	<=	1e-100
  COG2440	evalue	<=	1e-70
  COG0644	evalue	<=	1e-100
  COG1902	evalue	<=	1e-100
  COG2197	evalue	<=	1e-80
  COG0346	evalue	<=	1e-100
  COG1028	evalue	<=	1e-200
  COG2025	evalue	<=	1e-150
  COG2086	evalue	<=	1e-150
  29K0U	evalue	<=	1e-100
  COG1894	evalue	<=	1e-250
  COG1053	evalue	<=	1e-100

## Search
The search can be performed using the following command:

.. code-block::

  prokfunfind -d queries -f emap-search -g ./genome-list.tsv --outputprefix ./out/emap-search/emap

Based on this search it can be seen that the components of the function were detected
in both genomes:

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
  Detected function: Equol Gene Cluster in genome ./genomes//GTDB26128
  1 out of 1 essential components present
  3 out of 3 nonessential components present

What happens is similar to the issue seen in the domain-based
search, where we have non-specific hits to additional genes in the second genome.
You can check the tsv or GFF ouput in the `./out/emap-search/` directory to
confirm this by looking for larger clusters of putative hits on both genomes.
This does highlight one of the benefits of using the ProkFunFind search tool to
perform mixed searches using combinations of different approaches. A walkthrough
on how to set up and run those searches can be found in the `Mixed Search <6-mixedsearch.rst>`
tutorial section. 
