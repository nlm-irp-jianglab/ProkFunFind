Queries and Search Genomes
**************************

This tutorial section will cover how to format the queries and the
information for the genomes that you are going to search.

# Queries

The query format used by *ProkFunFind* is designed around the concept of what we
refer to as a *feature model*. This feature model is a collection of search
terms that are associated with a biological function of interest. These terms
are organized into a hierarchical structure that is used to represent the
relationships between different components of the biological system.

The tutorial materials for this section are located in the
./queries/ directory.

## Types of Search Terms
*ProkFunFind* supports multiple types of queries and additional support for new search approaches is actively being worked on.

You can perform searches using *ProkFunFind* with the following kinds of
search terms:

| Type of Search | Search Term |
|----------------|-------------|
| Protein Sequence | Amino Acid Sequence |
| Hidden Markov Model | Protein Profile HMMs |
| Protein Domains | Supported domains predicted by InterProScan (SEE: [InterProScan]). Including Pfam and TIGRFAM|
| Ortholog Groups | Kegg Orthology (KO) and Clusters of Orthologous Groups (COGs) |

[InterProScan]: https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html#included-analyses

Multiple queries can be associated with a single search term ID. For example
one component of your feature can be associated with multiple protein sequences,
a COG, ID, and a set of protein domains. These IDs are linked to the individual
search terms names through a two column table. This table can be seen in the
`search-term.tsv` file:

```
  enzyme1	Sequence1	blast
  enzyme1	PF07969	interproscan	1
  enzyme2	Profile1	hmmer
  enzyme2	COG0670	emapper
  enzyme2	K00554	kofamscan
```
This input table must have three columns, with an optional fourth column. The
first column specifies the search term ID, the second column specifies the
query associated with the search term, and the third column specifies the
search approach (blast for sequences, hmmer for profile HMMs, interproscan for
protein domains, emapper for COGs, or kofamscan for KOs). The last column can
be used to specify the domain precision *****


## Feature Model Definition
The search terms in *ProkFunFind* are organized in a hierarchical representation
of the biological function. A visual representation of this can be seen here:

![](./figures/feature-model.png)

At the highest level of this organization is the overall biological feature
that you want to search for. This would typically be something like a biological
pathway or an enzyme complex, but this structure can be used to represent any
collection of features.

The overall feature is broken up into subcomponents. These can be any groupings
that you want to use to organize your search terms. An example of a subcomponent
of a metabolic pathway would be an enzyme complex that catalyzes one of the
metabolic reactions. This format is also flexible, allowing for the Definition
of multiple levels of subcomponents.

The subcomponents of a feature are ultimately associated with one or more
search term IDs. When a search is performed using *ProkFunFind* the presence
or absence of a the feature is assessed based on the specified essentiality
and what search terms were detected in the search.

The feature model is provided to the *ProkFunFind* pipeline as a JSON formatted
file. An example of a feature model input for the toy example in the above
image can be seen in the `systems.json` file:

```
  {
      "name": "Pathway",
      "components": [
          {
              "name": "Reaction1",
              "presence": "essential",
              "components": [
                  {
                      "orthoID": "enzyme1",
                      "description": "Enzyme catalyzes Reaction1",
                      "presence": "essential"
                  }
              ]
          },
          {
              "name": "Reaction2",
              "presence": "nonessential",
              "components": [
                  {
                      "orthoID": "enzyme2",
                      "description": "Enzyme catalyzes Reaction2",
                      "presence": "essential"
                  }
              ]
          }
      ]
  }
```
This is just a toy example meant to show the format, but more complex
relationships can be represented as well. Examples of a more complex feature
model can be seen in the `systems-complex.json` file and additional examples
can also be seen in the other tutorial sections ({doc}`Sequence Searches <4-seqsearch>`
and {doc}`Annotation Searches <5-annotsearch>`).

# Search Configuration
The configuration file is where all of the search parameters are defined. This
central file should be named 'config.ini'. This file is broken up into different
sections where filtering thresholds and file naming patterns are defined. An
example of a config.ini file can be found in the `config.ini` file in the examples
directory.

```
  [main]
  cluster.tool   = DBSCAN
  system.file    = system.json
  faa_suffix     = .faa
  gff_suffix     = .gff3
  fna_suffix     = .fna
  search_terms = domain_precision.txt

  [DBSCAN]
  cluster.eps         = 4
  cluster.min_samples = 1.8

  [emapper]
  annot_suffix = .emapper.annotations
  evalue = 1e-3
  filter_file = hit_filter.tab

  [kofamscan]
  annot_suffix = .kofam.tsv
```

This example file is set up for running a search using EGGNog-mapper and
KOfamscan annotation results.

The 'main' section of the configuration file defines the names of the
feature model definition file in the 'system.file' property and the search
terms file in the 'search_terms' property. This section is also used to set
the clustering tool used to identify if the hits from the search are found in
any clusters within the genome. Currently on the DBSCAN algorithm is supported
for gene clustering. The last settings defined in the main section are the
'faa_suffix', 'gff_suffix', and 'fna_suffix' properties which are used to
specify the file extensions for the amino acid fasta files, gff files, and genome
fasta files respectively.

The 'DBSCAN' section is used to set the parameters used in the DBSCAN clustering
to determine if multiple genes are present in the same clusters in the genome.
The 'clsuter.eps' setting is used to set how far two observations can be while
still being considered to be in the same cluster. the 'cluster.min_samples'
parameter is used to determine how many genes must be in the same region for
them to be considered a cluster. See
`DBSCAN <https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html>`_
for more information on the DBSCAN implementation.

The other sections defined in the configuration file are search approach
specific. These sections are used to define the file extensions for the
annotation or query files as well as the filtering parameters for each search
approach. See the other search specific toturial sections for examples and the
'inputs' section of the documentation for a complete table of all settigns
allowed for each search approach.


# Search Space
The last component of the `ProkFunFind` approach is the genomes being searched.
The set of information needed for each genome depends on what kinds of searches
are being performed. At minimum each genome needs a genome fasta file, a GFF
file containing the predicted genes, and a protein fasta file of the predicted
protein sequences. With just this information searches can be performed using
BLAST or HMMER. To search using additional features, files containing the results
of running EGGNog-mapper, InterProScan, or KOfamscan also need to be present. See
the documentation sections and PFA tutorial section for more information on these
annotation formats and the ProkFunAnnotate pipeline that can be used to generate
them.
