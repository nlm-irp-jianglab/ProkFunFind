Sequence and Profile Searches
*****************************

This tutorial section will walk you through how to set up a search using BLAST or
HMMER within ProkFunFind. The materials for this section of the tutorial are
going to be in the `queries/blast-search/`, `queries/hmmer-search/`, and
the `genomes/` directories of the prokfunfind-tutorial repository. Pre-generated
output in the

# BLAST-based Searches
## Input
To perform a BLAST-based search with ProkFunFind you will need a protein
sequences for your search terms and need the genome fasta, protein fasta, and
GFF files for your search genomes.

For this tutorial the genome information is already prepared in the `genomes/`
directory. In this directory there are two genomes from the GTDB database for
*Adlercreutzia equolifaciens* and *Adlercreutzie celatus_A* named GTDB18040
and GTDB26128 respectively. Both of these genomes have pre-generated
annotation information. For this blast search the essential files are the
'.gff', '.fna', and '.faa' files.

The query for the search consists of all of the standard configuration file and
function defintion along with a protein fasta file containing the query sequences.
This protein fasta file needs to be indexed using the `makeblastdb` tool before
performing a search. These input files are in the 'queries/blast-search/'
directory.

## Equol Gene Cluster Query
The query for this tutorial is a set of features that can be used to identify
a gene cluster that encodes a biosynthesis pathway for the metabolite S-Equol
and multiple other genes that have been observed to be in the same gene cluster.
The function definition is broken up into four components with one of them being
essential for the function to be considered present.
