.. ProkFunFind

.. _design:


********************************
Designing a Function Definition
********************************
The function definition is a central feature of ProkFunFind. This
fucntion definition allows for the representation of complex biological
functions in a hierarchical structure that can account for essential
and non-essential components. This section of the documentation will 
contain a description of how to get from a biological concept to
a function defintion and some things that need to be considered along
the way. The structure and format of the YAML file will be described
briefly in this section, but more details are provided in the :ref:`Inputs<inputs>` 
documentation section. 

Designing a YAML Function Definition
#####################################

The design of a function defintion should start with a good understanding of
what the function is that you want to search for. This information can come from
literature, from the analysis of a specific organism that performs the function, or
even be based on a hypothesis about what kinds of protein domains or other features
would be involved in the function. 

For this example lets consider a relatively simple pathway for acetate production that 
is present in many organisms. This pathway consits of the 'pta' and ackA' genes which 
convert Acetyl-CoA to Acetyl-Phosphate and then convert Acetyl-phosphate to acetate, 
producing an ATP. 

Finding Information about a Function
*************************************
To define this function we first need to assess if there are any features that we can 
use to represent each gene, for example a representative sequence, HMM profile model, 
or orthology definition. 

For this function we can look at some commonly used databses to get an idea of what our
genes of interest look like and what features are associated with them. In this case lets 
at KEGG (https://www.genome.jp/kegg/)

====================  =================================================================================================================
Gene                   KEGG Link
====================  =================================================================================================================
pta                    https://www.genome.jp/entry/K13788+K00625+K15024+2.3.1.8+R00230
--------------------  -----------------------------------------------------------------------------------------------------------------
ackA                   https://www.genome.jp/entry/K00925+2.7.2.1+R00315
====================  =================================================================================================================

If we look at both of these pages we can get alot of information about our genes of interest. 
Notably we can find information on the COG and KO orthology groups that these genes are 
parts of. Orthology groups like these are good places to start with when searching for 
functions. In this case we can see that pta is associated with K13788, K15024, COG0280, and COG0857. 
ackA is associated with K00925 and COG0282. 

With this information in hand, we can start to write our function definition out. For this
this example we can start with defining our overall function, in this case we'll call it 
'Acetate Production'. 

.. code-block:: 

    ---
    name: Acetate Production
    components:


This overall function can consist of any number of sub-functions
allowing for the biological system to be flexibly represented. In this
case we only have two genes so we can either defined the genes directly
under the overall functions, or we can define it as part of a sub-function. 
If we were to define it directly under the overall function it would look
something like this: 

.. code-block:: 

    ---
    name: Acetate Production
    components: 
    - geneID: PTA 
    - geneID: ACKA 

For more complex systems though we may want to represent multiple sub-systems, 
or to leave open the possibility that we expand our function definition later on. 
For this purpose, we can define our pta-ackA genes as part of a sub-system, leaving
open the posisbiltiy of adding other sub-systems to the search in the future. This
kind of definition would look like this: 

.. code-block:: 

    ---
    name: Acetate Production
    components:
    - name: PTA-ACKA Pathway
      components:
      - geneID: PTA
      - geneID: ACKA

Customizing the Function Definition
*************************************

Now that we have the overall structure of our function defined, we can start to add 
in the search terms and other information needed for the search. The first thing to
define would be the essentiality of each component. This is defined under 'presence'
property for each 'name' and 'geneID' defined in the function except for the overall function name. 
The essentiality defined here determines how ProkFunFind reports on the overall results
of the search (i.e., was the function detected or not), but all hits will be reported
for each component weather it is essential or not. 

.. code-bock:: 

    ---
    name: Acetate Production
    components: 
    - name: PTA-ACKA Pathway
      presence: essential
      components: 
      - geneID: PTA
        presence: essential
      - geneID: ACKA
        presence: essential

With the essentility defined, we can then assign the search terms that we want to 
use for each gene. The search terms can be any of the supoorted features for 
ProkFunFind, in this case we are going to use a combination of KOs and COGs. 

    ---
    name: Acetate Production
    components: 
    - name: PTA-ACKA Pathway
      presence: essential
      components: 
      - geneID: PTA
        presence: essential
        terms: 
        - id: K13788
          method: kofam
        - id: K15024
          method: kofam
        - id: COG0280
          method: emapper
        - id: COG0857
          method: emapper
      - geneID: ACKA
        presence: essential
        terms: 
        - id: K00925
          method: kofam
        - id: COG0282
          method: emapper

In this case we can add the KOs,  with the search method being KofamScan and the COGs
with the search method being eggNOG-mapper. 

With the function being defined we can then define the search settings and default filtering thresholds
in the first section of the configuration file (see :ref:`Inputs<inputs>` for more information).

.. code-block::

---
    main:
      cluster_tool: DBSCAN
      faa_suffix: .faa
      gff_suffix: .gff
      fna_suffix: .fna
    kofamscan:
      annot_suffix: .kofam.tsv
      evalue: 1e-3
    emapper:
      annot_suffix: .emapper.annotations
      evalue: 1e-3
    DBSCAN:
      cluster_eps: 4
      cluster_min_samples: 1.8
    ---
    name: Acetate Production
    components: 
    - name: PTA-ACKA Pathway
      presence: essential
      components: 
      - geneID: PTA
        presence: essential
        terms: 
        - id: K13788
          method: kofam
        - id: K15024
          method: kofam
        - id: COG0280
          method: emapper
        - id: COG0857
          method: emapper
      - geneID: ACKA
        presence: essential
        terms: 
        - id: K00925
          method: kofam
        - id: COG0282
          method: emapper

These default settings designate what the annotation files are named and what the
default filtering thresholds are. With the current configuration file,
with default e-value filtering would be applied for each of the terms defined in
the function. These current settings would be a good place to start for most searches
but would likely lead to many spurrious hits with such a high e-value threshold
being used for the filtering. 

Refining search settings
*************************

Like with any search,  the thresholds used to filter hits are going to play a
significant role in the quality and reliability of the search results. 
In the example above a search performed with the default e-value filtering 
at less than or equal to 1e-3 would likely result in many spurrious hits. 
The process of refining a search and adjusting thresholds is usually an 
iterative process, and can vary greatly from gene to gene. 
The following sections will briefly describe some of the approaches that 
can be used to help refine the search and determine appropriate thresholds
to use. 

Comparison to a ground truth
-----------------------------
The easiest way to help refine your search results will be to compare
your results to a ground truth dataset. One way to do this would be to
compare your search to a collection of trait data like what is available
in Madin et al, 2020 (https://doi.org/10.1038/s41597-020-0497-4). If
your trait of interest is included in this or a similar trait database
this could be used to evaluate if your predictions of trait presence 
or absence agree with what is in the database. If you notice that
there are many disagreements between the prediction and ground truth 
then they can be used to evaluate if the hit filtering thresholds need to be
made more strict or more lenient. 

A ground truth comparison could also be done on a smaller scale, for
example looking a the search results for a smaller collection of 
well characterized genomes that you know include your genes of
interest or do not have them. These smaller comparisons will give
you a better idea of how many spurrious hits are being produced,
and the aproximate e-values that you might expect from a 
real hit to one of your genes of interest. This can help 
guide the subsequent refinement of the search settings. 


Evaluating E-value distributions
---------------------------------
When no ground truth is available, you can try to evalutate the quality of
the search by looking at the distirbution of e-values for all the hits
produced in a search. The e-values for searches are reported in the GFF output
from the ProkFunFind search results and can be parsed from those files for 
each genome and each hit. You can then look at the distribution of e-values 
and see if there is an obvious separation between different groups of hits. 
For example you might have a group of hits at an evalue of 1e-3 to 1e-10 and then 
another group of hits at 1e-100 to 1e-110. While you might not have a ground truth 
to compare to, you might be able to make an educated guess that the group of hits 
around 1e-3 to 1e-10 are likely off target hits, and that you could set a new 
e-value threshold at some intermediate value to help filter these out. 


Using a phylogenetic tree to deliniate potential clades
--------------------------------------------------------
Another common approach would be to build a phylogentic tree based on the
gene sequences of your potential hits. You can then evaluate if there are any
clear clades that are formed in the tree. You could then retroactively examine
the hits in each clade to determine if there is a clear difference in 
evalue or identity between the clades that would allow you to set a 
better threshold. 

Comparing multiple annotation approaches
-----------------------------------------
Another apporach would be to use a consensus approach to refine your search. 
This could rely on looking at multiple types of annotations, for example 
orthology definitions and protein domain hits. By checking that your hits 
from an initial search, also include a specific domain that is associated 
with your gene of interest, you might be able to refine your results and
identify a more specific subset of hits to examine. 