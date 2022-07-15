.. ProkFunFind

.. _outputs:


**************
Output dataset
**************

======================
Output format overview
======================

  ProkFunFind will output files with the same prefix. ``annot.gff``, ``json``, ``tsv``, ``pkl`` formatted files will be outputed every time.
  Other formatted file will be reported based on detection method.


blast.m6
========
  The blast output: query is protein sequences of genomes; database is the bait sequences; format is ``-m 6``

  ::
  
    $ cat output.prefix.blast.m6
    > GUT_GENOME143137_00182	ecoli_garR	52.203	295	134	3	2	292	3	294	3.80e-104	297
    > GUT_GENOME143137_00182	cclostridioforme_GarR	51.890	291	139	1	2	291	3	293	6.49e-104	296
    > GUT_GENOME143137_00419	ecoli_garD	30.208	384	244	7	5	368	122	501	2.39e-52	173
    > GUT_GENOME143137_00419	cclostridioforme_GarD	33.113	302	195	4	68	368	188	483	5.51e-51	169
    > GUT_GENOME143137_00648	cclostridioforme_gudA	28.041	296	205	3	82	372	46	338	1.34e-36	126
    > GUT_GENOME143137_00649	cclostridioforme_gudC	27.273	154	98	4	4	145	8	159	3.05e-06	35.0
    > GUT_GENOME143137_00650	cclostridioforme_gudB	41.837	490	277	4	13	497	17	503	3.93e-117	345
    > GUT_GENOME143137_00901	cclostridioforme_gudA	25.200	250	143	6	72	301	46	271	1.78e-12	57.4
    > GUT_GENOME143137_00903	cclostridioforme_gudB	34.647	482	311	3	17	496	19	498	2.31e-92	281
    > GUT_GENOME143137_00918	cclostridioforme_GarL	22.378	286	216	4	3	284	6	289	4.03e-21	81.3
  
  
annot.gff
=========
  The annotated gff file: genes matched with bait sequences and passed filter criteria were annotated with the first blast hit.
  This file can be imported to Geneious to manually curate.
  
  ::
  
    $ cat output.prefix.annot.gff
    > GUT_GENOME143137_1	GuFunFind	CDS	187675	188575	.	-	.	ID=GUT_GENOME143137_00182;Name=garR;Parent=Cl_0;Target=ecoli_garR 2 294;pct_identity=52.203;evalue=3.8e-104
    > GUT_GENOME143137_2	GuFunFind	CDS	38455	39622	.	+	.	ID=GUT_GENOME143137_00419;Name=garD;Parent=Cl_0;Target=ecoli_garD 121 501;pct_identity=30.208;evalue=2.39e-52
    > GUT_GENOME143137_2	GuFunFind	CDS	296187	297321	.	+	.	ID=GUT_GENOME143137_00648;Name=gudA;Parent=Cl_1;Target=cclostridioforme_gudA 45 338;pct_identity=28.041;evalue=1.34e-36
    > GUT_GENOME143137_2	GuFunFind	CDS	297794	299288	.	+	.	ID=GUT_GENOME143137_00650;Name=gudB;Parent=Cl_1;Target=cclostridioforme_gudB 16 503;pct_identity=41.837;evalue=3.93e-117
    > GUT_GENOME143137_3	GuFunFind	CDS	239793	240903	.	+	.	ID=GUT_GENOME143137_00901;Name=gudA;Parent=Cl_0;Target=cclostridioforme_gudA 45 271;pct_identity=25.2;evalue=1.78e-12
    > GUT_GENOME143137_3	GuFunFind	CDS	241398	242916	.	+	.	ID=GUT_GENOME143137_00903;Name=gudB;Parent=Cl_0;Target=cclostridioforme_gudB 18 498;pct_identity=34.647;evalue=2.31e-92
    > GUT_GENOME143137_3	GuFunFind	CDS	262245	263112	.	+	.	ID=GUT_GENOME143137_00918;Name=garL;Parent=Cl_1;Target=cclostridioforme_GarL 5 289;pct_identity=22.378;evalue=4.03e-21
    > GUT_GENOME143137_4	GuFunFind	CDS	87018	88551	.	-	.	ID=GUT_GENOME143137_01073;Name=gudB;Parent=Cl_0;Target=cclostridioforme_gudB 18 507;pct_identity=36.531;evalue=6.07e-95
    > GUT_GENOME143137_4	GuFunFind	CDS	89063	90152	.	-	.	ID=GUT_GENOME143137_01075;Name=gudA;Parent=Cl_0;Target=cclostridioforme_gudA 47 301;pct_identity=26.515;evalue=1.06e-15
    > GUT_GENOME143137_5	GuFunFind	CDS	36139	37639	.	-	.	ID=GUT_GENOME143137_01304;Name=gudB;Parent=Cl_0;Target=cclostridioforme_gudB 3 480;pct_identity=40.167;evalue=3.42e-120
  
  
json
====
  The json file is similar to the input ``system.json`` file with the "genes" added and the "completeness" of each subsystem.
  
  ::
  
    $ cat output.prefix.json
    > {
    >     "name": "Mucic_and_Saccharic_Acid",
    >     "components": [
    >         {
    >             "orthoID": "garD",
    >             "presence": "essential",
    >             "genes": [
    >                 "GUT_GENOME143137_00419",
    >                 "GUT_GENOME143137_02123",
    >                 "GUT_GENOME143137_02124",
  
  

tsv
====
  A tab separated file with three columns.
  
  ::
  
    $ cat output.prefix.tsv
    > Gene_Name	Cluster_ID	Functions
    > GUT_GENOME143137_00182	GUT_GENOME143137_1:Cl_0	Mucic_and_Saccharic_Acid/garR
    > GUT_GENOME143137_00419	GUT_GENOME143137_2:Cl_0	Mucic_and_Saccharic_Acid/garD
    > GUT_GENOME143137_00648	GUT_GENOME143137_2:Cl_1	Mucic_and_Saccharic_Acid/gudA
    > GUT_GENOME143137_00650	GUT_GENOME143137_2:Cl_1	Mucic_and_Saccharic_Acid/gudB
    > GUT_GENOME143137_00901	GUT_GENOME143137_3:Cl_0	Mucic_and_Saccharic_Acid/gudA
    > GUT_GENOME143137_00903	GUT_GENOME143137_3:Cl_0	Mucic_and_Saccharic_Acid/gudB
    > GUT_GENOME143137_00918	GUT_GENOME143137_3:Cl_1	Mucic_and_Saccharic_Acid/garL
    > GUT_GENOME143137_01073	GUT_GENOME143137_4:Cl_0	Mucic_and_Saccharic_Acid/gudB
    > GUT_GENOME143137_01075	GUT_GENOME143137_4:Cl_0	Mucic_and_Saccharic_Acid/gudA

  
  ==============   ==================================================================
  Name             Description
  ==============   ==================================================================
  Gene_Name        The name of gene
  --------------   ------------------------------------------------------------------
  Cluster_ID       Genes with the same Cluster_ID are within the same neighborhood
  --------------   ------------------------------------------------------------------
  Functions        The function in the system  
  ==============   ==================================================================

pkl
====
  pickle_ object of the output Genome object, which can be loaded to python for further analysis.
  
  .. _pickle: https://docs.python.org/3/library/pickle.html


==========================
Outputs in different mode
==========================

Single genome mode(``rep``)
===========================
  In this mode, it will output four(Interproscan) or five(blast) files with the same prefix.

Pangenome mode(``pan``)
=======================
  In this mode, it will output the results for all genomes belongs to the same pangenome.
  All the results will be in the same folder and outputs for the same genome will have the same prefix.

================
What to do next?
================

  1. Import the ``prefix.annot.gff`` to Geneious software and start the manually curation process
  2. Based on the curated result to change different parameters in ``detect.ini``, ``filter.ini`` and ``cluster.ini``
  3. Re-run ``prokfunfind.py`` to test if the change parameters can improve annotation
