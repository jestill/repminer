RepMiner
========

Description
-----------

The RepMiner package takes a graph theory approach to the classification and assembly of the repetitive fraction of genomic sequence data. Sequences analyzed by RepMiner can range from full length transposable elements in well characterized genomes to short length sequence reads resulting from low coverage sample sequence data. RepMiner makes use of transposable elements identified from model species to map the location of putative transposable elements onto homology based networks derived from comparing the sequences of the query genome to itself. Individual clusters representing Pseudo Assembly Networks (PANs) may be selected and assembled using the TGICL/CAP3 program.


Features
--------

Fully implemented features:

Automated All-by-All BLAST given a FASTA file for the query sequences of interest
Automated BLAST against a set of databases of transposable elements (TEs) from model species
Parsing of TE BLAST results into a classified set of transposable elements
Production of networks of similarity based on BLAST results that are suitable for visualization in the program Cytoscape
Automated assembly of manually selected Pseudo Assembly Networks using TGICL/CAP3
Automated comparison of transposable element assembly against a database of profile hidden Markov Models using HMMER
Parsing of all HMMER results into an easy to read tab delimited text file
Creation of an Apollo file for all assemblies resulting from TGICL/CAP3
Production of an HTML page summarizing all analyses for each of the PANS submitted to the program
