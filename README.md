# HiC_Genome_Assembler
Documentation is still in progres...
Code repository for a HiC based genome assembly algorithm I wrote during my time in the Baumann lab 
https://www.baumannlab.org/links.html<br>

## Algorithm overview

- This is an algorithm designed to assemble a near chromosome scale genome from HiC like data, and an already established primary genome assembly. The algorithm attempts to first group scaffolds to their respective chromosomes and then attempts to order and orient the grouped scaffolds with respect to one another<br><br>

- Our algorithm is dependent on the output of HiCPro, a well established bioinformatics software package that aligns HiC data to a primary genome assembly, removes artifacts, and normalizes the resulting data via a matrix balancing method known as ICED normalization. This particular approach to HiC data processing breaks the up the primary genome assembly into non-overlapping equal sized pieces. The size of these pieces is commonly referred to as the “resolution size”. Read pairs are then counted between any two pieces of the genome, and the resulting square matrix is then normalized to generate the input contact map for our assembly algorithm, which is broken up into two primary phases.<br><br>

- Briely, the first phase of the algorithm attempts to cluster the genomic loci, defined by the contact map, into discrete chromosome groups by contact values alone. In other words, scaffold information is not considered at this point of the algorithm (apart from the initial alignment step), but only contact values between pairs of loci. Moreover, this phase does not require a pre-defined number of chromosomes to cluster within the data. Instead, the number of chromosomes is found dynamically by the algorithms employed. By not incorporating scaffold information, and a pre-defined chromosome count, our process allows for an unbiased initial approach in clustering genomic loci into chromosome sized groups. After this process, scaffold information can then be leveraged to assess the accuracy of the algorithm as well as identify potential errors in the primary genome assembly itself.<br><br>

- The second phase of the algorithm then attempts to order and orient scaffolds within the found chromosome groupings. Amongst all the different ordering and orientation possibilities, we define the best permutation as the maximum value of a function that weights contact values that are close in linear space more heavily than values that are further away. We employ several different algorithms in a sequential fasion in an attempt to find the optimal ordering.


## Part1

Clusters contact map rows together via average clustering, and identifies matrix cut indices via HMM and modularity maximazation steps

## Part2

Orders and orients scaffolds with respect to on another on any given chromsome


## Part3

Attempts to correct potentially ambigously oriented scaffolds because of a resolution size that is too large. This step is optional and requires a restriction enzyme cut site file (that is also used by HiC-Pro) along with the valid pairs file produced by HiC-Pro that details read pair mappings

## Part4

Writes out the final assembled genome with new gaps introduced as a sequence of 100 "N" characters. Scaffolds that were unable to be grouped, are still written to this new version with their original sequence name intact.
