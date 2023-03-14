# HiC_Genome_Assembler
Code repository for a HiC based genome assembly algorithm I wrote during my time in the Baumann lab. Note, that documentation and code may change slightly between now and publication as bugs are encountered and fixed.
https://www.baumannlab.org/links.html<br>

## Algorithm overview

- This is an algorithm designed to assemble a near chromosome scale genome from HiC like data, and an already established primary genome assembly. The algorithm attempts to first group scaffolds to their respective chromosomes and then attempts to order and orient the grouped scaffolds with respect to one another<br><br>

- Our algorithm is dependent on the output of HiCPro, a well established bioinformatics software package that aligns HiC data to a primary genome assembly, removes artifacts, and normalizes the resulting data via a matrix balancing method known as ICED normalization. This particular approach to HiC data processing breaks the up the primary genome assembly into non-overlapping equal sized pieces. The size of these pieces is commonly referred to as the “resolution size”. Read pairs are then counted between any two pieces of the genome, and the resulting square matrix is then normalized to generate the input contact map for our assembly algorithm, which is broken up into two primary phases.<br><br>

- Briely, the first phase of the algorithm attempts to cluster the genomic loci, defined by the contact map, into discrete chromosome groups by contact values alone. In other words, scaffold information is not considered at this point of the algorithm (apart from the initial alignment step), but only contact values between pairs of loci. Moreover, this phase does not require a pre-defined number of chromosomes to cluster within the data. Instead, the number of chromosomes is found dynamically by the algorithms employed. By not incorporating scaffold information, and a pre-defined chromosome count, our process allows for an unbiased initial approach in clustering genomic loci into chromosome sized groups. After this process, scaffold information can then be leveraged to assess the accuracy of the algorithm as well as identify potential errors in the primary genome assembly itself.<br><br>

- The second phase of the algorithm then attempts to order and orient scaffolds within the found chromosome groupings. Amongst all the different ordering and orientation possibilities, we define the best permutation as the maximum value of a function that weights contact values that are close in linear space more heavily than values that are further away. We employ several different algorithms in a sequential fasion in an attempt to find the optimal ordering.


  ![high_levelOverview](https://user-images.githubusercontent.com/20343526/219607825-c99b6578-828b-4031-bc9a-7f730d115262.png)

<br><br>
## General usage and install instructions
Code is entirely python based and can be run by downloading the HIC_ASSEMBLER directory and installing the relevant packages provided in the **packageInstallCommands.txt** file. Additionally, it is recommended to create a new enviroment before running any installs to avoid complications.
Provided is a working example config file (hicAssembler_config_workingExample.txt), which is intended to provide the user with a template, but should be altered to fit ones needs. Once the config file has been updated with the relevant settings / parameters, one can run the following command(s) to run the program<br>

**To run the entire pipeline in one fell swoop, then run with the following command**
- python HIC_ASSEMBLER/run_hicAssembler.py -part1 -part2 -part3 -part4 -c file/path/to/hicAssembler_config.txt

**Note that -part3 is optional and depends on a validpair file produced by HiCpro as well as data produced by a rescriction enzyme (if your data is of newer hic like data types then this might not be applicable). In this case, set the finalOrderingsFile option in the config file to be the same as the chromosomeOrderFile**
- python HIC_ASSEMBLER/run_hicAssembler.py -part1 -part2 -part4 -c file/path/to/hicAssembler_config.txt

**Note that any part can be run independently, as long as the output from the previous parts has been generated. For example, to just run part2 then run the following command**
- python HIC_ASSEMBLER/run_hicAssembler.py -part2 -c file/path/to/hicAssembler_config.txt

## Part1

- Clusters contact map rows together via average clustering, and identifies matrix cut indices via hyper geometric cutting algorithm or HMM. Smaller groups are optionally resolved by modularity maximazation. Initial grouping assessment is also performed
  
  ![Screen Shot 2023-03-10 at 1 45 07 PM](https://user-images.githubusercontent.com/20343526/224424933-167ccb6d-a552-407f-8565-8aac53f63dad.png)


## Part2

- Orders and orients scaffolds with respect to on another on any given chromsome

  ![part2](https://user-images.githubusercontent.com/20343526/219610194-72c2aef0-8cde-45f9-8ecb-a73cba6ef967.png)



## Part3

- Attempts to correct potentially ambigously oriented scaffolds because of a resolution size that is too large. This step is optional and requires a restriction enzyme cut site file (that is also used by HiC-Pro) along with the valid pairs file produced by HiC-Pro that details read pair mappings

  ![part3](https://user-images.githubusercontent.com/20343526/219610278-6c374383-3f68-494b-816b-19f549fc4442.png)



## Part4

- Writes out the final assembled genome with new gaps introduced as a sequence of 100 "N" characters. Scaffolds that were unable to be grouped, are still written to this new version with their original sequence name intact.

  ![part4](https://user-images.githubusercontent.com/20343526/219610312-04df94eb-9a08-455b-b2fc-96bcbcb6628f.png)

