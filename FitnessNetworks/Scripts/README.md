## Fitness Network Scripts
This directory contains a number of scripts containing functions to generate/analyze genotypic fitness networks from SNP data encoded in 0/1/2 format, and for visualization. 

### NetworkSummaryFunctions-V2.R
This script contains a set of functions used to read in genotype data, resample and process into network format, and summarize aspects of the resultant networks shape, extract genotypic trajectories separating species (or any focal nodes), and summarize accessibility. 

### Access-A2M-A2P-AllSamp-V2.R 
This script uses the functions defined in the script above to sample genotypes, estimate networks, and summarize the shape of these networks, and the accessibility of interspecific genotypic trajectories. This version of the script used data for all hybrids (from both Martin & Wainwright 2013, as well as Martin & Gould 2020).

### Access-A2M-A2P-RTM2-V2.R
The same as the above, but only using hybrids from the second fitness experiment (Martin & Gould 2020).

### Summarize-SppDists-10snp-Nets.R
This script simply estimates networks from 10 SNPs sampled from all adaptive loci (irrespective of which species it swept in) to summarize overall network structure, irrespective fo fitness. 
