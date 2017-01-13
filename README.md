# HiC_spector

A matrix library for spectral analysis and reproducilbilty of Hi-C contact maps. It has been tested in Julia 0.4 and 0.5. Please refer to the beginning of the file HiC_spector.jl for the necessary packages.

Julia users can include the file HiC_spector.jl for their own analysis.


Several useful functions:

get_reproducibility
  - to calculate the reproducibility metric between 2 HiC contact maps

knight_ruiz
  - Knight Ruiz algorithm for matrix balancing

get_expect_vs_d_single_chr_v0
  - to find the average contact frequency as a function of genomic distance

get_compartment_A_B
  - to find A B compartment, using method described in Liberman et al. Science 2009

and a few functions for binning a genome, and reading HiC maps

Please refer to the file hic_spector_tutorial.jl for examples on using some of the functions. To get some contact maps for testing the code, please follow the instructions shown in the file data/readme_data.


Reference: 

Yan KK, Galip Gürkan Yardımcı, William S Noble and Gerstein M. HiC-Spector: a matrix library for spectral and reproducibility analysis of Hi-C contact maps. bioRxiv https://doi.org/10.1101/088922

http://biorxiv.org/content/early/2016/11/21/088922

