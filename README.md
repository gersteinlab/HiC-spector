# HiC-spector

A matrix library for spectral and reproducilbilty analysis of Hi-C contact maps. Several useful functions:

get_reproducibility
  - to calculate the reproducibility metric between 2 HiC contact maps
  
knight_ruiz
  - Knight Ruiz algorithm for matrix balancing
  
get_expect_vs_d_single_chr_v0
  - to find the average contact frequency as a function of genomic distance
  
get_compartment_A_B
  - to find A B compartment, using method described in Liberman et al. Science 2009
  
and a few functions for binning a genome, and reading HiC maps

<h3>Installation</h3> 
HiC-spector is mostly written in Julia. It has been tested in Julia 0.4 and 0.5. Both the Julia language and several packages have to be installed. Please refer to the beginning of the file HiC_spector.jl for the necessary packages. To get some contact maps for testing the code, please follow the instructions shown in the file data/readme_data.

There is a python script available for quantifying reproducibility. The python version is able to read files in .hic format. To do so, please download the python version of the tool straw (straw.py) developed by the Aiden lab (https://github.com/theaidenlab/straw).

<h3>Usage</h3>
Julia users can include the file HiC_spector.jl for their own analysis by simply using
> include("./HiC_spector.jl");

Please refer to the file hic_spector_tutorial.jl for examples on using some of the functions. 
The script run_reproducibility.jl is used to get the reproducibility score from a command-line interface. Usage:
> julia run_reproducibility.jl matrix_file1 matrix_file2 

Please refer to ./data/readme_data for the format of the input files

For non-julia users, one can use the python script run_reproducibility.py to obtain the reproducibility score. Usage:
> python run_reproducibility.py matrix_file1 matrix_file2

In addition to the text delimited input files, the python script can calculate reproducibility score by reading files in .hic format. Usage:
>

<h3>Aurthor/Support</h3>
Koon-Kiu Yan, koon-kiu.yan@yale.edu; Chengfei Yan, chengfei.yan@yale.edu

<h3>Reference</h3>
Yan KK, Galip Gürkan Yardımcı, William S Noble and Gerstein M. HiC-Spector: a matrix library for spectral and reproducibility analysis of Hi-C contact maps. bioRxiv https://doi.org/10.1101/088922
