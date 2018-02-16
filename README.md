# HiC-spector

A matrix library for spectral and reproducibility analysis of Hi-C contact maps. Several useful functions:

get_reproducibility
  - to calculate the reproducibility metric between 2 HiC contact maps
  
knight_ruiz
  - Knight Ruiz algorithm for matrix balancing
  
get_expect_vs_d_single_chr_v0
  - to find the average contact frequency as a function of genomic distance
  
get_compartment_A_B
  - to find A, B compartment, using method described in Liberman et al. Science 2009
  
and a few functions for binning a genome, and reading HiC maps

<h3>Installation</h3> 
HiC-spector is mostly written in Julia. It has been tested in Julia 0.4 and 0.5. Both the Julia language (http://julialang.org/) and the required packages have to be installed. Please refer to the beginning of the file HiC_spector.jl for the necessary packages. To get some contact maps for testing the code, please follow the instructions shown in the file data/readme_data.

There is a Python script available for quantifying reproducibility. The Python version can read files in genomic coordinates as well as the .hic format (https://github.com/theaidenlab/juicebox/wiki/Data). To do so, please download the Python version of the tool straw (straw.py) developed by the Aiden lab (https://github.com/theaidenlab/straw).

<h3>Usage</h3>
The script run_reproducibility.jl is used to get the reproducibility score from a command-line interface. Usage:
> julia run_reproducibility.jl matrix_file1 matrix_file2 

The input file here is a simple text delimited format with no header.

1 1 20

1 2 18

...

The first and second columns represent the row and column indices of a contact map, whereas the third column is the count. To represent a full matrix, only the upper-triangular component is required. Note that the index should begin with 1. 

Please use the files stored in the folder A549 mentioned in ./data/readme_data to test the script run_reproducibility.jl.

Julia users can include the file HiC_spector.jl for their own analysis by simply using
> include("./HiC_spector.jl");

Please refer to the file hic_spector_tutorial.jl for how to use some of the functions and how to read files in other formats.

For non-Julia users, one can use the Python script run_reproducibility.py to obtain the reproducibility score. Usage:
> python run_reproducibility.py -F matrix_file1 matrix_file2

If the matrix files are labeled in genomic coordinates of bins, USage:
> python run_reproducibility.py t matrix_file1 matrix_file2 40000

where 40000 is the bin size used in the two files

In addition to the text delimited input files, the Python script can calculate reproducibility score for contact maps stored in .hic format. Usage:
> python run_reproducibility.py -f hic_file1 hic_file2 chrid resolution

A script is provided in the tool straw (https://github.com/theaidenlab/straw/tree/master/python) for reading the headers (including chr id and the available resolutions) in .hic file. 

Regarding memory, given two contact maps of human chr1 binned in a bin-size of 10kb, the code works fine in a laptop (16GB memory) from our experience. 

<h3>Aurthor/Support</h3>
Koon-Kiu Yan, koonkiu.yan@gmail.com; Mark Gerstein, mark@gersteinlab.org

<h3>Reference</h3>
Yan KK, Galip Gürkan Yardımcı, William S Noble and Gerstein M. HiC-Spector: a matrix library for spectral and reproducibility analysis of Hi-C contact maps. Bioinformatics 22 March 2017. https://doi.org/10.1093/bioinformatics/btx152
