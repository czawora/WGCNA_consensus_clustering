
making a consensus network

---------------------------------------------------------------------------------
---------------------------------------------------------------------------------
---------------------------------------------------------------------------------

packages needed:

WGCNA
methods
edgeR - maybe??
RColorBrewer
data.table
pryr

---------------------------------------------------------------------------------
---------------------------------------------------------------------------------
---------------------------------------------------------------------------------

files included:

*********
##BEWARE hardcoded ibis server filepaths are present in these files, you should change them 
*********

make_subsamp_wgcna.py
subsamp_wgcna.R
seq_inidicator_mat.R
seq_adding_mat.R
merge_mats.R
consensus_cluster.R


Step 1 -- make subsample network bash files 
use: make_subsamp_wgcna.py 
	 subsamp_wgcna.R


###########inputs are: name of experiment directory
###########			list of gene names
###########			input count matrix name
###########			number of samplings
###########			boolean ("rand" = T, all else = F) describing if random parameters ###########         should be used,
###########         otherwise the default hardcoded parameters will be used

###########output: experiment directory is created in /cbcb/project-scratch/ZCL/wgcna/
###########		subdirectories sub_genes/, clusters/, ouput/, failed/, success/, 2many/, ###########     and bash/ are created
###########		file containing the number of times a gene was chosen and the number of ###########     subsamplings that were taken

###########	programs takes a subsample of gene names for the network and some parameters ########### for network creation and produces bash file to run the R script subsamp_wgcna.########### R in the bash/ subdirectory of the experiment directory

then run the bash files

Step 2 make indicator matrix (can be done while networks are running)
use: seq_inidicator_mat.R ## need to change hardcoded things

###########inputs:
###########exp_name - name of "experiment" created by make_subsamp_wgcna.py (need to ###########change hardcoded paths)
###########start and stop -- start and stop numbers of gene sample lists to create a ###########indicator matrix for. This allows for parallel running of seq_inidicator_matrix.###########R for a speed up

###########output: 
###########a gene-by-gene matrix counting the number of times each pair of genes was ###########sampled together for input

Step 3 make co-clustering matrix
use: seq_adding_mat.R ## need to change hardcoded things

###########same inputs and outputs as seq_inidicator_mat.R but now the output matrix is ###########the gene-by-gene matrix counting the number of times pairs clustered together 

Step 4 combine the matrices made in parallel into one indicator matrix and one co-clustering matrix
use: merge_mats.R

Step 5 divide the co-clustering matrix by the indicator matrix to make a consensus matrix that has values [0,1]. This can be used as a distance value now.

Step 6 use the 1 minus the consensus matrix for new distance unit to do WGCNA clustering
use consensus_cluster.R


