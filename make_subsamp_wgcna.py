

"""
inputs are: name of experiment directory
			list of gene names
			input count matrix name
			number of samplings
			boolean ("rand" = T, all else = F) describing if random parameters should be used,
			otherwise the default hardcoded parameters will be used

output: experiment directory is created in /cbcb/project-scratch/ZCL/wgcna/
		subdirectories sub_genes/, clusters/, ouput/, failed/, success/, 2many/, and bash/ are created
		file containing the number of times a gene was chosen and the number of subsamplings that were taken

	programs takes a subsample of gene names for the network and some parameters for network creation and produces bash file to run the R script subsamp_wgcna.R in the bash/ subdirectory of the experiment directory
"""


import sys
import os
import random
import subprocess

if "-ls" in sys.argv:

	subprocess.call(["ls", "-d", "*/"])

exp_dir = sys.argv[1]
count_mat = sys.argv[2]
n = int(sys.argv[3])
rand = sys.argv[4]

sub_samp_rate = 0.8

print("making directories")

####make experiment directory
try:

	os.mkdir("/cbcb/project2-scratch/ZCL/" + exp_dir)
	os.mkdir("/cbcb/project2-scratch/ZCL/" + exp_dir + "/clusters")
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir)
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/ouput")
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/failed")
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/success")
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/2many")
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/bash")
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/sub_genes")
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/clusters")
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/config")
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/adjmat")
	os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/indmat")

except:
	print("experiment name already used")
	exit(1)


####make gene combinations
gene_used = {}
gene_list = []
genes = open("/cbcb/lab/smount/ZCL/gene_list.txt")

print("reading gene file")
for g in genes:
	gene_list.append(g.strip("\n"))
	gene_used.update({g.strip("\n"):0})

gene_num = len(gene_list)
select_frac = int(gene_num * 0.8)

powers = [1,2,4,8,12,16]
minModuleSize = [40, 60, 90, 120, 150, 180, 210]
if count_mat[0] == "l":###for using lcm data, to restrict cluster numbers
	minModuleSize = [90,120,150,180,210]
merge_eigengene = [0,1]

cfg = open("/cbcb/project-scratch/ZCL/wgcna/consensus/"+  exp_dir +"/cluster_config.txt", 'w')
cfg.write("power\tminModSize\tmerge\n")
	
print("generating gene combinations")
for x in range(0, n ):

	samp = random.sample(gene_list, select_frac)
	print("sample size: " + str(len(samp)))

	f = open("/cbcb/project-scratch/ZCL/wgcna/consensus/"+exp_dir+ "/sub_genes/" + "sample" + str(x)  + ".txt", 'w')

	for s in samp:

		f.write(s + "\n")
		gene_used[s] += 1

	f.close()

	if rand == "rand":

		random.shuffle(powers)
		random.shuffle(minModuleSize)
		random.shuffle(merge_eigengene)
		p = powers[0]
		minM = minModuleSize[0]
		m = merge_eigengene[0]

	else:

		p = 2
		minM = 90
		merge = 0

	
	cfg.write(str(p) + "\t" + str(minM) + "\t" + str(m) + "\n")

	run = open("/cbcb/project-scratch/ZCL/wgcna/consensus/"+  exp_dir +"/bash/run" + str(x) + ".sh", 'w')

	run.write("#PBS -q throughput\n#PBS -l mem=36GB,walltime=12:00:00,ncpus=2\n")
	run.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/subsamp_wgcna.R " + str(p) + " " + str(minM) + " " + str(m) + " " + "sample" + str(x) + ".txt" +  " " + count_mat + " " + exp_dir + " " + str(x) + "\n")

	run.close()

cfg.close()

print("writing gene sample rates")
gene_file = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/gene_sample_rate.csv", 'w')

gene_file.write("gene\tselected\ttotal_sample\n")

for k in gene_used.keys():

	gene_file.write(k + "\t" + str(gene_used[k]) + "\t" + str(n) + "\n")

gene_file.close()















