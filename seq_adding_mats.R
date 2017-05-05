

#####
#####This script is for creating 1,0 indicator matrices for the clusters produced by subsamp_wgcna.R
####


library("data.table")
library(pryr)

args <- commandArgs(TRUE)

exp_name <- as.character(args[1])
start <- as.numeric(args[2])
stop <- as.numeric(args[3])

##create matrix from gene list
gene_f <- read.delim("/cbcb/lab/smount/ZCL/gene_list.txt", sep = "\t", header = FALSE)

genes <- as.character(gene_f$V1)


#############################################################
###construct master matrix as data.table
#############################################################

#make matrix
print("make matrix")
#con_mat <- matrix(0, nrow = length(genes), ncol = length(genes))
master_mat <- matrix(0, nrow = length(genes), ncol = length(genes))


#convert to data frame
print("make data frame")
master_df <- as.data.frame(master_mat)
rm(master_mat)

#set colnames
colnames(master_df) <- genes
master_df$g <- genes

print("make data table")
#convert to data table
master_dt <- data.table(master_df)
rm(master_df)

setkey(master_dt, g)

master_dt[, g := NULL]

#############################################################
###construct matrix as data.table
#############################################################

#make matrix
print("make matrix")
#con_mat <- matrix(0, nrow = length(genes), ncol = length(genes))
con_mat <- matrix(0, nrow = length(genes), ncol = length(genes))
print(dim(con_mat))


#convert to data frame
print("make data frame")
con_df <- as.data.frame(con_mat)
rm(con_mat)

#set colnames
colnames(con_df) <- genes
#con_df$g <- genes

print("make data table")
#convert to data table
con_dt <- data.table(con_df)
rm(con_df)

#setkey(con_dt, g)

#############################################################
#############################################################
#############################################################


for(runNum in seq(start,stop)){


	con_dt[,] = 0
	con_dt[, g := genes]
	setkey(con_dt, g)


	clust_dir <- "/cbcb/project2-scratch/ZCL/"


	#open cluster list
	clust_name <- paste(clust_dir, exp_name, "/clusters/",  runNum, ".csv", sep = "")

	clust <- read.csv(clust_name)

	for(m in unique(clust$group)){

		print(mem_used())
		if(m == "grey"){
			next
		}

		mod_genes <- as.character(clust[which(clust$group == m), "genes"])

		print(length(mod_genes))

		con_dt[mod_genes[1:length(mod_genes)], mod_genes[1:length(mod_genes)] := 1]

		# combos <- expand.grid(mod1 = mod_genes, mod2 = mod_genes)

	}

	con_dt[, g := NULL]

	master_dt <- master_dt + con_dt

}
# print(sum(con_dt[, lapply(.SD, sum, na.rm=TRUE), .SDcols=genes]))

fwrite(master_dt, file = paste( "/cbcb/project-scratch/ZCL/wgcna/consensus/", exp_name,"/adjmat/", start, "_", stop, ".csv" ,sep = ""))



