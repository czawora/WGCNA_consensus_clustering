
#setup
args<-commandArgs(TRUE)
norm_file_dir <- "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/subsets/"

#gather command line arguments
if (length(args) < 7){

  message("arguments are insufficient")
  quit(save = "no")

} else {
  softPower <- as.numeric(args[1])
  minModuleSize <- as.numeric(args[2])
  merge_eigengenes <- as.numeric(args[3]) #1 = T, else = F
  gene_sub_list_name <- as.character(args[4])
  norm_file <- as.character(args[5])
  save_dir <- as.character(args[6])
  runNum <- as.character(args[7]) ##just to keep track of process success e.g. run12
}

##working directory for output shoul 

filename <- paste(norm_file_dir, norm_file, sep = "")

#fail marker
file.create(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/",save_dir, "/failed/", runNum, sep =""))

setwd(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/",save_dir, sep = ""))

#redirect output
#sink(file = paste(getwd(), "/output/", runNum, ".txt", sep =""))

#print(paste("softPower: ", ,sep = ""))

##herheherhehrehrhehrehrheherh
suppressMessages(library("WGCNA"))
suppressMessages(library("methods"))
suppressMessages(library("edgeR"))
# suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
nthr <- enableWGCNAThreads()
print(nthr)

print("importing data")
selection <- read.csv(filename, header=TRUE)

#set rownames
rownames(selection) <- selection[,"X"]
#remove useless column
selection[,"X"] <- NULL
selection[,"X.1"] <- NULL

print(head(selection))

####
#open gene subsample file and subset names
####

sub_genes <- read.delim(paste("sub_genes/", gene_sub_list_name, sep = ""), sep = "\n", header = FALSE)
print(sub_genes$V1)
selection <- selection[as.character(sub_genes$V1), ]
print(head(selection))

print(dim(selection))
#transpose with only number to keep numeric
transpose <- t(selection)
#selection <- NULL

#"gene35204-v1.0-hybrid"

#remove non variance columsn
print("removing columns with very little variance")
transpose <- transpose[,apply(transpose, 2, var, na.rm=TRUE) > 0.05]
print(dim(transpose))

print("making adjacency matrix")
adjacency_res <- adjacency(transpose, power = softPower, type="signed")

dim(adjacency_res)


distanceMat <- 1 - adjacency_res


print("building gene Tree")
# Call the hierarchical clustering function
geneTree <- hclust(as.dist(distanceMat), method = "average")

# Module identification using dynamic tree cut:
print("cutting modules")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = distanceMat, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize, verbose= 4)

table(dynamicMods)

# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
l <- length(unique(dynamicColors))
l
table(dynamicColors)

if(l > 200){

        file.remove(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/",save_dir, "/failed/", runNum, sep =""))
        print("too many clusters")
        file.create(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/",save_dir, "/2many/", runNum, sep =""))

        quit(save = "no")
}

#write out clusters
gene_names <- rownames(adjacency_res)
clusters <- as.data.frame(cbind(genes= gene_names, group = dynamicColors))
clusters <- clusters[order(clusters$group), ]


if (merge_eigengenes == 0){
  write.csv(clusters , file = paste("/cbcb/project2-scratch/ZCL/", save_dir,"/clusters/", runNum, ".csv", sep =""), row.names = FALSE)
}

moduleColors <- dynamicColors
#if we want to merge cluster based on eigengenes
if (merge_eigengenes == 1){

  # Calculate eigengenes
  MEList <- moduleEigengenes(transpose, colors = dynamicColors)
  MEs <- MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss <- 1-cor(MEs)
  # Cluster module eigengenes
  METree <- hclust(as.dist(MEDiss), method = "average")

  MEDissThres <- 0.25

  # Call an automatic merging function
  merge <- mergeCloseModules(transpose, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  moduleColors <- mergedColors
  print(table(moduleColors))
  print(length(unique(moduleColors)))
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;


  #write out merged clusters
  gene_names <- rownames(adjacency_res)
  mclusters <- as.data.frame(cbind(genes= gene_names, group = mergedColors))
  mclusters <- mclusters[order(mclusters$group), ]
  write.csv(mclusters , file = paste("/cbcb/project2-scratch/ZCL/", save_dir,"/clusters/", runNum, ".csv", sep =""), row.names = FALSE)
}


print("complete")

#stop output redirection
#sink()

file.remove(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/",save_dir,"/failed/",runNum, sep =""))
file.create(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/",save_dir, "/success/",runNum,sep =""))

