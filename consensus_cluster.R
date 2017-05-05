
library(WGCNA)
library(data.table)

print("loading libraries")

#setup
args<-commandArgs(TRUE)

#gather command line arguments
if (length(args) < 3){

  message("arguments are insufficient")
  quit(save = "no")

} else {
  minModuleSize <- as.numeric(args[1])
  consensusMat_fname <- as.character(args[2])#includes full path
  save_dir <- as.character(args[3])
}

zz <- file(paste0(save_dir,"/log.log"), open="wt")
sink(zz, type = "message")

file.create(paste(save_dir,"/failed", sep =""))

setwd(paste(save_dir,"/..", sep = ""))


##working directory for output should

adjacency_res <- as.matrix(as.data.frame(fread(consensusMat_fname, header = TRUE)))
rownames(adjacency_res) <- colnames(adjacency_res)
dist_mat <- 1 - adjacency_res


#######loop through different module sizes and check on average consensus of cluster

print("building gene Tree")
# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dist_mat), method = "average")

# Module identification using dynamic tree cut:
print("cutting modules")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dist_mat, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize, verbose= 4)

table(dynamicMods)


write.csv(as.character(dynamicMods), file= paste0(save_dir,"/dynamicMods.csv"))
# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
l <- length(unique(dynamicColors))
l
table(dynamicColors)

if(l > 200){

        file.remove(paste(save_dir, "/failed", sep =""))
        print("too many clusters")
        file.create(paste(save_dir, "/2many",  sep =""))

        quit(save = "no")
}

#write out clusters

gene_color <- data.frame(gene = rownames(dist_mat), cluster = dynamicColors)

gene_color <- gene_color[which(gene_color$cluster != "grey"),]

gene_color <- gene_color[order(gene_color$cluster),]

gene_num <- data.frame(gene = gene_color$gene, cluster = as.numeric(gene_color$cluster))

gene_color$cluster <- as.character(gene_color$cluster)

write.csv(gene_color, file = paste0(save_dir, "/gene_cluster_color.csv"), row.names = FALSE)
write.csv(gene_num, file = paste0(save_dir, "/gene_cluster_num.csv"), row.names = FALSE)


file.remove(paste(save_dir,"/failed", sep =""))
file.create(paste(save_dir, "/success",sep =""))




