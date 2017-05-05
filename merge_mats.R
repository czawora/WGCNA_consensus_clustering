
args <- commandArgs(TRUE)

library(data.table)

mat_dir <- as.character(args[1])
first_mat_name <- as.character(args[2])


master <- fread(paste(mat_dir, first_mat_name, sep = ""), sep = ",", header = TRUE)

for(i in seq(3,length(args))){

	m <- as.character(args[i])

	next_mat <- fread(paste(mat_dir, m, sep = ""), sep = ",", header = TRUE)

	master <- next_mat + master

}

fwrite(master, file = paste(mat_dir, "master.csv", sep = ""))
