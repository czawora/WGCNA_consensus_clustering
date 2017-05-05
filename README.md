# WGCNA_consensus_clustering

A collection of files I used to implement consensus clustering on F. vesca gene expression data. Most files have hardcoded paths in them and need to be changed before any other use.

These files were written to take advantage of parallel processing on my local server and therefore are designed with many intermediary steps and files being generated along the way. WGCNA can be used to simply cluster genes, however I have also coupled it with an interative resampling method to try to establish consensus clusters. 

The steps I followed are borrowed from this paper 
http://dx.doi.org/10.1023/A:1023949509487
and using WGCNA for clustering

Of course you can also just cluster using WGCNA without going through the re-sampling method I used.

For details on each of the files, check the other README.txt file
