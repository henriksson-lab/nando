if(FALSE){
  install.packages("shiny")
}

library(shiny)

#Load packages early to ensure they exist
library(RSQLite)
library(Signac)
library(Seurat)
library(hdf5r)
library(stringr)
library(ggplot2)
library(patchwork)
library(sqldf)
library(reshape2)
library(cowplot)
library(ggforce) #for geom_circle
library(Matrix)


print("======= load nando ================ ")

source("../R/generics.R")
source("../R/includes.R")
source("../R/objects.R")
source("../R/plots.R")
source("../R/preprocessing.R")
source("../R/shinyexport.R")
source("../R/utils.R")




print("======= reading data to be cached in memory ================ ")

#Location to data files
nando_dir <- "/corgi/websites/tcellnet/expnando"

#TODO note, format changed
genecat <- ImportGeneCategories(nando_dir)  

shapregs_genelist <- ImportShapregsGenelist(nando_dir)

print("read adata") ############# comment out below for speedy testing
adata <- readRDS(file=file.path(nando_dir,"allcells.RDS"))  ### would be ideal if we did not have to do this. used for what?

available_clusterings <- colnames(ReadNandoClustering(nando_dir))
gene_categories <- ImportGeneCategories(nando_dir)
available_tf_in_ss <- rownames(gene_categories)[gene_categories[,1]=="IrreducibleTF"]

#Figure out which clusters in each clustering
available_clusters <- list()
nando_clustering <- ReadNandoClustering(nando_dir)
for(i in 1:ncol(nando_clustering)){
  available_clusters[[colnames(nando_clustering)[i]]] <- unique(nando_clustering[,i])
}

#Prepare h5 files for reading
h5hp <- PrepareProbabilityMatrixH5(file.path(nando_dir,"hp.h5"))
#h5ss <- PrepareProbabilityMatrixH5(file.path(nando_dir,"ss.h5"))



#website_dir <- "./"   #"/corgi/websites/tcellnet/shiny/"

print("========== global done ================")


