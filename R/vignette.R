
if(FALSE){

  
  ### The network fitting is able to use all available CPUs if set up correctly
  library(doParallel)
  totalCores = detectCores()
  cluster <- makeCluster(totalCores) 
  registerDoParallel(cluster)
  library(foreach)
  
  
  ### Experimental directory
  nando_dir <- "/corgi/websites/tcellnet/expnando"

  if(FALSE){
    library(Seurat)
    library(Nando)
    adata <- readRDS(file.path(nando_dir,"allcells.RDS"))
  }
  
  ################################################################################
  ################### Running Pando & first part of Nando data ###################
  ################################################################################
  
  ### First prepare a putative network. Here fitting a very limited set of genes for speed
  DefaultAssay(adata) <- "RNA"
  PrepareNandoDir(nando_dir, adata, use_genome=BSgenome.Hsapiens.UCSC.hg38, include_genes = VariableFeatures(adata)[1:10])
  #write.csv(data.frame(gene=VariableFeatures(adata)[1:10]), file.path(nando_dir,"genes.csv"))
  
  ### Fit all genes using Pando. This can be run in parallel using SLURM to speed it up (see separate vignette)
  ### note!!! Need our modified version of Pando to perform this
  RunPandoParallel(nando_dir)
  
  ### Select clusterings to separate explanations over. Names of clusters must be unique across clusterings
  use_clusterings <- data.frame(
    row.names=colnames(adata),
    
    #RNA based cell type annotation, and all cells together
    dice=adata$pred.dice,
    ALL="ALL",
    
    #Previous annotation, split over donors to assess replicability
    donor_dice=sprintf("%s %s",adata$donor_id, adata$pred.dice),
    donor_ALL=sprintf("%s %s",adata$donor_id, "ALL")
  )
  # use_clusterings$ALL <- "ALL"
  # use_clusterings <- data.frame(
  #   row.names=colnames(adata)
  # )
  # use_clusterings$ALL <- "ALL"
  SelectClusteringForNando(nando_dir, use_clusterings)
  
  ### Collect SHAPs per network. This can be run in parallel using SLURM to speed it up (see separate vignette)
  SummarizeShapByClusterParallel(nando_dir)
  
  
  ################################################################################
  ################### Analyzing Nando data #######################################
  ################################################################################
  
  ### We can now load all SHAP summaries as networks
  #nandonets <- LoadNandoNetworks(nando_dir)
  nandonets <- LoadNandoNetworks(nando_dir)#, keep_clusters = sprintf("donor%s",0:3))
  #hack. if just "donor0" then add "ALL"
  #hack. remove any non-dice for speed
  
  ### For testing, just pick a few
  if(FALSE){
    names(nandonets@nets) <- str_replace_all(names(nandonets@nets),"donordonor","donor")
    nandonets@nets <- nandonets@nets[1:4] #The first donors for the "ALL" network
  }
  
  ### How similar are the networks?
  netsim <- ComputeNandoNetworkSimilarity(nandonets)
  
  ### Compute steady states. This can be run using multiple CPUs if foreach set up
  nandonets <- ComputeSteadyState(nandonets)

  ### Compute hitting probabilities. This can be run using multiple CPUs if foreach set up
  nandonets <- ComputeHittingProbability(nandonets)  #acceptable speed... 20 30 min?
  
  ### It is now possible to check where you would end up if you start from steady state
  hpss <- ComputeHittingProbabilityFromSS(nandonets)
  PlotTopProbabilityMatrix(hpss, min.pmean = 1e-2)  #TODO exclude genes, ss
  #ComputeHittingProbabilityFromSS(nandonets@nets[[1]])
  
  saveRDS(nandonets, "/corgi/websites/tcellnet/finalout.johan/nandonets.rds")
  
  
  ### Plot top genes in steady state
  ss <- SteadyStateMatrix(nandonets)
  PlotTopProbabilityMatrix(ss, min.pmean = 1e-2)
  
  ### Compute a simplified network that is easier to visualize
  net <- nandonets@nets[[1]]
  net@ss <- net@ss[order(net@ss, decreasing = TRUE)]
  keep_genes <- names(net@ss)[1:20]
  tmat <- ComputeSimplifiedMatrix(net, keep_genes)
  plot(TransitionMatrixToIgraph(tmat))
  
  ### Save the network 
  ExportTransitionMatrixCSV(ComputeSimplifiedMatrix(net, keep_genes), "/corgi/websites/tcellnet/graphs/top20simplified.csv")
  ExportTransitionMatrixCSV(TransitionMatrix(net), "/corgi/websites/tcellnet/graphs/fullnetwork.edges.csv")
  write.csv(GeneCategories.NandoNetwork(net),"/corgi/websites/tcellnet/graphs/fullnetwork.genecat.csv", row.names = FALSE)

  ### Comparison of two steady states. If it changes, which edges does the probability flow over?
  flow <- MeltSparsematrix(ComputeSteadyStateChangeEdgeflow(nandonets@nets[[1]],nandonets@nets[[2]]))
  flow <- flow[order(flow$value),]
  head(flow)
  #TOX -> DACH1
  #Are the most variable connections also where most flow happens?
  
  #nandonet <- nandonets@nets[[1]]
  #nandonet@shap
  

  
  ### Store data for ShinyNando to visualize
  
  #TODO importance region
  
  
  ### Perform GO analysis 
  #TODO object for GO analysis? 
  
  
  ### Perturbation analysis of SS  

  
  ### Plotting of the graph
  plot(TransitionMatrixToIgraph(tmat_tf))

  
################################################################################
################### Walktrap analysis with Seurat ##############################
################################################################################

  ### Run a walktrap analysis
  walktraps <- ComputeNandoWalktrap(nandonets, 10)
  NumSteps(walktraps)  
  
  ### Show entropy
  PlotWalktrapEntropy(c("GATA3","C5"))
  
  ### Can also extract entropy of chosen genes and time points
  GetWalktrapEntropy(walktraps, c("GATA3","C5"))
  

    
  
  ### Extract Seurat-like object to compare genes by their walktrap
  library(Seurat)
  gdata <- NandoWalktrapToSeurat(walktraps)
  
  #Split donor & cell type 
  gdata$donor <- str_split_fixed(gdata$cluster," ",2)[,1]
  gdata$ct <- str_split_fixed(gdata$cluster," ",2)[,2]
  
  #Plot: Clustering vs cell type vs donor
  DimPlot(object = gdata, reduction = 'umap', label = TRUE)/
    DimPlot(object = gdata, reduction = 'umap', group.by = "ct", label = TRUE) 
  
  #Plot: Clustering vs cell type vs donor
  DimPlot(object = gdata, reduction = 'umap', group.by = "ct", label = TRUE) |
    DimPlot(object = gdata, reduction = 'umap', group.by = "donor", label = TRUE) 
  
  #Exclude ALL category
  DimPlot(object = gdata[,gdata$category!="ALL"], reduction = 'umap', label = TRUE)/
    DimPlot(object = gdata[,gdata$category!="ALL"], reduction = 'umap', group.by = "category", label = TRUE)
  
  #Plot: Some examples of genes
  DimPlot(object = gdata[,gdata$fromgene %in% unique(gdata$fromgene)[1:20]], reduction = 'umap', group.by = "fromgene", label = TRUE)

    
################################################################################
################### Seurat analysis of hitting probabilities ###################
################################################################################
  
  
  ### Extract Seurat-like object to compare genes by their hitting probabilities
  library(Seurat)
  gdata <- HittingProbabilitiesToSeurat(nandonets)
  
  #Split donor & cell type 
  gdata$donor <- str_split_fixed(gdata$cluster," ",2)[,1]
  gdata$ct <- str_split_fixed(gdata$cluster," ",2)[,2]
  unique(gdata$ct)

  #Plot: Clustering vs cell type vs donor
  DimPlot(object = gdata, reduction = 'umap', label = TRUE)/
    DimPlot(object = gdata, reduction = 'umap', group.by = "ct", label = TRUE) 
  
  #Plot: Clustering vs cell type vs donor
  DimPlot(object = gdata, reduction = 'umap', group.by = "ct", label = TRUE) |
    DimPlot(object = gdata, reduction = 'umap', group.by = "donor", label = TRUE) 

  #Plot: Some examples of genes
  DimPlot(object = gdata[,gdata$fromgene %in% unique(gdata$fromgene)[1:20]], reduction = 'umap', group.by = "fromgene", label = TRUE)
  


################################################################################
################### Shiny Nando ################################################
################################################################################

  #This exports everything such that it can be browsed with Shiny.
  #HP and SS must have been computed for this
  ExportShinyNando(nandonets)


}

