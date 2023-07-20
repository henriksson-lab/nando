library(markovchain)
library(Matrix)
library(stringr)
library(doParallel)
library(foreach)

# Use this to update the documentation
if(FALSE){
  setwd("/corgi/henriksson/jupyter/tcellpaper/Nando")
  devtools::document()
}

################################################################################
################### Running of Pando ###########################################
################################################################################

#' Prepare seurat object for analysis with Pando, and then Nando
#' 
#' @param adata Seurat object
#' @param nando_dir Location of all Nando files. Created if not existing
#' @param include_genes List of genes to fit; using highly variable genes if not set
#' @param use_genome The genome to use for the putative network
#' @return Nothing, creates directory and initial files
#' 
#' @export
PrepareNandoDir <- function(nando_dir, adata, include_genes=NULL, use_genome=BSgenome.Hsapiens.UCSC.hg38){

  if(!file.exists(nando_dir)){
    dir.create(nando_dir)
  }
  
  if(is.null(include_genes)){
    DefaultAssay(adata) <- "RNA"
    include_genes <- VariableFeatures(adata)
  }
  
  if(FALSE){
    variable_adata <- adata
    DefaultAssay(variable_adata) <- "RNA"
    variable_adata <- FindVariableFeatures(variable_adata, verbose = T, nfeatures=nfeatures)
    include_genes <- VariableFeatures(variable_adata)
  }
  
  ## Produce putative network. Note that we use the patched version
  adata <- fixed.initiate_grn.Seurat(adata, peak_assay="ATAC")
  adata <- find_motifs( # must run after initiate_grn
    adata,
    pfm = motifs,
    genome = use_genome
  )

  ### Save putative network + atac + rna data together for network fitting
  saveRDS(adata, file=file.path(nando_dir,"allcells.RDS"))

  ### Save a list of which genes to fit  
  write.csv(data.frame(gene=include_genes), file.path(nando_dir,"genes.csv"))
}




#' Run Pando to fit genes and store SHAP scores per cell
#' 
#' This function should be used in an R script run in a SLURM batch array.
#' It will use SLURM_ARRAY_TASK_COUNT and SLURM_ARRAY_TASK_ID to divide the
#' summarization between different nodes
#' 
#' @param nando_dir Location of all Nando files
#' @return Nothing, saves to disk immediately
#' 
#' @export
RunPandoParallel <- function(nando_dir){
  
  #Read the data + putative network
  adata <- readRDS(file=file.path(nando_dir,"allcells.RDS"))
  
  #Which genes to fit?
  list_genes <- read.csv(file.path(nando_dir,"genes.csv"))$gene
  
  #Figure out how to divide it all for SLURM
  num_tasks <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
  task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  if(is.na(num_tasks)){
    print("SLURM jobarray not detected; running in a single thread")
    num_tasks <- 1  
    task_id <- 0
  } else {
    print(paste("Divide slurm",num_tasks,task_id))
  }
  
  ## Figure out which files this particular job should do
  list_taskid <- (1:length(list_genes))%%num_tasks
  list_genes_to_process <- list_genes[list_taskid==task_id]
  print("Will process genes =====================")
  print(list_genes_to_process)
  print("========================================")
  
  ## Do the actual fitting
  adata_pando <- Nando::infer_grn(
    adata,
    genes = list_genes_to_process,
    parallel = FALSE,
    peak_to_gene_method = 'Signac', ############### todo new argument to tell where to store SHAPs
    verbose = 2,
    method = 'xgb',
    tf_cor = 0
  )
  
  print("Finished fitting genes")
}


################################################################################
################### Summarizing SHAP scores ####################################
################################################################################

#' Selection of clusterings over which SHAPs will be summarized
#' 
#' @param use_clusterings Data frame, each column with a clustering. Row names must be all the cell names from the Seurat object
#' 
#' @export
SelectClusteringForNando <- function(nando_dir, use_clusterings){
  write.csv(use_clusterings, file.path(nando_dir,"clusterings.csv"))
}


#' Helper function
#' 
SummarizeShapByCluster_summarizeClusterReg <- function(shap, for_category, onecat, curgene){
  #Summarize values for each group of cells
  foreach(onecat = unique(for_category), .combine = "rbind",  .verbose = F) %do% {
    shapcat <- colMeans(abs(shap[for_category==onecat,,drop=FALSE]))
    
    outdf <- data.frame(
      gene=curgene,
      category=onecat,
      #clustering=cur_clustering,
      region=names(shapcat),
      gain=shapcat
    )
    
    #This saves memory
    outdf <- outdf[outdf$gain!=0,,drop=FALSE]
    
    #Extract target gene
    pretarget <- str_split_fixed(outdf$region,":",2)
    outdf$tf     <- apply(pretarget,1,function(x) x[!str_starts(x,"chr")])
    outdf$region <- apply(pretarget,1,function(x) x[str_starts(x,"chr")])
    
    #Sum gain
    sqldf::sqldf("select category, region, tf, gene, sum(gain) as gain from outdf group by tf,gene,category,region")
  }
}

#' Helper function
SummarizeShapByCluster_summarizeClusterNoreg <- function(shap, for_category, onecat, curgene){
  foreach(onecat = unique(for_category), .combine = "rbind",  .verbose = F) %do% {
    shapcat <- colMeans(abs(shap[for_category==onecat,,drop=FALSE])) 
    
    outdf <- data.frame(
      gene=curgene,
      #clustering=cur_clustering,
      category=onecat,
      region=names(shapcat),
      gain=shapcat
    )
    
    #This saves memory
    outdf <- outdf[outdf$gain!=0,,drop=FALSE]  
    
    #Extract target gene
    pretarget <- str_split_fixed(outdf$region,":",2)
    outdf$tf <- apply(pretarget,1,function(x) x[!str_starts(x,"chr")])   #actually, this is the TF!
    
    #Sum gain over regions. Will save memory and speed up things later
    outdf <- sqldf::sqldf("select category, tf, gene, sum(gain) as gain from outdf group by tf,gene,category")
    outdf
  }
}



#' Compute summarized SHAP values for a number of clusters.
#' Note that cells can be part of more than one cluster.
#' 
#' This function is not meant to be used directly; most users should rather use a wrapper
#' that loops over all files
#' 
#' @param shap_dir Directory with SHAP per cell
#' @param list_files_to_process List of files to process
#' @param taskid Process ID when multiprocessing, or 0
#' @param summarize_shap_for_cluster Clusterings of cells
#' 
#' @export
SummarizeShapByCluster <- function(nando_dir, shap_dir, list_files_to_process, summarize_shap_for_cluster, taskid){
  
  #Create place to store summaries
  if(!file.exists(shap_dir)){
    dir.create(shap_dir)
  }
  
  #Loop over all files
  all_shap <- list()
  all_shap_reg <- list()
  for(f in list_files_to_process){
    
    curgene <- str_remove(str_remove(f,"shap ")," .RDS")
    print(curgene)
    
    shap <- readRDS(file.path(shap_dir,f))
    
    #Remove intercept and BIAS
    shap <- shap[,!(colnames(shap) %in% c("BIAS","(Intercept)")),drop=FALSE]

    #For each type of categorization
    all_shap_for_gene <- NULL
    all_shap_for_gene_reg <- NULL
    for(cur_clustering in colnames(summarize_shap_for_cluster)){
      for_category <- summarize_shap_for_cluster[,cur_clustering]

      if(nrow(shap)!=nrow(summarize_shap_for_cluster)){
        stop("Error; number of cells mismatch")
      }
      
      #Summarize values for each group of cells
      all_shap_for_gene <- rbind(
        all_shap_for_gene, 
        SummarizeShapByCluster_summarizeClusterNoreg(shap, for_category, onecat, curgene))
      all_shap_for_gene_reg <- rbind(
        all_shap_for_gene_reg, 
        SummarizeShapByCluster_summarizeClusterReg(shap, for_category, onecat, curgene))
      
    }
    all_shap[[curgene]] <- all_shap_for_gene
    all_shap_reg[[curgene]] <- all_shap_for_gene_reg
    
  }
  
  ## Set up and create output directories
  dir_shapsummary <- file.path(nando_dir,"summarizedshaps")
  dir_shapsummary_reg <- file.path(nando_dir,"summarizedshapregs")
  
  if(!file.exists(dir_shapsummary)){
    dir.create(dir_shapsummary)
  }

  if(!file.exists(dir_shapsummary_reg)){
    dir.create(dir_shapsummary_reg)
  }

  ## Save summarized SHAPs, region ignored
  df <- do.call("rbind", all_shap)
  rownames(df) <- NULL  
  saveRDS(df, file.path(dir_shapsummary,paste(task_id,".RDS")))

  ## Save summarized SHAPs, with region kept
  df <- do.call("rbind", all_shap_reg)
  rownames(df) <- NULL  
  saveRDS(df, file.path(dir_shapsummary_reg,paste(task_id,".RDS")))

}





#' Compute summarized SHAP values for a number of clusters.
#' Note that cells can be part of more than one cluster.
#' 
#' This function should be used in an R script run in a SLURM batch array.
#' It will use SLURM_ARRAY_TASK_COUNT and SLURM_ARRAY_TASK_ID to divide the
#' summarization between different nodes
#' 
#' @param nando_dir Location of all Nando files
#' @return Nothing
#' 
#' @export
SummarizeShapByClusterParallel <- function(nando_dir){

  shap_dir <- file.path(nando_dir, "shap")

  #Read list of clusterings  
  summarize_shap_for_cluster <- read.csv(file.path(nando_dir,"clusterings.csv"), row.names = "X")

  #Which files are there to process?
  list_shap_files <- list.files(shap_dir)
  list_shap_files <- list_shap_files[str_starts(list_shap_files, "shap ")]
  
  #Figure out how to divide it all for SLURM
  num_tasks <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
  task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  if(is.na(num_tasks)){
    print("SLURM jobarray not detected; running in a single thread")
    num_tasks <- 1  
    task_id <- 0
  } else {
    print(paste("Divide slurm",num_tasks,task_id))
  }
  
  ## Figure out which files this particular job should do
  shap_file_for_task <- (1:length(list_shap_files))%%num_tasks
  list_files_to_process <- list_shap_files[shap_file_for_task==task_id]
  print("Will process ===========================")
  print(list_files_to_process)
  print("========================================")
  
  ## Do the actual work
  SummarizeShapByCluster(
    nando_dir, shap_dir, list_files_to_process, summarize_shap_for_cluster, taskid)
  
  print("Finished summarizing shaps")
}





################################################################################
################### Network preparation ########################################
################################################################################

#' Make a transition matrix
#' 
#' @param onegrn A NandoNet
#' @param subsetGene,subsetTF Subset of genes to include. Default is all
#' @param includegenes Names of genes to forcefully add as row/columns in the matrix
#' @return Transition matrix, in sparse matrix format to save memory
#' 
#' @export
#' 
TransitionMatrix.NandoNetwork <- function(nandonet, subsetGene=NULL, subsetTF=NULL, includegenes=c()){
  TransitionMatrix(
    nandonet@shap, 
    subsetGene=subsetGene,
    subsetTF=subsetTF,
    includegenes=c())
}


#' Make a transition matrix from SHAP scores. Primarily for internal use
#' but can be used for various subsets of the SHAP entries in a network
#' 
#' @param onegrn Long format data frame of SHAPs
#' @param subsetGene,subsetTF Subset of genes to include. Default is all
#' @param includegenes Names of genes to forcefully add as row/columns in the matrix
#' @return Transition matrix, in sparse matrix format to save memory
#' 
#' @export
#' 
TransitionMatrix.data.frame <- function(onegrn, subsetGene=NULL, subsetTF=NULL, includegenes=c()){
  
  #Perform subsetting
  if(!is.null(subsetGene)){
    onegrn <- onegrn[onegrn$gene %in% subsetGene,]
  }
  if(!is.null(subsetGene)){
    onegrn <- onegrn[onegrn$tf %in% subsetTF,]
  }
  
  #Remove gain=0 connections as they just cause issues
  onegrn <- onegrn[onegrn$gain!=0,]
  
  #ignore regions; sum up gains
  onegrn <- sqldf::sqldf("select sum(gain) as gain, gene, tf from onegrn group by gene,tf")
  
  #and normalized per-gene gains are also log normal. however, p-values easiest to motivate on linear space
  onegrn_norm <- merge(sqldf::sqldf("select sum(gain) as sumgain, gene from onegrn group by gene"), onegrn)
  onegrn_norm$gain <- onegrn_norm$gain/onegrn_norm$sumgain
  
  #Figure out where all entries would sit in a sparse matrix
  allgenes <- sort(unique(c(onegrn$tf, onegrn$gene, includegenes)))
  numgenes <- length(allgenes)
  onegrn_norm <- merge(merge(onegrn_norm[,c("gene","tf","gain")],data.frame(i=1:numgenes, gene=allgenes)),data.frame(j=1:numgenes, tf=allgenes))
  
  #Create sparse matrix
  newmat <- sparseMatrix(
    onegrn_norm$i,
    onegrn_norm$j,
    x=onegrn_norm$gain, 
    dims=c(numgenes, numgenes), dimnames=list(allgenes,allgenes))

  #if there are no outgoing edges from a node, need to add one to itself at minimum
  #print(paste("number of absorbing states:", sum(rowSums(newmat)==0)))
  diag(newmat)[which(rowSums(newmat)==0)] <- 1

  newmat  
}

#' Helper function: Convert a sparse transition matrix into one for the R markovchain library
#' 
#' @param mat Transition matrix in sparse format
#' @return A markovchain object
#' 
transitionmatrix_to_mc <- function(mat){
  new("markovchain",transitionMatrix=as.matrix(mat),
               states=colnames(mat),
               name="mc")
}

#' Helper function: Prepare one network from a long list of SHAP.
#' Sets up gene categories
#' 
#' @param shaps Long format of shap scores
#' 
CreateOneNandoNetwork <- function(shaps){
  nandonet <- new("NandoNetwork", shap=shaps)
  nandonet <- compute_gene_classes(nandonet)
  nandonet
}


#' Load networks previously summarized by cluster
#'
#' @param nando_dir The directory holding per-cluster SHAP scores
#' @param keep_clusters Optional list of clusters to extract; default is all of them
#' @return A ListofNandoNetwork object
#' 
#' @export
LoadNandoNetworks <- function(nando_dir, keep_clusters=NULL){

  #Load all genes  
  print("Loading all SHAP scores...")
  shaps_dir <- file.path(nando_dir,"summarizedshaps")
  ct_shaps <- list()
  for(f in list.files(shaps_dir)){
    if(str_ends(f,"RDS")){
      print(f)
      ct_shaps[[f]] <- readRDS(file.path(shaps_dir,f))
    }
  }
  ct_shaps <- do.call("rbind", ct_shaps)
  rownames(ct_shaps) <- NULL

  #Pick all clusters by default
  if(is.null(keep_clusters)){
    keep_clusters <- unique(ct_shaps$category)
  } 

  #Organize networks by cluster
  print("Preparing network for each cell type...")
  list_nets <- foreach (cur_cat = keep_clusters, .combine=c, .verbose = F) %do% {
    
    #Could also go grab from each list directly here to skip rbind above
    
    print(cur_cat)
    CreateOneNandoNetwork(ct_shaps[ct_shaps$category==cur_cat,!(colnames(ct_shaps) %in% c("clustering","category"))])
    #ComputeHittingProbability.NandoNetwork(nandonets@nets[[x]])
  }
  names(list_nets) <- keep_clusters
  new("ListOfNandoNetwork", nets=list_nets)
}




################################################################################
################### Probabilistic interpretation ###############################
################################################################################


#'
#' Helper function: Computes what type category each gene is in the network.
#' In particular, finds the irreducible subset (assuming there is only one!)
#' 
#' @param nandonet One network
#' 
compute_gene_classes <- function(nandonet){

  shaps <- nandonet@shap
  
  #Figure out which are TFs and non-TFs
  nandonet@list_all_gene <- unique(c(shaps$tf, shaps$gene))
  nandonet@list_tfs <- sort(unique(shaps$tf))
  nandonet@list_nontfs <- sort(unique(shaps$gene[!(shaps$gene %in% nandonet@list_tfs)]))

  #Find distances between genes
  onegraph <- igraph::graph_from_edgelist(as.matrix(shaps[,c("gene","tf")]), directed = TRUE)
  d <- igraph::distances(onegraph, mode="out")
  d <- d[nandonet@list_tfs,nandonet@list_tfs]  #From TFs, how many jumps to other genes?

  #Irreducible TFs must be the ones being able to reach the big cluster
  num_reachable_from_tf <- rowSums(!is.infinite(d))    
  #nandonet@list_irreducible <- names(num_reachable_from_tf)[num_reachable_from_tf >= 0.5*max(num_reachable_from_tf)]
  almost_irreducible <- names(num_reachable_from_tf)[num_reachable_from_tf >= 0.5*max(num_reachable_from_tf)]
  
  #Rough approximation above. Some TFs pointing into this network will also be kept and most be removed in a second step.
  #These are any TFs that are unreachable from some TF
  d <- d[almost_irreducible,almost_irreducible]
  nandonet@list_irreducible <- colnames(d)[apply(!is.infinite(d),2,all)]
  
  #Best to actually validate the choice
  if(FALSE){
    shaps <- shaps[shaps$tf %in% nandonet@list_irreducible & shaps$gene %in% nandonet@list_irreducible,]
    onegraph <- igraph::graph_from_edgelist(as.matrix(shaps[,c("gene","tf")]), directed = TRUE)
    d <- igraph::distances(onegraph, mode="out")
    plot(onegraph)
    #d <- d[nandonet@list_irreducible,nandonet@list_irreducible]
  }
  
  nandonet
}


#' Compute hitting probabilities for a network.
#' 
#' A crucial trick is used internally to do this fast. HPs are first calculated for TFs alone, as
#' this requires a full linear equation to work out. Then nonTFs are calculated, and since these are
#' unreachable, their HPs can be worked out by simply weighing HPs from the TFs
#' 
#' @param nandonet One NandoNetwork
#' @return A NandoNetwork object with hitting probability set
#' 
#' @rdname ComputeHittingProbability
#' @export
#' @method ComputeHittingProbability ListOfNandoNetwork
ComputeHittingProbability.NandoNetwork <- function(nandonet){

  shaps <- nandonet@shap

  #First calculate this for all TFs. This will be a fairly dense but small matrix
  tmat_tf <- TransitionMatrix(
    shaps[shaps$gene %in% nandonet@list_tfs,], 
    includegenes = nandonet@list_tfs)
  tmat_tf.hp <- markovchain::hittingProbabilities(transitionmatrix_to_mc(tmat_tf)) #Slow part of this code. now really slow?
  
  #Prepare new hp matrix for all genes, and transfer computed hp's for TFs    #could be made sparse; but not affecting speed much
  full.hp <- matrix(0, length(nandonet@list_all_gene), length(nandonet@list_all_gene))
  colnames(full.hp) <- nandonet@list_all_gene
  rownames(full.hp) <- nandonet@list_all_gene
  full.hp[rownames(tmat_tf.hp), colnames(tmat_tf.hp)] <- tmat_tf.hp    
  
  #Pull out edges and calculate p values related to nonTFs
  nontf_grn <- shaps[shaps$gene %in% nandonet@list_nontfs,]
  nontf_grn <- sqldf::sqldf("select sum(gain) as gain, gene, tf from nontf_grn group by gene,tf")
  nontf_grn <- merge(sqldf::sqldf("select sum(gain) as sumgain, gene from nontf_grn group by gene"), nontf_grn)
  nontf_grn$p <- nontf_grn$gain/nontf_grn$sumgain
  nontf_grn
  
  #Get hybrid transition matrix, from genes to first upstream TF
  nontf_grn_mat <- reshape2::acast(nontf_grn, gene~tf, value.var = "p", fill=0)
  
  #Prepare to collect: upstream -> final destination
  firstupstream.hp <- tmat_tf.hp[colnames(nontf_grn_mat),,drop=FALSE]
  
  #For each nonTF, compute weighted hp of direct regulators
  full.hp[rownames(nontf_grn_mat), colnames(tmat_tf.hp)] <- nontf_grn_mat %*% firstupstream.hp
  
  if(!all(full.hp[,nandonet@list_nontfs]==0)){
    #Sanity check: nonTFs should not be reachable
    print("reachability issue")
    stop()
  }
  
  #90% sparse
  nandonet@hp <- as(full.hp, "sparseMatrix")
  nandonet
}




#' Compute hitting probabilities for each network.
#' This function is able to parallelize using foreach if set up properly
#' 
#' @param nandonets A ListOfNandoNetwork
#' @return A ListOfNandoNetwork object with hitting probability set
#' 
#' @rdname ComputeHittingProbability
#' @export
#' @method ComputeHittingProbability ListOfNandoNetwork
ComputeHittingProbability.ListOfNandoNetwork <- function(nandonets){

  outlist <- foreach (x = names(nandonets@nets), .combine=c, .verbose = F) %do% {
    print(x)
    ComputeHittingProbability.NandoNetwork(nandonets@nets[[x]])
  }
  names(outlist) <- names(nandonets@nets)
  nandonets@nets <- outlist
  nandonets
}





#' Compute steady state for a network. Only the irreducible TFs will be considered
#' 
#' @param nandonet One network
#' @return A NandoNetwork object with steady state set
#'
#' @rdname ComputeSteadyState
#' @export
#' @method ComputeSteadyState NandoNetwork
ComputeSteadyState.NandoNetwork <- function(nandonet){
  shaps <- nandonet@shap
  tmat_tf <- TransitionMatrix(shaps[
    shaps$gene %in% nandonet@list_irreducible &
      shaps$tf %in% nandonet@list_irreducible,])
  ss <- markovchain::steadyStates(transitionmatrix_to_mc(tmat_tf))
  
  nandonet@ss <- ss[1,nandonet@list_irreducible]
  nandonet
}



#' Compute steady state for each network.
#' This function is able to parallelize using foreach if set up properly
#' 
#' @param nandonets A list of networks
#' @return A ListOfNandoNetwork object with steady state set
#' 
#' @rdname ComputeSteadyState
#' @export
#' @method ComputeSteadyState ListOfNandoNetwork
ComputeSteadyState.ListOfNandoNetwork <- function(nandonets){
  
  outlist <- foreach (x = names(nandonets@nets), .combine=c, .verbose = F) %do% {
    print(x)
    ComputeSteadyState.NandoNetwork(nandonets@nets[[x]])
  }
  names(outlist) <- names(nandonets@nets)
  nandonets@nets <- outlist
  nandonets
}



################################################################################
################### Walktrap analysis ##########################################
################################################################################




#' Compute a series of P^t to figure out about the neighborhood of genes, for one network
#' 
#' @param nandonet One network
#' @param num_mul Calculate for P^1 ... P^num_mul
#' @param only_irreducible_TF Should only irreducible TFs be kept, and other TFs removed?
#' @return A NandoWalktrap object
#' 
#' @rdname ComputeNandoWalktrap
#' @export
#' @method ComputeNandoWalktrap NandoNetwork
ComputeNandoWalktrap.NandoNetwork <- function(nandonet, num_mul, only_irreducible_TF=TRUE){

  #Get the transition matrix  
  if(only_irreducible_TF){
    cur_tmat <- TransitionMatrix(nandonet@shap[nandonet@shap$tf %in% nandonet@list_irreducible,])
  } else {
    cur_tmat <- TransitionMatrix(nandonet)
  }

  entropy_over_steps <- matrix(nrow=nrow(cur_tmat), ncol=num_mul)
  rownames(entropy_over_steps) <- rownames(cur_tmat)
  
  #Compute distributions. p_ij, so each column is the distribution from one gene
  mat <- as(cur_tmat, "TsparseMatrix") #Sparse format that consists of triplets, easy to work with
  cumulative_mat <- mat
  totmat <- mat  
  totmat_meta <- data.frame(step=rep(1, ncol(mat)))  
  entropy_over_steps[,1] <- fast_entropy_sparsematrix(cumulative_mat)
  
  for(i in 2:num_mul){
    cumulative_mat <- cumulative_mat %*% mat
    entropy_over_steps[,i] <- fast_entropy_sparsematrix(cumulative_mat)
    totmat <- rbind(totmat, cumulative_mat)
    totmat_meta <- rbind(totmat_meta, 
                         data.frame(step=rep(i,ncol(mat))))
  }
  totmat_meta$fromgene <- rep(colnames(mat),num_mul)
  
  #Matrix of distributions, and metadata (fromgene, step) returned
  new("NandoWalktrap",mat=totmat, meta=totmat_meta, entropy=entropy_over_steps)
}



#' Compute a series of P^t to figure out about the neighborhood of genes, for a list of networks
#' 
#' @param nandonets List of network
#' @param num_mul Calculate for P^1 ... P^num_mul
#' @param only_irreducible_TF Should only irreducible TFs be kept, and other TFs removed?
#' @return A ListOfNandoWalktrap object
#' 
#' @rdname ComputeNandoWalktrap
#' @export
#' @method ComputeNandoWalktrap ListOfNandoNetwork
ComputeNandoWalktrap.ListOfNandoNetwork <- function(nandonets, num_mul, only_irreducible_TF=TRUE){
  
  outlist <- foreach (x = names(nandonets@nets), .combine=c, .verbose = F) %do% {
    print(x)
    ComputeNandoWalktrap.NandoNetwork(nandonets@nets[[x]], num_mul, only_irreducible_TF)
  }
  names(outlist) <- names(nandonets@nets)
  
  new("ListOfNandoWalktrap",nets=outlist)
}


  
#' Helper function: Compute entropy over a transition matrix
#' 
#' @param mat Matrix, which must be in format of TsparseMatrix
#' @return One-column matrix with the entropy from each gene
#' 
fast_entropy_sparsematrix <- function(mat){

  #Calculate entropy for each gene while keeping sparsity. x!=0 so no extra checking needed
  df <- data.frame(
    i=mat@i,
    pe=mat@x * log(mat@x)
  )
  df <- sqldf::sqldf("select -sum(pe) as entropy, i from df group by i")
  
  entropy_over_steps <- matrix(nrow=nrow(mat), ncol=1)
  rownames(entropy_over_steps) <- names(mat)
  
  entropy_over_steps[df$i+1, 1] <- df$entropy
  entropy_over_steps
}






#' Helper function to run basic Seurat initial processing
#' 
#' @param gdata A Seurat object
#' @return A Seurat object after PCA, UMAP and leiden clustering
#' 
BasicSeuratProcessingGdata <- function(gdata){
  #gdata <- Seurat::NormalizeData(object = gdata)
  gdata <- Seurat::FindVariableFeatures(gdata, selection.method="disp", nfeatures=nrow(gdata))
  gdata <- Seurat::ScaleData(gdata)
  gdata <- Seurat::RunPCA(gdata)#, npcs=min(10, length(gdata@assays$RNA@var.features)))
  gdata <- Seurat::FindNeighbors(gdata)#, k.param=1:min(dims, length(gdata@assays$RNA@var.features)))
  gdata <- Seurat::FindClusters(gdata)
  gdata <- Seurat::RunUMAP(gdata, dims = 1:20)
  gdata  
}


#' Convert Walktrap into Seurat to explore the cell state space using common single cell tools.
#' 
#' @param walktraps A ListOfNandoWalktraps
#' @param keep_steps Which steps to subset by; default is to keep all
#' @param keep_genes Which genes to subset by; default is to keep all
#' @return A Seurat object
#' 
#' @export
#' 
NandoWalktrapToSeurat <- function(walktraps, keep_steps=NULL, keep_genes=NULL){
  
  #Default is all genes
  if(is.null(keep_genes)){
    keep_genes <- rownames(walktraps@nets[[1]]@mat)
  }
  
  #Default is all steps
  if(is.null(keep_steps)){
    keep_steps <- unique(walktraps@nets[[1]]@meta$step)
  }

  #Gather all meta  
  tot_meta <- foreach(cur_ct = names(walktraps@nets), .combine = "rbind", .verbose = F) %do% {
    wt <- walktraps@nets[[cur_ct]]
    keep <- wt@meta$step %in% keep_steps & wt@meta$fromgene %in% keep_genes
    partmeta <- wt@meta[keep,,drop=FALSE]
    partmeta$cluster <- cur_ct
    partmeta
  }
  
  #Gather all matrices
  tot_mat <- foreach(cur_ct = names(walktraps@nets), .combine = "rbind", .verbose = F) %do% {
    wt <- walktraps@nets[[cur_ct]]
    keep <- wt@meta$step %in% keep_steps & wt@meta$fromgene %in% keep_genes
    wt@mat[keep,,drop=FALSE]
  }
  
  #Create seurat object
  gdata <- Seurat::CreateSeuratObject(counts = t(tot_mat))
  gdata@meta.data$fromgene <- tot_meta$fromgene
  gdata@meta.data$cluster <- tot_meta$cluster
  gdata@meta.data$step <- tot_meta$step
  
  # Cluster and compare all genes by neighbor similarity
  BasicSeuratProcessingGdata(gdata)
}






#' Gather entropies across all networks for chosen genes and steps
#' 
#' @param walktraps A ListOfNandoWalktraps
#' @param keep_genes Which genes to extract. Default is all (which is a lot!)
#' @param keep_steps Which steps to keep. Default is all
#' 
#' @return Molten data.frame with all values
#' 
#' @export 
GetWalktrapEntropy <- function(walktraps, keep_genes=NULL, keep_steps=NULL){
  
  if(is.null(keep_steps)){
    keep_steps <- 1:ncol(walktraps@nets[[1]]@entropy)
  }
  
  if(is.null(keep_genes)){
    keep_genes <- rownames(walktraps@nets[[1]]@entropy)
  }
  
  df <- foreach(cur_ct = names(walktraps@nets), .combine = "rbind", .verbose = F) %do% {
    wt <- walktraps@nets[[cur_ct]]
    df <- as.data.frame(wt@entropy[rownames(wt@entropy) %in% keep_genes,,drop=FALSE])
    colnames(df) <- 1:ncol(df)
    df <- df[,colnames(df) %in% keep_steps]
    df <- reshape2::melt(as.matrix(df))
    df$cluster <- cur_ct
    df
  }
  colnames(df) <- c("gene","step","entropy","cluster")
  df
}




################################################################################
################### Overlap analysis ###########################################
################################################################################


#' Compute network similarities
#'
#' @details
#' Similarities done using fuzzy logic, as our edges have weights [0,1] (if an edge not in list then p=0)
#' https://en.wikipedia.org/wiki/Fuzzy_set_operations  
#' discussion of |A| here: https://gvpress.com/journals/IJEIC/vol4_no1/2.pdf
#' 
#' Jaccard index: https://en.wikipedia.org/wiki/Jaccard_index
#' |A intersect B| / |A union B|
#'
#' SSI: https://en.wikipedia.org/wiki/Overlap_coefficient
#' |A intersect B| / min(|A|, |B|)
#'
#' @param nandonets ListOfNandoNetwork
#' @param method The similarity metric to apply
#' @return Matrix of similarity scores
#' 
#' @export
ComputeNandoNetworkSimilarity <- function(nandonets, method=c("jaccard","ssi")){
  method <- match.arg(method)

  #Obtain all network edges. Long format
  print("Gathering all networks in comparable fromat")
  tmat_long <- foreach (x = names(nandonets@nets), .combine=rbind, .verbose = F) %do% {
    print(x)
    out <- MeltSparsematrix(TransitionMatrix(nandonets@nets[[x]]))
    out$net <- x
    data.frame(
      net=out$net,
      p=out$value,
      edge=paste(out$row,out$col)
    )
  }
  tocompare <- reshape2::acast(tmat_long, edge~net, value.var="p", fill=0)
  
  
  #Allocate output matrix
  print("Doing comparisons")
  numnet <- length(nandonets@nets)
  out <- matrix(nrow=numnet, ncol=numnet)
  colnames(out) <- names(nandonets@nets)
  rownames(out) <- names(nandonets@nets)
  
  #Perform comparisons.
  #Note: pmin/pmax seem to do dense expansion. Comparison can be made in linear time and extremely fast,
  #but would likely need to be implemented in C
  
  #Potential optimization is to use the symmetry of the matrix. Compute half, add transposed. fill diagonal with 1
  
  #Fuzzy Jaccard index
  if(method=="jaccard"){
    for(i in 1:numnet){
      for(j in 1:numnet){
        print(paste(i,j))
        out[i,j] <- sum(pmin(tocompare[,i], tocompare[,j])) / sum(pmax(tocompare[,i], tocompare[,j]))
      }
    }
  }
  
  #Fuzzy SSI
  if(method=="ssi"){
    for(i in 1:numnet){
      for(j in 1:numnet){
        print(paste(i,j))
        out[i,j] <- sum(pmin(tocompare[,i],tocompare[,j])) / min(sum(tocompare[,i]),sum(tocompare[,j]))
        #Could precompute sums at the end to speed up marginally
      }
    }
  }  
  
  out
}



################################################################################
########################### Igraph #############################################
################################################################################



#' Convert a transition matrix into a weighted directed igraph object
#' 
#' @param tmat Transition matrix in sparse format
#' @return igraph object
#' 
#' @export
TransitionMatrixToIgraph <- function(tmat){
  sp <- as(tmat, "TsparseMatrix")
  onegraph <- igraph::graph_from_edgelist(
    as.matrix(data.frame(from=colnames(sp)[sp@i+1],to=colnames(sp)[sp@j+1])), 
    directed = TRUE)
  igraph::E(onegraph)$weight <- sp@x
  onegraph
}



################################################################################
################### GO analysis ################################################
################################################################################







#' Convert hitting probabilities into Seurat to explore the cell state space using common single cell tools.
#' 
#' @param walktraps A ListOfNandoNetwork
#' @param keep_genes Which genes to subset by; default is to keep all
#' @return A Seurat object
#' 
#' @export
#' 
HittingProbabilitiesToSeurat <- function(nandonets, keep_genes_from=NULL, keep_genes_to=NULL){
  
  #Default is all genes
  if(is.null(keep_genes_from)){
    keep_genes_from <- rownames(nandonets@nets[[1]]@hp)
  }
  if(is.null(keep_genes_to)){
    keep_genes_to <- colnames(nandonets@nets[[1]]@hp)
  }

  #Gather all meta
  tot_meta <- foreach(cur_ct = names(nandonets@nets), .combine = "rbind", .verbose = F) %do% {
    net <- nandonets@nets[[cur_ct]]
    meta <- data.frame(
      fromgene=rownames(net@hp)[rownames(net@hp) %in% keep_genes_from]
    )
    meta$cluster <- cur_ct
    meta
  }

  #Gather all matrices
  tot_mat <- foreach(cur_ct = names(nandonets@nets), .combine = "rbind", .verbose = F) %do% {
    net <- nandonets@nets[[cur_ct]]
    net@hp[
      rownames(net@hp) %in% keep_genes_from,
      colnames(net@hp) %in% keep_genes_to,
      drop=FALSE]
  }

  #Create seurat object
  gdata <- Seurat::CreateSeuratObject(counts = t(tot_mat))
  gdata@meta.data$fromgene <- tot_meta$fromgene
  gdata@meta.data$cluster <- tot_meta$cluster

  # Cluster and compare all genes by neighbor similarity
  BasicSeuratProcessingGdata(gdata)
}





################################################################################
################### Steady state perturbation analysis #########################
################################################################################



#' In moving from one state state pi_i to pi_j, which edges does the probability
#' mainly flow over?
#' 
#' @param netA,netB Two NandoNetworks to be compared
#' @return Matrix, from each gene showing how much flowed to another
#' 
ComputeSteadyStateChangeEdgeflow <- function(netA, netB){
  #netA <- nandonets@nets[[1]]
  #netB <- nandonets@nets[[2]]

  matA <- TransitionMatrix(netA, subsetGene=netA@list_irreducible, subsetTF=netA@list_irreducible)
  matB <- TransitionMatrix(netB, subsetGene=netB@list_irreducible, subsetTF=netB@list_irreducible)
  
  fromSS <- netA@ss  
  toSS <- netB@ss

  # Check some simplifying assumptions  
  if(!all(names(fromSS)==names(toSS))){
    stop("to implement unequal ss")
  }
  if(!all(colnames(matA)==names(fromSS))){
    stop("eep.")
  }

  # Integrate flow, moving from one steady state to another  
  equilibrium_flow <- toSS * matB  #not sure if this is correct
  accum_flow <- equilibrium_flow*0
  curstate <- fromSS
  for(i in 1:20){
    curflow <- curstate * matB   #not sure if this is correct
    delta_flow <- curflow - equilibrium_flow
    accum_flow <- accum_flow + delta_flow
    #print(sum(abs(delta_flow)))
    curstate <- curstate %*% matB
    curstate <- as.double(curstate) #ugly but matrix multiplication result in matrix; issuse with cycling in elementwise mul
  }
  
  accum_flow
}


#' Melt a sparse matrix to long format
MeltSparsematrix <- function(mat){
  out <- as(mat, "TsparseMatrix")
  data.frame(
    row=rownames(mat)[out@i+1],
    col=colnames(mat)[out@j+1],
    value=mat@x
  )
}




  