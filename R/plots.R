library(ggplot2)



# steal from website



################################################################################
################### Probability matrices #######################################
################################################################################



#' For each network, compute hitting probability, starting from steady state
#' 
#' @param nandonets A ListOfNandoNetworks
#' @return A sparse matrix of hitting probabilities; each row is one network
#' 
#' @export
#' 
ComputeHittingProbabilityFromSS <- function(nandonets){
  hp <- foreach(onecat = names(nandonets@nets), .combine = "rbind",  .verbose = F) %do% {
    nandonet <- nandonets@nets[[onecat]]
    
    #Get the steady state of same format as hitting probability matrix. Keep sparsity
    longss <- MatrixExtra::emptySparse(nrow=1, ncol=ncol(nandonet@hp), format="T")
    colnames(longss) <- colnames(nandonet@hp)
    longss[1,names(nandonet@ss)] <- nandonet@ss
    
    #ss hitting probability is simply a weighting using ss. Return a one-row matrix
    longss %*% nandonet@hp
  }
  rownames(hp) <- names(nandonets@nets)
  hp
}

#' For each network, extract steady states. Return as a matrix
#' 
#' @param nandonets A ListOfNandoNetworks
#' @return Matrix of steady states; each row is one network
#' 
#' @export
#' 
SteadyStateMatrix <- function(nandonets){
  ss <- foreach(onecat = names(nandonets@nets), .combine = "rbind",  .verbose = F) %do% {
    nandonet <- nandonets@nets[[onecat]]
    data.frame(
      gene=names(nandonet@ss),
      cluster=onecat,
      p=nandonet@ss)
  }
  reshape2::acast(ss, cluster~gene, value.var="p", fill = 0)
}


#' Plot a matrix consisting of probabilities
#' 
#' @param probs A matrix of probabilities; can be sparse or dense
#' @param min.pmean A lower cutoff on p values
#' @param dolog Show log10 p values or not
#' @return A ggplot2 object
#' 
#' @import ggplot2
#' @export
#' 
PlotTopProbabilityMatrix <- function(probs, min.pmean=1e-2, dolog=TRUE){
  #Order genes by average probability
  pmean <- colMeans(probs)
  newo <- order(pmean,decreasing = FALSE)
  probs <- probs[,newo]
  pmean <- pmean[newo]
  
  #Filter genes
  probs <- probs[,pmean>min.pmean]

  #Produce plot  
  long_probs <- reshape2::melt(probs)
  colnames(long_probs) <- c("cluster","gene","p")
  if(dolog){
    ggplot(long_probs, aes(cluster, gene, fill=log10(p))) + geom_tile() + labs(y="", x = "", fill = "Log10 p") #+ scale_y_reverse()
  } else {
    ggplot(long_probs, aes(cluster, gene, fill=p)) + geom_tile() + labs(y="", x = "", fill = "p") #+ scale_y_reverse()
  }
}





################################################################################
############################ Walktrap plots ####################################
################################################################################

#' Gather and plot entropy for chosen genes over time
#' 
#' @param for_genes Genes to show
#' @return A ggplot2 object
#' 
#' @import ggplot2
#' @export
PlotWalktrapEntropy <- function(for_genes){
  entr <- GetWalktrapEntropy(walktraps, for_genes)
  entr$cluster_gene <- paste(entr$cluster, entr$gene)
  ggplot(entr, aes(step, entropy, group=cluster_gene, color=gene)) + geom_line() 
}




# Plot of trajectories in seurat object



# plot similarity between networks
