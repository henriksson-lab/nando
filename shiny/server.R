library(RSQLite)
library(cowplot)


if(FALSE){
  #To run this app
  library(shiny)
  runApp(".")
}


################################################################################
########### Hitting probabilities of the full network ##########################
################################################################################


###########################
# Plot hitting probabilities from a given gene
PlotHittingProbabilityForGeneAsMatrixH5 <- function(genename, p_cutoff=0.1, use_clusters){
  
  #note; gene=tf, in the sense that we cannot tell what is a tf
  
  probs <- ReadProbabilityMatricesH5ForGene(h5hp, use_clusters, genename)
  
  # if(!donoraveraged){
  #   #Arrange in sensible order ... need donor be last?
  #   probs$ct <- str_split_fixed(probs$category," ",2)[,2]
  #   probs$donor <- str_split_fixed(probs$category," ",2)[,1]
  #   ct_table <- sqldf("select distinct ct, category from probs order by ct,donor")    ##unique(probs[,c("ct","category")])
  #   probs$category <- factor(probs$category, levels=ct_table$category)
  # }
  
  if(is.null(probs) || nrow(probs)==0){
    print(paste("no gene2",genename))
    return(paste("No data for gene ",genename))
  }
  
  meanp <- sqldf("select gene, avg(p) as meanp, max(p) as maxp from probs group by gene order by meanp desc")
  probs$gene <- factor(probs$gene, levels=rev(meanp$gene))
  
  ptot <- ggplot(probs[probs$gene %in% meanp$gene[meanp$maxp > p_cutoff],,drop=FALSE], aes(category, gene, fill=p)) + 
    geom_tile() +
    xlab("") + ylab("") +
    scale_x_discrete(guide = guide_axis(angle = 45))
  
  ptot
}
#plot_hp_for_gene_asmatrix("BCL6")





################################################################################
########### Functions to compute shap/p across regions #########################   
################################################################################


#' Compare SHAP p-values for two cell types at each location
#' 
#' @param shapregs A data frame with all p values
#' @param ct1,ct2 The two cell types to compare
#' @return delta p
#' 
CompareShapreg2CT <- function(shapregs,ct1,ct2){
  df1 <- shapregs[shapregs$category==ct1, c("tf","p","chrom","start","end")]
  df2 <- shapregs[shapregs$category==ct2, c("tf","p","chrom","start","end")]
  statecomp <- merge(
    df1,df2, 
    by=c("tf","chrom","start","end"), 
    all=TRUE)
  
  statecomp$p.x[is.na(statecomp$p.x)] <- 0
  statecomp$p.y[is.na(statecomp$p.y)] <- 0
  statecomp$deltap <- statecomp$p.x - statecomp$p.y
  statecomp
}


#' Plot differences between two cell types in terms of p
#' 
#' @param statecomp A previously made comparison
#' @param view_from,view_to Genomic range (integers; chromosome not needed)
#' @return A ggplot2 object
#' 
PlotShaprepComparisonTrack <- function(statecomp, view_from=NULL, view_to=NULL){
  #Show as points, up or down
  totp <- ggplot(statecomp, aes(start, deltap, color=tf, label=tf)) +
    xlab("") + ylab("delta p") +
    #geom_jitter() + 
    geom_text() + theme(legend.position="none")
  
  if(!is.null(view_from)){
    totp <- totp + xlim(view_from, view_to)
  }
  
  totp
}



################################################################################
########### Plotting of genome tracks ##########################################
################################################################################

# read clustering file to get all options. need function
#available_clusterings <- c("pred.dice","atac_clusters")


#' Plot genome tracks including SHAPs at enhancers
#' 
#' @param nando_dir Directory with Nando data
#' @param genename Which gene to show
#' @param use_clustering Which clustering to show tracks over
#' @param compare_ct1,compare_ct2 Two cell types to compare
#' @return A ggplot2 object
#' 
shiny_plot_genometracks <- function(nando_dir,
                                    genename, 
                                    use_clustering,
                                    compare_ct1, compare_ct2){
  
  #Load SHAPs for a gene and region
  shapregs <- ImportShapregsForGene(nando_dir, genename)
  
  if(is.null(shapregs) || nrow(shapregs)==0){
    print(paste("no gene",genename))
    return(paste("No data for gene ",genename))
  }
  
  #need to fill in 0s?
  
  #Compare two cell types for one gene
  statecomp <- CompareShapreg2CT(shapregs, compare_ct1, compare_ct2)
  
  #Adjust view location
  view_start <- min(statecomp$start)-500
  view_end <- max(statecomp$end)+500
  the_chrom <- statecomp$chrom[1]
  
  
  list_of_plots <- list()  
  
  #Produce adata-based plot tracks
  if(exists("adata")){
    #Potentially dangerous; affecting a global or not?
    Idents(adata) <- ReadNandoClustering(nando_dir)[,use_clustering]

    seurat_view_coord <- paste(the_chrom, view_start, view_end, sep="-")
    #Idents(adata) <- adata@meta.data$shiny_clustering  #Potentially dangerous; affecting a global or not?
    list_of_plots[["annotation"]] <- AnnotationPlot(object = adata, region = seurat_view_coord)
    list_of_plots[["coverage"]] <- CoveragePlot(object = adata, region = seurat_view_coord, annotation = TRUE, peaks = TRUE)
  } 
  
  #Produce SHAP tracks
  list_of_plots[["shap"]] <- PlotShaprepComparisonTrack(statecomp, view_start, view_end)
  
  #print(names(list_of_plots))
  
  #Assemble and return plots
  ptot <- cowplot::plot_grid(plotlist = list_of_plots, ncol = 1)
  ptot
}

if(FALSE){
  shiny_plot_genometracks(nando_dir, genename, use_clustering, compare_ct1 = "donor0 ALL", compare_ct2="donor1 ALL")
}




################################################################################
########### Steady state related network #######################################
################################################################################


###########################
# Plot core network steady state statistics
plot_corenetwork_ss_matrix <- function(themat, p_cutoff=0.01, include_tf=c()){
  #print(p_cutoff)
  #plot_corenetwork_somep_matrix(themat, p_cutoff, include_tf)
  PlotTopProbabilityMatrix(themat, p_cutoff, include_tf=include_tf)
}
#plot_corenetwork_ss_matrix(corenetwork_sump)



###########################
# Plot core network steady state statistics
plot_corenetwork_hp_matrix <- function(ssmat, hp.ssmat, p_cutoff=0.01, include_ss=TRUE, include_tf=c()){
  #print(head(ssmat))
  if(!include_ss){
    hp.ssmat <- hp.ssmat[,colnames(hp.ssmat) %in% include_tf | !(colnames(hp.ssmat) %in% colnames(ssmat)),drop=FALSE]
    #hp.ssmat <- hp.ssmat[hp.ssmat$tf %in% include_tf | !(hp.ssmat$tf %in% ssmat$tf),,drop=FALSE]
  }
  #plot_corenetwork_somep_matrix(hp.ssmat, p_cutoff, include_tf)
  PlotTopProbabilityMatrix(hp.ssmat, p_cutoff, include_tf=include_tf)
}
#plot_corenetwork_hp_matrix(corenetwork_sump, include_ss = FALSE)





################################################################################
########### The server declaration #############################################
################################################################################

server <- function(input, output, session) {

  ## Update option lists; long lists of options are best handled in this way  
  updateSelectizeInput(session, 'shapreg_gene', choices = shapregs_genelist, server = TRUE)
  updateSelectizeInput(session, 'geneinnet_gene', choices = shapregs_genelist, server = TRUE)
  updateSelectizeInput(session, 'network_gene', choices = available_tf_in_ss, server = TRUE)

  
  observeEvent(input$shapreg_useclustering, {
    print("update ct1 ct2")
    use_clustering <- input$shapreg_useclustering
    clusters_in_clustering <- available_clusters[[use_clustering]]
    updateSelectizeInput(session, 'shapreg_ct1', choices = clusters_in_clustering, server = TRUE)
    updateSelectizeInput(session, 'shapreg_ct2', choices = clusters_in_clustering, server = TRUE)
  })
  
  

  ################################################################################
  ########### Overall network browser ############################################
  ################################################################################

  output$networkPlotSS <- renderPlot(height=700, {
    
    #Should dynamically adapt plot size based on number of displayed genes
    include_tf <- input$network_gene
    if(is.null(include_tf)){
      include_tf <- c()
    }
    network_pcutoff <- 10**as.double(input$network_pcutoff)
    #corenetwork_sump <- get_corenetwork_summaryp(donoraveraged=TRUE) ## 666 no more

    use_clustering <- input$network_useclustering#"donor_ALL"
    
    ssmat <- ImportShinyNando.SteadyStates(nando_dir)
    clusters_in_clustering <- available_clusters[[use_clustering]]
    ssmat <- ssmat[rownames(ssmat) %in% clusters_in_clustering,]
    
    plot_corenetwork_ss_matrix(
      ssmat,
      include_tf = include_tf, 
      p_cutoff = network_pcutoff
    ) 
  })
  
  output$networkPlotHP <- renderPlot(height=700, {
      
    #TODO Should dynamically adapt plot size based on number of displayed genes
    include_tf <- input$network_gene
    if(is.null(include_tf)){
      include_tf <- c()
    }
    network_pcutoff <- 10**as.double(input$network_pcutoff)

    use_clustering <- input$network_useclustering
    clusters_in_clustering <- available_clusters[[use_clustering]]
    include_ss <- input$network_ss_in_hp
    
    ##SS matrix
    ssmat <- ImportShinyNando.SteadyStates(nando_dir)
    ssmat <- ssmat[rownames(ssmat) %in% clusters_in_clustering,,drop=FALSE]
    
    #HP.SS matrrix
    hp.ssmat <- ImportShinyNando.HP.SS(nando_dir)
    hp.ssmat <- hp.ssmat[rownames(hp.ssmat) %in% clusters_in_clustering,]
    
    plot_corenetwork_hp_matrix(
      ssmat,
      hp.ssmat, 
      include_ss = include_ss,
      include_tf = include_tf, 
      p_cutoff = network_pcutoff
    )
    

  })
  
  
  ################################################################################
  ########### Genome browser (shap per region viewer) ############################
  ################################################################################
  
  output$genomebrowserPlot <- renderPlot(height=1000, {
    
    if(is.na(input$shapreg_gene) | input$shapreg_gene==""){
      print("NA gene")
      return("NA gene")
    } else {
      
      shiny_plot_genometracks(
        nando_dir = nando_dir,
        genename = input$shapreg_gene, 
        compare_ct1 = input$shapreg_ct1,
        compare_ct2 = input$shapreg_ct2,
        use_clustering = input$shapreg_useclustering)
    }
    
  })
  
  
  
  ################################################################################
  ########### Gene in net browser ################################################
  ################################################################################
  
  #Direct neighbours
  output$geneinnetNeighbourPlot <- renderPlot(height=500, {

    
    if(is.na(input$geneinnet_gene) | input$geneinnet_gene==""){
      print("NA gene")
      return("NA gene")
    } else {

      use_clustering <- input$geneinnet_useclustering
      clusters_in_clustering <- available_clusters[[use_clustering]]
      
      network_pcutoff_direct <- 10**as.double(input$geneinnet_pcutoff_direct)
      PlotHittingProbabilityForGeneAsMatrixH5(
        genename = input$geneinnet_gene, 
        p_cutoff = network_pcutoff_direct, 
        use_clusters=clusters_in_clustering
      )  #max p rather than average?
      
    }
    
  })
  
  
  #Closest by proximity ........... ,  limited to core network?
  output$geneinneHittingPlot <- renderPlot(height=500, {
    
    use_clustering <- input$geneinnet_useclustering
    clusters_in_clustering <- available_clusters[[use_clustering]]
    
    if(is.na(input$geneinnet_gene) | input$geneinnet_gene==""){
      print("NA gene")
      return("NA gene")
    } else {

      network_pcutoff_hp <- 10**as.double(input$geneinnet_pcutoff_hp)
      PlotHittingProbabilityForGeneAsMatrixH5(
        use_clusters = clusters_in_clustering,
        genename = input$geneinnet_gene, 
        p_cutoff = network_pcutoff_hp
      )
    }
    
  })
  
}

