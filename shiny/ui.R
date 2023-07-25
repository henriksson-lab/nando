


################################################################################
########### Core network #######################################################
################################################################################
tab_network <- sidebarLayout(
  
  sidebarPanel(
    
    selectInput(
      inputId = "network_useclustering",
      label = "Clustering:",
      choices = names(available_clusters), 
      selected = names(available_clusters)[1]
    ),
    
    
    selectInput(
      inputId = "network_gene",
      label = "Include additional genes:",
      selectize = TRUE,
      multiple = TRUE,
      choices = c(""), 
      selected = NULL
    ),

    
    checkboxInput("network_ss_in_hp", "Show irreducible TFs in proximity matrix:", FALSE),
    
    sliderInput(
      inputId = "network_pcutoff", 
      label = "Log10 p cutoff:",
                min = -5, max = 0, value = -2,step = 0.1
    )#,
    
        
    # selectInput(
    #   inputId = "network_clustering",
    #   label = "Clustering:",
    #   choices = c("",available_clusterings),
    #   selected = "pred.dice"
    # )
    
    
  ),
  mainPanel(
    h3("Steady state (core network):"),
    plotOutput(outputId = "networkPlotSS", height = "700px"),
    br(),
    h3("Hitting probabilities from steady state (core network):"),
    plotOutput(outputId = "networkPlotHP")
  )
)

################################################################################
########### Genome browser #####################################################
################################################################################
tab_genomebrowser <- sidebarLayout(
  
  sidebarPanel(

    selectInput(
      inputId = "shapreg_useclustering",
      label = "Clustering:",
      choices = names(available_clusters), 
      selected = names(available_clusters)[1]
    ),
    
    selectInput(
      inputId = "shapreg_gene",
      label = "Gene:",
      selectize = TRUE,
      choices = c("CDKL5"), #shapregs_genelist,
      selected = "CDKL5"
    ),
    
    selectInput(
      inputId = "shapreg_ct1",
      label = "Compare cell type 1:",
      selectize = TRUE,
      choices = c(""),
      selected = ""
    ),
    
    selectInput(
      inputId = "shapreg_ct2",
      label = "Compare cell type 2:",
      selectize = TRUE,
      choices = c(""),
      selected = ""
    )
    
    
  ),
  mainPanel(
    plotOutput(outputId = "genomebrowserPlot")
  )
)

################################################################################
########### About page #########################################################
################################################################################
tab_about <- verbatimTextOutput("todo description of project here")



################################################################################
########### QC page ############################################################
################################################################################


################################################################################
########### Network for each gene ##############################################
################################################################################
#network<- read.csv('/home/mahogny/nando/cytoscape_edgeinfo.csv', header = T)
tab_geneinnet <- sidebarLayout(
  
  sidebarPanel(
    
    selectInput(
      inputId = "geneinnet_useclustering",
      label = "Clustering:",
      choices = names(available_clusters), 
      selected = names(available_clusters)[1]
    ),
    
    
    selectInput(
      inputId = "geneinnet_gene",
      label = "Gene:",
      selectize = TRUE,
      choices = c("CDKL5"), #shapregs_genelist,
      selected = "CDKL5"
    ),
    
    sliderInput(
      inputId = "geneinnet_pcutoff_direct", 
      label = "Log10 p cutoff (direct upstream):",
      min = -5, max = 0, value = -2,step = 0.1
    ),
    
    sliderInput(
      inputId = "geneinnet_pcutoff_hp", 
      label = "Log10 p cutoff (hitting probability):",
      min = -5, max = 0, value = -1,step = 0.1
    )
    
    #checkboxInput("geneinnet_showdonors", "Show donor averages (faster):", TRUE),
    
    
    
  ),
  mainPanel(
    h3("Direct upstream probabilities:"),
    plotOutput(outputId = "geneinnetNeighbourPlot"),
    
    h3("Top hitting probabilities:"),
    plotOutput(outputId = "geneinneHittingPlot"),
  )
)

################################################################################
########### Total page #########################################################
################################################################################


ui <- fluidPage(
  
  # App title 
  titlePanel("CD4 T cell network database"),
  
  
  tabsetPanel(type = "tabs",
    tabPanel("Gene in network", tab_geneinnet),
    tabPanel("Core Network", tab_network),
    tabPanel("Genome tracks", tab_genomebrowser),
    tabPanel("About", tab_about),
    tabPanel("QC", verbatimTextOutput("todo integrate sebastians code"))
  )
  
)



