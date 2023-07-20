#' @importFrom methods setClass
#' @importClassesFrom Matrix sparseMatrix
#' @importClassesFrom SeuratObject Seurat
NULL


# These let us specify NULL values initially
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("data.frameOrNULL",members=c("data.frame", "NULL"))
setClassUnion("characterOrNULL",members=c("character", "NULL"))
setClassUnion("sparseMatrixOrNULL",members=c("dgCMatrix", "NULL"))


################################################################################
################### NandoNetwork ###############################################
################################################################################


#' An S4 class to represent a list of networks
#'
#' @slot nets All networks
#' 
#' @export
#' 
ListOfNandoNetwork <- setClass("ListOfNandoNetwork", slots=list(
  nets="list"
  )
)

#' An S4 class to represent the model of a network.
#'
#' Note that if a gene is never listed in the TF column, it will be considered a nonTF.
#'
#' @slot list_all_gene All genes (TFs and nonTFs)
#' @slot list_tfs All TFs every mentioned
#' @slot list_irreducible TFs that are members of the irreducible net
#' @slot list_nontfs All nonTFs
#' 
#' @export
#' 
NandoNetwork <- setClass("NandoNetwork", slots=list(
  shap="data.frameOrNULL",
  
  list_all_gene="characterOrNULL",
  list_tfs="characterOrNULL",
  list_irreducible="characterOrNULL",
  list_nontfs="characterOrNULL",

  ss="numericOrNULL",
  hp="sparseMatrixOrNULL", 
  hp2ss="numericOrNULL"
  )
)


#' An S4 class to represent a probability matrix file. 
#' Provides speedy out-of-memory access.
#' 
#' @slot fname Name of the file
#' @slot totalcols All the column names ever mentioned
#' @slot totalrows All the row names ever mentioned
#' @slot list_cols Column names of each matrix
#' @slot list_rows Row names of each matrix
#' @slot list_names Names of each network
#' 
#' @export
#' 
setClass("ProbabilityMatrixH5", slots=list(
  
  fname="character",
  
  totalcols="array",
  totalrows="array",
  
  list_cols="list",
  list_rows="list",
  list_names="array"
  )
)


print.NandoNetwork <- function(nandonet){
  cat("NandoNetwork\n")
  cat(
    "#gene:", length(nandonet@list_all_gene),
    "  #TFs:", length(nandonet@list_tfs),
    "  #I.TFs:", length(nandonet@list_irreducible),
    "  #nonTFs:", length(nandonet@list_nontfs)
  )
}
setMethod('show', 'NandoNetwork', function(object) print(object))



print.ListOfNandoNetwork <- function(nandonets){
  cat("List of NandoNetworks\n")
  cat("nets: ",names(nandonets@nets))
}
setMethod('show', 'ListOfNandoNetwork', function(object) print(object))



names.ListOfNandoNetwork <- function(nandonets){
  names(nandonets@nets)
}





################################################################################
################### Walktrap ###################################################
################################################################################


#' An S4 class to represent a list of walktrap analyses
#'
#' @slot nets All walktraps
#' 
#' @export
#' 
setClass("ListOfNandoWalktrap", slots=list(
  nets="list"
  )
)

#' An S4 class to represent a walktrap analysis.
#' It holds P^t and entropy for each step
#' 
#' @slot mat All P stored concatenated
#' @slot meta For each row of P, what it corresponds to
#' @slot entropy Entropy for every step
#' 
#' @export
#' 
setClass("NandoWalktrap", slots=list(
  mat="sparseMatrix",
  meta="data.frame",
  entropy="matrix"
  )
)


print.ListOfNandoWalktrap <- function(nandonets){
  cat("List of NandoWalktrap\n")
  cat("List of nets: ",names(nandonets@nets),"\n")
  cat("Num steps: ",NumSteps(nandonets))
}
setMethod('show', 'ListOfNandoWalktrap', function(object) print(object))


names.ListOfNandoWalktrap <- function(wts){
  names(wts@nets)
}


NumSteps.ListOfNandoWalktrap <- function(wts){
  NumSteps(wts@nets[[1]])
}

NumSteps.NandoWalktrap <- function(wts){
  max(wts@meta$step)
}

