

################################################################################
########################### Export #############################################
################################################################################


#clusterings.csv  will be used to let users try different sets. especially donor average vs nonaverage



#' Store list of matrices 
#' 
#' @return Nothing
StoreProbabilityMatricesH5 <- function(fname, list_matrix){
  
  if(file.exists(fname)){
    file.remove(fname)
  }
  
  #Create file, write metadata
  rhdf5::h5createFile(fname)
  rhdf5::h5write(names(list_matrix), fname,"listnets")

  #Store "total" row and column, genes present anywhere
  total_col <- sort(unique(unlist(lapply(list_matrix, colnames))))
  total_row <- sort(unique(unlist(lapply(list_matrix, rownames))))
  rhdf5::h5write(total_col, fname, "totalcols")
  rhdf5::h5write(total_row, fname, "totalrows")
  
  #Store matrices
  for(i in 1:length(list_matrix)){
    print(paste(i,names(list_matrix)[i]))
    onemat <- as.matrix(list_matrix[[i]])  #can we avoid this?
    
    #The matrix body
    matname <- sprintf("mat_%s",i)
    rhdf5::h5createDataset(file = fname, dataset = matname, 
                           dims = dim(onemat), storage.mode = "double", 
                           chunk = c(2,ncol(onemat)), level = 6) ######### This makes queries for (i->j) very fast for a given i
    rhdf5::h5write(onemat, fname,matname)
    
    #The columns and rows
    rhdf5::h5write(colnames(onemat), fname, sprintf("cols_%s",i))
    rhdf5::h5write(rownames(onemat), fname, sprintf("rows_%s",i))
  }
}




#' Export HP data available to ShinyNando
#' 
#' @param nandonets A ListOfNandoNetwork
#' @return Nothing
ExportShinyNando.HP <- function(nando_dir, nandonets){
  list_mat <- lapply(nandonets@nets, function(net) net@hp)
  StoreProbabilityMatricesH5(file.path(nando_dir,"hp.h5"), list_mat)
}


#' Export SS data available to ShinyNando
#' 
#' @param nandonets A ListOfNandoNetwork
#' @return Nothing
ExportShinyNando.SteadyStates <- function(nandor_dir, nandonets){
  ssmat <- SteadyStateMatrix(nandonets)
  saveRDS(ssmat, file.path(nando_dir,"ss.RDS")) #Because it is small enough, best to just read into memory
}



#' Export all data available to ShinyNando
#' 
#' @param nandonets A ListOfNandoNetwork
#' @return Nothing
#' 
#' @export
ExportShinyNando <- function(nandor_dir, nandonets){
  ExportShinyNando.HP(nandor_dir, nandonets)
  ExportShinyNando.SteadyStates(nandor_dir, nandonets)
}



################################################################################
########################### Reading ############################################
################################################################################




#' Prepare a set of stored matrices for reading.
#' Read essentials to speed up later access
#' 
#' @param fname Name of the file
#' 
PrepareProbabilityMatrixH5 <- function(fname){

  list_names <- rhdf5::h5read(fname,"listnets")
  totalcols <- rhdf5::h5read(fname,"totalcols")
  totalrows <- rhdf5::h5read(fname,"totalrows")
  
  numnet <- length(list_names)
  list_cols <- list()
  list_rows <- list()
  for(i in 1:numnet){
    list_cols[[i]] <- rhdf5::h5read(fname,sprintf("cols_%s",i))
    list_rows[[i]] <- rhdf5::h5read(fname,sprintf("rows_%s",i))
  }
  new("ProbabilityMatrixH5", 
      fname=fname,
      totalcols=totalcols,
      totalrows=totalrows,
      list_cols=list_cols,
      list_rows=list_rows,
      list_names=list_names
  )  
}

h5hp <- PrepareProbabilityMatrixH5(file.path(nando_dir,"hp.h5"))
h5ss <- PrepareProbabilityMatrixH5(file.path(nando_dir,"ss.h5"))



StoredProbabilityMatrixH5

StoredProbabilityMatrixH5


