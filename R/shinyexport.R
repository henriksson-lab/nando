

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
ExportShinyNando.SteadyStates <- function(nando_dir, nandonets){
  ssmat <- SteadyStateMatrix(nandonets)
  saveRDS(ssmat, file.path(nando_dir,"ss.RDS")) #Because it is small enough, simplest to just read all into memory
}


#' Export HP.SS data available to ShinyNando
#' 
#' @param nandonets A ListOfNandoNetwork
#' @return Nothing
ExportShinyNando.HP.SS <- function(nando_dir, nandonets){
  hpss <- ComputeHittingProbabilityFromSS(nandonets)
  saveRDS(hpss, file.path(nando_dir,"hpss.RDS")) #Because it is small enough, simplest to just read all into memory
}



#' Import SS data available to ShinyNando
#' 
#' @param nandonets A ListOfNandoNetwork
#' @return Nothing
ImportShinyNando.SteadyStates <- function(nando_dir){
  readRDS(file.path(nando_dir,"ss.RDS"))
}


#' Import SS data available to ShinyNando
#' 
#' @param nandonets A ListOfNandoNetwork
#' @return Nothing
ImportShinyNando.HP.SS <- function(nando_dir){
  readRDS(file.path(nando_dir,"hpss.RDS"))
}





#' Export all data available to ShinyNando
#' 
#' @param nandonets A ListOfNandoNetwork
#' @return Nothing
#' 
#' @export
ExportShinyNando <- function(nando_dir, nandonets){
  ExportShinyNando.HP(nando_dir, nandonets)
  ExportShinyNando.SteadyStates(nando_dir, nandonets)
  ExportShinyNando.HP.SS(nando_dir, nandonets)
  ExportShapPerRegionSqlite(nando_dir)
  ExportShapSqlite(nando_dir)
  ExportGeneCategories(nando_dir, nandonets)
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


#' Read probability matrices for one gene
ReadProbabilityMatricesH5ForGene <- function(matrixfile, netnames, genename){  #this is hp, so gene from to all other genes?

  #matrixfile <- h5hp
  #netnames <- c("donor1","donor2")  #issue. not renamed in h5 yet. or store again?
  
  #genenames <- get_hp_genelist(donoraveraged)
  #netnames <- get_hp_netlist(donoraveraged)
  
  allhp <- NULL  
  for(this_net in netnames){
    
    net_index <- which(matrixfile@list_names==this_net)
    if(length(net_index)==0){
      stop(paste("Failed to find matrix", matrixfile@fname, this_net))
    }
    
    gene_index <- which(matrixfile@list_rows[[net_index]]==genename)
    if(length(gene_index)==0){
      stop(paste("Failed to find gene in matrix", matrixfile@fname, this_net, genename))
    }
    
    genenames <- matrixfile@list_cols[[net_index]]
    hp <- rhdf5::h5read(matrixfile@fname, sprintf("mat_%s",net_index), index=list(gene_index,1:length(genenames)))
    allhp <- rbind(
      allhp,
      data.frame(
        category=this_net,
        gene=genenames,
        p=as.double(hp)
      )
    )
  }
  allhp
}


if(FALSE){
  
  h5hp <- PrepareProbabilityMatrixH5(file.path(nando_dir,"hp.h5"))
  h5ss <- PrepareProbabilityMatrixH5(file.path(nando_dir,"ss.h5"))
  
  ReadProbabilityMatricesH5ForGene(h5hp, c("donor1 ALL","donor2 ALL"), "CD8A")
  
}



#' Store SHAPs per network, gene and region
#' 
#' @import RSQLite
#' @return Nothing
#' 
ExportShapSqlite <- function(nando_dir){
  sqlite_file_ct_shapregs <- file.path(nando_dir, "shaps.sqlite")
  
  ##Store in an sqlite file such that we can load only one gene of interest in the future
  if(file.exists(sqlite_file_ct_shapregs)){
    file.remove(sqlite_file_ct_shapregs)
  }
  
  con <- dbConnect(RSQLite::SQLite(), dbname = sqlite_file_ct_shapregs)
  dir_shapsummary <- file.path(nando_dir,"summarizedshaps")
  for(f in list.files(dir_shapsummary)){
    if(str_ends(f,"RDS")){
      print(f)
      #Read some genes
      shapregs <- readRDS(file.path(dir_shapsummary,f))
      
      shapregs$category <- str_replace_all(shapregs$category,"donordonor","donor") #hack, to remove
      
      
      #Normalize to get p-values
      shapregs <- merge(shapregs, sqldf::sqldf("select category, gene, sum(gain) as sumgain from shapregs group by gene,category"))
      shapregs$p <- shapregs$gain/shapregs$sumgain

      #Save in SQL
      shapregs <- shapregs[,c("tf","category","gene","p")]
      dbWriteTable(con, "shaps", as.data.frame(shapregs), append=TRUE)
    }
  }
  dbExecute(con, "CREATE INDEX shaps_gene ON shaps (gene);")
  
  #Store list of genes for fast query later  
  sql_result <- dbSendQuery(con, "SELECT distinct gene FROM shaps")
  genelist <- dbFetch(sql_result)
  dbClearResult(sql_result)
  dbWriteTable(con, "genelist", genelist)
  
  dbDisconnect(con)
}




#' Store SHAPs per network, gene and region
#' 
#' @import RSQLite
#' @return Nothing
#' 
ExportShapPerRegionSqlite <- function(nando_dir){
  sqlite_file_ct_shapregs <- file.path(nando_dir, "shapregs.sqlite")
  
  ##Store in an sqlite file such that we can load only one gene of interest in the future
  if(file.exists(sqlite_file_ct_shapregs)){
    file.remove(sqlite_file_ct_shapregs)
  }
  
  con <- dbConnect(RSQLite::SQLite(), dbname = sqlite_file_ct_shapregs)
  dir_shapregsummary <- file.path(nando_dir,"summarizedshapregs")
  for(f in list.files(dir_shapregsummary)){
    if(str_ends(f,"RDS")){
      print(f)
      #Read some genes
      shapregs <- readRDS(file.path(dir_shapregsummary,f))
      
      shapregs$category <- str_replace_all(shapregs$category,"donordonor","donor") #hack, to remove
      
      #Normalize to get p-values
      shapregs <- merge(shapregs, sqldf::sqldf("select category, gene, sum(gain) as sumgain from shapregs group by gene,category"))
      shapregs$p <- shapregs$gain/shapregs$sumgain
      #shapregs <- shapregs[,c("tf","category","gene","p","region")]
      
      #Expand coordinates into separate variables
      pos <- as.data.frame(str_split_fixed(shapregs$region,"_",3))
      colnames(pos) <- c("chrom","start","end")
      pos$start <- as.integer(pos$start)
      pos$end <- as.integer(pos$end)
      shapregs <- cbind(shapregs[,c("tf","p","category","gene")], pos)
      
      #Save in SQL
      shapregs <- shapregs[,c("tf","category","gene","p","chrom","start","end")]
      dbWriteTable(con, "shapregs", as.data.frame(shapregs), append=TRUE)
    }
  }
  dbExecute(con, "CREATE INDEX shapregs_gene ON shapregs (gene);")
  
  #Store list of genes for fast query later  
  sql_result <- dbSendQuery(con, "SELECT distinct gene FROM shapregs")
  genelist <- dbFetch(sql_result)
  dbClearResult(sql_result)
  dbWriteTable(con, "genelist", genelist)
  
  dbDisconnect(con)
}


#' Store category of each gene
#' @param nando_dir Nando directory
#' @param nandonets ListOfNandoNetwork
#' @return Nothing
ExportGeneCategories <- function(nando_dir, nandonets){
  df <- GeneCategoriesDF.ListOfNandoNetwork(nandonets)
  saveRDS(df, file.path(nando_dir,"cellcategory.RDS"))
}

#' Read category of each gene
#' @param nando_dir Nando directory
#' @return Data frame with categories per net
ImportGeneCategories <- function(nando_dir){
  readRDS(file.path(nando_dir,"cellcategory.RDS"))
}









#' Get all the upstream SHAP info for one gene, including region
#' 
#' 
#' 
ImportShapregsForGene <- function(nando_dir, genename){
  sqlite_file_ct_shapregs <- file.path(nando_dir, "shapregs.sqlite")
  con <- dbConnect(RSQLite::SQLite(), dbname = sqlite_file_ct_shapregs)
  iris_result <- dbSendQuery(con, "SELECT * FROM shapregs where [gene] = ?")
  dbBind(iris_result, list(genename))
  shapregs <- dbFetch(iris_result)
  dbClearResult(iris_result)
  dbDisconnect(con)
  shapregs
}

#' Get which genes are in this database
#' 
#' 
#' 
ImportShapregsGenelist <- function(nando_dir){
  sqlite_file_ct_shapregs <- file.path(nando_dir, "shapregs.sqlite")
  con <- dbConnect(RSQLite::SQLite(), dbname = sqlite_file_ct_shapregs)
  sql_result <- dbSendQuery(con, "SELECT * FROM genelist")
  genelist <- dbFetch(sql_result)
  dbClearResult(sql_result)
  dbDisconnect(con)
  genelist$gene
}



#' Get all the upstream SHAP info for one gene, including region
#' 
#' 
#' 
ImportShapsForGene <- function(nando_dir, genename){
  sqlite_file_ct_shapregs <- file.path(nando_dir, "shaps.sqlite")
  con <- dbConnect(RSQLite::SQLite(), dbname = sqlite_file_ct_shapregs)
  iris_result <- dbSendQuery(con, "SELECT * FROM shaps where [gene] = ?")
  dbBind(iris_result, list(genename))
  shapregs <- dbFetch(iris_result)
  dbClearResult(iris_result)
  dbDisconnect(con)
  shapregs
}

