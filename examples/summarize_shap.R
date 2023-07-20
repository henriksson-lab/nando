if(FALSE){
  install.packages("doParallel")
  devtools::install_github('henriksson-lab/Nando')
}



library(doParallel)
registerDoParallel(parallel::detectCores())
library(Nando)

print("number of cores")
print(parallel::detectCores())


nando_dir <- "/home/mahogny/mystore/nando2"

SummarizeShapByClusterParallel(nando_dir)
