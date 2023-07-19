ComputeHittingProbability <- function(object, ...){
  UseMethod(generic = 'ComputeHittingProbability', object = object)
}

ComputeSteadyState <- function(object, ...){
  UseMethod(generic = 'ComputeSteadyState', object = object)
}

ComputeNandoWalktrap <- function(object, ...){
  UseMethod(generic = 'ComputeNandoWalktrap', object = object)
}


NumSteps <- function(object, ...){
  UseMethod(generic = 'NumSteps', object = object)
}

TransitionMatrix <- function(object, ...){
  UseMethod(generic = 'TransitionMatrix', object = object)
}


