#' @export
ComputeHittingProbability <- function(object, ...){
  UseMethod(generic = 'ComputeHittingProbability', object = object)
}

#' @export
ComputeSteadyState <- function(object, ...){
  UseMethod(generic = 'ComputeSteadyState', object = object)
}

#' @export
ComputeNandoWalktrap <- function(object, ...){
  UseMethod(generic = 'ComputeNandoWalktrap', object = object)
}


#' @export
NumSteps <- function(object, ...){
  UseMethod(generic = 'NumSteps', object = object)
}

#' @export
TransitionMatrix <- function(object, ...){
  UseMethod(generic = 'TransitionMatrix', object = object)
}


#' @export
ComputeSimplifiedMatrix <- function(object, ...){
  UseMethod(generic = 'ComputeSimplifiedMatrix', object = object)
}

