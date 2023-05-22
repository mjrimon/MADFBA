#' Find biomass reaction in the model
#' @description Helper function to search for biomass reaction in available reactions of a model
#' 
#' @param model Object of class sybil::modelorg containging the genome sclae metabolic model
#' @param biomassRxn Biomass reaction name
#' @param biomassPartName Biomass partial reaction name
#' @return Reaction id for biomass reaction
#' @noRd
findBiomassReact <- function(model, biomassRxn= '', biomassPartName=c("biomass")){
  # we can use sybil::react_id(model) or model@react_id which one is preferred ???
  warn<-info<-''
      if (biomassRxn != '') {
        # find biomass reaction index
        # if the reaction is not found, we'll have an error
        biomassIdx <- which(sybil::react_id(model) %in% biomassRxn)
        if (length(biomassIdx) == 0) {
          # the reaction was not found
          cat(warn, 'The specified biomass reaction', biomassRxn, '\n')
          cat(warn, 'has not been found in the model: Ignoring it.\n')
          biomassRxn <- ''
        }
      }

      # if the biomass reaction has not been specified, we'll first look for
      # a reaction called 'biomass' and, if none is found, we'll use the
      # first objective function.
      #
      # NOTE that there may be more than one objective function: we'll
      # settle for the first, but this may not be the correct one!
      # e.g. if obj <- c("SECRETION", "BIOMASS"), then growth will be
      # tied to secretion, not to biomass!!!
      if (biomassRxn == '' || missing(biomassRxn)) {
        # Try to find the real Biomass function before accepting to use only
        # the objective function (which might even be multiple objectives)
        #
        # THIS IS VERY, VERY RISKY
        #       It might be that the biomass reaction is called something else
        #       and the one named *biomass* is not the real Biomass reaction
        #       or that there is more than one reaction called '*biomass*'
        #   and the proper one is not the first. Hopefully we will get it
        #       right
        #
        biomassIdx <- grep(biomassPartName, sybil::react_id(model), ignore.case=TRUE)[1]
        # if found use it
        if (length(biomassIdx) != 0) {
            biomassRxn <- sybil::react_id(model)[biomassIdx][1]
            cat(warn, "Setting Biomass reaction to the first reaction containing\n")
            cat('       BIOMASS in its name:', biomassRxn, '(', biomassIdx, ')\n')
        }    
        # if not, revert to objective function
        else {
            cat(warn, "Setting Biomass reaction to default:\n")
            cat(warn, " FIRST objective function!!!\n")
            biomassIdx <- which(sybil::obj_coef(model) != 0)[1]
            biomassRxn <- sybil::react_id(model)[biomassIdx][1]
        }
        cat(warn, 'IF THIS IS WRONG, USE PARAMETER "biomassRxn" TO SET IT\n\n')
      }
      cat(info, 'Biomass reaction:', biomassRxn, "(",biomassIdx, ')\n\n')
      return(biomassIdx)
}




#' Check status of lp optimization with the used solver
#' @description Check status of linear optimization depending on solver
#' 
#' @param opt_sol Object of class sybil::optsol containging the result of an optimization
#' @param solver Name of the solver used in the optimization
#' @return status OK?
#' @noRd
getLPstat <- function(opt_sol, solver = sybil::SYBIL_SETTINGS("SOLVER")){
  return(opt_sol$stat==solverStatOK(solver))
}

#' get value of ok status  for the solver
#' @description get ok status for the solver
#' 
#' @param solver Name of the solver 
#' @return Value of OK status for the solver
#' @noRd
solverStatOK <- function(solver = sybil::SYBIL_SETTINGS("SOLVER")){
  # values from function checkSolStat, file optObj_basicfunc.R, package sybil 
  switch(solver,
         glpkAPI =     {solve_ok <- 5},
         cplexAPI =    {solve_ok <- 1}, # tambien son validos 101 y 102, desconozco la diferencia
         clpAPI =      {solve_ok <- 0},
         lpSolveAPI =  {solve_ok <- 0},
         stop("Solver not suported!"))
  return(solve_ok)
}


#' showFH
#'
#' @param mod 
#' @param logRxns 
#' @param title 
#'
#' @noRd

showFH <- function(mod, logRxns, title){
  cat('\n\n', '  --------- ', title, ' ---------', '\n\n')
  
  h <-  '  Reacts |  biomass  '
  lb <- '  lowbnd |    --     '
  rsep<-'  ------ | ----------'
  for (i in 1:length(logRxns))  {
    h     <- sprintf("%s %s %10s", h, '|', logRxns[i])
    lb    <- sprintf("%s %s %10.3f", lb, '|', mod@lowbnd[mod@react_id %in% logRxns[i]])
    rsep    <- sprintf("%s %s %s", rsep, '|', '----------')
  }
  cat( h, '\n')
  cat(lb, '\n')
  cat(rsep, '\n')
}

#' showConcentrations
#'
#' @param mod 
#' @param opt_sol 
#' @param logRxns 
#'
#' @noRd
#' 
showConcentrations <- function (mod, opt_sol, logRxns){
  showFH(mod, logRxns, 'CONCENTRATIONS')
  
  for (i in 1:dim(opt_sol@concentrationMatrix)[2])  {
    c <- sprintf("%10.8f", opt_sol@biomassVec[i])
    for (j in 1:length(logRxns))  {
      c     <- sprintf("%s %s %10.3f", c, '|', opt_sol@concentrationMatrix[opt_sol@excRxnNames %in% logRxns[j] ,i])
    }
    cat(sprintf("%8s |", paste("[", i-1, "]", sep = "")), c, '\n')
  }
  
}

#' showFluxes
#'
#' @param mod 
#' @param opt_sol 
#' @param logRxns 
#'
#' @noRd
#' 
showFluxes <- function (mod, opt_sol, logRxns){
  showFH(mod, logRxns, 'FLUXES')
  
  if (is.null(opt_sol@all_fluxes)) return()
  for (i in 1:dim(opt_sol@all_fluxes)[2])  {
    f <- sprintf("%10.8f", opt_sol@biomassVec[i])
    for (j in 1:length(logRxns))  {
      f     <- sprintf("%s %s %10.3f", f, '|', opt_sol@all_fluxes[mod@react_id %in% logRxns[j] ,i])
    }
    cat(sprintf("%8s |", paste("[", i-1, "]", sep = "")), f, '\n')
  }
  
}

addAttribute <- function(mod, attr, val) {
    for (i in 1:length(mod)){
        prevAttr <- names(mod[[i]])
        mod[[i]] <- append(mod[[i]], val[i])
        names(mod[[i]]) <- c(prevAttr, attr)
    }
    return(mod)
}
