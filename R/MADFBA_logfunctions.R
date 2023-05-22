
# Log utility functions should be in the utility file
#' logSubstrateFH
#'
#' @param model 
#' @param biomass 
#' @param initConcentrations 
#'
#' @noRd
#'

logSubstrateFH <- function (model,
                            biomass,
                            initConcentrations
)
{
    if (MADFBA_SETTINGS("VERBOSE") < 2) return()
    maxLogFluxes <- ifelse(length(initConcentrations) > MADFBA_SETTINGS("MAXLOGFLUXES"),
                                  MADFBA_SETTINGS("MAXLOGFLUXES"),
                                  length(initConcentrations))
    substrate <- names(initConcentrations)[1:maxLogFluxes]
    h0 <-  '                     '
    h <-   '  Reacts |  biomass  '
    lb <-  '  lowbnd |    --     '
    c <-   sprintf("  %s | %10.8f", 'inicon', biomass)
    rsep <-'  ------ | ----------'
    h0 <- h <- lb <- c <- rsep <- ''
    for (i in 1:length(substrate))  {
        h0    <- sprintf("%s %s %10s", h0, '|', ' substrate')
        h     <- sprintf("%s %s %10s", h, '|', substrate[i])
        lbnd  <- ifelse(!is.logical(model@lowbnd[model@react_id %in% substrate[i]]), 
                                sprintf("%10.3f", model@lowbnd[model@react_id %in% substrate[i]]), 
                                '          ')
        lb    <- sprintf("%s %s %10s", lb, '|', lbnd)
        c     <- sprintf("%s %s %10.3f", c, '|', initConcentrations[i])
        rsep    <- sprintf("%s %s %s", rsep, '|', '----------')
        
    }  
    
    cat('\n\n')
    cat('                     ', h0, ' |', gsub(' substrate', '   flux   ', h0), '\n', sep = '')
    cat('  Reacts |  biomass  ', h, ' |', h, '\n', sep = '')
    cat('  lowbnd |    --     ',lb, ' |', lb, '\n', sep = '')
    cat(sprintf("  %s | %10.8f", 'inicon', biomass), c, ' |', c, '\n', sep = '')
    cat('  ------ | ----------', rsep, ' |', rsep, '\n', sep = '')
    
}


#' logFluxes
#'
#' @param biomass 
#' @param fluxes 
#' @param stepNo 
#'
#' @noRd
#' 
logFluxes <- function (biomass, fluxes, stepNo)
{
    if (MADFBA_SETTINGS("VERBOSE") <= 2) return()

    f <- sprintf("%10.8f", biomass)
    for (i in 1:length(fluxes))  {
        f     <- sprintf("%s %s %10.3f", f, '|', fluxes[i])
    }  
    
    cat(sprintf("%8s |", paste("[", stepNo, "]", sep = "")), f, '\n')

}

#' logConcentrations
#'
#' @param biomass 
#' @param concentrations 
#' @param stepNo 
#'
#' @noRd
#' 
logConcentrations <- function (biomass, concentrations, stepNo){
    
    if (MADFBA_SETTINGS("VERBOSE") <= 2) return()
    
    c <- sprintf("%10.8f", biomass)
    for (i in 1:length(concentrations))  {
        c     <- sprintf("%s %s %10.3f", c, '|', concentrations[i])
    }
    cat(sprintf("%8s |", paste("[", stepNo, "]", sep = "")), c, '\n')
    
}

#' logConcentFluxes
#'
#' @param biomass 
#' @param concentrations 
#' @param fluxes 
#' @param stepNo 
#'
#' @noRd
logConcentFluxes <- function (biomass, concentrations, fluxes, stepNo){
    
    if (MADFBA_SETTINGS("VERBOSE") < 2) return()
    maxLogFluxes <- ifelse(length(concentrations) > MADFBA_SETTINGS("MAXLOGFLUXES"),
                                  MADFBA_SETTINGS("MAXLOGFLUXES"),
                                  length(concentrations))

    c <- sprintf("%10.8f", biomass)
    f <- '|'
    for (i in 1:maxLogFluxes)  {
        c     <- sprintf("%s %s %10.3f", c, '|', concentrations[i])
        f     <- sprintf("%s %s %10.3f", f, '|', fluxes[i])
    }
    cat(sprintf("%8s |", paste("[", stepNo, "]", sep = "")), c, f, '\n', sep = ' ')

}

getPrefix <- function(prefix, depth = -2){
    if (prefix=='') prefix <- deparse(sys.call(depth)[[1]])
    return(prefix)
}
logDbg <- function(prefix, msg, ...){
    logMADFBA(level=9, prefix = prefix, 'DEBUG: ', msg, ...)
}
logInfo <- function (prefix, msg, ...){
    logMADFBA(level=3, prefix = prefix, 'INFO: ', msg, ...)
}
logInfo2 <- function (msg, ...){
    #if (prefix=='') 
    prefix <- getPrefix('')
    logMADFBA(level=2, prefix = prefix, 'INFO: ', msg, ...)
}
logWarn <- function (prefix, msg, ...){
    logMADFBA(level=1, prefix = prefix, 'WARNING: ', msg, ...)
}
logError <- function (prefix, msg, ...){
    logMADFBA(level=0, prefix = prefix, 'ERROR: ', msg, ...)
}

logMADFBA <- function (level=3, prefix='MADFBA', type ='INFO: ', msg, ...){
    if (MADFBA_SETTINGS("VERBOSE") >= level) {
        msg <- gettext(paste(c(prefix, type, msg, ...), collapse = " "))
        cat("# ", msg, "\n", sep = "")
    }
}

