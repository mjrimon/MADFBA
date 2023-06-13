#' Create a model list. 
#'
#' Create a model list. Each model is also a list of named elements with all variables needed to process 
#' the model. The model is an object of class sybil::modelorg
#' The params can be a unique value meaning all models (1:n) in the models params have the same value for 
#' example of initial biomass or a vector with a length equal to the number of models containing, in the
#' same order that the model, their respective value 
#'
#' @param models An object or a list of objects of sybil::modelorg
#' @param biomass Initial biomass of the models. One value -> same biomass for all models or a vector with each model biomass 
#' @param biomassRxn Biomass reaction for the models. One value or a list of values.
#' @param exclUptakeRxns Excluded reactions for the models. One value or a list of values.
#' @param deathrate Death rate of the models. One value or a list of values.

#'
#' @return A list of named lists with the models data
#'
#' @noRd
makeOrganisms <- function(models, biomass, biomassRxn, exclUptakeRxns, dynamicConstraints, deathrate){
    checkNumParams <- function(param, num, msg){
        logInfo(deparse(sys.call(0)[[1]]), 'Checking param', msg)
        if (length(param) == 1) {
            param <- rep(param[1], num)
        }
        stopifnot(length(param) == num)
        return(param)
    }
    # sys.call(0)  this function name
    # sys.call(-1) caller function name
    #logPrefix <- deparse(sys.call(0)[[1]])
    logInfo2('Checking models. All suplied models must be modelorg ...')
 
    # check models
    # We need, at least, one model. All models must be modelorg class
    if (!is.list(models)) models <- c(models)
    # Check that all models are modelorg class
    lapply(models, function(model){ stopifnot(is(model, "modelorg"))})
    models.num <- length(models)
   
    logInfo2('Models checked ok. Total models:', models.num)


    # check deathrate
    deathrate <- checkNumParams(deathrate, models.num, 'deathrate')

    # check biomass
    biomass <- checkNumParams(biomass, models.num, 'biomass')
 
    # check biomassRxn
    biomassRxn <- checkNumParams(biomassRxn, models.num, 'biomassRxn')
    
    # get index of biomass reaction of each model, 
    # probably findBiomassReact function is too verbose.
    biomassIdx <- mapply(findBiomassReact, models, biomassRxn)# SIMPLIFY = FALSE)

    # update biomassRxn to found biomass react
    biomassRxn <- mapply(function(model, bioIdx) {return(model@react_id[bioIdx])}, models, biomassIdx)
   
    # set default value for exclUptakeRxns. No more missing check needed
    # if (missing(exclUptakeRxns)){
    #     exclUptakeRxns <- c('EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)');
    # }
    if (!is.list(exclUptakeRxns)) exclUptakeRxns <- list(exclUptakeRxns)
    exclUptakeRxns <- checkNumParams(exclUptakeRxns, models.num, 'exclUptakeRxns')

    # check dynamicConstraints
    if (is.null(dynamicConstraints)) 
        dynamicConstraints <- lapply(models, function(model) {return(NULL)})
    else {
       if (is.data.frame(dynamicConstraints)) dynamicConstraints <- list(dynamicConstraints)
       dynamicConstraints <- checkNumParams(dynamicConstraints, models.num, 'dynamicConstraints')
    }

    # Get exchange reactions of each model
    excReact <- lapply(models, function(model) {return(sybil::findExchReact(model))})

    # keep track of all uptakes and secretions. All values starting in 0
    # negative values means uptakes and positive represent secretions
    exUptakeMatrix <- lapply(excReact, function(ex) {return(numeric(ex@react_num))})


    all_fluxes <- lapply(models, function(model) {return(numeric(model@react_num))})
    all_stat <- lapply(models, function(stat) {return(solverStatOK())})

    lpmod <- lapply(models, function(model) {return(sybil::sysBiolAlg(model, algorithm = "fba"))})

    organisims <- lapply(models, function(model) {return(list(model=model))})
    organisims <- addAttribute(organisims, 'deathrate', deathrate)
    organisims <- addAttribute(organisims, 'biomass', biomass)
    organisims <- addAttribute(organisims, 'biomassRxn', biomassRxn)
    organisims <- addAttribute(organisims, 'biomassIdx', biomassIdx)
    organisims <- addAttribute(organisims, 'excReact', excReact)
    organisims <- addAttribute(organisims, 'exclUptakeRxns', exclUptakeRxns)
    organisims <- addAttribute(organisims, 'dynamicConstraints', dynamicConstraints)
    organisims <- addAttribute(organisims, 'biomassVec', biomass)
    organisims <- addAttribute(organisims, 'exUptakeMatrix', exUptakeMatrix)
    organisims <- addAttribute(organisims, 'all_fluxes', all_fluxes)
    organisims <- addAttribute(organisims, 'all_stat', all_stat)
    organisims <- addAttribute(organisims, 'lpmod', lpmod)

    return(organisims)
}

addOrgsToMedium <- function(models, initConcentrations, contype) {

    # Exclude reactions with concentrations that will not be changed 
    exclUptakeRxnsUpd <- lapply(models, function(mod, subs) {
                                    exclUpdated <- mod$excReact@react_id %in% mod$exclUptakeRxns |
                                                     (mod$excReact@uptake & !mod$excReact@react_id %in% subs)
                                    return(mod$excReact[exclUpdated]@react_id)  
                                    }
                            ,names(initConcentrations)
                            )


    if (contype == 'legacy') {
        concentrations <- lapply(models, getconcentrations, initConcentrations)
    }else{
        concentrations <- lapply(models, getconcentrations, c()) #c('-1'= -1))
    }
    # We don't know right now all environment metabolites because it could be changes during simulations
    # therefore concentration matrix is initializaed to all exchange reactions with 0 value
    concentrationMatrix <- mapply(function(mod, concen, exc) {return(list(concen[mod$excReact@react_id]))}, models, concentrations)#, SIMPLIFY = FALSE);

    organisims <- addAttribute(models, 'exclUptakeRxnsUpd', exclUptakeRxnsUpd)
    organisims <- addAttribute(organisims, 'concentrations', concentrations)
    organisims <- addAttribute(organisims, 'concentrationMatrix', concentrationMatrix)
    organisims <- addAttribute(organisims, 'contype', contype)

    # Add indexes for later use
    excReactInd <- mapply(function(mod, exclUptakeRxnsUpd) {
                                    return(!mod$excReact@react_id %in% exclUptakeRxnsUpd)}
                        , models, exclUptakeRxnsUpd, SIMPLIFY = FALSE)

    excReactIdx <- mapply(function(mod, excReactInd) {
                                    return(mod$excReact[excReactInd]@react_pos)}
                        , models, excReactInd, SIMPLIFY = FALSE)

    excReactNames <- mapply(function(mod, excReactInd) {
                                    return(mod$excReact[excReactInd]@react_id)}
                        , models, excReactInd, SIMPLIFY = FALSE)
    
    substrateRxnsIdx <- lapply(models, function(mod, substrateRxns) {
                                    return(match(substrateRxns, mod$model@react_id))}
                        , names(initConcentrations))

    organisims <- addAttribute(organisims, 'excReactIdx', excReactIdx)
    organisims <- addAttribute(organisims, 'excReactNames', excReactNames)
    organisims <- addAttribute(organisims, 'substrateRxnsIdx', substrateRxnsIdx)

    return(organisims)
}


#' Get result data from a model simulation
#'
#' gets the simmulation result from one model like an optsol_dynamicFBA object or like a list
#'
#' @param mod List with one model data
#' @param timeVec  Simulation Time vector
#' @param fld fluxDistribution?
#' @param retOptSol optsol_dynamicFBA?
#' @param contype Concentrations summary type
#'
#' @return optsol_dynamicFBA object or named list lista, as retOptSol param, with simulation model result
#'
#' @noRd

getResult <- function(mod, timeVec, fld, retOptSol, contype) {
        
    if (contype=='legacy')
        concentmatrix <- mod$concentrationMatrix
    else
        concentmatrix <- mod$exUptakeMatrix
    
    # add names to the concentrationMatrix
    row.names(concentmatrix) <- mod$excReact@react_id;
    row.names(mod$all_fluxes) <- mod$model@react_id;

    # With multimodel, not all organism finish simulation at same time.
    # We only have valid data if state is OK, so filter timeVec
    time <- timeVec[mod$all_stat==solverStatOK()]

    # If the organism has defined 'deathrate' and is dying at some simulations steps
    # we have not ok solver state but we are applying death so filter timeVec
    if (length(mod$biomassVec) > length(time)) time <- timeVec[1:length(mod$biomassVec)]
    timeVec <- time

    stepNo <- length(timeVec)
    colnames(concentmatrix) <- timeVec
    # if we save all step fluxes
    if (fld) {
        # If finish by by extintion, with deathrate, we have one more column
        if (length(timeVec) == ncol(mod$all_fluxes))
            colnames(mod$all_fluxes) <- timeVec
        else
            colnames(mod$all_fluxes) <- timeVec[-1]
    }
    
    ## Preparing OUTPUT

    if (isTRUE(retOptSol)) {
        # return an optSol object
        if(is.null(mod$all_fluxes)) mod$all_fluxes=as.matrix(NA);

        return (MADFBA::optsol_dynamicFBA(solver = mod$lpmod@problem@solver,
                    method = mod$lpmod@problem@method,
                    nprob  = stepNo,
                    ncols  = mod$model@react_num,
                    nrows  = mod$model@met_num,
                    fld    = fld,
                    all_fluxes = mod$all_fluxes,
                    concmat= concentmatrix, 
                    exRxn  = rownames(concentmatrix), 
                    tmVec  = timeVec,  
                    bmVec  = mod$biomassVec,
                    contype = contype
                    )
            )
    }else{ # return a list instead of an optSol object
        return(optsol <- list(nprob  = stepNo,
                  ncols  = mod$model@react_num,
                  nrows  = mod$model@met_num,
                  all_fluxes = mod$all_fluxes,
                  all_stat= mod$all_stat,
                  concentrationMatrix=concentmatrix, 
                  excRxnNames=rownames(concentmatrix), 
                  timeVec= timeVec,  
                  biomassVec=mod$biomassVec,
                  mod = mod
                  )
           )
    }
}


#' In multimodel gets summary results of all models
#'
#' Gets the simulation summary results of all models like an optsol_dynamicFBA object or like a list
#'
#' @param mod List with one model data
#' @param mediumMatrix Environment concentration matrix of substrates in every simulation step
#' @param mediumBiomassVec Environment biomass vector in every simulation step
#' @param timeVec  Simulation Time vector
#' @param retOptSol optsol_dynamicFBA?
#'
#' @return optsol_dynamicFBA object or named list lista, as retOptSol param, with simulation model result
#'
#' @noRd

getResultSummary <- function(mods, mediumMatrix, mediumBiomassVec, timeVec, retOptSol) {

    stepNo <- length(timeVec)

    colnames(mediumMatrix) <- timeVec
    
    ## Preparing OUTPUT
    cm <- mediumMatrix
    dimcm <- dim(cm)
    if (isTRUE(retOptSol)) {
        # return an optSol object
        all_fluxes <- as.matrix(NA)
        if((length(mods) == 1) & (!is.null(mods[[1]]$all_fluxes))) all_fluxes=mods[[1]]$all_fluxes

        return (MADFBA::optsol_dynamicFBA(solver = sybil::SYBIL_SETTINGS("SOLVER"),
                    method = sybil::SYBIL_SETTINGS("METHOD"),
                    nprob  = stepNo,
                    ncols  = dimcm[2],
                    nrows  = dimcm[1],
                    fld    = FALSE,
                    all_fluxes = all_fluxes,
                    concmat= cm, #concentrationMatrix,
                    exRxn  = rownames(cm), #concentrationMatrix),
                    tmVec  = timeVec,  
                    bmVec  = mediumBiomassVec
                    )
            )
    }else{ # return a list instead of an optSol object
        return(optsol <- list(nprob  = stepNo,
                  ncols  = dimcm[2],
                  nrows  = dimcm[1],
                  all_fluxes = as.matrix(NA), # mod$all_fluxes,
                  all_stat= NA,           # posibly add an environment state?
                  concentrationMatrix=cm, 
                  excRxnNames=rownames(cm), 
                  timeVec= timeVec,  
                  biomassVec=mediumBiomassVec,
                  mod = mods
                  )
           )
    }
}

#' Compare 2 sybil::modelorg objects
#'
#' Compare two sybil::modelorg objects checking if they are equals
#'
#' @param A sybil::modelorg object to compare
#' @param B sybil::modelorg object to compare
#'
#' @return True if both models are equal else False
#'
#' @noRd
isEqualModel <- function(A, B){
    lapply(c(A, B), function(model) {if (!is(model, "modelorg")) return(FALSE)})
    m <- A
    Achecks = c(m@react_num, m@react_id, m@met_id, m@lowbnd, m@uppbnd, m@S@x)
    m <- B
    Bchecks = c(m@react_num, m@react_id, m@met_id, m@lowbnd, m@uppbnd, m@S@x)
    if (length(Achecks) != length(Bchecks)) return(FALSE)
    if (any(Achecks != Bchecks)) return(FALSE)
    return(TRUE)      
}

#' Compare reactions bounds from 2 sybil::modelorg objects
#'
#' Compare reactions bounds from two sybil::modelorg objects checking if they are equals
#'
#' @param A sybil::modelorg object to compare
#' @param B sybil::modelorg object to compare
#' @param react_pos model reaction indexes to compare
#'
#' @return True if both bounds are equal else False
#'
#' @noRd
isEqualModelReactBounds <- function(A, B, react_pos){
    lapply(c(A, B), function(model) {if (!is(model, "modelorg")) return(FALSE)})
    if (A@react_num != B@react_num) return(FALSE)
    if (A@react_num < react_pos) return(FALSE)
    return(A@lowbnd[react_pos] == B@lowbnd[react_pos] &&
           A@uppbnd[react_pos] == B@uppbnd[react_pos])      
}


#' Get initial environment concentrations for the model
#'
#' gets a named list with intial environ concentrations of all model reactions
#'
#' @param mod model
#' @param initConcentrations Named list with environ initial concentrations
#'
#' @return named list initial environ concentration for the model 
#'
#' @noRd
getconcentrations <- function(mod, initConcentrations){

    concentrations <- -mod$model@lowbnd
    names(concentrations) <- mod$model@react_id
    concentrations[names(initConcentrations)] <- initConcentrations
    # remove environ subtrate that are not present in the model
    concentrations <- concentrations[mod$model@react_id]
    return(concentrations)
}

#' Model substrate uptakes limit for the model
#'
#' Gets a vector with uptake limits for the model reactions
#'
#' @param model model to get uptake limits
#' @param excReact Model exchange reactions
#' @param concentrations Environment concentrations
#' @param biomass Environment biomass
#' @param timeStep Simulation time step
#'
#' @return Named vector with model reactions uptake limits
#'
#' @noRd
getUptakeBound <- function(model, excReact, concentrations, biomass, timeStep){

    # get union names from excReact and concentrations
    namesreacts <- unique(unlist(c(model@react_id, names(concentrations))))

    # Create vars modelub and concent with all the union reactions
    modelub <- numeric(length(namesreacts))
    names(modelub) <- namesreacts
    concent <- modelub

    # init vars
    if (length(modelub[model@react_id]) != length(model@lowbnd)) {
        logWarn('Length differs!!', 'length(modelub[model@react_id])', length(modelub[model@react_id])
               ,'length(model@lowbnd)', length(model@lowbnd))
    }
    modelub[model@react_id] <- -model@lowbnd
    concent[names(concentrations)] <- concentrations

    uptakeBound <- concent/(biomass*timeStep);
    uptakeBound <- ifelse(uptakeBound > modelub, modelub, uptakeBound)

    uptakeBound <- ifelse(abs(uptakeBound) < sybil::SYBIL_SETTINGS("TOLERANCE"),0,uptakeBound)

    return(uptakeBound[model@react_id])    

}

#' optimize Linear Programming problem
#'
#' optimize LP program to get max value of objective function
#'
#' @param mod model to optimize
#' @param excReactIdx Exchange reactions index to update uptake bounds
#' @param uptakeBound Uptake values for exchange reactions
#' @param method Method to use for optimization
#' @param algorithm Algorithm for optimization
#'
#' @return List object with the optimization results
#'
#' @noRd
optimizeMADFBA <- function(mod, excReactIdx, uptakeBound, method = 'FBA', algorithm = 'fba') {
    optsol.OK <- solverStatOK()
    mod$model@lowbnd[mod$excReactIdx]  = -uptakeBound[mod$excReactIdx]
    if (method == "directFBA") {
        sol = sybil::optimizeProb(mod$model, algorithm = algorithm, retOptSol=FALSE);
        ## check solution status (did we find a solution?)
        if ( (sol$ok == 0) && (sol$stat == optsol.OK) ){   ## checkSolStat
            # there is a solution, so check
            lfba.sol <- sybil::optimizeProb(mod$lpmod
                                            ,react=excReactIdx
                                            ,lb=mod$model@lowbnd[excReactIdx]
                                            ,ub=mod$model@uppbnd[excReactIdx]
            )
            
            if ((sol$obj - lfba.sol$obj) > sybil::SYBIL_SETTINGS("TOLERANCE")) cat("obj diff ERROR: ", sol$obj, lfba.sol$obj,  '\n') #stepNo, '\n')
            #flx <- (sol$fluxes - lfba.sol$fluxes)
            #if (length(flx[flx > sybil::SYBIL_SETTINGS("TOLERANCE")])>0) cat("fluxes diff ERROR: ", stepNo, flx[flx > sybil::SYBIL_SETTINGS("TOLERANCE")], '\n')
        }
    }
    else if (method == "FBA") {
        sol = sybil::optimizeProb(mod$lpmod
                                ,react=excReactIdx
                                ,lb=mod$model@lowbnd[excReactIdx]
                                ,ub=mod$model@uppbnd[excReactIdx]
                                )
    }         
    else if ((method == "MTF") || (method == "directMTF") || (method == "lpMTF")) {
        mtf.sol <- NULL
        # use MTF instead of FBA
        # choose one of these three
        if (method == "lpMTF") {
            # this is the two-step lpmod-based approach:
            # first solve FBA to get a value for the objective function(s)
            #sol = sybil::optimizeProb(lpmod);
            sol <- sybil::optimizeProb(mod$lpmod, react=excReactIdx
                                    ,lb=mod$model@lowbnd[excReactIdx]
                                    ,ub=mod$model@uppbnd[excReactIdx]
                                    )
            # then solve MTF with the computed objective
            #if ( length(sybil::checkSolStat(sol$stat, mod$lpmod@problem@solver))==0 ){
            if ( (sol$ok == 0) && (sol$stat == optsol.OK) )   ## checkSolStat
                mtf.sol <- sybil::optimizeProb(mod$model, algorithm="mtf", mtfobj=sol$obj)
            # else{
            #     cat(warn, 'lpMTF: No feasible solution - nutrients exhausted in step ', stepNo,  'for model', m, '->', mod$model@mod_name, '\n');
            #     #break; multimodel break by mods_stats
            # }
            
        } else if (method == "directMTF") {
            # first solve FBA to get a value for the objective function(s)
            sol <- sybil::optimizeProb(mod$model, algorithm="fba", retOptSol=FALSE)
            # then solve MTF with the computed objective
            if ( (sol$ok == 0) && (sol$stat == optsol.OK) )   ## checkSolStat
                mtf.sol <- sybil::optimizeProb(mod$model, algorithm="mtf", mtfobj=sol$obj)
            # else {
            #     cat('\n')
            #     cat(warn, 'directMTF: No feasible solution - nutrients exhausted in step ', stepNo,  'for model', m, '->', mod$model@mod_name, '\n');
            #     #break; multimodel break by mods_stats
            # }
        } else { # if (method == "MTF")
            # this one does an FBA internally and so does not need
            # a previous FBA calculation
            mtf.sol <- optimizeProb(mod$model, algorithm="mtf") 
            # should give a message saying FBA is computed
            # but it doesn't, although it seems to work 
            # correctly nevertheless, I need to trace the 
            # internal logic of Sybil.
        }
        # NOTE that we could also call MTF it with retOptSol=FALSE
        # as we did in FBA. This is only a matter of taste.

        if (!is.null(mtf.sol) && (mtf.sol@lp_ok == 0) && (mtf.sol@lp_stat == optsol.OK) ){   ## checkSolStat
            mtf.ok <- exit_code(checkOptSol(mtf.sol, onlywarn=FALSE))
            mtf.stat <- status_code(checkOptSol(mtf.sol))
            mtf.obj <- mod_obj(mtf.sol)
            mtf.fluxdist <- getFluxDist(mtf.sol)
        } else         #init mtf.sol & mtf.stat 
            mtf.ok <- mtf.stat <- mtf.obj <- mtf.fluxdist <- -1
        
        # convert to list
        sol <- list(
                    ok     = mtf.ok,
                    stat   = mtf.stat,
                    obj    = mtf.obj,
                    fluxes = unname(mtf.fluxdist)
                    #, exchRxnFluxes = mtf.flux.exchRxn
                    )
    } 
        
    return(sol)

}



#' Get Fluxes Variability analisys in every simulation step
#'
#' EXPERIMENTAL
#' Gets, in every step, the flux variability analisys of all reactions, saving
#' the results in files in the working directory
#'
#' @param mod model to get FVA
#' @param stepNo step simulation number
#'
#' @noRd
computeFVA <- function (mod, stepNo) {
    # Do FVA and save the output
    ranges <- fluxVar(mod$model, percentage=80, verboseMode=MADFBA_SETTINGS("VERBOSE"))
    print("\n\nBeing playful, aren't you? ;-)\n\n")

    tabFVAvalues <- paste(mod$model@mod_name, 'ADFBA-FVA-', stepNo, '.tab', sep='')
    cat(paste('\nSaving optimized output from FVA to "', tabFVAvalues,'"\n', sep=""))
    cat(paste('\nTotal length: ', length((ranges)), 
        '\nNo. of reactions: ', length(ranges)/2, '\n\n'))

    nreact=length(ranges)/2
    for (i in 1:nreact) {
    if (MADFBA_SETTINGS("VERBOSE") > 6) {
        print(paste('FVA,', stepNo, ',', i, ',', mod$model@react_id[i], ',', 
            lp_obj(ranges)[i], ',', lp_obj(ranges)[i+nreact]))
    }
    out <- paste(i, '       ', mod$model@react_id[i], '        ', 
            lp_obj(ranges)[i], '        ', lp_obj(ranges)[i+nreact])
    cat(out, file=tabFVAvalues, sep="\n", append=TRUE)
    }
    print(paste('--- END FVA', stepNo, '---\n'))
}

#' Apply death rate to the model
#'
#' Apply exponential death rate to the model when there is no growth posibility
#'
#' @param mod model to apply death rate
#' @param timeStep Step simulation time
#' @param fld all fluxes?
#' @param contype concentration type
#'
#' @return update model with death rate applied
#' 
#' @noRd
applyModDeath <- function (mod, timeStep, fld, contype) {
    death <- mod$deathrate * mod$biomass * timeStep
    mod$biomass <- mod$biomass - death
    mod$biomassVec <- c(mod$biomassVec,mod$biomass);
    # no growth -> no fluxes
    if (fld)
        mod$all_fluxes = cbind(mod$all_fluxes, numeric(mod$model@react_num))

    # Repeat last concentration
    mod$concentrationMatrix <- cbind(mod$concentrationMatrix, mod$concentrations[mod$excReact@react_pos])
    # There is no solution to LP problem therefore there is no fluxes, no uptakes 
    # in this step, add 0
    if (contype == 'totaluptake' & length(dim(mod$exUptakeMatrix)) > 1)
        mod$exUptakeMatrix <- cbind(mod$exUptakeMatrix, mod$exUptakeMatrix[, ncol(mod$exUptakeMatrix)])
    else
        mod$exUptakeMatrix <- cbind(mod$exUptakeMatrix, numeric(length(mod$excReact@react_num)))
        
    return(mod)
}


# from adynFBA.R
# Plot active exchange reactions normalized to 1 (this is useful to
# better visualize relationships between metabolite consumption rate
# changes)
#' Plot active exchange reactions
#'
#' Plot active exchange reactions normalizing values if selected
#'
#' @param dfs optosol_dynamicFBA object from MADFBA simulation
#' @param normalize option to normalize data
#'
#' @export

plotExchMets <- function(dfs, normalize=FALSE, main = 'Concentrations') {
    times <- dfs@timeVec
    concs <- as.matrix(dfs@concentrationMatrix)
    # if we want to plot only non-zero ExcRxns
    concs <- concs[apply(concs[, -1], 1, function(x) !all(abs(x)<SYBIL_SETTINGS('TOLERANCE'))),]

    # normalize the row values
    #concs <- t(apply(concs, 1, function(x)(x-min(x))/(max(x)-min(x))))
    if (normalize) {
        concs <- t(apply(concs, 1, function(x)(x)/(max(x))))
        main <- paste('Normalized', main)
    }


    # choose color palette
    nrxns <- dim(concs)[1]
    cl <- heat.colors(nrxns)
    cl <- terrain.colors(nrxns)
    cl <- topo.colors(nrxns)
    cl <- cm.colors(nrxns)
    cl <- rainbow(nrxns)
    col <- 0    

    # find plot limits
    xmin <- min(sapply(times, function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
    xmax <- max(sapply(times, function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
    ymin <- min(sapply(concs, function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
    ymax <- max(sapply(concs, function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
    # create inital empty plot
    plot(NULL, NULL, xlim=c(xmin, xmax), ylim=c(ymin,ymax), type="l", 
         main=main, ylab = "mmol",xlab='Time(hrs)')
    # plot the metabolites
    for (i in row.names(concs)) {
        col <- col + 1
        lines(concs[i, ], col=cl[col], type='l', lty=col)
        #points(concs[i, ], col=cl[col], pch=col, cex=0.5)
    }

    legend('topleft', legend=rownames(concs), col=cl, lwd=1, cex=0.5)
}

