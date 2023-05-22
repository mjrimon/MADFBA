#' @param dynamicConstraints    Values of dynamically changing reaction rates
#'                      Default: \code{NULL}
#
#' @param nutrientChanges       A facility to simulate changing concentrations
#'                      of nutrients: it may be either an array of delta values
#'                      to add at each time point, or a function that will
#'                      take the current concentrations and return the new ones (the old concentrations
#'                      will be substituted by the new ones, so ensure you return all).
#'                      Default: \code{NULL}
#' @noRd

# Check variables needed inside function
doNutrientChanges <- function (nutrientChanges, mods, concentrations, timeStep, stepNo, ...){
    if (is.function(nutrientChanges)) {
        # the change may be a function of current concentration or
        # of time
        concentrations <- nutrientChanges(concentrations, timeStep, stepNo)
    } else {
        ncNames <- names(nutrientChanges)
        # Check file format, v2 (new format) has fixed columns
        if (length(ncNames)>=3 & any(grepl("Time", ncNames, ignore.case = TRUE))
                               & any(grepl("Substrate", ncNames, ignore.case = TRUE))
                               & any(grepl("value", ncNames, ignore.case = TRUE))
           )
            concentrations <- doNutrientChanges_v2file(nutrientChanges, mods, concentrations, timeStep, stepNo, ...)
        else
            concentrations <- doNutrientChanges_v1file(nutrientChanges, mods, concentrations, timeStep, stepNo, ...)
    }
    return(concentrations)
}



doDynamicConstraints <- function (dynamicConstraints, model, concentrations, sol, timeStep, stepNo, ...){
    # First, modify the model to reflect time-specific constraints
    if (is.function(dynamicConstraints)) {
        # This should allow modifying the model depending on current
        # conditions, e.g. if we detect that time < 10h induce a lag
        # phase of growth, if we detect that some metabolite is low,
        # repress coupled uptake reactions or activate alternative
        # uptake mechanisms...
        # We also pass the last fluxes, so that inside the function
        # we may also detect (if desired) activation or repression of
        # a key reaction/route and trigger the corresponding changes.
        model <- dynamicConstraints(model, concentrations, sol$fluxes, timeStep, stepNo)
    } else {
        dcNames <- names(dynamicConstraints)
        # Check file format, v2 (new format) has fixed columns
        if (length(dcNames)>=4  
                               & any(grepl("Time", dcNames, ignore.case = TRUE))
                               & any(grepl("react_id", dcNames, ignore.case = TRUE))
                               & any(grepl("bound", dcNames, ignore.case = TRUE))
                               & any(grepl("value", dcNames, ignore.case = TRUE))
           )
            model <- doDynamicConstraints_v2file(dynamicConstraints, model, concentrations, sol, timeStep, stepNo, ...)
        else
            model <- doDynamicConstraints_v1file(dynamicConstraints, model, concentrations, sol, timeStep, stepNo, ...)
    } # if (is.function(dynamicConstraints))
    return(model)      
}


doNutrientChanges_v1file <- function (nutrientChanges, mods, concentrations, timeStep, stepNo, ...){
    # nutrientChanges MUST be a table of concentration deltas at
    # the specified time
    #  Time column must contains 'time' string in its name (no case sensitive)
    ncTimeInd <- grepl("Time", names(nutrientChanges), ignore.case = TRUE)
    ncTimeIdx <- which(ncTimeInd)
    if (length(ncTimeIdx) == 0) return(concentrations)  # NO time column, can't continue

    tn <- round(timeStep*stepNo, 3) # timeStep*stepNo  strange results for stepNo 239 & 240 with timeStep = 0.1
                                    # not same precision in t>=tn than t < (tn + timeStep) :(
                                    # miliHours

    ## Get rows for this time step. Time must be between tn and tn+1 (tn + timeStep, next time step)
    t <- nutrientChanges[,ncTimeIdx]            # Get all times with nutrient changes 
    rowncInd <- t>=tn & t < (tn + timeStep)     # Rows with nutrient changes for this timestep [tn, tn+1)

    # Is there any nutrient change to apply at this step? Sum all posible rows
    if (any(rowncInd)) {
        nc <- nutrientChanges[rowncInd, -ncTimeIdx]
        nc <- as.matrix(nc)
        delta <- colSums(nc) 
        if (any(delta != 0)) {
            concentrations[names(delta)] <- concentrations[names(delta)] + delta  # add/sustract concentration
            
            logPrefix <- deparse(sys.call(0)[[1]])   # sys.call(-1) for caller function name
            logInfo2('Applied Nutrient Changes, model ->', model@mod_id, 'Time step ->', tn, '/', tn+timeStep)
            logInfo(logPrefix, delta)
        }
    }
    return(concentrations)
}

doNutrientChanges_v2file <- function (nutrientChanges, mods, concentrations, timeStep, stepNo, ...){
    # nutrientChanges MUST be a table of concentration deltas at
    # the specified time
    # XXX mj XXX 
    # WE NEED TO KNOW PREVIOUS DATA FORMAT. if we want to keep compatibility
    # assume 3 columns table -> Time, Substrate, value
    # get time of current time step tn
    tn <- round(timeStep*stepNo, 3) # timeStep*stepNo  strange results for stepNo 239 & 240 with timeStep = 0.1
                                    # not same precision in t>=tn than t < (tn + timeStep) :(
                                    # miliHours

    ## Get rows for this time step. Time must be between tn and tn+1 (tn + timeStep, next time step)
    # Find Time column, case insentive, can specify units valid names eg. 'TIME', 'Time (h)' ...
    ncTimeInd <- grepl("Time", names(nutrientChanges), ignore.case = TRUE)
    ncTimeIdx <- which(ncTimeInd)
    if (length(ncTimeIdx) == 0) return(concentrations)  # NO time column, can't continue

    t <- nutrientChanges[,ncTimeIdx]            # Get all times with nutrient changes 
    rowncInd <- t>=tn & t < (tn + timeStep)     # Rows with nutrient changes for this timestep [tn, tn+1)

    # Is there any nutrient change to apply at this step?
    if (any(rowncInd)) {
        # Find value column, case insentive, can specify units valid names eg. 'VALUE', 'Value (mmol/l)' ...
        ncValueInd <- grepl("value", names(nutrientChanges), ignore.case = TRUE)
        ncValueIdx <- which(ncValueInd)

        # Find value column, case insentive, can specify units valid names eg. 'SUBSTRATE', 'substrate' ...
        ncSubstrateInd <- grepl("substrate", names(nutrientChanges), ignore.case = TRUE)
        ncSubstrateIdx <- which(ncSubstrateInd)

        delta <- nutrientChanges[rowncInd, ncValueIdx]             # o ncn[, 3]
        names(delta) <- nutrientChanges[rowncInd, ncSubstrateIdx]  # o ncn[, 2]
        concentrations[names(delta)] <- concentrations[names(delta)] + delta  # add/sustract concentration
        
        logInfo2('Applied Nutrient Changes. Time step ->', tn, '/', tn+timeStep)
        logInfo('', delta)
    }
    return(concentrations)
}


doDynamicConstraints_v1file <- function (dynamicConstraints, model, concentrations, sol, timeStep, stepNo, ...){
    old.model <- model
    ## Coerce exchange rates for this step to the values provided
    #
    #  Time column must contains 'time' string in its name (no case sensitive)
    dcNames <- names(dynamicConstraints)
    dcTimeInd <- grepl("Time", dcNames, ignore.case = TRUE)
    dcTimeIdx <- which(dcTimeInd)
    if (length(dcTimeIdx) == 0) return(model)  # NO time column, can't continue

    dcTimeIdx <- dcTimeIdx[1]   # In case more than 1 column choose the first one -> Reaction with string 'time' in name?
                                # Info needed to user?

    ## Get dynamicConstraints rows for this timeStep

    ## Get rows for this time step. Time must be between tn and tn+1 (tn + timeStep, next time step)
    tn <- round(timeStep*stepNo, 3)          # max precision miliHour, see NutrientChanges
    t <- dynamicConstraints[,dcTimeIdx]      # Get all times with nutrient changes 
    rowdcInd <- t>=tn & t < (tn + timeStep)  # Rows for this timestep [tn, tn+1)

    # Apply only last row ot this time step. It's possible, for testing, to specify a large timeStep and
    # find several rows, because of tab file definition the last row contains all reacts values, no need 
    # to apply all rows.
    if (any(rowdcInd)) {
        rowdcIdx <- max(which(rowdcInd))         # Apply only last dc row of this time step.
        # Remove time column 
        dcn <- dynamicConstraints[rowdcInd, -dcTimeIdx]

        # get reacts names
        reacts <- sub("(\\[|\\.)(low|upp)(\\]|\\.)", "", names(dcn))

        # Get upp col bounds
        dcubInd <- grepl("(\\[|\\.)upp(\\]|\\.)", names(dcn), ignore.case = TRUE)
        # Get low col bounds
        dclbInd <- grepl("(\\[|\\.)low(\\]|\\.)", names(dcn), ignore.case = TRUE)

        # add low&upp col bounds
        lbub <- !(dcubInd | dclbInd)
        dcubInd <- dcubInd | lbub
        dclbInd <- dclbInd | lbub

        # update dcn to list with react names
        dcn <- unlist(dcn)
        names(dcn) <- reacts

        # get lowbnd dynamic constrains and remove reactions not in model
        dclb <- dcn[dclbInd]
        dclb <- dclb[names(dclb) %in% model@react_id]

        # get uppbnd dynamic constrains and remove reactions not in model
        dcub <- dcn[dcubInd]
        dcub <- dcub[names(dcub) %in% model@react_id]

        # Apply constraints to model
        model@lowbnd[match(names(dclb), model@react_id)] <- dclb
        model@uppbnd[match(names(dcub), model@react_id)] <- dcub


        logPrefix <- deparse(sys.call(0)[[1]])   # sys.call(-1) for caller function name
        logInfo2('Applied Dynamic Constrains, model ->', model@mod_id, 'Time step ->', tn, '/', tn+timeStep)
        logInfo(logPrefix, 'lowbnd ->', names(dclb), 'old:', old.model@lowbnd[match(names(dclb), model@react_id)], 'new:', dclb)
        logInfo(logPrefix, 'uppbnd ->', names(dcub), 'old:', old.model@uppbnd[match(names(dcub), model@react_id)], 'new:', dcub)
    }
    return(model)

}


doDynamicConstraints_v2file <- function (dynamicConstraints, model, concentrations, sol, timeStep, stepNo, ...){
    old.model <- model
    ## Coerce exchange rates for this step to the values provided
    #
    # keepeing compatibility one file per model   -> field names:  "Time	react_id	bound	value"
    # option: same constraint file for all models -> field names:  "Time	mod_id	react_id	bound	value"

    ## Get dynamicConstraints rows for this model
    ## Get rows for this time step. Time must be between tn and tn+1 (tn + timeStep, next time step)
    tn <- round(timeStep*stepNo, 3)          # max precision miliHour, see NUtrientChanges

    dcNames <- names(dynamicConstraints)
    dcTimeInd <- grepl("Time", dcNames, ignore.case = TRUE)
    dcTimeIdx <- which(dcTimeInd)
    if (length(dcTimeIdx) == 0) return(model)  # NO time column, can't continue

    t <- dynamicConstraints[,dcTimeIdx]        # Get all times with nutrient changes 
    rowtnInd <- t>=tn & t < (tn + timeStep)    # Rows for this timestep [tn, tn+1)

    # Get rows with dynamicConstrains for this model at this time step

    # option: same constraint file for all models
    if (length(dcNames >= 5) & any(grepl("mod_id", dcNames, ignore.case = TRUE))) {
        dcModInd <- grepl("mod_id", dcNames, ignore.case = TRUE)
#        dcModIdx <- which(dcTimeInd)
        rowmodInd <- dynamicConstraints[, dcModInd] == model@mod_id
        rowdcInd <- rowmodInd & rowtnInd
    }else{
        rowdcInd <- rowtnInd
    }

    if (any(rowdcInd)) {
        # Find value column, case insentive, can specify units valid names eg. 'VALUE', 'Value (mmol/l)' ...
        ncValueInd <- grepl("value", names(dynamicConstraints), ignore.case = TRUE)
        ncValueIdx <- which(ncValueInd)

        dcn <- dynamicConstraints[rowdcInd, ]                      # Get rows with dc for this timestep and model
        dclbInd <- grepl("l", dcn[, 'bound'], ignore.case = TRUE)  # Get rows with lowbnd dynamic constraints
        dcubInd <- grepl("u", dcn[, 'bound'], ignore.case = TRUE)  # Get rows with uppbnd dynamic constraints

        # get lowbnd dynamic constrains and remove reactions not in model
        dclb <- dcn[dclbInd, ncValueIdx]
        names(dclb) <- dcn[dclbInd, 'react_id']
        dclb <- dclb[names(dclb) %in% model@react_id]

        # get uppbnd dynamic constrains and remove reactions not in model
        dcub <- dcn[dcubInd, ncValueIdx]
        names(dcub) <- dcn[dcubInd, 'react_id']
        dcub <- dcub[names(dcub) %in% model@react_id]


        model@lowbnd[model@react_id %in% names(dclb)] <- dclb
        model@uppbnd[model@react_id %in% names(dcub)] <- dcub


        logInfo2('Applied Dynamic Constrains, model ->', model@mod_id, 'Time step ->', tn, '/', tn+timeStep)
        logInfo('lowbnd ->', names(dclb), dclb)
        logInfo('uppbnd ->', names(dcub), dcub)
    }
    return(model)
}
