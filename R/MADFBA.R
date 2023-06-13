
#' Simulation Multimodel Adaptive Dynamic FBA 
#'
#' Runs a simulation of the evolution of the test 
#' with one or several models according to the specified parameters
#'
#' @param models A model or a list of models of type sybil::modelorg 
#' @param substrateRxns Name of the substrates of the reactions available in the environment
#' @param initConcentrations Initial concentration of environment substrates
#' @param initBiomass An integer or a vector with the initial biomass of the models 
#' @param timeStep Simulation time step (hours)
#' @param nSteps Maximum number of steps the simulation will take (may end earlier)
#' @param exclUptakeRxns Consumption reactions that are excluded from the simulation. 
#'                       If not specified, 'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)' are excluded by default.
#' @param retOptSol Sets the type of the result: TRUE an optsol_dynamicFBA object; FALSE a list
#' @param fld Indicates whether all flows of all simulation steps should be included in the result.
#' @param verboseMode Level of information presented
#' @param biomassRxn Name of the model growth reaction of the biomass, if not specified it is obtained by default.
#' @param dynamicConstraints Bound reactions limits for each model or callback function.
#' @param nutrientChanges Environment changes, table or callback function
#' @param method Method to use in the solver to optimize the problem
#' @param logModel Model for which to display the log during the simulation, default 0, the sum of all
#' @param contype Type of consumed/excreted concentration data per model to be saved at each step:
#'                    'legacy' apparent substrate concentration, initialized at initial substrate concentrations. 
#'                    'uptakerate' amount of substrate consumed/excreted at that instant.
#'                    'totaluptake' total amount of substrate consumed/excreted up to that instant.  
#' @param deathrate  models deathrate
#'                   If different from 0, default value, apply as death rate of the organism when running out of nutrients.
#' @param ... 
#'
#' @return An optsol_dynamicFBA object or a list with the result of a model 
#'            or a list of optsol_dynamicFBA objects or lists if the simulation is a multi-model simulation.
#' @export
#'
#' @examples
#'
#' # Load models data
#' data(Ec_core)
#' data(Slividans)
#' 
#' # set the models uptake bounds
#' Ec_core <- sybil::changeBounds(Ec_core, react = "EX_glc(e)", lb = -12)
#' Ec_core <- sybil::changeBounds(Ec_core, react = "EX_o2(e)",  lb = -10)
#' Ec_core <- sybil::changeBounds(Ec_core, react = "EX_ac(e)",  lb = -10)
#' 
#' Slividans <- sybil::changeBounds(Slividans, react = "EX_mnl(e)", lb = -7)
#' Slividans <- sybil::changeBounds(Slividans, react = "EX_glc(e)", lb = -5)
#' 
#' # define medium, substrate names and initial concentrations
#' init.source <- c("EX_ac(e)","EX_o2(e)","EX_glc(e)", 'EX_mnl(e)')
#' init.conc   <- c(10,1,28,4)

#' # set models and initial biomass for simulation
#' init.models <- c(Ec_core,Slividans)
#' init.bmass <- c(.01, .0075)
#' Ec_df_min  <- MADFBA(init.models,  #exclUptakeRxns = c(),  # Si no se especifican exclusiones por defecto exluye, al menos o2
#'                               substrateRxns      = init.source,
#'                               initConcentrations = init.conc,
#'                               initBiomass        = init.bmass,
#'                               timeStep=.1,nSteps=200,verbose=3,
#'                               method = "FBA",                   # "directMTF" #"FBA" #"MTF"  #"directFBA"
#'                               fld=TRUE,                         # All fluxes para poder representarlos luego si procede
#'                               retOptSol = TRUE
#'                               )


MADFBA <- function (models,                # 1 or list with n models.
                substrateRxns,             # substrates that exist in the environment, s substrates. Available for all models
                initConcentrations,        # Initial concentrations of the surrounding substrates. If a single value is specified for s substrates, they are all initialized with the same value.
                initBiomass,               # 1 or list with n biomass of each model. If there are n models and 1 biomass, the same biomass is applied to all the models
                timeStep,
                nSteps,
                exclUptakeRxns = c('EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'),  
                retOptSol = TRUE,
                fld = FALSE,
                verboseMode = 2, 
                biomassRxn = "",                # We will always need the reaction that defines the biomass, even if the objective is different, the biomass is necessary to calculate the growth.
                dynamicConstraints = NULL,      # Bound reactions limits for each model or callback function.
                nutrientChanges = NULL,         # Environment changes, table or callback function
                method = "FBA",                 
                logModel = 0,                   # model to display the log of during simulation, default 0, summary of all the models
                contype = 'legacy',             # type of consumed/excreted concentration data per model to be saved at each step:
                                                # 'legacy' apparent substrate concentration, initialized at initial substrate concentrations. 
                                                # 'uptakerate' amount of substrate consumed/excreted at that instant.
                                                # 'totaluptake' total amount of substrate consumed/excreted up to that instant.  
                deathrate = 0,
                ...){                    

    # aids for debugging/verbose messages
    funcName <- as.character(match.call()[[1]]) # get function name
    # or
    funcName <- deparse(sys.call(0)[[1]]) # get function name
    logPrefix <- paste('# mj //', funcName, '=> ')
    info <- paste(logPrefix, 'INFO: ')
    warn <- paste(logPrefix, 'WARNING: ')
    err  <- paste(logPrefix, 'ERROR: ')

    
    #optsol.OK <- 5 # only for glpk solver!!!! , use solverStatOK(solver)
    # Ok status code of linear optimization depending on solver
    optsol.OK <- solverStatOK()

    # use sybil logclass, log comments to screen and/or file
    # start logging
    # logObj <- sybil::sybilLog(filename = "",
    #                           loglevel = -1,
    #                           verblevel = verboseMode)
    MADFBA_SETTINGS("VERBOSE", verboseMode)

    ##--------------------------------------------------------------------------##
    # check prerequisites 
    #
    #  Environment, Sustrates (medium), method

    ##--------------------------------------------------------------------------##
    # initConcentrations must be same length that substrateRxns. If its lentgh is = 1
    # means all substrateRxns have same init concentration
    # 'border case', only 1 substrate,
    if (length(initConcentrations) == 1) {
        initConcentrations <- rep(initConcentrations[1], length(substrateRxns))
    }
    
    stopifnot(length(initConcentrations) == length(substrateRxns))

    # Name init concentrations values for use with model exchange reactions
    # substrateRxns and initConcentrations are the environment for all models
    # this handles unordered substrateRxns definition. 
    # Use > concentrations[names(initConcentrations)] <- initConcentrations

    names(initConcentrations) <- substrateRxns

    stopifnot(method %in% c('FBA', 'directFBA', 'MTF', 'lpMTF', 'directMTF'))
    
    stopifnot(contype %in% c('legacy', 'uptakerate', 'totaluptake'))


    # # Models

    mods <- makeOrganisms(models, initBiomass, biomassRxn, exclUptakeRxns, dynamicConstraints, deathrate)

    mod <- mods[[1]]

    # show log info
    if (logModel > length(mods)) 
        logModel <- length(mods)

    if (logModel==0) 
        logInfo2('Showing summary info table of all models')
    else 
        logInfo2('Showing info table of model', mods[[logModel]]$model@mod_id, '\n')

    logSubstrateFH(mod$model, mod$biomass, initConcentrations)

    # We prepare the 'possible medium'. It has to contain all the exchange reactions of all the models.
    # models. All secretions of metabolites are made to the medium and all consumptions are # obtained from the medium, which is common for all.
    # obtained from the medium, which is common for all. Possibly a function to obtain the medium
    
    mednames = list()
    mednames <- lapply(mods, function(mod) return(c(mednames, mod$excReact@react_id)))
    mednames <- unique(unlist(mednames))
    medium <- numeric(length(mednames))
    names(medium) <- mednames
    medium[names(initConcentrations)] <- initConcentrations

    mediumMatrix <- medium

    mods <- addOrgsToMedium(mods, initConcentrations, contype)


    mediumStepFluxes0 <- numeric(length(initConcentrations))
    names(mediumStepFluxes0) <- substrateRxns

    # add biomass to medium. commodity vars, we have all data like the sums of the models.
    # look for better solution
    # 
    mediumBiomass <- sum(sapply(mods, function(model) model$biomass))
    mediumBiomassVec <- mediumBiomass
    
    # In principle we do not store the flow matrix of the medium, as we can always access the sum.
    # Prepare a function to do the sum of matrices with different number of rows
    # (reactions of each model) and different number of columns (Steps where the 
    # model has been growing)
    
    timeVec <- 0

    for (stepNo in 1:nSteps){
        logInfo('\n','',type='')
        logInfo(logPrefix, 'Step', stepNo)
        # The changes of nutrients in the substrate are at the level of 'step', they belong to the
        # environment, and affect all models
        if ( ! missing(nutrientChanges) && ! is.null(nutrientChanges) )  {
                medium <- doNutrientChanges(nutrientChanges, mods, medium, timeStep, stepNo)
            } # if there are nutrientChanges
        
        mediumStepBiomass <- 0
        mediumStepFluxes <- mediumStepFluxes0
        for (m in 1: length(mods)) {
            mod <- mods[[m]]

            ##-----------------------------------------------------------------------------##
            # Do the simulation
            if ((! missing(dynamicConstraints)) && (! is.null(mod$dynamicConstraints))) {
                dc_model <- doDynamicConstraints(mod$dynamicConstraints, mod$model, medium, sol, timeStep, stepNo)
                # Have there been changes in the model? 
                # If the changes are in the objective function, recalculate the problem.
                # if so recalculate problem with new data
                if (!isEqualModel(mod$model, dc_model)) {
                    if (!isEqualModelReactBounds(mod$model, dc_model, mod$biomassRxn)) 
                        mod$lpmod <- sybil::sysBiolAlg(dc_model, algorithm = "fba")
                    mod$model <- dc_model
                }
            } # if there are dynamicConstraints

            uptakeBound <- getUptakeBound(mod$model, mod$excReact, medium, mod$biomass, timeStep)

            #mod$lpmod <- sybil::sysBiolAlg(mod$model, algorithm = "fba")
            sol <- optimizeMADFBA(mod, mod$excReactIdx, uptakeBound, method)

            ### j ###
            # UNDOCUMENTED FEATURE
            #       WORKS IN PROGRESS!!!    WORKS IN PROGRESS!!!
            #       Compute FVA at each time point
            #       THIS WILL SLOW DOWN CALCULATIONS HORRIBLY
            #       BUT WILL ALLOW US TO COLLECT ADDITIONAL INFORMATION FOR
            #       OUR STATISTICAL ANALYSES
            #       WE STILL NEED TO THINK OUT THE BEST WAY TO PRODUCE THE
            #       OUTPUT SO THAT IT IS AMENABLE FOR STATISTICAL ANALYSIS
            #       WITHOUT BREAKING COMPATIBILITY WITH sybil::dynamicFBA
            #
            #       Probably we should store them in a data frame with two
            #       matrices (low/high) and then save the data frame at the
            #       end of the calculation. Choice of time as rows or columns
            #       will depend on the statistics we finally use, but this
            #       still is a work in progress. Or maybe we should use only
            #       two long columns and add the time step as a factor, or
            #       use a list...
            #
            #       For now, we'll do with creating a series of output files.
            #
            if (verboseMode > 5) computeFVA(mod, stepNo)

            # add sol$ to all_stat vector
            mod$all_stat = c(mod$all_stat,sol$stat)
            # save mod values in case sol$stat != optsol.OK 
            #mods[[m]] <- mod

            # check stat. If this model has exhausted metabolites, go to the next model
            if ( (sol$ok != 0) || (sol$stat != optsol.OK)) {
                # Add last calculated biomass to accumulate for the medium
                # If we add death rate we will compute here 
                # XXX mj XXX Test death rate hard coded
                if (mod$deathrate > 0) {
                    mod <- applyModDeath(mod, timeStep, fld, contype)
                }
                # We add the remaining biomass to the biomass of the medium
                mediumStepBiomass <- mediumStepBiomass + mod$biomass
                # We can inform lack of nutrients the first time it occurs or at all steps
                # for all steps comment next line
                if (length(mod$all_stat[mod$all_stat!= optsol.OK]) == 1)
                    cat(warn, method, ': No feasible solution - nutrients exhausted in step ', stepNo,  'for model', m, '->', mod$model@mod_name, '\n');

                # show log model step if selected
                if (m==logModel) logConcentFluxes(mod$biomass, mod$concentrations[mod$substrateRxnsIdx], sol$fluxes[mod$substrateRxnsIdx], stepNo)        

                mods[[m]] <- mod
                next;
            }
            
            # mu is the growth rate, the flow rate of the biomass reaction or growth. 
            # Regardless of what the objective function is for calculating the outcome of the problem, 
            # we need to know how much the biomass grows so that we can calculate 
            # the metabolite that has been consumed or excreted,
            # with that data, calculate the metabolite consumed or excreted.
            # varma-palsson 1994 equations (6) and (7)
            mu <- sol$fluxes[mod$biomassIdx]
            logInfo(logPrefix, 'Model:', mod$model@mod_id, 'Objective at step', stepNo, '(', stepNo * timeStep, 'h ):', mu)

            # get uptake fluxes after FBA/MTF
            uptakeFlux = sol$fluxes[mod$excReactIdx];

            medium[mod$excReactNames]= medium[mod$excReactNames] - uptakeFlux/mu*mod$biomass*(1-exp(mu*timeStep));
            # No substrate can be negative 
            #medium[medium < sybil::SYBIL_SETTINGS("TOLERANCE")] <- 0
            medium[medium < 0] <- 0

            # keep track of model concentrations changes
            mod$concentrations[mod$excReactIdx]= mod$concentrations[mod$excReactIdx] - uptakeFlux/mu*mod$biomass*(1-exp(mu*timeStep));

            # save concentrations
            mod$concentrationMatrix <- cbind(mod$concentrationMatrix,mod$concentrations[mod$excReact@react_pos]);

            # kepp track of uptakes and secretions, no exluded reactions, in order to know which reacts
            # and quatity are needed to achive this step
            # These calculations can be done at the end of the simulation, with the result object, if we store
            # all fluxes at all steps. In this case this could be redundant. In this way we can give all
            # the necessary information about the simulation in one place. No further computations needed.
            #mod$exUptakeMatrix <- cbind(mod$exUptakeMatrix, - sol$fluxes[mod$excReact@react_pos]/mu*mod$biomass*(1-exp(mu*timeStep)))
            if ( contype == 'totaluptake' & length(dim(mod$exUptakeMatrix)) > 1)
                mod$exUptakeMatrix <- cbind(mod$exUptakeMatrix, mod$exUptakeMatrix[, ncol(mod$exUptakeMatrix)] - sol$fluxes[mod$excReact@react_pos]/mu*mod$biomass*(1-exp(mu*timeStep)))
            else
                mod$exUptakeMatrix <- cbind(mod$exUptakeMatrix, - sol$fluxes[mod$excReact@react_pos]/mu*mod$biomass*(1-exp(mu*timeStep)))

            # add the fluxes to the table of flux changes over time
            if(fld){
                if (stepNo == 1) {
                    mod$all_fluxes = sol$fluxes;
                }else{
                    mod$all_fluxes = cbind(mod$all_fluxes,sol$fluxes);
                }
            }

            # compute new biomass
            # once computed new concentrations we can compute new biomass
            # we need this code 'distribution' for use only one biomass variable
            mod$biomass = mod$biomass*exp(mu*timeStep);
            mod$biomassVec = c(mod$biomassVec,mod$biomass);

            logInfo(logPrefix, 'Model:', mod$model@mod_id, 'Biomass at t =', stepNo * timeStep, mod$biomass)

            # Calculate cumulative data from all models in this step
            # Add the biomass to the biomass of the medium
            mediumStepBiomass <- mediumStepBiomass + mod$biomass

            modStepFluxes <- sol$fluxes[mod$substrateRxnsIdx]
            names(modStepFluxes) <- mod$model@react_id[mod$substrateRxnsIdx]

            # As a multi-model, a model may not have all the metabolites (reactions) found in the substrate, so NA will appear.
            # in the substrate, so NA will appear. We omit them.
            modStepFluxes <- na.omit(modStepFluxes)
            
            # If we omit, delete some value, we must take it into account to update only the available values.
            # available values, for this purpose we use names(mStepFluxes)
            mediumStepFluxes[names(modStepFluxes)] <- mediumStepFluxes[names(modStepFluxes)] + modStepFluxes

            # show log model step if selected
            if (m==logModel) logConcentFluxes(mod$biomass, mod$concentrations[mod$substrateRxnsIdx], sol$fluxes[mod$substrateRxnsIdx], stepNo)        

            mods[[m]] <- mod
        }

        # Check the status of all models. If all have finished, you can't continue.
        # The code goes here so as not to add a new element to timeVec that could not be performed.
        # nor execute the possible change of nutrients that is at the beginning of the loop for a second round.
        # For compatibility of operation with the single model version.
        # If we apply the death rate of the organisms, we have to continue 
        # the simulation as long as there is biomass and obtain the death and accumulated data.
        
        # Is a variable needed to set the limit at which to decide that there is no biomass left?
        # for the moment MADFBA_SETTINGS('TOLERANCE')
        modsDeathSim <- lapply(mods, function(mod) return(mod$deathrate > 0 & mod$biomass > MADFBA_SETTINGS('EXTINCTION')))

        mods_stat <- lapply(mods, function(mod) return(tail(mod$all_stat, 1)))

        if (all(mods_stat != optsol.OK) & !any(modsDeathSim == TRUE)) {
            cat('Simulation terminated. No model can continue growing nor dying')
            break;
        }
        
        mediumMatrix <- cbind(mediumMatrix, medium)
        mediumBiomassVec <- c(mediumBiomassVec, mediumStepBiomass)

        if (logModel==0) 
            logConcentFluxes(mediumStepBiomass, medium[substrateRxns], mediumStepFluxes, stepNo)        
        
        # The simulation is not yet finished, we add a new element to the time vector for the next step.      
        timeVec = c(timeVec,stepNo*timeStep);
    }# end loop
    # close progress bar

    ##--------------------------------------------------------------------------##
    # Simulation completed... prepare RETURN output

    # We prepare the result with the data of the whole environment, the data of the evolution 
    # of the substrates and the total biomass growth of the system.
    resultsummary <- getResultSummary(mods, mediumMatrix, mediumBiomassVec, timeVec, retOptSol) 

    # Check if there is only one model, in which case the summary data, the object or the list is returned.
    # to maintain compatibility in the object returned by ADFBA. The summary data
    # of the system is the data of the single model.
    if (length(mods) == 1) return(resultsummary)

    # There is more than one model. We prepare the data for each model and generate 
    # a list with the results of all the models separately plus one more result
    # with the summary of all the models.

    result <- lapply(mods, getResult, timeVec, fld, retOptSol, contype)

    return(append(result, resultsummary))
 
}

