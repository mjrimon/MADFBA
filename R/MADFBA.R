
#' Simulacion Multimodelo Adaptativo Dinámico FBA
#'
#' Ejecuta una simulación de la evolución del ensayo con uno o varios modelos
#' de acuerdo a los parámetros especificados
#'
#' @param models Un modelo o un lista de modelos de tipo sybil::modelorg 
#' @param substrateRxns Nombre de los substratos de las reacciones disponibles en el entorno
#' @param initConcentrations Concentracion inicial de los substratos del entorno
#' @param initBiomass Un entero o un vector con la biomasa inicial de los modelos 
#' @param timeStep Paso de tiempo de simulación (horas)
#' @param nSteps Numero de pasos que tendrá la simulación como máximo (puede terminar antes)
#' @param exclUptakeRxns Reacciones de consumo que se excluye de la simulacion. si no se especifica, por defecto se excluyen 'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'
#' @param retOptSol Establece el tipo del resultado: TRUE  un objeto optsol_dynamicFBA; FALSE una lista
#' @param fld Indica si se deben incluir en el resultado todos los flujos de todos los pasos de la simulación
#' @param verboseMode Nivel de infrmación presentada
#' @param biomassRxn Nombre de la reacción de crecimiento del modelo, de la biomasa, si no se especifica se obtiene por defecto
#' @param dynamicConstraints 
#' @param nutrientChanges 
#' @param method 
#' @param logModel Modelo del que mostrar el log durante la simulación, por defecto el 0, la suma de todos
#' @param contype Tipo de datos de concentraciones consumidas/excretadas por modelo que se guardarán en cada paso:
#'                    'legacy' concentración de substratos aparente, incializado a concentraciones iniciales de substrato. 
#'                    'uptakerate'  cantidad de substrato consumido/excretado en ese instante.
#'                    'totaluptake' cantidad total de substrato consumido/excretado hasta ese instante.  
#' @param deathrate  deathrate de los modelos
#'                   Si es distinto de 0, valor por defecto, aplicar como tasa de muerte del organismo al quedarse sin nutrientes
#' @param ... 
#'
#' @return Un objeto optsol_dynamicFBA o una lista con el resultado de un modelo o una lista de objetos optsol_dynamicFBA o listas si la simulación es multimodelo
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


MADFBA <- function (models,                # 1 o lista con n modelos. Hapy parametros referentes a los modelos que habrá que definir aquí
                substrateRxns,             # substratos que existen en el entorno, s substratos. Disponibles para todos los modelos
                initConcentrations,        # Concentraciones iniciales de los substratos del entorno. Si se especifica un unico valor para s substratos, se inicializan todos con el mismo valor
                initBiomass,               # 1 o lista con n biomasa de cada modelo. Si hay n modelos y 1 biomasa, se aplica la misma biomasa a todos los modelos, avisar
                timeStep,
                nSteps,
                exclUptakeRxns = c('EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'),  # Valor por defecto, no mas missing checks
                retOptSol = TRUE,
                fld = FALSE,
                verboseMode = 2, 
                biomassRxn = "",                # Siempre necesitaremos la reacción que define la biomasa, aunque el objetivo sea otro, la biomasa es necesaria para calcular el crecimiento
                dynamicConstraints = NULL,      # Referidas a los modelos. Funcion o lista para un modelo. Si hay n modelos... ¿Varias funciones? ¿lista de listas? repensar..
                nutrientChanges = NULL,         # Cambio del entorno
                method = "FBA",                 # Si hay n modelos, ¿Todos el mismo método?, ¿cada modelo un metodo?, ¿Tiene sentido biológico?
                logModel = 0,                   # modelo del que mostrar el log durante la simulación, por defecto el 0, resumen de todos
                contype = 'legacy',             # tipo de datos de concentraciones consumidas/excretadas por modelo que se guardarán en cada paso:
                                                #       'legacy' concentración de substratos aparente, incializado a concentraciones iniciales de substrato. 
                                                #       'uptakerate'  cantidad de substrato consumido/excretado en ese instante.
                                                #       'totaluptake' cantidad total de substrato consumido/excretado hasta ese instante.  
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
    
    # Podrían ser parámetros de configuración...
    stopifnot(method %in% c('FBA', 'directFBA', 'MTF', 'lpMTF', 'directMTF'))
    
    stopifnot(contype %in% c('legacy', 'uptakerate', 'totaluptake'))


    # XXX mj XXX
    # A partir de aquí todo es por modelo, hay que hacerlo para cada modelo.

    # # Modelos

#    mods <- makeOrganisms(models, initBiomass, biomassRxn, exclUptakeRxns, initConcentrations, contype, timeStep, logObj)
    mods <- makeOrganisms(models, initBiomass, biomassRxn, exclUptakeRxns, dynamicConstraints, deathrate)

    #mods <- addOrgsToMedium(mods, initConcentrations, contype)
    

    mod <- mods[[1]]

    # show log info
    if (logModel > length(mods)) 
        logModel <- length(mods)

    if (logModel==0) 
        logInfo2('Showing summary info table of all models')
    else 
        logInfo2('Showing info table of model', mods[[logModel]]$model@mod_id, '\n')

    logSubstrateFH(mod$model, mod$biomass, initConcentrations)

    # Preparamos el 'posible medio'. Tiene que contener todas las reacciones de intercambio de todos
    # los modelos. Todas las secreciones de metabolitos se realizan al medio y todos los consumos se
    # obtienen del medio, el cual es comun para todos. Posiblemente una funcion para obtener el medio

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
    
    # En principio no guardamos la matriz de flujos del medio, como siempre podemos acceder
    # a la suma. Prepara funcion para hacer la suma de matrices con distinto numero de filas
    # (reacciones de cada modelo) y distinto numero de columnas (Steps donde ha estado 
    # creciendo el modelo)
    # mediumFluxesMatrix <- NA
    
    timeVec <- 0

    for (stepNo in 1:nSteps){
        logInfo('\n','',type='')
        logInfo(logPrefix, 'Step', stepNo)
        # Los cambios de nutrientes en el sustrato son a nivel de 'paso', pertenecen al
        # medio, al entorno, y afectan a todos los modelos,
        # ver que hay que recalcular si se produce un cambio de  nutrientes
        if ( ! missing(nutrientChanges) && ! is.null(nutrientChanges) )  {
                medium <- doNutrientChanges(nutrientChanges, mods, medium, timeStep, stepNo)
            } # if there are nutrientChanges
        
        mediumStepBiomass <- 0
        mediumStepFluxes <- mediumStepFluxes0
        for (m in 1: length(mods)) {
            mod <- mods[[m]]

            ##-----------------------------------------------------------------------------##
            # Do the simulation

            ###########  XXX mj XXX  
            # Estas 'constraints' son por modelo, ver como se puede enfocar
            # Se puede llamar en cada paso, para todos los modelos, o solo el modelo que se quiere
            # modificar. 
            # if ((! missing(dynamicConstraints)) && (! is.null(dynamicConstraints))) {
            #     dc_model <- doDynamicConstraints(dynamicConstraints, mod$model, medium, sol, timeStep, stepNo)
            if ((! missing(dynamicConstraints)) && (! is.null(mod$dynamicConstraints))) {
                dc_model <- doDynamicConstraints(mod$dynamicConstraints, mod$model, medium, sol, timeStep, stepNo)
                #Ha habido cambios en el modelo? 
                # Si los cambios son en la funcion objetivo hay que recalcular el problema
                # es asi recalcular problema con los nuevos datos
                if (!isEqualModel(mod$model, dc_model)) {
                    if (!isEqualModelReactBounds(mod$model, dc_model, mod$biomassRxn))
                        mod$lpmod <- sybil::sysBiolAlg(dc_model, algorithm = "fba")
                    mod$model <- dc_model
                }
            } # if there are dynamicConstraints

            # probably a function for low bounds calcs, constrains and environment evolution

            # the model limits have changed
            #       so we need to upfate originalUptake to reflect the new limits
            ##        originalUptake[excReactInd] <- -lowbnd(model)[excReactInd]
            # this is where we do the update referred to earlier in the loop
            #modified <- lowbnd(model) != lowbnd(old.model)
            #originalUptake[modified & excReactInd] = -lowbnd(model)[modified & excReactInd]
            


            #uptakeBound <- getUptakeBound(mod$model, mod$concentrations, mod$biomass, timeStep)
            uptakeBound <- getUptakeBound(mod$model, mod$excReact, medium, mod$biomass, timeStep)
            # Las exclusiones de sustratos se realizaron al eliminarlos de excReactInd. Comprobar si hay alguna manera mas 'limpia'
            
            # Recalcular siempre???
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
            #  XXX mj XXX 
            # OJO!!!! se genera un archivo por cada modelo en cada paso, ademas hay que tener en cuenta 
            # que los modelos tienen miles de reacciones que hay que calcular en cada paso y almacenar en el fichero
            # puede ser una operación muuuyyyy costosa en tiempo
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
                # Sumamos la biomasa que queda a la biomasa del medio
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
            
            # if the objective function is not Biomass, this is plainly wrong!
            #mu_obj =  sol$obj;  ##objvalue sol.f
            # XXX j XXX any of these two alternatives should fix the issue
            #mu_bmidx <- sol$fluxes[mod$biomassIdx]
            #mu_bmrxn <- sol$fluxes[mod$model@react_id == mod$biomassRxn]
            #mu <- mu_bmidx

            # mu es el growth rate, el flujo de la reaccion de biomasa o crecimiento. 
            # Independientemente de cual sea la funcion objetivo para calcular el resultado del problema,
            # necesitamos cuanto crece la biomasa para, con ese dato, calcular el metabolito se ha consumido o excretado
            # varma-palsson 1994 equaciones (6) y (7)
            mu <- sol$fluxes[mod$biomassIdx]
            logInfo(logPrefix, 'Model:', mod$model@mod_id, 'Objective at step', stepNo, '(', stepNo * timeStep, 'h ):', mu)

            # get uptake fluxes after FBA/MTF
            uptakeFlux = sol$fluxes[mod$excReactIdx];
            # Update concentrations
            # XXX mj XXX we need previous biomass (Xo) to compute concentrations according
            # to function (7) of varma-palsson 1994
            # therefore new biomass must be computed later
            # keep track of medium changes
            medium[mod$excReactNames]= medium[mod$excReactNames] - uptakeFlux/mu*mod$biomass*(1-exp(mu*timeStep));
            # No substrate can be negative 
            #medium[medium < sybil::SYBIL_SETTINGS("TOLERANCE")] <- 0
            medium[medium < 0] <- 0

            # keep track of model concentrations changes
            mod$concentrations[mod$excReactIdx]= mod$concentrations[mod$excReactIdx] - uptakeFlux/mu*mod$biomass*(1-exp(mu*timeStep));
            # next line commented to keep adfba compatibility
            #mod$concentrations[mod$concentrations < sybil::SYBIL_SETTINGS("TOLERANCE")] <- 0;
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

            # Calcular datos acumulados de todos los modelos en este paso
            # Sumamos la biomasa a la biomasa del medio
            mediumStepBiomass <- mediumStepBiomass + mod$biomass

            modStepFluxes <- sol$fluxes[mod$substrateRxnsIdx]
            names(modStepFluxes) <- mod$model@react_id[mod$substrateRxnsIdx]

            # Al ser multimodelo puede que un modelo no tenga todos los metabolitos (reacciones) que se encuentran
            # en el substrato, por lo cual aparecera NA. Los omitimos.
            modStepFluxes <- na.omit(modStepFluxes)
            
            # Si omitimos, eliminamos algún valor, hay que tenerlo en cuenta para actualizar solo los
            # valores disponibles, para ello utilizamos names(mStepFluxes)
            mediumStepFluxes[names(modStepFluxes)] <- mediumStepFluxes[names(modStepFluxes)] + modStepFluxes

            # show log model step if selected
            if (m==logModel) logConcentFluxes(mod$biomass, mod$concentrations[mod$substrateRxnsIdx], sol$fluxes[mod$substrateRxnsIdx], stepNo)        

            mods[[m]] <- mod
        }

        # Check el estado de todos los modelos. Si todos han terminado no se puede continuar.
        # El codigo va aquí para no añadir un nuevo elemento a timeVec que no se podría realizar
        # ni ejecutar el posible cambio de nutrientes que está al inicio del bulcle para una segunda vuelta
        # Por compatibilidad de funcionamiento con la versión monomodelo.
        # Si aplicamos la tasa de muerte de los organismos hay que continuar con
        # la simulación mientras haya bioamsa y obtener los datos de defunción y acumulados

        # Es necesaria una variable para establecer el limite a partir del cual decidir que no queda biomasa?
        # de momento SYBIL_SETTINGS('TOLERANCE')
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
        
        # Aún no se ha terminado la simulación, añadimos un nuevo elemento al vector de tiempo para el siguiente paso      
        timeVec = c(timeVec,stepNo*timeStep);
    }# end loop
    # close progress bar

    ##--------------------------------------------------------------------------##
    # Simulation completed... prepare RETURN output

    # Preparamos el resultado con los datos de todo el entorno, los datos de evolucion de los substratos y
    # de crecimiento total de la biomasa del sistema
    resultsummary <- getResultSummary(mods, mediumMatrix, mediumBiomassVec, timeVec, retOptSol) 

    # Chequear si solo hay un modelo, en cuyo caso se devuelve los datos resumen, el objeto o la lista
    # para mantener la compatibilidad en el objeto devuelto por ADFBA. Los datos resumen
    # del sistema son los datos del unico modelo
    if (length(mods) == 1) return(resultsummary)

    # Hay mas de un modelo. Preparamos los datos de cada modelo y generamos una lista con los
    # resultados de todos los modelos por separado mas un resutado mas con el resumen de todo 
    # el sistema

    result <- lapply(mods, getResult, timeVec, fld, retOptSol, contype)

    return(append(result, resultsummary))
 
}

