# NOTE: this function comes verbatim from package 'sybilDynFBA'
#        It is included here because it is not an exported function of
#        sybilDynFBA, but we need to use it to return a sybilDynFBA
#        compatible optsol object.
#
#        This source code belongs, therefore to the authors of 'sybilDynFBA'
#
# 1/6- Definition: to be added to file AllClasses.R

#------------------------------------------------------------------------------#
#                  definition of the class optsol_dynamicFBA                   #
#------------------------------------------------------------------------------#

#' Title
#'
#' @exportClass optsol_dynamicFBA
#' @import sybil methods
#' @importFrom graphics par layout points lines legend
#' @importFrom stats spline
#' @rdname optsol_dynamicFBA
#'
#' @slot concentrationMatrix matrix. 
#' @slot excRxnNames character. 
#' @slot timeVec numeric. 
#' @slot biomassVec numeric. 
#' @slot all_fluxes matrix. 
#'
#'
#' @noRd
setClass("optsol_dynamicFBA",
      representation(
        concentrationMatrix="matrix",  # Matrix of extracellular metabolite concentrations
        excRxnNames="character",       # Names of exchange reactions for the EC metabolites
        timeVec="numeric",             # Vector of time points
        biomassVec="numeric",          # Vector of biomass values
        all_fluxes="matrix"            # Matrix of all fluxes at all steps   24/7/2015
      ),
      contains = "optsol_optimizeProb",
      package = "sybil"
)

#showClass("optsol_dynamicFBA")


#-------------------------------------------------------------------------------#

# 2/6-Constructor: to be added to file AllClasses-constructors.R
# optsol_dynamicFBAClass
#' Title
#'
#' @param solver 
#' @param method 
#' @param nprob 
#' @param ncols 
#' @param nrows 
#' @param fld 
#' @param concmat 
#' @param exRxn 
#' @param tmVec 
#' @param bmVec 
#' @param all_fluxes 
#' @param ... 
#'
#' @export
#'
#' @noRd
optsol_dynamicFBA <- function(solver, method, nprob,
                #lpdir,
                ncols, nrows, 
                #objf,
                fld,concmat,exRxn,tmVec,bmVec,all_fluxes, ...) {
    if (missing(solver) || 
        missing(method) ||
        missing(nprob)  ||
        #missing(lpdir)  ||
        missing(ncols)  ||
        missing(nrows)  ||
        #missing(objf)   ||
        missing(fld)    ||
        missing(bmVec) ||
        missing(tmVec)  ||
        missing(all_fluxes)
       ) {
        stop("Not enough arguments for creating an object of class optsol_dynamicFBA!")
    }

    if (fld == TRUE) {
        fldist <- sybil::fluxDistribution(all_fluxes, ncols, nprob)
    }
    else {
        fldist <- sybil::fluxDistribution(NA)
    }

 
    new("optsol_dynamicFBA", 
        solver       = as.character(solver),
        method       = as.character(method),
        num_of_prob  = as.integer(nprob),
        lp_num_cols  = as.integer(ncols),
        lp_num_rows  = as.integer(nrows),
        lp_obj       = numeric(nprob),
        lp_ok        = integer(nprob),
        lp_stat      = integer(nprob),
        #lp_dir       = as.character(lpdir),
        #obj_function = as.character(objf),
        fluxdist     = fldist,
     #   num_of_steps_executed=nsteps,
        concentrationMatrix=concmat,
        excRxnNames=exRxn,
        timeVec= tmVec,  
        biomassVec = bmVec,
        all_fluxes = as.matrix(all_fluxes)
       )
}


#-----------------------------------------------------------------------------#

# 3/6-dynamicFBA: dynamicFBA.R
##              3.1 Check
##              3.2 optimizeProb -> get OptObj
##              3.3 Set Bounds
##              3.4 Call SimpleFBA
##              3.5 Store Solution
##              3.6 Adjust OUTPUT to class optsol_dynamicFBA

#-----------------------------------------------------------------------------#

# 4/6-  accessors: to be added to file optsol_dynamicFBA-accessors.R

#-----------------------------------------------------------------------------#

# 5/6-plot: to be added to file plot-methods.R

# optsol_dynamicFBAClass
# #' @export
# setGeneric("plot", function(x,y,
#                             ylim=50,
#                             xlab = "",
#                             ylab = "Value",
#                             type = "p",
#                             pch = 20,
#                             col = "black",             
#                             #               collower, colupper, pchupper, pchlower,
#                             #               dottedline = TRUE,
#                             plotRxns=NULL,
#                             baseline = 0,
#                             main = "Cell density",
#                             mainConcent = "Concentrations",
#                             ...) {standardGeneric("plot")})
#' @export
setMethod("plot", signature("optsol_dynamicFBA","missing"),
          function(x,y,
                   ylim=50,
                   xlab = "",
                   ylab = "Value",
                   type = "p",
                   pch = 20,
                   col = "black",             
    #               collower, colupper, pchupper, pchlower,
    #               dottedline = TRUE,
                   plotRxns=NULL,
                   baseline = 0,
                   main = "Cell density",
                   mainConcent = "Concentrations",
                   legend.pos = "topleft",
                    ...) {
                
                if(missing(plotRxns)){
                   plot(x@timeVec,x@biomassVec,main=main,xlab='Time',ylab=ylab, type="l");
                   points(x@timeVec,x@biomassVec, col = "red",lwd=2);
                   }
                else {  
                        # XXX mj XXX
                        # multimodelo
                        # Actualizar plotRxns para elimnar reacciones NO DISPONIBLES en el modelo
                        plotRxns <- x@excRxnNames[x@excRxnNames %in% plotRxns]
                        # continuamos con el codigo
                        def.par <- par(no.readonly = TRUE);
                        layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE))
                        #layout.show(2);
                        # first plot biomass
                         plot(x@timeVec,x@biomassVec, col = 1
                                            ,main=main,xlab='Time(hrs)',ylab="X(g/l)",type="l",lwd=2);
                         #lines(spline(x@timeVec,x@biomassVec, n = 201, method = "natural"), col = 2, type="c",lwd=2);
                         points(x@timeVec,x@biomassVec, col = "red",lwd=2);

                        # plot concentrations
                        ##plot(x@timeVec,2*x@biomassVec,main='Biomass',xlab='Time',ylab=ylab);
                        ## define min/max ()plot(x@timeVec,
                        ymin <- min(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
                        ymax <- max(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
                        for ( i in 1:length(plotRxns) ){
                            plotInd=(x@excRxnNames %in% plotRxns[i]);
                            #print( x@concentrationMatrix[plotInd]);
                            # multimodelo
                            # Si hay una reacción que que no se encuentra en este
                            # modelo, continuar a la siguiente reacción
                            # mas facil, actualizar plotRxns para que solo incluya las reacciones del modelo
                            # Así tambien actualizará en legend y no mostrará reacciones no disponibles.
                            if (!any(plotInd)) next
                            if(i==1){
                                plot(x@timeVec, x@concentrationMatrix[plotInd],
                                        type="l", col =i, main=mainConcent, ylab = "mmol",xlab='Time(hrs)'
                                        ,ylim=c(ymin,ymax));
                                #lines(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"), type="c", col =i);
                            }
                            else{
                                lines(x@timeVec, x@concentrationMatrix[plotInd], col =i, type='l');
                                #lines(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"), col =i, type='c');
                            }
                        }
                        legend(legend.pos, plotRxns, col=1:length(plotRxns), lty=1, lwd=1, cex=0.8);
                }
                  
                #if (!missing(plotRxns)){                   }
          }
)


#-----------------------------------------------------------------------------#

# 6/6-Documentation: optsol_dynamicFBA.Rd


#New version-----------------------------
