#  zzz.R
# Inspired by sybil with help, as always, from:
#       https://adv-r.hadley.nz/index.html
#       stackoverflow (en este caso https://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r)
#       https://www.rdocumentation.org/
#

.MADFBAenv <- new.env(parent = emptyenv())

.onLoad <- function(lib, pkg) {

    # -------------------------------------------------------------- #
    # settings in MADFBA

    .MADFBAenv$settings <- list(
        MADFBA_VERSION      = "0.1",
        EXTINCTION          = 1e-3,
        MAXLOGFLUXES        = 4,
        VERBOSE             = 2
    )


    # -------------------------------------------------------------- #
    # other conf vars

    # .MADFBAenv$config <- 

}


#' @export
MADFBA_SETTINGS <- function(parm, value, ...) {

    if ( (missing(parm)) && (missing(value)) ) {
       return(.MADFBAenv$settings)
    }

    if (missing(value)) {
        if (!parm %in% names(.MADFBAenv$settings)) {
            stop("unknown parameter ", sQuote(parm))
        }
        return(.MADFBAenv$settings[[parm]])
    }
    
    if (length(parm) != 1) {
        stop("arguments 'parm' and 'value' must have a length of 1")
    }
    
    switch(parm,
        "MADFBA_VERSION" = {
        	stop("this value must not be set by the user!")
        },
    
        "VERBOSE" = {
            .MADFBAenv$settings[["VERBOSE"]] <- as.numeric(value)
        },
    
        "EXTINCTION" = {
            .MADFBAenv$settings[["EXTINCTION"]] <- as.numeric(value)
        },
    
        "MAXLOGFLUXES" = {
            .MADFBAenv$settings[["MAXLOGFLUXES"]] <- as.numeric(value)
        },
    
        {
            stop("unknown parameter: ", sQuote(parm))
        }
    )
}

