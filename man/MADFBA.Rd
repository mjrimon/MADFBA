% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MADFBA.R
\name{MADFBA}
\alias{MADFBA}
\title{Simulacion Multimodelo Adaptativo Dinámico FBA}
\usage{
MADFBA(
  models,
  substrateRxns,
  initConcentrations,
  initBiomass,
  timeStep,
  nSteps,
  exclUptakeRxns = c("EX_co2(e)", "EX_o2(e)", "EX_h2o(e)", "EX_h(e)"),
  retOptSol = TRUE,
  fld = FALSE,
  verboseMode = 2,
  biomassRxn = "",
  dynamicConstraints = NULL,
  nutrientChanges = NULL,
  method = "FBA",
  logModel = 0,
  contype = "legacy",
  deathrate = 0,
  ...
)
}
\arguments{
\item{models}{Un modelo o un lista de modelos de tipo sybil::modelorg}

\item{substrateRxns}{Nombre de los substratos de las reacciones disponibles en el entorno}

\item{initConcentrations}{Concentracion inicial de los substratos del entorno}

\item{initBiomass}{Un entero o un vector con la biomasa inicial de los modelos}

\item{timeStep}{Paso de tiempo de simulación (horas)}

\item{nSteps}{Numero de pasos que tendrá la simulación como máximo (puede terminar antes)}

\item{exclUptakeRxns}{Reacciones de consumo que se excluye de la simulacion. si no se especifica, por defecto se excluyen 'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'}

\item{retOptSol}{Establece el tipo del resultado: TRUE  un objeto optsol_dynamicFBA; FALSE una lista}

\item{fld}{Indica si se deben incluir en el resultado todos los flujos de todos los pasos de la simulación}

\item{verboseMode}{Nivel de infrmación presentada}

\item{biomassRxn}{Nombre de la reacción de crecimiento del modelo, de la biomasa, si no se especifica se obtiene por defecto}

\item{dynamicConstraints}{}

\item{nutrientChanges}{}

\item{method}{}

\item{logModel}{Modelo del que mostrar el log durante la simulación, por defecto el 0, la suma de todos}

\item{contype}{Tipo de datos de concentraciones consumidas/excretadas por modelo que se guardarán en cada paso:
'legacy' concentración de substratos aparente, incializado a concentraciones iniciales de substrato.
'uptakerate'  cantidad de substrato consumido/excretado en ese instante.
'totaluptake' cantidad total de substrato consumido/excretado hasta ese instante.}

\item{deathrate}{deathrate de los modelos
Si es distinto de 0, valor por defecto, aplicar como tasa de muerte del organismo al quedarse sin nutrientes}

\item{...}{}
}
\value{
Un objeto optsol_dynamicFBA o una lista con el resultado de un modelo o una lista de objetos optsol_dynamicFBA o listas si la simulación es multimodelo
}
\description{
Ejecuta una simulación de la evolución del ensayo con uno o varios modelos
de acuerdo a los parámetros especificados
}
\examples{

# Load models data
data(Ec_core)
data(Slividans)

# set the models uptake bounds
Ec_core <- sybil::changeBounds(Ec_core, react = "EX_glc(e)", lb = -12)
Ec_core <- sybil::changeBounds(Ec_core, react = "EX_o2(e)",  lb = -10)
Ec_core <- sybil::changeBounds(Ec_core, react = "EX_ac(e)",  lb = -10)

Slividans <- sybil::changeBounds(Slividans, react = "EX_mnl(e)", lb = -7)
Slividans <- sybil::changeBounds(Slividans, react = "EX_glc(e)", lb = -5)

# define medium, substrate names and initial concentrations
init.source <- c("EX_ac(e)","EX_o2(e)","EX_glc(e)", 'EX_mnl(e)')
init.conc   <- c(10,1,28,4)
# set models and initial biomass for simulation
init.models <- c(Ec_core,Slividans)
init.bmass <- c(.01, .0075)
Ec_df_min  <- MADFBA(init.models,  #exclUptakeRxns = c(),  # Si no se especifican exclusiones por defecto exluye, al menos o2
                              substrateRxns      = init.source,
                              initConcentrations = init.conc,
                              initBiomass        = init.bmass,
                              timeStep=.1,nSteps=200,verbose=3,
                              method = "FBA",                   # "directMTF" #"FBA" #"MTF"  #"directFBA"
                              fld=TRUE,                         # All fluxes para poder representarlos luego si procede
                              retOptSol = TRUE
                              )
}
