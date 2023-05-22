# MADFBA
The Multi-model Adaptive Dynamic FBA package is a free and open source software package that has been developed as part of a master's thesis from the National Center of Biotechnology-CSIC and works from an R environment. It is dependent on the Sybil package and its function is to perform dynamic, multi-model and adaptive flow balance analysis (FBA), i.e. it allows the user to modify the environment at any time during the simulation.

Inside the MADFBA package, in the "R" folder you can find the 8 scripts that compose the MADFBA algorithm.
### MADFBA.R
This is the main script, in charge of initializing all the necessary variables and executing the selected algorithm within the simulation. 
It is composed of the main function "MADFBA", whose parameters are mostly similar to those of the ADFBA algorithm developed at the CNB scientific computing service. Two new parameters to mention are "contype" and "deathrate", which will be defined at the beginning. The first one refers to the type of plots that will be presented at the end and there are three options: "legacy", "uptakerate" and "totaluptakerate". The second is the death rate, by definition it will be set to 0.
The main loop within this script is defined by a pre-specified maximum number of execution steps (nSteps), although it may terminate earlier if the nutrients in the medium are exhausted. In addition, the time that these steps last (timeSteps) is specified, the user can define it as he wants, although it is usually in hours. 
The way the loop works is as follows: over a series of steps and as long as the system can evolve, it will check for changes in the medium and update them. Then, the variables will be initialized and the dynamic constraints, if any, will be applied to each of the models. It is now when optimization occurs with the selected algorithm and solver. After this, a check is made to see if the model being run has exhausted all its metabolites, if so, it moves on to the next one and, in addition, the death rate is applied if it has been specified. If the substrate has not been exhausted, the growth rates are calculated, updated and the model states and fluxes, biomass, consumptions, etc. are saved. Finally, a summary of the whole procedure is saved.
### MADFBA_functions.R
Script with the main functions of the package defined.
### MADFBA_logfunctions.R
Script with the package's log functions, i.e. functions that present information.
### MADFBA_utils.R
Script with utility functions of the package.
### matlab.R
Script inspired by a function of the BacArena program that allows to convert matlab formatted models to modelorg so that they can be read by sybil.
### Optsol_dynamicFBA.R
Class taken from the sybilDynFBA package to obtain an optosol object compatible with sybilDynFBA. The purpose of this script is to be able to use the MADFBA package without having to depend on the sybilDynFBA package, which caused errors. It is not essential, but we want to maintain compatibility with the ADFBA package. Some modifications have been made to the method that displays the graphs.
### settings.R
Creation of the global variable MADFBA_SETTINGS inspired by SYBIL_SETTINGS.
### doDynamic.R
In this script the functions that produce variations during the execution of the algorithm are presented, that is, they add or remove nutrients from the medium or cause changes in the limits of the different reactions of the model.
The data folder contains text documents with examples of the tables used for nutrient changes and constraints and that can be used to play with the algorithm. It is very important that the researcher, when creating the tables and their modifications, write the titles of each column as shown in these examples.
