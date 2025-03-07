library(nidaarem)
source("SimulationCode/LassoFunctions.R", chdir=TRUE)


for(k in 1:nreps) {
   fname <- paste("SimulationResults/Lasso/run1_n", tot_sampsize, ".RData", sep="")

   save(aa, file=fname)
}