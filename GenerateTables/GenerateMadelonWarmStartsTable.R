## Code to Generate the results of penalized logistic regression
## study with the madelon data and all 10 penalty values
## using "warm starts"

setwd("~/Documents/nidaarem_reproduce")
load("SimulationResults/LogisticLasso/MadelonFullPath.RData")


nmethods <- 9
MadelonWMResultsTable <- matrix(NA, nrow=nmethods, ncol=2)
rownames(MadelonWMResultsTable) <- c("PGD", "Nesterov", "Nesterov w restarts", "SQUAREM",
                                       "DAAREM (fm)", "DAAREM (cm)", "NIDAAREM (fm)", "NIDAAREM (cm)", "MPE")
colnames(MadelonWMResultsTable) <- c("Number of Iterations", "Timing")


  ## Record PGD results
  MadelonWMResultsTable[1,1] <- sum(PGDNI)
  MadelonWMResultsTable[1,2] <- sum(PGDTime)
 
  ## Record Nesterov results
  MadelonWMResultsTable[2,1] <- sum(NESTNI)
  MadelonWMResultsTable[2,2] <- sum(NESTTime)
 
  ## Record Nesterov with restarts results
  MadelonWMResultsTable[3,1] <- sum(RENESTNI)
  MadelonWMResultsTable[3,2] <- sum(RENESTTime)
 
  ## Record SQUAREM results
  MadelonWMResultsTable[4,1] <- sum(SQNI)
  MadelonWMResultsTable[4,2] <- sum(SQTime)
 
  ## Record DAAREM (with epsilon monotonicity) results
  MadelonWMResultsTable[5,1] <- sum(DAAREMNI)
  MadelonWMResultsTable[5,2] <- sum(DAAREMTime)
 
  ## Record DAAREM (with cyclic monotonicity) results
  MadelonWMResultsTable[6,1] <- sum(DAARAMNI)
  MadelonWMResultsTable[6,2] <- sum(DAARAMTime)
 
  ## Record NIDAAREM (with epsilon monotonicity) results
  MadelonWMResultsTable[7,1] <- sum(NIDAAREMNI)
  MadelonWMResultsTable[7,2] <- sum(NIDAAREMTime)
 
  ## Record NIDAAREM (with cyclic monotonicity) results
  MadelonWMResultsTable[8,1] <- sum(NIDAARAMNI)
  MadelonWMResultsTable[8,2] <- sum(NIDAARAMTime)
 
  ## Record MPE Results
  MadelonWMResultsTable[9,1] <- sum(MPENI)
  MadelonWMResultsTable[9,2] <- sum(MPETime)


save(MadelonWMResultsTable, file="Tables/MadelonWarmStartsResults.RData")

## To generate a table for Latex use the following code:
library(xtable)
xtable(MadelonWMResultsTable, digits=c(0, 0, 1))









