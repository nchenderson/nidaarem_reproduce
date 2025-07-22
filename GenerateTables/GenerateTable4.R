## Code to Generate the results of matrix completion simulation
## where matrix completion is performed over a series
## of different penalty terms using "warm starts".

setwd("~/Documents/nidaarem_reproduce")
load("SimulationResults/MatrixCompletion/MatrixCompleteWarmStarts.RData")


MCWarmStartsTable <- matrix(NA, nrow=8, ncol=2)
rownames(MCWarmStartsTable) <- c("PGD", "Nesterov", "Nesterov w restarts", "SQUAREM",
                              "DAAREM (cm)", "NIDAAREM (cm)", "RNIDAAREM (cm)", "SoftImpute")
colnames(MCWarmStartsTable) <- c("Total Iterations", "Total Time")

## Record PGD results
MCWarmStartsTable[1,1] <- sum(PGDNI)
MCWarmStartsTable[1,2] <- sum(PGDTime)

## Record Nesterov results
MCWarmStartsTable[2,1] <- sum(NESTNI)
MCWarmStartsTable[2,2] <- sum(NESTTime)

## Record Nesterov with restarts results
MCWarmStartsTable[3,1] <- sum(RENESTNI)
MCWarmStartsTable[3,2] <- sum(RENESTTime)

## Record SQUAREM results
MCWarmStartsTable[4,1] <- sum(SQNI)
MCWarmStartsTable[4,2] <- sum(SQTime)

## Record DAAREM (with cyclic monotonicity) results
MCWarmStartsTable[5,1] <- sum(DAAREMNI)
MCWarmStartsTable[5,2] <- sum(DAAREMTime)

## Record NIDAAREM (with cyclic monotonicity) results
MCWarmStartsTable[6,1] <- sum(NIDAARAMNI)
MCWarmStartsTable[6,2] <- sum(NIDAARAMTime)

## Record NIDAAREM (with cyclic monotonicity and residual 
## -based monotonicity monitoring) results
MCWarmStartsTable[7,1] <- sum(RNIDAARAMNI)
MCWarmStartsTable[7,2] <- sum(RNIDAARAMTime)

## Record SoftImpute results
MCWarmStartsTable[8,2] <- sum(SITime)


#save(MCWarmStartsTable, file="Tables/MatrixCompletionWarmStartsTable.RData")

## To generate a table for Latex use the following code:
#library(xtable)
#xtable(MCWarmStartsTable, digits=c(0,0,1))









