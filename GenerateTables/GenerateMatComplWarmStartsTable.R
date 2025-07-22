## Code to Generate the table of the results of matrix completion
## using the movielens data and all 10 penalty values
## using "warm starts"

setwd("~/Documents/nidaarem_reproduce")
load("SimulationResults/MatrixCompletion/MatrixCompleteWarmStarts.RData")

nmethods <- 8
MatrixComplWMResultsTable <- matrix(NA, nrow=nmethods, ncol=2)
rownames(MatrixComplWMResultsTable) <- c("PGD", "Nesterov", "Nesterov w restarts", "SQUAREM",
                                     "DAAREM (fm)", "NIDAAREM (cm)", "NIDAAREM (cm-resid)", "Soft-Impute")
colnames(MatrixComplWMResultsTable) <- c("Number of Iterations", "Timing")

## Record PGD results
MatrixComplWMResultsTable[1,1] <- sum(PGDNI)
MatrixComplWMResultsTable[1,2] <- sum(PGDTime)

## Record Nesterov results
MatrixComplWMResultsTable[2,1] <- sum(NESTNI)
MatrixComplWMResultsTable[2,2] <- sum(NESTTime)

## Record Nesterov with restarts results
MatrixComplWMResultsTable[3,1] <- sum(RENESTNI)
MatrixComplWMResultsTable[3,2] <- sum(RENESTTime)

## Record SQUAREM results
MatrixComplWMResultsTable[4,1] <- sum(SQNI)
MatrixComplWMResultsTable[4,2] <- sum(SQTime)

## Record DAAREM (with epsilon monotonicity) results
MatrixComplWMResultsTable[5,1] <- sum(DAAREMNI)
MatrixComplWMResultsTable[5,2] <- sum(DAAREMTime)

## Record NIDAAREM (with cyclic monotonicity) results
MatrixComplWMResultsTable[6,1] <- sum(NIDAARAMNI)
MatrixComplWMResultsTable[6,2] <- sum(NIDAARAMTime)

## Record RINIDAAREM (with residual monitoring) results
MatrixComplWMResultsTable[7,1] <- sum(RNIDAARAMNI)
MatrixComplWMResultsTable[7,2] <- sum(RNIDAARAMTime)

## Record RINIDAAREM (with residual monitoring) results
MatrixComplWMResultsTable[8,2] <- sum(SITime)


save(MatrixComplWMResultsTable, file="Tables/MatrixCompletionWarmStartsResults.RData")

## To generate a table for Latex use the following code:
library(xtable)
xtable(MatrixComplWMResultsTable, digits=c(0, 0, 1))









