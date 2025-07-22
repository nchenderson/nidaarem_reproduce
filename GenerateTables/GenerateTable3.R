## Code to Generate the results of penalized logistic regression
## study with the madelon data and two fixed values of lambda

setwd("~/Documents/nidaarem_reproduce")
load("SimulationResults/LogisticLasso/MadelonFixedLambdas.RData")


nmethods <- 9
MadelonResultsTable <- matrix(NA, nrow=nmethods*2, ncol=7)
rownames(MadelonResultsTable) <- rep(c("PGD", "Nesterov", "Nesterov w restarts", "SQUAREM",
                              "DAAREM (fm)", "DAAREM (cm)", "NIDAAREM (fm)", "NIDAAREM (cm)", "MPE"), 2)
colnames(MadelonResultsTable) <- c("Lambda", "Mean Iterations", "Median Iterations", "SD Iterations",
                                   "Mean Timing", "Median Timing", "Mean Objfn Value")

MadelonResultsTable[,1] <- rep(c(67.6, 8.7), each=nmethods)

for(k in 0:1) {

## Record PGD results
MadelonResultsTable[k*nmethods + 1,2] <- mean(PGDNI[k+1,])
MadelonResultsTable[k*nmethods + 1,3] <- median(PGDNI[k+1,])
MadelonResultsTable[k*nmethods + 1,4] <- sd(PGDNI[k+1,])
MadelonResultsTable[k*nmethods + 1,5] <- mean(PGDTime[k+1,])
MadelonResultsTable[k*nmethods + 1,6] <- median(PGDTime[k+1,])
MadelonResultsTable[k*nmethods + 1,7] <- mean(PGDObj[k+1,])

## Record Nesterov results
MadelonResultsTable[k*nmethods + 2,2] <- mean(NESTNI[k+1,])
MadelonResultsTable[k*nmethods + 2,3] <- median(NESTNI[k+1,])
MadelonResultsTable[k*nmethods + 2,4] <- sd(NESTNI[k+1,])
MadelonResultsTable[k*nmethods + 2,5] <- mean(NESTTime[k+1,])
MadelonResultsTable[k*nmethods + 2,6] <- median(NESTTime[k+1,])
MadelonResultsTable[k*nmethods + 2,7] <- mean(NESTObj[k+1,])

## Record Nesterov with restarts results
MadelonResultsTable[k*nmethods + 3,2] <- mean(RENESTNI[k+1,])
MadelonResultsTable[k*nmethods + 3,3] <- median(RENESTNI[k+1,])
MadelonResultsTable[k*nmethods + 3,4] <- sd(RENESTNI[k+1,])
MadelonResultsTable[k*nmethods + 3,5] <- mean(RENESTTime[k+1,])
MadelonResultsTable[k*nmethods + 3,6] <- median(RENESTTime[k+1,])
MadelonResultsTable[k*nmethods + 3,7] <- mean(RENESTObj[k+1,])

## Record SQUAREM results
MadelonResultsTable[k*nmethods + 4,2] <- mean(SQNI[k+1,])
MadelonResultsTable[k*nmethods + 4,3] <- median(SQNI[k+1,])
MadelonResultsTable[k*nmethods + 4,4] <- sd(SQNI[k+1,])
MadelonResultsTable[k*nmethods + 4,5] <- mean(SQTime[k+1,])
MadelonResultsTable[k*nmethods + 4,6] <- median(SQTime[k+1,])
MadelonResultsTable[k*nmethods + 4,7] <- mean(SQObj[k+1,])

## Record DAAREM (with epsilon monotonicity) results
MadelonResultsTable[k*nmethods + 5,2] <- mean(DAAREMNI[k+1,])
MadelonResultsTable[k*nmethods + 5,3] <- median(DAAREMNI[k+1,])
MadelonResultsTable[k*nmethods + 5,4] <- sd(DAAREMNI[k+1,])
MadelonResultsTable[k*nmethods + 5,5] <- mean(DAAREMTime[k+1,])
MadelonResultsTable[k*nmethods + 5,6] <- median(DAAREMTime[k+1,])
MadelonResultsTable[k*nmethods + 5,7] <- mean(DAAREMObj[k+1,])

## Record DAAREM (with cyclic monotonicity) results
MadelonResultsTable[k*nmethods + 6,2] <- mean(DAARAMNI[k+1,])
MadelonResultsTable[k*nmethods + 6,3] <- median(DAARAMNI[k+1,])
MadelonResultsTable[k*nmethods + 6,4] <- sd(DAARAMNI[k+1,])
MadelonResultsTable[k*nmethods + 6,5] <- mean(DAARAMTime[k+1,])
MadelonResultsTable[k*nmethods + 6,6] <- median(DAARAMTime[k+1,])
MadelonResultsTable[k*nmethods + 6,7] <- mean(DAARAMObj[k+1,])

## Record NIDAAREM (with epsilon monotonicity) results
MadelonResultsTable[k*nmethods + 7,2] <- mean(NIDAAREMNI[k+1,])
MadelonResultsTable[k*nmethods + 7,3] <- median(NIDAAREMNI[k+1,])
MadelonResultsTable[k*nmethods + 7,4] <- sd(NIDAAREMNI[k+1,])
MadelonResultsTable[k*nmethods + 7,5] <- mean(NIDAAREMTime[k+1,])
MadelonResultsTable[k*nmethods + 7,6] <- median(NIDAAREMTime[k+1,])
MadelonResultsTable[k*nmethods + 7,7] <- mean(NIDAAREMObj[k+1,])

## Record NIDAAREM (with cyclic monotonicity) results
MadelonResultsTable[k*nmethods + 8,2] <- mean(NIDAARAMNI[k+1,])
MadelonResultsTable[k*nmethods + 8,3] <- median(NIDAARAMNI[k+1,])
MadelonResultsTable[k*nmethods + 8,4] <- sd(NIDAARAMNI[k+1,])
MadelonResultsTable[k*nmethods + 8,5] <- mean(NIDAARAMTime[k+1,])
MadelonResultsTable[k*nmethods + 8,6] <- median(NIDAARAMTime[k+1,])
MadelonResultsTable[k*nmethods + 8,7] <- mean(NIDAARAMObj[k+1,])

## Record MPE Results
MadelonResultsTable[k*nmethods + 9,2] <- mean(MPENI[k+1,])
MadelonResultsTable[k*nmethods + 9,3] <- median(MPENI[k+1,])
MadelonResultsTable[k*nmethods + 9,4] <- sd(MPENI[k+1,])
MadelonResultsTable[k*nmethods + 9,5] <- mean(MPETime[k+1,])
MadelonResultsTable[k*nmethods + 9,6] <- median(MPETime[k+1,])
MadelonResultsTable[k*nmethods + 9,7] <- mean(MPEObj[k+1,])

}
## Multiply last column by negative 1 to reflect the minimization
## version of this optimization problem
MadelonResultsTable[,7] <- (-1)*MadelonResultsTable[,7]

save(MadelonResultsTable, file="Tables/MadelonTwoFixedLambdas.RData")

## To generate a table for Latex use the following code:
#library(xtable)
#xtable(MadelonResultsTable[,-1], digits=c(rep(1, 4), 2, 2, 6))









