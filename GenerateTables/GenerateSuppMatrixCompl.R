## Code to Generate the results of matrix completion simulation
## study where the penalty term is set to lambda = 20.
## Results from this simulation study will be in the supplement

setwd("~/Documents/nidaarem_reproduce")
load("SimulationResults/MatrixCompletion/MatrixCompleteSimFixedLambda200.RData")


MCResultsTable200 <- matrix(NA, nrow=8, ncol=6)
rownames(MCResultsTable200) <- c("PGD", "Nesterov", "Nesterov w restarts", "SQUAREM",
                                "DAAREM (cm)", "NIDAAREM (cm)", "RNIDAAREM (cm)", "SoftImpute")
colnames(MCResultsTable200) <- c("Mean Iterations", "Median Iterations", "SD Iterations",
                                "Mean Timing", "Median Timing", "Mean Objfn Value")

## Record PGD results
MCResultsTable200[1,1] <- mean(PGDNI)
MCResultsTable200[1,2] <- median(PGDNI)
MCResultsTable200[1,3] <- sd(PGDNI)
MCResultsTable200[1,4] <- mean(PGDTime)
MCResultsTable200[1,5] <- median(PGDTime)
MCResultsTable200[1,6] <- mean(PGDObj)

## Record Nesterov results
MCResultsTable200[2,1] <- mean(NESTNI)
MCResultsTable200[2,2] <- median(NESTNI)
MCResultsTable200[2,3] <- sd(NESTNI)
MCResultsTable200[2,4] <- mean(NESTTime)
MCResultsTable200[2,5] <- median(NESTTime)
MCResultsTable200[2,6] <- mean(NESTObj)

## Record Nesterov with restarts results
MCResultsTable200[3,1] <- mean(RENESTNI)
MCResultsTable200[3,2] <- median(RENESTNI)
MCResultsTable200[3,3] <- sd(RENESTNI)
MCResultsTable200[3,4] <- mean(RENESTTime)
MCResultsTable200[3,5] <- median(RENESTTime)
MCResultsTable200[3,6] <- mean(RENESTObj)

## Record SQUAREM results
MCResultsTable200[4,1] <- mean(SQNI)
MCResultsTable200[4,2] <- median(SQNI)
MCResultsTable200[4,3] <- sd(SQNI)
MCResultsTable200[4,4] <- mean(SQTime)
MCResultsTable200[4,5] <- median(SQTime)
MCResultsTable200[4,6] <- mean(SQObj)

## Record DAAREM (with cyclic monotonicity) results
MCResultsTable200[5,1] <- mean(DAARAMNI)
MCResultsTable200[5,2] <- median(DAARAMNI)
MCResultsTable200[5,3] <- sd(DAARAMNI)
MCResultsTable200[5,4] <- mean(DAARAMTime)
MCResultsTable200[5,5] <- median(DAARAMTime)
MCResultsTable200[5,6] <- mean(DAARAMObj)

## Record NIDAAREM (with cyclic monotonicity) results
MCResultsTable200[6,1] <- mean(NIDAARAMNI)
MCResultsTable200[6,2] <- median(NIDAARAMNI)
MCResultsTable200[6,3] <- sd(NIDAARAMNI)
MCResultsTable200[6,4] <- mean(NIDAARAMTime)
MCResultsTable200[6,5] <- median(NIDAARAMTime)
MCResultsTable200[6,6] <- mean(NIDAARAMObj)

## Record NIDAAREM (with cyclic monotonicity and residual 
## -based monotonicity monitoring) results
MCResultsTable200[7,1] <- mean(RNIDAARAMNI)
MCResultsTable200[7,2] <- median(RNIDAARAMNI)
MCResultsTable200[7,3] <- sd(RNIDAARAMNI)
MCResultsTable200[7,4] <- mean(RNIDAARAMTime)
MCResultsTable200[7,5] <- median(RNIDAARAMTime)
MCResultsTable200[7,6] <- mean(RNIDAARAMObj)

## Record SoftImpute results
MCResultsTable200[8,4] <- mean(SITime)
MCResultsTable200[8,5] <- median(SITime)
MCResultsTable200[8,6] <- mean(SIObj)


save(MCResultsTable200, file="Tables/MatrixCompletionLambda200SimTable.RData")

## To generate a table for Latex use the following code:
library(xtable)
xtable(MCResultsTable200[rownames(MCResultsTable200)!="MPE",], digits=c(rep(1, 6), 6))









