## Code to Generate the results of matrix completion simulation
## study where the penalty term is set to lambda = 20

setwd("~/Documents/nidaarem_reproduce")
load("SimulationResults/MatrixCompletion/MatrixCompleteSimFixedLambda20.RData")


MCResultsTable20 <- matrix(NA, nrow=8, ncol=6)
rownames(MCResultsTable20) <- c("PGD", "Nesterov", "Nesterov w restarts", "SQUAREM",
                                "DAAREM (cm)", "NIDAAREM (cm)", "RNIDAAREM (cm)", "SoftImpute")
colnames(MCResultsTable20) <- c("Mean Iterations", "Median Iterations", "SD Iterations",
                               "Mean Timing", "Median Timing", "Mean Objfn Value")

## Record PGD results
MCResultsTable20[1,1] <- mean(PGDNI)
MCResultsTable20[1,2] <- median(PGDNI)
MCResultsTable20[1,3] <- sd(PGDNI)
MCResultsTable20[1,4] <- mean(PGDTime)
MCResultsTable20[1,5] <- median(PGDTime)
MCResultsTable20[1,6] <- mean(PGDObj)

## Record Nesterov results
MCResultsTable20[2,1] <- mean(NESTNI)
MCResultsTable20[2,2] <- median(NESTNI)
MCResultsTable20[2,3] <- sd(NESTNI)
MCResultsTable20[2,4] <- mean(NESTTime)
MCResultsTable20[2,5] <- median(NESTTime)
MCResultsTable20[2,6] <- mean(NESTObj)

## Record Nesterov with restarts results
MCResultsTable20[3,1] <- mean(RENESTNI)
MCResultsTable20[3,2] <- median(RENESTNI)
MCResultsTable20[3,3] <- sd(RENESTNI)
MCResultsTable20[3,4] <- mean(RENESTTime)
MCResultsTable20[3,5] <- median(RENESTTime)
MCResultsTable20[3,6] <- mean(RENESTObj)

## Record SQUAREM results
MCResultsTable20[4,1] <- mean(SQNI)
MCResultsTable20[4,2] <- median(SQNI)
MCResultsTable20[4,3] <- sd(SQNI)
MCResultsTable20[4,4] <- mean(SQTime)
MCResultsTable20[4,5] <- median(SQTime)
MCResultsTable20[4,6] <- mean(SQObj)

## Record DAAREM (with cyclic monotonicity) results
MCResultsTable20[5,1] <- mean(DAARAMNI)
MCResultsTable20[5,2] <- median(DAARAMNI)
MCResultsTable20[5,3] <- sd(DAARAMNI)
MCResultsTable20[5,4] <- mean(DAARAMTime)
MCResultsTable20[5,5] <- median(DAARAMTime)
MCResultsTable20[5,6] <- mean(DAARAMObj)

## Record NIDAAREM (with cyclic monotonicity) results
MCResultsTable20[6,1] <- mean(NIDAARAMNI)
MCResultsTable20[6,2] <- median(NIDAARAMNI)
MCResultsTable20[6,3] <- sd(NIDAARAMNI)
MCResultsTable20[6,4] <- mean(NIDAARAMTime)
MCResultsTable20[6,5] <- median(NIDAARAMTime)
MCResultsTable20[6,6] <- mean(NIDAARAMObj)

## Record NIDAAREM (with cyclic monotonicity and residual 
## -based monotonicity monitoring) results
MCResultsTable20[7,1] <- mean(RNIDAARAMNI)
MCResultsTable20[7,2] <- median(RNIDAARAMNI)
MCResultsTable20[7,3] <- sd(RNIDAARAMNI)
MCResultsTable20[7,4] <- mean(RNIDAARAMTime)
MCResultsTable20[7,5] <- median(RNIDAARAMTime)
MCResultsTable20[7,6] <- mean(RNIDAARAMObj)

## Record SoftImpute results
MCResultsTable20[8,4] <- mean(SITime)
MCResultsTable20[8,5] <- median(SITime)
MCResultsTable20[8,6] <- mean(SIObj)


save(MCResultsTable20, file="Tables/MatrixCompletionLambda20SimTable.RData")

## To generate a table for Latex use the following code:
library(xtable)
xtable(MCResultsTable20[rownames(MCResultsTable20)!="MPE",], digits=c(rep(1, 6), 6))









