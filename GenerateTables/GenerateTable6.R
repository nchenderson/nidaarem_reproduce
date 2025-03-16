## Code to Generate the results of the constrained QP simulation study

setwd("~/Documents/nidaarem_reproduce")
load("SimulationResults/ConstrainedQPSimResults.RData")


CQPResultsTable <- matrix(NA, nrow=8, ncol=6)
rownames(CQPResultsTable) <- c("PGD", "Nesterov", "Nesterov w restarts", "SQUAREM",
                               "MPE", "DAAREM (cm)", "NIDAAREM", "NIDAAREM (cm)")
colnames(CQPResultsTable) <- c("Mean Iterations", "Median Iterations", "SD Iterations",
                               "Mean Timing", "Median Timing", "Mean Objfn Value")

## Record PGD results
CQPResultsTable[1,1] <- mean(PGDNI)
CQPResultsTable[1,2] <- median(PGDNI)
CQPResultsTable[1,3] <- sd(PGDNI)
CQPResultsTable[1,4] <- mean(PGDTime)
CQPResultsTable[1,5] <- median(PGDTime)
CQPResultsTable[1,6] <- mean(PGDObj)

## Record Nesterov results
CQPResultsTable[2,1] <- mean(NESTNI)
CQPResultsTable[2,2] <- median(NESTNI)
CQPResultsTable[2,3] <- sd(NESTNI)
CQPResultsTable[2,4] <- mean(NESTTime)
CQPResultsTable[2,5] <- median(NESTTime)
CQPResultsTable[2,6] <- mean(NESTObj)

## Record Nesterov with restarts results
CQPResultsTable[3,1] <- mean(RENESTNI)
CQPResultsTable[3,2] <- median(RENESTNI)
CQPResultsTable[3,3] <- sd(RENESTNI)
CQPResultsTable[3,4] <- mean(RENESTTime)
CQPResultsTable[3,5] <- median(RENESTTime)
CQPResultsTable[3,6] <- mean(RENESTObj)

## Record SQUAREM results
CQPResultsTable[4,1] <- mean(SQNI)
CQPResultsTable[4,2] <- median(SQNI)
CQPResultsTable[4,3] <- sd(SQNI)
CQPResultsTable[4,4] <- mean(SQTime)
CQPResultsTable[4,5] <- median(SQTime)
CQPResultsTable[4,6] <- mean(SQObj)

## Record MPE results
CQPResultsTable[5,1] <- mean(MPENI)
CQPResultsTable[5,2] <- median(MPENI)
CQPResultsTable[5,3] <- sd(MPENI)
CQPResultsTable[5,4] <- mean(MPETime)
CQPResultsTable[5,5] <- median(MPETime)
CQPResultsTable[5,6] <- mean(MPEObj)

## Record DAAREM (with cyclic monotonicity) results
CQPResultsTable[6,1] <- mean(DAARAMNI)
CQPResultsTable[6,2] <- median(DAARAMNI)
CQPResultsTable[6,3] <- sd(DAARAMNI)
CQPResultsTable[6,4] <- mean(DAARAMTime)
CQPResultsTable[6,5] <- median(DAARAMTime)
CQPResultsTable[6,6] <- mean(DAARAMObj)

## Record NIDAAREM results
CQPResultsTable[7,1] <- mean(NIDAAREMNI)
CQPResultsTable[7,2] <- median(NIDAAREMNI)
CQPResultsTable[7,3] <- sd(NIDAAREMNI)
CQPResultsTable[7,4] <- mean(NIDAAREMTime)
CQPResultsTable[7,5] <- median(NIDAAREMTime)
CQPResultsTable[7,6] <- mean(NIDAAREMObj)

## Record NIDAAREM (with cyclic monotonicity) results
CQPResultsTable[8,1] <- mean(NIDAARAMNI)
CQPResultsTable[8,2] <- median(NIDAARAMNI)
CQPResultsTable[8,3] <- sd(NIDAARAMNI)
CQPResultsTable[8,4] <- mean(NIDAARAMTime)
CQPResultsTable[8,5] <- median(NIDAARAMTime)
CQPResultsTable[8,6] <- mean(NIDAARAMObj)


#save(CQPResultsTable, file="Tables/ConstrQuadProgramSimTable.RData")

## library(xtable)
## xtable(CQPResultsTable)





  



