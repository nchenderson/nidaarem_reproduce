
### Note that we only need to use a subset of the rows in each of the results.
setwd("~/Documents/nidaarem_reproduce")
#load("SimulationResults/Lasso/LassoRho8.RData")

nmethods <- 11
LassoResultsTable <- matrix(NA, nrow=2*nmethods, ncol=6)
rownames(LassoResultsTable) <- rep(c("PGD", "Nesterov", "Nesterov w restarts", 
                                "DAAREM (fm)", "DAAREM (am)", "NIDAAREM (fm)", "NIDAAREM (am)",
                                "RNIDAAREM (am)", "SDAAREM (am)", "SNIDAAREM (am)", "SQUAREM"), 2)
colnames(LassoResultsTable) <- c("Mean Iterations", "Median Iterations", "SD Iterations",
                                "Mean Timing", "Median Timing", "Mean Objfn Value")

design_setting <- 6
for(k in c(0, nmethods)) {
  
  if(k==0) {
      load("SimulationResults/Lasso/LassoRho8.RData")
  } else if(k==nmethods) {
      load("SimulationResults/Lasso/LassoRho95.RData")
  }
  ## Record PGD results
  LassoResultsTable[k+1,1] <- mean(PGDNI[,design_setting])
  LassoResultsTable[k+1,2] <- median(PGDNI[,design_setting])
  LassoResultsTable[k+1,3] <- sd(PGDNI[,design_setting])
  LassoResultsTable[k+1,4] <- mean(PGDTime[,design_setting])
  LassoResultsTable[k+1,5] <- median(PGDTime[,design_setting])
  LassoResultsTable[k+1,6] <- mean(PGDObj[,design_setting])

  ## Record Nesterov results
  LassoResultsTable[k+2,1] <- mean(NESTNI[,design_setting])
  LassoResultsTable[k+2,2] <- median(NESTNI[,design_setting])
  LassoResultsTable[k+2,3] <- sd(NESTNI[,design_setting])
  LassoResultsTable[k+2,4] <- mean(NESTTime[,design_setting])
  LassoResultsTable[k+2,5] <- median(NESTTime[,design_setting])
  LassoResultsTable[k+2,6] <- mean(NESTObj[,design_setting])

  ## Record Nesterov with restarts results
  LassoResultsTable[k+3,1] <- mean(RENESTNI[,design_setting])
  LassoResultsTable[k+3,2] <- median(RENESTNI[,design_setting])
  LassoResultsTable[k+3,3] <- sd(RENESTNI[,design_setting])
  LassoResultsTable[k+3,4] <- mean(RENESTTime[,design_setting])
  LassoResultsTable[k+3,5] <- median(RENESTTime[,design_setting])
  LassoResultsTable[k+3,6] <- mean(RENESTObj[,design_setting])

  ## Record DAAREM (with fixed monotonicity) results
  LassoResultsTable[k+4,1] <- mean(DAAREMNI[,design_setting])
  LassoResultsTable[k+4,2] <- median(DAAREMNI[,design_setting])
  LassoResultsTable[k+4,3] <- sd(DAAREMNI[,design_setting])
  LassoResultsTable[k+4,4] <- mean(DAAREMTime[,design_setting])
  LassoResultsTable[k+4,5] <- median(DAAREMTime[,design_setting])
  LassoResultsTable[k+4,6] <- mean(DAAREMObj[,design_setting])

  ## Record DAAREM (with alternating monotonicity) results
  LassoResultsTable[k+5,1] <- mean(DAARAMNI[,design_setting])
  LassoResultsTable[k+5,2] <- median(DAARAMNI[,design_setting])
  LassoResultsTable[k+5,3] <- sd(DAARAMNI[,design_setting])
  LassoResultsTable[k+5,4] <- mean(DAARAMTime[,design_setting])
  LassoResultsTable[k+5,5] <- median(DAARAMTime[,design_setting])
  LassoResultsTable[k+5,6] <- mean(DAARAMObj[,design_setting])

  ## Record NIDAAREM (with fixed monotonicity) results
  LassoResultsTable[k+6,1] <- mean(NIDAAREMNI[,design_setting])
  LassoResultsTable[k+6,2] <- median(NIDAAREMNI[,design_setting])
  LassoResultsTable[k+6,3] <- sd(NIDAAREMNI[,design_setting])
  LassoResultsTable[k+6,4] <- mean(NIDAAREMTime[,design_setting])
  LassoResultsTable[k+6,5] <- median(NIDAAREMTime[,design_setting])
  LassoResultsTable[k+6,6] <- mean(NIDAAREMObj[,design_setting])

  ## Record NIDAAREM (with alternating monotonicity) results
  LassoResultsTable[k+7,1] <- mean(NIDAARAMNI[,design_setting])
  LassoResultsTable[k+7,2] <- median(NIDAARAMNI[,design_setting])
  LassoResultsTable[k+7,3] <- sd(NIDAARAMNI[,design_setting])
  LassoResultsTable[k+7,4] <- mean(NIDAARAMTime[,design_setting])
  LassoResultsTable[k+7,5] <- median(NIDAARAMTime[,design_setting])
  LassoResultsTable[k+7,6] <- mean(NIDAARAMObj[,design_setting])

  ## Record NIDAAREM (with alternating monotonicity and residual 
  ## -based monotonicity monitoring) results
  LassoResultsTable[k+8,1] <- mean(RNIDAARAMNI[,design_setting])
  LassoResultsTable[k+8,2] <- median(RNIDAARAMNI[,design_setting])
  LassoResultsTable[k+8,3] <- sd(RNIDAARAMNI[,design_setting])
  LassoResultsTable[k+8,4] <- mean(RNIDAARAMTime[,design_setting])
  LassoResultsTable[k+8,5] <- median(RNIDAARAMTime[,design_setting])
  LassoResultsTable[k+8,6] <- mean(RNIDAARAMObj[,design_setting])

  ## Subsetted DAAREM with alternating monotonicity
  LassoResultsTable[k+9,1] <- mean(SDAARAMNI[,design_setting])
  LassoResultsTable[k+9,2] <- median(SDAARAMNI[,design_setting])
  LassoResultsTable[k+9,3] <- sd(SDAARAMNI[,design_setting])
  LassoResultsTable[k+9,4] <- mean(SDAARAMTime[,design_setting])
  LassoResultsTable[k+9,5] <- median(SDAARAMTime[,design_setting])
  LassoResultsTable[k+9,6] <- mean(SDAARAMObj[,design_setting])

  ## Subsetted NIDAAREM with alternating monotonicity
  LassoResultsTable[k+10,1] <- mean(SNIDAARAMNI[,design_setting])
  LassoResultsTable[k+10,2] <- median(SNIDAARAMNI[,design_setting])
  LassoResultsTable[k+10,3] <- sd(SNIDAARAMNI[,design_setting])
  LassoResultsTable[k+10,4] <- mean(SNIDAARAMTime[,design_setting])
  LassoResultsTable[k+10,5] <- median(SNIDAARAMTime[,design_setting])
  LassoResultsTable[k+10,6] <- mean(SNIDAARAMObj[,design_setting])

  ## Record SQUAREM results
  LassoResultsTable[k+11,1] <- mean(SQNI[,design_setting])
  LassoResultsTable[k+11,2] <- median(SQNI[,design_setting])
  LassoResultsTable[k+11,3] <- sd(SQNI[,design_setting])
  LassoResultsTable[k+11,4] <- mean(SQTime[,design_setting])
  LassoResultsTable[k+11,5] <- median(SQTime[,design_setting])
  LassoResultsTable[k+11,6] <- mean(SQObj[,design_setting])
}
save(LassoResultsTable, file="Tables/Table1LassoSimResults.RData")

library(xtable)
xtable(LassoResultsTable, digits=c(rep(1, 6), 6))




