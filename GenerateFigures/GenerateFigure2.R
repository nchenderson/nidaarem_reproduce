load("SimulationResults/Lasso/LassoRho8.RData")

k1 <- 5


nreps <- 50
num_iters <- c(PGDNI[,k1], NESTNI[,k1], RENESTNI[,k1], SQNI[,k1],
               DAAREMNI[,k1], NIDAAREMNI[,k1], NIDAARAMNI[,k1], RNIDAARAMNI[,k1], 
               SDAARAMNI[,k1], SNIDAARAMNI[,k1])
method <- rep(c("PGD", "Nesterov", "Nesterov (restarts)", "SQUAREM", 
                "DAAREM", "NIDAAREM (fm)", "NIDAAREM (cm)", "NIDAAREM (cm-resid)",
                "SDAAREM (cm)", "SNIDAAREM (cm)"), each=nreps)

NumIterDat <- data.frame(method=method, numiters=num_iters)
NumIterDat$method <- factor(NumIterDat$method, levels=c("PGD", "Nesterov", "Nesterov (restarts)", "SQUAREM", 
                                                        "DAAREM", "NIDAAREM (fm)", "NIDAAREM (cm)", "NIDAAREM (cm-resid)",
                                                        "SDAAREM (cm)", "SNIDAAREM (cm)"))

boxplot(num_iters ~ method, data=NumIterDat, las=1, cex.axis=0.5, log="y")

