#!/usr/bin/env Rscript

library("ape")
library("BactDating")

#two arguments the prefix and the mcmc chain length
args = commandArgs(trailingOnly=TRUE)

#mcmc chain length
nbIts <- as.numeric(args[2])
cluster_tree <-loadGubbins(paste("bacdacter_files/",args[1], sep=""))
dates <- read.csv(paste("bacdacter_files/", args[1], "_dates.csv", sep=""), header = TRUE, sep =",", stringsAsFactors = FALSE)
rownames(dates) <- dates$ID
cluster_dates <- dates[cluster_tree$tip.label,2]
results <- bactdate(tree=cluster_tree, date=cluster_dates, initMu = NA, initAlpha = NA, initSigma = NA,
                            updateMu = T, updateAlpha = T, updateSigma = T, updateRoot = T,
                            nbIts = as.numeric(args[2]), thin = ceiling(10000), useCoalPrior = T,
                            model = "mixedgamma", useRec = T, showProgress = F)

jpeg(paste("bacdacter_out/",args,"_bactdate_plot_r2.jpg", sep=""))
plot(results,"trace")
dev.off()

save(results, file = paste("bacdacter_out/",args[1], "_bactdate.results_r2.RData", sep=""))
