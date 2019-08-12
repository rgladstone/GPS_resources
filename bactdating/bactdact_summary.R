library(BactDating)
library(ape)
library(coda)

setwd("/Users/rg9/rg9_documents/GPS/MGEN_bioresource/bactdating/")

#function to convert bactdate mcmc to coda format
as.mcmc.resBactDating <- function(x,burnin=0.1) {
  record=x$record
  record=record[max(1,round(nrow(record)*burnin)):nrow(record),]
  mat=cbind(
    record[,'mu'],
    record[,'sigma'],
    record[,'alpha'])
  colnames(mat)<-c('mu','sigma','alpha')
  return(coda::as.mcmc(mat))
}

#function to compare models 
modelcompare_out = function(res1,res2) {
  dic1=res1$dic
  dic2=res2$dic
  dif=dic2-dic1
  if (dif>10){
    res_mod <- 'Model 1 is definitely better'
  }
  if (dif>5 && dif<10){
    res_mod <- 'Model 1 is slightly better'
  }
  if (abs(dif)<5) {
    res_mod <- 'The difference is not significant'
  }
  if (dif< -5 && dif>-10) {
    res_mod <- 'Model 1 is slightly worse'
  }
  if (dif< -10) {
    res_mod <- 'Model 1 is definitely worse'
  }
return(c(dic1,dic2,res_mod))
}

#loop through lineages and extract stats to assess dating 

GPSCs <- c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,30,31,32,33,34,37,38,39,40,41,43,47,48,50,51,52,53,54,55,56,57,58,61,62,65,67,68,70,72,76,77,78,79,80,81,90,91,93,94,97,103,105,108,115,117,131)
out <- matrix(nrow=0,ncol=11)

for(i in GPSCs){
    #1st iter load and convert
    tmp <- load(file = paste("GPSC", i, "_bactdate.results.RData", sep = ""))
    results_r1 <- get(tmp)
    coda_res1 <- as.mcmc.resBactDating(results_r1, burnin = 0.1)
    if(dim(coda_res1)[1]==9001){
      coda_res1 <- as.mcmc(coda_res1[seq(1,9001,10),])
    }
    #write out tree
    write.tree(results_r1$tree, file =  paste("GPSC", i, "_BD.tre", sep = ""))
    
    #2nd iter load and convert
    tmp <- load(file = paste("GPSC", i, "_bactdate.results_r2.RData", sep = ""))
    results_r2 <- get(tmp)
    coda_res2 <- as.mcmc.resBactDating(results_r2, burnin = 0.1)
    if(dim(coda_res2)[1]==9001){
      coda_res2 <- as.mcmc(coda_res2[seq(1,9001,10),])
    }
    #3rd iter load and convert
    tmp <- load(file = paste("GPSC", i, "_bactdate.results_r3.RData", sep = ""))
    results_r3 <- get(tmp)
    coda_res3 <- as.mcmc.resBactDating(results_r3, burnin = 0.1)
    if(dim(coda_res3)[1]==9001){
      coda_res3 <- as.mcmc(coda_res3[seq(1,9001,10),])
    }
    #combine iterations
    coda_iters <- mcmc.list(coda_res1,coda_res2,coda_res3)
    
    #rand dates iter load
    filename <- paste("GPSC", i, "_bactdate.results_rand_dates.RData", sep = "")
    if(file.exists(filename)){
      tmp <- load(file = paste("GPSC", i, "_bactdate.results_rand_dates.RData", sep = ""))
      results_rand_dates <- get(tmp)
    } else {
      results_rand_dates <- "NA"
    }
    
    #extract and store stats
    gpsc_out <- matrix(nrow = 1, ncol = 11)
    gpsc_out[1,1] <- paste("GPSC",i,sep = "")
    gpsc_out[1,2] <- results_r1$rootprob
    loop <- 3
    for (nam in c('likelihood','prior','mu','sigma','alpha')) {
      v=results_r1$record[,nam]
      v=v[(1+length(v)/2):length(v)]
      v=sort(v)
      vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])
      gpsc_out[1,loop] <- as.character(sprintf('%.2e [%.2e;%.2e],',vals[1],vals[2],vals[3]))
      loop <- loop+1
    }
    v=results_r1$record[,Ntip(results_r1$tree)+1]
    v=sort(v[(1+length(v)/2):length(v)])
    vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])
    gpsc_out[1,8] <- as.character(sprintf('%.2f [%.2f;%.2f]',vals[1],vals[2],vals[3]))
    gpsc_out[1,9] <- as.character(sprintf('%.2f; %.2f; %.2f',effectiveSize(coda_res1)[1],effectiveSize(coda_res1)[2],effectiveSize(coda_res1)[3]))
    if(results_rand_dates=="NA"){
      model_test1 <- "NA;NA;NA"
    } else {
      model_test1 <- modelcompare_out(results_r1,results_rand_dates)
    }
    gpsc_out[1,10] <- as.character(sprintf('%s; %s; %s',model_test1[1],model_test1[2],model_test1[3]))
        converg_test <-gelman.diag(coda_iters, confidence = 0.95, transform=FALSE, autoburnin=FALSE,multivariate=TRUE)
    gpsc_out[1,11] <- as.character(sprintf('%.2f; %.2f; %.2f',converg_test$psrf[4],converg_test$psrf[5],converg_test$psrf[6]))
    out <- rbind(out, gpsc_out)

}

colnames(out) <- c('GPSC','rootprob','likelihood','prior','mu','sigma','alpha','Root date', 'ESS(mu; sigma; alpha)','DIC(dates;random;verdict)', 'gelman(mu;sigma,alpha)')

write.csv(out,file="bactdate.replicates_burn10_results.csv",row.names=FALSE)