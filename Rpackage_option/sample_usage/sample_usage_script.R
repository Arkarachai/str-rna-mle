#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
library(combinat)
library(snow)
library(str.rna.mle)
 

###############
## read data ##
###############

filename=args[1]
oridataset=read.table(filename,colClasses = "character")    
columnname=(c('locus','DNA','motif','class','functional','seq1','seq2'))
colnames(oridataset)=columnname
transformdataset=list()
for (line in 1:nrow(oridataset)){
    individualdataset=list(as.numeric(oridataset$DNA[line]),oridataset$motif[line]
    ,stringToNumList(oridataset$seq1[line]),stringToNumList(oridataset$seq2[line]))
    
    transformdataset[[length(transformdataset)+1]] <- individualdataset
}

############################
## optimizing using optim ##
############################

    tempep1=10**runif(1,-9,-0.3)
    tempQ1=runif(1,0,1)
    tempep2=10**runif(1,-9,-0.3)
    tempQ2=runif(1,0,1)
    initialparameter=c(tempep1,tempQ1,tempep2,tempQ2)
    print(c("initialparameter",initialparameter))
    optimresult=optim(initialparameter,likelihood_lumpMLE,NULL,method="L-BFGS-B",
    lower=c(1e-9,0,1e-9,0),
    upper=c(0.5,1,0.5,1),
    control=list(ndeps=c(0.00001,0.01,0.00001,0.01),
    fnscale=-1),
    tempdataset=transformdataset,bin_size=2,N_core=5)

print('done')
