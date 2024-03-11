# 2. Step: Simulate aDMM to evaluate the impact of magnitude on parameters recovered via DDM
#Packages#
#rm(list = ls())
install.packages("plyr")
install.packages("foreach")
install.packages("psych")

library(plyr)
library(foreach)
library(psych)
set.seed(51091)
#Change to path
setwd("C:/Users/Asus/OneDrive - University of Aberdeen/Desktop/Aberdeen_Uni/cap/HELP_ME_UNDERSTAND_MODELLING/Krajbich_Simulations_aDDM_2")
fixOptions<-500



# Low OV and VD conditions
lowValVD2 <- c(1, 3)      # OV = 4    # VD = 2 low
midValVD4 <- c(2, 6)      # OV = 8    # VD = 4 high
modValVD4 <- c(4, 8)      # OV = 12    # VD = 4 high
highValVD2 <- c(7, 9)     # OV = 16    # VD = 2 low

# Constructing valMat matrix
valMat <- matrix(0, nrow = 4, ncol = 4)
valMat[, 1] <- c(1, 2, 3, 4)
valMat[, 2] <- c(2, 4, 4, 2)
valMat[, 3] <- c(lowValVD2[1], midValVD4[1], modValVD4[1], highValVD2[1])
valMat[, 4] <- c(lowValVD2[2], midValVD4[2], modValVD4[2], highValVD2[2])



experimentN=8                     #8 hyperparameter setups (CSV Files), according to Smith & Krajbich (2019)

#Generate hyper-params based on ranges from Smith & Krajbich (2019) Psych Science
#Note that EXP 2 uses the same hyper-params as EXP1, but with greater variability
if (experimentN %in% c(1,2)){
  hypparams<-read.csv("Exp1_hyperparams.csv")
} else {
  dmin<-0.0002
  dmax<-0.0004
  sigmin<-0.02
  sigmax<-0.03 
  thetamin<-.1
  thetamax<-.9
  ndtmin<-350
  ndtmax<-650
  dsigmin<-0
  dsigmax<-0
  hypparams<-matrix(0,nrow=1,ncol=5)
  hypparams[,1]<-runif(1,dmin,dmax)
  hypparams[,2]<-runif(1,sigmin,sigmax)
  hypparams[,3]<-runif(1,thetamin,thetamax)
  hypparams[,4]<-round(runif(1,ndtmin,ndtmax))
  hypparams[,5]<-runif(1,dsigmin,dsigmax)
  hypparamsDF<-data.frame(hypparams)
  colnames(hypparamsDF)<-c("d","sigma","theta","ndt","dsig")
  write.csv(hypparamsDF,
            file = paste0("Exp", experimentN, 
                          "_hyperparams.csv"),row.names = F)
}


SubjNs = c(1:40)
nparams = length(SubjNs)
ntrial = 100
#Experiment 1 used a narrow range of values
if (experimentN==1){
  subjParams<-matrix(0,nrow=nparams,ncol=6)
  subjParams[,1]<-SubjNs
  subjParams[,2]<-rnorm(nparams,hypparams[,1], .00001) #d
  subjParams[,3]<-rnorm(nparams,hypparams[,2], .001) #sigma
  subjParams[,4]<-rnorm(nparams,hypparams[,3], .01) #theta
  subjParams[,5]<-round(rnorm(nparams,hypparams[,4], 10)) #ndt
} else if(experimentN %in% c(2,7,8)){
  subjParams<-matrix(0,nrow=nparams,ncol=6)
  subjParams[,1]<-SubjNs
  subjParams[,2]<-rnorm(nparams,hypparams[,1],.0001) #.00001) #d
  subjParams[,3]<-rnorm(nparams,hypparams[,2], .01) #.001) #sigma
  subjParams[,4]<-rnorm(nparams,hypparams[,3], .1) #.01) #theta
  subjParams[,5]<-round(rnorm(nparams,hypparams[,4], 100)) #10)) #ndt
} else {
  #Exp 3 tests if d is more variable
  #Exp 4 tests if only sigma is more variable
  #Exp 5 tests if theta is more variable
  #Exp 6 tests if ndt is more variable
  subjParams<-matrix(0,nrow=nparams,ncol=6)
  subjParams[,1]<-SubjNs
  subjParams[,2]<-abs(rnorm(nparams,hypparams[,1], 
                            ifelse(experimentN==3,.0001,.00001))) #.00001) #d
  subjParams[,3]<-rnorm(nparams,hypparams[,2], 
                        ifelse(experimentN==4,.01, .001)) #.001) #sigma
  subjParams[,4]<-rnorm(nparams,hypparams[,3], 
                        ifelse(experimentN==5,.1,.01)) #.01) #theta
  subjParams[,5]<-round(rnorm(nparams,hypparams[,4],
                              ifelse(experimentN==6, 100, 10))) #10)) #ndt
}
#subjParams <- read.csv("Exp1_Subj_Params.csv")

colnames(subjParams) <- c("SubjN","d","sigma","theta","ndt","dsigma")
write.csv(subjParams,paste0("Exp",experimentN,"_Subj_Params.csv"),row.names = F)


experimentN <- 8
subjParams <- read.csv(paste0("Exp", experimentN, "_Subj_Params.csv"))

sapply(fixOptions, function(x) {
  thisFix <- x
  valDifferences <- unique(valMat[,2])
  
  sapply(valDifferences, function(y){
    vDiff <- y
    valOV <- unique(valMat[,1])
    
    sapply(valOV, function(z){
      # Find the correct pairing based on the values in columns 1 and 2
      if (z == 1){
        vDiff <- 2
      } else if (z == 2){
        vDiff <- 4
      } else if (z == 3) {
        vDiff <- 4
      } else if (z == 4) {
        vDiff <- 2
      }
      
      theseVals <- valMat[valMat[,1] == z & valMat[,2] == vDiff, 3:4]
      nVals <- length(theseVals)
      value_min <- min(theseVals)
      value_max <- max(theseVals)
    
    #Function and Data#
    source("sim_addm_fun_ctrlValueDiff_v18.R") #was v18
    
    value_diff <- vDiff
    value_range <- theseVals
    fix_length <- thisFix
    
    params <- matrix(0, nrow = nparams, ncol = 11)
    params[, 1] <- subjParams[, 1]
    params[, 2] <- ntrial
    params[, 3] <- subjParams[, 2] #d
    params[, 4] <- subjParams[, 3] #sigma
    params[, 5] <- subjParams[, 4] #theta
    params[, 6] <- subjParams[, 5] #ndt
    params[, 7] <- 0               #dsigma
    params[, 8] <- value_diff
    params[, 9] <- fix_length
    params[, 10] <- value_min
    params[, 11] <- value_max
    
    simData <- matrix(0, nrow = params, ncol = 11)
    iterations <- nrow(params)
    
    simData <- foreach(i = 1:iterations, .combine = rbind) %do% {
      sim_addm_fun_ctrlValueDiff_v18(params[i, 1], params[i, 2], params[i, 3], params[i, 4],
                                     params[i, 5], params[i, 6], params[i, 7], params[i, 8],
                                     params[i, 9], params[i, 10], params[i, 11])
    }
    
    valClass <- ifelse(mean(value_range) > 7, "High", 
                       ifelse(mean(value_range) < 3, "Low", 
                              ifelse(mean(value_range) > 3 & mean(value_range) < 5, "LowMid", 
                                     ifelse(mean(value_range) > 5 & mean(value_range) < 7 & length(value_range) < 4, "HighMid", "Mix"))))
    
    if (vDiff < 1 | vDiff > 0) 
      vDiffC <- paste0("exp(-", vDiff * 10, ")")
    if (vDiff == 0 | vDiff >= 1) 
      vDiffC <- vDiff
    
    write.csv(simData, 
              paste0("EXP", experimentN, "_simADDM_", x, "ms_", nparams, "Subj_100Trials_Dif", 
                     vDiffC, "_", valClass, "Conds_v18.csv"))
    })
  })
})
