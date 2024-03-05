# 3. Step: Combing the aDDM simulations from Simulation_ADDM_Vero file, according to Krajbich et al. (2011)
# the goal is to figure out whether fitting the aDDM on simulated data can give me a good idea of the effect and sample size I will need, to investigate the overall value effect
#rm(list = ls())

library(plyr);
library(dplyr);
library(gtools)
setwd("C:/Users/Asus/OneDrive - University of Aberdeen/Desktop/Aberdeen_Uni/cap/HELP_ME_UNDERSTAND_MODELLING/Krajbich_Simulations_aDDM_2")

experimentN = 8                                                                              #insert experiment number, make sure to have a folder with 'EXP' and the number in your path, insert all important param/hyperparam files
#Combine condition-level files
experimentN= 8

df500_f <- list.files(path = paste0(getwd(),"/EXP",experimentN),pattern = paste0("500ms"))    #is a function in R used to list files and directories in a specified directory
df500_f <- mixedsort(df500_f)
#Read in each file, but first add trial numbers
df500_f2 <- lapply(df500_f, function(x){
  f <- read.csv(paste0(getwd(),"/EXP",experimentN,"/",x), header = TRUE)
  # Data/",x), header = TRUE)
  f$trial <- f$newtrial <- 0                                                                  # adds 2 new columns to the data frame f, trial and newtrial, initializes them with 0
  f2_list<-lapply(unique(f$subject),function(s){                                              # to each sub. in the frame, function x is applied,
    tmpdf <- f[f$subject==s,]
    tmpdf$newtrial<- ifelse(tmpdf$rt==lag(tmpdf$rt) & 
                              tmpdf$choice==lag(tmpdf$choice) &
                              tmpdf$nfix==lag(tmpdf$nfix),0,1)
    tmpdf$newtrial[1] <- 1
    count=0
    for(i in 1:nrow(tmpdf)){
      if(tmpdf$newtrial[i]==1) count <- count+1
      tmpdf$trial[i] <- count
    }
    return(tmpdf)
  })
  f2<- do.call(rbind.data.frame,f2_list)
  return(f2)
})
df500_f3 <- do.call(rbind.data.frame, df500_f2)

#Exclude extra fixations to export to HDDM
df500_f3$VD <- round(df500_f3$leftval - df500_f3$rightval,2)
df500_f3$absVD <- abs(df500_f3$VD)
df500_f3$OV <- round(df500_f3$leftval + df500_f3$rightval,2)
df500<-ddply(df500_f3, .(subject, absVD, OV, trial, choice, rt, leftval, rightval), summarise,
             VD = unique(VD),
             fixDurLeft = sum(fixdur[roi==1]),
             fixDurRight = sum(fixdur[roi==2]),
             fixDurTotal = sum(fixdur))
df500$EL <- ((df500$fixDurLeft)/(df500$fixDurTotal))
df500$ER <- ((df500$fixDurRight)/(df500$fixDurTotal))
df500$x1 <- (df500$leftval*df500$EL)-(df500$rightval*df500$ER)
df500$x2 <- (df500$leftval*df500$ER)-(df500$rightval*df500$EL)
write.csv(df500,paste0("addmSims_MathPsych_EXP",experimentN,"_HDDM_500ms.csv"), row.names = F)


