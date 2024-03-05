#this code expects a data file that has already been run through the "data_prep_for_addm" script and saved as "dat"; nsim is number of simulations, d is drift parameter, sigma is within trial variability in drift, theta is attention parameter, ndt is non-decision time in ms, dsigma is across trial variability in drift. It produces a simulated dataset called "sims".
#Remove left bias (Force half to be Left, the other Right)
#Uniform fixations (25 ms, 100 ms, or 750ms)
#Fixation distribution is uniform (constrained)
#L/R values vary by vlue_diff

sim_addm_fun_ctrlValueDiff_v18<-function(subjn,ntrial,d,sigma,theta,ndt,dsigma,value_diff,fix_length,value_min,value_max,nsim=1){
  # random number generator library
  library(RcppZiggurat)
  library(plyr)
  
  value_range<-range(value_min,value_max)
  
  #ntrial<-nrow(dat) #number of unique trials in the dataset (get this from the data)
  #maxfix<-(ncol(dat)-4)/2 #maximum number of fixations in a given trial
  maxfix<-50
  #preallocate the simulation vectors
  subject <- subjn #subject number
  simSize<-nsim*ntrial #nrow of simulated data
  rt<-rep(0,simSize) #reaction time vector
  choice<-rep(-1,simSize) #choice vector
  fixnum<-rep(-1,simSize) #fixation number vector
  revfixnum<-rep(-1,simSize) #reverse fixation number vector
  fixdur<-rep(-1,simSize) #fixation duration vector
  nfix<-rep(-1,simSize) #total number of fixations vector
  roi<-rep(-1,simSize) #roi vector
  leftval<-rep(-1,simSize) #left value vector
  rightval<-rep(-1,simSize) #right value vector
  
  #New code shuffles all the trial level information and allocates them into
  #a new matrix whose size matches the size of the to-be-simulated data
  
  #Create functions that bootstrap values and puts them in new matrix
  
  #fixations
  #temp_fixlength_mat<-dat[,(5+maxfix):ncol(dat)] #fixation length data
  
  fixationsArray<-round(runif(simSize*maxfix,(fix_length-(fix_length/2)),(fix_length+(fix_length/2))))
  #fixationsArray<-fix_length
  #fixlength_mat<-matrix(0,ncol=maxfix,nrow=ntrial) #
  #each fix length is drawn from an uniform distributoon 0 - 2x (x is mean fixlength)
  # nFixations<-c(1:(maxfix*simSize)) #Number of fixations
  # fixationsArray<-matrix(0,nrow=1,ncol=(maxfix*simSize))
  
  #convert array to matrix
  fixlength_mat<-matrix(fixationsArray,ncol=maxfix,nrow=simSize)  
  
  
  #temp_roi_mat<-dat[,5:(4+maxfix)]
  roi_mat<-matrix(0,nrow = simSize,ncol=maxfix) #ROI data
  LeftROI<-rep_len(c(1,2),maxfix) #ROI starting Left
  RightROI<-rep_len(c(2,1),maxfix) #ROI starting Right
  is.whole <- function(x) is.numeric(x) && floor(x)==x 
  #LeftRight<-rep_len(0,simSize)
  
  #Force ROI to left (even trials) and right (odd trials)
  for (i in 1:(simSize)){
    #LeftRight[i]<-rbinom(n = 1, size = 1, p=.5)
    #if (LeftRight[i]==1) roi_mat[i,]<- LeftROI
    #if (LeftRight[i]==0) roi_mat[i,]<- RightROI
    if (is.whole(i/2)==T) roi_mat[i,]<-LeftROI
    if (is.whole(i/2)==F) roi_mat[i,]<-RightROI
  }
  #temp_value_mat<-(temp_roi_mat==1)*dat[,3]+(temp_roi_mat==2)*dat[,4] #value data
  # lvalue_mat<-dat[,3]*matrix(1,nrow=(dat),ncol=maxfix) #lvalue data
  #rvalue_mat<-dat[,4]*matrix(1,nrow=(dat),ncol=maxfix) #rvalue data
  #row0 = c(1:simSize) #number of trials to be simulated
  #rowOrdr<-sapply(row0, function(x){
  #  bootStrapRow<-sample(c(1:nrow(dat)),1000,replace=T)
  #  thisRow<-sample(bootStrapRow,1)
  # return(thisRow)})
  #value_mat<-temp_value_mat[rowOrdr,]
  # lvalue_mat<-temp_lvalue_mat[rowOrdr,]
  # rvalue_mat<-temp_rvalue_mat[rowOrdr,]
  
  #if more than one value, randomly sample a value from the value range 
  #But, if max value, second value is less; if min value, second value is more
  if (length(value_range)==1) {
    val1 = value_range
    val2 = value_range}
  if (length(value_range)>1){
    val1<-sample(value_range,size=(simSize),replace = T) #lvalue data
    val2<-sapply(c(1:length(val1)),function(x){
      v1<-val1[x]
      rand<-sample(c(0,1),1)
      if (v1==min(value_range)){v2<-v1+value_diff} 
      if (v1==max(value_range)){v2<-v1-value_diff}
      if(v1>min(value_range) & v1<max(value_range)) {v2<-ifelse(rand==1,v1+value_diff,v1-value_diff)}
      return(v2)
    })}
  lvalue_mat<-matrix(val1,nrow=(simSize),ncol=maxfix)
  rvalue_mat<-matrix(val2,nrow=(simSize),ncol=maxfix)
  
  drift_mat<-d*((roi_mat==1)*(lvalue_mat-theta*rvalue_mat)+(roi_mat==2)*(theta*lvalue_mat-rvalue_mat))
  
  counter<-1
  #loop over unique trials
  for (k in 1:ntrial){
    nroi<-roi_mat[k,]
    nfixlength<-fixlength_mat[k,]
    ndrift<-drift_mat[k,]
    ntemp = cumsum(nfixlength)  
    index = rep(0,max(ntemp))
    index[ntemp[1:length(ntemp)-1]+1] = 1
    index[1] = 1
    index = cumsum(index)
    bigdrift = ndrift[index]
    
    #presample vector of noise to sample from
    maxt<-max(ntemp) #maximum possible RT (in ms)
    noises<-zrnorm(maxt*nsim)*sigma #pre-assemble vector of noise to sample from 
    
    #presample across-trial noise
    dnoises<-zrnorm(nsim)*dsigma
    
    j<-0
    
    #loop over number of simulations per trial
    for (i in 1:nsim){
      trt<-Inf
      j<-j+1
      noise<-noises[((j-1)*maxt+1):((j)*maxt)] #sample noise realizations
      evidence<-bigdrift+noise+dnoises[i] #add drift to noise
      #evidence<-noise
      RDV<-cumsum(evidence)#sum up the evidence
      absRDV<-abs(RDV)
      trt<-min(which(absRDV>=1))#find the first time that the evidence exceeds a magnitude of 1
      ifelse(trt==Inf,trt<-maxt,trt<-trt)#set maxt as the finishing time if time runs out
      
      #save the trial data
      fix.num<-min(which(cumsum(nfixlength)>=trt))
      nfix[counter:(counter+fix.num-1)]<-fix.num
      choice[counter:(counter+fix.num-1)]<-ceiling(RDV[trt]/1000000) #was the RDV + or -?
      rt[counter:(counter+fix.num-1)]<-trt+ndt #add non-decision time to the RT
      fixnum[counter:(counter+fix.num-1)]<-(1:fix.num) 
      revfixnum[counter:(counter+fix.num-1)]<-(fix.num:1)
      roi[counter:(counter+fix.num-1)]<-nroi[1:fix.num]
      leftval[counter:(counter+fix.num-1)]<-lvalue_mat[k,1]
      rightval[counter:(counter+fix.num-1)]<-rvalue_mat[k,1]
      #fix the last fixation time
      tfix.length = nfixlength[1:fix.num]
      erro = sum(tfix.length)-trt
      tfix.length[fix.num] = tfix.length[fix.num]-erro
      fixdur[counter:(counter+fix.num-1)]<-tfix.length
      
      counter<-counter+fix.num
    }
  }
  sims<-data.frame(subject,choice,rt,fixnum,nfix,revfixnum,roi,leftval,rightval,fixdur)	
  sims<-sims[sims$choice>=0,]	
  sims<-sims[sims$rt<(maxt+ndt),] #remove trials that did not finish
  rm(noises,absRDV,bigdrift,choice,counter,drift_mat,erro,evidence,fix.num,fixdur,fixlength_mat,fixnum,i,index,j,k,leftval,lvalue_mat,maxfix,maxt,ndrift,nfix,nfixlength,noise,nroi,nsim,ntemp,ntrial,RDV,revfixnum,rightval,roi,roi_mat,rt,rvalue_mat,tfix.length,trt)
  
  return(sims)
}

