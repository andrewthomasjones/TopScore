#source("https://bioconductor.org/biocLite.R")
#biocLite("tspair")

library(tcltk)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(rbenchmark)
library(ggplot2)
library(tspair)
library(ktspair)
library(rbenchmark)
library(boot)

#existing stuff
data(tspdata)
#testing
dat_test <- dat[1:50, c(1:5, 45:50)]
grp_test <- grp[c(1:5, 45:50)]
dat_test2 <- dat[1:50, c(20:25, 40:45)]

tsp1 <- tspcalc(dat,grp)
# tspplot(tsp1)
# 
# #sig test
# out <- tspsig(dat,grp,B=50,seed=12355)
# 
# #predict
# predict(tsp1,eSet2)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
Rcpp::sourceCpp('/windowsDrive/Users/Andy/Dropbox/Work/TSP/TSP/Tsp_C++.cpp')

###################################################
# tie breaker differs from Tan et al 2005. Unsure why.
ktsp<-function(dat, grp, k=1, weights = rep(1,n)){
    
  #convert to numbers, keep labels for later
  numGrp <- as.numeric(factor(grp))
  numGrp2 <- 2*(as.numeric(factor(grp))-1.5)
  namesGrp <- levels(factor(grp))
    
  n<-length(grp)
  #add something here to sort columns / labels into blocks
  

  
  #send into C  
  pairs <- getData(dat, numGrp,  k, weights)
  # r indexing
  pairs [,1:2] <-pairs [,1:2]+1 
  colnames(pairs)= c("1st","2nd","delta", "gamma")
  return(list(sample_size=n, n_pairs=k, names_group=namesGrp, labels_group=numGrp2, pairs = pairs))
  
}


###################################################
#basic clasifier
classify <-function(dat, tsp_classif){
  
  N<-dim(dat)[2]
  
  vote_dir <- 2*((tsp_classif$pairs[,3]>0) - 0.5) 
  votes <- array(0,c(tsp_classif$n_pairs,N))
  
  for(k in 1:tsp_classif$n_pairs){
    
    votes[k,]<-2*((dat[tsp_classif$pairs[k,1],] < dat[tsp_classif$pairs[k,2],])- 0.5) 
    
    
  }
  
  directed_votes <- (vote_dir)%*%votes
  C = 2*( (directed_votes > 0) - 0.5) 
  
  
  return(C)
  
}
classifyAda <-function(dat, tsp_classif){
  
  N<-dim(dat)[2]
  
  vote_dir <- 2*((tsp_classif$pairs[,3]>0) - 0.5)*tsp_classif$alpha 
  
  votes <- array(0,c(tsp_classif$n_pairs,N))
  
  for(k in 1:tsp_classif$n_pairs){
    
    votes[k,]<-2*((dat[tsp_classif$pairs[k,1],] < dat[tsp_classif$pairs[k,2],])- 0.5) 
    
    
  }
  
  directed_votes <- (vote_dir)%*%votes
  C = 2*( (directed_votes > 0) - 0.5) 
  
  
  return(C)
  
} 

###################################################
#bagging loop
bagger <- function(dat,grp,M){
  #convert to numbers, keep labels for later
  numGrp <- as.numeric(factor(grp))
  numGrp2 <- 2*(as.numeric(factor(grp))-1.5)
  namesGrp <- levels(factor(grp))
  
  #initialize
  n<-length(grp)
  pairs<-array(0, c(M,4))
  
  #bootstrap resample
  for(m in 1:M){
    samp <-sample(n)
    pairs[m,] <- getData(dat[,samp], numGrp[samp],  1, rep(1,n))
    
  }
  pairs[,1:2]<-pairs[,1:2]+1
  
  return(list(sample_size=n, n_pairs=M, names_group=namesGrp, labels_group=numGrp2, pairs = pairs))
}

##################################################
#calculate alpha
alphacalc <- function(error){
  alpha <- log((1-error)/error)
  return(alpha)
} 

#recalculate weights
newweight <-function(weight, alpha, classes, grp){
  w<- weight*exp(alpha*(classes!=grp))
  if(all(is.na(w))){
    n<-length(grp)
    w<-rep(1/n,n)
  }
  return(w)
  
}



#boosting loop
booster <- function(dat,grp,M=5){
  
  if(M==0){
    #put something for non-fixed iterations, stopping conditions
  }
  
  #convert to numbers, keep labels for later
  numGrp <- as.numeric(factor(grp))
  numGrp2 <- 2*(as.numeric(factor(grp))-1.5)
  namesGrp <- levels(factor(grp))
  
  #initialize
  n<-length(grp)
  pairs<-array(0, c(M,4))
  
  weight=rep(1/n,n)
  weightM<-matrix(data=0,nrow=M,ncol=n)
  errM<-matrix(data=0,nrow=M,ncol=1)
  alphaM<-matrix(data=0,nrow=M,ncol=1)
  classesM<-matrix(data=0,nrow=M,ncol=n)
  pairM<-list()
  errs<-list()
  confuseM<-array(data=0, dim=c(M,2,2))
  
  
  for(m in 1:M){
    weightM[m,] <- weight 
    pairM[[m]] <- ktsp( dat,numGrp, 1, weight)
    classesM[m,] <- classify(dat, pairM[[m]])
    confuseM[m,,]<-table(numGrp2 ,classesM[m,])
    errM[m] <- (confuseM[m,2,1]+confuseM[m,1,2] ) /sum(confuseM[m,,]) 
    alphaM[m] <- alphacalc(errM[m])
    weight <- newweight(weight, alphaM[m], classesM[m,], numGrp2)
    pairs[m,] <- pairM[[m]]$pairs[1,]
    
    if(m>1){
      temp_clasif = list(sample_size=n, n_pairs=m, names_group=namesGrp, labels_group=numGrp2, pairs = pairs[1:m,], alpha = unlist(alphaM[1:m]))
      tab <- table(numGrp2,classifyAda(dat,temp_clasif))
      errs[m]<-(tab[2,1] + tab[1,2])/sum(tab)
    }else{
      errs[1] =0
    }
  }
  
  return(list(sample_size=n, n_pairs=M, names_group=namesGrp, labels_group=numGrp2, pairs = pairs, err =unlist(errs[1:M]), alpha = unlist(alphaM[1:M])))
}



#boosting loop
kbooster <- function(dat,grp,M=5, k=3){
  
  if(M==0){
    #put something for non-fixed iterations, stopping conditions
  }
  
  #convert to numbers, keep labels for later
  numGrp <- as.numeric(factor(grp))
  numGrp2 <- 2*(as.numeric(factor(grp))-1.5)
  namesGrp <- levels(factor(grp))
  
  #initialize
  n<-length(grp)
  pairs<-array(0, c(M*k,4))
  
  weight=rep(1/n,n)
  weightM<-matrix(data=0,nrow=M,ncol=n)
  errM<-matrix(data=0,nrow=M,ncol=1)
  alphaM<-matrix(data=0,nrow=M*k,ncol=1)
  classesM<-matrix(data=0,nrow=M,ncol=n)
  pairM<-list()
  errs<-list()
  confuseM<-array(data=0, dim=c(M,2,2))
  
  i=1
  for(m in 1:M){
    weightM[m,] <- weight 
    pairM[[m]] <- ktsp( dat,numGrp, k, weight)
    classesM[m,] <- classify(dat, pairM[[m]])
    confuseM[m,,]<-table(numGrp2 ,classesM[m,])
    errM[m] <- (confuseM[m,2,1]+confuseM[m,1,2] ) /sum(confuseM[m,,]) 
    alphaM[i:(i+k-1),1] <- rep(alphacalc(errM[m]),k)
    weight <- newweight(weight, alphaM[m], classesM[m,], numGrp2)
    #print(pairM[[m]]$pairs[1:k,])
    #print(pairs[i:(i+k-1),])
    pairs[i:(i+k-1),] <- pairM[[m]]$pairs[1:k,]
    i<-i+k
    if(m>1){
      temp_clasif = list(sample_size=n, n_pairs=m, names_group=namesGrp, labels_group=numGrp2, pairs = pairs[1:m,], alpha = unlist(alphaM[1:m]))
      tab <- table(numGrp2,classifyAda(dat,temp_clasif))
      errs[m]<-(tab[2,1] + tab[1,2])/sum(tab)
    }else{
      errs[1] =0
    }
  }
  
  return(list(sample_size=n, n_pairs=M*k, names_group=namesGrp, labels_group=numGrp2, pairs = pairs, err =unlist(errs[1:M]), alpha = unlist(alphaM[1:(M*k)])))
}
















#################################################################
boosterCI<- function(data, indices){
  M=5
  basic<-ktsp(data[,indices], grp[indices], M)
  basicclassify <- classify(data[,indices], basic)
  basictab<-table(basic$labels_group, basicclassify)
  
  #errors, which is which
  err0<- (basictab[3]/(basictab[1]+basictab[3]))
  err1<-(basictab[2]/(basictab[2]+basictab[4]))
  err2<-(basictab[2]+basictab[3])/sum(basictab)
  return(c(err0,err1,err2))
}
##########################################################################

#testing
tsp1<-ktsp(dat,grp ,7)
tsp1
tsp2<-bagger(dat, grp , 10)
tsp3<-booster(dat, grp, 20)

tspk<-kbooster(dat, grp, 2)

#tsp3<-bagger(dat, grp ,100)
cl1<-classify(dat,tsp1)
cl2<-classify(dat,tsp2)
cl3<-classifyAda(dat,tspk)

table(tsp1$labels_group ,cl1)
table(tsp2$labels_group ,cl2)
table(tsp3$labels_group ,cl3)

boot(data=t(dat), statistic=boosterCI, R=100)

##################################################
#LOOCV
err<-0
numGrp2 <- 2*(as.numeric(factor(grp))-1.5)
for(i in 1:length(grp)){
  print(paste(i, " of ",length(grp) ))
  full<-kbooster(dat[,-i], grp[-i], 100,3)
  iclass<-classify(dat[,c(i,1)], full)[1]
  err<-err+as.numeric(iclass==numGrp2[i])
}
err/length(grp)
#################################################



#TO DO
# -  put boost/bag wholely in C
# -  classifier in C to speed up boost
# -   disjoint pairs all methods
# -  k-fold CV
# -  CV into C
# - warning message zero error case for boost
# - method chooser / clean interface
# - profile all code for speed









