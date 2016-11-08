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
#existing stuff
data(tspdata)

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
getTSPClass<-function(dat, grp, k=1, weights = rep(1,n)){
  
  #convert to numbers, keep labels for later
  numGrp <- as.numeric(factor(grp))
  numGrp2 <- 2*(as.numeric(factor(grp))-1.5)
  namesGrp <- levels(factor(grp))
  
  n<-length(grp)
  #add something here to sort columns / labels into blocks
  
  #send into C  
  pairs <- getData(numGrp, dat, k, weights)
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
###################################################

# tie breaker differs from Tan et al 2005. Unsure why.
dat_test <- dat[1:50, c(1:5, 45:50)]
tsp1<-getTSPClass(dat_test,grp[c(1:5, 45:50)],5)
tsp2<-getTSPClass(dat,grp,3)

#quickTSP(numGrp, dat, rep(1,n))

dat_test2 <- dat[1:50, c(20:25, 40:45)]

classify(dat_test2,tsp1)
classify(dat,tsp2)

#TO DO
#--
# - boost
# - bag
# - tree
# - method chooser
# - Bootstrap
# - CV loop
# - profile all code for speed







