#source("http://bioconductor.org/biocLite.R")
#biocLite("tspair")

#labels in grp  1 x N
#levels in dat  N x G
library(tspair)
library(bigmemory)
data(tspdata)

#prob. function with weights
#p <- function(ilevels, jlevels, weight=rep(1,dim(dat)[2])){
#  pij <-sum(weight*(ilevels<jlevels))/length(ilevels)
#  return(pij)
#}


##weighted pr function, uses p. defaults to unweighted if none supplied
#wpr <- function(dat, N, weight=rep(1,N)){
#  pMat <- matrix(data=0,nrow=N,ncol=N)
#   for(i in 1:N){
#      
#      for(j in i:N){
#        pMat[i,j] <- sum(weight*(dat[i,] < dat[j,]))/length(dat[i,])
#      }
#  }
#  #pMat=(t(pMat)+pMat)
#  return(pMat)  #return symmetric prob matrix i.e. all Pij
#}


#weighted pr function, uses p. defaults to unweighted if none supplied
wpr <- function(dat,weight){
  N<-dim(dat)[1]
  pMat <- matrix(data=0,nrow=N,ncol=N)
   for(i in 1:N){
      for(j in i:N){
        pMat[i,j] <- sum(weight*(dat[i,] < dat[j,]))/length(dat[i,])
      }
  }
  #pMat=(t(pMat)+pMat)
  return(pMat)  #return symmetric prob matrix i.e. all Pij
}

pr <- function(dat,i,j,weight){
     return(sum(weight*(dat[i,] < dat[j,])))
}


thing <- function(dat,i,weight){
      N<-dim(dat)[1]
      p <- matrix(data=0,nrow=1,ncol=N)
      p <- mapply(pr, j=i:N, MoreArgs=list(dat=dat, weight=weight, i=i))
      return(p/N)
}

#weighted pr function, uses p. defaults to unweighted if none supplied, partially vectorized
wpr2 <- function(dat,weight){
  N<-dim(dat)[1]
  pMat <- matrix(data=0,nrow=N,ncol=N)
  pMat <- mapply(thing, i=1:N, MoreArgs=list(dat=dat, weight=weight),SIMPLIFY = TRUE)
  return(pMat)  #return symmetric prob matrix i.e. all Pij
}


#function takes a maxtrix of expression values and a vector of group labels
toppair <- function(dat,label, weight){

   
  
  #subset data based on labels
  C1 <- dat[,label==1]
  C2 <- dat[,label==-1]
  
  weight1 <- dat[,label==1]
  weight2 <- dat[,label==-1]
  
  #probs for each group
  prC1<-wpr(C1,weight1)
  prC2<-wpr(C2,weight2)
  
  #delta (prob difference between classes)
  delta<-abs(prC1-prC2)
  
  deltamax <- which(delta==max(delta),arr.ind=T)

  return(list(deltamax, delta, prC1, prC2))
}

lessthan<-function(X,Y,dat,m){return(as.numeric(dat[X,m]<dat[Y,m]))}

tspandy<- function(dat, grp, weight){
  #1 and 2 as factor for group labels
  grp<-as.factor(as.numeric(factor(grp)))
  
#split into group
  C1<-dat[,which(grp==1)]
  C2<-dat[,which(grp==2)]

#split into group
  W1<-weight[which(grp==1)]
  W2<-weight[which(grp==2)]
  
  #group sizes
  M1<-length(subset(grp, grp==1))
  M2<-length(subset(grp, grp==2))
  #number of observations
  M<-M1+M2
  #number of variables
  N<-dim(dat)[1]
   
  
  pr1<- array(0, dim = c(M1, N, N))
  pr2<- array(0, dim = c(M2, N, N))
  
  for(m in 1:M1){
    pr1[m,,]<-outer(1:N, 1:N, "lessthan", C1, m)  
  }
  
  for(m in 1:M2){
    pr2[m,,]<-outer(1:N, 1:N, "lessthan", C2, m)  
  }
  
  #probs for each class
  prob1<-colSums((W1*pr1),1)/M1
  prob2<-colSums((W2*pr2),1)/M2
  
  Delta<-abs(prob1-prob2)
    
  pairs<-which(Delta == max(Delta), arr.ind = TRUE) 
  
	
  return(list(pairs,Delta, prob1,prob2))
}

tspandyB<- function(dat, grp, weight){
  #1 and 2 as factor for group labels
  grp<-as.factor(as.numeric(factor(grp)))
  
#split into group
  C1<-dat[,grp==1]
  C2<-dat[,grp==2]

#split into group
  W1<-weight[grp==1]
  W2<-weight[grp==2]
  
  #group sizes
  M1<-length(subset(grp, grp==1))
  M2<-length(subset(grp, grp==2))
  #number of observations
  M<-M1+M2
  #number of variables
  N<-dim(dat)[1]
   

	Pr1list<-list()
	for(i in 1:M1){
 		Pr1list[[i]]<-filebacked.big.matrix(N, N, init=0, type="double", backingpath="C:/Users/s4076911/Documents/R/", backingfile=paste("pr1",i,".bin", sep=""), descriptorfile=paste("pr1",i,".desc", sep=""))
	}
	
	Pr2list<-list()
	for(i in 1:M2){
 		Pr1list[[i]]<-filebacked.big.matrix(N, N, init=0, type="double", backingpath="C:/Users/s4076911/Documents/R/", backingfile=paste("pr2",i,".bin", sep=""), descriptorfile=paste("pr2",i,".desc", sep=""))
	}
 
      
  for(m in 1:M1){
    Pr1list[[m]]<-outer(1:N, 1:N, "lessthan", C1, m)  
  }
  
  for(m in 1:M2){
    Pr2list[[m]]<-outer(1:N, 1:N, "lessthan", C2, m)  
  }
  
  #probs for each class
  prob1<-big.matrix(N, N, init=0, type="double")
  prob2<-big.matrix(N, N, init=0, type="double")
  Delta<-big.matrix(N, N, init=0, type="double")
  
  prob1<-
  prob2<-colSums((W2*pr2),1)/M2
  
  Delta<-abs(prob1-prob2)
    
  pairs<-which(Delta == max(Delta), arr.ind = TRUE) 
  
	
  return(list(pairs,Delta, prob1,prob2))
}




#bagging loop
bagger <- function(dat,grp,M=10){
      
  #sets labels as -1, +1
  label=as.numeric(factor(grp))-1
  label[label==0]=-1
  
  #initial weights, all 1/N
  N<-dim(dat)[2]
  errM<-matrix(data=0,nrow=M,ncol=1)
  classesM<-matrix(data=0,nrow=M,ncol=N)
  pairM<-vector("list", M)
  confuseM<-array(data=0, dim=c(M,2,2))
  alphaM<-matrix(data=0,nrow=M,ncol=1)
  
  for(m in 1:M){
    samp<-sample(N,N, replace =T)
    pairM[[m]] <- tspandy(dat[,samp],label[samp],rep(1,N))
    classesM[m,] <- classify(dat, pairM[[m]], pairM[[m]][3], pairM[[m]][4])
    #errM[m] <- error(dat, label, pairM[[m]], pairM[[m]][3], pairM[[m]][4], weight)
    confuseM[m,,]<-table(label,classesM[m,])
    errM[m] <- (confuseM[m,2,1]+confuseM[m,1,2] ) /sum(confuseM[m,,]) 
    
  }
  
  return(list(alphaM, errM, classesM, confuseM, pairM))
}







#boosting loop
boost <- function(dat,grp,M=5){

  if(M==0){
   #put something for non-fixed iterations, stopping conditions
  }

  #sets labels as -1, +1
  label=as.numeric(factor(grp))-1
  label[label==0]=-1
   
  #initial weights, all 1/N
  N<-dim(dat)[2]
  weight=rep(1/N,N)
  weightM<-matrix(data=0,nrow=M,ncol=N)
  errM<-matrix(data=0,nrow=M,ncol=1)
  alphaM<-matrix(data=0,nrow=M,ncol=1)
  classesM<-matrix(data=0,nrow=M,ncol=N)
  pairM<-vector("list", M)
  confuseM<-array(data=0, dim=c(M,2,2))
  
  
  for(m in 1:M){
    weightM[m,]<-weight 
    pairM[[m]] <- tspandy(dat,label,weight)
    classesM[m,] <- classify(dat, pairM[[m]], pairM[[m]][3], pairM[[m]][4])
    #errM[m] <- error(dat, label, pairM[[m]], pairM[[m]][3], pairM[[m]][4], weight)
    confuseM[m,,]<-table(label,classesM[m,])
	errM[m] <- (confuseM[m,2,1]+confuseM[m,1,2] ) /sum(confuseM[m,,]) 
	alphaM[m] <- alphacalc(errM[m])
	weight <- newweight(weight, alphaM[m], classesM[m,], label)
  }

  return(list(alphaM, errM, classesM, confuseM, pairM))
}


#boosting loop
boostB <- function(dat,grp,M=0){

  if(M==0){
   #put something for non-fixed iterations, stopping conditions
  }

  #sets labels as -1, +1
  label=as.numeric(factor(grp))-1
  label[label==0]=-1
   
  #initial weights, all 1/N
  N<-dim(dat)[2]
  weight=rep(1/N,N)
  weightM<-big.matrix(M, N, init=0, type="double")
  errM<-matrix(data=0,nrow=M,ncol=1)
  alphaM<-matrix(data=0,nrow=M,ncol=1)
  classesM<-big.matrix(M, N, init=0,type="double")
  pairM<-vector("list", M)
  confuseM<-array(data=0, dim=c(M,2,2))
  
  
  for(m in 1:M){
    weightM[m,]<-weight 
    pairM[[m]] <- tspandy(dat,label,weight)
    classesM[m,] <- classify(dat, pairM[[m]], pairM[[m]][3], pairM[[m]][4])
    #errM[m] <- error(dat, label, pairM[[m]], pairM[[m]][3], pairM[[m]][4], weight)
    confuseM[m,,]<-table(label,classesM[m,])
	errM[m] <- (confuseM[m,2,1]+confuseM[m,1,2] ) /sum(confuseM[m,,]) 
	alphaM[m] <- alphacalc(errM[m])
	weight <- newweight(weight, alphaM[m], classesM[m,], label)
  }

  return(list(alphaM, errM, classesM, confuseM, pairM))
}


#simple weighted error calculation
error <-function(dat, grp, pair, pr1, pr2, weight){
  E<-0
  E <- sum(weight*(classify(dat,pair,pr1,pr2) != grp))/sum(weight)
  return(E)
}

#simple unweighted error calculation
error2 <-function(dat, grp, pair, pr1, pr2, weight){
  E<-0
  E <- sum((classify(dat,pair,pr1,pr2) != grp))/sum(weight)
  return(E)
}


#basic clasifier, quick anyway, could clean up code up clear what's happening this way.
#this works
classify <-function(dat, pair, pr1, pr2){
  pair<-unlist(pair)
  
	N<-dim(dat)[2]

	C<-rep(0,N)



  pij1 <-  pr1[[1]][pair[1],pair[2]]
	
  pij2 <-  pr2[[1]][pair[1],pair[2]]
	



    for(n in 1:N){

          if( pij1 >= pij2){
              if(dat[pair[1],n] < dat[pair[2],n]){
                  C[n] <- -1
              }
              if(dat[pair[1],n] >= dat[pair[2],n]){
                  C[n] <- 1
              }
          }
          
          if( pij1 < pij2){
              if(dat[pair[1],n] < dat[pair[2],n]){
                  C[n] <- 1
              }
              if(dat[pair[1],n] >= dat[pair[2],n]){
                  C[n] <- -1
              }
          }
          
          
          
    }
        
       
  return(C)
   
}  
classifysingle <-function(dat, pair, pr1, pr2){
  pair<-unlist(pair)
 

	C<-0



  pij1 <-  pr1[[1]][pair[1],pair[2]]
	
  pij2 <-  pr2[[1]][pair[1],pair[2]]




          if( pij1 >= pij2){
              if(dat[pair[1]] < dat[pair[2]]){
                  C <- -1
              }
              if(dat[pair[1]] >= dat[pair[2]]){
                  C <- 1
              }
          }
          
          if( pij1 < pij2){
              if(dat[pair[1]] < dat[pair[2]]){
                  C <- 1
              }
              if(dat[pair[1]] >= dat[pair[2]]){
                  C <- -1
              }
          }
          
          
          
    
        
       
  return(C)
   
} 



           

#calculate alpha
alphacalc <- function(error){
    alpha <- log((1-error)/error)
    return(alpha)
} 

#recalculate weights
newweight <-function(weight, alpha, classes, grp){
  w<- weight*exp(alpha*(classes!=grp))    
  return(w)

}

#bigclassify(test[[1]],test[[3]])
#full boosted classifier
bigclassify <- function(all){
  bclass<- sign(colSums(as.vector(all[[1]])*as.matrix(all[[3]]),1))
  return(bclass)
}


bigclassify2 <- function(all, dat){
print(all[[1]])
classes<-array(data=0,dim=length(all[[1]]))
for(i in length(all[[1]])){
	classes[i]<-classifysingle(dat, all[[5]], all[[5]][[i]][3], all[[5]][[i]][4])
}	
bclass<- sign(as.vector(all[[1]])%*%classes)
return(bclass)
}





# ptm<- proc.time()
# test1<-wpr(dat)
# proc.time() - ptm
# 
# ptm<- proc.time()
# test2<-wpr2(dat)
# proc.time() - ptm



# ptm<- proc.time()
# test<-boost(dat, grp, 100)
# proc.time() - ptm


#test<-boost(dat, grp, 5)

library(boot)
basicCI<- function(data, indices){
basic<-tspandy(data[,indices], label[indices], rep(1/(length(label[indices])),length(label[indices])))
basicclassify <- classify(data[,indices], basic, basic[3], basic[4])
basictab<-table(label,basicclassify)
err0<-1-(basictab[1]/(basictab[1]+basictab[3]))
err1<-(basictab[2]/(basictab[2]+basictab[4]))
err2<-(basictab[2]+basictab[3])/sum(basictab)
return(c(err0,err1,err2))
}

#basic1<-boot(data=t(dat), statistic=basicCI, R=100)

fullCI<- function(data, indices){
full<-boost(data[,indices], grp[indices],25)
fullclass<-bigclassify(full)
fulltab<-table(label,fullclass)
err0<-1-(fulltab[1]/(fulltab[1]+fulltab[3]))
err1<-(fulltab[2]/(fulltab[2]+fulltab[4]))
err2<-(fulltab[2]+fulltab[3])/sum(fulltab)
return(c(err0,err1,err2))
}
 

#


full1<-boot(data=t(dat2), statistic=fullCI, R=10)
full2<-boot(data=t(dat2), statistic=basicCI, R=10)

#LOOCV
err<-0
label=as.numeric(factor(grp2))-1
label[label==0]=-1
for(i in 1:length(grp2)){
	full<-boost(dat[1:5000,-i], grp[-i], 10)
	iclass<-bigclassify(full)
	err<-err+as.numeric(iclass==label[i])
}
err





#save(err, file = "C:\\Users\\s4076911\\Dropbox\\Work\\TSP\\test8.RData")

#save.image("C:\\Users\\s4076911\\Dropbox\\Work\\TSP\\test8.RData")

#LOOCV for base is 0.36






# BIG<-read.big.matrix("C:\\Users\\s4076911\\Dropbox\\Work\\TSP\\cancer1.csv", type="single",header=T,backingfile="cancer.bin", descriptorfile ="cancer.desc")
# BIG<- BIG[,-1]
# 
# 
# grp1<-as.vector(data1[dim(data1)[1],2:21])
# full<-boost(data2[1:100,], grp1, 1)
# 
# 
# normalize <- function(x){
# 	return((x-mean(x))/sd(x))
# }
# 
# BIG2<- apply(BIG, 2, normalize)

thing1<-boost(dat2[1:5000,],grp2)



