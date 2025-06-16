library(hdi)

###############################define some functions
localDiagnostic<-function (y, h) {    #see original SaRa algorithm 
  yy = c(rep(0, h - 1), y, rep(0, h))
  n = length(y)
  z = rep(0, n)
  for(i in 1:n){
    z[i]=sum(yy[i:(h+i-1)])/h-sum(yy[(h+i):(2*h-1+i)])/h
  }
  return(z)
}

localMax<-function (y, span = 5) {   #see original SaRa algorithm
  if (length(y) < span * 2 + 1) 
    return(NULL)
  n = length(y)
  index = NULL
  for (i in (span + 1):(n - span)) {
    if (y[i] == max(y[(i - span):(i + span)])) 
      index = c(index, i)
  }
  return(index)
}

refine<-function(y, candidate, itermax=10, delta=3){          #Step2
  #y: data; candidate: change-points candidates in Step 1
  k1=length(candidate)
  psi=c(0,candidate,length(y)+1)
  for(iter in 1:itermax){
    r=sample(2:(k1+1),size=k1,replace=FALSE)
    for(j in 1:k1){
      a=r[j]
      b=min(psi[a+1]-psi[a]-2, psi[a]-psi[a-1]-2, delta)
      candid2=(psi[a]-b):(psi[a]+b)
      s=1:(2*b+1)
      for(l in 1:(2*b+1)){
        leftmean=mean(y[(psi[a-1]+1):candid2[l]])
        rightmean=mean(y[(candid2[l]+1):(psi[a+1]-1)])
        s[l]=(leftmean-rightmean)^2*length(y[(psi[a-1]+1):candid2[l]])*length(y[(candid2[l]+1):(psi[a+1]-1)])
      }
      psi[a]=candid2[which.max(s)]
    }
  }
  return(psi[2:(k1+1)])
}

global<-function(y, candidate){          #exhaustive search in Step2
  #y: data; candidate: change-points candidates in Step 1
  k=length(candidate)
  psi=c(0,candidate,length(y)+1)
  for(j in 2:(k+1)){
    candid2=(psi[j-1]+2):(psi[j+1]-2)
    s=1:length(candid2)
    for(l in 1:length(candid2)){
      leftmean=mean(y[(psi[j-1]+1):candid2[l]])
      rightmean=mean(y[(candid2[l]+1):(psi[j+1]-1)])
      s[l]=(leftmean-rightmean)^2*length(y[(psi[j-1]+1):candid2[l]])*length(y[(candid2[l]+1):(psi[j+1]-1)])
    }
    psi[j]=candid2[which.max(s)]
  }
  return(psi[2:(k+1)])
}

test<-function(y,candidate2,alpha){          #Step3
  #candidate2:  refined change points candidates in Step 2
  #alpha:  significance level
  n=length(y)
  X1=matrix(1,nrow=n,ncol=n)
  for(i in 1:(n-1))
    X1[i,(i+1):n]=0
  X=X1[,2:n]
  yy=y-mean(y)     #centralized y
  XX=NULL           #centralized x
  for(i in 1:length(candidate2)){
    XX=cbind(XX,X[,candidate2[i]]-mean(X[,candidate2[i]]))
  }
  fit=ridge.proj(XX, yy, standardize=FALSE, lambda=1/n)        #ridge projection method
  candidate2[which(fit$pval.corr <alpha)]
}

###############################Additional Simulation
set.seed(18)
n=500
x=1:n
signal=(x>=0.3*n)-1*(x>=0.4*n)+1*(x>=0.8*n)-1*(x>=0.9*n)
lambda=0.45
h1=15
loc=NULL       #Step 1
Loc=NULL       #Step 2
LOC=NULL      #exhaustive search
time=NULL     #running time for Step 2
Time=NULL     #running time for exhaustive search
for(i in 1:200){
  y=signal+rnorm(n,mean=0,sd=0.5)
  D=abs(localDiagnostic(y, h=h1))
  index=localMax(D,span=h1)
  candidate=index[which(D[index]>lambda)]
  time1<- system.time({candidate2=refine(y,candidate)})
  time2<- system.time({candidate3=global(y,candidate)})
  time=c(time, time1[1])
  Time=c(Time, time2[1])
  loc=c(loc,candidate[which.min(abs(candidate-0.3*n))], candidate[which.min(abs(candidate-0.4*n))], candidate[which.min(abs(candidate-0.8*n))], candidate[which.min(abs(candidate-0.9*n))])
  Loc=c(Loc, candidate2[which.min(abs(candidate2-0.3*n))], candidate2[which.min(abs(candidate2-0.4*n))], candidate2[which.min(abs(candidate2-0.8*n))], candidate2[which.min(abs(candidate2-0.9*n))])
  LOC=c(LOC,candidate3[which.min(abs(candidate3-0.3*n))], candidate3[which.min(abs(candidate3-0.4*n))], candidate3[which.min(abs(candidate3-0.8*n))], candidate3[which.min(abs(candidate3-0.9*n))])
}
mean(time)
mean(Time)
par(mfrow=c(1,3))
hist(loc, breaks=80 ,xlab="locations", main="Step 1")
hist(Loc, breaks=80, xlab="locations", main="Step 2")
hist(LOC, breaks=80, xlab="locations", main="exhaustive search")


######omit Step 2
set.seed(54)
n=3000
J=50
x=1:n
tau=sort(sample(1:149, size=J, replace = FALSE))*20
h=c(8,10,12)
lambda=0.3*0.95^(1:4)
num=matrix(0,nr=3,ncol=4)       #average K in Step 1
Num=matrix(0,nr=3,ncol=4)       #average K in Step 3
FDP=matrix(0,nr=3,ncol=4)        #average FDP
TR=matrix(0,nr=3,ncol=4)         #average true detection rate
beta=sample(c(-1,1), size=J, replace = TRUE)
signal=0
for(i in 1:J){
  signal<-signal+beta[i]*(x>tau[i])
}
for(k in 1:3){
  h1=h[k]
  for(j in 1:4){
    n1=rep(0,200)
    n2=rep(0,200)
    fp=rep(0,200)
    tr=rep(0,200)
    for(i in 1:200){
      y=signal+rnorm(n,mean=0,sd=0.3) 
      D=abs(localDiagnostic(y, h=h1))
      index=localMax(D,span=h1)
      candidate=index[which(D[index]>lambda[j])]
      estimate=test(y,candidate,0.05)
      for(l in 1:length(estimate)){
        if(min(abs(estimate[l]-tau))>1)
          fp[i]=fp[i]+1
      }
      fp[i]=fp[i]/length(estimate)
      for(l in 1:J){
        if(min(abs(estimate-tau[l]))<=1)
          tr[i]=tr[i]+1/J
      }
      n1[i]=length(candidate)
      n2[i]=length(estimate)
    }
    FDP[k,j]=mean(fp)
    TR[k,j]=mean(tr)
    num[k,j]=mean(n1)
    Num[k,j]=mean(n2)
  }
}

#####################################Table 1
cbind(Num[1,],TR[1,],FDP[1,],Num[2,],TR[2,],FDP[2,],Num[3,],TR[3,],FDP[3,])
