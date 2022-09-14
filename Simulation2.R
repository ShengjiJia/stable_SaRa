library(hdi)
library(DNAcopy)
library(tilingArray)

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

refine<-function(y,candidate){          #Step2
  #y: response or data
  #candidate: change points candidates after screening in Step 1
  k1=length(candidate)
  psi=c(0,candidate,length(y)+1)
  for(iter in 1:10){
    r=sample(2:(k1+1),size=k1,replace=FALSE)
    for(j in 1:k1){
      a=r[j]
      b=min(psi[a+1]-psi[a]-2,psi[a]-psi[a-1]-2,5)
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

###############################Simulation2
set.seed(12)
n=3000
J=50
x=1:n
X1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1))
  X1[i,(i+1):n]=0
X=X1[,2:n]
tau=sort(sample(1:149, size=J, replace = FALSE))*20
beta=sample(c(-0.75,-0.5,-0.25,0.25,0.5,0.75), size=J, replace = TRUE)
signal=0
for(i in 1:J){
  signal<-signal+beta[i]*(x>tau[i])
}
weaktau=tau[which(abs(beta)==0.25)]         #change points with weak signals
Num=matrix(0,nr=5,ncol=4)         #number of change points for 4 cases and 5 methods
tr=matrix(0,nr=5,ncol=4)              #true detection rate for 4 cases and 5 methods

#####case 1 i.i.d. errors
num=matrix(0,nr=5,ncol=200)       #final estimated number  
TR=matrix(0,nr=5,ncol=200)          #true detection rate
weakTR=matrix(0,nr=5,ncol=200)  #true detection rate of weak signals
for(i in 1:200){
  y=signal+rnorm(n,mean=0,sd=0.2) 
  ######CBS
  CBS=DNAcopy::segment(CNA(y, rep(1,n), 1:n))
  num[1,i]=length(CBS$output[,4])-1
  estimate1=CBS$output[1:(length(CBS$output[,4])-1),4]
  TR[1,i]=0 
  for(l in 1:J){
    if(min(abs(estimate1-tau[l]))<=3)
      TR[1,i]=TR[1,i]+1/J
  }
  weakTR[1,i]=0 
  for(l in 1:length(weaktau)){
    if(min(abs(estimate1-weaktau[l]))<=3)
      weakTR[1,i]=weakTR[1,i]+1/length(weaktau)
  }
  ######DP
  DP=tilingArray::segment(y, maxseg=100, maxk=n/10)
  num[2,i]=which.max(logLik(DP, penalty="BIC"))-1
  estimate2=DP@breakpoints[[which.max(logLik(DP, penalty="BIC"))]][,"estimate"]
  TR[2,i]=0 
  for(l in 1:J){
    if(min(abs(estimate2-tau[l]))<=3)
      TR[2,i]=TR[2,i]+1/J
  }
  weakTR[2,i]=0 
  for(l in 1:length(weaktau)){
    if(min(abs(estimate2-weaktau[l]))<=3)
      weakTR[2,i]=weakTR[2,i]+1/length(weaktau)
  }
  ######SaRa-BIC
  h1=10
  D=abs(localDiagnostic(y, h=h1))
  index=localMax(D,span=h1)
  lambda=2*sqrt(2/h1)*mad(diff(y))/sqrt(2)        #use robust estimator of sigma
  candidate=index[which(D[index]>lambda)]
  candidate1=candidate[order(D[candidate],decreasing=TRUE)]          #ranking
  bic=rep(0,length(candidate1))
  mbic=rep(0,length(candidate1))
  for(j in 1:length(candidate1)){
    model=lm(y~X[,candidate1[1:j]])
    bic[j]=n*log(summary(model)$sigma)+j*log(n)
    mbic[j]=n*log(summary(model)$sigma)+1.5*j*log(n)+0.5*sum(log(diff(c(0,sort(candidate1[1:j]),n))/n))
  }
  j1=which.min(bic)
  estimate3=candidate1[1:j1]
  num[3,i]=length(estimate3)
  TR[3,i]=0 
  for(l in 1:J){
    if(min(abs(estimate3-tau[l]))<=3)
      TR[3,i]=TR[3,i]+1/J
  }
  weakTR[3,i]=0 
  for(l in 1:length(weaktau)){
    if(min(abs(estimate3-weaktau[l]))<=3)
      weakTR[3,i]=weakTR[3,i]+1/length(weaktau)
  }
  ######SaRa-mBIC
  j2=which.min(mbic)
  estimate4=candidate1[1:j2]
  num[4,i]=length(estimate4)
  TR[4,i]=0 
  for(l in 1:J){
    if(min(abs(estimate4-tau[l]))<=3)
      TR[4,i]=TR[4,i]+1/J
  }
  weakTR[4,i]=0 
  for(l in 1:length(weaktau)){
    if(min(abs(estimate4-weaktau[l]))<=3)
      weakTR[4,i]=weakTR[4,i]+1/length(weaktau)
  }
  ######proposed s-SaRa
  candidate2=refine(y,candidate)
  estimate5=test(y,candidate2,0.05)
  num[5,i]=length(estimate5)
  TR[5,i]=0 
  for(l in 1:J){
    if(min(abs(estimate5-tau[l]))<=3)
      TR[5,i]=TR[5,i]+1/J
  }
  weakTR[5,i]=0 
  for(l in 1:length(weaktau)){
    if(min(abs(estimate5-weaktau[l]))<=3)
      weakTR[5,i]=weakTR[5,i]+1/length(weaktau)
  }
}
Num[,1]=apply(num,1,mean)
tr[,1]=apply(TR,1,mean)

#######case 2 AR(1) errors
num=matrix(0,nr=5,ncol=200)      #final estimated number  
TR=matrix(0,nr=5,ncol=200)         #true detection rate
for(i in 1:200){
  y=signal+arima.sim(n=(length(x)+100), list(ar = 0.2), sd = 0.2)[101:(length(x)+100)]
  ######CBS
  CBS=DNAcopy::segment(CNA(y, rep(1,n), 1:n))
  num[1,i]=length(CBS$output[,4])-1
  estimate1=CBS$output[1:(length(CBS$output[,4])-1),4]
  TR[1,i]=0 
  for(l in 1:J){
    if(min(abs(estimate1-tau[l]))<=3)
      TR[1,i]=TR[1,i]+1/J
  }
  ######DP
  DP=tilingArray::segment(y, maxseg=100, maxk=n/10)
  num[2,i]=which.max(logLik(DP, penalty="BIC"))-1
  estimate2=DP@breakpoints[[which.max(logLik(DP, penalty="BIC"))]][,"estimate"]
  TR[2,i]=0 
  for(l in 1:J){
    if(min(abs(estimate2-tau[l]))<=3)
      TR[2,i]=TR[2,i]+1/J
  }
  ######SaRa-BIC
  h1=10
  D=abs(localDiagnostic(y, h=h1))
  index=localMax(D,span=h1)
  lambda=2*sqrt(2/h1)*mad(diff(y))/sqrt(2)        #use robust estimator of sigma
  candidate=index[which(D[index]>lambda)]
  candidate1=candidate[order(D[candidate],decreasing=TRUE)]          #ranking
  bic=rep(0,length(candidate1))
  mbic=rep(0,length(candidate1))
  for(j in 1:length(candidate1)){
    model=lm(y~X[,candidate1[1:j]])
    bic[j]=n*log(summary(model)$sigma)+j*log(n)
    mbic[j]=n*log(summary(model)$sigma)+1.5*j*log(n)+0.5*sum(log(diff(c(0,sort(candidate1[1:j]),n))/n))
  }
  j1=which.min(bic)
  estimate3=candidate1[1:j1]
  num[3,i]=length(estimate3)
  TR[3,i]=0 
  for(l in 1:J){
    if(min(abs(estimate3-tau[l]))<=3)
      TR[3,i]=TR[3,i]+1/J
  }
  ######SaRa-mBIC
  j2=which.min(mbic)
  estimate4=candidate1[1:j2]
  num[4,i]=length(estimate4)
  TR[4,i]=0 
  for(l in 1:J){
    if(min(abs(estimate4-tau[l]))<=3)
      TR[4,i]=TR[4,i]+1/J
  }
  ######proposed s-SaRa
  candidate2=refine(y,candidate)
  estimate5=test(y,candidate2,0.05)
  num[5,i]=length(estimate5)
  TR[5,i]=0 
  for(l in 1:J){
    if(min(abs(estimate5-tau[l]))<=3)
      TR[5,i]=TR[5,i]+1/J
  }
}
Num[,2]=apply(num,1,mean)
tr[,2]=apply(TR,1,mean)

#######case 3 heteroscedastic errors
num=matrix(0,nr=5,ncol=200)      #final estimated number  
TR=matrix(0,nr=5,ncol=200)         #true detection rate
for(i in 1:200){
  y=signal+rnorm(n,mean=0,sd=0.1*(1+2*(1:n)/n))
  ######CBS
  CBS=DNAcopy::segment(CNA(y, rep(1,n), 1:n))
  num[1,i]=length(CBS$output[,4])-1
  estimate1=CBS$output[1:(length(CBS$output[,4])-1),4]
  TR[1,i]=0 
  for(l in 1:J){
    if(min(abs(estimate1-tau[l]))<=3)
      TR[1,i]=TR[1,i]+1/J
  }
  ######DP
  DP=tilingArray::segment(y, maxseg=100, maxk=n/10)
  num[2,i]=which.max(logLik(DP, penalty="BIC"))-1
  estimate2=DP@breakpoints[[which.max(logLik(DP, penalty="BIC"))]][,"estimate"]
  TR[2,i]=0 
  for(l in 1:J){
    if(min(abs(estimate2-tau[l]))<=3)
      TR[2,i]=TR[2,i]+1/J
  }
  ######SaRa-BIC
  h1=10
  D=abs(localDiagnostic(y, h=h1))
  index=localMax(D,span=h1)
  lambda=2*sqrt(2/h1)*mad(diff(y))/sqrt(2)        #use robust estimator of sigma
  candidate=index[which(D[index]>lambda)]
  candidate1=candidate[order(D[candidate],decreasing=TRUE)]          #ranking
  bic=rep(0,length(candidate1))
  mbic=rep(0,length(candidate1))
  for(j in 1:length(candidate1)){
    model=lm(y~X[,candidate1[1:j]])
    bic[j]=n*log(summary(model)$sigma)+j*log(n)
    mbic[j]=n*log(summary(model)$sigma)+1.5*j*log(n)+0.5*sum(log(diff(c(0,sort(candidate1[1:j]),n))/n))
  }
  j1=which.min(bic)
  estimate3=candidate1[1:j1]
  num[3,i]=length(estimate3)
  TR[3,i]=0 
  for(l in 1:J){
    if(min(abs(estimate3-tau[l]))<=3)
      TR[3,i]=TR[3,i]+1/J
  }
  ######SaRa-mBIC
  j2=which.min(mbic)
  estimate4=candidate1[1:j2]
  num[4,i]=length(estimate4)
  TR[4,i]=0 
  for(l in 1:J){
    if(min(abs(estimate4-tau[l]))<=3)
      TR[4,i]=TR[4,i]+1/J
  }
  ######proposed s-SaRa
  candidate2=refine(y,candidate)
  estimate5=test(y,candidate2,0.05)
  num[5,i]=length(estimate5)
  TR[5,i]=0 
  for(l in 1:J){
    if(min(abs(estimate5-tau[l]))<=3)
      TR[5,i]=TR[5,i]+1/J
  }
}
Num[,3]=apply(num,1,mean)
tr[,3]=apply(TR,1,mean)

#######case 4 t-distribution
num=matrix(0,nr=5,ncol=200)      #final estimated number  
TR=matrix(0,nr=5,ncol=200)         #true detection rate
for(i in 1:200){
  y=signal+0.1*rt(n,df=3)
  ######CBS
  CBS=DNAcopy::segment(CNA(y, rep(1,n), 1:n))
  num[1,i]=length(CBS$output[,4])-1
  estimate1=CBS$output[1:(length(CBS$output[,4])-1),4]
  TR[1,i]=0 
  for(l in 1:J){
    if(min(abs(estimate1-tau[l]))<=3)
      TR[1,i]=TR[1,i]+1/J
  }
  ######DP
  DP=tilingArray::segment(y, maxseg=100, maxk=n/10)
  num[2,i]=which.max(logLik(DP, penalty="BIC"))-1
  estimate2=DP@breakpoints[[which.max(logLik(DP, penalty="BIC"))]][,"estimate"]
  TR[2,i]=0 
  for(l in 1:J){
    if(min(abs(estimate2-tau[l]))<=3)
      TR[2,i]=TR[2,i]+1/J
  }
  ######SaRa-BIC
  h1=10
  D=abs(localDiagnostic(y, h=h1))
  index=localMax(D,span=h1)
  lambda=2*sqrt(2/h1)*mad(diff(y))/sqrt(2)        #use robust estimator of sigma
  candidate=index[which(D[index]>lambda)]
  candidate1=candidate[order(D[candidate],decreasing=TRUE)]          #ranking
  bic=rep(0,length(candidate1))
  mbic=rep(0,length(candidate1))
  for(j in 1:length(candidate1)){
    model=lm(y~X[,candidate1[1:j]])
    bic[j]=n*log(summary(model)$sigma)+j*log(n)
    mbic[j]=n*log(summary(model)$sigma)+1.5*j*log(n)+0.5*sum(log(diff(c(0,sort(candidate1[1:j]),n))/n))
  }
  j1=which.min(bic)
  estimate3=candidate1[1:j1]
  num[3,i]=length(estimate3)
  TR[3,i]=0 
  for(l in 1:J){
    if(min(abs(estimate3-tau[l]))<=3)
      TR[3,i]=TR[3,i]+1/J
  }
  ######SaRa-mBIC
  j2=which.min(mbic)
  estimate4=candidate1[1:j2]
  num[4,i]=length(estimate4)
  TR[4,i]=0 
  for(l in 1:J){
    if(min(abs(estimate4-tau[l]))<=3)
      TR[4,i]=TR[4,i]+1/J
  }
  ######proposed s-SaRa
  candidate2=refine(y,candidate)
  estimate5=test(y,candidate2,0.05)
  num[5,i]=length(estimate5)
  TR[5,i]=0 
  for(l in 1:J){
    if(min(abs(estimate5-tau[l]))<=3)
      TR[5,i]=TR[5,i]+1/J
  }
}
Num[,4]=apply(num,1,mean)
tr[,4]=apply(TR,1,mean)

#################################Table 2
cbind(Num[,1],tr[,1],Num[,2],tr[,2],Num[,3],tr[,3],Num[,4],tr[,4])

#################################true detecion rate (TR) for weak signals for case 1
apply(weakTR,1,mean)[3:5]