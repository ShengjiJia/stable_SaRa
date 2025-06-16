library(hdi)
library(DNAcopy)
library(mosum)
library(cumSeg)

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


#######Real Data
SNPdata <- read.delim("~/Research/Projects/change points(stable)/SNPdata.txt")
y=na.omit(SNPdata$X99HI0700A.Log.R.Ratio[which(SNPdata$Chr==20)])
y=y[-which(y<=-1.5)]
n=length(y)
y=y[1:n]
x=1:n
X1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1))
  X1[i,(i+1):n]=0
X=X1[,2:n]
###CSB
set.seed(12)
CBS=DNAcopy::segment(CNA(y, rep(1,n), 1:n))
num1=length(CBS$output[,4])-1
loc1=CBS$output[1:num1,4]
loc11=c(0,loc1,n)
estimate1=NULL
for(i in 1:(num1+1)){
  m=mean(y[(loc11[i]+1):loc11[i+1]])
  estimate1=c(estimate1,rep(m, loc11[i+1]-loc11[i]))
}
par(mfrow=c(2,3))
plot(y,xlab="locations",main="CBS",pch=20,col=8)
lines(estimate1,col=2,lwd=2)
###MOSUM
MOSUM=mosum(y, G=10)
loc2=MOSUM$cpts
num2=length(loc2)
loc22=c(0,loc2,n)
estimate2=NULL
for(i in 1:(num2+1)){
  m=mean(y[(loc22[i]+1):loc22[i+1]])
  estimate2=c(estimate2,rep(m, loc22[i+1]-loc22[i]))
}
plot(y,xlab="locations",main="MOSUM",pch=20,col=8)
lines(estimate2,col=2,lwd=2)
###SaRa-BIC
h1=20
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
loc3=sort(candidate1[1:j1])
num3=length(loc3)
loc33=c(0,loc3,n)
estimate3=NULL
for(i in 1:(num3+1)){
  m=mean(y[(loc33[i]+1):loc33[i+1]])
  estimate3=c(estimate3,rep(m, loc33[i+1]-loc33[i]))
}
plot(y,xlab="locations",main="SaRa-BIC",pch=20,col=8)
lines(estimate3,col=2,lwd=2)
###SaRa-mBIC
j2=which.min(mbic)
loc4=sort(candidate1[1:j2])
num4=length(loc4)
loc44=c(0,loc4,n)
estimate4=NULL
for(i in 1:(num4+1)){
  m=mean(y[(loc44[i]+1):loc44[i+1]])
  estimate4=c(estimate4,rep(m, loc44[i+1]-loc44[i]))
}
plot(y,xlab="locations",main="SaRa-mBIC",pch=20,col=8)
lines(estimate4,col=2,lwd=2)
###s-SaRa
candidate2=refine(y,candidate)
loc5=test(y,candidate2,0.05)
num5=length(loc5)
loc55=c(0,loc5,n)
estimate5=NULL
for(i in 1:(num5+1)){
  m=mean(y[(loc55[i]+1):loc55[i+1]])
  estimate5=c(estimate5,rep(m, loc55[i+1]-loc55[i]))
}
plot(y,xlab="locations",main="s-SaRa",pch=20,col=8)
lines(estimate5,col=2,lwd=2)
qqnorm(y-estimate5)
