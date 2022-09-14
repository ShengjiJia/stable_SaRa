#######################Simulation 0 (Figure 1 in Introduction)
#######################define some functions 
#Calculate the value for local diagnostic function
localDiagnostic<-function (y, h) { 
  yy = c(rep(0, h - 1), y, rep(0, h))
  n = length(y)
  z = rep(0, n)
  for(i in 1:n){
    z[i]=sum(yy[i:(h+i-1)])/h-sum(yy[(h+i):(2*h-1+i)])/h
  }
  return(z)
}

#Get the local maximizers of local diagnostic function
localMax<-function (y, span = 5) {  
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


#######################Simulation 0 
set.seed(9)
n=1000
J=20
x=1:n
tau=sort(sample(1:49, size=J, replace = FALSE))*20
h=c(6,8,10)
lambda=0.35*0.9^((1:5)-1)
num=matrix(0,nr=3,ncol=5)       #average number of candidates
beta=sample(c(-0.6,-0.4,0.4,0.6), size=J, replace = TRUE)
signal=0
for(i in 1:J){
  signal<-signal+beta[i]*(x>tau[i])
}
for(i in 1:3){
  h1=h[i]
  for(j in 1:5){
    n1=rep(0,200)
    for(k in 1:200){
      y=signal+rnorm(n,mean=0,sd=0.25) 
      D=abs(localDiagnostic(y, h=h1))
      index=localMax(D,span=h1)
      n1[k]=sum(D[index]>lambda[j])
    }
    num[i,j]=mean(n1)
  }
}

############################Figure 1
par(mfrow=c(1,2))
plot(x,signal,type="l", lwd=2, xlab="locations", ylab="signal", main="true signal")
plot(x=lambda, y=num[1, ], pch=1, col=3, type="b", ylim=c(10,50), xlab = "lambda",ylab ="K",main="estimated K")
lines(x=lambda, y=num[2, ], pch=2, col=2, type="b")
lines(x=lambda, y=num[3, ], pch=3, col=1, type="b")
legend("topright", legend=c("h=6","h=8","h=10"), pch=c(1,2,3), col=c(3,2,1), bty="n")