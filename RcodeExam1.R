##This is R Code for Statistical Computing

## Problem 1

## Load MASS Library
library (MASS)
attach (Boston)

x<-Boston$medv
cv<-sd(x)/mean(x);cv  

N=2000
n<-length(x);n
cvb<-numeric(N)

##Bootstrap estimate of standard error of cv

for (k in 1:N) 
{
i<-sample(1:n, size=n, replace=TRUE)  

xbar<-mean(x[i])
s<-sd(x[i])
cvb[k]<-s/xbar
}

se.cv<-sd(cvb);se.cv 
hist(cvb,prob=TRUE)
m<-2000
cvm<-numeric(m)

for(k in 1:m)
{
  i<-sample(1:n, size=n, replace=TRUE)
  xbar<-mean(x[i])
  s<-sd(x[i])
  cvm[k]<-s/xbar
}

mean.cvm<-mean(cvm);mean.cvm ##Answer 1-1.

bias<-mean(cvm)-cv;bias

##Bootstrap Percentile confidence interval

alfa<-0.05
ci<-quantile (cvm,c(alfa/2,1-alfa/2),type=1);ci     ## Answer for 1-2


##Calculate Standard error for T Distribution

error<-qt (.975,df=n-1)*s/sqrt(n);error

cit.upper<-mean(cvm)+error;cit.upper
cit.bot<-mean(cvm)-error;cit.bot     ## Answer for 1-3



##Calculate Standard error for BCa Bootstrap confidence interval

for (k in 1:p)
{
  i<-sample (1:n, size=n, replace=TRUE)
  xbar<-mean(x[i])
  minus<-(x[i]-xbar)^2
  
  
  
}

library(boot)

boot(data=x, statistic=c, R=1000, conf=0.95)

boot.ci(cvm,type="bca")

