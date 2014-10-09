install.packages ("rmarkdown")
install.packages ("knitr")
##This is R Code for Statistical Computing

## Problem 1

# Clear working environment
rm(list=ls())

## Load MASS Library
library (MASS)
attach (Boston)

x<-Boston$medv
cv<-sd(x)/mean(x);cv ##Answer 1-1.

N=2000
n<-length(x);n

##Bootstrap estimate of standard error of cv
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

sd(cvm)

mean.cvm<-mean(cvm);mean.cvm 

bias<-mean(cvm)-cv;bias

##Bootstrap Percentile confidence interval

alfa<-0.05
ci<-quantile (cvm,c(alfa/2,1-alfa/2),type=1);ci     ## Answer for 1-2


##Calculate Standard error for T Distribution

error<-qt (.975,df=n-1)*s/sqrt(n);error

##these code does not compatible with it self, this is along with BOB's code
cit.upper<-mean(cvm)+qt (.975,df=n-1)*std_err_theta/sqrt(n) ;cit.upper
cit.bot<-mean(cvm)-qt (.975,df=n-1)*std_err_theta/sqrt(n) ;cit.bot

cit.ci=c(cit.bot,cit.upper);cit.ci
## Answer for 1-3



##Calculate Standard error for BCa Bootstrap confidence interval

for (k in 1:p)
{
  i<-sample (1:n, size=n, replace=TRUE)
  xbar<-mean(x[i])
  minus<-(x[i]-xbar)^2
}

library(boot)


f<-function(data,i)
{
  i<-sample(1:n, size=n, replace=TRUE)
  xbar<-mean(x[i])
  s<-sd(x[i])
  cvm<-s/xbar
  
}

bca.err<-boot(data=x, statistic=f , R=5000)

boot.ci(bca.err,type="bca")  ##Answer 1-4.


## Problem 2 --------------------------------------------------

##Load Library(DAAG) as indicated in example 7.18 page 210 in Rizzo's book.

install.packages("DAAG")
library (DAAG)

##Create 10 equally size folds

install.packages("cvTools")

library(cvTools) #run the above line if you don't have this library

folds.k<- 10 #the number of folds

folds<- cvFolds(NROW(ironslag), K=folds.k)

ironslag$holdoutpred_lin <- rep(0,nrow(ironslag))
ironslag$holdoutpred_quad <- rep(0,nrow(ironslag))
ironslag$holdoutpred_exp <- rep(0,nrow(ironslag))
ironslag$holdoutpred_log <- rep(0,nrow(ironslag))

## Calculate MSE for Linear
for(i in 1:k){
  train <- ironslag[folds$subsets[folds$which != i], ] #Set the training set
  validation <- ironslag[folds$subsets[folds$which == i], ] #Set the validation set
  new_lin <- lm(magnetic ~ chemical,data=train) #Set linear model 
  new_quad <- lm(magnetic ~ chemical + I(chemical^2),data=train) #Set Quadratic model 
  new_exp <- lm(log(magnetic) ~ chemical,data=train) #Set Exponential model 
  new_log <- lm(log(magnetic) ~ log(chemical),data=train) #Set Log-Log model 
  
  newpred_lin <- predict(new_lin,newdata=validation) #Get the predicitons for the validation set
  newpred_quad <- predict(new_quad,newdata=validation) #Get the predicitons for the validation set
  newpred_exp <- predict(new_exp,newdata=validation) #Get the predicitons for the validation set
  newpred_log <- predict(new_log,newdata=validation) #Get the predicitons for the validation set
  
  ironslag[folds$subsets[folds$which == i], ]$holdoutpred_lin <- newpred_lin #Put the hold out prediction in the data set for later use
  ironslag[folds$subsets[folds$which == i], ]$holdoutpred_quad <- newpred_quad #Put the hold out prediction in the data set for later use
  ironslag[folds$subsets[folds$which == i], ]$holdoutpred_exp <- newpred_exp #Put the hold out prediction in the data set for later use
  ironslag[folds$subsets[folds$which == i], ]$holdoutpred_log <- newpred_log #Put the hold out prediction in the data set for later use
}
ironslag$holdoutpred_exp #Oringinal Dataset with Predictions
ironslag$holdoutpred_quad #Oringinal Dataset with Predictions

mse.lin<-mean((ironslag$holdoutpred_lin-ironslag$magnetic)^2) ##Calculate MSE for Linear Model
mse.quad<-mean((ironslag$holdoutpred_quad-ironslag$magnetic)^2) ##Calculate MSE for Linear Model
mse.exp<-mean((ironslag$holdoutpred_exp-ironslag$magnetic)^2) ##Calculate MSE for Linear Model
mse.log<-mean((ironslag$holdoutpred_log-ironslag$magnetic)^2) ##Calculate MSE for Linear Model

mse<-c(mse.lin, mse.quad,mse.exp,mse.log);mse



