# Clear working environment
rm(list=ls())

# Options for document compilation
knitr::opts_chunk$set(warning=FALSE, message=FALSE, comment=NA, fig.width=4, fig.height=3)
library(MASS) # load namespace `MASS`
data(Boston) # the object "Boston" now appears in the global environment
s <- sd(Boston$medv) # sample std dev
xbar <- mean(Boston$medv) # arithmetic mean
mycv <- s/xbar # coefficient of variation
mycv
B <- 1000 # we will do 1000 bootstrap replicates
n <- length(Boston$medv) # number of obs to be in each replicate
set.seed(8675309) # set random seed for repeatability
# make B replicate samples and place them in a list object
my_replicates <- lapply(1:B, function(x) sample(Boston$medv, size=n, replace=T))
theta_hats <- sapply(my_replicates, function(x) sd(x)/mean(x)) # coef of variation for each replicate
theta_hat_bar <- mean(theta_hats)
theta_hat_bar
theta_hat_errors <- theta_hats-theta_hat_bar
theta_hat_sq_errors <- theta_hat_errors^2
std_err_theta <- sqrt( (1/(B-1)) * sum(theta_hat_sq_errors) )
std_err_theta
alpha <- 0.05
z <- qnorm(alpha/2, lower.tail=F)
lower <- mycv-z*std_err_theta
upper <- mycv+z*std_err_theta
c(lower, upper)
# define function to take bootstrap standard error
bootstrap_se <- function(x, num_replicates, estimate_fun){
  # sample size
  n <- length(x)
  # pre-allocate empty vector for bootstrap estimates
  boot_ests <- vector(mode='numeric', length=num_replicates)
  # take bootstrap samples and find estimate of each
  for(i in 1:num_replicates){
    boot_ests[i] <- estimate_fun(sample(x, size=n, replace=T))
  }
  # find mean of bootstrap estimates
  boot_est_mean <- mean(boot_ests)
  # find squared error of each bootstrap estimate
  boot_sq_errs <- (boot_ests-boot_est_mean)^2
  # find std error
  boot_std_err <- sqrt( (1/(num_replicates-1)) * sum(boot_sq_errs) )
  
  return(boot_std_err)
}

# define coef_var function to pass to the estimate_fun argument of bootstrap_se()
coef_var <- function(x) sd(x)/mean(x)

# apply bootstrap_se() to each bootstrap sample in my_replicates
se_hats <- sapply(my_replicates, bootstrap_se, num_replicates=200, estimate_fun=coef_var)
my_ts <- (theta_hats-mycv)/(se_hats)
t_lower <- quantile(my_ts, (1-alpha/2))
t_upper <- quantile(my_ts, alpha/2)
c(t_lower, t_upper)
lower <- mycv-t_lower*std_err_theta
upper <- mycv-t_upper*std_err_theta
c(lower, upper)
z_0hat <- qnorm(sum(theta_hats>mycv)/B)
z_0hat
loo_theta_hats <- vector(mode='numeric', length=n)
for(i in 1:n){
  loo_theta_hats[i] <- coef_var(Boston$medv[-i])
}
loo_theta_hats_mean <- mean(loo_theta_hats)
loo_thete_hat_errs <- loo_theta_hats_mean-loo_theta_hats
numerator <- sum((loo_thete_hat_errs)^3)
denominator <- 6*sum(loo_thete_hat_errs^2)^(3/2)
ahat <- numerator/denominator
ahat
tmp <- z_0hat+qnorm(alpha/2)
tmp2 <- z_0hat+((tmp)/(1-(ahat*tmp)))
alpha_1 <- pnorm(tmp2)
alpha_1
tmp <- z_0hat+qnorm((1-(alpha/2)))
tmp2 <- z_0hat+((tmp)/(1-(ahat*tmp)))
alpha_2 <- pnorm(tmp2)
alpha_2
quantile(theta_hats, c(alpha_1, alpha_2))
library(DAAG) # load library for ironslag data set
data(ironslag) # data now exists in the global work space
k <- 10 # for 10-fold cross-validation
n <- nrow(ironslag) # sample size
seg_size <- floor(n/k) # approx size of each segment
leftovers <- n%%k
seg_assignments <- rep(1:k, each=seg_size) # create segment assignments
seg_assignments <- c(seg_assignments, 1:leftovers) # add "leftover" observations
set.seed(8675309) # set random seed
seg_assignments <- sample(seg_assignments) # randomly permute segment assignments
train <- ironslag[seg_assignments!=1, ]
lin_mod <- lm(magnetic~chemical, data=train) # Fit linear model using training data
quad_mod <- lm(magnetic~chemical+I(chemical^2), data=train ) # Fit quadratic model using training data
exp_mod <- lm(log(magnetic)~chemical, data=train) # Fit exponential model using training data
log_mod <- lm(log(magnetic)~log(chemical), data=train) # Fit log-log model using training data
validation <- ironslag[seg_assignments==1, ]
lin_pred <- predict(lin_mod, newdata=validation) # use linear model to predict validation set
quad_pred <- predict(quad_mod, newdata=validation) # use quadratic model to predict validation set
exp_pred <- exp(predict(exp_mod, newdata=validation)) # remember to exponentiate predictions
log_pred <- exp(predict(log_mod, newdata=validation)) # remember to exponentiate predictions
lin_errors <- (validation$magnetic-lin_pred)
quad_errors <- (validation$magnetic-quad_pred)
exp_errors <- (validation$magnetic-exp_pred)
log_errors <- (validation$magnetic-log_pred)

lin_MSE <- mean(lin_errors^2)
quad_MSE <- mean(quad_errors^2)
exp_MSE <- mean(exp_errors^2)
log_MSE <- mean(log_errors^2)

c(lin_MSE, quad_MSE, exp_MSE, log_MSE)
mse_df <- data.frame(lin_MSE=NA, quad_MSE=NA, exp_MSE=NA, log_MSE=NA)
for(i in 1:k){
  train <- ironslag[seg_assignments!=i, ] # define training set i
  
  # fit models for training set i
  lin_mod <- lm(magnetic~chemical, data=train) # linear model
  quad_mod <- lm(magnetic~chemical+I(chemical^2), data=train) # quadratic model 
  exp_mod <- lm(log(magnetic)~chemical, data=train) # exponential model
  log_mod <- lm(log(magnetic)~log(chemical), data=train) # log-log model
  
  validation <- ironslag[seg_assignments==i, ] # define validation set i
  
  # make predictions for validation set i
  lin_pred <- predict(lin_mod, newdata=validation) # linear model predictions
  quad_pred <- predict(quad_mod, newdata=validation) # quadratic model predictions
  exp_pred <- exp(predict(exp_mod, newdata=validation)) # exponential model predictions
  log_pred <- exp(predict(log_mod, newdata=validation)) # log-log model predictions
  
  # calculate prediction errors for validation set i
  lin_errors <- (validation$magnetic-lin_pred)
  quad_errors <- (validation$magnetic-quad_pred)
  exp_errors <- (validation$magnetic-exp_pred)
  log_errors <- (validation$magnetic-log_pred)
  
  # calculate MSE for validation set i
  lin_MSE <- mean(lin_errors^2)
  quad_MSE <- mean(quad_errors^2)
  exp_MSE <- mean(exp_errors^2)
  log_MSE <- mean(log_errors^2)
  
  # store calculated MSEs in mse_df
  mse_df[i, ] <- c(lin_MSE, quad_MSE, exp_MSE, log_MSE)
}
colMeans(mse_df)
# load library `bootstrap` for the scor data set
library(bootstrap)
data(scor) # call the data set into the global environment
myfunc <- function(mysamp){
  sigma_hat <- cov(mysamp) # MLE of Sigma
  evals <- eigen(sigma_hat)$values # eigenvalues of Sigma hat
  lambda_1 <- max(evals) # largest eigenvalue
  lambda_sum <- sum(evals) # sum of all eigenvalues
  theta_hat <- lambda_1/lambda_sum # proportion of variance explained by first principal component
  return(theta_hat)
}
my_theta_hat <- myfunc(scor) # find theta hat using all x's
my_theta_hat
theta_hat_vec <- vector(mode='numeric', length=nrow(scor)) # pre-allocate empty vector for theta hats
for(i in 1:nrow(scor)){
  samp_x_i <- scor[-i, ] # leave out observation i 
  theta_hat_i <- myfunc(samp_x_i) # calculate theta hat without obs i, using `myfunc`
  theta_hat_vec[i] <- theta_hat_i # assign computed theta hat to ith element in "theta_hat_vec"
}
theta_hat_bar <- mean(theta_hat_vec) # mean of all theta hat i's
theta_hat_bar
n <- nrow(scor) # sample size
jacknife_bias <- (n-1)*(theta_hat_bar-my_theta_hat) # jacknife est of bias
jacknife_bias
theta_hat_errors <- theta_hat_vec-theta_hat_bar # error of each theta hat i
theta_hat_sq_errors <- theta_hat_errors^2 # squared errors from raw errors
jacknife_se <- sqrt( ((n-1)/n) * sum(theta_hat_sq_errors) ) # compute jacknife se
jacknife_se