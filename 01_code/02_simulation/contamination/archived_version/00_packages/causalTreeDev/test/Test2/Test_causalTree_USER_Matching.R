# This program tests out the features of the causalTree package
#'Splitting Rule 1, Matching, Outlier'

#install.packages("devtools")
 library(devtools)
#install.packages("rpart", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("rpart.plot", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("reshape2", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("plyr", dependencies=TRUE, repos='http://cran.us.r-project.org')
library(rpart)
library(rpart.plot)
# install_github("susanathey/causalTree",force=TRUE)
#library(causalTree)
library(reshape2)
library(plyr)

set.seed(124)
# Generate data 
# parameters for data generating
p <- 20 # number of total covariates
pt <- 4 # number of covariates affecting treatment effects
py <- 4 # number of covariates affecting outcomes but not treatment effects
asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
n <- 3000 # total size of the dataset

propens <- .5 #treatment probability
sig = .01
treatsize <- .5 # treatment effect size
levsize <- 1

# draw W
w <- rbinom(n, 1, propens)
w_ <- 1-w


# data with outliers?
outliers = F
# True 
# False
# outlier scenario according to fdaoutlier (https://cran.r-project.org/web/packages/fdaoutlier/vignettes/simulation_models.html)
outlier_scenario = 3
# 1: outliers are shifted from the ‘normal’ non-outlying observations
# 2: non-persistent magnitude outliers, i.e., the outliers are magnitude outliers for only a portion of the domain of the functional data
# 3: outliers that are magnitude outliers for a part of the domain
# 4: outliers defined on the reversed interval of the main model
# 5: shape outliers with a different covariance structure from that of the main model
# 6: shape outliers that have a different shape for a portion of the domain
# 7: pure shape outliers that are periodic
# 8: pure shape outliers that are periodic
# 9: Periodic functions with outliers of different amplitude
# dist for rng if no outliers
dist_scenario_no_outliers = 'exp'
# norm
# exp
# cauchy
# lnorm
# unif



draw_X <- function(outliers=F, n, p, outlier_scenario=1, dist_scenario_no_outliers='norm'){
  if(outliers==F){
    
    dist_name <- paste('r', dist_scenario_no_outliers, sep='')
    rng_name <- match.fun(dist_name) 
    X <- matrix(rng_name(n * p), nrow = n, ncol = p)
    cat("Data without outliers. You've used ", dist_name, " for data generation.")
    
    
  } else {
    
    library(fdaoutlier)
    
    model_name <- paste('simulation_model', outlier_scenario, sep='')
    func_name <- match.fun(model_name) 
    outlier_sim <- func_name(n = n,
                             p = p,
                             outlier_rate = .15,
                             seed = 50, plot = F)
    X <- as.matrix(outlier_sim$data)
    cat("Data with outliers. You've used ", model_name, " of 'fdaoutlier' for outlier generation.")
    
  }
  
  return(X)
}

X <- draw_X(outliers, n, p, outlier_scenario, dist_scenario_no_outliers='unif')

# draw X

X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# generate treatment effects as function of X
if (p<pt+py) print("error: p>=pt+py required")
tau <- 0
for (iii in 1:pt) {
  tau <- tau + treatsize*pmax(X[,iii],array(0,n))*(2*X[,iii])
}

# generate average value of outcomes
mu <- treatsize*rowSums(X[,1:pt])+levsize*rowSums(X[,(pt+1):(pt+py)])

# generate outcomes as function of treatment status, mu, tau, and noise
y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)
y_ <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)

# create formulas for estimation
# if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
# create a formula such as f with the list of x variables separated by "+"
f <- ""
nextx <- ""
if (p>1) {
  for (ii in 1:(p-1)) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
    f <- paste(f, nextx, "+", sep="")
  }
  f <- paste(f, "x", ii+1, sep="")
} else if (p==1) {
  f <- "x1"
}

for (ii in 1:p) {
  nextx <- paste("x",ii, sep="")
  if (ii==1) {name <- nextx}
  if (ii>1) {name <- c(name, nextx)}
}

name <- c( name,  "y", "w", "tau_true")

tau_true <- (1-2*w)*(y_ - y)

ntr <- round(.333*n)
nest <- round(.333*n)
ntest <- n - ntr - nest

# set global parameters
minsize.temp = 25
split.Bucket.temp = F
bucketNum.temp = 5
bucketMax.temp = 100

X<-data.frame(X)
for (tmp1 in 1:ncol(X)){
  xtmp<-X[,tmp1]
  unxtmp<-unique(xtmp)
  if(length(unxtmp)<bucketMax.temp) #convert to factor
    X[,tmp1]<-factor(X[,tmp1])
}

dataTrain <- data.frame(X[1:ntr,], y[1:ntr], w[1:ntr], tau_true[1:ntr])
dataEst <- data.frame(X[(ntr+1):(ntr+nest),], y[(ntr+1):(ntr+nest)], w[(ntr+1):(ntr+nest)], tau_true[(ntr+1):(ntr+nest)])
dataTest <- data.frame(X[(ntr+nest+1):n,], y[(ntr+nest+1):n], w[(ntr+nest+1):n], tau_true[(ntr+nest+1):n])


###
names(dataTrain)=name
names(dataEst)=name
names(dataTest)=name

tree_honest_prune_list <- vector(mode="list", length=5)
tree_dishonest_prune_list <- vector(mode="list", length=5)

# preselect cross-validation groups to remove randomness in comparing methods
xvalvec = sample(5, nrow(dataTrain), replace=TRUE)
#xvalvec=5


# Do causal tree estimation
split.Rule.temp = "CT" #CT
cv.option.temp = "matching" #CT
split.Honest.temp = T
cv.Honest.temp = T
split.alpha.temp = .5
cv.alpha.temp = .5



#This function is a wrapper for honest causal tree
tree <- causalTree::honest.causalTree(as.formula(paste("y~",paste(f))),
                  data=dataTrain, treatment=dataTrain$w,
                  est_data=dataEst, est_treatment=dataEst$w,
                  split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp,
                  bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp,
                  split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
#You can still prune as usual; the cptable is the one from training the tree
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp)

# save the results
tree_honest_CT <- tree
tree_honest_CT_prune <- tree_prune

#manually convert honest to dishonest: honest causaltree re-estimation using training data itself
tree <- causalTree::honest.causalTree(as.formula(paste("y~",paste(f))),
                          data=dataTrain, treatment=dataTrain$w,
                          est_data=dataTrain, est_treatment=dataTrain$w,
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp,
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp,
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)


# get the dishonest version--estimated leaf effects on training sample
tree <- causalTree::causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                   split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                   bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                   split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 
tree_dishonest_CT_prune <- tree_prune

# we can get the honest tree manually, by estimating the leaf effects on a new sample
tree_honest_CT_prune2 <- causalTree::estimate.causalTree(object=tree_dishonest_CT_prune,data=dataEst, treatment=dataEst$w)

print(tree_honest_CT_prune)

#add to list
tree_dishonest_prune_list[[1]] <- tree_dishonest_CT_prune
tree_honest_prune_list[[1]] <- tree_honest_CT_prune



# Other honest splitting rules, cv

#tstats
split.Rule.temp = "tstats"
cv.option.temp = "matching"
split.Honest.temp = F
cv.Honest.temp = T
split.alpha.temp = 1

tree <- causalTree::honest.causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                          est_data=dataEst, est_treatment=dataEst$w,
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
# prune
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_honest_tstats <- tree
tree_honest_tstats_prune <- tree_prune
tree_honest_prune_list[[2]] <- tree_honest_tstats_prune

#honest user
split.Rule.temp = "user"
cv.option.temp = "matching"
split.Honest.temp = F
cv.Honest.temp = T
split.alpha.temp = 1
?causalTree

tree <- causalTree::honest.causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                          est_data=dataEst, est_treatment=dataEst$w,
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
#You can still prune as usual; the cptable is the one from training the tree
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_honest_USER <- tree
tree_honest_USER_prune <- tree_prune

tree_honest_prune_list[[3]] <- tree_honest_USER_prune

#honest TOT
split.Rule.temp = "TOT"
cv.option.temp = "matching"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = 1

tree <- causalTree::honest.causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                          est_data=dataEst, est_treatment=dataEst$w,
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
#You can still prune as usual; the cptable is the one from training the tree
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_honest_TOT <- tree
tree_honest_TOT_prune <- tree_prune
tree_honest_prune_list[[4]] <- tree_honest_TOT_prune
#fit

#honest fit

split.Rule.temp = "fit"
cv.option.temp = "matching"
split.Honest.temp = T
cv.Honest.temp = T
split.alpha.temp = .5

tree <- causalTree::honest.causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                          est_data=dataEst, est_treatment=dataEst$w,
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
#You can still prune as usual; the cptable is the one from training the tree
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_honest_fit <- tree
tree_honest_fit_prune <- tree_prune

tree_honest_prune_list[[5]] <- tree_honest_fit_prune

##################################
#Dishonest approaches


# CT
split.Rule.temp = "CT"
cv.option.temp = "matching"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = 1

tree <- causalTree::causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                   split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                   bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                   split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)

opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_dishonest_CT <- tree
tree_dishonest_CT_prune <- tree_prune
tree_dishonest_prune_list[[1]] <- tree_dishonest_CT_prune

#tstats
split.Rule.temp = "tstats"
cv.option.temp = "matching"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = 1

tree <- causalTree::causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                   split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                   bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                   split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)

opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_dishonest_tstats <- tree
tree_dishonest_tstats_prune <- tree_prune

tree_dishonest_prune_list[[2]] <- tree_dishonest_tstats_prune

#user
split.Rule.temp = "user"
cv.option.temp = "matching"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = 1

tree <- causalTree::causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)

opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_dishonest_user <- tree
tree_dishonest_user_prune <- tree_prune
tree_dishonest_prune_list[[3]] <- tree_dishonest_user_prune

#TOT
split.Rule.temp = "TOT"
cv.option.temp = "matching"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = 1

tree <- causalTree::causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                   split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                   bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                   split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
#You can still prune as usual; the cptable is the one from training the tree
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_dishonest_TOT <- tree
tree_dishonest_TOT_prune <- tree_prune
tree_dishonest_prune_list[[4]] <- tree_dishonest_TOT_prune

#fit
split.Rule.temp = "fit"
cv.option.temp = "matching"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = .5

tree <- causalTree::causalTree(as.formula(paste("y~",paste(f))), 
                   data=dataTrain, treatment=dataTrain$w, 
                   split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                   bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                   split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)

opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_dishonest_fit <- tree
tree_dishonest_fit_prune <- tree_prune
tree_dishonest_prune_list[[5]] <- tree_dishonest_fit_prune

# performance on test set
# if modifying for a real dataset, comment out the parts with infeasible MSE since we won't know that
MSEtau_infeas_honest <- array(0,5)
MSEtau_infeas_dishonest <- array(0,5)
MSEtau_TOT_honest <- array(0,5)
MSEtau_TOT_dishonest <- array(0,5)
MSEy_honest <- array(0,5)
MSEy_dishonest <- array(0,5)

# 1. alternative performance measure: MAPE (Mean Absolute Percentage Error)
MAPEtau_infeas_honest <- array(0,5)
MAPEtau_infeas_dishonest <- array(0,5)
MAPEtau_TOT_honest <- array(0,5)
MAPEtau_TOT_dishonest <- array(0,5)
MAPEy_honest <- array(0,5)
MAPEy_dishonest <- array(0,5)

# 2. alternative performance measure: MAE (Mean Absolute Error)
MAEtau_infeas_honest <- array(0,5)
MAEtau_infeas_dishonest <- array(0,5)
MAEtau_TOT_honest <- array(0,5)
MAEtau_TOT_dishonest <- array(0,5)
MAEy_honest <- array(0,5)
MAEy_dishonest <- array(0,5)

ystar = dataTest$y*(dataTest$w-propens)/(propens*(1-propens))
for (i in 1:5) {
  # use the predictions from honest trees on the Test set
  predicthonest = predict(tree_honest_prune_list[[i]],newdata=dataTest,type="vector")
  # use the predictions from dishonest trees on Test set
  predictdishonest = predict(tree_dishonest_prune_list[[i]],newdata=dataTest,type="vector")
  # this is simulated data, so we know the true tau for each individual.  This is infeasible MSE
  MSEtau_infeas_honest[[i]] = mean((predicthonest - dataTest$tau_true)^2, na.rm=T)
  MSEtau_infeas_dishonest[[i]] = mean((predictdishonest - dataTest$tau_true)^2, na.rm=T)
  # in a typical dataset, can use TOT loss function
  MSEtau_TOT_honest[[i]] = mean((predicthonest-ystar)^2, na.rm=T)
  MSEtau_TOT_dishonest[[i]] = mean((predictdishonest-ystar)^2, na.rm=T)
  ######### MAPE Criterion 
  MAPEtau_infeas_honest[[i]] = mean(abs((predicthonest - dataTest$tau_true)/dataTest$tau_true), na.rm=T)
  MAPEtau_infeas_dishonest[[i]] = mean(abs((predictdishonest - dataTest$tau_true)/dataTest$tau_true), na.rm=T)
  # in a typical dataset, can use TOT loss function
  MAPEtau_TOT_honest[[i]] = mean(abs((predicthonest - ystar)/ystar), na.rm=T)
  MAPEtau_TOT_dishonest[[i]] = mean(abs((predictdishonest - ystar)/ystar), na.rm=T)
  ########MAE
  MAEtau_infeas_honest[[i]] = mean(abs((predicthonest - dataTest$tau_true)), na.rm=T)
  MAEtau_infeas_dishonest[[i]] = mean(abs((predictdishonest - dataTest$tau_true)), na.rm=T)
  MAEtau_TOT_honest[[i]] = mean(abs((predicthonest - ystar)), na.rm=T)
  MAEtau_TOT_dishonest[[i]] = mean(abs((predictdishonest - ystar)), na.rm=T)
  
  # look at how the different partitions do predicting outcomes within a leaf
  dataTest$leaff <- as.factor(round(predicthonest,4))
  
  #honest estimation of predicted y
  dataEst$leaff <- as.factor(round(predict(tree_honest_prune_list[[i]], newdata=dataEst,type="vector"),4))
  yPredHonestTable <- melt(tapply(dataEst$y, list(dataEst$leaff, dataEst$w), mean), varnames=c("leaff","w"))
  yPredHonestTable <- rename(yPredHonestTable,replace=c("value"="ypredhon"))
  dataTest <- merge(dataTest,yPredHonestTable, by.x=c("leaff", "w"), by.y=c("leaff","w"))
  MSEy_honest[[i]] = mean((dataTest$ypredhon-dataTest$y)^2, na.rm=T)
  MAPEy_honest[[i]] = mean(abs((dataTest$ypredhon-dataTest$y)/dataTest$y), na.rm=T)
  MAEy_honest[[i]] = mean(abs((dataTest$ypredhon-dataTest$y)), na.rm=T)


  #dishonest estimation of predicted y
  dataTest$leaffd <- as.factor(round(predictdishonest,4))
  dataTrain$leaffd <- as.factor(round(predict(tree_dishonest_prune_list[[i]], newdata=dataTrain,type="vector"),4))
  yPredDishonestTable <- melt(tapply(dataTrain$y, list(dataTrain$leaffd, dataTrain$w), mean), varnames=c("leaffd","w"))
  yPredDishonestTable <- rename(yPredDishonestTable,replace=c("value"="ypreddishon"))
  dataTest <- merge(dataTest,yPredDishonestTable, by.x=c("leaffd", "w"), by.y=c("leaffd","w"))
  MSEy_dishonest[[i]] = mean((dataTest$ypreddishon-dataTest$y)^2, na.rm=T)
  MAPEy_dishonest[[i]] = mean(abs((dataTest$ypreddishon-dataTest$y)/dataTest$y), na.rm=T)
  MAEy_dishonest[[i]] = mean(abs((dataTest$ypreddishon-dataTest$y)), na.rm=T)
  
  dataTest <- dataTest[, !(names(dataTest) %in% c("ypredhon","ypreddishon"))]
}

print("infeasible MSE(tau)'s for honest CT, tstats, user, TOT, fit")
print(MSEtau_infeas_honest)
print("infeasible MSE(tau)'s for dishonest CT, tstats, user, TOT, fit")
print(MSEtau_infeas_dishonest)
print("TOT MSE(tau)'s for honest CT, tstats, user, TOT, fit")
print(MSEtau_TOT_honest)
print("TOT MSE(tau)'s for dishonest CT, tstats, user, TOT, fit")
print(MSEtau_TOT_dishonest)
print("MSE(y)'s for honest CT, tstats, user, TOT, fit")
print(MSEy_honest)
print("MSE(y)'s for dishonest CT, tstats, user, TOT, fit")
print(MSEy_dishonest)
############### MAPE
print("infeasible MAPE(tau)'s for honest CT, tstats, user, TOT, fit")
print(MAPEtau_infeas_honest)
print("infeasible MAPE(tau)'s for dishonest CT, tstats, user, TOT, fit")
print(MAPEtau_infeas_dishonest)
print("TOT MAPE(tau)'s for honest CT, tstats, user, TOT, fit")
print(MAPEtau_TOT_honest)
print("TOT MAPE(tau)'s for dishonest CT, tstats, user, TOT, fit")
print(MAPEtau_TOT_dishonest)
print("MAPE(y)'s for honest CT, tstats, user, TOT, fit")
print(MAPEy_honest)
print("MAPE(y)'s for dishonest CT, tstats, user, TOT, fit")
print(MAPEy_dishonest)

################

print("infeasible MAE(tau)'s for honest CT, tstats, user, TOT, fit")
print(MAEtau_infeas_honest)
print("infeasible MAE(tau)'s for dishonest CT, tstats, user, TOT, fit")
print(MAEtau_infeas_dishonest)
print("TOT MAE(tau)'s for honest CT, tstats, user, TOT, fit")
print(MAEtau_TOT_honest)
print("TOT MAE(tau)'s for dishonest CT, tstats, user, TOT, fit")
print(MAEtau_TOT_dishonest)
print("MAE(y)'s for honest CT, tstats, user, TOT, fit")
print(MAEy_honest)
print("MAE(y)'s for dishonest CT, tstats, user, TOT, fit")
print(MAEy_dishonest)

########################

# easy trick to get coefficients and standard errors

dataTrain$leaves <- predict(tree_dishonest_CT_prune, newdata=dataTrain, type = 'vector')
dataEst$leaves <- predict(tree_dishonest_CT_prune, newdata=dataEst, type = 'vector')
dataTest$leaves <- predict(tree_dishonest_CT_prune, newdata=dataTest, type = 'vector')


dataTrain$leavesf <- factor(round(dataTrain$leaves,4))
dataEst$leavesf <- factor(round(dataEst$leaves,4))
dataTest$leavesf <- factor(round(dataTest$leaves,4))

# run regressions with indicators for the leaves interacted with the treatment indicator
if (length(levels(dataTrain$leavesf)) == 1){
  
  modelTrain <- lm(y~w, data=dataTrain)
  modelEst <- lm(y~w, data=dataEst)
  modelTest <- lm(y~w, data=dataTest)
  
  summary(modelTrain)
  summary(modelEst)
  summary(modelTest)
  
} else{
  
  modelTrain <- lm(y~-1+leavesf+leavesf*w-w, data=dataTrain)
  modelEst <- lm(y~-1+leavesf+leavesf*w-w, data=dataEst)
  modelTest <- lm(y~-1+leavesf+leavesf*w-w, data=dataTest)
  
  print("Leaf names match estimated treatment effects on training set")
  print(summary(modelTrain))
  print("Estimated treatment effects on estimation set typically more moderate than training set")
  print(summary(modelEst))
  print("Estimated treatment effects on test set typically more moderate than training set")
  print(summary(modelTest))
  
  
  # extract the coefficient vectors which are the leaf treatment effects
  coefnumh <- length(coef(modelEst))
  coefnuml <- length(coef(modelEst))/2 + 1
  
  Train.coeftr <- coef(modelTrain)[coefnuml:coefnumh]
  Est.coeftr <- coef(modelEst)[coefnuml:coefnumh]
  Test.coeftr <- coef(modelTest)[coefnuml:coefnumh]
  
  # calculate leaf probabilities
 
  leafprobEst <- tapply(dataEst$y,list(dataEst$leavesf),length)
  leafprobTrain <- tapply(dataTrain$y,list(dataTrain$leavesf),length)
  leafprobTest <- tapply(dataTest$y,list(dataTest$leavesf),length)
  leafprob <- (leafprobEst + leafprobTrain + leafprobTest)/(nrow(dataEst) + nrow(dataTrain) + 
                                                              nrow(dataTest))
  
  #calculate variance of estimated treatment effects--typically this is higher in the training set, since there is overfitting there
  Train.coefvar <- sum(leafprob * Train.coeftr^2)-(sum(leafprob*Train.coeftr)^2)
  Est.coefvar <- sum(leafprob * Est.coeftr^2)-(sum(leafprob*Est.coeftr)^2)
  Test.coefvar <- sum(leafprob * Test.coeftr^2)-(sum(leafprob*Test.coeftr)^2)
  
  print("Variance of estimated treatment effects: Train, Estimation, Test Sets")
  print("Typically train has higher var--more extreme estimates--due to overfitting")
  print(c(Train.coefvar,Est.coefvar,Test.coefvar))
}      

