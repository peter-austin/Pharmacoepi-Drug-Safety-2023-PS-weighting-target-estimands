# Computes the true estimand for ATE/ATT/MW/OW in the different scenarios.
# Assumes an interaction between treatment and covariates in the outcome model.
# Assume a continuous outcome. There are two continuous linear predictors:
# 1) one for treatment-selection model.
# 2) one for outcome model.

# NOTE: THIS CODE IS FOR ILLUSTRATIVE PURPOSES ONLY AND COMES WITH
#       ABSOLUTELY NO WARRANTY.
library(mvtnorm)
library(rms)

treat.caliper <- 200
# Calipers for grid search for parameters.

N.pop <- 1000000
# Size of super-population.

for (diff.means in 1:1){
# Regression coefficient for treatment effect on continuous outcome.

for (prop.treat in seq(0.1,0.9,0.1)){
# Prevalence of treatment.

for (cstat in c(0.65,0.75,0.85)){
# c-statistic of propensity score model.

for (LP.corr in c(0.2,0.4,0.6,0.8)){
# Correlation between linear predictors for treatment and outcome.
# Don't allow zero colleration, since this would induce zero confounding.

for (beta.interact in c(0,0.2,0.4,0.6,0.8,1)){
# Coefficient for interaction between LP.outcome and treatment status.

  set.seed(25112021)

################################################################################
# Generate two continuous linear predictors (covariates). One for 
# treatment-selection model, one for outcome model.
# We assume that they are correlated.
################################################################################

x <- rmvnorm(N.pop,mean=c(0,0),sigma=matrix(c(1,LP.corr,LP.corr,1),nrow=2,
  ncol=2))

LP.treat <- x[,1]
LP.outcome <- x[,2]

remove(x)

X.treat <- cbind(1,LP.treat)
X.outcome <- cbind(1,LP.outcome)

################################################################################
# Generate a binary treatment variable with the given prevalence of treatment.
################################################################################

B.treat <- qnorm(cstat)*sqrt(2)
# Initial estimate of regression coefficient using BMC MRM paper.

cstat.emp <- 0

iter.cstat <- 0

scalar <- 1
B.treat.modified <- B.treat

while (abs(cstat.emp - cstat) > 0.005){

  iter.cstat <- iter.cstat + 1

  B.treat.modified <- scalar * B.treat.modified
  # Modified coefficients to increase/decrease c-statistic of PS model.

  treat.function <- function(b0.treat){
  # Function for determining prevalence of treatment with a given intercept.

    beta.treat.modified <- c(b0.treat,B.treat.modified)
    # Set the intercept of the treat model to the given value.

    XB <- X.treat %*% beta.treat.modified

    p.treat <- exp(XB)/(1 + exp(XB))
    Y <- rbinom(N.pop,1,p.treat)

    return(list(Y,beta.treat.modified))

    remove(Y,beta.treat.modified,XB,p.treat)

  }
  # End of function for estimating prevalence of treatment.

# Use a bisection approach to determine the intercept that results in the
# desired prevalence of treatment.
int.low <- -treat.caliper
int.high <- treat.caliper

iter.treat <- 0
treat.prev <- 1

while(abs(treat.prev - prop.treat) > 0.0001){
  iter.treat <- iter.treat + 1
  set.seed(iter.treat)

  int.mid <- (int.low + int.high)/2

  treat.fun.results <- treat.function(b0.treat=int.mid)
  treat <- treat.fun.results[[1]]
  treat.prev <- mean(treat)

  beta.treat.final <- treat.fun.results[[2]]
  remove(treat.fun.results)

  if (treat.prev < prop.treat) int.low <- int.mid else
    int.high <- int.mid
  
  cat(iter.treat,int.low,int.mid,int.high,treat.prev,prop.treat,
    file="iter.treat.out",fill=T,append=T)
}

beta0.treat <- int.mid
# Intercept for treatment-selection model.

remove(cstat.emp,int.low,int.mid,int.high,treat.prev,treat.function)

################################################################################
# Estimate propensity score model in super-population
################################################################################

psm.pop <- lrm(treat ~ LP.treat)

cstat.emp <- psm.pop$stats["C"]

remove(psm.pop)

if (cstat.emp < cstat) scalar <- 0.95
   scalar <- 1.05

cat(iter.cstat,iter.treat,B.treat.modified,
  cstat,cstat.emp,prop.treat,mean(treat),
  file="estimand.iter.out",fill=T,append=T)

}

remove(iter.treat,scalar,iter.cstat,B.treat.modified)

################################################################################
# Generate a continuous outcome
################################################################################

R2 <- 0.25
# Target R2 for outcomes model

XB <- 0 + 1*LP.outcome

sigma2 <- (var(XB) - R2*var(XB))/R2

Ycontinuous0 <- 0 + 1*LP.outcome + rnorm(N.pop,0,sqrt(sigma2)) 

Ycontinuous1 <- 0 + 1*LP.outcome + diff.means*1 + beta.interact*LP.outcome +
  rnorm(N.pop,0,sqrt(sigma2))
# Allow for interaction between LP and treatment.

Y.continuous <- treat*Ycontinuous1 + (1-treat)*Ycontinuous0

remove(R2,XB,sigma2,LP.treat,LP.outcome)

##############################################################################
# Create true propensity score function. 
# Compute the true effect size for MW and OW.
##############################################################################

XB <- X.treat %*% beta.treat.final

ps <- exp(XB)/(1 + exp(XB))

remove(B.treat,beta0.treat,beta.treat.final,XB)

u <- pmin(ps,1-ps)
entropy <- -ps*log(ps) - (1-ps)*log(1-ps)

wt.ate <- (treat/ps) + (1-treat)/(1-ps)
wt.att <- treat + (1-treat)*ps/(1-ps)
wt.mw <- ((treat*u)/ps) + ((1-treat)*u)/(1-ps)
wt.ow <- treat*(1-ps) + (1-treat)*ps
wt.ew <- treat*entropy/ps + (1-treat)*entropy/(1-ps)

remove(u,entropy,ps)

Y.continuous.mean.0.ate <- weighted.mean(Y.continuous[treat==0],wt.ate[treat==0])
Y.continuous.mean.1.ate <- weighted.mean(Y.continuous[treat==1],wt.ate[treat==1])
Y.continuous.mean.0.att <- weighted.mean(Y.continuous[treat==0],wt.att[treat==0])
Y.continuous.mean.1.att <- weighted.mean(Y.continuous[treat==1],wt.att[treat==1])
Y.continuous.mean.0.mw <- weighted.mean(Y.continuous[treat==0],wt.mw[treat==0])
Y.continuous.mean.1.mw <- weighted.mean(Y.continuous[treat==1],wt.mw[treat==1])
Y.continuous.mean.0.ow <- weighted.mean(Y.continuous[treat==0],wt.ow[treat==0])
Y.continuous.mean.1.ow <- weighted.mean(Y.continuous[treat==1],wt.ow[treat==1])
Y.continuous.mean.0.ew <- weighted.mean(Y.continuous[treat==0],wt.ew[treat==0])
Y.continuous.mean.1.ew <- weighted.mean(Y.continuous[treat==1],wt.ew[treat==1])

diff.means.ate <- Y.continuous.mean.1.ate - Y.continuous.mean.0.ate
diff.means.att <- Y.continuous.mean.1.att - Y.continuous.mean.0.att
diff.means.mw <- Y.continuous.mean.1.mw - Y.continuous.mean.0.mw
diff.means.ow <- Y.continuous.mean.1.ow - Y.continuous.mean.0.ow
diff.means.ew <- Y.continuous.mean.1.ew - Y.continuous.mean.0.ew

remove(Y.continuous.mean.0.ate,Y.continuous.mean.1.ate,
  Y.continuous.mean.0.att,Y.continuous.mean.1.att,
  Y.continuous.mean.0.mw,Y.continuous.mean.1.mw,
  Y.continuous.mean.0.ow,Y.continuous.mean.1.ow,
  Y.continuous.mean.0.ew,Y.continuous.mean.1.ew,
  wt.ate,wt.att,wt.mw,wt.ow,wt.ew)

cat(diff.means,prop.treat,cstat,LP.corr,beta.interact,diff.means.ate,
  "DiffMeans","ATE",file="estimand.out",fill=T,append=T)
cat(diff.means,prop.treat,cstat,LP.corr,beta.interact,diff.means.att,
  "DiffMeans","ATT",file="estimand.out",fill=T,append=T)
cat(diff.means,prop.treat,cstat,LP.corr,beta.interact,diff.means.mw,
  "DiffMeans","MW",file="estimand.out",fill=T,append=T)
cat(diff.means,prop.treat,cstat,LP.corr,beta.interact,diff.means.ow,
  "DiffMeans","OW",file="estimand.out",fill=T,append=T)
cat(diff.means,prop.treat,cstat,LP.corr,beta.interact,diff.means.ew,
  "DiffMeans","EW",file="estimand.out",fill=T,append=T)

remove(diff.means.ate,diff.means.att,diff.means.mw,diff.means.ow,
  diff.means.ew,X.outcome,X.treat,Y.continuous,Ycontinuous0,
  Ycontinuous1,treat)

}
}
}
}
}