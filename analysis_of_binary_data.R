setwd('C:/Users/celib/Desktop/DENMARK/DTU/1-SEMESTER/Statistical Modelling/projects/survival_data')

library(numDeriv)
library(matlib)

### ------------------------------------ ###
### Part 1 - Analysis of the binary data ###
### ------------------------------------ ###

## See 4.3. in the book

# Read data
logistic <- read.table('logistic.txt', header = TRUE, sep = '\t')


# Looks like:
# ATZ     AIDS_yes      n
# Yes     25            170
# No      44            168


# For the analysis, think:
#           AZT           control           total
# AIDS      x=25          y=44              t=69
# no AIDS   m-x=145       n-y=124           N-t=269
# total     m=170         n=168             N=338


### ------------------------------------------------ ###
### ------------------------------------------------ ###

## Fit the Binomial distribution to the data (consider all data as coming from the same population):
# --> success = individual develops AIDS

# N (number of observations) = total number of individuals (338)
N <- logistic$n[1] + logistic$n[2]    
# t (number of successes) = total number of AIDS cases (69)
t <- logistic$AIDS_yes[1] + logistic$AIDS_yes[2]    


#observed proportion of AIDS in the population (0.2041)
obs_prop <- t/N     


### ------------------------------------------------ ###
### ------------------------------------------------ ###

## Fit the Binomial separately to the two distributions:

# Group 1: AZT treated
# m (number of observations) = number of individuals treated with AZT (170)
m <- logistic$n[1]
# x (number of successes) = number of AIDS cases, of those treated with AZT (25)
x <- logistic$AIDS_yes[1]

#observed proportion of AIDS in AZT-treated group (0.1471)
obs_prop_AZT <- x/m


# Group 2: Control
# n (number of observations) = number of control individuals (168)
n <- logistic$n[2]
# y (number of successes) = number of AIDS cases, of the control (44)
y <- logistic$AIDS_yes[2]

#observed proportion of AIDS in Group 2 (0.2619)
obs_prop_ct <- y/n



### ------------------------------------------------ ###
### ------------------------------------------------ ###

# To visualize the data as a random binomial distribution with the given parameters:


# Full population
binom <- dbinom(x = 0:N, prob = obs_prop, size = N)
binom_norm <- binom/max(binom)
plot(x = 0:N, y = binom_norm, type = "h", col = "purple", xlab = "Individuals with AIDS",
     ylab = "Probability", main = "Probability mass function")

# AZT-treated group
binom_AZT <- dbinom(x = 0:N, prob = obs_prop_AZT, size = N)
binom_AZT_norm <- binom_AZT/max(binom_AZT)
lines(x = 0:N, y = binom_AZT_norm, type = "h", col = "violet")

# Control group
binom_ct <- dbinom(x = 0:N, prob = obs_prop_ct, size = N)
binom_ct_norm <- binom_ct/max(binom_ct)
lines(x = 0:N, y = binom_ct_norm, type = "h", col = "lightblue")

# Add a legend
legend("topright", legend=c("Full population", "AZT-treated", "Control"),
       fill=c("purple", "violet", "lightblue"), cex=0.8)


### ------------------------------------------------ ###
### ------------------------------------------------ ###

## Do these groups have equal success probabilities?

# Chi squared test on proportions
prop.test(c(x,y), c(m, n))

# with a p-val of 0.01, we can reject the null hypothesis that the proportions are equal


### ------------------------------------------------ ###
### ------------------------------------------------ ###

## Estimate parameters in the model (p0 probability of AIDS in control group, 
## p1 probability of AIDS in treatment group) and report a confidence interval 
## for the parameter describing the difference (beta1), compare with the result above.

# We will do this with a likelihood analysis

# We consider two distinct binomial distributions, one for the AZT-treated group and 
# another one for the control
#
# X ~ B(obs_prop_ct, n)
# Y ~ B(obs_prop_AZT, m)
#
# In order to build one single model that takes into consideration both distributions, 
# we parametrize: find a parameter that explains the relationship between obs_prop_ct and
# obs_prop_AZT (beta1) and put all the other variation in a nuissance parameter (beta0). 
# With this, build a model.
#
# We choose the log-odds ratio as our theta 
# beta1 = (obs_prop_AZT/(1-obs_prop_AZT)) / (obs_prop_ct/(1-obs_prop_ct))
#
# In small samples, the likelihood of the log-odds ratio is more regular than the 
# likelihood of other parameters, hence our choice.


## Negative log-Likelihood function (see section 4.3. of the book: comparing two proportions)

# My optimizer only takes one parameter (or a vector of parameters), therefore I use the variable
# params, which is a vector of beta1 and beta0, in that order(params[1] is beta1 and params[2] is beta0)

# We have to make the likelihood negative because nlminb finds the minimum and we want the maximum

jllikfun <- function(params){
  a <- exp(params[1]*x)
  b <- exp(params[2]*(x+y))
  c <- (1+exp(params[1]+params[2]))**(-m)
  d <- (1+exp(params[2]))**(-n)
  return(-log(a*b*c*d))
  }

# Optimization
fitted <- nlminb(jllikfun, start=c(-0.5, -1.0))   
beta1_hat <- fitted$par[1]    # -0.7218
beta0_hat <- fitted$par[2]      # -1.0361
llik <- fitted$objective


print("Parameters (beta1 and beta0)")
beta1_hat
beta0_hat


## We can calculate confidence intervals two ways: 

## Wald confidence intervals (quadratic approximation - from Information matrix) 
# s.e (beta1) = sqrt(I(beta1)) and s.e (beta0) = sqrt(I(beta0)), and I = H^-1
hess <- hessian(jllikfun, c(beta1_hat,beta0_hat))
I <- inv(hess)

se_beta1 = sqrt(I[1,1]) # 0.2787263
se_beta0 = sqrt(I[2,2]) # 0.1754759

CI_beta1 <- beta1_hat + c(-1,1) * qnorm(0.975) * se_beta1 # qnorm gives the quantile function and incorporates the confidence level in the interval
CI_beta0 <- beta0_hat + c(-1,1) * qnorm(0.975) * se_beta0

print("Wald confidence intervals (beta1 and beta0)")
CI_beta1
CI_beta0



## Likelihood-based confidence intervals (plotting the likelihood)
# Graphically: plotting the profile likelihood

# I'm defining the likelihood function instead of using the log-likelihood that I defined before 
# but it would be okay to look at the log-likelihood as well. Just take the exponent out of 
# exp(-0.5 * qchisq(0.95, df = 1)) when plotting the confidence threshold in the profile likelihood

jlikfun <- function(beta1,beta0){            
  a <- exp(beta1*x)
  b <- exp(beta0*(x+y))
  c <- (1+exp(beta1+beta0))**(-m)
  d <- (1+exp(beta0))**(-n)
  return(a*b*c*d)
}

# To plot the profile likelihood we need to evaluate the likelihood function at a range of values 
# of the parameter we are interested in (beta1) and fix the other one (beta0 = beta0_hat = MLE = -1.036)
beta1 <- seq(-3,2, by = 0.01)
proflik = sapply(beta1, jlikfun, beta0=beta0_hat)
matplot(beta1, proflik/max(proflik), type="l", ylab="Profile Likelihood (beta0=MLE)")
lines(range(beta1), exp(-0.5 * qchisq(0.95, df = 1)) * c(1,1), col=2, lty=2)


# Numerically: finding the intersections with uniroot (I can't get this to work)



## INTERPRETATION
# Parameter beta1 gives us info about the log ratio between both proportions 
# (deaths/AIDS cases in AZT-treated and control groups)
# 
# If the CI of beta1 includes 0, then we cannot reject the null hypothesis that both 
# proportions are the same (log1=0)
# 
# In this case, we can reject it.



#--------------------------------#
#### ASSIGNMENT 2
#--------------------------------#
## Fit a logistic regression model for the binary outcome AIDS="yes" 
## versus AIDS="no" with the explanatory variable treatment with AZT (Yes, No). 
## Present the odds ratio for the effect of AZT on AIDS with 95% confidence interval 
## and interpret the result in words.

## I need my data to be a dataframe with 2 columns: 
## one with values 0/1 depending on treatment, the other with 0/1 depending on AIDS/noAIDS

treatment = c(rep(1,170), rep(0,168))
outcome = c(rep(1,25), rep(0,145), rep(1,44), rep(0,124))

data = data.frame(treatment,outcome)

y <- data[ ,2]
X <- cbind(1,data[ ,1])

fit = glm(y~-1+X,family=binomial)   # means "regress x on y, but leave out the intercept"
summ = summary(fit)

# We find the parameter estimates and SEs in the glm and summary objects
# Parameter estimates
X1_hat = fit$coefficients[1]    # X1 is beta0 (nuissance, same as beta0 before)
X2_hat = fit$coefficients[2]    # X2 is beta1 (treatment effect, same as beta1 before)

# SEs
X1_se = summ$coefficients[1,2]
X2_se = summ$coefficients[2,2]

print('Parameter estimates:')
print(X1_hat)
print(X2_hat)

print('S.E.s of parameter estimates:')
print(X1_se)
print(X2_se)


## ODDS RATIO (and its 95% CI) of the effect of AZT on AIDS

# The odds ratio (OR) is a measure of how strongly an event is associated 
# with exposure. The odds ratio is a ratio of two sets of odds: the odds of 
# the event occurring in an exposed group versus the odds of the event occurring 
# in a non-exposed group.

# Since we assigned y to something else, reassign to number of successes
x = 25
y = 44
m = 170
n = 168

# Odds in the AZT-treated group = AZT-treated with AIDS/AZT-treated without AIDS
O_AZT = x/(m-x)
# Odds in the control group = control with AIDS/control without AIDS
O_ct = y/(m-y)
# Odds ratio = odds in the AZT-treated group/odds in the control group - 0.49
OR = O_AZT/O_ct

# CI - (0.2860442, 0.8522099)
CI_OR <- exp(log(OR) + c(-1,1) * qnorm(0.975) * sqrt(1/x+1/(m-x)+1/y+1/(m-y)))


## Interpretation
# If the confidence interval for the odds ratio includes the number 1 then the calculated odds ratio 
# would not be considered statistically significant. This can be seen from the interpretation of the odds ratio. 

# An odds ratio greater than 1 implies there are greater odds of the event happening in the exposed versus 
# the non-exposed group. 
# An odds ratio of less than 1 implies the odds of the event happening in the exposed group 
# are less than in the non-exposed group. 
# An odds ratio of exactly 1 means the odds of the event happening are the exact same in the exposed 
# versus the non-exposed group. 

# Thus, if the confidence interval includes 1 (eg, [0.01, 2], [0.99, 1.01], or [0.99, 100]), 
# then the expected true population odds ratio may be above or below 1, so it is uncertain whether 
# the exposure increases or decreases the odds of the event happening with our specified level of confidence.

# In this case, the interval does not contain 1, thus the AZT has a true effect on the development of the disease.
# More specifically, the chances of AZT-treated HIV patients developing AIDS is 28-84% of those in the control. 




#----------------------------------------------------------#
## Test the hypothesis of no effect of AZT on AIDS using:
##    - The likelihood ratio test
##    - The Wald test
##    - The score test (lecture 9)
#----------------------------------------------------------#

## LIKELIHOOD RATIO TEST
# H0: beta1 = 0
# H1: beta1 != 0

# We have to optimize the likelihood again when beta1 is fixed to 0, in order to 
# find the likelihood of H0 

jllikfun0 <- function(beta0){
  a <- exp(0*x)
  b <- exp(beta0*(x+y))
  c <- (1+exp(0+beta0))**(-m)
  d <- (1+exp(beta0))**(-n)
  return(-log(a*b*c*d))
}

fitted0 <- nlminb(jllikfun0, start=c(-1.0))
beta0_hat_0 <- fitted0$par
llik0 <- fitted0$objective

# Does -2 * log(likelihood of H0/likelihood of H1) follow a chi2 distribution?
# Equivalently: 2 * (log-likelihood_H0(beta1=0)) - log-likelihood_H1(beta1=beta1_hat))
# llik0 is the log-likelihood of the null hypothesis
# llik is the log-likelihood of the alternative hypothesis

Q <- 2 * (llik0-llik) ## test statistic

#print('p-value for the LRT')
1 - pchisq(Q, df = 1) # p-value 0.00848

#----------------------------------------------------------#

## WALD TEST
W <- beta1_hat**2/se_beta1**2
1-pchisq(W, df=1) # p-value 0.009611126

#----------------------------------------------------------#

## SCORE TEST (from observed information)
# We have to evaluate the likelihood with 2 parameters, where beta1=0 
# (but beta0 is optimized in the one-parameter likelihood function)

# Score vector at H0
S <- jacobian(jllikfun, c(0,beta0_hat_0))
# Information matrix at H0
I <- hessian(jllikfun, c(0,beta0_hat_0))

z.obs <- S %*% solve(I) %*% t(S)
p.score.obs <- 1-pchisq(z.obs,df=1) ## p-value 0.008816727

#----------------------------------------------------------#

