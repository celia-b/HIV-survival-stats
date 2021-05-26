setwd('C:/Users/celib/Desktop/DENMARK/DTU/1-SEMESTER/Statistical Modelling/projects/survival_data')

library(numDeriv)
library(matlib)
library(dplyr)
library(survival)
library(stats)

### ------------------------------------------- ###
### Part 2 - Analysis of the survival time data ###
### ------------------------------------------- ###


actg <- read.table('actg320.txt', header = TRUE, sep = '\t')


### ------------------------------------------------ ###
### ------------------------------------------------ ###

## How many patients got AIDS or died in the two treatment groups?
## What is the proportion of patients that got AIDS or died in the two
## group? Other relevant number that could be calculated?

# Data summary (number of patients by treatment and by event)
tr_event_count <- actg %>% group_by(tx) %>% count(event)

# Total number of patients in the study
N <- tr_event_count$n[1]+tr_event_count$n[2]+tr_event_count$n[3]+tr_event_count$n[4]

# Patients on treatment 0 (old)
t0 <- tr_event_count$n[1]+tr_event_count$n[2]
# Patients on treatment 1 (new)
t1 <- tr_event_count$n[3]+tr_event_count$n[4]

# Total patients who leave the study/the study ends and they haven't developed AIDS (event 0)
e0 <- tr_event_count$n[1]+tr_event_count$n[3]
# Total patients who die/develop AIDS (event 1)
e1 <- tr_event_count$n[2]+tr_event_count$n[4]

# Number of patients w/AIDS or dead (event 1) in each treatment
t0e1 <- tr_event_count$n[2]   #old treatment --> 63 deaths/AIDS cases
t1e1 <- tr_event_count$n[4]   #new treatment --> 33 deaths/AIDS cases


# Proportion of patients w/AIDS or dead in full population
prop <- e1/N      # 0.0834

# Proportion of patients w/AIDS or dead in treatment 0 (old)
prop0 <- t0e1/t0  # 0.1092

# Proportion of patients w/AIDS or dead in treatment 1 (new)
prop1 <- t1e1/t1  # 0.0575


### ------------------------------------------------ ###
### ------------------------------------------------ ###

## Fit an exponential distribution, using numerical methods, to the time
## of event (time) in the data set, remember to take into account that
## some of the data is censored (i.e. we only know that the time to the
## event is longer that the reported time).

# To deal with the censor data we divide our data based on the event. 
# For data with event=1, we use the probability density function, 
# as we know exactly the time of diagnosis/death.
# For data with event=0, we use the cumulative probability function, 
# as we only know that the diagnosis/death will happen after the given time. 

# Getting my data to look like rat example (lecture 6)
timevector <- dplyr::pull(actg, time)
eventvector <- dplyr::pull(actg, event)
treatmentvector <- dplyr::pull(actg, tx)

data <- list(time=as.double(timevector), event=as.double(eventvector), treatment=as.double(treatmentvector))


# A quick look at data (ignoring censoring)
# ecdf is "Empirical Cumulative Distribution Function"
plot(ecdf(data$time[data$treatment==0]),xlim = range(actg$time), main = "Empirical Cumulative Distribution Function")  # control treatment
lines(ecdf(data$time[data$treatment==1]),col=2)                 # new treatment



#-----------------------------------------------------------------#
## likelihood estimation: Exponential model, disregarding treatment
#-----------------------------------------------------------------#

# I'm looking at theta (average number of days until event) instead of lambda 
# (rate of event happening) because I like it better, but changing rate=1/theta 
# to rate=lambda we can optimize lambda instead

ll_exp <- function(theta, delta, data){
  - sum(dexp(data[delta==1], rate=1/theta, log=TRUE)) -       # non-censored data = probability density function (dexp)
    sum(log(1-pexp(data[delta==0],rate=1/theta)))             # censored data = cumulative probability function (pexp; in this case 1-pexp because the diagnosis/death will happen after the recorded time)
}

## Eg. estimating theta for control treatment from global model:
print('Control, based on global model ($minimum = mean number of days to AIDS/death):')
optimize(ll_exp, c(0,3000),
         delta = data$event[data$treatment==0],
         data = data$time[data$treatment==0])

## Eg. estimating theta for new treatment from global model:
print('New treatment, based on global model ($minimum = mean number of days to AIDS/death):')
optimize(ll_exp,c(0,5000),
         delta = data$event[data$treatment==1],
         data = data$time[data$treatment==1])

## Common mean for both treatments
print('Common mean for both treatments ($minimum = mean number of days to AIDS/death):')
(OPT0 <- optimize(ll_exp,c(0,3000),
                  delta = data$event,
                  data = data$time))

#------------------------------------------------------------------------#
## likelihood estimation: Exponential model, taking treatment into account
#------------------------------------------------------------------------#

ll_exp2 <- function(theta, delta, data, treatment){
  ll_exp(theta[1], delta[treatment==0], data[treatment==0]) +
    ll_exp(theta[2], delta[treatment==1], data[treatment==1])
}

print('$par[1] = Control, $par[2] = New treatment')
(OPT1 <- nlminb(c(2000,4000),ll_exp2,
                delta = data$event,
                data = data$time, 
                treatment = data$treatment))

#------------------------------------------------------------------------#
## likelihood comparison: all data vs. separated by treatment
#------------------------------------------------------------------------#

## Parameters
# All data (theta: mean number of days to AIDS/death):
OPT0$minimum 
# All data (lambda: rate of AIDS development/death):
1/OPT0$minimum 

# Groups separated (theta: mean number of days to AIDS/death):
OPT1$par 
# Groups separated (lambda: rate of AIDS development/death):
1/OPT1$par #groups separated

## Negative log-likelihood for each model
OPT0$objective #groups together
OPT1$objective # groups separated

# AICs
2 * OPT0$objective + 2*1 #groups together
2 * OPT1$objective + 2*2 # groups separated - not sure if 2*2 or 2*1, I think 2*2 bc 2 parameters

## AIC for the second model seems slightly better, but is it significant?
## Do Likelihood ratio test
chi2 <- - 2 * ((OPT1$objective)-OPT0$objective) ## test statistics
1 - pchisq(chi2, df = 1) ## p-value = 0.000828
## Hence difference!!! The model taking into account groups is better


## -----------------------------------------------------------------------------------#
## likelihood estimation: Exponential model, one parameter indicates treatment effect
#-------------------------------------------------------------------------------------#

## Formulate a model where one parameter indicate the treatment effect, 
## find the MLE and compare with the result above. You should also find the 
## Wald confidence interval for the treatment parameter


# theta[1] is beta0 and theta[2] is beta1
ll_exp3 <- function(theta, delta, data, treatment){
  ll_exp(theta[1], delta[treatment==0], data[treatment==0]) +
    ll_exp(theta[1]+theta[2], delta[treatment==1], data[treatment==1]) #theta[2] is the difference between the groups
}

# $par[1] = beta0(nuissance), $par[2] = beta1(difference, what the treatment adds to the baseline)
(OPT2 <- nlminb(c(250,0),ll_exp3,
                delta = data$event,
                data = data$time, 
                treatment = data$treatment))


#------------------------------------------------------------------------------------------------------#
## likelihood comparison: separated by treatment vs. single parameters estimates difference
#------------------------------------------------------------------------------------------------------#
## Parameters
# Separated by treatment (mean)
OPT1$par
# Separated by treatment (rate)
1/OPT1$par
# Single parameter estimates difference (mean)
OPT2$par
# Single parameter estimates difference (rate)
1/OPT2$par

## Negative log-likelihood
# Separated by treatment
OPT1$objective
# Single parameter estimates difference
OPT2$objective

# AICs
2 * OPT1$objective + 2*2 # by groups
2 * OPT2$objective + 2*2 # one parameter estimates difference

## Is the new model better than the last?
## Do Likelihood ratio test
chi2 <- - 2 * ((OPT2$objective)-OPT1$objective) ## test statistics
print('p-value for the LRT')
1 - pchisq(chi2, df = 1) ## p-value = 0.884
## Hence no difference (obviously, because it is essentially the same model but formulated in a different way)


#-----------------------------------------------#
## Parameter uncertainty
#-----------------------------------------------#
# Find the Wald confidence interval for the treatment parameter in the model above

H <- hessian(ll_exp3, OPT2$par, delta = data$event,
             data = data$time, treatment = data$treatment)
# S.E. for beta0 and beta1
sqrt(diag(solve(H)))

## CI for the difference between groups (theta2)
CI.diff <- OPT2$par[2] + c(-1,1) * qnorm(0.975) * sqrt(diag(solve(H)))[2] #qnorm gives the quantile function
# Wald CI for beta1 (the treatment effect parameter)
CI.diff



## Derive the theoretical results for the models above, including the standard error 
## estimates, use this to formulate and implement the profile likelihood function 
## for the treatment parameter.

### Likelihood-based CI
## Graphically, by profiling the parameter (Here I show the log-likelihood, when before I showed the likelihood. Should probably make it consistent, although results are the same)
pll <- function(theta2, delta, data, treatment){
  f.tmp <- function(theta1, theta2, delta, data, treatment){
    ll_exp3(c(theta1,theta2), delta, data, treatment)
  }
  optimize(f.tmp,c(0,4000),theta2 = theta2, delta = delta,
           data = data, treatment = treatment)$objective
}

theta2 <- seq(100,5000)
# This takes a while to plot because the interval is very big
matplot(theta2, -(sapply(theta2,pll,delta=data$event, 
                         data = data$time, treatment = data$treatment) -
                    OPT1$objective),type="l",ylab="Profile log-likelihood (beta0=MLE)")
lines(range(theta2), - c(1,1)*qchisq(0.95,df=1)/2,col=2,lty=2)
rug(CI.diff)


## Numerically: intersect between profile likelihood and red line with uniroot() --> see week1 exercises


#-----------------------------------------------------#
#-----------------------------------------------------#
#### ASSIGNMENT 2 ####
#-----------------------------------------------------#
#-----------------------------------------------------#

## Descriptive statistics (see above)

# Follow-up time

# Control treatment:
median(actg$time[actg$tx==0])
# New treatment:
median(actg$time[actg$tx==1])

#-----------------------------------------------------#

## Plot the survival functions in the two treatment groups.
## Which group seems to be doing best?
surv_by_tx = survfit(Surv(time, event == 1) ~ tx, conf.type="log-log", data=actg)
surv_by_tx  # NA values?
plot(surv_by_tx, las=1, xlab='Days since admission', ylab='Estimated survival probability',
     main='Survival functions', col=c('black', 'red'), lwd=2)
legend("bottomleft", col=c('black', 'red'), c('Control', 'New Treatment'), lwd=2)

# The group with the new treatment has higher estimated survival probability 
# throughout the study, especially as time passes

#-----------------------------------------------------#

## Plot the cumulative incidence functions for the two groups.
## Which plot would you prefer?
plot(surv_by_tx, fun=function(x) {1- x}, conf.int=F, las=1, ylim=c(0,0.15), xlab="Days since admission",
     ylab="Estimated Failure Probability", main='Cumulative incidence functions', col=c('black', 'red'), lwd=2)
legend("bottomright", col=c('black', 'red'), c("Control","New Treatment"), lwd=2)

# This plot is more readable than the other one, so we prefer this one.

#-----------------------------------------------------#

## Compare the survival in the two treatment groups using a log-rank test.
survdiff(Surv(time, event == 1) ~ tx, data=actg)
# With a p-val = 0.001, we can reject the null hypothesis that the two curves are the same

#-----------------------------------------------------#
## Parametric survival models 
## Fit parametric survival models containing treatment (tx) and CD4 count (cd4) 
## as explanatory variables. 
## Try using the exponential, Weibull and log-logistic models, 
## which one gave the best fit (and why)?

# Exponential model
mod_exp <- survreg(Surv(time, event) ~ tx + cd4, data=actg, dist="exponential")
summary(mod_exp)

# Weibull model
mod_weib <- survreg(Surv(time, event) ~ tx + cd4, data=actg, dist="weibull")
summary(mod_weib)

# Log-logistic model
mod_loglog <- survreg(Surv(time, event) ~ tx + cd4, data=actg, dist="loglogistic")
summary(mod_loglog)

# Comparing models
AIC_exp <- -2*mod_exp$loglik[2]+ 2*mod_exp$df; AIC_exp
AIC_weib <- -2*mod_weib$loglik[2]+ 2*mod_weib$df; AIC_weib
AIC_loglog <- -2*mod_loglog$loglik[2]+ 2*mod_loglog$df; AIC_loglog

# The log-logistic model has the lowest AIC and is therefore chosen.
# Is there a biological sense to this?

#-----------------------------------------------------#

## Using the survival model you chose, make a table of estimates and their 
## 95% confidence intervals.
confint(mod_loglog)

## Using your model compute the time ratio for the treatment effect. 
## Similarly, compute the time ratio for the effect of increasing the CD4 count 
## with 50. In both cases unceartainty evaluation (e.g. confidence intervals) 
## should be included. Interpret the results in words.
round(exp(cbind(TR = coef(mod_loglog), confint(mod_loglog))), 2)
round(exp(cbind(TR = coef(mod_loglog)*50, confint(mod_loglog)*50)), 2)


## INTERPRETATION

# The time ratio is 2.32 and the confidence interval is [1.41, 4.10]. 
# Receiving the treatment increases the median survival time by 2.32 times.

# Increasing the CD4 count by 50 increases the median survival time with a factor of 2.83.


#-----------------------------------------------------#

## Assess the goodness of fit of this model using a plot based on the Cox Snell residuals.
# This is wrong - the scale parameter should somehow be included
actg$CoxSnell <- actg$time*exp(-mod_loglog$linear.predictors)
surv <- survfit(Surv(CoxSnell, event==1)~1 , data=actg)
plot(surv$time, -log(surv$surv), main='Cox Snell residuals plot')
abline(a=0, b=1)

#-----------------------------------------------------#

## Give a graphical presentation of your model

xrange <- range(actg$time)
t <- seq(xrange[1],xrange[2],length=500)

coef <- mod_loglog$coefficients
z1 <- (log(t)-(coef[1]+coef[3]*43))/mod_loglog$scale # cd4=43, tx=0
z2 <- (log(t)-(coef[1]+coef[3]*86))/mod_loglog$scale # cd4=86, tx=0
z3 <- (log(t)-(coef[1]+coef[2]+coef[3]*43))/mod_loglog$scale # cd4=43, tx=1
z4 <- (log(t)-(coef[1]+coef[2]+coef[3]*86))/mod_loglog$scale # cd4=86, tx=1

S1 <- (1+exp(z1))^-1
S2 <- (1+exp(z2))^-1
S3 <- (1+exp(z3))^-1
S4 <- (1+exp(z4))^-1

yrange <- range(1-S1)

plot(xrange, yrange, type="n", xlab="Days since admission",
     ylab="1-survival", las=1, ylim=c(0,0.1)) 

lines(t, 1-S1, type="l", col=1, lty=2, lwd=2)
lines(t, 1-S2, type="l", col=1, lwd=2)
lines(t, 1-S3, type="l", col=2, lty=2, lwd=2)
lines(t, 1-S4, type="l", col=2, lwd=2)

legend(x="bottomright", col=c(1,1,2,2), lty=c(2,1,2,1), lwd=2,
       c('cd4=43, tx=0', 'cd4=86, tx=0', 'cd4=43, tx=1','cd4=86, tx=1' ))


