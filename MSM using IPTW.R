# we are working through a marginial structural model
# using inverse probablity of treatment weighting
# we use the same r heart categorization


# import required packages
library(tableone)
library(ipw)  # used to get IPTW
library(sandwich) # for robust variance estimation
library(survey)  # it is defined for survey, we use for calculating weights
library(ggplot2)

#read in data
load(url("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.sav"))

# viewing the data
View(rhc)



# converting variable to numeric
# one hot encoding

ARF <- as.numeric(rhc$cat1=='ARF')
CHF <- as.numeric(rhc$cat1=='CHF')
Cirr <- as.numeric(rhc$cat1=='Cirrhosis')
colcan <- as.numeric(rhc$cat1 =='Colon Cancer')
Coma <- as.numeric(rhc$cat1=='Coma')
COPD <- as.numeric(rhc$cat1=='COPD')
lungcan <- as.numeric(rhc$cat1=='Lung Cancer')
MOSF <- as.numeric(rhc$cat1=='MOSF w/Malignancy')
sepsis <- as.numeric(rhc$cat1=='MOSF w/Sepsis')
female <- as.numeric(rhc$sex=='Female')
died <- as.numeric(rhc$death=='Yes')
age<-rhc$age
treatment <- as.numeric(rhc$swang1=='RHC')
meanbp1 <- rhc$meanbp1
aps <- rhc$aps1


# we will setup a set of confounders and treatment variables

mydata <- cbind(ARF, CHF, Cirr, colcan, Coma,lungcan,
                MOSF, sepsis, age, female, meanbp1, treatment,
                died,aps)

mydata <-data.frame(mydata)

# let's get the propensity score model
psmodel <- glm(treatment ~ age+female+meanbp1+ARF+CHF+
                 Cirr+colcan+Coma+lungcan+
                 MOSF+sepsis+aps,
               family= binomial(link='logit'))

# let get all the propensity scores
ps <- predict(psmodel, type='response') #type:predited value(proba)

summary(psmodel)


# create weights
weights <-ifelse(treatment==1, 1/ps, 1/(1-ps))

# apply weight to data
# this is where the magic of creating pseudo-population happends

weighted_data <- svydesign(ids = ~1, data=mydata, weights = ~ weights)

# define X variables
xvars <-c("ARF","CHF", "Cirr", "colcan", "lungcan",
          "MOSF", "sepsis", "age", "female", "aps","meanbp1","Coma")




# weighted table 1
# this is where the magic happends of creating table one
weightedtable <- svyCreateTableOne(vars=xvars, strata='treatment',
                                   data=weighted_data, test = FALSE)
# show table with SMD
print(weightedtable, smd=TRUE)


# MSM
# causal relative risk (weighted GLM)
# log is used for peforming relative risk
glm.obj <- glm(died~treatment, weights = weights,
               family=binomial(link=log))
# my goal is to see how outcome is related to treatment
# get a summary of glm
summary(glm.obj)
betaiptw <- coef(glm.obj)
# to take weighting into account use asymptotic weighting (sandwitch)
SE <- sqrt(diag(vcovHC(glm.obj, type='HC0'))) # to account for fact that
# we have inflated the population, we account for it

# getting point estimate and CI for relative risk
causalrr <- exp(betaiptw[2])
lcl <- exp(betaiptw[2] - 1.96 * SE[2])
ucl <- exp(betaiptw[2] + 1.96 * SE[2])
c(lcl, causalrr, ucl)

# now let's get the causal risk difference
# we use identity link function for it
glm.obj <- glm(died~treatment, weights = weights,
               family=binomial(link='identity'))
# my goal is to see how outcome is related to treatment
# get a summary of glm
summary(glm.obj)
betaiptw <- coef(glm.obj)
# to take weighting into account use asymptotic weighting (sandwitch)
SE <- sqrt(diag(vcovHC(glm.obj, type='HC0'))) # to account for fact that
# we have inflated the population, we account for it

causalrr <- (betaiptw[2])
lcl <- (betaiptw[2] - 1.96 * SE[2])
ucl <- (betaiptw[2] + 1.96 * SE[2])
c(lcl, causalrr, ucl)


# we also can simply use the ipw package to do all we did
weightmodel <- ipwpoint(exposure = treatment,
                        family='binomial',
                        link='logit',
                        denominator=~age+female+meanbp1+ARF+
                          CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis,
                        data=mydata)
# denminator of weight=propensity score

# it calculates the weights
# calculates the propensity score
# calculate the MGM

summary(weightmodel$ipw.weights)

# we can get the plot of weights
ipwplot(weights = weightmodel$ipw.weights,
        logscale = FALSE,
        main='weights',
        xlim=c(0,22))

# now that you got th weights
# you can fit the marginal structural model

msm <- (svyglm(died~treatment, 
              design=svydesign(~1, weights=~weightmodel$ipw.weights,
                               data=mydata)))
# this performs the sandwitch for you
# it only needs weights from ipw

coef(msm)
confint(msm)
# treatment   0.02333029 0.07976873
# with ipw and survey we got the same results

# let's see what to do what to do for truncation
truncweight <- replace(weight, weight>10, 10)

# there is another way
#ipwpoint has a trucated option, 
# trunc=0.01 this is a percentile (1-99 th percentile)