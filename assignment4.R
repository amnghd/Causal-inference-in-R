library(tableone)
library(Matching)
library(ipw)
library(survey)
library(MatchIt)
data(lalonde)

# Q1
xvars = c('age', 'educ', 'black', 'hispan', 'married', 'nodegree',
          're74', 're75')

# let's get the propensity score model
psmodel <- glm(treat~.-re78, data=lalonde,
               family= binomial(link='logit'))

# let get all the propensity scores
ps <- predict(psmodel, type='response') #type:predited value(proba)

print(min(ps))
print(max(ps))


#Q2
# standardized difference for weighted coefficient
# I will use ipw
weightmodel <- ipwpoint(exposure = re78,
                        family='binomial',
                        link='logit',
                        denominator=~age+educ+black+hispan+married+
                        nodegree+re74+re75,
                        data=lalonde)
# denminator of weight=propensity score
summary(weightmodel$ipw.weights)
# it calculates the weights



#develop the weighted data
weighted_data <- svydesign(ids = ~1, data=lalonde,
                           weights = ~ weightmodel$ipw.weights)


# weighted table 1
# this is where the magic happends of creating table one
weightedtable <- svyCreateTableOne(vars=xvars, strata='re78',
                                   data=weighted_data, test = FALSE)
# show table with SMD
print(weightedtable, smd=TRUE)