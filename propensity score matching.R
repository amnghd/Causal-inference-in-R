library(tableone) # for matching comparison
library(MatchIt) # propensity score mathcing with plots
library(Matching) # carries out the matching
library(gtools)

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

# new dataset
# we will setup a set of confounders and treatment variables

mydata <- cbind(ARF, CHF, Cirr, colcan, lungcan,
                MOSF, sepsis, age, female, treatment,
                died,aps)

mydata <-data.frame(mydata)


# let's first perform propensity score calculaton using logistic regression

psmodel <- glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+
                 MOSF+sepsis+age+female+meanbp1+aps,
               family = binomial(), data=mydata) # binomial informs using logistic regression

# show coeficients
summary(psmodel)

# create propensity score
pscore <- psmodel$fitted.values

# now let's take a look at the distribution of treatment
# there are a lot of overlap


# now use matchit for propensity score matching

m.out <- matchit(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+
                 MOSF+sepsis+age+female+meanbp1+aps,
               data=mydata,
               method='nearest')

summary(m.out)


plot(m.out, type='jitter')
plot(m.out, type='hist')


# just doing another matching, taking caliper into account
# doing greedy matching
psmatch <- Match(Tr=mydata$treatment, M=1, X=logit(pscore), replace=FALSE)
matched <-mydata[unlist(psmatch[c("index.treated","index.control")]),]
xvars <-c("ARF","CHF", "Cirr", "colcan", "lungcan",
          "MOSF", "sepsis", "age", "female", "aps")
matchedtab1 <- CreateTableOne(vars=xvars, strata = "treatment",
                              data=matched, test=FALSE)
print(matchedtab1, smd=TRUE)

# ttest comparison
# now let's analyze the outcome
y_trt <- matched$died[matched$treatment==1]
y_con <- matched$died[matched$treatment==0]
print("performing t test ...")
t.test(y_trt, y_con, paired=TRUE)
# we are not happy with the match
# we see alot of SMD >0.1


# how about using caliper
psmatch <- Match(Tr=mydata$treatment, M=1, X=logit(pscore), replace=FALSE, caliper=0.2)
matched <-mydata[unlist(psmatch[c("index.treated","index.control")]),]
xvars <-c("ARF","CHF", "Cirr", "colcan", "lungcan",
          "MOSF", "sepsis", "age", "female", "aps")
matchedtab1 <- CreateTableOne(vars=xvars, strata = "treatment",
                              data=matched, test=FALSE)
print(matchedtab1, smd=TRUE)

# ttest comparison
# now let's analyze the outcome
y_trt <- matched$died[matched$treatment==1]
y_con <- matched$died[matched$treatment==0]
print("performing t test ...")
t.test(y_trt, y_con, paired=TRUE)



