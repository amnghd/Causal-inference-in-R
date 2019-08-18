#install.packages("tableone")
#install.packages("Matching")

#load packages
library(tableone)  # useful for performing match analysis
# it allows checing the balanacing
library(Matching) # carries out the matching

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


# new dataset
# we will setup a set of confounders and treatment variables

mydata <- cbind(ARF, CHF, Cirr, colcan, lungcan,
                MOSF, sepsis, age, female, treatment,
                died)

mydata <-data.frame(mydata)

# list of covariates
# we are just going with a smaller set of confounders
xvars <-c("ARF","CHF", "Cirr", "colcan", "lungcan",
          "MOSF", "sepsis", "age", "female")

# let's create a table one prematch
table1 <- CreateTableOne(vars=xvars, strata='treatment', data=mydata, test=FALSE)
# print the results including the standardized mean differences (SMD)
print("pre matched table 1")
print(table1,smd=TRUE)

# now let's perform matching
print("matching ...")
greedymatch <- Match(Tr=treatment, M=1, X=mydata[xvars])
# Tr is the treatment, M=1 pair match or many to one
# above generates indices for control and treatment slices of main data
indices = unlist(greedymatch[c("index.treated", "index.control")])
matched <- mydata[indices,]

# now let's get the table one on matched data
print("table1 for matching")
table2 <- CreateTableOne(vars=xvars, strata='treatment', data=matched, test=FALSE)
# print the results including the standardized mean differences (SMD)
print("pre matched table 1")
print(table2,smd=TRUE)


# now let's analyze the outcome
y_trt <- matched$died[matched$treatment==1]
y_con <- matched$died[matched$treatment==0]

# pairwise differences
diffy <- y_trt - y_con

#performing the t-test
# null hypothesis: there is no differences
# alt hypothesis: there is some differences
print("performing t test ...")
t.test(diffy)


print("similar results for paired test ...")
t.test(y_trt, y_con, paired=TRUE)



# let's perform similar test using mcnemar test
print(table(y_trt, y_con))
ex <- matrix(c(1414, 750, 583, 423), 2,2)
mcnemar.test(ex)