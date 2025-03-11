# we assume that both files, that is experiments.R and models.R are located in
# the same directory

source("models.R")

##########################################################################
# An example of running the experiment.
##########################################################################
example <- function(dataFile = "data/new_data_1.csv", B = 2000){ # use 2000 bootstrap samples to estimate p-value.
  # load data from disk
  d <- fLoadData(dataFile)
  # run the bootstrap likelihood ratio test; it also produces p-value according
  # to the Wilks theorem; this call will take about 200 seconds.
  res <- LRTest(disoutcome ~ age + comasc + log(comadur), # formula describing the model 
                data = d, # this is the 
                # function for fitting a model under H_0, we use bclm function, but any other could be used
                h0ModelFunc = function(formula, data) poclm(formula, data, restricted = TRUE),
                # function for fitting a model under H_1, we use bclm function, but any other could be used
                h1ModelFunc = function(formula, data) soblm(formula, data, restricted = FALSE),
                B = B) # use 2000 bootstrap samples to estimate p-value.
  res
}
# example(B = 500)
#########Results for disoutcome ~ age + comasc + log(comadur)##########
#H_0: soblm H1: poclm; p-value 0.162
#H_0: poclm H1: soblm; p-value 0.814
#H_0: bclm  H1: poclm; p-value 0.142
#H_0: polcm H1: bclm;  p-value 0.722
#H_0: soblm H1: bclm;  p-value 0.386
#H_0: bclm  H1: soblm: p-value 0.636

##########################################################################
##   Good-of-Fit-Test
##########################################################################
#bclm baseline category model
#poclm proportional odds model
#soblm sequence of binomial models

###################Goodness-Fit-of-Test############################
example1 <- function(dataFile = "data/new_data_1.csv", B = 2000){ # use 2000 bootstrap samples to estimate p-value.
  # load data from disk
  d <- fLoadData(dataFile)
  # run the bootstrap likelihood ratio test; it also produces p-value according
  # to the Wilks theorem; this call will take about 200 seconds.
  res <- LRTest(disoutcome ~ age + comasc + log(comadur)+conv+log(lactate)+log(platelet)+log(HRP2+1), # formula describing the model 
                data = d, # this is the 
                # function for fitting a model under H_0, we use bclm function, but any other could be used
                h0ModelFunc = function(formula, data) poclm(formula, data, restricted = TRUE),
                # function for fitting a model under H_1, we use bclm function, but any other could be used
                h1ModelFunc = function(formula, data) soblm(formula, data, restricted = FALSE),
                B = B) # use 2000 bootstrap samples to estimate p-value.
  res
}
#example1(B = 1000)
###Updated covariates with HRP2##### 
#H_0: soblm H1: poclm; p-value 0.280
#H_0: poclm H1: soblm; p-value 0.318
#H_0: bclm  H1: poclm; p-value 0.372
#H_0: polcm H1: bclm;  p-value 0.284
#H_0: soblm H1: bclm;  p-value 0.186
#H_0: bclm  H1: soblm: p-value 0.802

# with HRP2 #H_0: poclm constraint H_1: bclm baseline category model p-value: 0.073
#H_0: poclm constraint H_1: sequence of binomial models p-value: 0.070


##############################################
example2 <- function(dataFile = "data/new_data_1.csv", B = 2000){ # use 2000 bootstrap samples to estimate p-value
  # load data from disk
  d <- fLoadData(dataFile)
  # run the bootstrap likelihood ratio test; it also produces p-value according
  # to the Wilks theorem; this call will take about 200 seconds.
  res <- LRTest(disoutcome ~ age + comasc + log(comadur)+conv+log(lactate)+log(platelet), # formula describing the model 
                data = d, # this is the 
                # function for fitting a model under H_0, we use bclm function, but any other could be used
                h0ModelFunc = function(formula, data) poclm(formula, data, restricted = TRUE),
                # function for fitting a model under H_1, we use bclm function, but any other could be used
                h1ModelFunc = function(formula, data) soblm(formula, data, restricted = FALSE),
                B = B) # use 2000 bootstrap samples to estimate p-value.
  res
}
#example2(B = 500)

###Updated covariates without HRP2##### 
#H_0: soblm H1: poclm; p-value 0.344
#H_0: poclm H1: soblm; p-value 0.236
#H_0: bclm  H1: poclm; p-value 0.448
#H_0: poclm H1: bclm;  p-value 0.166
#H_0: soblm H1: bclm;  p-value 0.154
#H_0: bclm  H1: soblm: p-value 0.882

#without HRP2
#H_0: poclm constraint H_1: sequence of binomial models p-value: 0.042
#H_0: poclm constraint H_1: bclm baseline category model p-value 0.018

