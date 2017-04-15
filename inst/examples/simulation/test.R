library(AZRsim)
library(AZRaux)
library(microbenchmark)

rm(list = ls())

#########################
# Simple PK model with single dose
#########################

model       <- AZRmodel("modelPKtest.txt")
simtime     <- seq(0,200,1)

dosingTable <- data.frame(TIME=0,DOSE=400,DURATION=2,INPUT=1)

results     <- AZRsimulate(model,simtime=simtime,dosingTable=dosingTable)

plot(results$TIME,results$Cc,type="l")

#########################
# Simple PK model with more complex dosing
#########################

model        <- AZRmodel("modelPKtest.txt")
simtime      <- seq(0,200,1)

dosing_INPUT1 <- data.frame(
  TIME     = c(12.5,seq(24,120,24)),
  DOSE     = 120,
  DURATION = 2,
  INPUT    = 1,
  LAGTIME  = 2)

dosing_INPUT2 <- data.frame(
  TIME     = 12,
  DOSE     = 300,
  DURATION = 120,
  INPUT    = 2,
  LAGTIME  = 0)

dosingTable <- rbind(dosing_INPUT1,dosing_INPUT2)

results     <- AZRsimulate(model,simtime=simtime,dosingTable=dosingTable)

plot(results$TIME,results$Cc,type="l")

#########################
# Assess simulation speed for complex dosing model
#########################

microbenchmark(
  AZRsimulate(model,dosingTable=dosingTable)
  ,times=1000
)


################
# Population simulation
################

model <- AZRmodel("modelPKtest.txt")

dosing_INPUT1 <- data.frame(
  TIME     = c(12.5,seq(24,120,24)),
  DOSE     = 120,
  DURATION = 2,
  INPUT    = 1,
  LAGTIME  = 2)

dosing_INPUT2 <- data.frame(
  TIME     = 12,
  DOSE     = 300,
  DURATION = 120,
  INPUT    = 2,
  LAGTIME  = 0)

dosingTable <- rbind(dosing_INPUT1,dosing_INPUT2)

# Determine individual parameters by sampling
Nsubjects     <- 1000
center        <- c(ka=0.3, CL=0.4, Vc=12)
centerTrans   <- log(center)
indivParam    <- exp(t(replicate(n=Nsubjects,centerTrans + rnorm(mean=0,sd=0.5,n=length(centerTrans)))))

microbenchmark(
  results <- AZRsimpop(model,ncores=1,simtime=seq(0,200),parameterTable=indivParam,dosingTable=dosingTable)
  ,times=1)

microbenchmark(
  results <- AZRsimpop(model,ncores=8,simtime=seq(0,200),parameterTable=indivParam,dosingTable=dosingTable)
  ,times=1)
