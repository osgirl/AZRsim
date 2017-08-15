

# install.packages(c("Rcpp", "tidyverse", "shiny", "foreach", "doParallel"))
# install.packages("microbenchmark","mrgsolve")
# install.packages("microbenchmark","mrgsolve")
# path<-"J:/Projects/Programmers/Rpackages/azrtest.7z"


#.libPaths( c( .libPaths("J:/R/R-3.3.3/library"), "J:/Projects/Programmers/Rpackages/Library") )
#.libPaths()[2]



library(AZRsim)
#library(AZRaux)
library(microbenchmark)
library(mrgsolve)
library(tidyverse)
library(knitr)


rm(list = ls())

#########################
# Simple PK model with single dose
#########################


#setwd("J:/Projects/Programmers/Rpackages/AZRsim-master/inst/examples/simulation")
list.files()



model       <- create_model("modelPKtest.txt")
simtime     <- seq(0,200,1)
dosingTable <- data.frame(TIME=0,DOSE=400,DURATION=0,INPUT=1)
results     <- simulate(model,simtime=simtime,dosing_table = dosingTable)
plot(results$TIME,results$Cc,type="l")
shiny_plot(results)



#########################
# Simple PK model with more complex dosing
#########################

model       <- create_model("modelPKtest.txt")
simtime      <- seq(0,200,1)

dosing_INPUT1 <- data.frame(
  TIME     = c(12.5,seq(24,120,24)),
  DOSE     = 120,
  DURATION = 2,
  INPUT    = 1,
  LAGTIME  = 10)

dosing_INPUT2 <- data.frame(
  TIME     = 12,
  DOSE     = 300,
  DURATION = 120,
  INPUT    = 2,
  LAGTIME  = 0)


dosingTable <- rbind(dosing_INPUT1,dosing_INPUT2)

results     <- simulate(model,simtime=simtime,dosing_table=dosingTable)

plot(results$TIME,results$Cc,type="l")

# Not run:
model <- create_model(system.file("examples/NovakTyson.txt", package="AZRsim"))
x <- simulate(model,400)
shiny_plot(x)


#########################
# Assess simulation speed for complex dosing model
#########################

microbenchmark(
  #AZRsimulate(model,dosingTable=dosingTable)
  simulate(model,dosing_table = dosingTable)
  ,times=1000
)


################
# Population simulation
################

model       <- create_model("modelPKtest.txt")

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
  results <- AZRsimpop(model,ncores=1,simtime=seq(0,200),parameterTable=indivParam,dosing_table=dosingTable)
  ,times=1)

microbenchmark(
  results <- AZRsimpop(model,ncores=8,simtime=seq(0,200),parameterTable=indivParam,dosingTable=dosingTable)
  ,times=1)

# RxODE Test Code Segment

ode <- " C2 = centr/V2;
      d/dt(depot) = -KA * depot;
      d/dt(centr) = KA * depot - CL *C2;
        "
library(RxODE)
m1 <- RxODE(ode)
theta <- c(KA = 0.2, CL = 0.4, V2 = 12)
ev <- eventTable()
ev$add.dosing(dose=300, dosing.to = 2, rate = 300/120, start.time = 0)
ev$add.dosing(dose=120, nbr.doses = 4, dosing.interval = 24, rate = 60, start.time = 2)
ev$add.sampling(time = seq(0,200,1))
x <- m1$run(theta, ev)
dosing <- ev$get.dosing()
nsub     <- 1000
#theta.all <- cbind(kA=0.3, CL = 0.4*exp(rnorm(nsub,0,.5)),V2=12)
centerTrans2   <- log(theta)
theta.all    <- exp(t(replicate(n=nsub,centerTrans2 + rnorm(mean=0,sd=0.5,n=length(centerTrans2)))))
cp.all = matrix(NA, 4, nsub)

# mrgsolve has both analytical solutions and odes
# for such a simple example as henning provided they should be similar in results
# code files can be defined as regular files (preferred), but using inline here as so simple
ana_code <- '
$PARAM CL = 0.4, V = 12, KA=0.2
$CMT GUT CENT
$GLOBAL
#define CP (CENT/V)
$PKMODEL ncmt=1, depot=TRUE
$CAPTURE CP
'

ode_code <- '
$PARAM CL = 0.4, V = 12, KA=0.2
$CMT GUT CENT
$GLOBAL
#define CP (CENT/V)
$ODE
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - CP*CL/V;
$CAPTURE CP'

# compile all models
ana_mod <- mcode("ana", ana_code)
ode_mod <- mcode("ode", ode_code)

# mrgsolve dosing takes a nonmem-style dataframe, though can also use the `ev()` for simple dosing.
# For example `ev(amt = 100, cmt = 1, addl = 5, ii = 24)`. In this case we will construct
# an entire dataset.

dosing_df <- map_df(1:1000, function(id) {
  # create code that is a reasonable approximation of hennings - I'm not sure what he is demoing with the
  # lag time exactly
  # this code and his code both give doses into two compartments
  bind_rows(
    data_frame(ID = 1, time = 0, cmt = 2, amt=300, addl = 0, ii = 0, rate = 300/120, evid = 1),
    data_frame(ID = 1, time = 0.5, cmt = 1, amt=120, addl = 5, ii = 24, rate = 60, evid = 1)
  ) %>% mutate(ID = id) })

model <- create_model("modelPKtest.txt")

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

# mrgsolve supports nonmem conventions of `$OMEGA` and `$SIGMA`, however
# these must be simulated outside in azrsim. Luckily, mrgsolve also allows explicitly overriding values
# for each ID, so we can use the `indivParam` by passing it to `idata_set`

mrg_params <- indivParam %>% as_data_frame %>% rename(KA = ka, V = Vc) %>% mutate(ID = 1:n())


# lets wrap everything in a function just so names are easy to read in microbenchmark output
analytical <- function() {
  ana_mod %>%
    data_set(dosing_df)%>%
    idata_set(mrg_params) %>%
    mrgsim(end = 200)
}

ode_mrgsolve <- function() {
  ode_mod %>% data_set(dosing_df) %>%
    idata_set(mrg_params) %>%
    mrgsim(end = 200)
}
result_rxode = matrix(data = NA,1,4)

ode_rxode <- function(){
  for(i in 1:nsub)
  {
    theta = theta.all[i,]
    x <- m1$run(theta, ev)
    result_rxode <- rbind(result_rxode, x)

  }
}

ode_azrsim <- function() {
  AZRsimpop(model,ncores=1,simtime=seq(0,200),parameterTable=indivParam,dosing_table=dosingTable)
}

res <- microbenchmark::microbenchmark(
  analytical(),
  ode_mrgsolve(),
  ode_rxode(),
  ode_azrsim(),
  times = 10L
)


res

ggplot2::autoplot(res)














