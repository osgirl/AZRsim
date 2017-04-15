# basic simulation tests

library(deSolve)
context("SHO Test")

test_that("AZRmodel runs", {
  model <- AZRmodel(system.file("examples/sho.txt", package="AZRsim"))
  AZRmodel_sho <- AZRsimulate(model,100)
})

test_that("AZRmodel and deSolve return the same simulation for a SHO model", {
  model <- AZRmodel(system.file("examples/sho.txt", package="AZRsim"))
  AZRmodel_sho <- AZRsimulate(model,100)
  sho <- function(t, y, parms) {
    with(as.list(c(y,parms)), {
      dy2_dt = y2
      dy1_dt = - y1 - theta * y2
      return(list(c(dy2_dt, dy1_dt)))
    })
  }
  pars <- c("theta" = 0.15)
  yini <- c("y1" = 1, "y2" = 0)
  deSolve_sho <- deSolve::ode(y = yini, times = AZRmodel_sho$TIME, func = sho, parms = pars)
  expect_equal(AZRmodel_sho$y1, deSolve_sho[,2], tolerance = 1e-3)
})
