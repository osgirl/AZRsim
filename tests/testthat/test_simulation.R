# basic simulation tests

library(deSolve)
context("SHO Test")

test_that("create_model runs", {
  model <- create_model(system.file("examples/sho.txt", package="AZRsim"))
  sho <- simulate(model,100)
})

test_that("create_model and deSolve return the same simulation for a SHO model", {
  model <- create_model(system.file("examples/sho.txt", package="AZRsim"))
  azr_sho <- simulate(model,100)
  sho <- function(t, y, parms) {
    with(as.list(c(y,parms)), {
      dy2_dt = y2
      dy1_dt = - y1 - theta * y2
      return(list(c(dy2_dt, dy1_dt)))
    })
  }
  pars <- c("theta" = 0.15)
  yini <- c("y1" = 1, "y2" = 0)
  deSolve_sho <- deSolve::ode(y = yini, times = azr_sho$TIME, func = sho, parms = pars)
  expect_equal(azr_sho$y1, deSolve_sho[,2], tolerance = 1e-3)
})
