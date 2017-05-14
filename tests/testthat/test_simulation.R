# basic simulation tests

library(deSolve)
context("SHO Tests")

sho <- function(t, y, parms) {
  with(as.list(c(y,parms)), {
    dy2_dt = y2
    dy1_dt = - y1 - theta * y2
    return(list(c(dy2_dt, dy1_dt)))
  })
}

test_that("create_model and simulate run", {
  model <- create_model(system.file("examples/sho.txt", package="AZRsim"))
  azr_sho <- simulate(model,100)
})

model <- create_model(system.file("examples/sho.txt", package="AZRsim"))

test_that("simulate `simtime` argument works", {
  time_steps <- seq(1, 50, by=0.1)
  azr_sho <- simulate(model, time_steps)
  expect_equal(azr_sho$TIME, time_steps)
})

test_that("simulate `IC` argument works", {
  azr_sho <- simulate(model, seq(0, 20, by=0.1), IC = c("y1" = 0, "y2" = 3))
  pars <- c("theta" = 0.15)
  yini <- c("y1" = 0, "y2" = 3)
  deSolve_sho <- deSolve::ode(y = yini, times = azr_sho$TIME, func = sho, parms = pars)
  expect_equal(azr_sho$y1, deSolve_sho[,2], tolerance = 1e-3)
})

test_that("simulate `parameters` argument works", {
  azr_sho <- simulate(model, seq(0, 20, by=0.1), parameters = c("theta" = 0.22))
  pars <- c("theta" = 0.22)
  yini <- c("y1" = 1, "y2" = 0)
  deSolve_sho <- deSolve::ode(y = yini, times = azr_sho$TIME, func = sho, parms = pars)
  expect_equal(azr_sho$y1, deSolve_sho[,2], tolerance = 1e-3)
})

test_that("simulate `outputs` argument works", {
  azr_sho <- simulate(model, seq(0, 20, by=0.1), outputs = "y1")
  expect_equal(names(azr_sho), c("TIME", "y1"))
})

test_that("simulate `verbose = TRUE` argument works", {
  azr_sho <- simulate(model, seq(0, 20, by=0.1), verbose = TRUE)
})

test_that("simulate `opt_method_stiff = TRUE` argument works", {
  azr_sho <- simulate(model, seq(0, 20, by=0.1), opt_method_stiff = TRUE)
})

test_that("create_model and deSolve return the same simulation for a SHO model", {
  model <- create_model(system.file("examples/sho.txt", package="AZRsim"))
  azr_sho <- simulate(model,100)
  pars <- c("theta" = 0.15)
  yini <- c("y1" = 1, "y2" = 0)
  deSolve_sho <- deSolve::ode(y = yini, times = azr_sho$TIME, func = sho, parms = pars)
  expect_equal(azr_sho$y1, deSolve_sho[,2], tolerance = 1e-3)
})
