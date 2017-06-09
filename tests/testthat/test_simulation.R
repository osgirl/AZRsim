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

test_that("a dosing table in simulate.azrmod and an event data frame in deSolve return the same simulation", {
          model <- create_model(system.file("examples/one_cpt_dt.txt", package="AZRsim"))
          dt <- data.frame("TIME" = seq(1,9, by = 1),
                           "DOSE" = 40,
                           "DURATION" = 0,
                           "INPUT" = 1,
                           "LAGTIME" = 0,
                           stringsAsFactors = FALSE)
          one_cpt_dt_simulation <- simulate(model, 10, dosingTable = dt)

          one_cpt <- function(t, y, p) {
            with(as.list(c(y,p)), {
            dydt = -theta * y
            return(list(dydt))
            })
          }
          eventdf <- data.frame(var = "y", time = dt$TIME, value = dt$DOSE, method = "add")
          deSolve_one_cpt <- deSolve::ode(y = c(y=0), times = one_cpt_dt_simulation$TIME,
                                          func = one_cpt, parms = c(theta = 0.6), events = list(data = eventdf))
          expect_equal(one_cpt_dt_simulation[,"y"], deSolve_one_cpt[,"y"], tol = 0.001)
          })
