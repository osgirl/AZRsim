# test that model blocks work appropriately

library(deSolve)

context("test that 'MODEL FUNCTIONS' works")

test_that("custom AZRsim function and deSolve return the same simulation for a SHO model", {
  model <- create_model(system.file("examples/sho_func.txt", package="AZRsim"))
  azr_sho <- simulate(model,100)
  sho <- function(t, y, parms) {
    with(as.list(c(y,parms)), {
      dy2_dt = y2
      dy1_dt = - y1 - theta * (y2 * 0.01)
      return(list(c(dy2_dt, dy1_dt)))
    })
  }
  pars <- c("theta" = 0.15)
  yini <- c("y1" = 1, "y2" = 0)
  deSolve_sho <- deSolve::ode(y = yini, times = azr_sho$TIME, func = sho, parms = pars)
  expect_equal(azr_sho$y1, deSolve_sho[,2], tolerance = 1e-3)
})

context("test that 'MODEL REACTIONS/VARIABLES' works")

test_that("custom AZRsim reactions/variables and deSolve return the same simulation for a SHO model", {
  model_react <- create_model(system.file("examples/sho_react.txt", package="AZRsim"))
  model_var <- create_model(system.file("examples/sho_var.txt", package="AZRsim"))
  azr_sho_react <- simulate(model_react,100)
  azr_sho_var <- simulate(model_var,100)
  sho <- function(t, y, parms) {
    with(as.list(c(y,parms)), {
      dy2_dt = y2
      dy1_dt = - y1 - theta * y2
      return(list(c(dy2_dt, dy1_dt)))
    })
  }
  pars <- c("theta" = 0.15)
  yini <- c("y1" = 1, "y2" = 0)
  deSolve_sho <- deSolve::ode(y = yini, times = azr_sho_react$TIME, func = sho, parms = pars)
  expect_equal(azr_sho_react$y1, deSolve_sho[,2], tolerance = 1e-3)
  expect_equal(azr_sho_var$y1, deSolve_sho[,2], tolerance = 1e-3)
})

context("test that dosing works")

one_cpt <- function(t, y, parms) {
  with(as.list(c(y,parms)), {
    dy_dt = -theta * y
    return(list(c(dy_dt)))
  })
}
pars <- c("theta" = 0.6)
yini <- c("y" = 0)
injectevents <- data.frame(var = "y",
                           time = 1:9,
                           value = 40,
                           method = "add")
deSolve_one_cpt <- deSolve::ode(y = yini, times = seq(0, 10, by=0.01), func = one_cpt, parms = pars, events = list(data = injectevents))
drop_obs <- seq(101, 901, by = 100)

test_that("AZRsim events and deSolve events return the same discrete dosing simulation", {
  model <- create_model(system.file("examples/one_cpt.txt", package="AZRsim"))
  azr_one_cpt <- simulate(model, seq(0, 10, by=0.01))
  expect_equal(azr_one_cpt$y[-drop_obs], deSolve_one_cpt[-drop_obs,2], tolerance = 1e-3)
})

test_that("AZRsim dosing and deSolve events return the same dosing table simulation", {
  model <- create_model(system.file("examples/one_cpt_dt.txt", package="AZRsim"))
  dt <- data.frame("TIME" = seq(1,9, by = 1),
                   "DOSE" = 40,
                   "DURATION" = 0,
                   "INPUT" = 1,
                   "LAGTIME" = 0,
                   stringsAsFactors = FALSE)
  azr_one_cpt_dt <- simulate(model, seq(0, 10, by=0.01), dosing_table = dt, output = c("y"))
  expect_equal(azr_one_cpt_dt$y[-drop_obs], deSolve_one_cpt[-drop_obs,2], tolerance = 1e-3)
})
