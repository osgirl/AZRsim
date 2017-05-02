# test that functions return objects of the correct class

context("azrmod/azrsim")

test_that("create_model returns an object of class 'azrmod', simulate returns an object of class 'azrsim'", {
  model <- create_model(system.file("examples/sho.txt", package="AZRsim"))
  sho <- simulate(model,100)
  expect_equal(class(model), "azrmod")
  expect_equal(class(sho), "azrsim")
})
