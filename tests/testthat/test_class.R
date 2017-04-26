# test that functions return objects of the correct class

context("azrmod")

test_that("create_model returns an object of class 'azrmod'", {
  model <- create_model(system.file("examples/sho.txt", package="AZRsim"))
  expect_equal(class(model), "azrmod")
})

# context("simulate")
#
# test_that("simulate returns an object of class 'azrsim", {
#
# })
