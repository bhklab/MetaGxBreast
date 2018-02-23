
library(MetaGxBreast)

context("Checking loadBreastEsets")


test_that("ensure datasets and duplicates are properly loaded from the hub and package", {
  esetsAndDuplicates = MetaGxBreast::loadBreastEsets()
  esets = esetsAndDuplicates$esets
  duplicates = esetsAndDuplicates$duplicates
  expect_equal(duplicates[[1]], "DFHCC2.DFHCC2_REF12rep")
  expect_equal(class(esets[[1]])[1], "ExpressionSet")
  expect_equal(class(esets[[1]]@assayData$exprs), "matrix")
  
})
