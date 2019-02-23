library(MetaGxBreast)

context("Checking loadBreastEsets")


test_that("ensure datasets and duplicates are properly loaded from the hub and package", {
  esetsAndDuplicates = MetaGxBreast::loadBreastEsets()
  esets = esetsAndDuplicates$esets
  #still class warning regardless of replacing with is
  expect_equal(is(esets[[1]])[1], "ExpressionSet")
  expect_equal(is(esets[[1]]@assayData$exprs)[1], "matrix")
})
