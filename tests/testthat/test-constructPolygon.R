test_that("returns a polygon", {
  X <- c(0,1000,0,1000)
  Y <- c(0,0,1000,1000)
  
  poly_sf <- constructPolygon(X, Y, 4268)
  
  expect_s3_class(poly_sf[[1]], "POLYGON")
})

test_that("returned polygon can be plotted", {
  X <- c(0,1000,0,1000)
  Y <- c(0,0,1000,1000)

  poly_sf <- constructPolygon(X, Y, 4268)
  
  expect_error(plot(poly_sf), NA)
})

test_that("warning if bad crs provided", {
  X <- c(0,1000,0,1000)
  Y <- c(0,0,1000,1000)

  expect_warning(constructPolygon(X, Y, "this is not valid"))
})
