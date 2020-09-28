test_that("returns a polygon", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  poly_sf <- polyCon(points, 4268)
  
  expect_s3_class(poly_sf$geometry[[1]], "POLYGON")
})

test_that("returned polygon can be plotted", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  poly_sf <- polyCon(points, 4268)
  
  expect_success(plot(poly_sf))
})

test_that("errors if bad crs provided", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  expect_error(polyCon(points, "this is not valid"))
})
