test_that("error returned for wrong column names", {
  points <- data.frame(
    lat=c(0,1000,0,1000),
    long=c(0,0,1000,1000)
  )
  
  expect_error(eoo(points))
})

test_that("returns a number by default", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  expect_type(eoo(points), "double")
})

test_that("calculates km-squared correctly", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  expect_equal(eoo(points, returnV="S"), 1)
})

test_that("returns sf class for polygon mode", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  expect_s3_class(eoo(points, returnV="SF"), "sf")
})

test_that("returns an sf with a polygon geometry", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  eoo_sf <- eoo(points, returnV="SF")
  
  expect_s3_class(eoo_sf$geometry[[1]], "POLYGON")
})

test_that("returned sf can be plotted", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  eoo_sf <- eoo(points, returnV="SF")
  
  expect_success(plot(eoo_sf))
})