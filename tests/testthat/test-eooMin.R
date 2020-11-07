test_that("returns a number by default", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  min_eoo <- eooMin(points, defaultRadius=100)
  
  expect_type(min_eoo, "double")
})

test_that("calculates correct value for square", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  # error radius of ~353 gives a 500x500 m square
  min_eoo <- eooMin(points, defaultRadius=353.5534)
  
  # should return eoo in km2, so 0.25 km2
  expect_equal(min_eoo, 0.25)
})

test_that("returns an sf with a polygon geometry", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  min_eoo_geom <- eooMin(points, defaultRadius=100, returnV="SF")
  
  expect_s3_class(min_eoo_geom[[1]], "POLYGON")
})

test_that("returns a list with 2 values (min and max)", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  min_eoo <- eooMin(points, defaultRadius=100, returnV="EX")
  
  expect_equal(length(min_eoo), 2)
})

test_that("max value is greater than min value", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  min_eoo <- eooMin(points, defaultRadius=100, returnV="EX")
  
  expect_gte(min_eoo$max, min_eoo$min)
})
