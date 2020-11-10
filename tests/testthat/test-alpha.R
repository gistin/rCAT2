test_that("CRS in m returns area in km^2", {
  set.seed(420)
  
  points <- ptsSquare(19,0.1)
  names(points) <- c("lat","long")
  points <- simProjWiz(points)
  
  area <- aHullMean(points)
  units <- units(area)
  expect_equal(units$numerator, c("km", "km"))
})

test_that("CRS in long/lat returns area in km^2", {
  set.seed(420)
  
  points <- ptsSquare(19, 0.1)
  attr(points, "crs") <- 4326
  
  area <- aHullMean(points)
  units <- units(area)
  expect_equal(units$numerator, c("km", "km"))
})

test_that("Missing CRS returns area in km^2", {
  set.seed(420)
  
  points <- ptsSquare(19, 0.1)
  names(points) <- c("lat", "long")
  points <- simProjWiz(points)
  attr(points, "crs") <- NULL
  
  area <- aHullMean(points)
  units <- units(area)
  expect_equal(units$numerator, c("km", "km"))
})







