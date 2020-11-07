test_that("returns a polygon", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  poly_sf <- constructPolygon(points, 4268)
  
  expect_s3_class(poly_sf[[1]], "POLYGON")
})

test_that("returned polygon can be plotted", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  poly_sf <- constructPolygon(points, 4268)
  
  expect_error(plot(poly_sf), NA)
})

test_that("warning if bad crs provided", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  expect_warning(constructPolygon(points, "this is not valid"))
})
