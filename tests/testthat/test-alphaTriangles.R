test_that("CRS in long/lat raises warning", {
  set.seed(420)
  
  points <- ptsSquare(19, 0.1)
  
  expect_warning(alphaTriangles(points$X, points$Y, crs=4326))
})

test_that("Missing CRS raises warning", {
  set.seed(420)
  
  points <- ptsSquare(19, 0.1)
  
  expect_warning(alphaTriangles(points$X, points$Y, crs=NULL), 
                 "(M|m)issing CRS")
})