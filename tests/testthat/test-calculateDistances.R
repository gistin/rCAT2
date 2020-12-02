test_that("returns numbers", {
  points <- ptsSquare(19,0.1)
  names(points) <- c("lat","long")
  points <- simProjWiz(points)
  
  crs <- attr(points, "crs")
  triangles <- alphaTriangles(points$X, points$Y, crs)
  
  distances <- calculateDistances(triangles[[1]])
  expect_type(distances, "double")
})

test_that("returns 3 values", {
  points <- ptsSquare(19,0.1)
  names(points) <- c("lat","long")
  points <- simProjWiz(points)
  
  crs <- attr(points, "crs")
  triangles <- alphaTriangles(points$X, points$Y, crs)
  
  distances <- calculateDistances(triangles[[1]])
  expect_length(distances, 3)
})
