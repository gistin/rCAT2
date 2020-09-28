test_that("zero error radius returns original edge point", {
  expect_equal(l_c_intercepts(c(-1,-1), c(1,1), 0), c(1,1,1,1))
})

test_that("intersecting midpoint returns midpoint", {
  expect_equal(l_c_intercepts(c(-1,0), c(0,0), 1), c(-1,0,1,0))
})

test_that("distance between returned points is same as diameter", {
  midpoint <- c(0,0)
  point <- c(4,0)
  radius <- 2
  
  intersection <- l_c_intercepts(midpoint, point, radius)
  expect_equal(intersection[3] - intersection[1], 2*radius)
})