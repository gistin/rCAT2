test_that("returns number by default", {
  points <- data.frame(
    X=runif(10, min=0, max=1000),
    Y=runif(10, min=0, max=1000)
  )
  
  expect_type(aooFixedRotation(points), "double")
})

test_that("correct rotation for a square", {
  points <- data.frame(
    X=c(0,900,0,900),
    Y=c(0,0,900,900)
  )
  
  aoo_info <- aooFixedRotation(points, returnV="E", cellsize=1000)
  
  expect_equal(aoo_info$rotation, 0)
})
