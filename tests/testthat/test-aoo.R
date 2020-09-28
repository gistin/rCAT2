test_that("returns a number by default", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  aoo_value <- aoo(points)
  
  expect_type(aoo_value, "double")
})

test_that("all values in same cell works", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  expect_equal(aoo(points), 4)
})

test_that("all points in individual cells works", {
  points <- data.frame(
    X=c(0,1000,0,1000),
    Y=c(0,0,1000,1000)
  )
  
  expect_equal(aoo(points, cellsize=100), 0.04)
})

test_that("returns full cell info", {
  points <- data.frame(
    X=c(200,1100,200,1100),
    Y=c(200,200,1100,1100)
  )
  
  aoo_info <- aoo(points, cellsize=100, returnV="E")
  
  correct_info <- list(area=0.04, nocells=4, rotation=0, xshift=0, yshift=0)
  
  expect_equal(aoo_info, correct_info)
})
