test_that("warns for lots of points", {
  points <- data.frame(
    X=runif(100, min=0, max=1000),
    Y=runif(100, min=0, max=1000)
  )

  expect_warning(aooFixedGrid(points))
})

test_that("returns number by default", {
  points <- data.frame(
    X=runif(10, min=0, max=1000),
    Y=runif(10, min=0, max=1000)
  )
  
  expect_type(aooFixedGrid(points), "double")
})

test_that("gets correct size of offset points", {
  points <- data.frame(
    X=c(200,1100,200,1100),
    Y=c(200,200,1100,1100)
  )
  
  expect_equal(aooFixedGrid(points, cellsize=1000), 1)
})

test_that("returns full shift info", {
  points <- data.frame(
    X=c(200,1100,200,1100),
    Y=c(200,200,1100,1100)
  )
  
  aoo_info <- aooFixedGrid(points, cellsize=1000, returnV="E")
  
  correct_info <- list(area=1, nocells=1, rotation=0, xshift=200, yshift=200)
  
  expect_equal(aoo_info, correct_info)
})

test_that("returns sf polygon", {
  points <- data.frame(
    X=c(200,1100,200,1100),
    Y=c(200,200,1100,1100)
  )
  
  aoo_poly <- aooFixedGrid(points, cellsize=1000, returnV="SF")
  
  expect_s3_class(aoo_poly$geometry[[1]], "POLYGON")
})

test_that("full results dataframe has right number of rows", {
  points <- data.frame(
    X=c(200,1100,200,1100),
    Y=c(200,200,1100,1100)
  )
  
  aoo_info <- aooFixedGrid(points, cellsize=1000, returnV="ALL")
  
  expect_equal(nrow(aoo_info), nrow(points)^2)
})

