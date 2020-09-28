test_that("rotation gives expected results for unit square", {
  square <- data.frame(
    x=c(0,1,1,0),
    y=c(0,0,1,1)
  )
  
  rotated_square <- data.frame(
    x=c(0,1/sqrt(2),sqrt(2),1/sqrt(2)),
    y=c(0,-1/sqrt(2),0,1/sqrt(2))
  )
  
  expect_equal(rotateP(square, deg2rad(45)), rotated_square)
})


test_that("returned dataframe has same column names", {
  square <- data.frame(
    X=c(0,1,1,0),
    Y=c(0,0,1,1)
  )
  
  rotated <- rotateP(square, deg2rad(45))
  expect_equal(colnames(rotated), colnames(square))
})
