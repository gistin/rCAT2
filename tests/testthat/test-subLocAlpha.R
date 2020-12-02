test_that("returns an integer by default", {
  set.seed(420)
  
  points <- ptsSquare(19,0.1)
  names(points) <- c("lat","long")
  points <- simProjWiz(points)
  
  subpops <- subLocAlpha(points)
  
  expect_type(subpops, "integer")
})
