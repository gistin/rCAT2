test_that("returns correct length data frame", {
  cell_coords <- data.frame(
    x=c(0,900,0,900),
    y=c(0,0,900,900)
  )
  
  cells <- buildCellPolys(cell_coords, 100)
  expect_length(cells$id, nrow(cell_coords) * 5)
})
