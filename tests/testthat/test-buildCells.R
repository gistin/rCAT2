test_that("returns polygons", {
  cell_coords <- data.frame(
    x=c(0,900,0,900),
    y=c(0,0,900,900)
  )
  
  cells <- buildCells(cell_coords, 100, 0, 0, 0, crs=4628)
  expect_s3_class(cells[[1]], "POLYGON")
})

test_that("cells have the specified CRS", {
  cell_coords <- data.frame(
    x=c(0,900,0,900),
    y=c(0,0,900,900)
  )
  
  cells <- buildCells(cell_coords, 100, 0, 0, 0, crs=4628)
  expect_equal(st_crs(cells)$epsg, 4628)
})

test_that("cells can be plotted", {
  cell_coords <- data.frame(
    x=c(0,900,0,900),
    y=c(0,0,900,900)
  )
  
  cells <- buildCells(cell_coords, 100, 0, 0, 0, crs=4628)
  expect_error(plot(cells), NA)
})
