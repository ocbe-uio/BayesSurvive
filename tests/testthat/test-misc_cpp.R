test_that("list_to_matrix works", {
  r_list <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9))
  result_matrix <- list_to_matrix(r_list)
  expect_equal(result_matrix, matrix(1:9, nrow = 3, ncol = 3))
})

test_that("list_to_cube works", {
  r_list <- list(matrix(1:6, 3), matrix(7:12, 3), matrix(13:18, 3))
  result_cube <- list_to_cube(r_list)
  expect_equal(result_cube, array(1:18, dim = c(3, 2, 3)))
})

test_that("list_to_vector works", {
  r_list <- as.list(1:5)
  result_vector <- list_to_vector(r_list)
  expect_equal(result_vector, matrix(1:5))
})
