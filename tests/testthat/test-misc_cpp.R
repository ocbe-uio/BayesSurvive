test_that("list_to_matrix works", {
  r_list <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9))
  result_matrix <- list_to_matrix(r_list)
  expect_equal(result_matrix, matrix(1:9, nrow = 3, ncol = 3))
})
