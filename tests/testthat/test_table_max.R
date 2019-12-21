context("table_max")

test_that("consistency of output size", {

  df = data.frame("PRCP" = c(1, 2, 0, 1, 0, 0, 3, 3, 0),
              "YEAR" = c(1998, 1998, 1998, 1999, 1999, 1999, 2000, 2000, 2000))
  output <- table_max(df, Nt = 3, reshuffle_days = FALSE)

  expect_equal( length(output$Xi), length(output$Fi))

})
