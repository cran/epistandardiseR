test_that("direct_age_standardise returns expected output", {

age_test_data <- data.frame(
  variable = rep(c("Gloucester"), each = 19),
  AgeBand = rep(c("0-4", "5-9","10-14","15-19", "20-24", "25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"),
                times = 1),
  Observed = c( 3, 3, 1, 0, 1, 10, 10, 23, 36, 68, 146, 271, 348, 401, 434, 481, 323, 229, 356),
  Population = c( 636,  558,   609, 543, 384, 594, 669,  843,  660,  576, 606, 522, 426, 273, 300, 288, 153, 123, 450))

result  <- direct_age_standardise(
  data = age_test_data,
  var = "variable",
  ageband = "AgeBand",
  observed = "Observed",
  population = "Population",
  multiplier = 1000)

expect_equal(nrow(result), 1)

expect_named(
  result,
  c(
    "variable",
    "Observed",
    "Population",
    "CrudeRate",
    "AgeStandardisedRate",
    "LowerCI",
    "UpperCI"
  )
)

expect_equal(result$Observed, c(3144))
expect_equal(result$Population, c(9213))

expect_equal(result$CrudeRate,  round((3144 / 9213) * 1000, 2))
expect_true(all(result$LowerCI <= result$AgeStandardisedRate))
expect_true(all(result$UpperCI >= result$AgeStandardisedRate))
})
