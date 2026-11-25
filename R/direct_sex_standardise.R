#' Direct Sex Standardisation with Confidence Intervals
#'
#' Computes directly sex-standardised rates using a 50:50 male/female
#' reference population. Confidence intervals are calculated using
#' Byar's approximation. Lower confidence interval limits below zero
#' are reported as zero.
#'
#' @param data Dataframe containing observed counts and population.
#' @param var One or more grouping variable names supplied as a character
#'   vector, e.g. c("Region", "FinancialYear").
#' @param sex Column containing sex categories.
#' @param observed Column containing observed counts.
#' @param population Column containing denominator population.
#' @param multiplier Scale for rates, default is 1,000.
#' @param ref_population Reference population containing Sex and
#'   StdPopulation columns. Defaults to `std_sex`.
#'
#' @return A tibble containing the grouping variables, sex-standardised
#'   rate, and 95% confidence intervals.
#'
#' @examples
#' example_data <- data.frame(
#'   variable = c("A", "A", "B", "B"),
#'   sex = c("F", "M", "F", "M"),
#'   observed = c(10, 20, 15, 25),
#'   population = c(1000, 1200, 1100, 1300)
#' )
#'
#' direct_sex_standardise(
#'   data = example_data,
#'   var = "variable",
#'   sex = "sex",
#'   observed = "observed",
#'   population = "population"
#' )
#'
#' @export

direct_sex_standardise <- function(
    data,
    var,
    sex,
    observed,
    population,
    multiplier = 1000,
    ref_population = std_sex) {

  # ------------------------------------------------------------
  # Check grouping variables
  # ------------------------------------------------------------

  if (!is.character(var)) {
    stop("`var` must be a character vector of column names, e.g. ","c('PCN', 'FinancialYear').", call. = FALSE)}

  # ------------------------------------------------------------
  # Convert single-column arguments to symbols
  # ------------------------------------------------------------

  if (!is.character(sex) || length(sex) != 1) {
    sex <- rlang::as_string(rlang::ensym(sex))}

  if (!is.character(observed) || length(observed) != 1) {
    observed <- rlang::as_string(rlang::ensym(observed))}

  if (!is.character(population) || length(population) != 1) {
    population <- rlang::as_string(rlang::ensym(population))}

  # ------------------------------------------------------------
  # Convert column names to symbols
  # ------------------------------------------------------------

  var_syms <- rlang::syms(var)
  sex_sym <- rlang::sym(sex)
  observed_sym <- rlang::sym(observed)
  population_sym <- rlang::sym(population)

  # ------------------------------------------------------------
  # Check reference population
  # ------------------------------------------------------------

  if (!all(c("Sex", "StdPopulation") %in% names(ref_population))) {
    stop("`ref_population` must contain columns named `Sex` and ","`StdPopulation`.", call. = FALSE)}

  # ------------------------------------------------------------
  # Standardisation
  # ------------------------------------------------------------

  out <- data %>%

    # Aggregate observed and population by grouping variables and sex
    dplyr::group_by(dplyr::across(!!!var_syms),!!sex_sym) %>%

    dplyr::summarise(
      obs = sum(!!observed_sym, na.rm = TRUE),
      pop = sum(!!population_sym, na.rm = TRUE),
      .groups = "drop") %>%

    # Calculate sex-specific crude rate
    dplyr::mutate(
      rate = dplyr::if_else(
        pop > 0,
        obs / pop,
        0
      )
    ) %>%

    # Add standard population weights
    dplyr::left_join(
      ref_population,
      by = setNames("Sex", sex)
    ) %>%

    # Calculate expected cases
    dplyr::mutate(
      expected = rate * StdPopulation
    ) %>%

    # Aggregate across sex
    dplyr::group_by(dplyr::across(!!!var_syms)) %>%

    dplyr::summarise(
      dss_rate =
        sum(expected, na.rm = TRUE) /
        sum(StdPopulation, na.rm = TRUE),

      total_observed =
        sum(obs, na.rm = TRUE),

      total_population =
        sum(pop, na.rm = TRUE),

      crude_rate =
        ifelse(
          total_population > 0,
          total_observed / total_population,
          0
        ),

      .groups = "drop"
    ) %>%

    # Byar 95% confidence intervals
    dplyr::mutate(

      lower = dplyr::if_else(
        total_observed == 0,
        0,
        dss_rate *
          (
            1 -
              (1 / (9 * total_observed)) -
              (1.96 / (3 * sqrt(total_observed)))
          )^3
      ),

      upper =
        dss_rate *
        (
          1 -
            (1 / (9 * (total_observed + 1))) +
            (1.96 / (3 * sqrt(total_observed + 1)))
        )^3,

      # Ensure lower CI cannot be negative
      lower = pmax(lower, 0),

      # Apply rate multiplier and rounding
      dss_rate = round(multiplier * dss_rate, 2),
      lower = round(multiplier * lower, 2),
      upper = round(multiplier * upper, 2)
    ) %>%

    dplyr::select(
      dplyr::all_of(var),
      Observed = total_observed,
      Population = total_population,
      CrudeRate = crude_rate,
      SexStandardisedRate = dss_rate,
      LowerCI = lower,
      UpperCI = upper)

  return(out)
}
