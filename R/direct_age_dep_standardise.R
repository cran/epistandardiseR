
#' Direct Age and Deprivation Standardisation with Confidence Intervals
#'
#' Computes directly age and deprivation standardised rates with 95%
#' confidence intervals using the 2013 European Standard Population
#' with an equal 10% distribution across deprivation deciles.
#' Confidence intervals are calculated using Byar's method. Lower
#' confidence interval limits below zero are reported as zero.
#'
#' @param data Dataframe containing observed counts and population.
#' @param var One or more grouping variables (quoted or unquoted).
#' @param ageband Column specifying age bands.
#' @param deprivation Column specifying deprivation deciles.
#' @param observed Column containing observed counts.
#' @param population Column containing denominator population.
#' @param multiplier Scale for rates (default = 1,000).
#'
#' @return A tibble containing the grouping variables, total observed
#'   events, total population, crude rate, age- and deprivation-
#'   standardised rate, and 95% confidence intervals.
#'
#' @examples
#' example_data <- data.frame(
#'   variable = c("A", "A", "A", "A"),
#'   ageband = c("0-4", "0-4", "5-9", "5-9"),
#'   deprivation = c(1, 2, 1, 2),
#'   observed = c(10, 12, 20, 22),
#'   population = c(1000, 1100, 1200, 1300)
#' )
#'
#' direct_age_dep_standardise(
#'   data = example_data,
#'   var = "variable",
#'   ageband = "ageband",
#'   deprivation = "deprivation",
#'   observed = "observed",
#'   population = "population"
#' )
#'
#' @export

direct_age_dep_standardise <- function(
    data,
    var,
    ageband,
    deprivation,
    observed,
    population,
    multiplier = 1000) {

  # ------------------------------------------------------------
  # Helper to convert quoted or unquoted grouping variables
  # ------------------------------------------------------------

  sym_flex <- function(x) {
    if (is.character(x)) {
      rlang::syms(x)
    } else {
      rlang::syms(as.character(rlang::enexpr(x)))
    }
  }

  var_syms <- sym_flex(var)

  ageband_sym <- rlang::ensym(ageband)
  deprivation_sym <- rlang::ensym(deprivation)
  observed_sym <- rlang::ensym(observed)
  population_sym <- rlang::ensym(population)

  ageband_name <- rlang::as_string(ageband_sym)
  deprivation_name <- rlang::as_string(deprivation_sym)

  # ------------------------------------------------------------
  # Check reference population
  # ------------------------------------------------------------

  if (!all(c("AgeBand", "DepDecile", "StdPopulation") %in%
           names(std_age_dep))) {

    stop(
      "`std_age_dep` must contain columns named `AgeBand`, ",
      "`DepDecile`, and `StdPopulation`.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # Compute observed counts and population by age × deprivation
  # ------------------------------------------------------------

  out <- data %>%
    dplyr::group_by(
      !!!var_syms,
      !!ageband_sym,
      !!deprivation_sym
    ) %>%
    dplyr::summarise(
      obs = sum(!!observed_sym, na.rm = TRUE),
      pop = sum(!!population_sym, na.rm = TRUE),
      .groups = "drop"
    ) %>%

    # ----------------------------------------------------------
  # Ensure all reference age × deprivation combinations exist
  # ----------------------------------------------------------

  dplyr::right_join(
    std_age_dep,
    by = setNames(
      c("AgeBand", "DepDecile"),
      c(ageband_name, deprivation_name)
    )
  ) %>%

    dplyr::mutate(
      obs = dplyr::coalesce(obs, 0),
      pop = dplyr::coalesce(pop, 0)
    ) %>%

    # ----------------------------------------------------------
  # Calculate age × deprivation specific rates
  # ----------------------------------------------------------

  dplyr::mutate(
    rate = dplyr::if_else(
      pop > 0,
      obs / pop,
      0
    )
  ) %>%

    # ----------------------------------------------------------
  # Calculate expected events and variance components
  # ----------------------------------------------------------

  dplyr::mutate(
    expected = rate * StdPopulation,

    var_component = dplyr::if_else(
      pop > 0,
      (StdPopulation^2) * (obs / (pop^2)),
      0
    )
  ) %>%

    # ----------------------------------------------------------
  # Aggregate to requested grouping level
  # ----------------------------------------------------------

  dplyr::group_by(!!!var_syms) %>%
    dplyr::summarise(
      total_observed = sum(obs, na.rm = TRUE),

      total_population = sum(pop, na.rm = TRUE),

      crude_rate = dplyr::if_else(
        sum(pop, na.rm = TRUE) > 0,
        sum(obs, na.rm = TRUE) /
          sum(pop, na.rm = TRUE),
        0
      ),

      dss_rate =
        sum(expected, na.rm = TRUE) /
        sum(StdPopulation, na.rm = TRUE),

      variance =
        sum(var_component, na.rm = TRUE) /
        (sum(StdPopulation, na.rm = TRUE)^2),

      .groups = "drop"
    ) %>%

    # ----------------------------------------------------------
  # Byar 95% confidence intervals
  # ----------------------------------------------------------

  dplyr::mutate(

    lower = ifelse(
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
    lower = pmax(0, lower),

    # --------------------------------------------------------
    # Apply multiplier and rounding
    # --------------------------------------------------------

    crude_rate =
      round(multiplier * crude_rate, 2),

    dss_rate =
      round(multiplier * dss_rate, 2),

    lower =
      round(multiplier * lower, 2),

    upper =
      round(multiplier * upper, 2)
  ) %>%

    # ----------------------------------------------------------
  # Select output
  # ----------------------------------------------------------

  dplyr::select(
    !!!var_syms,
    Observed = total_observed,
    Population = total_population,
    CrudeRate = crude_rate,
    AgeDeprivationStandardisedRate = dss_rate,
    LowerCI = lower,
    UpperCI = upper)

  return(out)
}

