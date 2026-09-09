#' Direct Age Standardisation with Confidence Intervals
#'
#' Computes directly age-standardised rates with 95% confidence intervals
#' using the 2013 European Standard Population and Byar's method.
#'
#' @param data Dataframe with observed values and population.
#' @param var One or more grouping variables (quoted or unquoted).
#' @param ageband Column specifying age bands.
#' @param observed Column with observed counts.
#' @param population Column with denominator population.
#' @param multiplier Scale for rates (default = 1,000).
#'
#' @return A tibble with directly standardised rates and 95% confidence
#'   intervals.
#'
#' @examples
#' example_data <- data.frame(
#'   variable = c("A", "A", "B", "B"),
#'   ageband = c("0-4", "5-9", "0-4", "5-9"),
#'   observed = c(10, 20, 15, 25),
#'   population = c(1000, 1200, 1100, 1300)
#' )
#'
#' direct_age_standardise(
#'   data = example_data,
#'   var = "variable",
#'   ageband = "ageband",
#'   observed = "observed",
#'   population = "population"
#' )
#'
#' @export

direct_age_standardise <- function(
    data,
    var,
    ageband,
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
  observed_sym <- rlang::ensym(observed)
  population_sym <- rlang::ensym(population)

  # ------------------------------------------------------------
  # Calculate observed counts and population by age band
  # ------------------------------------------------------------

  out <- data %>%
    dplyr::group_by(
      !!!var_syms,
      !!ageband_sym
    ) %>%
    dplyr::summarise(
      obs = sum(!!observed_sym, na.rm = TRUE),
      pop = sum(!!population_sym, na.rm = TRUE),
      .groups = "drop"
    ) %>%

    # ----------------------------------------------------------
  # Ensure all reference age bands are present
  # ----------------------------------------------------------

  dplyr::right_join(
    ref_population %>%
      dplyr::select(AgeBand),
    by = setNames(
      "AgeBand",
      rlang::as_string(ageband_sym)
    )
  ) %>%

    dplyr::mutate(
      obs = dplyr::coalesce(obs, 0),
      pop = dplyr::coalesce(pop, 0)
    ) %>%

    # ----------------------------------------------------------
  # Calculate age-specific rates
  # ----------------------------------------------------------

  dplyr::mutate(
    rate = dplyr::if_else(
      pop > 0,
      obs / pop,
      0
    )
  ) %>%

    # ----------------------------------------------------------
  # Join the 2013 European Standard Population
  # ----------------------------------------------------------

  dplyr::left_join(
    ref_population,
    by = setNames(
      "AgeBand",
      rlang::as_string(ageband_sym)
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
  # Aggregate to the requested grouping level
  # ----------------------------------------------------------

  dplyr::group_by(!!!var_syms) %>%
    dplyr::summarise(
      total_observed = sum(obs, na.rm = TRUE),

      total_population = sum(pop, na.rm = TRUE),

      crude_rate = dplyr::if_else(
        sum(pop, na.rm = TRUE) > 0,
        sum(obs, na.rm = TRUE) / sum(pop, na.rm = TRUE),
        0
      ),

      das_rate =
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
      das_rate *
        (
          1 -
            (1 / (9 * total_observed)) -
            (1.96 / (3 * sqrt(total_observed)))
        )^3
    ),

    upper =
      das_rate *
      (
        1 -
          (1 / (9 * (total_observed + 1))) +
          (1.96 / (3 * sqrt(total_observed + 1)))
      )^3,

    # Lower confidence limits cannot be negative
    lower = pmax(0, lower),

    # --------------------------------------------------------
    # Apply multiplier and rounding
    # --------------------------------------------------------

    crude_rate = round(multiplier * crude_rate, 2),

    dss_rate = round(multiplier * das_rate, 2),

    lower = round(multiplier * lower, 2),

    upper = round(multiplier * upper, 2)
  ) %>%

    # ----------------------------------------------------------
  # Select output
  # ----------------------------------------------------------

  dplyr::select(
    !!!var_syms,
    total_observed,
    total_population,
    crude_rate,
    dss_rate,
    lower,
    upper
  )

  return(out)
}
