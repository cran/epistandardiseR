#' Direct Deprivation Standardisation with Confidence Intervals
#'
#' Computes directly deprivation-standardised rates with 95% confidence
#' intervals using a standard national deprivation profile with equal
#' population in each deprivation decile and Byar’s method.
#'
#' @param data Data frame containing observed counts and population denominators
#' @param var One or more grouping variables (quoted or unquoted)
#' @param dep_decile Column specifying deprivation deciles (1–10)
#' @param observed Column with observed counts
#' @param population Column with denominator population
#' @param multiplier Scale for rates (default = 1,000)
#'
#' @return A tibble with deprivation-standardised rates and 95% confidence intervals
#'
#' @examples
#' example_data <- data.frame(
#'   group = c("A", "A", "B", "B"),
#'   dep_decile = factor(c(1, 2, 1, 2)),
#'   observed = c(10, 20, 15, 25),
#'   population = c(1000, 1200, 1100, 1300)
#' )
#'
#' direct_dep_standardise(
#'   data = example_data,
#'   var = "group",
#'   dep_decile = "dep_decile",
#'   observed = "observed",
#'   population = "population"
#' )
#'
#' @export
direct_dep_standardise <- function(
    data,
    var,
    dep_decile,
    observed,
    population,
    multiplier = 1000
) {

  # Helper to convert quoted or unquoted names to symbols
  sym_flex <- function(x) {
    if (is.character(x)) rlang::syms(x)
    else rlang::syms(as.character(rlang::enexpr(x)))
  }

  var_syms <- sym_flex(var)
  dep_decile_sym <- rlang::ensym(dep_decile)
  observed_sym <- rlang::ensym(observed)
  population_sym <- rlang::ensym(population)

  # Convert deprivation decile to factor if numeric
  if (is.numeric(data[[rlang::as_string(dep_decile_sym)]])) {
    data[[rlang::as_string(dep_decile_sym)]] <- factor(
      data[[rlang::as_string(dep_decile_sym)]],
      levels = sort(unique(data[[rlang::as_string(dep_decile_sym)]]))
    )
  }

  out <- data %>%
    dplyr::group_by(!!!var_syms, !!dep_decile_sym) %>%
    dplyr::summarise(
      obs = sum(!!observed_sym, na.rm = TRUE),
      pop = sum(!!population_sym, na.rm = TRUE),
      .groups = "drop"
    ) %>%

    # ---- Ensure all deprivation deciles are present ----
  dplyr::right_join(
    std_dep_decile %>% dplyr::select(dep_decile),
    by = setNames("dep_decile", rlang::as_string(dep_decile_sym))
  ) %>%
    dplyr::mutate(
      obs = tidyr::replace_na(obs, 0),
      pop = tidyr::replace_na(pop, 0)
    ) %>%

    # ---- Safe rate calculation ----
  dplyr::mutate(
    rate = dplyr::if_else(pop > 0, obs / pop, 0)
  ) %>%

    # ---- Join standard deprivation population ----
  dplyr::left_join(
    std_dep_decile,
    by = setNames("dep_decile", rlang::as_string(dep_decile_sym))
  ) %>%

    dplyr::mutate(
      expected = rate * std_population,
      var_component = (std_population^2) * (obs / (pop^2))
    ) %>%

    dplyr::group_by(!!!var_syms) %>%
    dplyr::summarise(
      total_observed = sum(obs, na.rm = TRUE),
      total_population = sum(pop, na.rm = TRUE),
      crude_rate = sum(obs, na.rm = TRUE) / sum(pop, na.rm = TRUE),
      dds_rate = sum(expected, na.rm = TRUE) / sum(std_population, na.rm = TRUE),
      variance = sum(var_component, na.rm = TRUE) /
        (sum(std_population, na.rm = TRUE)^2),
      .groups = "drop"
    ) %>%

    dplyr::mutate(
      # ---- Byar confidence intervals ----
      lower = ifelse(
        total_observed == 0, 0,
        dds_rate * ((1 - (1 / (9 * total_observed)) -
                       (1.96 / (3 * sqrt(total_observed))))^3)
      ),
      upper = dds_rate *
        ((1 - (1 / (9 * (total_observed + 1))) +
            (1.96 / (3 * sqrt(total_observed + 1))))^3),
      lower = pmax(0, lower),

      # ---- Apply multiplier and rounding ----
      crude_rate = round(multiplier * crude_rate, 2),
      dds_rate = round(multiplier * dds_rate, 2),
      lower = round(multiplier * lower, 2),
      upper = round(multiplier * upper, 2)
    ) %>%

    dplyr::select(
      !!!var_syms,
      total_observed,
      total_population,
      crude_rate,
      dds_rate,
      lower,
      upper
    )

  return(out)
}
