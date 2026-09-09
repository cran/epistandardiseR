#' Standard Age and Sex Reference Population
#'
#' A reference population derived from the European Standard Population,
#' assigning equal weight to females and males within each age band.
#' This 50:50 distribution can be used for direct age and sex standardisation.
#'
#' @format A data frame with age and sex-specific standard population
#'   weights:
#' \describe{
#'   \item{AgeBand}{Age band from the European Standard Population}
#'   \item{Sex}{Sex category}
#'   \item{StdPopulation}{Standard population weight}
#' }
#'
#' @source Derived from the European Standard Population using a 50:50
#'   female/male distribution.
#'
#' @keywords datasets
"std_age_sex"
