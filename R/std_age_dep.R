#' Standard Age and Deprivation Reference Population
#'
#' A reference population derived from the 2013 European Standard
#' Population, assigning equal weight to each deprivation decile within
#' each age band. This reference population can be used for direct age
#' and deprivation standardisation.
#'
#' @format A data frame with age- and deprivation-specific standard
#'   population weights:
#' \describe{
#'   \item{AgeBand}{Age band from the European Standard Population}
#'   \item{DepDecile}{Deprivation decile}
#'   \item{StdPopulation}{Standard population weight}
#' }
#'
#' @source Derived from the 2013 European Standard Population using an
#'   equal 10% distribution across deprivation deciles.
#'
#' @keywords datasets
"std_age_dep"
