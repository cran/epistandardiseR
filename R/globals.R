# Declare global variables used inside dplyr pipelines
utils::globalVariables(c(
  # direct age standardisation
  "obs",
  "pop",
  "rate",
  "expected",
  "StdPopulation",
  "var_component",
  "das_rate",
  "dss_rate",
  "lower",
  "upper",
  "total_obs",
  "ref_population",
  "std_age_sex",
  "std_age_dep",

  # additional DSR variables
  "AgeBand",
  "total_observed",
  "crude_rate",
  "total_population",
  "std_sex",

  # deprivation standardisation
  "std_dep_decile",
  "std_population",
  "dds_rate"
))
