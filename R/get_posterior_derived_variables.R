#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET POSTERIOR DERIVED VARIABLES ARRANGED BY VARIABLE OR BY SAMPLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_posterior_derived_variables

# Called from run-model-v1.R

# Author: Robyn Forrest
# Date created:  August 2, 2024
# Last Modified: August 2, 2024

# This function handles getting out the posterior samples of derived model variables
#  i.e., REPORT() objects from model

#  See readme at https://github.com/kaskr/tmbstan:
#    "What if you want a posterior for derived quantities in the report? Just
#    loop through each posterior sample (row) and call the report function"

# Version 1: posteriors_by_variable
# For reporting, graphs etc
# This is a list where each element is a variable of interest
# Make a list for putting posterior estimates of biomass,
#    recruits and q
# Note that parameters, rec devs and logf are
#    already reported in the mc object
#   (but have added Ft to REPORT to simplify projection model)

# Version 2: posteriors_by_sample
#  for projections
# This is a list where each element is a posterior sample
#  containing all the variables needed for the proj model
#  Gets passed to projection model using purrr::map2_df()

# arguments:
# mc = tmbstan output from fitmcmc in run_model_v1.R
# type = sets whether to return posteriors by variable (type=byvariable)
#  or posteriors by sample (type=bysample)

# Returns:
# if type=byvariable:
#   A list with length 6 (biomass, numbers, recruits, surv, Ft, q)
#   Each list item contains a dataframe with nrow=number of posterior samples (minus warmup)
# if type=bysample:
#   A list with length = number of posterior samples (minus warmup)
#   Each list item contains a list with variables needed for projections
#   Also contains posteriors for the objective function components
# Note the burn in (warmup) samples have already been removed
