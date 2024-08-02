#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN PROJECTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run_projections

# Author: Robyn Forrest
# Date created:  August 2, 2024
# Last Modified: August 2, 2024

# This function handles calling the projection model under the range
#  of future tacs and calls the reference points functions

# arguments:
# posteriors = posterior outputs from the MCMC, includes derived variables

# Returns:
# A list with length ntac (number of projected tacs)
# Each list item contains a dataframe with all the projection variables
#   for each posterior sample
# Note the burn in (warmup) samples have already been removed

# NOTE: as for all the function scripts, this function assumes that dat, ctl and pfc
#  are already in the global space

run_projections <- function(posteriors){


}
