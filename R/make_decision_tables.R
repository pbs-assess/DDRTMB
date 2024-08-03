#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAKE DECISION TABLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make_decision_table
# Called from run-model-v1.R

# Author: Robyn Forrest
# Date created:  August 2, 2024
# Last Modified: August 2, 2024

# This function takes the outputs of the projection model and calculates probabilities
# of various metrics for each TAC

# arguments:
# proj_outputs = the list object returned by run_projections()

# Returns:
# A decision table

make_decision_table <- function(proj_outputs){

  dtable <- as.data.frame(matrix(NA,
                              ncol = 6,
                              nrow = length(tac)))
  tac <- pfc$num.tac

  for(t in seq_along(tac)){
    d <- as.data.frame(model$mcmccalcs$proj.dat)
    d <- d[d$TAC == tac[t],]
    dat[t, 1] <- f(tac[t], 0)
    dat[t, 2] <- f(mean(d$B2022B2021 < 1), 2)
    dat[t, 3] <- f(mean(d$F2021F2020 > 1), 2)
    dat[t, 4] <- f(mean(d$B2022Bmin < 1), 2)
    dat[t, 5] <- f(mean(d$B2022BAvgS < 1), 2)
    dat[t, 6] <- f(mean(d$F2021FAvgS > 1), 2)
  }

  col.names <- c("2021 Catch (mt)",
                 "P(B2022 < B2021)",
                 "P(F2021 > F2020)",
                 "P(B2022 < LRP)",
                 "P(B2022 < USR)",
                 "P(F2021 > LRR)")
}

