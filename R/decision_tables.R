decision.table <- function(proj_outputs){

  dtable <- as.data.frame(matrix(NA,
                              ncol = 6,
                              nrow = length(tac)))
  col.names <- c("2021 Catch (mt)",
                   "P(B2022 < B2021)",
                   "P(F2021 > F2020)",
                   "P(B2022 < LRP)",
                   "P(B2022 < USR)",
                   "P(F2021 > LRR)")

  tac <- model$proj$tac.vec
  if(!is.na(tac.vec[1])){
    tac <- tac.vec[tac.vec %in% tac]
  }
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

  if(make.lt.gt){
    dat <- mutate_at(dat, -1,
                     function(x) gsub('0.00', '<0.01', x))
    dat <- mutate_at(dat, -1,
                     function(x) gsub('1.00', '>0.99', x))
  }



}

