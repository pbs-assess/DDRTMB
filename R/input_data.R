# TODO: Add descriptions of inputs

#' pcod.dat iscam input file, filtered to only include inputs relevant to delay-difference model
#'
#' Data file used for input for iscam in the Pcod 2020 assessment
#'
#' @format ## `pcod2020dat`
#' A list with 25 items:
#' #' \describe{
#'   \item{gearNames}{Names of the separate gears considered in the model, includes commercial and survey}
#'   \item{syr}{Start year of data}
#'   \item{nyr}{End year of data}
#'   \item{sage}{Start age -- might not use this, just use kage}
#'   \item{nage}{Maximum age}
#'   \item{ngear}{Number of gears, includes commercial and survey}
#'   \item{alloc}{Allocation of catch between gears, CHECK: Catch allocation: 0 < alloc <= 1 if commercial, 0 if survey}
#'   \item{linf}{Growth, von Bertalanffy Linf}
#'   \item{k}{Growth, von Bertalanffy k}
#'   \item{to}{Growth, von Bertalanffy t0}
#'   \item{lwscal}{Length-weight "a" (scale) parameter}
#'   \item{lwpow}{Length-weight "b" (power) parameter}
#'   \item{kage}{Knife-edge age-at-recruitment -- should be same as sage}
#'   \item{alpha.g}{Growth, Ford-Walford alpha}
#'   \item{rho.g}{Growth, Ford-Walford rho}
#'   \item{wk}{Weight at age of recruitment}
#'   \item{nctobs}{Number of catch observations}
#'   \item{catch}{Catch data matrix, nctobs rows, 4 columns: year, gear, type, value (tonnes). Type: 1 = catch in weight; 2 = catch in numbers}
#'   \item{nit}{Number of survey indexes}
#'   \item{nitnobs}{Number of observations for each index of abundance. A vector of length nit.}
#'   \item{survtype}{Survey type, a vector of length nit. Type: 1 = survey is proportional to numbers, 2 = survey is proportional to biomass, 3 = survey is proportional to spawning biomass (e.g., a spawn survey)}
#'   \item{indices}{List of length nit of matrices nitnobs[i] rows x 5 columns: iyr = Survey year, omits missing years, it = survey value, gear = gear index from ngear (includes commercial and surveys), wt = 1/CV (used as precision to multiplicatively weight observations), timing = the fraction of total mortality that has occurred prior to survey. Usually zero.}
#'   \item{nmeanwt}{Number of mean weight series}
#'   \item{nmeanwtobs}{Number of observations for mean weight of the catch}
#'   \item{meanwtdata}{Annual commercial mean weights (kg)}
#' }
#' @source <https://github.com/pbs-assess/pacific-cod-2020>
"pcod2020dat"

#' pcod.ctl iscam input file for pcod 2020 assessment
#'
#' Control file used for input for iscam in the Pcod 2020 assessment
#'
#' @format ## `pcod2020ctl`
#' A list with 8 items:
#' \describe{
#'   \item{num.params}{Number of leading parameters in the model. Does not include catchability parameters (q)}
#'   \item{params}{A matrix $num.params rows x 7 columns.
#'
#'    Parameters (rows):
#'    \itemize{
#'      \item \code{log_ro} log unfished recruitment.
#'      \item \code{steepness} steepness.
#'      \item \code{log_m} log natural mortality.
#'      \item \code{log_avrec} log average recruitment.
#'      \item \code{log_recinit} log average of initial recruitments to fill first year.
#'      \item \code{rho} EIV: fraction of the total variance associated with observation error.
#'      \item \code{kappa (varphi)} EIV: total precision (inverse of variance) of the total error.
#'    }
#'     Settings (columns):
#'    \itemize{
#'      \item \code{ival} starting value
#'      \item \code{lb} lower bound (not sure if this can be used in tmb)
#'      \item \code{ub} upper bound (not sure if this can be used in tmb)
#'      \item \code{phz} phase of estimation (in ADMB). Won't need this in tmb but need to check how to fix a par
#'      \item \code{prior} prior type:
#'     \itemize{
#'        \item 0 = uniform      (0,0),
#'        \item 1 =  normal       (p1=mu,p2=sig),
#'        \item 2 =  lognormal    (p1=log(mu),p2=sig),
#'        \item 3 =  beta         (p1=alpha,p2=beta),
#'        \item 4 =  gamma        (p1=alpha,p2=beta)
#'    }
#'      \item \code{p1} prior distribution parameters, see \code{prior}
#'      \item \code{p2} prior distribution parameters, see \code{prior}}
#'    }
#'   \item{num.indices}{Number of indexes of abundance (surveys only). DEPRECATE THIS AND USE DAT$NIT}
#'   \item{surv.q}{Controls for estimating catachability coefficients. Matrix 3 x num.indices (i.e., 1 column per index).
#'      Settings (rows):
#'   \itemize{
#'    \item 1. Prior type:
#'   \itemize{
#'      \item  0 = uniform
#'      \item  1 = normal on log q
#'      \item  2 = random walk in q
#'    }
#'    \item 2. Prior mean
#'    \item 3. Prior sd}
#'    }
#'   \item{fit.mean.weight}{Switch for fitting to annual mean weights. 0=FALSE, 1=TRUE.}
#'   \item{num.mean.weight}{Number of mean weight series}
#'   \item{weight.sig}{SD for mean weight likelihood}
#'   \item{misc}{13 miscellaneous controls:
#'   \itemize{
#'    \item    1.  verbose ADMB output (0=off, 1=on)
#'    \item    2.  recruitment model (1=beverton-holt, 2=ricker)
#'    \item    3.  std in observed catches in first phase.
#'    \item    4.  std in observed catches in last phase.
#'    \item    5.  assume unfished equilibrium in first year (0=FALSE, 1=TRUE, 2 = At equilibrium with F in first year)
#'    \item    6.  Maternal effects multiplier
#'    \item    7.  Mean fishing mortality for regularizing the estimates of Ft
#'    \item    8.  std in mean fishing mortality in first phase
#'    \item    9.  std in mean fishing mortality in last phase
#'    \item    10. phase for estimating m_deviations (use -1 to turn off mdevs)
#'    \item    11. std in deviations for natural mortality
#'    \item    12. number of estimated nodes for deviations in natural mortality
#'    \item    13. fraction of total mortality that takes place prior to spawning}
#' }
#' }
#' @source <https://github.com/pbs-assess/pacific-cod-2020>
"pcod2020ctl"







#' pcod.pfc iscam input file for pcod 2020 assessment
#'
#' projection file used for input for iscam in the Pcod 2020 assessment
#'
#' @format ## `pcod2020pfc`
#' A list with 5 items:
#' \describe{
#'   \item{num.tac}{}
#'   \item{tac.vec}{}
#'   \item{num.ctl.options}{}
#'   \item{ctl.options}{}
#'   \item{eof}{end of file check}
#' }
#' @source <https://github.com/pbs-assess/pacific-cod-2020>
"pcod2020pfc"
