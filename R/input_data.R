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
#' A list with 12 items:
#' \describe{
#'   \item{num.params}{Number of leading parameters in the model. Does not include catchability parameters (q)}
#'   \item{params}{A matrix $num.params rows x 7 columns.
#'          Parameters (rows): 1. log_ro = log unfished recruitment,
#'          2. steepness = steepness,
#'          3. log_m = log natural mortality,
#'          4. log_avrec = log average recruitment,
#'          5. log_recinit = log average of initial recruitments to fill first year,
#'          6. rho = EIV: fraction of the total variance associated with observation error,
#'          7. kappa (varphi) = EIV: total precision (inverse of variance) of the total error.
#'          Settings (columns):
#'          1. ival = starting value
#'          2. lb = lower bound (not sure if this can be used in tmb)
#'          3. ub = upper bound (not sure if this can be used in tmb)
#'          4. phz = phase of estimation (in ADMB). Won't need this in tmb but need to check how to fix a par
#'          5. prior = prior type:
#'           0 uniform      (0,0)
#'           1 normal       (p1=mu,p2=sig)
#'           2 lognormal    (p1=log(mu),p2=sig)
#'           3 beta         (p1=alpha,p2=beta)
#'           4 gamma        (p1=alpha,p2=beta)
#'          6. p1 = prior distribution parameters, see above
#'          7. p2 = prior distribution parameters, see above}
#'   \item{num.indices}{Number of indexes of abundance (surveys only). DEPRECATE THIS AND USE DAT$NIT}
#'   \item{surv.q}{Controls for estimating catachability coefficients. Matrix 3 x num.indices (i.e., 1 column per index).
#'           Settings (rows):
#'           1. Prior type:
#'            0 uniform
#'            1 normal on log q
#'            2 random walk in q
#'          2. Prior mean
#'          3. Prior sd}
#'   \item{fit.mean.weight}{Switch for fitting to annual mean weights. 0=FALSE, 1=TRUE.}
#'   \item{num.mean.weight}{Number of mean weight series}
#'   \item{weight.sig}{SD for mean weight likelihood}
#'   \item{misc}{13 miscellaneous controls:
#'         1  -verbose ADMB output (0=off, 1=on)
#'        2  -recruitment model (1=beverton-holt, 2=ricker)
#'        3  -std in observed catches in first phase.
#'        4  -std in observed catches in last phase.
#'        5  -Assume unfished equilibrium in first year (0=FALSE, 1=TRUE, 2 = AT EQUILIBRIUM WITH FISHING MORTALITY IN SYR - IMPLEMENTED ONLY IN DELAY DIFF MODEL)
#'        6  -Maternal effects multiplier
#'        7  -Mean fishing mortality for regularizing the estimates of Ft
#'        8  -std in mean fishing mortality in first phase
#'        9  -std in mean fishing mortality in last phase
#'        10 -phase for estimating m_deviations (use -1 to turn off mdevs)
#'        11 -std in deviations for natural mortality
#'        12 -number of estimated nodes for deviations in natural mortality
#'        13 -fraction of total mortality that takes place prior to spawning}
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
