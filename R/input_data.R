# TODO: Add descriptions of inputs

#' pcod.dat iscam input file, filtered to only include inputs relevant to delay-difference model
#'
#' Data file used for input for iscam in the Pcod 2020 assessment
#'
#' @format ## `pcod2020dat`
#' A list with xx items:
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
#'   \item{catch}{Catch data matrix, nctobs rows, 4 columns: year, gear, type, value. Type: 1 = catch in weight; 2 = catch in numbers}
#'   \item{nit}{Number of survey indexes}
#'   \item{nitnobs}{Number of observations for each index of abundance. A vector of length nit.}
#'   \item{survtype}{Survey type, a vector of length nit. Type: 1 = survey is proportional to numbers, 2 = survey is proportional to biomass, 3 = survey is proportional to spawning biomass (e.g., a spawn survey)}
#'   \item{indices}{List of length nit of matrices nitnobs[i] rows x 5 columns: iyr, it, gear, wt, timing}
#'   \item{nmeanwt}{Number of mean weight series (if none, must be set to 1 - CAN NOW CHANGE THIS REQUIREMENT)}
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
#'   \item{num.params}{number of parameters in the model}
#'   \item{params}{matrix of parameters initial guesses, bounds, estimation phases and prior parameters}
#'   \item{age.size}{matrix of age composition controls -add more description }
#'   \item{sel}{selectivity control matrices, 10 rows and ncols=pcod2020dat$ngear}
#'   \item{start.yr.time.block}{}
#'   \item{num.indices}{number of indexes of abundance}
#'   \item{surv.q}{controls for estimating catachability coefficients}
#'   \item{fit.mean.weight}{}
#'   \item{num.mean.weight.cv}{}
#'   \item{weight.sig}{}
#'   \item{misc}{miscellaneous controls}
#'   \item{eof}{end of file check}
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
