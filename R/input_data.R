#' pcod.dat iscam input file
#'
#' Data file used for input for iscam in the Pcod 2020 assessment
#'
#' @format ## `pcod2020dat`
#' A list with 44 items:
#' \describe{
#'   \item{hasGearNames}{}
#'   \item{gearNames}{Names of the separate gears considered in the model}
#'   \item{hasAgeCompN}{Logical}
#'   \item{narea}{number of areas in the assessment}
#'   \item{ngroup}{number of groups in the assessment}
#'   \item{nsex}{number of sex in the assessment}
#'   \item{syr}{start year for the assessment model}
#'   \item{nyr}{end year for the assessment model}
#'   \item{sage}{start age of assessment}
#'   \item{nage}{maximum age in the assessment}
#'   \item{ngear}{number of gears in the assessment}
#'   \item{alloc}{allocation of catch between gears}
#'   \item{linf}{Von Bertallanffy Linf growth parameter estimate}
#'   \item{k}{Von Bertallanffy k growth parameter estimate}
#'   \item{to}{Von Bertallanffy t0 parameter estimate}
#'   \item{lwscal}{length-weight "a" (scale) parameters}
#'   \item{lwpow}{length-weight "b" (power) parameters}
#'   \item{age50}{}
#'   \item{sd50}{}
#'   \item{usemat}{}
#'   \item{matvec}{}
#'   \item{dd.kage}{}
#'   \item{dd.alpha.g}{}
#'   \item{dd.rho.g}{}
#'   \item{dd.wk}{}
#'   \item{nctobs}{}
#'   \item{catch}{catch data matrix, 65 rows and 7 columns}
#'   \item{nit}{number of survey indexes}
#'   \item{nitnobs}{number of observations for each index of abundance}
#'   \item{survtype}{}
#'   \item{indices}{list of length nit including matrices with index of abundance information}
#'   \item{nagears}{number of gears with age data}
#'   \item{nagearsvec}{}
#'   \item{nagearssage}{start age in catch at age data}
#'   \item{nagearsnage}{start age in catch at age data}
#'   \item{eff}{effective sample size for age composition}
#'   \item{agecompflag}{}
#'   \item{agearsN}{}
#'   \item{nwttab}{}
#'   \item{nwtobs}{}
#'   \item{nmeanwt}{}
#'   \item{nmeanwtobs}{number of observations for mean weight of the catch}
#'   \item{meanwtdata}{mean weight of the catch data}
#'   \item{eof}{end of file check}
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
