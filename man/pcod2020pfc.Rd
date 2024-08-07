% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input_data.R
\docType{data}
\name{pcod2020pfc}
\alias{pcod2020pfc}
\title{pcod.pfc iscam input file for pcod 2020 assessment}
\format{
\subsection{\code{pcod2020pfc}}{

A list with 4 items:
\describe{
\item{num.tac}{Number of fixed TAC values to be used in projections}
\item{tac.vec}{The fixed TAC values to loop over in projections}
\item{num.ctl.options}{Number of control options for projections}
\item{ctl.options}{A vector with 7 control options:
\itemize{
\item \code{syrmeanm} = Start year for calculating average M
\item \code{nyrmeanm} = End year for calculating average M
\item \code{syrmeanrecproj} = Start year for average recruitment period in projections
\item \code{nyrmeanrecproj} = End year for average recruitment period in projections
\item \code{shortcntrlpts} = Use short time series for averaging historical reference points (syr-\code{shortcntrlpts})
\item \code{longcntrlpts} = Use long time series for averaging historical reference points (syr-\code{longcntrlpts})
\item \code{bmin} = Year of bmin for historical Limit Reference Point (LRP), i.e., year of "minimum biomass from which the stock recovered to above average"
}
}
}
}
}
\source{
\url{https://github.com/pbs-assess/pacific-cod-2020}
}
\usage{
pcod2020pfc
}
\description{
projection file used for input for iscam in the Pcod 2020 assessment
}
\keyword{datasets}
