% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImmatureorMature.R
\name{mature}
\alias{mature}
\title{Assign the sampled individuals as immature or mature according with the maturity stage}
\usage{
mature(dataBio, param.stage_mature)
}
\arguments{
\item{dataBio}{Dataframe with the varable maturity stage included, this should a numeric variable}

\item{param.stage_mature}{This is the first maturity stage of the mature fish}
}
\value{
A dataframe with the variable "mat_stg" included, with only 0s and 1s. The 0 if the individual is immature and as 1 if the individual is mature.
}
\description{
Assign the sampled individuals as immature or mature according with the maturity stage
}
\examples{
mature(dataBio)
}
