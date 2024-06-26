% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Maturity_ogive.R
\name{Maturity_ogive}
\alias{Maturity_ogive}
\title{Maturity ogive parameters estimation (L25, L50, L75)}
\usage{
Maturity_ogive(data)
}
\arguments{
\item{data}{Dataframe with information from the simulations, containning thefollowing variables: length, maturity (0- immature; 1  - mature); ID_sim (identification of the number of the simulation run); type (number of individuals selected in the current simulation).}
}
\value{
L25 (length at which 25% of the individuals are mature), L50 (length at which 50% of the individuals are mature), L75 (length at which 75% of the individuals are mature)
}
\description{
Maturity ogive parameters estimation (L25, L50, L75)
}
\examples{
Maturity_ogive(simulation_biodata)

}
