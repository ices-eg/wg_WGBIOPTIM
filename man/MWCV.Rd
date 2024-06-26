% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MWCV.R
\name{MWCV}
\alias{MWCV}
\title{MWCV (mean weigthed CV) calculation}
\usage{
MWCV(df1, variable)
}
\arguments{
\item{df1}{data frame of sampled data by length class in CA format (RDB) (individual measurements)}

\item{variable}{"lenCls"}
}
\value{
MWCV
}
\description{
MWCV (mean weigthed CV) calculation
}
\examples{
MWCV(example_samples,"lenCls")
}
