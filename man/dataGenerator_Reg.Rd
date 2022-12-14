% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataGenerator.R
\name{dataGenerator_Reg}
\alias{dataGenerator_Reg}
\title{dataGenerator_Reg}
\usage{
dataGenerator_Reg(chpts = 100, A = 2, B = -1, meansX = 0, sdX = 1, sdNoise = 1)
}
\arguments{
\item{chpts}{a vector of increasing change-point indices}

\item{A}{vector of regression coefficients A in A*x+B simple regression model}

\item{B}{vector of regression coefficients B in A*x+B simple regression model}

\item{meansX}{vector of mean values for x values generated by a Gaussian model}

\item{sdX}{vector of standard deviation values for x values generated by a Gaussian model}

\item{sdNoise}{standard deviation of the Gaussian noise}
}
\value{
a dataframe with time series x and y for change-points in regression of type y = A*x + B + noise
}
\description{
Generating data for changes in bivariate independent time series
}
\examples{
myData <- dataGenerator_Reg(chpts = c(40,90), A = c(2,-1),  B = c(-1,2), meansX = c(1,2))
}
