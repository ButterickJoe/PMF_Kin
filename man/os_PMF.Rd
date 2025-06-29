% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/os_PMF.R
\name{os_PMF}
\alias{os_PMF}
\title{Title Unconditional PMF for Older sisters}
\usage{
os_PMF(y, s1, U_mat, F_mat, Q)
}
\arguments{
\item{y}{Scalar. age of Focal at present}

\item{s1}{Scalar. age of Focal's older sister at present}

\item{U_mat}{Matrix. probabilities of surviving from age class i to i+1 in the i-th entry of sub-diagonal}

\item{F_mat}{Matrix. fertility rate of age class i in the i-th entry of top-row}

\item{Q}{Scalar. Life-time maximum kin-number}
}
\value{
Vector. PMF for the number-distribution of older sisters of age s1 when Focal is y
}
\description{
Title Unconditional PMF for Older sisters
}
