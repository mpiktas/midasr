% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/midas_nlpr_methods.R
\name{plot_midas_coef.midas_nlpr}
\alias{plot_midas_coef.midas_nlpr}
\title{Plot MIDAS coefficients}
\usage{
\method{plot_midas_coef}{midas_nlpr}(
  x,
  term_name = NULL,
  title = NULL,
  compare = NULL,
  normalize = FALSE,
  ...
)
}
\arguments{
\item{x}{\code{midas_r} object}

\item{term_name}{the term name for which the coefficients are plotted. Default is \code{NULL}, which selects the first MIDAS term}

\item{title}{the title string of the graph. The default is \code{NULL} for the default title.}

\item{compare}{the parameters for weight function to compare with the model, default is NULL}

\item{normalize}{logical, if FALSE use the weight from the model, if TRUE, set the normalization coefficient of the weight function to 1.}

\item{...}{not used}
}
\value{
a data frame with restricted MIDAS coefficients, unrestricted MIDAS coefficients and lower and upper confidence interval limits. The data
frame is returned invisibly.
}
\description{
Plots MIDAS coefficients of a MIDAS regression for a selected term.
}
\details{
Plots MIDAS coefficients of a selected MIDAS regression term together with corresponding MIDAS coefficients and their confidence intervals
of unrestricted MIDAS regression
}
\author{
Virmantas Kvedaras, Vaidotas Zemlys
}
