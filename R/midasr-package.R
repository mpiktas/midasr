##' Package for estimating and testing MIDAS regression.
##' 
##' The main feature of this package is function \code{\link{hAh.test}}
##' which performs test whether coefficients of MIDAS regression have certain functional form
##' 
##' @name midasr-package
##' @aliases midasr
##' @docType package
##' @title Estimating and testing MIDAS regression
##' @author Virmantas Kvedaras \email{virmantas.kvedaras@@mif.vu.lt}, Vaidotas Zemlys (maintainer) \email{zemlys@@gmail.com}
##' @keywords package
NULL

##' US monthly unemployment rate
##'
##' The monthly unemployment rate for United States from 1948 to 2011.
##' 
##' @name USunempr
##' @docType data
##' @format A \code{\link{ts}} object.
##' @source \href{http://data.bls.gov/timeseries/LNS14000000}{U.S. Bureau of Labor Statistics}
##' @keywords datasets
NULL


##' US annual gross domestic product in billions of chained 2005 dollars
##'
##' The annual gross domestic product in billions of chained 2005 dollars for US from 1948 to 2011.
##' 
##' @name USrealgdp
##' @docType data
##' @format A \code{\link{ts}} object.
##' @source \href{http://www.bea.gov/national/xls/gdplev.xls}{U.S. Department of Commerce, Bureau of Economic Analysis}
##' @keywords datasets
NULL

##' Realized volatility of S&P500 index
##'
##' Realized volatility of S&P500(Live) index of the period 2000 01 03 - 2012 05 22
##'
##' @name rvsp500
##' @docType data
##' @format A \code{data.frame} object with two columns. First column contains date id, and the second the realized volatility for S&P500 index.
##' @source \href{http://realized.oxford-man.ox.ac.uk/media/1366/oxfordmanrealizedvolatilityindices.zip}{Heber, Gerd, Asger Lunde, Neil Shephard and Kevin Sheppard (2009) "Oxford-Man Institute's realized library", Oxford-Man Institute, University of Oxfors}
##' @keywords datasets
NULL
