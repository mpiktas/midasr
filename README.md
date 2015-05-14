[![Coverage Status](https://img.shields.io/coveralls/mpiktas/midasr.svg)](https://coveralls.io/r/mpiktas/midasr?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/midasr)](http://cran.r-project.org/web/packages/midasr)
# midasr

The **midasr** R package provides econometric methods for working with mixed frequency data. The package provides tools for estimating time series MIDAS regression, where response and explanatory variables are of different frequency, e.g. quarterly vs monthly. The fitted regression model can be tested for adequacy and then used for forecasting. More specifically, the following main functions are available:

  - ```midas_r``` -- MIDAS regression estimation using NLS
  - ```mls``` -- time series embedding to lower frequency, flexible function for specifying MIDAS models
  - ```hAh.test``` and  ```hAhr.test``` -- adequacy testing of MIDAS regression
  - ```forecast``` -- forecasting MIDAS regression
  - ```midasr_ic_table``` -- lag selection using information criteria
  - ```average_forecast``` -- calculate weighted forecast combination
  - ```select_and_forecast``` -- perform model selection and then use the selected model for forecasting.

The package provides the usual methods for generic functions which can be used on fitted MIDAS regression object: ```summary```, ```coef```, ```residuals```, ```deviance```, ```fitted```, ```predict```, ```logLik```. It also
has additional methods for estimating robust standard errors: ```estfun``` and ```bread```. 

The package also provides all the popular MIDAS regression restrictions such as normalized Almon exponential, normalized beta and etc. 

The package development was influenced by features of the [MIDAS Matlab toolbox][3] created by Eric Ghysels.

The package has the project [webpage][1] and you can follow its development on  [github][2]. 

The detailed description of the package features can be found in the [User guide][4].  All of the code examples in the user guide and some additional examples together with the user guide .Rnw file can be found in the midasr-user-guide github [repository][5]. 

# Development
To install the development version of midasr, it's easiest to use the `devtools` package:

    # install.packages("devtools")
    library(devtools)
    install_github("midasr","mpiktas")

[1]: http://mpiktas.github.com/midasr
[2]: http://github.com/mpiktas/midasr
[3]: http://www.unc.edu/~eghysels/
[4]: https://github.com/mpiktas/midasr-user-guide/raw/master/midasr-user-guide.pdf
[5]: https://github.com/mpiktas/midasr-user-guide/
