# reginv
 Regression Inversion for Parametric Inference

Tools to use inversion for accurate small-sample inference about a target parameter.
`reginv` is a generic function which, given functions defining a simulation model and how to construct a statistic from simulated data, will use inversion to estimate a confidence interval for a parameter of interest using quantile regression or a related technique.

The main application to date is estimating extinction time from the fossil record using the `cutt` distribution, as in `reginv_cutt`.

### Installation
Install the **development** version of `reginv`
from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("eco-stats/reginv")

library(reginv)
```

### Spot a bug?

Thanks for finding the bug! We would appreciate it if you can pop over
to our [Issues page](https://github.com/eco-stats/reginv/issues) and
describe how to reproduce the bug!
