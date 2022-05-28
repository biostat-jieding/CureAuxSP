

### This file contains all the documentation of datasets ###


#' toydata1_smc - a toy simulated dataset from the logistic/porprotional hazards mixture cure model
#'
#' @description toydata1 is a toy simulated dataset from additive risk model.
#' The censoring rate is approximately 30%. There are two covariates and the
#' true corresponding covariate effects are 0.5 and 0.8, respectively.
#'
#' @format A data frame with 300 rows and 4 variables, including:
#' \describe{
#' \item{yobs}{the observed failure time}
#' \item{delta}{the censoring indicator}
#' \item{X1}{covariate 1 comes from U(0,1) distribution}
#' \item{X2}{covariate 2 comes from b(1,0.5) distribution}
#' }
#'
#' @examples
#' data("toydata_smc")
#' dim(toydata_smc)
#' summary(toydata)
#' 1-mean(toydata1[,'delta'])  # censoring rate
#'
#'
#' @keywords datasets
#'
"toydata_smc"



