
#==========================================================================#
# Calculate Subgroup survival rates using KM
#==========================================================================#

#' @title Subgroup survival rates
#'
#' @description Calculate Subgroup survival rates basedon the Kaplan-Meier estimation procedure
#'
#' @aliases St.Sub.KM
#'
#' @param tstar time points that the survival rates will be estimated at.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param G a matrix used to indicate which subgroups he/she belongs to for each of these subjects.
#'
#' @export St.Sub.KM
St.Sub.KM <- function(tstar,yobs,delta,G){

  K <- nrow(G)
  tSP <- rep(0, K)
  for(k in 1:K){
    idx <- (G[k,]==1)
    fit.k <- summary(survival::survfit(survival::Surv(yobs, delta) ~ 1, subset=idx))
    tSP[k] <- min( c(1,fit.k$surv)[ c(0,fit.k$time) <= tstar ] )
  }
  return( tSP )
}

#==== Sub-function: disturbing function used in the bootstrap procedure ====# (insert the elements of this function into the main body)
#' @title Disturbing function used in the bootstrap procedure
#'
#' @description Disturb the auxiliary in the bootstrap procedure
#' @aliases AuxSP.disturb.KM
#'
#' @param nboot specifies the number of bootstrap sampling. The default \code{nboot = 100}.
#' @param V.KM the estimated variance-covariance matrix of the auxiliary information based on Kaplan-Meier estimator
#' @param auxinfo a matrix that collects the auxiliary information in a pre-specified way.
#'
#' @export AuxSP.disturb.KM
AuxSP.disturb.KM <- function(nboot,V.KM,auxinfo){

  ScaleP <- diag(1/sqrt(auxinfo[,'M']))
  disturbs <- mvtnorm::rmvnorm(nboot,mean=rep(0,nrow(auxinfo)),
                               sigma=ScaleP %*% V.KM %*% ScaleP)
  return(disturbs)

}






#==========================================================================#
# Semi-parametric mixture cure model with auxiliary information - using CV (Pen/MA) ####
#==========================================================================#


#==== The function for fitting mixture cure model with homogeneous auxiliary information (t*-year survival probability) ====#
#' @title Semiparametric mixture cure model with auxiliary subgroup survival information using AIS method
#'
#' @description Fit the semi-parametric mixture cure model with auxiliary subgroup survival probability information using AIS method.
#' The test statistic for evaluating the homogeneity assumption will also be calculated.
#' @aliases SMC.AuxSP.AIS
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}.
#'   The \code{response} is a \code{Surv} object with right censoring.
#'   It is used to specify the covariate (risk factor) effects on the failure time of uncured subjects.
#'   See the documentation for \code{survreg} and \code{Surv} in R package \code{survival} for details.
#'   The expression to the right of the "~" specifies the effect of covariates on the failure time of uncured patients.
#' @param cureform indicator function a formula expression, of the form \code{cureform ~ predictors}.
#'   It is used to specify the effects of covariates on the cure rate.
#'   A covariate may be used in both \code{formula} and \code{cureform}.
#' @param data a data frame in which to interpret the variables named in the \code{formula} and the \code{cureform}.
#' @param aux indicates the historical summary data. It is a list of lists, and each sub-list represents auxiliary information from a study.
#'   In each study, we have several time points and each time point contains the following four elements
#'   \code{M} the sample size of the external study;
#'   \code{tstar} the time point that the auxiliary information was calculated at;
#'   \code{gfun} a function used to identify the subgroup;
#'   \code{sprob} auxiliary subgroup survival rates for each subgroup at the current time point.
#' @param latency specifies the model used in latency part.
#'   It can be \code{ph} which represents the proportional hazards model, or \code{aft} which represents the accelerated failure time model.
#' @param stdz a logical value. If it is \code{TRUE}, all the covariates are standardized. Otherwise, non-transformed covariates are used.
#' @param nboot specifies the number of bootstrap sampling. The default \code{nboot = 100}.
#'
#' @details
#'    This is a function used to fit the semiparametric mixture cure model with auxilary subgrop survival information.
#'    The method used here is the control variate technique.
#'    The test statistic for evaluating the homogeneity assumption will also be calculated.
#'
#' @return
#'   An object of class \code{SMC.AuxSP} is returned. It can be examined by \code{print.SMC.AuxSP()}.
#'   Note that the fitting results of corresponding mixture cure model without taking auxiliary information into consideration will be contained in the returned object.
#'
#' @examples
#' ### Example 1. Fit the PHMC model with auxiliary information for the simulated toy data.
#'
# library(CureAuxSP)
# data("toydata_smc")
#' set.seed(1)
#'
#' ## split the dataset into two parts: one for internal study and the other for external study
#' sdata.internal.select <- c(1:5000) %in% sample(1:5000,500,replace=FALSE)
#' sdata.internal <- toydata_smc[sdata.internal.select,]
#' sdata.external <- toydata_smc[!sdata.internal.select,]
#'
#' ## prepare our auxiliary information: derived from the external study defined above
#' gfunc = function(X,Z){
#'   rbind(  ( X[,1] >= -1.0 & X[,1] < -0.5 & X[,2] == 0 ),
#'           ( X[,1] >= -0.5 & X[,1] <  0.0 & X[,2] == 0 ),
#'           ( X[,1] >=  0.0 & X[,1] <  0.5 & X[,2] == 0 ),
#'           ( X[,1] >=  0.5 & X[,1] <  1.0 & X[,2] == 0 ),
#'           ( X[,1] >= -1.0 & X[,1] < -0.5 & X[,2] == 1 ),
#'           ( X[,1] >= -0.5 & X[,1] <  0.0 & X[,2] == 1 ),
#'           ( X[,1] >=  0.0 & X[,1] <  0.5 & X[,2] == 1 ),
#'           ( X[,1] >=  0.5 & X[,1] <  1.0 & X[,2] == 1 )
#'   ) * 1
#' }
#' sprob <- St.Sub.KM(
#'   tstar=1,yobs=sdata.external[,1],delta=sdata.external[,2],
#'   G=gfunc(sdata.external[,-c(1,2)],sdata.external[,-c(1,2)])
#' )
#' aux.homo <- list(
#'   study1 = list(
#'     time1 = list(
#'       M = 4500, tstar = 1,
#'       gfunc = gfunc, sprob = sprob
#'     )))
#' aux.hetero <- list(
#'   study1 = list(
#'     time1 = list(
#'       M = 4500, tstar = 1,
#'       gfunc = gfunc, sprob = sprob + c(0.2,-0.1,0,0,0,0,0,0)
#'     )))
#'
#' ## fit the model with homogeneous auxiliary information
#' toydata.SMC.AuxSP.ph <-  SMC.AuxSP.AIS(
#'   formula = Surv(yobs,delta)~X1+X2, cureform = ~X1, data=sdata.internal,
#'   aux = aux.homo, latency = "ph", stdz=FALSE, nboot = 100
#' )
#' print.SMC.AuxSP(toydata.SMC.AuxSP.ph)
#'
#' ## fit the model with heterogeneous auxiliary information
#' # without adjusting the heterogeneity (AIS)
#' toydata.SMC.AuxSP.ph.homo <-  SMC.AuxSP.AIS(
#'   formula = Surv(yobs,delta)~X1+X2, cureform = ~X1, data=sdata.internal,
#'   aux = aux.hetero, latency = "ph", stdz=FALSE, nboot = 100
#' )
#' print.SMC.AuxSP(toydata.SMC.AuxSP.ph.homo)
#' # use penalization techniques (PAIS)
#' toydata.SMC.AuxSP.ph.hetero_pen <-  SMC.AuxSP.PAIS(
#'   formula = Surv(yobs,delta)~X1+X2, cureform = ~X1, data=sdata.internal,
#'   aux = aux.hetero, latency = "ph", stdz=FALSE, nboot = 100
#' )
#' print.SMC.AuxSP(toydata.SMC.AuxSP.ph.hetero_pen)
#'
#' @export SMC.AuxSP.AIS
SMC.AuxSP.AIS <- function(formula, cureform, data, aux, latency = c("ph","aft")[1],
                          stdz=FALSE, nboot = 100
){

  ## basic (whether these inputs are valid enough to make the function works well)
  call <- match.call()
  Y <- model.response(model.frame(formula,data))
  if (!inherits(Y, "Surv")){stop("Response must be a survival object")}

  ## prepare - necesary vectors and matrices + extract vars from formulas
  mf.latency <- model.frame(formula,data)
  mf.incidence <- model.frame(cureform,data)
  X <-  model.matrix(attr(mf.latency, "terms"), mf.latency)[,-1,drop=F]
  Z <-  model.matrix(attr(mf.incidence, "terms"), mf.incidence)[,-1,drop=F]
  vars.name.response <- all.vars(formula)[c(1,2)]
  vars.name.latency <- sapply(colnames(X),function(x){for(xx in c(" ","(",")",":")){x <- sub(xx,"_",x,fixed=T)};x},USE.NAMES=F)
  vars.name.incidence <- sapply(colnames(Z),function(x){for(xx in c(" ","(",")",":")){x <- sub(xx,"_",x,fixed=T)};x},USE.NAMES=F)
  vars.name.covariate <- union(vars.name.latency,vars.name.incidence)
  colnames(X) <- vars.name.latency; colnames(Z) <- vars.name.incidence
  yobs  <- data[,vars.name.response[1]]
  delta <- data[,vars.name.response[2]]

  ## fit the model under homogeneous asumption
  res <- SMC.AuxSP.AIS.fit(
    yobs=yobs,delta=delta,X=X,Z=Z,aux=aux,latency=latency,stdz=stdz,
    methods.suffix=c("ori","homo"),tuning=NULL,auxsort=NULL,nboot=nboot
  )

  ## tidy the results
  fit <- list()
  fit$call <- call
  fit$info <- list(
    hetero=FALSE,
    latency=latency
  )
  fit$res.noaux <- res$res.ori
  fit$test.chisq <- res$res.homo$test.chisq
  fit$incidence <- res$res.homo$incidence
  fit$latency   <- res$res.homo$latency
  class(fit) <- c("SMC.AuxSP")

  ## return
  return(fit)

}


#==== The function for fitting mixture cure model with heterogeneous auxiliary information (t*-year survival probability) ====#
#' @title Semiparametric mixture cure model with auxiliary subgroup survival information using PAIS method
#'
#' @description Fit the semi-parametric mixture cure model with heterogeneous auxiliary subgroup survival probability information using PAIS method.
#' The test statistic for evaluating the homogeneity assumption will also be calculated.
#'
#' @aliases SMC.AuxSP.PAIS
#'
#' @inherit SMC.AuxSP.AIS
#'
#' @details
#'    This is a function used to fit the semiparametric mixture cure model with auxilary subgrop survival information.
#'    The method used here is the control variate technique.
#'    The test statistic for evaluating the homogeneity assumption will also be calculated.
#'
#' @return
#'   An object of class \code{SMC.AuxSP} is returned. It can be examined by \code{print.SMC.AuxSP()}.
#'   Note that the fitting results of corresponding mixture cure model without taking auxiliary information into consideration will be contained in the returned object.
#'
#' @examples
#'   See examples presented in \code{?SMC.AuxSP.AIS}
#'
#' @export SMC.AuxSP.PAIS
SMC.AuxSP.PAIS <- function(formula, cureform, data, aux=NULL, latency = c("ph","aft")[1],
                           stdz=FALSE, nboot = 100
){

  ## basic (whether these inputs are valid enough to make the function works well)
  call <- match.call()
  Y <- model.response(model.frame(formula,data))
  if (!inherits(Y, "Surv")){stop("Response must be a survival object")}

  ## prepare - necesary vectors and matrices + extract vars from formulas
  mf.latency <- model.frame(formula,data)
  mf.incidence <- model.frame(cureform,data)
  X <-  model.matrix(attr(mf.latency, "terms"), mf.latency)[,-1,drop=F]
  Z <-  model.matrix(attr(mf.incidence, "terms"), mf.incidence)[,-1,drop=F]
  vars.name.response <- all.vars(formula)[c(1,2)]
  vars.name.latency <- sapply(colnames(X),function(x){for(xx in c(" ","(",")",":")){x <- sub(xx,"_",x,fixed=T)};x},USE.NAMES=F)
  vars.name.incidence <- sapply(colnames(Z),function(x){for(xx in c(" ","(",")",":")){x <- sub(xx,"_",x,fixed=T)};x},USE.NAMES=F)
  vars.name.covariate <- union(vars.name.latency,vars.name.incidence)
  colnames(X) <- vars.name.latency; colnames(Z) <- vars.name.incidence
  yobs  <- data[,vars.name.response[1]]
  delta <- data[,vars.name.response[2]]

  ## fit the model under homogeneous asumption
  res <- SMC.AuxSP.AIS.fit(
    yobs=yobs,delta=delta,X=X,Z=Z,aux=aux,latency=latency,stdz=stdz,
    methods.suffix=c("ori","homo","hetero_pen"),tuning=NULL,auxsort=NULL,nboot=nboot
  )

  ## tidy the results
  fit <- list()
  fit$call <- call
  fit$info <- list(
    hetero=TRUE,
    method="PAIS",
    latency=latency
  )
  fit$res.noaux <- res$res.ori
  fit$test.chisq <- res$res.homo$test.chisq
  fit$incidence  <- res$res.hetero_pen$incidence
  fit$latency    <- res$res.hetero_pen$latency
  fit$tau        <- res$res.hetero_pen$tau

  ## return
  return(fit)

}




#==== The function for fitting mixture cure model with heterogeneous auxiliary information (t*-year survival probability) ====#
#' @title Semiparametric mixture cure model with auxiliary subgroup survival information using MAIS method
#'
#' @description Fit the semi-parametric mixture cure model with heterogeneous auxiliary subgroup survival probability information using PAIS method.
#' The test statistic for evaluating the homogeneity assumption will also be calculated.
#'
#' @aliases SMC.AuxSP.MAIS
#'
#' @inherit SMC.AuxSP.AIS
#' @param method specify the method that will used to accomodate heterogeneous auxiliary information.
#'   The available methods include:
#'    \code{MAIS} indicates the model averaging auxiliary information synthesis (MAIS) method (with all possible candidate models);
#'    \code{MAIS_Nest} indicates the model averaging auxiliary information synthesis (MAIS) method (with nested candidate models sorted by PAIS).
#' @details
#'    This is a function used to fit the semiparametric mixture cure model with auxilary subgrop survival information.
#'    The method used here is the control variate technique.
#'    The test statistic for evaluating the homogeneity assumption will also be calculated.
#'
#' @return
#'   An object of class \code{SMC.AuxSP} is returned. It can be examined by \code{print.SMC.AuxSP()}.
#'   Note that the fitting results of corresponding mixture cure model without taking auxiliary information into consideration will be contained in the returned object.
#'
#' @examples
#'   See examples presented in \code{?SMC.AuxSP.AIS}
#'
#' @export SMC.AuxSP.MAIS
SMC.AuxSP.MAIS <- function(formula, cureform, data, aux=NULL, latency = c("ph","aft")[1],
                           tuning = NULL, auxsort = NULL,stdz = FALSE, nboot = 100
){

  # This function has not been finished yet!

  if(fit$method=="MA1"){
    cat("\nEstimated Nuiance Parameters Tau:\n")
    print(fit$tau)
    cat("\nCandidate Models With Highest Weights:\n")
    printnum <- min(max(1,round((ncol(fit$candidates)-1)/2)),nrow(fit$candidates))
    print(fit$candidates[order(fit$candidates[,1],decreasing = T)[1:printnum],])
  }else if(fit$method=="MA2"){
    cat("\nEstimated Nuiance Parameters Tau:\n")
    print(fit$tau)
    cat("\nCandidate Models With Highest Weights:\n")
    printnum <- min(max(1,round((ncol(fit$candidates)-1)/2)),nrow(fit$candidates))
    print(fit$candidates[order(fit$candidates[,1],decreasing = T)[1:printnum],])
  }
}


#==== The summary function for fitting mixture cure model with auxiliary information ====#
#' @title Print SMC.AuxSP object
#'
#' @description Output of \code{SMC.AuxSP} object.
#'
#' @aliases print.SMC.AuxSP
#'
#' @param fit an object of SMC.AuxSP
#'
#' @export print.SMC.AuxSP
print.SMC.AuxSP <- function(fit){

  ## print information for model fitting
  if (!is.null(cl <- fit$call)) {
    cat("Call:\n")
    dput(cl)
  }

  ## print the fitting results for regression coefficients
  cat("\nCure Probability Model:\n")
  print(fit$incidence)
  cat("\nFailure Time Distribution Model:\n")
  print(fit$latency)

  ## print test statistics for evaluating homgeneity
  cat("\nTest Statistics For Evaluating Homogeneity:\n")
  print(fit$test.chisq)

  ## If PAIS, print etimated Nuisance Parameters
  if(fit$info$hetero==TRUE){
    cat("\nEstimated Nuiance Parameters Tau:\n")
    print(fit$tau)
  }

  ## finish
  invisible(fit)

}






#==== The main function for fitting the model ====#
#' @title Semiparametric mixture cure model with auxilary subgrop survival information
#'
#' @description Fit the semiparametric mixture cure model with auxilary subgrop survival information.
#'
#' @aliases SMC.AuxSP.AIS.fit
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param Z a matrix of covariates corresponding to incidence part.
#' @param aux indicates the historical summary data. It is a list of lists, and each sub-list represents auxiliary information from a study.
#'   In each study, we have several time points and each time point contains the following four elements
#'   \code{M} the sample size of the external study;
#'   \code{tstar} the time point that the auxiliary information was calculated at;
#'   \code{gfun} a function used to identify the subgroup;
#'   \code{sprob} auxiliary subgroup survival rates for each subgroup at the current time point.
#' @param latency specifies the model used in latency part. It can be \code{ph} which represents the proportional hazards model, or \code{aft} which represents the accelerated failure time model.
#' @param stdz standardize the covariates. Default is FALSE.
#' @param method.suffix specify the method that will used to accomodate heterogeneous auxiliary information. Note that it is valid only when \code{hetero=TRUE}. The available methods include:
#'    \code{hetero_pen} indicates the penalized auxiliary information synthesis (PAIS) method;
#'    \code{hetero_ma} indicates the model averaging auxiliary information synthesis (MAIS) method (with all possible candidate models);
#'    \code{hetero_ma_nestpen} indicates the model averaging auxiliary information synthesis (MAIS) method (with nested candidate models sorted by PAIS).
#'    \code{hetero_ma_nestfix} indicates the model averaging auxiliary information synthesis (MAIS) method (with nested candidate models defined before).
#' @param tuning a parameter used in the model averaging auxiliary information systhesis method. Default is log(length(yobs))
#' @param auxsort the pre knowledge of the sort of these available auxiliary information. Default is NULL.
#' @param nboot specifies the number of bootstrap sampling. The default \code{nboot = 100}.
#'
#' @export SMC.AuxSP.AIS.fit
SMC.AuxSP.AIS.fit <- function(yobs,delta,X,Z,aux,latency="ph",stdz=FALSE,
                              methods.suffix=c("ori","homo","hetero_pen","hetero_ma","hetero_ma_nestpen","hetero_ma_nestfix"),
                              tuning=NULL,auxsort=NULL,nboot=100,tau.true=NULL
){

  # tau.true is used to produce oracle homo estimator
  if(!is.null(tau.true)){methods.suffix <- c(methods.suffix,"oracle")} # only for simulatio purpose

  #-------------------------------------------------------#
  ### Specify the dimension
  pX <- ncol(X)   # number of parameters in internal data for latency part
  pZ <- ncol(Z)+1
  p <- pX+pZ # total number of parameters
  bet.names <- colnames(X)
  gam.names <- c('Intercept',colnames(Z))
  N <- length(yobs)   # sample size in internal data
  K <- length(aux) # number of external studies
  if(is.null(tuning)){tuning <- log(length(yobs))}

  #-------------------------------------------------------#
  ### prepare the auxinfo matrix
  tstar.unique <- unique(unlist(lapply(aux,function(studyi){sapply(studyi,function(timej){timej$tstar})})))
  auxinfo <- as.numeric()
  for(istudy in 1:length(aux)){ # istudy <- 2; itime <- 1
    for(itime in 1:length(aux[[istudy]])){
      aux.c <- aux[[istudy]][[itime]]
      K.c <- length(aux[[istudy]][[itime]]$sprob)
      auxinfo <- rbind(
        auxinfo,
        cbind(rep(istudy,K.c),
              rep(aux.c$tstar,K.c),
              rep(which(tstar.unique==aux.c$tstar),K.c),
              aux.c$sprob,
              aux.c$M,
              aux.c$gfunc(X,Z)))
    }
  }
  colnames(auxinfo) <- c('study','tstar','tstaridx','sprob','M',paste('ind',1:N,sep="")) # rename this matrix
  if(stdz==TRUE){
    X <- apply(X,2,function(x){(x-mean(x))/sd(x)})
    Z <- apply(Z,2,function(x){(x-mean(x))/sd(x)})
  }

  #-------------------------------------------------------#
  ### Specify functions
  if(latency=="ph"){
    SMC.fit <- SMC.PH.fit
    SMC.AuxSP.EE <- SMC.AuxSP.EE.PH
  }else if(latency=="aft"){
    SMC.fit <- SMC.AFT.fit
    SMC.AuxSP.EE <- SMC.AuxSP.EE.AFT
  }

  #-------------------------------------------------------#
  ### Ests for original estimators (without auxiliary information)
  smcfit <- SMC.fit(yobs,delta,X,Z)
  info.ori <- list(convergence=smcfit$convergence)
  bet.ori <- smcfit$latency[,1]
  gam.ori <- smcfit$incidence[,1]
  dif <- SMC.AuxSP.EE(bet.ori,gam.ori,yobs,delta,X,Z,smcfit$w,auxinfo); dif
  pai <- apply(auxinfo[,-c(1:5)],1,mean)

  #-------------------------------------------------------#
  ### bootstrap proceduce for preparation
  bet.all <- array(NA,dim=c(nboot,pX))
  gam.all <- array(NA,dim=c(nboot,pZ))
  dif.all <- sprob.KM.all <- array(NA,dim=c(nboot,nrow(auxinfo)))
  pai.all <- array(NA,dim=c(nboot,nrow(auxinfo)))
  iboot <- 1; nwrong <- 0
  while(iboot <= nboot){
    if(nwrong > (nboot*0.4)){stop("An error!")}
    ### get boot sample idx
    idx <- sample(1:N,N,replace=TRUE)
    ### Fit the current bootstrap
    smctry <- try({
      smcfit.iboot <- SMC.fit(yobs[idx],delta[idx],X[idx,,drop=F],Z[idx,,drop=F])
    }, silent = T) # smcfit.iboot$convergence==FALSE
    if( is(smctry, "try-error") == TRUE){nwrong<-nwrong+1;next} # | smcfit.iboot$convergence==FALSE
    ### Get the value for beta and gam
    w.iboot <- smcfit.iboot$w
    bet.iboot <- smcfit.iboot$latency[,1]; gam.iboot <- smcfit.iboot$incidence[,1]
    bet.all[iboot,] <- bet.iboot
    gam.all[iboot,] <- gam.iboot
    ### get the initial estimate of the probs (to get consistent estimator of 0s)
    sprob.KM.all[iboot,] <- sapply(1:nrow(auxinfo),function(iaux){
      St.Sub.KM(auxinfo[iaux,'tstar'],yobs[idx],delta[idx],auxinfo[iaux,idx+5,drop=F])})
    dif.all[iboot,] <- SMC.AuxSP.EE(bet.iboot,gam.iboot,yobs[idx],delta[idx],X[idx,,drop=F],
                                    Z[idx,,drop=F],w.iboot,auxinfo[,c(1:5,idx+5)])
    pai.all[iboot,] <- apply(auxinfo[,idx+5],1,mean)
    ### prepare for next iteration
    iboot <- iboot + 1
  }
  (bet.ori.se <- apply(bet.all,2,sd))
  (gam.ori.se <- apply(gam.all,2,sd))
  V.KM <- cov(sprob.KM.all)*N
  Sigma.all <- cov(cbind(bet.all,gam.all,dif.all))*N
  dis.all <- pai.all * AuxSP.disturb.KM(nboot,V.KM,auxinfo) # for bootstraping


  #-------------------------------------------------------#
  ### obtain matrices: Sigma.all+V.KM+c(bet.ori,gam.ori)+dif+N+auxinfo
  # get some matrices
  Sigma <- Sigma.all[1:p,1:p]
  Gam <- Sigma.all[(p+1):(p+nrow(auxinfo)),1:p,drop=F]
  V <- Sigma.all[(p+1):(p+nrow(auxinfo)),(p+1):(p+nrow(auxinfo)),drop=F]
  # reformulate V (to avoid singularity)
  V.new <- V
  for(k in 1:K){
    idxk <- (auxinfo[,'study']==k)
    VVk <- diag(pai[idxk])%*%V.KM[idxk,idxk]%*%diag(pai[idxk])
    V.new[idxk,idxk] <- V.new[idxk,idxk] + VVk*(N/auxinfo[idxk,'M'][1])
  }
  V.inv <- MASS::ginv(V.new,tol=1e-6) # eigen(V.new)$values solve(V.new)  V %*% V.inv %*% V
  A <- Sigma-t(Gam)%*%V.inv%*%Gam
  A.inv <- solve(A)
  W <- rbind(cbind(A.inv,-A.inv%*%t(Gam)%*%V.inv), # W0 <- solve(Sigma.all)
             cbind(-V.inv%*%Gam%*%A.inv,V.inv+V.inv%*%Gam%*%A.inv%*%t(Gam)%*%V.inv))
  W.inv <- MASS::ginv(W,tol=1e-6) # Sigma.all (a new version)


  #-------------------------------------------------------#
  ### estimators with auxiliary information - homo
  # obtain estimates
  if("homo" %in% methods.suffix){
    the.homo <- as.vector(c(bet.ori,gam.ori)-t(V.inv%*%Gam)%*%dif);the.homo
    the.homo.se  <- sqrt( diag((A)/N) ); the.homo.se # sqrt( diag((Sigma-2*t(Gam)%*%V.inv%*%Gam+t(lam)%*%V%*%lam)/N) )
    bet.homo <- the.homo[1:pX]; gam.homo <- the.homo[-c(1:pX)]
    bet.homo.se <- the.homo.se[1:pX]; gam.homo.se <- the.homo.se[-c(1:pX)]
    # do test: Chisq tests
    test.chisq.value <- as.vector(N*t(dif)%*%V.inv%*%dif)
    test.chisq.dfreedom <- nrow(auxinfo)
    test.chisq.pvalue <- 1-pchisq(test.chisq.value,test.chisq.dfreedom) # smaller is worse
    test.chisq <- data.frame(value= test.chisq.value,df=test.chisq.dfreedom,pvalue=test.chisq.pvalue)
    info.homo <- list(test.chisq=test.chisq)
  }

  #-------------------------------------------------------#
  ### estimators with auxiliary information - oracle homo !!!!
  if("oracle" %in% methods.suffix){
    idx.oracle <- which(tau.true==0)
    Sigma2 <- Sigma
    Gam2 <- Gam[idx.oracle,,drop=F]
    V.new2 <- V.new[idx.oracle,idx.oracle,drop=F]
    V.inv2 <- MASS::ginv(V.new2,tol=1e-6) # eigen(V.new2)$values solve(V.new2)  V2 %*% V.inv2 %*% V2
    the.oracle <- as.vector(c(bet.ori,gam.ori)-t(V.inv2%*%Gam2)%*%dif[idx.oracle]);the.oracle
    A2 <- Sigma2-t(Gam2)%*%V.inv2%*%Gam2
    the.oracle.se  <- sqrt( diag((A2)/N) ); the.oracle.se # sqrt( diag((Sigma-2*t(Gam)%*%V.inv%*%Gam+t(lam)%*%V%*%lam)/N) )
    bet.oracle <- the.oracle[1:pX]; gam.oracle <- the.oracle[-c(1:pX)]
    bet.oracle.se <- the.oracle.se[1:pX]; gam.oracle.se <- the.oracle.se[-c(1:pX)]
    info.oracle <- list()
  }

  #-------------------------------------------------------#
  ### estimators with auxiliary information - hetero - penalty
  if("hetero_pen" %in% methods.suffix){
    # transform necessary matrices: A+Gam+V.inv+auxinfo+nu
    V.inv.SVD <- svd(V.inv)
    V.inv.root <- V.inv.SVD$u%*%diag(sqrt(V.inv.SVD$d))%*%t(V.inv.SVD$v)
    # solve adaptive lasso using lars
    sol.penfit <-   SMC.AuxSP.AIS.Pen.fit(bet=bet.ori,gam=gam.ori,dif=dif,pai=pai,V.inv.root=V.inv.root,
                                          Gam=Gam,V.inv=V.inv,auxinfo=auxinfo)
    # obtain final estimates
    tau.hetero_pen <- sol.penfit$tau; tau.hetero_pen
    the.hetero_pen <-  sol.penfit$the; the.hetero_pen
    tau.path <- sol.penfit$tau.path
    # bootstrap type estimator of variance
    the.hetero_pen.all <- array(NA,dim=c(nboot,pX+pZ))
    tau.hetero_pen.all <- array(NA,dim=c(nboot,nrow(auxinfo)))
    tau.path.all <- rep(list(list()),nboot)
    for(iboot in 1:nboot){
      sol.penfit.iboot <- SMC.AuxSP.AIS.Pen.fit(
        bet=bet.all[iboot,],gam=gam.all[iboot,],dif=dif.all[iboot,]-dis.all[iboot,],pai=pai.all[iboot,],
        V.inv.root=V.inv.root,Gam=Gam,V.inv=V.inv,auxinfo=auxinfo)
      # obtain final estimates
      tau.hetero_pen.all[iboot,] <- sol.penfit.iboot$tau
      the.hetero_pen.all[iboot,] <- sol.penfit.iboot$the
      tau.path.all[[iboot]]      <- sol.penfit.iboot$tau.path
    }
    (tau.hetero_pen.se <- apply(tau.hetero_pen.all,2,sd))
    (the.hetero_pen.se <- apply(the.hetero_pen.all,2,sd))
    bet.hetero_pen <- the.hetero_pen[1:pX]; gam.hetero_pen <- the.hetero_pen[-c(1:pX)]
    bet.hetero_pen.se <- the.hetero_pen.se[1:pX]; gam.hetero_pen.se <- the.hetero_pen.se[-c(1:pX)]
    info.hetero_pen <- list(tau=tau.hetero_pen)
  }

  #-------------------------------------------------------#
  ### estimators with auxiliary information - hetero - model averaging - all candidate models
  if("hetero_ma" %in% methods.suffix){
    # prepare indices for each submodel: use all the 2^M combinations as candidate models
    indices <- matrix(0,nrow=2^nrow(auxinfo),ncol=nrow(auxinfo))
    idxxx <- 2;i0 <- 0
    for(Km in 1:nrow(auxinfo)){
      strm <- paste(
        paste("for(i",1:Km," in ",paste("(i",0:(Km-1),"+1)",sep=""),":(nrow(auxinfo)-",(Km-1):0,")){",sep="",collapse=""),
        paste("indices[idxxx,c(",paste("i",1:Km,sep="",collapse=","),")]<-1",sep=""),
        ";idxxx<-idxxx+1",
        paste(rep("}",Km),collapse=""),
        sep=""
      )
      eval(parse(text=strm))
    }
    # the estimator
    sol.mafit <- SMC.AuxSP.AIS.MA.fit(indices=indices,bet=bet.ori,gam=gam.ori,dif=dif,pai=pai,
                                      W=W,auxinfo=auxinfo,tuning=tuning)
    thetau.hetero_ma <- sol.mafit$thetau
    ws.hetero_ma <- sol.mafit$ws
    bet.hetero_ma <- thetau.hetero_ma[1:pX]
    gam.hetero_ma <- thetau.hetero_ma[(pX+1):(pX+pZ)]
    tau.hetero_ma <- thetau.hetero_ma[-c(1:(pX+pZ))]
    # obtain variance (bootstrap) SE: the corresponding bootstrap standard error
    thetau.hetero_ma.all <- array(NA,dim=c(nboot,pX+pZ+nrow(auxinfo)))
    for(iboot in 1:nboot){
      sol.mafit.iboot <- SMC.AuxSP.AIS.MA.fit(
        indices=indices,bet=bet.all[iboot,],gam=gam.all[iboot,],dif=dif.all[iboot,]-dis.all[iboot,],
        pai=pai.all[iboot,],W=W,auxinfo=auxinfo,tuning=tuning)
      thetau.hetero_ma.all[iboot,] <- sol.mafit.iboot$thetau
    }
    (thetau.hetero_ma.se <- apply(thetau.hetero_ma.all,2,sd))
    bet.hetero_ma.se <- thetau.hetero_ma.se[1:pX]
    gam.hetero_ma.se <- thetau.hetero_ma.se[(pX+1):(pX+pZ)]
    tau.hetero_ma.se <- thetau.hetero_ma.se[-c(1:(pX+pZ))]
    dimnames(indices) <- list(paste("Candidate",1:nrow(indices),sep=""),paste("aux",1:nrow(auxinfo),sep=""))
    candidates.hetero_ma <- data.frame(weight=round(ws.hetero_ma,4),indices)
    info.hetero_ma <- list(tau=tau.hetero_ma,candidates=candidates.hetero_ma)
  }

  #-------------------------------------------------------#
  ### estimators with auxiliary information - hetero - model averaging - nested candidate models - sort by fixed pattern (pre-defined)
  if("hetero_ma_nestfix" %in% methods.suffix){
    # prepare indices for each submodel: nested - fixed ordering pre-defined
    if(is.null(auxsort)){auxsort <- 1:nrow(auxinfo)}
    indices <- rbind(0,apply(diag(length(auxsort))[,auxsort],2,cummax))
    # the estimator
    sol.mafit <- SMC.AuxSP.AIS.MA.fit(indices=indices,bet=bet.ori,gam=gam.ori,dif=dif,pai=pai,
                                      W=W,auxinfo=auxinfo,tuning=tuning)
    thetau.hetero_ma_nestfix <- sol.mafit$thetau
    ws.hetero_ma_nestfix <- sol.mafit$ws
    bet.hetero_ma_nestfix <- thetau.hetero_ma_nestfix[1:pX]
    gam.hetero_ma_nestfix <- thetau.hetero_ma_nestfix[(pX+1):(pX+pZ)]
    tau.hetero_ma_nestfix <- thetau.hetero_ma_nestfix[-c(1:(pX+pZ))]
    # obtain variance (SE): the corresponding bootstrap standard error
    thetau.hetero_ma_nestfix.all <- array(NA,dim=c(nboot,pX+pZ+nrow(auxinfo)))
    for(iboot in 1:nboot){
      sol.mafit.iboot <- SMC.AuxSP.AIS.MA.fit(
        indices=indices,bet=bet.all[iboot,],gam=gam.all[iboot,],dif=dif.all[iboot,]-dis.all[iboot,],
        pai=pai.all[iboot,],W=W,auxinfo=auxinfo,tuning=tuning)
      thetau.hetero_ma_nestfix.all[iboot,] <- sol.mafit.iboot$thetau
    }
    (thetau.hetero_ma_nestfix.se <- apply(thetau.hetero_ma_nestfix.all,2,sd))
    bet.hetero_ma_nestfix.se <- thetau.hetero_ma_nestfix.se[1:pX]
    gam.hetero_ma_nestfix.se <- thetau.hetero_ma_nestfix.se[(pX+1):(pX+pZ)]
    tau.hetero_ma_nestfix.se <- thetau.hetero_ma_nestfix.se[-c(1:(pX+pZ))]
    dimnames(indices) <- list(paste("Candidate",1:nrow(indices),sep=""),paste("aux",1:nrow(auxinfo),sep=""))
    candidates.hetero_ma_nestfix <- data.frame(weight=round(ws.hetero_ma_nestfix,4),indices)
    info.hetero_ma_nestfix <- list(tau=tau.hetero_ma_nestfix,candidates=candidates.hetero_ma_nestfix)
  }

  #-------------------------------------------------------#
  ### estimators with auxiliary information - hetero - model averaging - nested candidate models - sort by penalty
  if("hetero_ma_nestpen" %in% methods.suffix){
    # prepare indices for each submodel: nested - penalty
    indices <- as.matrix(t(1*(tau.path!=0)))
    # the estimator
    sol.mafit <- SMC.AuxSP.AIS.MA.fit(indices=indices,bet=bet.ori,gam=gam.ori,dif=dif,pai=pai,
                                      W=W,auxinfo=auxinfo,tuning=tuning)
    thetau.hetero_ma_nestpen <- sol.mafit$thetau
    ws.hetero_ma_nestpen <- sol.mafit$ws
    bet.hetero_ma_nestpen <- thetau.hetero_ma_nestpen[1:pX]
    gam.hetero_ma_nestpen <- thetau.hetero_ma_nestpen[(pX+1):(pX+pZ)]
    tau.hetero_ma_nestpen <- thetau.hetero_ma_nestpen[-c(1:(pX+pZ))]
    # obtain variance (SE): the corresponding bootstrap standard error
    thetau.hetero_ma_nestpen.all <- array(NA,dim=c(nboot,pX+pZ+nrow(auxinfo)))
    for(iboot in 1:nboot){
      # prepare indices for each submodel:  sort the auxiliary information using penalty method
      indices.iboot <- as.matrix(t(1*(tau.path.all[[iboot]]!=0)))
      sol.mafit.iboot <- SMC.AuxSP.AIS.MA.fit(
        indices=indices.iboot,bet=bet.all[iboot,],gam=gam.all[iboot,],dif=dif.all[iboot,]-dis.all[iboot,],
        pai=pai.all[iboot,],W=W,auxinfo=auxinfo,tuning=tuning)
      thetau.hetero_ma_nestpen.all[iboot,] <- sol.mafit.iboot$thetau
    }
    (thetau.hetero_ma_nestpen.se <- apply(thetau.hetero_ma_nestpen.all,2,sd))
    bet.hetero_ma_nestpen.se <- thetau.hetero_ma_nestpen.se[1:pX]
    gam.hetero_ma_nestpen.se <- thetau.hetero_ma_nestpen.se[(pX+1):(pX+pZ)]
    tau.hetero_ma_nestpen.se <- thetau.hetero_ma_nestpen.se[-c(1:(pX+pZ))]
    dimnames(indices) <- list(paste("Candidate",1:nrow(indices),sep=""),paste("aux",1:nrow(auxinfo),sep=""))
    candidates.hetero_ma_nestpen <- data.frame(weight=round(ws.hetero_ma_nestpen,4),indices)
    info.hetero_ma_nestpen <- list(tau=tau.hetero_ma_nestpen,candidates=candidates.hetero_ma_nestpen)
  }


  #-------------------------------------------------------#
  ### do inference and combine results into pre-specified style
  out0 <- list()
  for(suffix in methods.suffix){
    # for latency part
    bet.c <- get(paste(c("bet"),suffix,sep="."))
    bet.c.se <- get(paste(c("bet"),suffix,"se",sep="."))
    zvalue.bet.c <- bet.c/bet.c.se
    pvalue.bet.c <- 2*(1-pnorm(abs(zvalue.bet.c)))
    # for incidence part
    gam.c <- get(paste(c("gam"),suffix,sep="."))
    gam.c.se <- get(paste(c("gam"),suffix,"se",sep="."))
    zvalue.gam.c <- gam.c/gam.c.se
    pvalue.gam.c <- 2*(1-pnorm(abs(zvalue.gam.c)))
    # combine inference results
    latency.c <- data.frame(Est=bet.c,SE=bet.c.se,zvalue=zvalue.bet.c,pvalue=pvalue.bet.c,row.names=bet.names)
    incidence.c <- data.frame(Est=gam.c,SE=gam.c.se,zvalue=zvalue.gam.c,pvalue=pvalue.gam.c,row.names=gam.names)
    # combine all results
    info.c <- get(paste(c("info"),suffix,sep="."))
    out0.c <- c(list(latency = latency.c,incidence = incidence.c),info.c)
    out0 <- c(out0,list(out0.c))
  }
  names(out0) <- paste("res",methods.suffix,sep=".")
  out0

  #-------------------------------------------------------#
  ### extract output values
  out <- c(list(model=paste("Semi-parametric mixcure " ,"logistic/",latency," model",sep="")),
           out0)
  return(out)

}



#==== Auxiliary information's Estimating Equations at different time points (for ph) ====#
#' @title Estimating equation for the auxiliary subgroup survival rates under the proportional hazards mixture cure model
#'
#' @description Calculate the estimating equation for the auxiliary subgroup survival rates under the proportional hazards mixture cure model
#'
#' @aliases SMC.AuxSP.EE.PH
#'
#' @param bet unknown parameters corresponding to latency part.
#' @param gam unknown parameters corresponding to incidence part.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param Z a matrix of covariates corresponding to incidence part.
#' @param w conditional probability of the individual remaining uncured.
#' @param auxinfo a matrix that collects the auxiliary information in a pre-specified way.
#'
#' @export SMC.AuxSP.EE.PH
SMC.AuxSP.EE.PH <- function(bet,gam,yobs,delta,X,Z,w,auxinfo){

  # the estimating equations in individual levels
  Stx <- SMC.PH.Stx(yobs,delta,X,Z,bet,gam,w,cross=TRUE,tm=auxinfo[,"tstar"])$Stx
  Psi <- (t(Stx)-auxinfo[,'sprob'])*auxinfo[,-c(1:5)]

  # output
  return(apply(Psi,1,mean))
}

#==== Auxiliary information's Estimating Equations at different time points (for aft) ====#
#' @title Estimating equation for the auxiliary subgroup survival rates under the accelerated failure time mixture cure model
#'
#' @description Calculate the estimating equation for the auxiliary subgroup survival rates under the accelerated failure time mixture cure model
#'
#' @aliases SMC.AuxSP.EE.AFT
#'
#' @param bet unknown parameters corresponding to latency part.
#' @param gam unknown parameters corresponding to incidence part.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param Z a matrix of covariates corresponding to incidence part.
#' @param w conditional probability of the individual remaining uncured.
#' @param auxinfo a matrix that collects the auxiliary information in a pre-specified way.
#'
#' @export SMC.AuxSP.EE.AFT
SMC.AuxSP.EE.AFT <- function(bet,gam,yobs,delta,X,Z,w,auxinfo){
  # Auxiliary Information's Estimating Equations at different time points

  # the estimating equations in individual levels
  Stx <- SMC.AFT.Stx(yobs,delta,X,Z,bet,gam,w,cross=TRUE,tm=auxinfo[,"tstar"])$Stx
  Psi <- (t(Stx)-auxinfo[,'sprob'])*auxinfo[,-c(1:5)]

  # output
  return(apply(Psi,1,mean))
}


#==== Sub-function: for fitting penaltized estimator ====#
#' @title Penalized auxiliary information synthesis method under the mixture cure model
#'
#' @description Fit the penalized auxiliary information synthesis method under the mixture cure model
#'
#' @aliases SMC.AuxSP.AIS.Pen.fit
#'
#' @param indices the indices for candidate models (one row indicate one candidate model).
#' @param bet unknown parameters corresponding to latency part.
#' @param gam unknown parameters corresponding to incidence part.
#' @param dif the calculated values of estimating equations related to auxiliary information.
#' @param pai the frequencies of each subgroup.
#' @param V.inv.root the square root of the inverse of the variance-covariance matrix of the estimated coefficients based on the internal data only
#' @param Gam correlation between the estimated coefficients based on the internal data only and the estimating equations related to auxiliary information.
#' @param V.inv the inverse of the variance-covariance matrix of the estimated coefficients based on the internal data only
#' @param auxinfo a matrix that collects the auxiliary information in a pre-specified way.
#'
#' @export SMC.AuxSP.AIS.Pen.fit
SMC.AuxSP.AIS.Pen.fit <- function(bet,gam,dif,pai,V.inv.root,Gam,V.inv,auxinfo){

  N <- ncol(auxinfo[,-c(1:5)])
  y.tilde <- as.vector( V.inv.root%*%dif )
  X.tilde <- V.inv.root%*%diag(pai)

  # solve adaptive lasso using lars
  w <- (1/abs(dif))*(pai)
  X.tilde.star <- t(t(X.tilde)/w)
  sol.lars <- lars::lars(X.tilde.star,y.tilde,trace=FALSE,normalize=FALSE,intercept=FALSE)
  tau.path <- t(as.matrix(sol.lars$beta))/w # each
  tau.path.RSS <- apply(X.tilde %*% tau.path - y.tilde,2,function(x){sum(x^2)})
  tau.path.Card <- apply(tau.path,2,function(x){sum(x!=0)})
  IC.all <- as.numeric( tau.path.RSS + log(N)/N * tau.path.Card ) # BIC type criterion
  min_IC.idx <- which.min( IC.all  )
  # output: return final estimates
  tau <- tau.path[,min_IC.idx]
  return(list(
    the=as.vector(c(bet,gam)-t(V.inv%*%Gam)%*%(dif-pai*tau)),
    tau=tau,
    tau.path=tau.path
  ))

}



#==== Sub-function: for fitting model averaging estimator ====#
#' @title Model averaging auxiliary information synthesis method under the mixture cure model
#'
#' @description Fit the model averaging auxiliary information synthesis method under the mixture cure model
#'
#' @aliases SMC.AuxSP.AIS.MA.fit
#'
#' @param indices the indices for candidate models (one row indicate one candidate model).
#' @param bet unknown parameters corresponding to latency part.
#' @param gam unknown parameters corresponding to incidence part.
#' @param dif the calculated values of estimating equations related to auxiliary information.
#' @param pai the frequencies of each subgroup.
#' @param W a matrix
#' @param auxinfo a matrix that collects the auxiliary information in a pre-specified way.
#' @param tuning a parameter used in the model averaging auxiliary information systhesis method. Default is log(length(yobs))
#'
#' @export SMC.AuxSP.AIS.MA.fit
SMC.AuxSP.AIS.MA.fit <- function(indices,bet,gam,dif,pai,W,auxinfo,tuning){

  # some preparation
  pX <- length(bet); pZ <- length(gam); N <- ncol(auxinfo[,-c(1:5)])
  # prepare matrices
  betgamtau.all <- Gpre <- array(NA,dim=c(nrow(indices),pX+pZ+nrow(auxinfo)))
  for(m in 1:nrow(indices)){
    Paim <- diag(nrow(auxinfo))[indices[m,]==1,,drop=F]
    Paim.aug <- as.matrix(Matrix::bdiag(diag(pX+pZ),Paim))
    # est for a single estimator from a candidate model
    Qm <- as.matrix(Matrix::bdiag(diag(pX+pZ),Paim%*%diag(pai)))
    Pm <- solve(Qm%*%W%*%t(Qm))%*%Qm%*%W
    betgamtau.m <- as.vector(Pm%*%c(bet,gam,dif))
    # record elemtns
    Gpre[m,] <- as.vector( t(Qm) %*% betgamtau.m - c(bet,gam,dif) )
    betgamtau.all[m,] <- as.vector( t(Paim.aug)%*%betgamtau.m )
  }
  G <- Gpre%*%W%*%t(Gpre)
  g <- apply(indices,1,sum)+pX+pZ

  # calculate weights and final estimator
  # tuning <- log(N) # 2-AIC; log(N)-BIC
  ipop.sol <- kernlab::ipop(c=(tuning/2)*g/N,H=G,A=t(rep(1,nrow(indices))),b=1,r=0,l=rep(0,nrow(indices)),u=rep(1,nrow(indices)))
  ws <- kernlab::primal(ipop.sol)
  thetau <- apply(betgamtau.all * ws,2,sum)

  # output
  return(list(
    thetau=thetau,
    ws=ws
  ))

}





