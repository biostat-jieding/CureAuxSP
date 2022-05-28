#==========================================================================#
# Semi-parametric cox mixture cure model  -- from smcure ####
#==========================================================================#

#==== The main function for fitting model ====#
#' @title Semiparametric proportional hazards mixture cure model
#'
#' @description Fit semiparametric proportional hazards (PH) mixture cure model by the EM algorithm.
#'
#' @aliases SMC.PH.fit
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param Z a matrix of covariates corresponding to incidence part.
#' @param incidence specifies the link in incidence part. The "logit", "probit" or complementary loglog ("cloglog") links are available. By default incidence = c("logit").
#' @param Var a logical value. If it is \code{TRUE}, the program returns Standard errors based on the bootstrap method is used. By default, \code{Var = FALSE}.
#' @param nboot specifies the number of bootstrap sampling. The default \code{nboot = 100}.
#' @param em.maxit specifies the maximum iteration number. If the convergence criterion is not met, the EM iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{em.maxit = 100}.
#' @param em.tol tolerance for convergence. The default is \code{em.tol = 1e-5}. Iteration stops once the relative change in deviance is less than \code{em.tol}.
#'
#' @export SMC.PH.fit
SMC.PH.fit <- function(yobs,delta,X,Z,incidence=c("logit"),Var=FALSE,nboot=100,
                       em.maxit=100,em.tol=1e-5){

  ### Specify the dimension and prepare data
  N <- length(yobs) # number of individuals
  pX <- ncol(X)
  pZ <- ncol(Z)+1
  ZI <- as.matrix(cbind(rep(1,N),Z)) # [dataframe: for incidence part]

  ### do EM algorithm
  # obtain initial bet, gam and baseline nonparametric part
  bet.old <- survival::coxph(survival::Surv(yobs,delta)~X,subset=(delta==1),method="breslow")$coef
  gam.old <- eval(parse(text=paste("glm", "(", "delta~Z",",family = quasibinomial(link='",incidence,"'",")",")",sep = "")))$coef
  Sutx.old <- SMC.PH.Sutx(yobs,delta,X,bet.old,delta,cross=FALSE,tm=NULL)$Sutx
  numit <- 1
  repeat{

    # calculate w first
    expZgam <- exp(ZI%*%gam.old)
    uncureprob <- expZgam/(1+expZgam)
    w <- as.vector(delta+(1-delta)*(uncureprob*Sutx.old)/(1-uncureprob+uncureprob*Sutx.old))

    # update gam
    gam <- eval(parse(text=paste("glm","(","w~Z",",family = quasibinomial(link='",incidence,"'", ")",")", sep = "")))$coef
    # update bet
    bet <- survival::coxph(survival::Surv(yobs,delta)~X+offset(log(w)),subset=(w!=0),method="breslow")$coef
    # bet <- coxph(Surv(yobs[w>0],delta[w>0])~X[w>0,,drop=F],weights=w[w>0],method="breslow")$coef
    # update Sutx
    Sutx <- SMC.PH.Sutx(yobs,delta,X,bet,w,cross=FALSE,tm=NULL)$Sutx

    ### update or stop
    convergence.value <- max(c(abs(bet.old-bet)/abs(bet),abs(gam.old-gam)/abs(gam)))
    if(convergence.value >= em.tol & numit < em.maxit){
      bet.old <- bet; gam.old <- gam
      Sutx.old <- Sutx
      numit <- numit + 1
    }else{
      break
    }

  } # end em reps
  convergence <- (numit<em.maxit & convergence.value < em.tol)

  ### extract output values
  out <- list(
    latency = data.frame(Est=bet,row.names=colnames(X)),
    incidence = data.frame(Est=gam,row.names=c('Intercept',colnames(Z))),
    convergence=convergence,numit=numit,
    w=w
  )
  return(out)

}


#==== Survival function for uncured patients (with PH latency) ====#
#' @title Survival function for uncured patients (with PH latency)
#'
#' @description Calculate the survival function for uncured patients (with PH latency) at various time points.
#'
#' @aliases SMC.PH.Sutx
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param bet unknown parameters corresponding to latency part, or the fitted values of regression coefficients in the latency part returned by the function \code{\link{SMC.PH.fit()}}.
#' @param w conditional probability of the individual remaining uncured.
#' @param Var a logical value. If it is \code{TRUE}, the program returns Standard errors based on the bootstrap method is used. By default, \code{Var = FALSE}.
#' @param cross a logical value. The default is \code{FALSE}.
#' @param tm the time points that the survival function will be calculated at. The default is \code{NULL}.
#'
#' @export SMC.PH.Sutx
SMC.PH.Sutx <- function(yobs,delta,X,bet,w,cross=FALSE,tm=NULL){ # for uncured

  # preparation: calculate hazards
  expXbet <- as.vector(exp(X%*%bet))
  sumRisk <- sapply(yobs,function(yobsi){sum((yobs>=yobsi)*w*expXbet)})
  ht <- ifelse(delta==1,delta/sumRisk,0)
  # calculate
  if(cross){
    # calculate hazards
    Ht <- sapply(tm,function(tmi){sum((yobs<=tmi)*ht)})
    Ht[tm>max(yobs[delta==1])] <- Inf; Ht[tm<min(yobs[delta==1])] <- 0
    # calculate Survivals
    St.baseline <- exp(-Ht)
    Sutx <- outer(expXbet,St.baseline,function(x,y){y^x})
  }else{
    # calculate hazards
    Ht <- sapply(yobs,function(yobsi){sum((yobs<=yobsi)*ht)})
    Ht[yobs>max(yobs[delta==1])] <- Inf; Ht[yobs<min(yobs[delta==1])] <- 0
    # calculate Survivals
    St.baseline <- exp(-Ht)
    Sutx <- St.baseline^expXbet
  }
  # output
  list(Sutx = Sutx)

}


#==== Survival function for the whole population (with PH latency) ====#
#' @title Survival function for the whole population (with PH latency)
#'
#' @description Calculate the survival function for the whole population (with PH latency) at various time points.
#'
#' @aliases SMC.PH.Stx
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param bet unknown parameters corresponding to latency part,   or the fitted values of regression coefficients in the latency part returned by the function \code{\link{SMC.PH.fit()}}.
#' @param gam unknown parameters corresponding to incidence part, or the fitted values of regression coefficients in the incidence part returned by the function \code{\link{SMC.PH.fit()}}.
#' @param w conditional probability of the individual remaining uncured.
#' @param Var a logical value. If it is \code{TRUE}, the program returns Standard errors based on the bootstrap method is used. By default, \code{Var = FALSE}.
#' @param cross a logical value. The default is \code{FALSE}.
#' @param tm the time points that the survival function will be calculated at. The default is \code{NULL}.
#'
#' @export SMC.PH.Stx
SMC.PH.Stx <- function(yobs,delta,X,Z,bet,gam,w,cross=FALSE,tm=NULL){ # for all

  # calculate uncureprob
  expZgam <- exp(cbind(1,Z)%*%gam)
  uncureprob <- as.vector(expZgam/(1+expZgam))
  Sutx <- SMC.PH.Sutx(yobs,delta,X,bet,w,cross,tm)$Sutx
  Stx <- Sutx*uncureprob+1-uncureprob
  # output
  return(list(Stx=Stx))

}





#==========================================================================#
# Semi-parametric aft mixture cure model  -- from smcure ####
#==========================================================================#

#==== The main function for fitting model ====#
#' @title Semiparametric accelerated failure time mixture cure model
#'
#' @description Fit semiparametric accelerated failure time (AFT) mixture cure model by the EM algorithm.
#'
#' @aliases SMC.AFT.fit
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param Z a matrix of covariates corresponding to incidence part.
#' @param incidence specifies the link in incidence part. The "logit", "probit" or complementary loglog ("cloglog") links are available. By default incidence = c("logit").
#' @param intercept a logical value. If \code{TRUE}, the intercept term in the latency part will be included.
#' @param Var a logical value. If it is \code{TRUE}, the program returns Standard errors based on the bootstrap method is used. By default, \code{Var = FALSE}.
#' @param nboot specifies the number of bootstrap sampling. The default \code{nboot = 100}.
#' @param em.maxit specifies the maximum iteration number. If the convergence criterion is not met, the EM iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{em.maxit = 100}.
#' @param em.tol tolerance for convergence. The default is \code{em.tol = 1e-5}. Iteration stops once the relative change in deviance is less than \code{em.tol}.
#'
#' @export SMC.AFT.fit
SMC.AFT.fit <- function(yobs,delta,X,Z,incidence=c("logit"),intercept=FALSE,
                        Var=FALSE,nboot=100,em.maxit=100,em.tol=1e-5){

  ### Specify the dimension and prepare data
  N <- length(yobs) # number of individuals
  pX <- ifelse(intercept==FALSE,ncol(X),ncol(X)+1)
  pZ <- ncol(Z)+1
  XI <- as.matrix(cbind(rep(1,N),X))
  ZI <- as.matrix(cbind(rep(1,N),Z)) # [dataframe: for incidence part]

  ### do EM algorithm
  # obtain initial bet, gam and baseline nonparametric part
  bet.old <- survival::survreg(survival::Surv(yobs,delta)~X)$coef # from survival package
  if(intercept==FALSE){bet.old <- bet.old[-1]}
  gam.old <- eval(parse(text=paste("glm", "(", "delta~Z",",family = quasibinomial(link='",incidence,"'",")",")",sep = "")))$coef
  Sutx.old <- SMC.AFT.Sutx(yobs,delta,X,bet.old,delta,cross=FALSE,tm=NULL,intercept)$Sutx
  numit <- 1
  repeat{

    # calculate w first
    expZgam <- exp(ZI%*%gam.old)
    uncureprob <- expZgam/(1+expZgam)
    w <- as.vector(delta+(1-delta)*(uncureprob*Sutx.old)/(1-uncureprob+uncureprob*Sutx.old))

    # update gam
    gam <- eval(parse(text=paste("glm","(","w~Z",",family = quasibinomial(link='",incidence,"'", ")",")", sep = "")))$coef
    # update bet
    bet <- optim(par=rep(0,pX),fn=SMC.AFT.rank,method="Nelder-Mead",
                 control=list(maxit=500,reltol=1e-04),
                 yobs=yobs,delta=delta,X=X,w=w,intercept=intercept)$par
    # update Sutx
    Sutx <- SMC.AFT.Sutx(yobs,delta,X,bet,w,cross=FALSE,tm=NULL,intercept)$Sutx

    ### update or stop
    # convergence.value <- max(max(abs(c(gam-gam.old,bet-bet.old))),max(abs(c(Sutx-Sutx.old))))
    convergence.value <- max(c(abs(bet.old-bet)/abs(bet),abs(gam.old-gam)/abs(gam)))
    if(convergence.value >= em.tol & numit < em.maxit){
      bet.old <- bet; gam.old <- gam
      Sutx.old <- Sutx
      numit <- numit + 1
    }else{
      break
    }

  } # end em reps
  convergence <- (numit<em.maxit & convergence.value < em.tol)

  ### extract values
  names.latency <- c('Intercept',colnames(X))
  if(intercept==FALSE){names.latency <- names.latency[-1]}
  names.incidence <- c('Intercept',colnames(Z))
  out <- list(
    latency = data.frame(Est=bet,row.names=names.latency),
    incidence = data.frame(Est=gam,row.names=names.incidence),
    convergence = convergence,numit=numit,
    w=w
  )
  return(out)

}

#==== Survival function for uncured patients (with AFT latency) ====#
#' @title Survival function for uncured patients (with AFT latency)
#'
#' @description Calculate the survival function for uncured patients (with AFT latency) at various time points.
#'
#' @aliases SMC.AFT.Sutx
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param bet unknown parameters corresponding to latency part, or the fitted values of regression coefficients in the latency part returned by the function \code{\link{SMC.AFT.fit()}}.
#' @param w conditional probability of the individual remaining uncured.
#' @param Var a logical value. If it is \code{TRUE}, the program returns Standard errors based on the bootstrap method is used. By default, \code{Var = FALSE}.
#' @param cross a logical value. The default is \code{FALSE}.
#' @param tm the time points that the survival function will be calculated at. The default is \code{NULL}.
#' @param intercept a logical value. If \code{TRUE}, the intercept term in the latency part will be included.
#'
#' @export SMC.AFT.Sutx
SMC.AFT.Sutx <- function(yobs,delta,X,bet,w,cross=FALSE,tm=NULL,intercept=FALSE){ # for uncured

  # prepare error
  if(intercept==FALSE){Xbet<-as.vector(X%*%bet)}else{Xbet<-as.vector(cbind(1,X)%*%bet)}
  error <- log(yobs)-Xbet
  # preparation: calculate hazards
  sumRisk <- sapply(error,function(errori){sum((error>=errori)*w)})
  ht <- ifelse(delta==1,delta/sumRisk,0)
  # calculate
  if(cross){
    # prepare errors in tm cross
    error2 <- outer(Xbet,log(tm),function(x,y){y-x})
    # calculate hazards
    Ht <- sapply(error2,function(errori){sum((error<=errori)*ht)})
    Ht[error>max(error[delta==1])] <- Inf; Ht[error<min(error[delta==1])] <- 0
    Ht <- matrix(Ht,ncol=length(tm))
    # calculate Survivals
    Sutx <- exp(-Ht)
  }else{
    # calculate hazards
    Ht <- sapply(error,function(errori){sum((error<=errori)*ht)})
    Ht[error>max(error[delta==1])] <- Inf; Ht[error<min(error[delta==1])] <- 0
    # calculate Survivals
    Sutx <- exp(-Ht)
  }
  # output
  list(Sutx = Sutx)

}


#==== Survival function for the whole population (with AFT latency) ====#
#' @title Survival function for the whole population (with PH latency)
#'
#' @description Calculate the survival function for the whole population (with AFT latency) at various time points.
#'
#' @aliases SMC.AFT.Stx
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param bet unknown parameters corresponding to latency part,   or the fitted values of regression coefficients in the latency part returned by the function \code{\link{SMC.AFT.fit()}}.
#' @param gam unknown parameters corresponding to incidence part, or the fitted values of regression coefficients in the incidence part returned by the function \code{\link{SMC.AFT.fit()}}.
#' @param w conditional probability of the individual remaining uncured.
#' @param Var a logical value. If it is \code{TRUE}, the program returns Standard errors based on the bootstrap method is used. By default, \code{Var = FALSE}.
#' @param cross a logical value. The default is \code{FALSE}.
#' @param tm the time points that the survival function will be calculated at. The default is \code{NULL}.
#'
#' @export SMC.AFT.Stx
SMC.AFT.Stx <- function(yobs,delta,X,Z,bet,gam,w,cross=FALSE,tm=NULL,intercept=FALSE){ # for all

  # calculate uncureprob
  expZgam <- exp(cbind(1,Z)%*%gam)
  uncureprob <- as.vector(expZgam/(1+expZgam))
  Sutx <- SMC.AFT.Sutx(yobs,delta,X,bet,w,cross,tm,intercept)$Sutx
  Stx <- Sutx*uncureprob+1-uncureprob
  # output
  return(list(Stx=Stx))

}

#==== The rank function (with AFT latency) ====#
#' @title The rank function (with AFT latency)
#'
#' @description Rank estimating equation used in the M-step of the EM algorithm for the AFT mixture cure model.
#'
#' @aliases SMC.AFT.rank
#'
#' @param bet unknown parameters corresponding to latency part, or the fitted values of regression coefficients in the latency part returned by the function \code{\link{SMC.AFT.fit()}}.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates corresponding to latency part.
#' @param w conditional probability of the individual remaining uncured.
#' @param intercept a logical value. If \code{TRUE}, the intercept term in the latency part will be included.
#'
#' @export SMC.AFT.rank
SMC.AFT.rank <- function(bet,yobs,delta,X,w,intercept=FALSE){

  if(intercept==FALSE){
    error <- as.vector(log(yobs)-X%*%bet)
  }else{
    error <- as.vector(log(yobs)-cbind(1,X)%*%bet)
  }
  LossGehan <- mean(sapply(1:length(yobs),function(i){sum((error[i]<error)*abs(error-error[i])*w*delta[i])}))
  return(LossGehan)

}

