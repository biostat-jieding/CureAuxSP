# CureAuxSP
This is an R package used to fit **Cure** models with **Aux**iliary information in type of **S**ubgroup survival **P**robabilities.
- The underlying methods are based on the paper titled "**Efficient auxiliary information synthesis for cure rate model**", which has been published in *Journal of the Royal Statistial Society Series C (Applied Statistics)* with DOI: 10.1093/jrsssc/qlad106.
- We have also wrote a paper titled "**CureAuxSP: An R package for estimating mixture cure models with auxiliary survival probabilities**" that introduces the detailed information about this package and this paper has been published in *Computer Methods and Programs in Biomedicine* with DOI: 10.1016/j.cmpb.2024.108212.
- It can also be downloaded as a standard R package from https://cran.r-project.org/web/packages/CureAuxSP/index.html.

## Package description and included main functions

We focus on cure models: (Semiparametric) **PH mixture cure model** or **AFT mixture cure model** *with or without* **auxiliary subgroup survival probabilities**.

Our provided main functions are (we refer to their help pages for more details):
- **SMC.AuxSP**: fit the models in various ways with syntax
```R
SMC.AuxSP(formula, cureform, sdata, aux = NULL, hetero = FALSE, N = Inf, latency = "PH", nboot = 400)
```
- **print.SMC.AuxSP**: print outputted results from SMC.AuxSP() with syntax
```R
print.SMC.AuxSP(object)
```

## Two numerical illustrations

### An example using a simulated dataset is shown below:
```R
## library
library(survival)
library(CureAuxSP)

## generate both the internal dataset of interest and the external dataset

# - the internal dataset
set.seed(1)
sdata.internal <- sdata.SMC(n = 300)
head(sdata.internal)

# - the external dataset 
set.seed(1)
sdata.external <- sdata.SMC(n = 10000)

## prepare the auxiliary information based on the external dataset

# - define two functions for subgroup splitting
gfunc.t1 <- function(X,Z=NULL){
  rbind((X[,1] <  0 & X[,2] == 0), (X[,1] >= 0 & X[,2] == 0),
        (X[,1] <  0 & X[,2] == 1), (X[,1] >= 0 & X[,2] == 1))}
gfunc.t2 <- function(X,Z=NULL){rbind((X[,2] == 0), (X[,2] == 1))}

# - calculate subgroup survival rates
sprob.t1 <- Probs.Sub(tstar = 1, sdata = sdata.external,
                      G = gfunc.t1(X = sdata.external[,-c(1,2)]))
sprob.t2 <- Probs.Sub(tstar = 2, sdata = sdata.external, 
                      G = gfunc.t2(X = sdata.external[,-c(1,2)]))
cat("Information at t* = 1:", sprob.t1, "\nInformation at t* = 2:", sprob.t2)

# - prepare the set that collects information about auxiliary data
aux <- list(
  time1 = list(tstar = 1, gfunc = gfunc.t1, sprob = c(0.73,0.70,0.88,0.83)),
  time2 = list(tstar = 2, gfunc = gfunc.t2, sprob = c(0.62,0.76)-0.20)
)

## fit the model without auxiliary information
set.seed(1)
sol.PHMC <-  SMC.AuxSP(
  formula = Surv(yobs,delta) ~ X1 + X2, cureform = ~ X1, 
  sdata = sdata.internal, aux = NULL, latency = "PH"
)
print.SMC.AuxSP(object = sol.PHMC)

## fit the model with auxiliary information

# - ignore heterogeneity
set.seed(1)
sol.PHMC.Homo <-  SMC.AuxSP(
  formula = Surv(yobs,delta) ~ X1 + X2, cureform = ~ X1, 
  sdata = sdata.internal, aux = aux, hetero = FALSE, latency = "PH"
)
print.SMC.AuxSP(object = sol.PHMC.Homo)

# - consider heterogeneity
set.seed(1)
sol.PHMC.Hetero <-  SMC.AuxSP(
  formula = Surv(yobs,delta) ~ X1 + X2, cureform = ~ X1, 
  sdata = sdata.internal, aux = aux, hetero = TRUE, latency = "PH"
)
print.SMC.AuxSP(object = sol.PHMC.Hetero)
```

### An example using a real dataset from TCGA program is shown below:
```R
## library
library(survival)
library(CureAuxSP)

## prepare the breast cancer dataset

# - download clinical data from the TCGA website
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical", file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

# - a preparation
sdata.pre <- data.frame(
  yobs   = ifelse(!is.na(clinical[,'days_to_death']),clinical[,'days_to_death'],clinical[,'days_to_last_followup'])/365,
  delta  = ifelse(!is.na(clinical[,'days_to_death']),1,0),
  ER     = ifelse(clinical[,'breast_carcinoma_estrogen_receptor_status']=='Positive',1,ifelse(clinical[,'breast_carcinoma_estrogen_receptor_status']=='Negative',0,NA)),
  Age    = clinical[,'age_at_initial_pathologic_diagnosis'],
  Race   = ifelse(clinical[,'race_list']=='BLACK OR AFRICAN AMERICAN','black',ifelse(clinical[,'race_list']=='WHITE','white','other')),
  Gender = ifelse(clinical[,'gender']=='FEMALE','Female',ifelse(clinical[,'gender']=='MALE','Male',NA)),
  Stage  = sapply(clinical[,'stage_event_pathologic_stage'],function(x,pattern='Stage X|Stage IV|Stage [I]*'){ifelse(grepl(pattern,x),regmatches(x,regexpr(pattern,x)),NA)},USE.NAMES = FALSE)
)

# - extract covariates and remove undesiable subjects and NA 
sdata.TCGA <- na.omit(
  sdata.pre[
    sdata.pre[,'yobs'] > 0
    & sdata.pre[,'Age'] <= 75
    & sdata.pre[,'Gender'] == "Female"
    & sdata.pre[,'Race'] %in% c('white')
    & sdata.pre[,'Stage'] %in% c('Stage I','Stage II','Stage III'),
    c('yobs','delta','Age','ER'),
  ]
)
rownames(sdata.TCGA) <- NULL

# - summary statistics of the internal dataset
summary(sdata.TCGA)

## plot a figure to show the existence of a cure fraction
# pdf("Figure_KM_TCGA_BRCA.pdf",width=8.88,height=6.66); {
plot(
  survival::survfit(survival::Surv(yobs, delta) ~ 1, data = sdata.TCGA), 
  conf.int = T, mark.time = TRUE, lwd = 2,
  ylab = "Survival Probability", xlab = "Survival Time (in Years)", 
  xlim = c(0,25), ylim = c(0,1)
)
# }; dev.off()

## fit the model without auxiliary information

# - rescale the Age variable
Age.Min <- min(sdata.TCGA$Age); Age.Max <- max(sdata.TCGA$Age)
sdata.TCGA$Age <- (sdata.TCGA$Age-Age.Min)/(Age.Max-Age.Min)

# - fit the model
set.seed(1)
sol.PHMC <-  SMC.AuxSP(
  formula = Surv(yobs,delta) ~ Age + ER, cureform = ~ Age + ER, 
  sdata = sdata.TCGA, aux = NULL, latency = "PH"
)
print.SMC.AuxSP(object = sol.PHMC)

## fit the model with auxiliary information

# - prepare the auxiliary information
Age.Cut <- c(0,(c(40,50,60)-Age.Min)/(Age.Max-Age.Min),1)
gfunc.t1 <- function(X,Z){
  rbind((X[,1] >= Age.Cut[1] & X[,1] <  Age.Cut[2] & X[,2] == 1),
        (X[,1] >= Age.Cut[2] & X[,1] <  Age.Cut[3] & X[,2] == 1),
        (X[,1] >= Age.Cut[3] & X[,1] <  Age.Cut[4] & X[,2] == 1),
        (X[,1] >= Age.Cut[4] & X[,1] <= Age.Cut[5] & X[,2] == 1),
        (X[,1] >= Age.Cut[1] & X[,1] <  Age.Cut[2] & X[,2] == 0),
        (X[,1] >= Age.Cut[2] & X[,1] <  Age.Cut[3] & X[,2] == 0),
        (X[,1] >= Age.Cut[3] & X[,1] <  Age.Cut[4] & X[,2] == 0),
        (X[,1] >= Age.Cut[4] & X[,1] <= Age.Cut[5] & X[,2] == 0))}
gfunc.t2 <- function(X,Z){rbind((X[,2] == 1), (X[,2] == 0))}
aux <- list( 
  time1 = list(tstar = 5, gfunc = gfunc.t1,
               sprob = c(0.810,0.935,0.925,0.950,0.695,0.780,0.830,0.850)),
  time2 = list(tstar = 10, gfunc = gfunc.t2, sprob = c(0.825,0.705))
)

# - ignore heterogeneity
set.seed(1)
sol.PHMC.Homo <-  SMC.AuxSP(
  formula = Surv(yobs,delta) ~ Age + ER, cureform = ~ Age + ER, 
  sdata = sdata.TCGA, aux = aux, hetero = FALSE, N = 1910, latency = "PH"
)
print.SMC.AuxSP(object = sol.PHMC.Homo)

# - consider heterogeneity
set.seed(1)
sol.PHMC.Hetero <-  SMC.AuxSP(
  formula = Surv(yobs,delta) ~ Age + ER, cureform = ~ Age + ER, 
  sdata = sdata.TCGA, aux = aux, hetero = TRUE, N = 1910, latency = "PH"
)
print.SMC.AuxSP(object = sol.PHMC.Hetero)
```


