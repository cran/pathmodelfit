#' Compute fit indices for the path component of latent variable structural equation models.
#'
#' \code{pathmodelfit} computes fit indices for evaluating the path component of latent variable structural equation models. Available fit indices include RMSEA-P and NSCI-P originally presented and evaluated by Williams and O'Boyle (2011) and demonstrated by O'Boyle and Williams (2011) and Williams, O'Boyle, & Yu, (2019). Also included are fit indices described by Hancock and Mueller (2011).
#'
#' @param lavaanoutput A \code{lavaan} \code{sem} object.
#' @return A vector with RMSEA-P, a p-value for the chi-square test comparing the theoretical and saturated model, a 90 percent confidence interval for RMSEA-P, NSCI-P, and SRMRs, RMSEAs, TLIs, and CFIs.
#' @import lavaan
#' @export
#' @references
#' Hancock, G. R., & Mueller, R. O. (2011). The reliability paradox in assessing structural relations within covariance structure models. Educational and Psychological Measurement, 71(2), 306-324.
#'
#' McNeish, D., & Hancock, G. R. (2018). The effect of measurement quality on targeted structural model fit indices: A comment on Lance, Beck, Fan, and Carter (2016). Psychological Methods, 23(1), 184–190. https://doi.org/10.1037/met0000157
#'
#' O'Boyle, E. H., Jr., & Williams, L. J. (2011). Decomposing model fit: Measurement vs. theory in organizational research using latent variables. Journal of Applied Psychology, 96(1), 1–12. https://doi.org/10.1037/a0020539
#'
#' Williams, L. J., & O’Boyle, E. H. (2011). The myth of global fit indices and alternatives for assessing latent variable relations. Organizational Research Methods, 14, 350-369.
#'
#' Williams, L. J., O’Boyle, E. H., & Yu, J. (2020). Condition 9 and 10 tests of model confirmation: A review of James, Mulaik, and Brett (1982) and contemporary alternatives. Organizational Research Methods, 23, 1, 6-29.
#'
#' @examples
#' library(lavaan)
#'
#' model4 <- '
#' Ldrrew =~ LdrrewI1 + LdrrewI2 + LdrrewI3
#' Jobcom =~ JobcomI1 + JobcomI2 + JobcomI3
#' Jobsat =~ JobsatI1 + JobsatI2 + JobsatI3
#' Orgcom =~ OrgcomI1 + OrgcomI2 + OrgcomI3
#' Jobsat ~ Ldrrew + Jobcom
#' Orgcom ~ Jobsat'
#'
#' data(mediationVC)
#'
#' fit <- sem(model4, sample.cov = mediationVC, sample.nobs = 232)
#' pathmodelfit(fit)
#'

pathmodelfit<-function(lavaanoutput){
  #extracting data features
  fit<-lavaanoutput
  N<-fit@loglik$ntotal
  model<-eval(fit@call$model)
  #setting structural paths from fit to 0
  structuralpaths<-which(fit@ParTable$op=='~')
  ptab0<-fit@ParTable
  ptab0$est<-NULL
  ptab0$se<-NULL
  ptab0$free[structuralpaths]<-0

  #updating lavaan call to use ptab0 model
  tmpnullcall<-as.list(fit@call)
  tmpnullcall$model <- as.name('ptab0')
  nullcall <- as.call(tmpnullcall)
  fit0 <- eval(nullcall)

  #saturated structural model
  eqpos<-gregexpr("\n",model)[[1]]
  neqs<-length(eqpos)
  totchars<-nchar(model)
  eqpos<-c(eqpos[1:neqs],totchars)
  modeqs<-character(neqs)
  measurement<-logical(neqs)
  covariance<-logical(neqs)

  for(i in 1:neqs){
    tmp<-substr(model,start=eqpos[i],stop=eqpos[i+1])
    measurement[i]<-grepl('=~',tmp)
    covariance[i]<-grepl('~~',tmp)
    modeqs[i]<-tmp
  }
  structural<-measurement==F & covariance==F
  modsat<-paste0(modeqs[!structural],collapse = '')

  tmpsatcall<-as.list(fit@call)
  tmpsatcall$model <- as.name('modsat')
  satcall <- as.call(tmpsatcall)
  fit1 <- eval(satcall)

  #compute RMSEA-P
  X2ss<-fit1@test[[1]]$stat
  dfss<-fit1@test[[1]]$df

  X2t<-fit@test[[1]]$stat
  dft<-fit@test[[1]]$df
  RMSEAP<-sqrt(((X2t-X2ss)-(dft-dfss))/((dft-dfss)*(N-1)))

  #RMSEAP CI
  #computations from Williams & O'Boyle excel spreadsheet
  B2<-N
  E2<-X2t
  C2<-X2ss
  F2<-dft
  D2<-dfss
  G2<-E2-C2
  H2<-F2-D2
  N2<-G2-H2
  O2<-H2+N2
  P2<-2*H2+4*N2
  Q2<-8*H2+24*N2
  R2<-1-((O2*Q2)/(3*(P2^2)))
  S2<-(P2/(2*O2))-(Q2/(4*P2))
  T2<-1+(R2/O2)*(((P2*(R2-1))/(2*O2))+S2)
  U2<-(R2^2)*P2/(O2^2)
  V2<-sqrt(U2)
  W2<-T2-(V2*1.645)
  X2<-T2+(V2*1.645)
  Y2<-O2*(W2^(1/R2))-S2
  Z2<-O2*(X2^(1/R2))-S2
  AA2<-Y2-H2
  AB2<-Z2-H2
  L2<-sqrt(AA2/(H2*B2-1))
  M2<-sqrt(AB2/(H2*B2-1))

  #computing p-value for rmsea-p model comparison
  pvalue <- anova(fit,fit1)$"Pr(>Chisq)"[2]

  #compute NSCI-P, which is CFI-P
  X2sn<-fit0@test[[1]]$stat
  dfsn<-fit0@test[[1]]$df
  NSCIP<-((X2sn-X2t)-(dfsn-dft))/((X2sn-X2ss)-(dfsn-dfss))
  output<-c(RMSEAP,pvalue,L2,M2,NSCIP)
  labs<-c('RMSEA-P','Theo. vs. Sat. p-value','RMSEA-P 90% lower bound','RMSEA-P 90% upper bound','NSCI-P')

  #compute Hancock and Mueller
  #extracting the implied latent variable vc matrix.
  #this will be a list if there are multiple groups
  impliedlvvcmatrix<-lavTech(fit1,what='cov.lv')
  latenteqs<-(fit1@ParTable$op=="=~")
  lvnames<-unique(fit1@ParTable$lhs[latenteqs])
  for (i in seq_along(impliedlvvcmatrix)){
    colnames(impliedlvvcmatrix[[i]]) <- lvnames
    rownames(impliedlvvcmatrix[[i]]) <- lvnames
  }
  modlvstructure<-paste0(modeqs[structural],collapse = '')

  tmpimpliedcall<-as.list(fit@call)
  tmpimpliedcall$model <- as.name('modlvstructure')

  #check whether 'data' or 'sample.cov' is used
  name_eq_data<-(names(tmpimpliedcall) %in% 'data')
  check<-max(name_eq_data)

  #if a data frame was specified, change 'data' to 'sample.cov' to input implied vc
  if(check==1){
    names(tmpimpliedcall)[which(name_eq_data==T)] <- 'sample.cov'
  }
  tmpimpliedcall$sample.cov <- as.name('impliedlvvcmatrix')

  #specifysample size, which may be a list for multi group application
  Samplesize<-fit@Data@nobs
  tmpimpliedcall$sample.nobs <- as.name('Samplesize')
  impliedcall <- as.call(tmpimpliedcall)
  fitwithimplied <- eval(impliedcall)

  Hancock<-fitMeasures(fitwithimplied,fit.measures = c('srmr','rmsea','tli','cfi'))
  names(Hancock)<-paste0(names(Hancock),'.s')

  Est<-c(output,Hancock)
  foutput<-as.data.frame(Est)
  rownames(foutput)<-c(labs,names(Hancock))
  print(foutput,digits=4)
}
