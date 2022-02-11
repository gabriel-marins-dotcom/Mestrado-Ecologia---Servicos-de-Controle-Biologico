#PMA.R
#================================================================================
#DEFINE CLASSES
#================================================================================
setOldClass("negbin")
setClass("aphidmodel", representation(model="negbin",parambio="data.frame", GoFglobal="data.frame"))
setClass("control", representation(peakaphid="data.frame",CAD="data.frame"))
setClass("pma", representation(aphidensity="data.frame", mrate="data.frame",gamma="data.frame",aphidmodel="aphidmodel",control="control"))
#================================================================================
#DEFINE FUNCTIONS
#================================================================================
#ANALYSIS
#--------------------------------------------------------------------------------
setGeneric("regtobiol", function(paramreg) {standardGeneric("regtobiol")})
setMethod("regtobiol", signature(paramreg="data.frame"), function(paramreg)
#Description: convert regression parameters into biological parameters, in the Costamagna model
#Input:   paramreg    data.frame with one row, and columns named b0, b1, b2 and N0
#Output:  parambiol   data.frame with one row, and columns named t0, r0, a and N0
{
  b0=paramreg$b0
  b1=paramreg$b1
  b2=paramreg$b2
  N0=paramreg$N0
  
  t0=min(Re(polyroot(c(b0-log(N0),b1,b2))))
  r0=b1+2*b2*t0
  a=-2*b2/r0
  parambiol=data.frame(N0=N0,t0=t0,r0=r0,a=a)
  return(parambiol)  
}) 
#--------------------------------------------------------------------------------
setGeneric("bioltoreg", function(parambiol) {standardGeneric("bioltoreg")})
setMethod("bioltoreg", signature(parambiol="data.frame"), function(parambiol)
#Description: convert regression parameters into biological parameters, in the Costamagna model
#Input:   parambiol   data.frame with one row, and columns named t0, r0, a and N0
#Output:  paramreg    data.frame with one row, and columns named b0, b1, b2 and N0
{
  t0=parambiol$t0
  r0=parambiol$r0
  a=parambiol$a
  N0=parambiol$N0
  
  b0=log(N0)-r0*t0*(1+1/2*a*t0)
  b1=r0*(1+a*t0)
  b2=-1/2*a*r0
  paramreg=data.frame(N0=N0,b0=b0,b1=b1,b2=b2)
  return( paramreg)
})
#--------------------------------------------------------------------------------
setGeneric("propagE_Costamagna", function(varb, param) {standardGeneric("propagE_Costamagna")})
setMethod("propagE_Costamagna", signature(varb="data.frame",param="data.frame"), function(varb,param)
#Description: propagation of error from regression parameters to those of biological parameters, for the Costamagna model
#Input:   varb      variances of regresssion parameters; data.frame with one row, and columns named b0, b1, b2
#         param     values of biological parameters (also includes tmax and Nmax); data.frame with one row, and columns named t0, r0, a, Nmax, tmax
#Output:  var       variance of biological parameters; data.frame with one row, and columns named t0, r0, a, Nmax, tmax
{
  #Extract variance of regression parameters and biological parameter values
  var.b0 = varb$b0
  var.b1 = varb$b1
  var.b2 = varb$b2
  
  t0 = param$t0
  r0 = param$r0
  a = param$a
  Nmax = param$Nmax
  tmax = param$tmax
  
  #Partial derivatives
  t0_b0 = -1/r0
  t0_b1 = -t0/r0
  t0_b2 = -t0^2/r0
  
  r0_b0 = a
  r0_b1 = 1+a*t0
  r0_b2 = 2*t0+a*t0^2
  
  a_b0 = -a^2/r0
  a_b1 = -a/r0*(1+a*t0)
  a_b2 = -1/r0*(2+2*a*t0+a^2*t0^2)
  
  Nmax_b0 = Nmax*(1/(2*a)*r0_b0-r0/(4*a^2)*a_b0)
  Nmax_b1 = Nmax*(1/(2*a)*r0_b1-r0/(4*a^2)*a_b1)
  Nmax_b2 = Nmax*(1/(2*a)*r0_b2-r0/(4*a^2)*a_b2)

  tmax_b0 = t0_b0-1/a^2*a_b0
  tmax_b1 = t0_b1-1/a^2*a_b1
  tmax_b2 = t0_b2-1/a^2*a_b2
  
  #Propagation of uncertainty
  var.t0 = (t0_b0)^2*var.b0 + (t0_b1)^2*var.b1 + (t0_b2)^2*var.b2
  var.r0 = (r0_b0)^2*var.b0 + (r0_b1)^2*var.b1 + (r0_b2)^2*var.b2
  var.a = (a_b0)^2*var.b0 + (a_b1)^2*var.b1 + (a_b2)^2*var.b2

  var.Nmax = (Nmax_b0)^2*var.b0 + (Nmax_b1)^2*var.b1 + (Nmax_b2)^2*var.b2
  var.tmax = (tmax_b0)^2*var.b0 + (tmax_b1)^2*var.b1 + (tmax_b2)^2*var.b2
  
  var=data.frame(t0=var.t0,r0=var.r0,a=var.a,Nmax=var.Nmax,tmax=var.tmax)
  return(var)
}) 
#--------------------------------------------------------------------------------
setGeneric("Costamagna", function(time,parambiol) {standardGeneric("Costamagna")})
setMethod("Costamagna", signature(time="numeric",parambiol="data.frame"), function(time,parambiol)
#Description: return expected aphid density from the Costamagna model (in absence of parasitoids), provided biological parameters (parambiol)
#Input:   time        numeric, time series at which evaluate density 
#         parambiol   data.frame with one row and columns named t0, r0, a and N0
#Output:  aphid       data.frame with one columns named time and mu
  
{
  paramreg=bioltoreg(parambiol)
  b0=paramreg$b0
  b1=paramreg$b1
  b2=paramreg$b2
  N=exp(b0+b1*time+b2*time^2) #N0 is already included in b0
  aphid=data.frame(time=time,mu=N)
  return(aphid)
}) 
#--------------------------------------------------------------------------------
setGeneric("compute_gamma", function(mrate,aphid,time.out) {standardGeneric("compute_gamma")})
setMethod("compute_gamma", signature(mrate="data.frame",aphid="data.frame",time.out="numeric"), function(mrate,aphid,time.out) 
#Description: compute gamma provided average mrate and aphid densities. In case aphid density is zero, there is no change in gamma even if mrate is greater than zero.
#Input:   mrate     average mummification rate per plant per day
#         aphid     average observed aphid density per plant
#         time.out  time series at which to calculate gamma; min and max also define the integral interval
#Output:  gamma  a dataframe with column named time (corresponding to time.out) and gamma
{
  zeroind=which(aphid$mu==0,arr.ind=TRUE) 
  r=mrate$mu/aphid$mu
  r[zeroind]=0 #r is set to zero when aphid density is zero (mrate is then assumed also to be zero)

  t=mrate$time
  I=integrate.data(t,r,x.out=time.out)$int
  
  #Compute gamma
  gamma=exp(-I)
  df=data.frame(time=time.out,gamma=gamma)
  return(df)
})  
#--------------------------------------------------------------------------------
setGeneric("param_aphid", function(modelN,N0=1,refdate=NULL) {standardGeneric("param_aphid")})
setMethod("param_aphid", signature(), function(modelN,N0=1,refdate=NULL)
#Description: extract regression parameters from a model, convert to the biological parameters of Costamagna and compute SE. In case a reference date is provided, also express t0 and tmax in dates. 
#Input:   modelN  negbin model arising from a regression using glmNB_aphid
#         N0      numeric, Known N0 used to define t0
#         refdate Date (optional), date associated to t = 0
#Output:  data.frame with parameters and SE (t0,SE.t0,r0,SE.r0,a,SE.a,tmax,SE.tmax,Nmax,SE.Nmax; all numeric) as columns. If refdate is provided, also return  refdate,tmax.date,t0.date as Date
{
  #Coefficient, variance and df from the regression
  b0=coef(modelN)[[1]] 
  b1=coef(modelN)[[2]]
  b2=coef(modelN)[[3]]
  
  var.b0=summary(modelN)$coefficients[1,2]^2
  var.b1=summary(modelN)$coefficients[2,2]^2
  var.b2=summary(modelN)$coefficients[3,2]^2
  
  #Biological parameters and maximum
  paramreg=data.frame(b0=b0,b1=b1,b2=b2,N0=N0)
  parambio=regtobiol(paramreg)
  t0=parambio$t0
  r0=parambio$r0
  a=parambio$a
  tmax=t0+1/a
  Nmax=N0*exp(r0/(2*a))
  
  #Variance of model parameters and maximum
  param=data.frame(t0=t0,r0=r0,a=a,Nmax=Nmax,tmax=tmax)
  varb=data.frame(b0=var.b0,b1=var.b1,b2=var.b2)
  var.param=propagE_Costamagna(varb=varb,param=param)
  
  var.t0=var.param$t0
  var.r0=var.param$r0
  var.a=var.param$a
  var.tmax=var.param$tmax
  var.Nmax=var.param$Nmax
  
  SE.t0=sqrt(var.t0)
  SE.r0=sqrt(var.r0)
  SE.a=sqrt(var.a)
  SE.tmax=sqrt(var.tmax)
  SE.Nmax=sqrt(var.Nmax)
  
  #Defining the dataframe of parameter values for output
  df=data.frame(t0=t0,SE.t0=SE.t0,r0=r0,SE.r0=SE.r0,a=a,SE.a=SE.a,tmax=tmax,SE.tmax=SE.tmax,Nmax=Nmax,SE.Nmax=SE.Nmax)
  
  #Conversion of tmax and t0 into dates if refdate is provided
  if(!(is.null(refdate)))
  {
    tmax.date=as.Date(refdate, format="%Y-%m-%d")+tmax
    t0.date=as.Date(refdate, format="%Y-%m-%d")+t0
    
    df$refdate=format(refdate, "%b-%d")
    df$tmax.date=format(tmax.date, "%b-%d")
    df$t0.date=format(t0.date, "%b-%d")
  }
  
  return(df)
})
#--------------------------------------------------------------------------------
setGeneric("glmNB_aphid", function(aphid,gamma,refdate=NULL) {standardGeneric("glmNB_aphid")})
#Description: 
#Input:   aphid       individual aphid record per plant; data.frame with time and Np as columns
#         gamma       gamma at same time than provided by aphids (values not repeated); data.frame with time and gamma as columns
#Output:  object of class aphidmodel
#         @model      negbin model; Np is obtained by using proper gamma; N is obtained by setting gamma to 1
#         @parambio   biological parameters associated to the model
#         @GoFglobal  Goodness of fit of the model
setMethod("glmNB_aphid", signature(aphid="data.frame", gamma="data.frame"), function(aphid,gamma,refdate) 
{
  ##Join accurate gamma to aphid
  aphidgamma=merge.float(aphid,gamma,by="time",digits=5)
  
  #Fit the model describing the observed Np (in presence of the parasitoid), i.e. gamma*N, with N (in absence of the parasitoid) described by the model of Costamagna 2008.
  model=glm.nb(Np~1+I(time)+I(time^2)+offset(log(gamma)),data=aphidgamma)
  theta=model$theta

  ##Fit the the null model considering gamma
  model.null=glm(Np~1+offset(log(gamma)),data=aphidgamma,family=negative.binomial(theta=theta))
  
  ##Fit the the full model considering gamma
  degree = as.factor(1:length(aphidgamma$Np))
  model.sat=glm(Np~degree+offset(log(gamma)),data=aphidgamma,family=negative.binomial(theta=theta),control=list(maxit=200))
  
  ##Calculate R2 of the model (R2.McF = McFadden R2 on the whole data set; R2.mean is the adjusted R2 on mean densities)
  GoFglobal=GoF_global(model, model.null, model.sat)

  ##Calculate parameters
  parambio=param_aphid(model,N0=1,refdate=refdate)
  
  #Generate a aphidmodel object to export
  aphidmodel=new("aphidmodel",model=model,parambio=parambio,GoFglobal=GoFglobal)
  
  return(aphidmodel)
})
#--------------------------------------------------------------------------------
setGeneric("pma", function(aphid,mrate,...) {standardGeneric("pma")})
setMethod("pma", signature(aphid="data.frame",mrate="data.frame"), function(aphid,mrate,thresholdperiod=NULL,refdate=NULL,correction=NULL)
#Description: perform the propagated mortality analysis (PMA) from aphid density records and average mrate
#Input:   aphid             individual aphid record per plant; data.frame with time and Np as columns
#         mrate             gamma at same time than provided by aphids (values not repeated); data.frame with time and gamma as columns
#         thresholdperiod   (optional) Range = c(tE1,tE2) over which threshold hold; used for CAD estimates of control
#         refdate=NULL      (optional) Date, corresponding to t = 0; t0 and tmax are expressed also as Date if provided
#         correction=NULL   (optional) numeric, constant correction for mummification rate. Correction = 1/s, with m_overall=1/s*m_measured.
#Output: object of class pma, components of the analysis can be obtain using @mrate(data.frame), @gamma(data.frame), @aphidmodel (class aphidmodel), @control (class control)
{
  #Generate an object of class pma
  pma = new("pma")
  
  #Summarize aphid density and mrate
  pma@aphidensity=ddply(aphid, "time", summarise, mu = mean(Np), sd=sd(Np), se=sem(Np))
  pma@mrate=mrate

  if(is.null(correction)==FALSE){mrate=correction*mrate}
  
  #compute gamma
  pma@gamma=compute_gamma(pma@mrate,pma@aphidensity,time.out=pma@aphidensity$time)

  #Compute aphid model
  pma@aphidmodel=glmNB_aphid(aphid=aphid,gamma=pma@gamma,refdate=refdate)
  
  #Compute the estimate of control
  ##At peak density
  tmax=pma@aphidmodel@parambio$tmax
  Nmax=pma@aphidmodel@parambio$Nmax
  gammapeak=compute_gamma(pma@mrate,aphid=pma@aphidensity,time.out=c(min(pma@aphidensity$time),tmax))[2,2]
  if(is.null(correction)==FALSE){gammapeak=gammapeak^correction}
  
  peakaphid=data.frame(t=tmax,N=Nmax,gamma=gammapeak)
  if(is.null(refdate)==FALSE)
  {
    peakaphid$t.date=pma@aphidmodel@parambio$tmax.date
  }
  pma@control@peakaphid=peakaphid
  
  ##As cumulative aphid density (CAD)
  if(is.null(thresholdperiod)==FALSE)
  {
    #Compute CAD in the presence of parasitoids  
    time=pma@aphidensity$time
    Np=pma@aphidensity$mu
    CADp=integrate.data(time,Np,x.out=thresholdperiod)$int[2]
    
    #Compute CAD in absence of parasitoids
    time=pma@aphidensity$time
    N=Np/pma@gamma$gamma
    CAD=integrate.data(time,N,x.out=thresholdperiod)$int[2]
    
    rCAD=CADp/CAD
    df=data.frame(tE1=thresholdperiod[1], tE2=thresholdperiod[2],CAD=CAD, rCAD=rCAD)
    if(is.null(refdate)==FALSE)
    {
      tE1.date=as.Date(refdate, format="%y-%m-%d")+thresholdperiod[1]
      tE2.date=as.Date(refdate, format="%y-%m-%d")+thresholdperiod[2]
      df$tE1.date=format(tE1.date, "%b-%d")
      df$tE2.date=format(tE2.date, "%b-%d")
    }
    pma@control@CAD=df
  }
  return(pma)
})
#--------------------------------------------------------------------------------
#PLOTS
#--------------------------------------------------------------------------------
setGeneric("plot_aphid", function(pma) {standardGeneric("plot_aphid")})
setMethod("plot_aphid", signature(pma="pma"), function(pma)
#Description: plot N (solid line) and Np (dotted line) predictions from pma, as well as average and SD of experimental Np. If thresholdperiod is provided, the grey zone represent CAD over this period.
#Input: object of class pma
#output: graph (see description)
{
  #Assign objects to be used
  model=pma@aphidmodel@model
  Npobs=pma@aphidensity
  mrate=pma@mrate

  #Estimate the fit of N and Np 
  ti=min(Npobs$time)
  tf=max(Npobs$time)
  t.fit=seq(ti,tf,length.out=100)
  
  gamma.fit=compute_gamma(mrate,Npobs,t.fit)$gamma
  Np.fit=predict(model,data.frame(time=t.fit,gamma=gamma.fit),type="response",se.fit = FALSE)
  N.fit=Np.fit/gamma.fit
  fit=data.frame(t=t.fit,N=N.fit,Np=Np.fit)
  
  #Plot experimental data (mu,sd) as well as the prediction (fit)
  p=ggplot() +
    geom_point(data = Npobs, aes(x = time, y = mu),size=2) +
    geom_errorbar(data = Npobs, aes(x = time, ymin = mu - sd, ymax = mu + sd), width = 0.4) +
    geom_line(data=fit,aes(x = t, y = N)) +
    geom_point(data=fit,aes(x = t, y = Np),shape=16,size=0.5) +
    
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  #Draw threshold period if CAD is provided
  if(dim(pma@control@CAD)[1]!=0)
  {
    ti=min(Npobs$time)
    tE1=pma@control@CAD$tE1
    tE2=pma@control@CAD$tE2
    t.fit=seq(tE1,tE2,length.out=100)
    tiE1=c(ti,tE1)
    gamma0=compute_gamma(mrate,Npobs,tiE1)$gamma[2]
    gamma.fit=gamma0*compute_gamma(mrate,Npobs,t.fit)$gamma
    Np.fit=predict(model,data.frame(time=t.fit,gamma=gamma.fit),type="response",se.fit = FALSE)
    N.fit=Np.fit/gamma.fit
    
    N.fit=c(0,N.fit,0)
    t.fit=c(t.fit[1],t.fit,t.fit[100])
    fit=data.frame(t=t.fit,N=N.fit)
    p$layers <- c(geom_polygon(data=fit,aes(x = t, y = N), fill = gray(0.8)), p$layers)
  }
  return(p)
})
#--------------------------------------------------------------------------------
setGeneric("plot_mrate", function(pma) {standardGeneric("plot_mrate")})
setMethod("plot_mrate", signature(pma="pma"), function(pma) 
#Description: plot mummification rate from a pma object
#Input: object of class pma
#output: graph (see description)
{
  p=ggplot() + geom_point(data = pma@mrate, aes(x = time, y = mu),colour=1,size=2) +
    geom_line(data = pma@mrate, aes(x = time, y = mu)) +
    
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
})
#--------------------------------------------------------------------------------
setGeneric("plot_gamma", function(pma) {standardGeneric("plot_gamma")})
setMethod("plot_gamma", signature(pma="pma"), function(pma) 
#Description: plot gamma from a pma object
#Input: object of class pma
#output: graph (see description)
{
  p=ggplot() +
    geom_point(data = pma@gamma, aes(x = time, y = gamma),colour=1,size=2) +
    geom_line(data = pma@gamma, aes(x = time, y = gamma)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
})
#--------------------------------------------------------------------------------
#UTILS GENERAL
#--------------------------------------------------------------------------------
merge.float = function(x,y,by,digits=5)
#Description: merge two dataframes with a row made of floats, by rounding first at the desired number of digits
{
  #Rounding time for float merging
  x[by]=signif(x[by], digits = digits)
  y[by]=signif(y[by], digits = digits)
  z=merge(x,y,by=by,all=TRUE)
  return(z)
}
#--------------------------------------------------------------------------------
integrate.data = function(x,y,x.out=x)
#Description: fit a cubic spline to data and integrate. 
{
  #Define the cubic spline interpolating data points (x,y)
  S=splinefun(x,y,method="natural")
  
  #Define a function to numerically integrate the spline S
  integrate.value=function(x,lower,upper){integrate(x,lower=lower,upper=upper)$value}
  integrate.vector=Vectorize(FUN=integrate.value,vectorize.args = c("lower","upper"))
  
  #Define x values on which evaluate the integral
  lower=rep(min(x.out),length(x.out))
  upper=x.out
  I=integrate.vector(S,lower=lower,upper=upper)
  return(list(x=x.out,int=I))
}
#--------------------------------------------------------------------------------
#UTILS STATISTICS
#--------------------------------------------------------------------------------
sem = function(x)
#Description: compute the standard error of the mean (SEM) of the vector x
#Argument:  x   vector
#Return:    SEM numeric
{
  SEM=sd(x)/sqrt(length(x))
  return(SEM)
}
#--------------------------------------------------------------------------------
GoF_global = function(model,model.null,model.sat)
#Description:   Compute the r2 dev from Cameron 1996 as well as a Likelihood ratio test (LR.test.p) with H0 assuming equality between the model and the null model
#Arguments:     model       model from glm.nb
#               model.null  null model from glm.nb
#Return:        Dataframe containing R2DEV and the p value from Likelihood ratio test (LR.test.p) 
{
  LR.test.p=anova(model.null,model,test = "Chisq")[["Pr(>Chi)"]][2]
  R2.DEV=1-(model.sat$deviance-model$deviance)/(model.sat$deviance-model.null$deviance)
  return(data.frame(R2.DEV=as.numeric(R2.DEV),LR.test.p=LR.test.p))
}
#--------------------------------------------------------------------------------
