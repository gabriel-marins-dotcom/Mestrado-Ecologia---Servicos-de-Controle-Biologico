#Examples.R
#================================================================================
#1. NOTES TO THE USERS
#================================================================================
#REQUIREMENTS
#All functions required to perform the propagated mortality analysis are contains in PMA.R (source PMA.R at the beginning of your file).
#The code requires the plyr, MASS and ggplot2 packages to work. 

#EXAMPLES
#In the first part of the example, a dataset with the proper format is generated to perform the analysis. We suggest you to start with this simulated data before analysing your own.
#In the second part, the propagated mortality analysis (PMA) is performed on the later dataset. 
#Two examples are provided, the first one is related to the basics of the function pma(). 
#The second example brings the use of optional parameters including: (1) a threshold period (required to calculate CAD estimate of control), (2) a correction of mummification rates to compensate for biases, and (3) the use of a reference date associated to t = 0.

#We hope this code would facilitate your use of the PMA!
#================================================================================
#2. SOURCES AND PACKAGES
#--------------------------------------------------------------------------------
require(plyr)
require(MASS)
require(ggplot2)

source("PMA.R")
#================================================================================
#3. GENERATE DATA FOR THE EXAMPLES
#================================================================================
#Generate dataset for aphid population dynamics in absence of the studied parasitoids
##Define biological parameters of the Costamagna model 
parambiol=data.frame(N0=1,t0=25,r0=0.25,a=0.03)

##Define the expected population without parasitoids (aphid density following the Costamagna model)
time=seq(0,100,4)
Nmean=Costamagna(time=time, parambiol=parambiol)

#Incorporate the impact of parasitoids
##Define parasitoid mumification rate at the same dates than aphids. We chose here to define control as 1 percent of the expected aphid density per plant per day.
mrate=data.frame(time=Nmean$time,mu=0.01*Nmean$mu)

##Calculate gamma from mrate
gamma=compute_gamma(mrate,Nmean,time.out=Nmean$time)

##Define averages for aphid density in presence of the parasitoid (Np), i.e. observed densities
Npmean=data.frame(time=Nmean$time,mu=gamma$gamma*Nmean$mu)

##Generate a random dataset for the observed aphid density Np, using 25 replicates per dates
numberofdates=length(Npmean$time)
numberofrep=25
time=rep.int(Npmean$time,numberofrep)
Np=rnegbin(n=numberofrep*numberofdates, mu = Npmean$mu, theta=1) # Alternatively with Poisson Distribution: Np=rpois(n=numberofrep*numberofdates,lambda=Npmean$mu), Poisson distribution arises when 1/theta = 0
Npind=data.frame(time=time, Np=Np)
#================================================================================
#4. PERFORMING THE ANALYSIS
#================================================================================
#Example 1: basic use of PMA
##Consider aphid density record (Npind) and the mummification rate (mrate), both known for the same time sequence
Npind #dataframe with columns named time (numeric, usually in days) and Np (numeric, observed aphid density on a plant); one row per individual plant, not the average
mrate #dataframe with columns named time (numeric, usually in days) and mu (numeric, average mummification events per plant per day)

#PMA can be computed from
Analysis=pma(aphid=Npind,mrate=mrate)

##Components of the analysis can be obtained using @, pma() is S4 programmed:
Analysis@aphidensity
Analysis@mrate
Analysis@gamma
Analysis@aphidmodel
Analysis@aphidmodel@model
Analysis@aphidmodel@parambio
Analysis@aphidmodel@GoFglobal
Analysis@control
Analysis@control@peakaphid
Analysis@control@CAD

##Plots can be directly called from the pma object
plot_mrate(Analysis)
plot_gamma(Analysis)
plot_aphid(Analysis)
#--------------------------------------------------------------------------------
#Example 2: optional parameters of PMA
##There are 3 optional parameters that can be passed to pma(), all are set to NULL by default.

##1.THRESHOLD PERIOD (required for CAD estimate of control)
##Define threshold period limits (tE1 and tE2). Define the range over which to calculate the CAD estimates of control.
##When aphids impact on the crop start at their establishment, tE1 might be set as the onset of the experiment tS (here t = 0), but one may want to choose another time.
thresholdperiod=c(0,70) #Alternative example: thresholdperiod=c(40,70)

##2.CORRECTION
##If we underestimated (or overestimated) mummification rate m_measured = s * m_overall, with s a constant. 
##For example, if we underestimated mrate by 25% (s = 0.75)
mrate_under=0.75*mrate

##The correction could be brough directly on mrate before the pma, or specified when calling the pma by using the correction parameter (correction = 1/s)(note: mrate from the pma output would also be corrected)
correction=1/0.75

##3.REFERENCE DATE
##Define a reference date corresponding to t = 0. Automatically convert t0 and tmax into dates.
refdate=as.Date("2017-05-01", format="%Y-%m-%d")

##The call to pma() is as follows:
Analysis=pma(aphid=Npind,mrate=mrate_under,correction=correction,thresholdperiod=thresholdperiod,refdate=refdate)
Analysis

##Plots are called in the same way, but the plot of aphid density represent CAD (in absence of parasitoids) over the threshold period (grey area)
plot_mrate(Analysis)
plot_gamma(Analysis)
plot_aphid(Analysis)
#--------------------------------------------------------------------------------
