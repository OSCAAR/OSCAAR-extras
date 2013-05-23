import numpy as np
from matplotlib import pyplot as plt
import oscaar
from oscaar import transiter

sampleData = oscaar.load('sampleOutput/oscaarDataBase.pkl')

#Fitting to the Data
NormFlux=sampleData.lightCurve
times=sampleData.times

#Orbital Parameters for Observed Exoplanet
Period = 1.580 #days
Rp = 0.11
aRstar = 14.0
inc = 89.0
dt = times[2]-times[1]
ecc = 0.0
arg_per = 0.0

#Running the fit 
fit,success=transiter.fittransit(NormFlux,Rp,aRstar,inc,dt,Period)

#Pulling output parameters from model parameters
Flux,Rp,aRstarFIT,incFIT,aRstar,midtrantime=transiter.output_params(times,NormFlux,fit,success,Period,ecc,arg_per)

