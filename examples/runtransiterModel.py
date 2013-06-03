import numpy as np
from matplotlib import pyplot as plt
import oscaar
from scipy import optimize
from numpy.random import shuffle
import transiterFit
from oscaar.extras.knownSystemParameters import returnSystemParams

#Load in the data from OSCAAR output. 
#sampleData = oscaar.load('sampleOutput/oscaarDataBase.pkl')
#NormFlux=sampleData.lightCurve
#timeObs=sampleData.times
#flux_error=sampleData.lightCurveError

#Planetary Parameters to create fake data. 
fk_RpRs=0.12
fk_aRs=12.7
fk_per=1.58
fk_inc=89.5
fk_t0=2454344.30867
fk_gam1=0.23
fk_gam2=0.50
fk_ecc=0.0
fk_argper=0.0
stddev=0.0015 #Standard Deviation of data. 

#Creating fake dataset, including flux, time, and uncertainty
timeObs,NormFlux=transiterFit.fake_data(stddev,fk_RpRs,fk_aRs,fk_per,fk_inc,fk_t0,fk_gam1,fk_gam2,fk_ecc,fk_argper)
flux_error=stddev*np.ones(np.size(timeObs))

#Initial Guess for Planetary Parameters by pulling them from the exoplanet website.
#For this example we're using GJ 1214 b
RpOverRs= 0.12
aOverRs = 13.0
period = 1.58
inclination = 89.5
epoch = 2454344.30867
gamma1 = 0.23
gamma2 = 0.30
longPericenter = 0.0

#Get parameters from exoplanet.org 
planet = 'GJ 1214 b'
#[p,ap,P,inc,eccentricity] = returnSystemParams.transiterParams(planet)

#RpOverRs = returnSystemParams.exoplanetDB[planet]['PER']
#aOverRs =  ap #returnSystemParams.exoplanetDB[planet]['PER']
#period = P
#inclination = inc
#epoch = timeObs[np.size(timeObs)/2] #Guess at the middle of the dataset
#gamma1 = 0.23 #For now have to look these up. 
#gamma2 = 0.30
#eccentricity = eccentricity
#longPericenter = 0.0

#Run the intial fit with input parameters stated above. 
fit,success = transiterFit.run_LMfit(timeObs,NormFlux,flux_error,RpOverRs,aOverRs,inclination,epoch,
                                     fitLimbDark='quadratic',
                                     plotting=True)

#Run MC fit to estimate uncertainty
n_iter = 800 #Number of random datasets to fit to. 
transiterFit.run_MCfit(n_iter,timeObs,NormFlux,flux_error,fit,success,plotting=True)

