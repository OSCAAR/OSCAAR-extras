import numpy as np
from matplotlib import pyplot as plt
import oscaar
from scipy import optimize
from numpy.random import shuffle

#Load in the data from OSCAAR output. 
sampleData = oscaar.load('sampleOutput/oscaarDataBase.pkl')
NormFlux=sampleData.lightCurve
timeObs=sampleData.times
flux_error=sampleData.lightCurveError

#Orbital Parameters for GJ1208
#RpOverRs = 0.117		## R_p/R_s = ratio of planetary radius to stellar radius
#aOverRs = 14.7                  ## a/R_s = ratio of semimajor axis to stellar radius
#period = 1.58			## Period [days]
#inclination = 89.0		## Inclination [degrees]
#epoch = timeObs[np.size(timeObs)/2]+0.001 # Mid transit time [Julian date], guess is at middle of the dataset. 
#gamma1 = 0.23			## Linear limb-darkening coefficient
#gamma2 = 0.30			## Quadratic limb-darkening coefficient
#eccentricity = 0.0		## Eccentricity
#longPericenter = 0.0		## Longitude of pericenter

## Initial Guesses for Planetary Parameters  
RpOverRs = 0.117
aOverRs = 14.7
period = 1.58
inclination = 89.5
epoch = timeObs[np.size(timeObs)/2]
gamma1 = 0.23
gamma2 = 0.30
eccentricity = 0.0
longPericenter = 0.0

#Parameters for making fake data for a fake planet
fk_RpRs=0.17
fk_aRs=10.7
fk_per=1.58
fk_inc=89.9
fk_t0=timeObs[np.size(timeObs)/2]
fk_gam1=0.23
fk_gam2=0.30
fk_ecc=0.0
fk_argper=0.0
stddev=0.001 #Standard Deviation of data. 

#Creating fake dataset, including flux, time, and uncertainty
#timeObs,NormFlux=fake_data(stddev,fk_RpRs,fk_aRs,fk_per,fk_inc,fk_t0,fk_gam1,fk_gam2,fk_ecc,fk_argper)
#flux_error=stddev*np.ones(np.size(timeObs))

modelParams = [RpOverRs,aOverRs,period,inclination,gamma1,gamma2,eccentricity,longPericenter,epoch]

#Look at the value of the impact parameter and make sure it's reasonable.
b = aOverRs * np.cos(inclination*np.pi/180.)
print "Impact Parameter: ",b  
print ""
print "Care should be taken to make sure the input parameters reflect archival results."
