import numpy as np
from matplotlib import pyplot as plt
import oscaar
from oscaar import transiter
from scipy import optimize
from numpy.random import shuffle
#import random

sampleData = oscaar.load('sampleOutput/oscaarDataBase.pkl')

#Fitting to the Data
NormFlux=sampleData.lightCurve
timeObs=sampleData.times
flux_error=sampleData.lightCurveError

#Orbital Parameters for GJ1208
## Define system parameters for planet GJ 
RpOverRs = 0.13			## R_p/R_s = ratio of planetary radius to stellar radius
aOverRs = 15.0                  ## a/R_s = ratio of semimajor axis to stellar radius
period = 1.58			## Period [days]
inclination = 88.0		## Inclination [degrees]
epoch = timeObs[np.size(timeObs)/2]	## Mid transit time [Julian date], guess is at middle of the dataset. 
gamma1 = 0.23			## Linear limb-darkening coefficient
gamma2 = 0.45			## Quadratic limb-darkening coefficient
eccentricity = 0.0		## Eccentricity
longPericenter = 0.0		## Longitude of pericenter

modelParams = [RpOverRs,aOverRs,period,inclination,gamma1,gamma2,eccentricity,longPericenter,epoch]

#Look at the value of the impact parameter and make sure it's reasonable.
b = aOverRs * np.cos(inclination*np.pi/180.)
print "Impact Parameter: ",b
print "If the value of the impact parameter is close to or greater than 1, then you are inputting an inital guess of a grazing transit. This may lead the fitting function to a unaccurate fit. I suggest intial guesses of b = 0.5 or so should give enough room to find a more accurate answer."  

#Duration of Transit
durationObs = timeObs[np.size(timeObs)-1]-timeObs[0]

#Running the fit
fit,success=optimize.curve_fit(oscaar.transitModel.occultquad,xdata=timeObs,
                               ydata=NormFlux,
                               p0=(RpOverRs,aOverRs,inclination,epoch),
                               sigma=sampleData.lightCurveError,
                               maxfev=10000000,
                               xtol=2e-15,
                               ftol=2e-16)

#Look at the results
for i in range(0,np.size(fit)):
    print fit[i],np.sqrt(success[i][i])

#Visually check to see if it's reasonable
plt.plot(timeObs,oscaar.transitModel.occultquad(timeObs,fit[0],fit[1],fit[2],fit[3]))
plt.plot(timeObs,NormFlux,'o')
plt.show()
plt.close()

#Run Prayer-Bead or Random Markov Chain to estimate uncertainties.
def shuffle_func(x):
    shuffle(x)
    return x

#Generate residuals
modelOut  = oscaar.transitModel.occultquad(timeObs,fit[0],fit[1],fit[2],fit[3])
residuals = NormFlux - modelOut

RpFit,aRsFit,incFit,epochFit = fit[0],fit[1],fit[2],fit[3]

#Generate random datasets based on residuals from inital fit. 
n_sets = 1000
Rp,aRs,inc,mid=[],[],[],[]
MCset,randSet = [],[]
for i in range(0,n_sets):
    MCset = shuffle_func(residuals)
    randSet = MCset + modelOut
    fit,success=optimize.curve_fit(oscaar.transitModel.occultquad,xdata=timeObs,
                               ydata=randSet,
                               p0=(RpFit,aRsFit,incFit,epochFit),
                               maxfev=10000000,
                               sigma=flux_error,
                               xtol=2e-15,
                               ftol=2e-16)
    
    #Save output parameters from fit
    Rp.append(fit[0])
    aRs.append(fit[1])
    inc.append(fit[2])
    mid.append(fit[3])
    
    #Plotting output fits for a visual check
    plt.plot(timeObs,oscaar.transitModel.occultquad(timeObs,fit[0],fit[1],fit[2],fit[3]))

plt.errorbar(timeObs,NormFlux,yerr=flux_error,linestyle='None',marker='.',label="Data")
plt.plot(timeObs,oscaar.transitModel.occultquad(timeObs,RpFit,aRsFit,incFit,epochFit),lw=3.0,color='k',label="Inital Fit")
plt.title('Results from Random MC Fits')
plt.xlabel('JD (days)')
plt.ylabel('Normalized Flux')
plt.legend()
plt.show()

print "Planetary to Stellar Radius: ",np.mean(Rp),np.std(Rp)
print "Semi-major Axis to Stellar Radius: ",np.mean(aRs),np.std(aRs)
print "Inclination of Orbit: ",np.mean(inc),np.std(inc)
print "Mid-Transit Time [JD]: ",np.mean(mid),np.std(mid)
