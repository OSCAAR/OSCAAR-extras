

import oscaar
import numpy as np


dataBankPath = "/Users/bmorris/git/OSCAAR-extras/examples/sampleOutput/oscaarDataBase.pkl"
initParams = [0.11,14.1,1.580400,90.0,0.23,0.3,0.0,0.0,2456427.9425593214]
Nsteps = 1e5#1e6
initBeta = np.zeros([4]) + 0.005    ## length = number of free params
idealAcceptanceRate = 0.30
saveInterval = 1e3     ## Must be a number such that Nsteps % saveInterval == 0
burnFraction = 0.20
mcmcinstance = oscaar.fitting.mcmcfit(dataBankPath,initParams,initBeta,Nsteps,\
                saveInterval,idealAcceptanceRate,burnFraction)

mcmcinstance.run(updatepkl=True,plots=True)
mcmcinstance.plot()