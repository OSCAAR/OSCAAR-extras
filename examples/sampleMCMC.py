'''
Created on May 7, 2013

@author: bmmorris
'''
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import oscaar

## If you've run sampleMCMC.py before and want to load the results of that
## run, set loadvals = True. To run a fresh chain, set loadvals = False.
loadvals = False 

if loadvals == False:
    ##############################
    # Generate sample transit light curve with random Gaussian noise
    x = np.arange(0.1,0.4,2e-3)
    noiseFactor = 5.0e-3      ## Set how much noise is added
    modelParams = [0.11,4.1,2.2,83.1,0.23,0.3,0.0,0.0,0.25] ## Set simulated planet's parameters
    y = oscaar.occultquad(x,modelParams) + np.random.normal(0,1,x.shape)*noiseFactor
    sigma_y = np.zeros_like(x) + noiseFactor
    Npoints = len(x)

    plt.errorbar(x,y,yerr=sigma_y)     ## Plot the simulated light curve
    plt.show()

    ##############################
    # Fit with MCMC
    Nsteps = int(0.75e6)     ## Number of MCMC "links in the chain" / number of steps
    saveInterval = 1e3      ## interval between MCMC links that get saved for uncertainty estimation
    inputParams = np.array([0.12,3.9,83.0,0.26],dtype=np.float64)	## Initial fit parameters
    ## Parameters: [Rp/Rs, a/Rs, inclination, mid-transit time]

    def occult4params(t,params):
        '''Allow 4 parameters to vary freely, keep the others fixed at the values assigned below'''
        p,ap,i,t0 = params
        return oscaar.occultquad(t,[p,ap,2.2,i,0.23,0.3,0.0,0.0,t0])

    ## Set the initial "beta" vector, discussed in Ford (2005), Section 5
    beta = np.zeros_like(inputParams,dtype=np.float64) + np.float64(0.05)
    print 'initial betas:',beta
    ## Tweak the beta vector until it will produce the desired acceptance rate
    beta = oscaar.mcmc.optimizeBeta(x,y,sigma_y,inputParams,occult4params,beta,idealAcceptanceRate=0.30)
    print 'Optimized betas:',beta

    ## Run the markov chain monte carlo algorithm, return the best fit parameters,
    ## the links in each chain, and the acceptance rate for the chain overall
    bestp, allparams, acceptanceRate = oscaar.mcmc.mcmc(x,y,sigma_y,inputParams,occult4params,Nsteps,beta,saveInterval,verbose=True)
    print 'acceptanceRate =',acceptanceRate
    np.save('bestp.npy',bestp)
    np.save('allparams.npy',allparams)  ## Save the links in each chain

if loadvals:
    ## If loading old values, set some other variables that would have been set if not loading old values
    allparams = np.load('allparams.npy')  ## Save the links in each chain
    bestp = np.load('bestp.npy')
    x = np.arange(0.1,0.4,5e-3)
    noiseFactor = 5e-3      ## Set how much noise is added
    modelParams = [0.07,4.1,2.2,83.1,0.23,0.3,0.0,0.0,0.25] ## Set simulated planet's parameters
    y = oscaar.occultquad(x,modelParams) + np.random.normal(0,1,x.shape)*noiseFactor
    sigma_y = np.zeros_like(x) + noiseFactor
    Npoints = len(x)
    def occult4params(t,params):
        '''Allow 4 parameters to vary freely, keep the others fixed at the values assigned below'''
        p,ap,i,t0 = params
        return oscaar.occultquad(t,[p,ap,2.2,i,0.23,0.3,0.0,0.0,t0])
##############################
# Prepare figures
fig = plt.figure()
ax1 = fig.add_subplot(331)
ax2 = fig.add_subplot(332)
ax3 = fig.add_subplot(333)
ax4 = fig.add_subplot(334)
ax5 = fig.add_subplot(335)
ax6 = fig.add_subplot(336)
ax7 = fig.add_subplot(337)
ax8 = fig.add_subplot(338)
ax9 = fig.add_subplot(339)
yfit = occult4params(x,bestp)
ax1.errorbar(x,y,yerr=sigma_y,fmt='o-')
ax1.plot(x,yfit,'r')
ax1.set_title("Fit with MCMC")

##############################
# Plot traces and histograms of mcmc params
p = allparams[0,:]
ap = allparams[1,:]
i = allparams[2,:]
t0 = allparams[3,:]
abscissa = np.arange(len(allparams[0,:]))   ## Make x-axis for trace plots
burnFraction = 0.20     ## "burn" or ignore the first 20% of the chains

ax2.plot(abscissa,p,'k.')
ax2.set_title('p trace')
ax2.axvline(ymin=0,ymax=1,x=burnFraction*len(abscissa),linestyle=':')

ax3.plot(abscissa,ap,'k.')
ax3.set_title('ap trace')
ax3.axvline(ymin=0,ymax=1,x=burnFraction*len(abscissa),linestyle=':')

ax4.plot(abscissa,i,'k.')
ax4.set_title('i trace')
ax4.axvline(ymin=0,ymax=1,x=burnFraction*len(abscissa),linestyle=':')

ax5.plot(abscissa,t0,'k.')
ax5.set_title('t0 trace')
ax5.axvline(ymin=0,ymax=1,x=burnFraction*len(abscissa),linestyle=':')

def histplot(parameter,axis,title,bestFitParameter):
    postburn = parameter[burnFraction*len(parameter):len(parameter)]    ## Burn beginning of chain
    Nbins = 15              ## Plot histograms with 15 bins
    n, bins, patches = axis.hist(postburn, Nbins, normed=0, facecolor='white')  ## Generate histogram
    plus,minus = oscaar.mcmc.get_uncertainties(postburn,bestFitParameter)   ## Calculate uncertainties on best fit parameter
    axis.axvline(ymin=0,ymax=1,x=bestFitParameter+plus,ls=':',color='r')    ## Plot vertical lines representing uncertainties
    axis.axvline(ymin=0,ymax=1,x=bestFitParameter-minus,ls=':',color='r')        
    axis.set_title(title)
## Plot the histograms
histplot(p,ax6,'p',bestp[0])
histplot(ap,ax7,'ap',bestp[1])
histplot(i,ax8,'i',bestp[2])
histplot(t0,ax9,'t0',bestp[3])

plt.savefig("mcmc_results.png",bbox_inches='tight')     ## Save plot
plt.show()