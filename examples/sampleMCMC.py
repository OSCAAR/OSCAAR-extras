'''
Created on May 7, 2013

@author: bmmorris
'''
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import oscaar

##############################
# Generate sample transit light curve with random Gaussian noise
x = np.arange(0.1,0.4,5e-3)
noiseFactor = 5e-4      ## Set how much noise is added
modelParams = [0.07,4.1,2.2,83.1,0.23,0.3,0.0,0.0,0.25] ## Set simulated planet's parameters
y = oscaar.occultquad(x,modelParams) + np.random.normal(0,1,x.shape)*noiseFactor
sigma_y = np.zeros_like(x) + noiseFactor
Npoints = len(x)

#plt.errorbar(x,y,yerr=sigma_y)     ## Plot the simulated light curve
#plt.show()

##############################
# Fit with MCMC

Nsteps = int(1e6) ## Number of MCMC "links in the chain" / number of steps
saveInterval = 1e3      ## interval between MCMC links that get saved for uncertainty estimation
inputParams = np.array([0.08,3.9,82.0,0.26],dtype=np.float64)	## Initial fit parameters
                    ## [Rp/Rs, a/Rs, inclination, mid-transit time]

def occult4params(t,params):
    '''Allow 4 parameters to vary freely, keep the others fixed at the values assigned below'''
	p,ap,i,t0 = params
	return oscaar.occultquad(t,[p,ap,2.2,i,0.23,0.3,0.0,0.0,t0])

## Set the initial "beta" vector, discussed in Ford (2005), Section 5
beta = np.zeros_like(inputParams,dtype=np.float64) + np.float64(0.05)
print 'initial betas:',beta
## Tweak the beta vector until it will produce the desired acceptance rate
beta = oscaar.mcmc.optimizeBeta(x,y,sigma_y,inputParams,occult4params,beta)
print 'Optimized betas:',beta

## Run the markov chain monte carlo algorithm, return the best fit parameters,
## the links in each chain, and hte acceptance rate for the chain overall
bestp, allparams, acceptanceRate = oscaar.mcmc.mcmc(x,y,sigma_y,inputParams,occult4params,Nsteps,beta,saveInterval,verbose=True)
print 'acceptanceRate =',acceptanceRate

np.save('allparams.npy',allparams)  ## Save the links in each chain

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
yfit = occult2params(x,bestp)
ax1.errorbar(x,y,yerr=sigma_y,fmt='o-')
ax1.plot(x,yfit,'r')
ax1.set_title("Fit with MCMC")

##############################
# Plot traces and histograms of mcmc params
p = allparams[0,:]
ap = allparams[1,:]
i = allparams[2,:]
t0 = allparams[3,:]
abscissa = np.arange(len(allparams[0,:]))
Nbins = 15
burnFraction = 0.25

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

def histplot(parameter,axis,title):
    n, bins, patches = axis.hist(parameter[burnFraction*len(parameter):len(parameter)], Nbins, normed=0, facecolor='white')
    axis.set_title(title)
histplot(p,ax6,'p')
histplot(ap,ax7,'ap')
histplot(i,ax8,'i')
histplot(t0,ax9,'t0')

#n, bins, patches = ax6.hist(ints[(1-burnFraction)*len(slopes):len(slopes)], Nbins, normed=0, facecolor='white')
#ax6.set_title('Intercepts')

print "MCMC best fit:",bestp

plt.savefig("mcmc_results.png",bbox_inches='tight')     ## Save plot
plt.show()