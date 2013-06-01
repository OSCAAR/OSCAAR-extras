import numpy as np
from matplotlib import pyplot as plt
import oscaar
from scipy import optimize
from numpy.random import shuffle

#Work below is from Nolan Matthews while using some of Brett Morris's
#code, particularly the transitModel.py,simulatedLightCurve.py, and
#modelLightCurve.py scripts.

#The script below allows two functionalities. There are some functions
#defined that allow the user to create a fake dataset. These were used
#to develop a random MC fitting routine which relies on the scipy
#function optimize.curve_fit. There is a tutorial on how to use the functions
#which can found at XXX. The optimize.curve_fit tool uses a least square 
#Levenburg-Marquadt (sp?) algorithm. Caution should be taken to the inital
#guesses for the parameters as least sq. LM fitting typically will find
#a local minimum. One can find archival results using XXX.

#Make Fake Datasets using random number generator to test fitting function
def fake_data(stddev,RpRs,aRs,per,inc,midtrantime,gamma1,gamma2,ecc,argper):
    
    #Define Times (in days) centered at mid-transit time. 
    expTime = 45./(3600*24.) #Set to be 1-minute, somewhat typical for observing.
    Nimages = 200. #Needs to be long enough to cover entire transit w/ some baseline. 
    times = np.arange(start=midtrantime-expTime*Nimages/2.,
                      stop=midtrantime+expTime*Nimages/2.,
                      step=expTime)

    #Creates Gaussian Distributed Data using numpy.random.normal function based on standard dev. 
    
    #Commented out stuff - talk to Brett about input parameters for model.
    #modelParams = [RpRs,aRs,per,inc,gamma1,gamma2,ecc,argper,midtrantime]
    #fake_data = oscaar.occultquad(times,modelParams) + np.random(scale=stddev,size=np.size(times))
    
    #Uses alternate input parameters setup for occultquad.
    perfect_data = oscaar.transitModel.occultquad(times,RpRs,aRs,inc,midtrantime)
    random_dist = np.random.normal(scale=stddev,size=np.size(times))
    fk_data = perfect_data + random_dist
    
    return times,fk_data

#Runs the intial fit using the LM least sq. algorithm. 
def run_LMfit(timeObs,NormFlux,flux_error,plotting=True):

    fit,success=optimize.curve_fit(oscaar.transitModel.occultquad,
                                   xdata=timeObs,
                                   ydata=NormFlux,
                                   p0=(RpOverRs,aOverRs,inclination,epoch),
                                   sigma=flux_error,
                                   maxfev=100000,
                                   xtol=2e-15,
                                   ftol=2e-16)

    #Check for Convergence    
    if type(success) != np.ndarray:
        print "The inital fit was not able to converge. Check to see if the input parameters are accurate."
        print ""
        
    #If Convergence is True, look at the results to double check.
    else:
        print "Results from the inital fit w/ uncertainties based on the sq. root of the covariance matrix"
        params = ["Rp/Rs","a/Rs","inc","Mid-Tran Time"]
        for i in range(0,np.size(fit)):
            print params[i],fit[i],"+/-",np.sqrt(success[i][i])
            
    #Visually check to see if it's reasonable
    if plotting == True:
        plt.plot(timeObs,oscaar.transitModel.occultquad(timeObs,fit[0],fit[1],fit[2],fit[3]))
        plt.plot(timeObs,NormFlux,'o')
        plt.show()
        plt.close()
    
    return fit,success
    
#Run Prayer-Bead or Random Markov Chain to estimate uncertainties.

#Shuffle Function, had to be modded for data type reasons.
def shuffle_func(x):
    shuffle(x)
    return x

#Function that allows one to determine model uncertainties using a random
#Monte Carlo method. 
def run_MCfit(n_iter,timeObs,NormFlux,flux_error,fit,success,plotting=True):

    modelOut  = oscaar.transitModel.occultquad(timeObs,fit[0],fit[1],fit[2],fit[3])
    residuals = NormFlux - modelOut

    RpFit,aRsFit,incFit,epochFit = fit[0],fit[1],fit[2],fit[3]

    #Generate random datasets based on residuals from inital fit. 
    n_sets = n_iter
    Rp,aRs,inc,mid=[],[],[],[]
    for i in range(0,n_sets):
    
        #Randomly shuffling both data/uncertainties together 
        MCset,randSet,SigSet = [],[],[]
        index_shuf = range(len(residuals))
        shuffle(index_shuf)
        for i in index_shuf:
            MCset.append(residuals[i])
            SigSet.append(flux_error[i])
        
        #Generate random dataset and fit to the function.
        randSet = MCset + modelOut
        fit,success=optimize.curve_fit(oscaar.transitModel.occultquad,xdata=timeObs,
                                   ydata=randSet,
                                   p0=(RpFit,aRsFit,incFit,epochFit),
                                   maxfev=10000,
                                   sigma=SigSet,
                                   xtol=2e-15,
                                   ftol=2e-16)
        
        #Save output parameters from fit
        Rp.append(fit[0])
        aRs.append(fit[1])
        inc.append(fit[2])
        mid.append(fit[3])
        
        if plotting == True: 
            plt.plot(timeObs,oscaar.transitModel.occultquad(timeObs,fit[0],fit[1],fit[2],fit[3]))
        
    #Visually compare MC fits to inital fit and observational data.
    if plotting == True:
        plt.errorbar(timeObs,NormFlux,yerr=flux_error,linestyle='None',marker='.',label="Data")
        plt.plot(timeObs,oscaar.transitModel.occultquad(timeObs,RpFit,aRsFit,incFit,epochFit),lw=3.0,color='k',label="Inital Fit")
        plt.title('Results from Random MC Fits')
        plt.xlabel('JD (days)')
        plt.ylabel('Normalized Flux')
        plt.legend()
        plt.show()
        plt.close()
        
    print "Results from random MC fit . . . . . "     
    print "Planetary to Stellar Radius: ",np.mean(Rp),"+/-",np.std(Rp)
    print "Semi-major Axis to Stellar Radius: ",np.mean(aRs),"+/-",np.std(aRs)
    print "Inclination of Orbit: ",np.mean(inc),"+/-",np.std(inc)
    print "Mid-Transit Time [JD]: ",np.mean(mid),"+/-",np.std(mid)
    
    