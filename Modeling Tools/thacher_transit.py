import numpy as np
import matplotlib.pyplot as plt
import sys,os,re,pdb,time,gc,glob
import constants as c
import scipy as sp
from scipy.stats.kde import gaussian_kde
from length import length
from astropy.stats import mad_std
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation
import batman,gc,time,emcee
from scipy.optimize import minimize
import phot_pipe as pp
import pandas as pd

# Target information
targname = 'TIC398733009.01'
obstype = 'Transit'
rastr = '02:08:40.74'
decstr = '-09:02:45.8'
date = '20191019'

# Path and date information
paths = pp.get_paths(targname=targname,obstype=obstype)
dates = pp.get_dates(targname=targname)

# Coordinates 
coords = SkyCoord(rastr,decstr,unit=(u.hour,u.deg))
ra  = coords.ra.deg
dec = coords.dec.deg

# Ephemeris
ingress = 2458775.821
midtrans = 2458775.841
egress = 2458775.861
period = 0.818605

# Data from data file
data = pd.read_csv(paths['output']+date+'/AAS_398733009-01_20191019_measurements.csv')

# Wrangle data. Use plotting to help
t = data['BJD_TDB'].values
f = data['rel_flux_T1'].values
inds, = np.where(f < 0.1)
t = t[inds][:-31]
f = f[inds][:-31]
plt.ion()
plt.figure(1)
plt.clf()
plt.plot(t,f,'ko')
plt.plot(t[0:8],f[0:8],'ro')
plt.plot(t[-30:],f[-30:],'ro')

# Cast in terms of expected mid-transit time
tfit = np.append(t[0:8],t[-30:])
ffit = np.append(f[0:8],f[-30:])
tmed  = np.median(tfit)

# Fit out of transit trend
fit,cov = np.polyfit(tfit-tmed,ffit,2,cov=True)
tplt = np.linspace(np.min(tfit-tmed),np.max(tfit-tmed),1000)
fplt = np.polyval(fit,tplt)
plt.plot(tplt+tmed,fplt)

# Renormalize
fnorm = np.polyval(fit,t-tmed)
f = f/fnorm

# Estimate errors from STD of out of transit data
ferr = np.ones(len(t))*np.std(f[-30:],ddof=1)

# Exposure time (for smoothing)
exp_time = 180./(86400.)

# Initial transit parameters (might have to guess with some)
param_info = {'t0':midtrans,'P':period ,'rprs':np.sqrt(0.074) ,
              'aors':6.8 , 'inc':84.0 ,'ecc':0.0, 'omega':0.0,
              'q1':0.3 ,'q2':0.3 ,'logf':-10.0,
              'exp_time':exp_time, 'supersample_factor':7,
              'limb_model':'quadratic'}

# Parameter "vector"
theta = [param_info['t0'],param_info['P'],param_info['rprs'],param_info['aors'],
         param_info['inc'],param_info['ecc'],param_info['omega'],param_info['q1'],param_info['q2'],
         param_info['logf']]

# Create a model from the initial parameters and see how it matches up. 
tmodel = np.linspace(np.min(t),np.max(t),100000)
model = compute_trans(theta, tmodel, param_info)
plt.ion()
plt.figure(2)
plt.clf()
plt.errorbar(t,f,yerr=ferr,fmt='o',zorder=1)
plt.plot(tmodel,model,'k-',zorder=2)
plt.savefig(paths['output']+date+'/Preliminary_fit.png',dpi=300)


'''
targname = 'TIC125442121.01'
obstype = 'Transit'
rastr = '23 26 19.13'
decstr = '+38 33 10.7'
date = '20191115'

paths = pp.get_paths(targname=targname,obstype=obstype)
dates = pp.get_dates(targname=targname)

coords = SkyCoord(rastr,decstr,unit=(u.hour,u.deg))
ra  = coords.ra.deg
dec = coords.dec.deg

ingress = 2458802.66963
midtrans = 2458802.74463
egress = 2458802.81963

data = pd.read_csv(paths['output']+date+'/TIC125442121-01_20191115.csv')

t = data['BJD_TDB'].values
f = data['rel_flux_T1'].values
ferr = np.ones(len(t))*np.std(f[-75:],ddof=1)
ferr /= np.median(f[-70:])
f /= np.median(f[-70:])

exp_time = 20./(86400.)


param_info = {'t0':midtrans,'P':3.01963 ,'rprs':np.sqrt(0.013178085938) ,
              'aors':6.8 , 'inc':90.0 ,'ecc':0.0, 'omega':0.0,
              'q1':0.3 ,'q2':0.3 ,'logf':-10.0,
              'exp_time':exp_time, 'supersample_factor':7,
              'limb_model':'quadratic'}

theta = [param_info['t0'],param_info['P'],param_info['rprs'],param_info['aors'],
         param_info['inc'],param_info['ecc'],param_info['omega'],param_info['q1'],param_info['q2'],
         param_info['logf']]

tmodel = np.linspace(np.min(t),np.max(t),100000)
model = compute_trans(theta, tmodel, param_info)
plt.ion()
plt.figure(1)
plt.clf()
plt.errorbar(t,f,yerr=ferr,fmt='o',zorder=1)
plt.plot(tmodel,model,'k-',zorder=2)
plt.savefig(paths['output']+date+'/Preliminary_fit.png',dpi=300)
'''


'''
targname = 'TIC125442121.01'
obstype = 'Transit'
rastr = '23 26 19.13'
decstr = '+38 33 10.7'
date = '20191112'

paths = pp.get_paths(targname=targname,obstype=obstype)
dates = pp.get_dates(targname=targname)

coords = SkyCoord(rastr,decstr,unit=(u.hour,u.deg))
ra  = coords.ra.deg
dec = coords.dec.deg

ingress = 2458799.650
midtrans = 2458799.725
egress = 2458799.800

data = pd.read_csv(paths['output']+date+'/Final_LC_r.csv')

t = data['BJD'].values
f = data['flux'].values
ferr = data['flux_err'].values

exp_time = 20./(86400.)


param_info = {'t0':2458799.725,'P':3.01963 ,'rprs':np.sqrt(0.013178085938) ,
              'aors':6.8 , 'inc':90.0 ,'ecc':0.0, 'omega':0.0,
              'q1':0.3 ,'q2':0.3 ,'logf':-10.0,
              'exp_time':exp_time, 'supersample_factor':7,
              'limb_model':'quadratic'}


theta = [param_info['t0'],param_info['P'],param_info['rprs'],param_info['aors'],
         param_info['inc'],param_info['ecc'],param_info['omega'],param_info['q1'],param_info['q2'],
         param_info['logf']]

tmodel = np.linspace(np.min(t),np.max(t),100000)
model = compute_trans(theta, tmodel, param_info)
plt.clf()
plt.errorbar(t,f,yerr=ferr,fmt='o',zorder=1)
plt.plot(tmodel,model,'k-',zorder=2)
plt.savefig(paths['output']+date+'/Preliminary_fit.png',dpi=300)
'''

'''
targname = 'TESS'
obstype = 'Tabby'
rastr = '20 06 15.45'
decstr = '+44 27 24.79'
date = 'TESS'

paths = pp.get_paths(targname=targname,obstype=obstype)
#dates = pp.get_dates(targname=targname)

coords = SkyCoord(rastr,decstr,unit=(u.hour,u.deg))
ra  = coords.ra.deg
dec = coords.dec.deg

ingress = 2458729.8
midtrans = 2458730.15
egress = 2458730.55

data = pd.read_csv(paths['output']+'TESS_2min_S15.csv')

t = data['time'].values + 2457000.0
f = data['flux'].values
ferr = data['flux_err'].values
ferr /= np.median(f)
f /= np.median(f)

exp_time = 120./(86400.)

P = 1.31*365.25 # Estimate from 20h duration and stellar density of 0.36 solar
# Other parameters estimated from Winn 2010 eq 19

param_info = {'t0':midtrans,'P':P,'rprs':np.sqrt(0.01) ,
              'aors':200.0, 'inc':90.0 ,'ecc':0.0, 'omega':0.0,
              'q1':0.3 ,'q2':0.3 ,'logf':-10.0,
              'exp_time':exp_time, 'supersample_factor':7,
              'limb_model':'quadratic'}


theta = [param_info['t0'],param_info['P'],param_info['rprs'],param_info['aors'],
         param_info['inc'],param_info['ecc'],param_info['omega'],param_info['q1'],param_info['q2'],
         param_info['logf']]

tmodel = np.linspace(np.min(t),np.max(t),100000)
model = compute_trans(theta, tmodel, param_info)
plt.clf()
plt.errorbar(t,f,yerr=ferr,fmt='o',zorder=1)
plt.plot(tmodel,model,'k-',zorder=2)
plt.xlim(ingress-0.3,egress+0.3)
plt.savefig(paths['output']+'Preliminary_fit.png',dpi=300)
'''

# Use only data right around transit
inds, = np.where((t >= ingress-0.3)&(t <= egress+0.3))
t = t[inds]
f = f[inds]
ferr = ferr[inds]



def compute_trans(theta, t, param_info):

    """
    ----------------------------------------------------------------------
    compute_trans:
    --------------
    Function to compute transit curve for given transit parameters. Returns
    model flux. 

    """

    t0,P,rprs,aors,inc,ecc,omega,q1,q2,logf = theta
    
    params = batman.TransitParams()       #object to store transit parameters
    params.t0 = t0                        #time of inferior conjunction
    params.per = P                        #orbital period
    params.rp = rprs                      #planet radius (in units of stellar radii)
    params.a = aors                       #semi-major axis (in units of stellar radii)
    params.inc = inc                      #orbital inclination (in degrees)
    params.ecc = ecc                      #eccentricity
    params.w = omega                      #longitude of periastron (in degrees)
    params.limb_dark = param_info['limb_model']         #limb darkening model
    u1,u2 = qtou(q1,q2,limb=param_info['limb_model'])
    params.u = [0.5, 0.5]                 #limb darkening coefficients [u1, u2, u3, u4]
    
    m = batman.TransitModel(params, t,
                            supersample_factor = param_info['supersample_factor'],
                            exp_time = param_info['exp_time'])    #initializes model

    model = m.light_curve(params)

    return model



def residuals(theta, t, f, param_info):
    """ 
    ----------------------------------------------------------------------
    resuduals:
    ----------
    Calculate residuals given an input model.
    ----------------------------------------------------------------------
    """

    t0,P,rprs,aors,inc,ecc,omega,q1,q2,logf = theta
    
    model = compute_trans(theta,t, param_info)

    resid = f-model
    
    return resid



def log_prior(theta):
    '''
    Return prior on model parameters
    '''

    t0,P,rprs,aors,inc,ecc,omega,q1,q2,logf = theta
 
    if rprs > 0.5 or rprs < 0:
        return -np.inf

    b = aors*np.cos(np.radians(inc))*(1-ecc**2)/(1+ecc*np.sin(np.radians(omega)))
    
    if b > (1 + rprs) or b < 0:
        return -np.inf

    if aors < 0:
        return -np.inf

    if ecc > 1 or ecc < 0:
        return -np.inf

    if q1<0 or q1>1 or q2<0 or q2>1:
        return -np.inf

    if P < 0:
        return -np.inf

    if inc > 90.0 or inc < 0:
        return -np.inf
    
    return 0.0



def log_likelihood_all(theta, t, f, ferr, param_info):
    
    model = compute_trans(theta,t,param_info)

    sigma2 = ferr ** 2 + np.exp(2 * logf) # This was a fractional error in DFM example!

    loglike = -0.5 * np.sum((f - model) ** 2 / sigma2 + np.log(sigma2))
    logp = log_prior(theta)
    
    return loglike + logp


def log_likelihood_standard(theta, t, f, ferr, param_info):
    
    t0,rprs,aors,inc,q1,q2,logf = theta

    P = param_info['P']
    ecc = param_info['ecc']
    omega = param_info['omega']
    
    theta = t0,P,rprs,aors,inc,ecc,omega,q1,q2,logf
    
    model = compute_trans(theta,t,param_info)
    
    sigma2 = ferr ** 2 + np.exp(2 * logf) # This was a fractional error in DFM example!

    #probs = 1./np.sqrt(2*np.pi*sigma2)*np.exp(-0.5*((f-model)**2/sigma2))
    
    loglike =  -0.5 * np.sum((f - model) ** 2 / sigma2 + np.log(2*np.pi*sigma2))

    logp = log_prior(theta)

    return loglike + logp


def log_likelihood_nolimb(theta, t, f, ferr, param_info):
    
    t0,rprs,aors,inc,logf = theta
    
    P = param_info['P']
    ecc = param_info['ecc']
    omega = param_info['omega']
    q1 = param_info['q1']
    q2 = param_info['q2']

    theta = t0,P,rprs,aors,inc,ecc,omega,q1,q2,logf
    
    model = compute_trans(theta,t,param_info)
    sigma2 = ferr ** 2 + np.exp(2 * logf) # This was a fractional error in DFM example!
    loglike = -0.5 * np.sum((f - model) ** 2 / sigma2 + np.log(sigma2))
    logp = log_prior(theta)

    return loglike + prior




def fit_single(t,f,ferr,nwalkers=250,burnsteps=5000,mcmcsteps=5000,clobber=False,
               param_info={'t0':None,'P':None,'rprs':None,'aors':None,'inc':None,
                           'ecc':None,'omega':None,'q1':None,'q2':None,'logf':-10.,
                           'FIXED':[]}):

    """
    fit_single:
    -----------
    Fit a single transit signal with specified mcmc parameters return 
    chains and log likelihood.

    """

    if len(param_info['FIXED']) > 0:
        if 'q1' in param_info['FIXED']:
            loglike = log_likelihood_nolimb
            initial = [param_info['t0'],param_info['rprs'],param_info['aors'],
                       param_info['inc'],param_info['logf']]
        else:
            loglike = log_likelihood_standard
            initial = [param_info['t0'],param_info['rprs'],param_info['aors'],
                       param_info['inc'],param_info['q1'],param_info['q2'],param_info['logf']]
    else:
        loglike = log_likelihood_all
        initial = [param_info['t0'],param_info['P'],param_info['rprs'],param_info['aors'],
                   param_info['inc'],param_info['ecc'],param_info['omega'],
                   param_info['q1'],param_info['q2'],param_info['logf']]


    '''
    This used to get inital starting parameters for MCMC
    nll = lambda *args: -loglike(*args)

    initial = np.array(initial) 
    soln = minimize(nll, initial, args=(t, f, ferr,param_info),method='L-BFGS-B')
    m_ml, b_ml, log_f_ml = soln.x
    
    print("Maximum likelihood estimates:")
    print("m = {0:.3f}".format(m_ml))
    print("b = {0:.3f}".format(b_ml))
    print("f = {0:.3f}".format(np.exp(log_f_ml)))

    plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
    plt.plot(x0, m_true * x0 + b_true, "k", alpha=0.3, lw=3, label="truth")
    plt.plot(x0, np.dot(np.vander(x0, 2), w), "--k", label="LS")
    plt.plot(x0, np.dot(np.vander(x0, 2), [m_ml, b_ml]), ":k", label="ML")
    plt.legend(fontsize=14)
    plt.xlim(0, 10)
    plt.xlabel("x")
    plt.ylabel("y");
    '''


    nw = nwalkers
    bs = burnsteps
    mcs = mcmcsteps
    ndim = len(initial)
    
    directory = paths['output']+date+'/MCMC/'
    if not os.path.exists(directory):
        os.makedirs(directory)

# Do not redo MCMC unless clobber flag is set
    done = os.path.exists(directory+targname+'*steps_rprs.txt')
    if done == True and clobber == False:
        print("MCMC run already completed")
        return False,False

    print("")
    print("----------------------------------------------------------------------")
    os.system('date')
    print("Starting MCMC fitting routine for "+targname)


# Set up MCMC sampler
    print("... initializing emcee sampler")
    tstart = time.time()
    sampler = emcee.EnsembleSampler(nw, ndim, loglike, args=(t, f, ferr,param_info))

    twomin = 2./1440.
    onesec = 1./86400.

#    q1 = (sdata[0,3]+sdata[0,4])**2
#    q2 = sdata[0,3]/(2.0*(sdata[0,3]+sdata[0,4]))

# Initial chain values
    print("... deriving starting values for chains")
    p0_1 = np.random.uniform(initial[0]-twomin,initial[0]+twomin,nw)     # t0
    p0_2 = np.random.uniform(0.9*initial[1],1.1*initial[1], nw)          # Rp/R* 
    p0_3 = np.random.uniform(0.9*initial[2],1.1*initial[2], nw)          # a/R*
    p0_4 = np.random.uniform(89.5, 89.9, nw)                                 # inc
    p0_5 = np.random.normal(initial[4],0.01,nw)                          # q1
    p0_6 = np.random.normal(initial[5],0.01,nw)                          # q2
    p0_7 = np.random.uniform(-11,-10,nw)                                 # logf
    p0 = np.array([p0_1,p0_2,p0_3,p0_4,p0_5,p0_6,p0_7]).T
    variables =["t0","Rp/R*","a/R*","inc","q1","q2","logf"]

    # Run burn-in
    print("")
    print("Running burn-in with "+str(bs)+" steps and "+str(nw)+" walkers")
    pos, prob, state = sampler.run_mcmc(p0, bs, progress=True)


    # Plot chains
    samples = sampler.get_chain()
    labels = variables
    fig, axes = plt.subplots(len(labels), figsize=(10, 7), sharex=True)
    for i in range(len(labels)):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
    axes[-1].set_xlabel("step number");
    plt.savefig(directory+'burnin_chains.png',dpi=300)
    
    # Calculate G-R scale factor for each variable
    Rs = GR_test(sampler.chain,variables=variables)

    tau = sampler.get_autocorr_time(discard=500)    
    for var in np.arange(ndim):
        acout = "Autocorrelation time for "+variables[var]+" = {0:0.3f}"
        print(acout.format(tau[var]))

    print("")
    afout = "Mean acceptance fraction: {0:0.3f}"
    print afout.format(np.mean(sampler.acceptance_fraction))


    # Save burn in stats
    burn = np.append(Rs,tau)
    burn = np.append(burn,np.mean(sampler.acceptance_fraction))
    np.savetxt(directory+targname+'_burnstats.txt',burn)

    mcs = np.int(np.round(np.max(tau)*50,decimals=-2))
    
# Reset sampler and run MCMC for reals
    print "... resetting sampler and running MCMC with "+str(mcs)+" steps"
    sampler.reset()
    posf, probf, statef = sampler.run_mcmc(pos, mcs, progress=True)

 # Calculate G-R scale factor for each variable
    Rs = GR_test(sampler.chain,variables=variables)

# Autocorrelation times
    tau = sampler.get_autocorr_time(discard=500)    
    for var in np.arange(ndim):
        acout = "Autocorrelation time for "+variables[var]+" = {0:0.3f}"
        print acout.format(tau[var])

    afout = "Final mean acceptance fraction: {0:0.3f}"
    print afout.format(np.mean(sampler.acceptance_fraction))

    stats = np.append(Rs,tau)
    stats = np.append(stats,np.mean(sampler.acceptance_fraction))
    np.savetxt(directory+targname+'_fit_finalstats.txt',stats)

    # Plot chains
    samples = sampler.get_chain()
    labels = variables
    fig, axes = plt.subplots(len(labels), figsize=(10, 7), sharex=True)
    for i in range(len(labels)):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
    axes[-1].set_xlabel("step number");
    plt.savefig(directory+'final_chains.png',dpi=300)

    
    print "Writing MCMC chains to disk"
    tdist  = sampler.flatchain[:,0]
    rdist  = sampler.flatchain[:,1]
    adist  = sampler.flatchain[:,2]
    idist  = sampler.flatchain[:,3]
    q1dist  = sampler.flatchain[:,4]
    q2dist  = sampler.flatchain[:,5]
    logfdist= sampler.flatchain[:,6]

    np.savetxt(directory+targname+'_'+str(np.int(nw))+'walkers_'+str(np.int(mcs))+'steps_t0.txt',tdist)
    np.savetxt(directory+targname+'_'+str(np.int(nw))+'walkers_'+str(np.int(mcs))+'steps_rprs.txt',rdist)
    np.savetxt(directory+targname+'_'+str(np.int(nw))+'walkers_'+str(np.int(mcs))+'steps_aors.txt',adist)
    np.savetxt(directory+targname+'_'+str(np.int(nw))+'walkers_'+str(np.int(mcs))+'steps_inc.txt',idist)
    np.savetxt(directory+targname+'_'+str(np.int(nw))+'walkers_'+str(np.int(mcs))+'steps_q1.txt',q1dist)
    np.savetxt(directory+targname+'_'+str(np.int(nw))+'walkers_'+str(np.int(mcs))+'steps_q2.txt',q2dist)
    np.savetxt(directory+targname+'_'+str(np.int(nw))+'walkers_'+str(np.int(mcs))+'steps_logf.txt',logfdist)
        
    lp = sampler.lnprobability.flatten()
    np.savetxt(directory+targname+'_'+str(np.int(nw))+'walkers_'+str(np.int(mcs))+'steps_logprob.txt',lp)

    chains = sampler.flatchain

    return chains,lp

def done_in(tmaster):
    import time
    import numpy as np

    t = time.time()
    hour = (t - tmaster)/3600.
    if np.floor(hour) == 1:
        hunit = "hour"
    else:
        hunit = "hours"

    minute = (hour - np.floor(hour))*60.
    if np.floor(minute) == 1:
        munit = "minute"
    else:
        munit = "minutes"

    sec = (minute - np.floor(minute))*60.


    if np.floor(hour) == 0 and np.floor(minute) == 0:
        tout = "done in {0:.2f} seconds"
        out = tout.format(sec)
#        print tout.format(sec)
    elif np.floor(hour) == 0:
        tout = "done in {0:.0f} "+munit+" {1:.2f} seconds"
        out = tout.format(np.floor(minute),sec)
#        print tout.format(np.floor(minute),sec)
    else:
        tout = "done in {0:.0f} "+hunit+" {1:.0f} "+munit+" {2:.2f} seconds"
        out = tout.format(np.floor(hour),np.floor(minute),sec)
#        print tout.format(np.floor(hour),np.floor(minute),sec)

    print " "

    return out



def bin_lc(x,y,nbins=100):

    """
    ----------------------------------------------------------------------    
    bin_lc:
    -------
    Utility to bin data and return standard deviation in each bin
    
    For visual aid in plots, mostly

    example:
    --------
    tbin,fbin,errbin = bin_lc(ttm,flux,nbins=200)

    """

    n, I = np.histogram(x, bins=nbins)
    sy, I = np.histogram(x, bins=nbins, weights=y)
    sy2, I = np.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)

    binvals = (I[1:] + I[:-1])/2
    yvals = mean
    yerr = std/np.sqrt(len(std))
    
    return binvals,yvals,yerr



def plot_model(modelparams,short=False,tag='',markersize=5,smallmark=2,
               nbins=100,errorbars=False,pdf=False):

    """
    ----------------------------------------------------------------------
    plot_model:
    -----------
    Plot transit model given model params.

    ----------------------------------------------------------------------
    """


    plt.rcdefaults()

# Check for output directory   
    directory = path+'MCMC/'
    if not os.path.exists(directory):
        os.makedirs(directory)

# Fold data on input parameters
    tplt = foldtime(t,period=modelparams[4],t0=pdata[0,3]+modelparams[3])
    fplt = flux
    eplt = e_flux

# Bin data (for visual aid)
    tfit, yvals, yerror = bin_lc(tplt,flux,nbins=nbins)

    plt.figure(2,figsize=(11,8.5),dpi=300)
    plt.subplot(2,1,1)
    plt.plot(tplt,(fplt-1)*1e6,'.',color='gray',markersize=smallmark)
    if errorbars:
        plt.errorbar(tfit,(yvals-1)*1e6,yerr=yerror*1e6,fmt='o',color='blue',
                     markersize=markersize)
    else:
        plt.plot(tfit,(yvals-1)*1e6,'bo',markersize=markersize)

# Center transit in plot
    wid1 = np.max(tplt)
    wid2 = abs(np.min(tplt))
    if wid1 < wid2:
        wid = wid1
    else:
        wid = wid2
    plt.xlim(np.array([-1,1])*wid)

    ldc = modelparams[5]

# Get model, raw (unsmoothed) model, and tfull

    tmodel,model,rawmodel = compute_trans(modelparams[0],modelparams[1],modelparams[2],\
                                              0.0,modelparams[4],\
                                              ldc,unsmooth=True)

#    tfull = tottofull(modelparams[1],modelparams[0],modelparams[2],modelparams[4])
 
    sig = rb.std((fplt-1)*1e6)
    med = np.median((fplt-1)*1e6)
#    min = max(np.min((model-1)*1e6),-4*sig)
# For deep transits
    min = np.min((model-1)*1e6)
#    yrange = np.array([min-3*sig,med+15*sig])
    yrange = np.array([min-4*sig,med+12*sig])
    plt.ylim(yrange) 
    locs,labels = plt.yticks()
    plt.yticks(locs, map(lambda x: "%.0f" % x, locs))
    plt.ylabel("ppm")

    plt.plot(tmodel,(rawmodel-1)*1e6,'g')
    plt.plot(tmodel,(model-1)*1e6,'r')
    plt.title(name+" Transit Fit")

    res = residuals(modelparams)

    dof = len(res) - ndim - 1.
    chisq = np.sum((res/eplt)**2)/dof

    rhoi = sdata[0,0]/(4./3.*np.pi*sdata[0,1]**3)

    aors = np.sqrt( ((1+modelparams[0])**2 - modelparams[2]**2)/
                    ((np.sin(np.pi*modelparams[1]/modelparams[4]))**2) + modelparams[2]**2)    

#    aors    = 2.0 * np.sqrt(modelparams[0]) * modelparams[4] / \
#        (np.pi*np.sqrt(ttot**2 - tfull**2))
    
    rhostar =  3.0*np.pi/( c.G * (modelparams[4]*24.*3600.)**2 ) * aors**3
    

    plt.annotate(r'$P$ = %.7f d' % modelparams[4], [0.5,0.87],horizontalalignment='center',
                 xycoords='figure fraction',fontsize='large')

    plt.annotate(r'$\chi^2_r$ = %.5f' % chisq, [0.87,0.85],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
    plt.annotate(r'$\rho_*$ (e=0) = %.3f' % rhostar, [0.87,0.81],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
    plt.annotate(r'$\rho_*$ (orig) = %.3f' % rhoi, [0.87,0.77],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
    val =  modelparams[1]*24.
    plt.annotate(r'$\tau_{\rm tot}$ = %.4f h' % val, [0.87,0.73],horizontalalignment='right',
                  xycoords='figure fraction',fontsize='large')
    val = (modelparams[4]-pdata[0,4])*24.*3600.
    plt.annotate(r'$\Delta P$ = %.3f s' % val, [0.15,0.85],
                 xycoords='figure fraction',fontsize='large')
    plt.annotate(r'$R_p/R_*$ = %.5f' % modelparams[0], [0.15,0.81],
                  xycoords='figure fraction',fontsize='large')   
    plt.annotate('b = %.2f' % modelparams[2], [0.15,0.77],
                  xycoords='figure fraction',fontsize='large')
    t0out = (modelparams[3])*24.*60.*60.
    plt.annotate(r'$\Delta t_0$ = %.3f s' % t0out, [0.15,0.73],
                  xycoords='figure fraction',fontsize='large')

   # Limb darkening parameters
    if limb == 'nlin':
        a1 = ldc[0]
        a2 = ldc[1]
        a3 = ldc[2]
        a4 = ldc[3]
        # need annotation for non-linear LD fits
    else:        
        q1in = ldc[0]
        q2in = ldc[1]
        u1, u2 = qtou(q1in,q2in,limb=limb)

        u1out = u1
        plt.annotate(r'$u_1$ = %.2f' % u1out, [0.15,0.59],
                     xycoords='figure fraction',fontsize='large')
        u2out = u2
        plt.annotate(r'$u_2$ = %.2f' % u2out, [0.15,0.55],
                     xycoords='figure fraction',fontsize='large')

    plt.subplot(2,1,2)
    s = np.argsort(tplt)
    plt.plot(tplt[s],res*1e6,'.',markersize=smallmark,color='gray')
    tres, yres, yreserr = bin_lc(tplt[s],res,nbins=nbins)
    if errorbars:
        plt.errorbar(tres,yres*1e6,yerr=yreserr*1e6,fmt='o',color='blue',markersize=markersize)
    else:
        plt.plot(tres,yres*1e6,'bo',markersize=markersize)

    plt.xlim(np.array([-1,1])*wid)
    sig = rb.std(res*1e6)
    med = np.median(res*1e6)
    plt.ylim(np.array([-5*sig,5*sig]))
    plt.axhline(y=0,color='r')
    locs,labels = plt.yticks()
    plt.yticks(locs, map(lambda x: "%.0f" % x, locs))
    plt.title("Residuals")
    plt.xlabel("Time from Mid Transit (days)")
    plt.ylabel("ppm")

    if pdf:
        ftype = '.pdf'
    else:
        ftype = '.png'

    plt.savefig(directory+targname+stag+ctag+rtag+lptag+ltag+'fit'+tag+ftype, dpi=300)
    plt.clf()

    gc.collect()
    return



def plot_final(modelparams,dispparams=False,short=False,tag='',ms=10,sm=8,
               errorbars=False,pdf=False):

    """
    ----------------------------------------------------------------------
    plot_final:
    -----------
    Plot transit model given model params.

    ----------------------------------------------------------------------
    """

    import gc
    from matplotlib import gridspec
    import matplotlib as mpl
    plt.rcdefaults()

    titletxt = 'Long Cadence'
    alpha = 0.5
    ysigmin = 1.5
    ysigmax = 3
    ms = 10
    sm = 8
    intfac = 1.0

    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['xtick.major.size'] = 5
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['ytick.major.size'] = 5
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams['xtick.labelsize'] = 18
    mpl.rcParams['ytick.labelsize'] = 18
    fs = 20
    lw = 1.5

# Check for output directory   
    directory = paths['output']+date+'/MCMC/'
    if not os.path.exists(directory):
        os.makedirs(directory)

# Bin data (for visual aid)
    #tfit, yvals, yerror = bin_lc(tplt,flux,nbins=nbins)

    plt.figure(2,figsize=(11,8.5),dpi=300)
    gs = gridspec.GridSpec(3, 1,wspace=0)
    ax1 = plt.subplot(gs[0:2, 0])    

    tplt = (t-modelparams[3])*24
    ax1.plot(tplt,(f-1)*1e3,'o',color='black',markersize=sm,alpha=0.75)
    #if errorbars:
    #    ax1.errorbar(t,(f-1)*1e3,yerr=yerror*1e6,fmt='o',color='black',
    #                 markersize=ms)
    #else:
    #    ax1.plot(tfit,(yvals-1)*1e6,'o',color='black',markersize=ms)

# Center transit in plot
    #ax1.set_xlim(np.array([-1,1])*wid)


# Center transit in plot
    wid1 = np.max(tplt)
    wid2 = abs(np.min(tplt))
    if wid1 < wid2:
        wid = wid1
    else:
        wid = wid2

    nbins = 2.0*wid*86400/(integration*intfac)

# Get model, raw (unsmoothed) model, and tfull
    theta = [modelparams[3],param_info['P'],modelparams[0],modelparams[1],modelparams[2],
             param_info['ecc'],param_info['omega'],modelparams[4],modelparams[5],modelparams[6]]

    tmodel = np.linspace(np.min(t),np.max(t),1000)
    model = compute_trans(theta, tmodel, param_info)
    tmodelplt = (tmodel-modelparams[3])*24
    
    ax1.plot(tmodelplt,(model-1)*1e3,'r',linewidth=lw)

    sig = mad_std((f-1)*1e3)
    med = np.median((f-1)*1e3)
#    min = max(np.min((model-1)*1e6),-4*sig)
# For deep transits
    mmin = np.min((model-1)*1e3)
#    yrange = np.array([min-3*sig,med+15*sig])
    yrange = np.array([mmin-ysigmin*sig,med+ysigmax*sig])
    ax1.set_ylim(yrange)
    ax1.set_xlim((np.min(tplt),np.max(tplt)))
    locs = ax1.yaxis.get_majorticklocs()
    ax1.yaxis.set_ticks(locs, map(lambda x: "%.0f" % x, locs))
    ax1.set_ylabel("ppt",fontsize=fs)
    ax1.set_xticklabels(())
    ax1.xaxis.set_tick_params(length=10, width=lw, labelsize=fs)
    ax1.yaxis.set_tick_params(length=10, width=lw, labelsize=fs)

    
    res = residuals(theta, t, f, param_info)

    ndim = len(modelparams)
    dof = len(res) - ndim - 1.
    chisq = np.sum((res/ferr)**2)/dof

    
#    rprs     = dispparams[0] 
#    duration = dispparams[1]
#    impact   = dispparams[2]
#    period   = dispparams[4]
#    aors     = np.sqrt( ((1+rprs)**2 - impact**2)/
#                        ((np.sin(np.pi*duration/period))**2) + impact**2)    
    
#    rhostar =  3.0*np.pi/( c.G * (period*24.*3600.)**2 ) * aors**3
    
#    ax1.annotate(r'$P$ = %.5f d' % period, [0.55,0.86],horizontalalignment='center',
#                 xycoords='figure fraction',fontsize=fs-1)
#    ax1.annotate(r'$\chi^2_r$ = %.3f' % chisq, [0.9,0.83],horizontalalignment='right',
#                 xycoords='figure fraction',fontsize=fs-1)
#    ax1.annotate(r'$\rho_*$ (e=0) = %.2f g/cc' % rhostar, [0.9,0.78],horizontalalignment='right',
#                 xycoords='figure fraction',fontsize=fs-1)
#    ax1.annotate(r'$\rho_{*,0}$ = %.2f g/cc' % rhoi, [0.9,0.73],horizontalalignment='right',
#                 xycoords='figure fraction',fontsize=fs-1)
#    ax1.annotate(r'$R_p/R_*$ = %.3f' % dispparams[0], [0.2,0.83],
#                  xycoords='figure fraction',fontsize=fs-1)   
#    val =  duration*24.
#    ax1.annotate(r'$\tau_{\rm tot}$ = %.2f h' % val, [0.2,0.78],
#                  xycoords='figure fraction',fontsize=fs-1)
#    t0out = (dispparams[3]) + pdata[0,3] # + 2454833.0
#    ax1.annotate(r'$t_0$ = %.4f BKJD' % t0out, [0.2,0.73],
#                  xycoords='figure fraction',fontsize=fs-1)

    ax2 = plt.subplot(gs[2, 0])
    ax2.plot(tplt,res*1e3,'o',markersize=sm,color='black',alpha=0.75)

    #tres, yres, yreserr = bin_lc(tplt[s],res,nbins=nbins)
    #if errorbars:
    #    ax2.errorbar(tres,yres*1e6,yerr=yreserr*1e6,fmt='o',color='black',markersize=ms)
    #else:
    #    ax2.plot(tres,yres*1e6,'o',color='black',markersize=ms)
    ax2.set_xlim((np.min(tplt),np.max(tplt)))

    sig = mad_std(res*1e3)
    med = np.median(res*1e3)
    ax2.set_ylim(np.array([-7*sig,7*sig]))
    ax2.axhline(y=0,color='r',linewidth=lw)
    locs = ax2.yaxis.get_majorticklocs()
    ax2.yaxis.set_ticks(locs, map(lambda x: "%.0f" % x, locs))
#    ax2.set_title("Residuals")
    ax2.set_xlabel("Time from Mid-Transit (hours)",fontsize=fs)
    ax2.set_ylabel("ppt",fontsize=fs)
    ax2.xaxis.set_tick_params(length=10, width=lw, labelsize=fs)
    ax2.yaxis.set_tick_params(length=10, width=lw, labelsize=fs)

    ax1.set_title('Transit fit for '+targname,fontsize=fs+2)

    plt.subplots_adjust(hspace=0.15,left=0.15,right=0.95,top=0.92)

    
    if pdf:
        ftype='.pdf'
    else:
        ftype='.png'
    plt.savefig(directory+targname+'_finalfit'+ftype, dpi=300)

#    plt.rcdefaults()
    plt.subplots_adjust(hspace=0.2,left=0.125,right=0.9,top=0.9)
    plt.clf()


    gc.collect()
    return




def qtou(q1,q2,limb='quadratic'):
    '''
    Transform between Kipping q parameters and the Claret u parameters for limb darkening
    '''
    if limb == 'quadratic':
        u1 =  2*np.sqrt(q1)*q2
        u2 =  np.sqrt(q1)*(1-2*q2)
    if limb == 'square-root':
        u1 = np.sqrt(q1)*(1-2*q2)
        u2 = 2*np.sqrt(q1)*q2

    return u1, u2


def utoq(u1,u2,limb='quadratic'):
    '''
    Transform between the Claret u parameters and the Kipping q parameters for limb darkening
    '''
    if limb == 'quadratic':
        q1 = (u1+u2)**2
        q2 = u1/(2.0*(u1+u2))
    if limb == 'square-root':
        q1 = (u1+u2)**2
        q2 = u2/(2.0*(u1+u2))

    return q1, q2


def get_limb_qs(Mstar=0.5,Rstar=0.5,Tstar=3800.0,limb='quad',network=None):
    '''
    Given the mass, radius and effective temperature of a star, estimate limb
    darkening from the Claret grid

    In cgs!
    '''

    import constants as c

    Ms = Mstar*c.Msun
    Rs = Rstar*c.Rsun
    loggstar = np.log10( c.G * Ms / Rs**2. )
 
    
    if limb == 'nlin':
        a,b,c,d = get_limb_coeff(Tstar,loggstar,network=network,limb=limb,interp='linear')
        return a,b,c,d
    else:
        a,b = get_limb_coeff(Tstar,loggstar,network=network,limb=limb,interp='linear')
        q1,q2 = utoq(a,b,limb=limb)
        return q1, q2
 

def get_qdists(sdata,sz=100,errfac=1):
    from scipy.stats.kde import gaussian_kde
    global q1pdf_func
    global q2pdf_func

    Mstar = sdata[0,0]/c.Msun
    eMstar = sdata[1,0]/c.Msun
    Rstar = sdata[0,1]/c.Rsun
    eRstar = sdata[1,1]/c.Rsun
    Tstar = sdata[0,2]
    eTstar = sdata[1,2]
    
    print "... using error factor of "+str(errfac)
    Ms = np.random.normal(Mstar,eMstar*errfac,sz)
    Rs = np.random.normal(Rstar,eRstar*errfac,sz)
    Ts = np.random.normal(Tstar,eTstar*errfac,sz)

    q1v = []
    q2v = []

    print "... generating LD distribution of "+str(sz)+" values"
    for i in range(sz):
        q1,q2 = get_limb_qs(Mstar=max(Ms[i],0.001),Rstar=max(Rs[i],0.001),Tstar=max(Ts[i],100),limb=limb)
        q1v = np.append(q1v,q1)
        q2v = np.append(q2v,q2)
        
    q1v = q1v[~np.isnan(q1v)]
    q2v = q2v[~np.isnan(q2v)]

    vals  = np.linspace(0,1,10000)
    q1kde = gaussian_kde(q1v)
    q1pdf = q1kde(vals)
    q1pdf_func = sp.interpolate.interp1d(vals,q1pdf,kind='nearest')

    q2kde = gaussian_kde(q2v)
    q2pdf = q2kde(vals)
    q2pdf_func = sp.interpolate.interp1d(vals,q2pdf,kind='nearest')


    return q1v, q2v



def get_qvals(q1v,q2v,nsamp=100):
    from scipy.stats.kde import gaussian_kde
    global q1pdf_func
    global q2pdf_func

    vals  = np.linspace(0,1,1000)
    q1kde = gaussian_kde(q1v)
    q1pdf = q1kde(vals)
    q1pdf_func = sp.interpolate.interp1d(vals,q1pdf,kind='linear')
    q1c   = np.cumsum(q1pdf)/np.sum(q1pdf)
    q1func = sp.interpolate.interp1d(q1c,vals,kind='linear')
    q1samp = q1func(np.random.uniform(0,1,nsamp))

    q2kde = gaussian_kde(q2v)
    q2pdf = q2kde(vals)
    q2pdf_func = sp.interpolate.interp1d(vals,q2pdf,kind='linear')
    q2c   = np.cumsum(q2pdf)/np.sum(q2pdf)
    q2func = sp.interpolate.interp1d(q2c,vals,kind='linear')
    q2samp = q2func(np.random.uniform(0,1,nsamp))

    return q1samp, q2samp




def GR_test(chains,variables=False):

    """
    ----------------------------------------------------------------------
    GR_test:
    --------
    Compute the Gelman-Rubin scale factor for each variable given input
    flat chains
    ----------------------------------------------------------------------
    """

    nwalkers = np.float(np.shape(chains)[0])
    nsamples = np.float(np.shape(chains)[1])
    ndims    = np.shape(chains)[2]
    Rs = np.zeros(ndims)
    for var in np.arange(ndims):
        psi = chains[:,:,var]
        psichainmean = np.mean(psi,axis=1)
        psimean = np.mean(psi)

        B = nsamples/(nwalkers-1.0)*np.sum((psichainmean - psimean)**2)

        s2j = np.zeros(np.int(nwalkers))
        for j in range(np.int(nwalkers)):
            s2j[j] = 1.0/(nsamples-1.0)*np.sum((psi[j,:] - psichainmean[j])**2)

        W = np.mean(s2j)

        varbarplus = (nsamples-1.0)/nsamples * W + 1/nsamples * B

        R = np.sqrt(varbarplus/W)

        if len(variables) == ndims:
            out = "Gelman Rubin scale factor for "+variables[var]+" = {0:0.3f}"
            print(out.format(R))

        Rs[var] = R

    return Rs

    

def bestvals(targname,param_info,date,bindiv=50,thin=False,
             frac=0.003,nbins=100,rpmax=1,durmax=10,sigfac=5.0,pdf=False):

    """
    ----------------------------------------------------------------------
    bestvals:
    ---------
    Find the best values from the 1-d posterior pdfs return best values 
    and the posterior pdf for rp/rs
    ----------------------------------------------------------------------
    """
    from statsmodels.nonparametric.kde import KDEUnivariate as KDE_U
    from scipy.stats.kde import gaussian_kde

    plt.rcdefaults()

    directory = paths['output']+date+'/MCMC/'

    print 'Importing MCMC chains...'
    chainfile = glob.glob(directory+targname+'*t0.txt')    
    tdist = np.loadtxt(chainfile[0])
    tdist = (tdist - param_info['t0'])*86400
    chainfile = glob.glob(directory+targname+'*rprs.txt')
    rdist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*aors.txt')
    adist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*inc.txt')
    idist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*q1.txt')
    q1dist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*q2.txt')
    q2dist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*logf.txt')    
    logfdist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*logprob.txt')    
    lp = np.loadtxt(chainfile[0])

#  Get maximum likelihood values
    maxlike = np.max(lp)
    imax = np.argmax(lp) 
    t0val = np.float(tdist[imax])
    rprsval = np.float(rdist[imax])
    aorsval = np.float(adist[imax])
    incval  = np.float(idist[imax])
    q1val   = np.float(q1dist[imax])
    q2val   = np.float(q2dist[imax])
    logfval = np.float(logfdist[imax])
        
    if thin:
        tdist = tdist[0::thin]
        rdist = rdist[0::thin]
        adist = adist[0::thin]
        idist = idist[0::thin]
        q1dist = q1dist[0::thin]
        q2dist = q2dist[0::thin]
        logfdist = logfdist[0::thin]
        lp     = lp[0::thin]

    nsamp = len(tdist)


    print('')
    print('Best fit parameters for '+targname)
    plt.figure(4,figsize=(8.5,11),dpi=300)    
   

    # Rp/Rstar
    sigsamp = 5.0
    rprsmin = np.min(rdist)
    rprsmax = np.max(rdist)
    rprssig = mad_std(rdist)
    bw = rprssig/sigsamp
    rprsmed = np.median(rdist)
    minval = max(rprsmed - sigfac*rprssig,rprsmin)
    maxval = min(rprsmed + sigfac*rprssig,rprsmax)
    if rprsval < 1.0e-6:
        print "It looks like there is no planet!!!"
        rprsmode = (maxval-minval)/2.0
        rprshi = maxval
        rprslo = minval
        nb = 100
        rprss = np.array([minval,maxval])
        rprspdf = np.array([1,1])
    else:
        rprss = np.linspace(rprsmin-2*np.std(rdist),rprsmax+2*np.std(rdist),100)
        rprs_kde = gaussian_kde(rdist)
        rprspdf = rprs_kde(rprss)
        rprsmode = rprss[np.argmax(rprspdf)]
        ps = np.percentile(rdist, [15.87, 50, 84.13])
        rprshi = ps[2]
        rprslo = ps[0]
        rprsmed = ps[1]
        nb = np.int(np.round(np.ceil((rprsmax-rprsmin) / (np.abs(maxval-minval)/bindiv))))

    rprsout = 'Rp/R*: max = {0:.5f}, med = {1:.5f}, mode = {2:.5f}, 1 sig int = {3:.5f}'
    print(rprsout.format(rprsval, rprsmed, rprsmode, rprshi-rprslo))
    plt.subplot(4,1,1)
    plt.hist(rdist,bins=nb,normed=True)
#    plt.plot(rprss,rprspdf,color='c')
#    rprsmin = np.float(np.sort(rprsdist)[np.round(frac*nsamp)])
#    rprsmax = min(np.float(np.sort(rprsdist)[np.round((1-frac)*nsamp)]),rpmax)
    plt.xlim([minval,maxval])
#    plt.axvline(x=rprsval,color='r',linestyle='--')
#    plt.axvline(x=rprslo,color='c',linestyle='--')
#    plt.axvline(x=rprsmed,color='c',linestyle='--')
#    plt.axvline(x=rprshi,color='c',linestyle='--')
    plt.xlabel(r'$R_p/R_{*}$')
    plt.ylabel(r'$dP/d(R_p/R_{*})$')
    plt.title('Parameter Distributions for '+targname)
    

    # a/R*                          
    aorsmin = np.min(adist)
    aorsmax = np.max(adist)
    aorssig = mad_std(adist)
    bw = aorssig/sigsamp
    aorsmed = np.median(adist)
    minval = max(aorsmed - sigfac*aorssig,aorsmin)
    maxval = min(aorsmed + sigfac*aorssig,aorsmax)
    aorss = np.linspace(aorsmin-2*np.std(adist),aorsmax+2*np.std(adist),100)
    aors_kde = gaussian_kde(adist)
    aorspdf = aors_kde(aorss)
    aorsmode = aorss[np.argmax(aorspdf)]
    ps = np.percentile(adist, [15.87, 50, 84.13])
    aorshi = ps[2]
    aorslo = ps[0]
    aorsmed = ps[1]
    nb = np.int(np.round(np.ceil((aorsmax-aorsmin) / (np.abs(maxval-minval)/bindiv))))

    out = 'a/R*: max = {0:.4f}, med = {1:.4f}, mode = {2:.4f}, 1 sig int = {3:.4f}'
    print(out.format(aorsval, aorsmed, aorsmode, aorshi-aorslo))

    plt.subplot(4,1,2)    
    plt.hist(adist,bins=nb,normed=True)
#    plt.plot(ds,dpdf,color='c')
#    plt.axvline(x=dval,color='r',linestyle='--')
#    plt.axvline(x=dlo,color='c',linestyle='--')
#    plt.axvline(x=dmed,color='c',linestyle='--')
#    plt.axvline(x=dhi,color='c',linestyle='--')
#    dmin = np.float(np.sort(ddist)[np.round(frac*nsamp)])
#    dmax = min(np.float(np.sort(ddist)[np.round((1-frac)*nsamp)]),durmax)
    plt.xlim([minval,maxval])
    plt.xlabel(r'a/R*')
    plt.ylabel(r'$dP/d(a/R_{*})$')


   # Inclination
    incmin = np.min(idist)
    incmax = np.max(idist)
    incsig = mad_std(idist)
    bw = incsig/sigsamp
    incmed = np.median(idist)
    minval = max(incmed - sigfac*incsig,incmin)
    maxval = min(incmed + sigfac*incsig,incmax)
    incs = np.linspace(incmin-2*np.std(idist),incmax+2*np.std(idist),100)
    inc_kde = gaussian_kde(idist)
    incpdf = inc_kde(incs)
    incmode = incs[np.argmax(incpdf)]
    ps = np.percentile(idist, [15.87, 50, 84.13])
    inchi = ps[2]
    inclo = ps[0]
    incmed = ps[1]
    nb = np.int(np.round(np.ceil((incmax-incmin) / (np.abs(maxval-minval)/bindiv))))

    out = 'inclination: max = {0:.4f}, med = {1:.4f}, mode = {2:.4f}, 1 sig int = {3:.4f}'
    print(out.format(incval, incmed, incmode, inchi-inclo))


    plt.subplot(4,1,3)
    plt.hist(idist,bins=nb,normed=True)
#    plt.plot(bss,incpdf,color='c')
#    bmin = np.float(np.sort(bdist)[np.round(frac*nsamp)])
#    bmax = np.float(np.sort(bdist)[np.round((1-frac)*nsamp)])
    plt.xlim([minval,maxval])
#    plt.axvline(x=bval,color='r',linestyle='--')
#    plt.axvline(x=blo,color='c',linestyle='--')
#    plt.axvline(x=bmed,color='c',linestyle='--')
#    plt.axvline(x=bhi,color='c',linestyle='--')
    plt.xlabel('Inclination (degrees)')
    plt.ylabel(r'$dP/di$')


   # t0
    t0min = np.min(tdist)
    t0max = np.max(tdist)
    t0sig = mad_std(tdist)
    bw = t0sig/sigsamp
    t0med = np.median(tdist)
    minval = max(t0med - sigfac*t0sig,t0min)
    maxval = min(t0med + sigfac*t0sig,t0max)
    t0s = np.linspace(t0min-2*np.std(tdist),t0max+2*np.std(tdist),100)
    t0_kde = gaussian_kde(tdist)
    t0pdf = t0_kde(t0s)
    t0mode = t0s[np.argmax(t0pdf)]
    ps = np.percentile(tdist, [15.87, 50, 84.13])
    t0hi = ps[2]
    t0lo = ps[0]
    t0med = ps[1]
    nb = np.int(np.round(np.ceil((t0max-t0min) / (np.abs(maxval-minval)/bindiv))))

    out = 't0: max = {0:.4f}, med = {1:.4f}, mode = {2:.4f}, 1 sig int = {3:.4f}'
    print(out.format(t0val, t0med, t0mode, t0hi-t0lo))


    plt.subplot(4,1,4)
    plt.hist(tdist,bins=nb,normed=True)
#    plt.plot(t0s,t0pdf,color='c')
#    bmin = np.float(np.sort(bdist)[np.round(frac*nsamp)])
#    bmax = np.float(np.sort(bdist)[np.round((1-frac)*nsamp)])
    plt.xlim([minval,maxval])
#    plt.axvline(x=bval,color='r',linestyle='--')
#    plt.axvline(x=blo,color='c',linestyle='--')
#    plt.axvline(x=bmed,color='c',linestyle='--')
#    plt.axvline(x=bhi,color='c',linestyle='--')
    plt.xlabel(r'$\Delta$ Mid-Transit Time (s)')
    plt.ylabel(r'$dP/dt$')
    plt.annotate(r'$t_0$ = %.6f d' % param_info['t0'], xy=(0.97,0.8),
                 ha="right",xycoords='axes fraction',fontsize='large')

    plt.subplots_adjust(hspace=0.4)

    plt.savefig(directory+targname+'_fit_params1.png', dpi=300)
    plt.clf()



# Second set of parameters

    plt.figure(5,figsize=(8.5,11),dpi=300)    

 # q1
    q1min = np.min(q1dist)
    q1max = np.max(q1dist)
    q1sig = mad_std(q1dist)
    bw = q1sig/sigsamp
    q1med = np.median(q1dist)
    minval = max(q1med - sigfac*q1sig,q1min)
    maxval = min(q1med + sigfac*q1sig,q1max)
    q1s = np.linspace(q1min-2*np.std(q1dist),q1max+2*np.std(q1dist),100)
    q1_kde = gaussian_kde(q1dist)
    q1pdf = q1_kde(q1s)
    q1mode = q1s[np.argmax(q1pdf)]
    ps = np.percentile(q1dist, [15.87, 50, 84.13])
    q1hi = ps[2]
    q1lo = ps[0]
    q1med = ps[1]
    nb = np.int(np.round(np.ceil((q1max-q1min) / (np.abs(maxval-minval)/bindiv))))

    out = 'Kipping q1: max = {0:.4f}, med = {1:.4f}, mode = {2:.4f}, 1 sig int = {3:.4f}'
    print(out.format(q1val, q1med, q1mode, q1hi-q1lo))

    plt.subplot(2,1,1)
    
    plt.hist(q1dist,bins=nb,normed=True)
 #    plt.plot(q1s,q1pdf,color='c')
    plt.xlim([0,1])
#    plt.axvline(x=q1val,color='r',linestyle='--')
#    plt.axvline(x=q1lo,color='c',linestyle='--')
#    plt.axvline(x=q1med,color='c',linestyle='--')
#    plt.axvline(x=q1hi,color='c',linestyle='--')
    plt.xlabel(r'$q_1$')
    plt.ylabel(r'$dP/q_1$')

# q2
    q2min = np.min(q2dist)
    q2max = np.max(q2dist)
    q2sig = mad_std(q2dist)
    bw = q2sig/sigsamp
    q2med = np.median(q2dist)
    minval = max(q2med - sigfac*q2sig,q2min)
    maxval = min(q2med + sigfac*q2sig,q2max)
    q2s = np.linspace(q2min-2*np.std(q2dist),q2max+2*np.std(q2dist),100)
    q2_kde = gaussian_kde(q2dist)
    q2pdf = q2_kde(q2s)
    q2mode = q2s[np.argmax(q2pdf)]
    ps = np.percentile(q2dist, [15.87, 50, 84.13])
    q2hi = ps[2]
    q2lo = ps[0]
    q2med = ps[1]
    nb = np.int(np.round(np.ceil((q2max-q2min) / (np.abs(maxval-minval)/bindiv))))

    out = 'Kipping q2: max = {0:.4f}, med = {1:.4f}, mode = {2:.4f}, 1 sig int = {3:.4f}'
    print(out.format(q2val, q2med, q2mode, q2hi-q2lo))

    plt.subplot(2,1,2)
    
    plt.hist(q2dist,bins=nb,normed=True)
 #    plt.plot(q2s,q2pdf,color='c')
    plt.xlim([0,1])
#    plt.axvline(x=q2val,color='r',linestyle='--')
#    plt.axvline(x=q2lo,color='c',linestyle='--')
#    plt.axvline(x=q2med,color='c',linestyle='--')
#    plt.axvline(x=q2hi,color='c',linestyle='--')
    plt.xlabel(r'$q_2$')
    plt.ylabel(r'$dP/q_2$')
    
    plt.subplots_adjust(hspace=0.4)

    plt.savefig(directory+targname+'_fit_params2.png', dpi=300)
    plt.clf()




'''
# calculate limb darkening parameters
    if limb == 'quad' or limb == 'sqrt':
        u1, u2 = qtou(q1val,q2val,limb=limb)
        ldc = [u1,u2]
    elif limb == 'nlin':
        ldc = [q1val,q2val,q3val,q4val]
    else:
        pdb.set_trace()

    plot_limb_curves(ldc=ldc,limbmodel=param_info['limb_model'],write=True)
'''
 
# vals
# 0 = rprsval
# 1 = aorsval
# 2 = inclination 
# 3 = offset time in seconds
# 4 = Kipping q1 val
# 5 = Kipping q2 val
# 6 = Log f (noise factor)

    vals = [rprsval,aorsval,incval,t0val/86400. + param_info['t0'],q1val,q2val,logfval]
    meds = [rprsmed,aorsmed,incmed,t0med/86400. + param_info['t0'],q1med,q2med,np.median(logfdist)]
    modes = [rprsmode,aorsmode,incmode,t0mode/86400. + param_info['t0'],q1mode,q2mode,np.median(logfdist)]
    onesig = [rprshi-rprslo,aorshi-aorslo,(t0hi-t0lo)inchi-inclo,(thi-tlo)/(24.*3600.),(phi-plo)/(24.*3600.),q1hi-q1lo,q2hi-q2lo]

    bestvals = [[vals],[meds],[modes],[onesig]]

    fmodelvals = [vals[0],vals[1],vals[2],vals[3]-pdata[0,3],vals[4],[vals[5],vals[6]]]

    plot_model(fmodelvals,tag='_MCMCfit',pdf=pdf)

    plot_final(fmodelvals,tag='_MCMCfit',pdf=pdf)


    outstr = name+ ' %.5f  %.5f  %.5f  %.5f  %.4f  %.4f  %.4f  %.4f  %.2f  %.2f  %.2f  %.2f  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.0f' % (vals[0],meds[0],modes[0],onesig[0],vals[1],meds[1],modes[1],onesig[1],vals[2],meds[2],modes[2],onesig[2],vals[3],meds[3],modes[3],onesig[3],vals[4],meds[4],modes[4],onesig[4],vals[5],meds[5],modes[5],onesig[5],vals[6],meds[6],modes[6],onesig[6],sampfac)

    f = open(directory+targname+'_fitparams.txt','w')
    f.write(outstr+'\n')
    f.closed

    # collect garbage
    gc.collect()
    return bestvals




def distparams(dist):

    vals = np.linspace(np.min(dist)*0.5,np.max(dist)*1.5,1000)
    kde = gaussian_kde(dist)
    pdf = b_kde(vals)
    dist_c = np.cumsum(pdf)/np.sum(pdf)
    func = sp.interpolate.interp1d(dist_c,vals,kind='linear')
    lo = np.float(func(math.erfc(1./np.sqrt(2))))
    hi = np.float(func(math.erf(1./np.sqrt(2))))
    med = np.float(func(0.5))
    mode = vals[np.argmax(pdf)]

    disthi = np.linspace(.684,.999,100)
    distlo = distthi-0.6827
    disthis = func(distthi)
    distlos = func(disttlo)
    
    interval = np.min(rhohis-rholos)

    return med,mode,interval,lo,hi



def triangle_plot(modelparams,targname=None,thin=False,sigfac=3.0,bindiv=75,sigsamp=5.0,pdf=False):
    import numpy as np
    import scipy as sp
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import robust as rb
    import sys,math,pdb,time,gc
    from scipy.stats.kde import gaussian_kde
    import matplotlib as mpl
    plt.rcdefaults()

    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['xtick.major.size'] = 5
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['ytick.major.size'] = 5
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    fs = 14
    lw = 1.5

    tmaster = time.time()
    print("Reading in MCMC chains")
    '''
    if chains is not False:
        rprsdist = chains[:,0]
        ddist = chains[:,1]*24.
        bdist = chains[:,2]
        tdist = chains[:,3]*24.*60.*60.
        pdist = (chains[:,4]-pdata[0,4])*24.*3600.
        q1dist = chains[:,5]
        q2dist = chains[:,6]
    else:
    '''
    print('... importing MCMC chains')
    chainfile = glob.glob(directory+targname+'*t0.txt')    
    tdist = np.loadtxt(chainfile[0])
    tdist = (tdist - param_info['t0'])*86400
    chainfile = glob.glob(directory+targname+'*rprs.txt')
    rdist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*aors.txt')
    adist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*inc.txt')
    idist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*q1.txt')
    q1dist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*q2.txt')
    q2dist = np.loadtxt(chainfile[0])
    chainfile = glob.glob(directory+targname+'*logf.txt')    
    logfdist = np.loadtxt(chainfile[0])

    #if lp is False:
    chainfile = glob.glob(directory+targname+'*logprob.txt')    
    lp = np.loadtxt(chainfile[0])

    print('... done importing chains!')
    print(done_in(tmaster))

#  Get maximum likelihood values
    print("Determining maximum likelihood values")
    maxlike = np.max(lp)
    imax = np.argmax(lp) 
    t0val = np.float(tdist[imax])
    rprsval = np.float(rdist[imax])
    aorsval = np.float(adist[imax])
    incval  = np.float(idist[imax])
    q1val   = np.float(q1dist[imax])
    q2val   = np.float(q2dist[imax])
    logfval = np.float(logfdist[imax])
        
    if thin:
        tdist = tdist[0::thin]
        rdist = rdist[0::thin]
        adist = adist[0::thin]
        idist = idist[0::thin]
        q1dist = q1dist[0::thin]
        q2dist = q2dist[0::thin]
        logfdist = logfdist[0::thin]
        lp     = lp[0::thin]

    nsamp = len(tdist)

    print("")
    print("Starting grid of posteriors...")
    plt.figure(6,figsize=(8.5,8.5))
    plt.clf()
    nx = 6
    ny = 6

    gs = gridspec.GridSpec(nx,ny,wspace=0.1,hspace=0.1)
    print("")
    print("... top plot of first column")
    tcol = time.time()
    top_plot(rdist,gs[0,0],val=rprsval,sigfac=sigfac,bindiv=bindiv)
    print(done_in(tcol))
    t = time.time()
    print("... first column plot")
    column_plot(rdist,adist,gs[1,0],val1=rprsval,val2=aorsval,ylabel=r'$a/R_\star$',sigfac=sigfac)
    print(done_in(t))
    column_plot(rdist,idist,gs[2,0],val1=rprsval,val2=incval,ylabel=r'$i$ (deg)',sigfac=sigfac,sigsamp=sigsamp)
    column_plot(rdist,tdist,gs[3,0],val1=rprsval,val2=t0val,ylabel=r'$t_0$ (s)',sigfac=sigfac,sigsamp=sigsamp)
    column_plot(rdist,q1dist,gs[4,0],val1=rprsval,val2=q1val,ylabel=r'$q_1$',sigfac=sigfac,sigsamp=sigsamp)
    corner_plot(rdist,q2dist,gs[5,0],val1=rprsval,val2=q2val,xlabel=r'$R_p/R_\star$',ylabel=r'$q_2$',
                sigfac=sigfac,sigsamp=sigsamp)
    print("First column: ")
    print(done_in(tcol))

    print("... starting second column")
    t2 = time.time()
    top_plot(adist,gs[1,1],val=aorsval,sigfac=sigfac,bindiv=bindiv)    
    middle_plot(adist,idist,gs[2,1],val1=aorsval,val2=incval,sigfac=sigfac,sigsamp=sigsamp)
    middle_plot(adist,tdist,gs[3,1],val1=aorsval,val2=t0val,sigfac=sigfac,sigsamp=sigsamp)
    middle_plot(adist,q1dist,gs[4,1],val1=aorsval,val2=q1val,sigfac=sigfac,sigsamp=sigsamp)
    row_plot(adist,q2dist,gs[5,1],val1=aorsval,val2=q2val,xlabel=r'$a/R_\star$',sigfac=sigfac,sigsamp=sigsamp)
    print(done_in(t2))

    print("... starting third column")
    t3 = time.time()
    top_plot(idist,gs[2,2],val=incval,sigfac=sigfac,bindiv=bindiv)    
    middle_plot(idist,tdist,gs[3,2],val1=incval,val2=t0val,sigfac=sigfac,sigsamp=sigsamp)
    middle_plot(idist,q1dist,gs[4,2],val1=incval,val2=q1val,sigfac=sigfac,sigsamp=sigsamp)
    row_plot(idist,q2dist,gs[5,2],val1=incval,val2=q2val,xlabel=r'$i$ (deg)',sigfac=sigfac,sigsamp=sigsamp)
    print(done_in(t3))

    print("... starting fourth column")
    t4 = time.time()
    top_plot(tdist,gs[3,3],val=t0val,sigfac=sigfac,bindiv=bindiv)    
    middle_plot(tdist,q1dist,gs[4,3],val1=t0val,val2=q1val,sigfac=sigfac,sigsamp=sigsamp)
    row_plot(tdist,q2dist,gs[5,3],val1=t0val,val2=q2val,xlabel=r'$t_0$ (s)',sigfac=sigfac,sigsamp=sigsamp)
    print(done_in(t4))

    print("... starting fifth column")
    t5 = time.time()
    top_plot(q1dist,gs[4,4],val=q1val,sigfac=sigfac,bindiv=bindiv)    
    row_plot(q1dist,q2dist,gs[5,4],val1=q1val,val2=q2val,xlabel=r'$q_1$',sigfac=sigfac,sigsamp=sigsamp)
    print(done_in(t5))

    print("... starting the last plot")
    t6 = time.time()
    top_plot(q2dist,gs[5,5],val=q2val,xlabel=r'$q_2$',sigfac=sigfac,bindiv=bindiv)    
    print(done_in(t6))

    plt.suptitle(targname+" Fit Posterior Distributions",fontsize=fs)

    print "Saving output figures"
    if pdf:
        ftype='.pdf'
    else:
        ftype='.png'
    plt.savefig(directory+targname+'_triangle'+ftype,dpi=300)
    plt.savefig(directory+targname+'_triangle.eps', format='eps',dpi=300)

    print("Procedure finished!")
    print(done_in(tmaster))

    gc.collect()

    mpl.rcdefaults()

    return



def top_plot(dist,position,val=False,sigfac=4.0,minval=False,maxval=False,bindiv=20,aspect=1,xlabel=False):
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from scipy.stats.kde import gaussian_kde
    from scipy.interpolate import interp1d
    import time

#    len = np.size(dist)
#    min = np.float(np.sort(dist)[np.round(frac*len)])
#    max = np.float(np.sort(dist)[np.round((1.-frac)*len)])
#    dists = np.linspace(np.min(dist)*0.5,np.max(dist)*1.5,1000)
#    kde = gaussian_kde(dist)
#    pdf = kde(dists)
#    cumdist = np.cumsum(pdf)/np.sum(pdf)
#    func = interp1d(cumdist,dists,kind='linear')
#    lo = np.float(func(math.erfc(1./np.sqrt(2))))
#    hi = np.float(func(math.erf(1./np.sqrt(2))))
    med = np.median(dist)
#    pdb.set_trace()
    sig = mad_std(dist)
    if sig < 1.0e-5 or np.abs(med) < 1.0e-10:
        nb = 10
    else:
        datamin = np.min(dist)
        datamax = np.max(dist)
        if not minval:
            minval = max(med - sigfac*sig,datamin)
        if not maxval:
            maxval = min(med + sigfac*sig,datamax)
        nb = np.int(np.round(np.ceil( (np.max(dist)-np.min(dist))*bindiv / np.abs(maxval-minval))))

    ax = plt.subplot(position)
    plt.hist(dist,bins=nb,normed=True,color='black')
    if not xlabel: 
        ax.set_xticklabels(())
    ax.set_yticklabels(())
    ax.set_xlim(min(minval,maxval),max(minval,maxval))
    xlimits = ax.get_xlim()
    ylimits = ax.get_ylim()
#    print "x range for top plot:"
#    print xlimits
    ax.set_aspect(abs((xlimits[1]-xlimits[0])/(ylimits[1]-ylimits[0]))/aspect)
    if val:
        pass
#        plt.axvline(x=val,color='w',linestyle='--',linewidth=2)
    if xlabel:
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(8) 
            tick.label.set_rotation('vertical')
        ax.set_xlabel(xlabel,fontsize=12)
 
    plt.gca().yaxis.set_major_locator(plt.NullLocator())

    return


def column_plot(dist1,dist2,position,val1=False,val2=False,sigfac=2.0,sigsamp=5.0,
                min1=False,max1=False,min2=False,max2=False,ylabel=None):
    from statsmodels.nonparametric.kernel_density import KDEMultivariate as KDE

#    len1 = np.size(dist1)
#    min1 = np.float(np.sort(dist1)[np.round(frac*len1)])
#    max1 = np.float(np.sort(dist1)[np.round((1.-frac)*len1)])
    med1 = np.median(dist1)
    sig1 = mad_std(dist1)
    datamin1 = np.min(dist1)
    datamax1 = np.max(dist1)
    if not min1:
        min1 = max(med1 - sigfac*sig1,datamin1)
    if not max1:
        max1 = min(med1 + sigfac*sig1,datamax1)

#    len2 = np.size(dist2)
#    min2 = np.float(np.sort(dist2)[np.round(frac*len2)])
#    max2 = np.float(np.sort(dist2)[np.round((1.-frac)*len2)])
    med2 = np.median(dist2)
    sig2 = mad_std(dist2)
    datamin2 = np.min(dist2)
    datamax2 = np.max(dist2)
    if not min2:
        min2 = max(med2 - sigfac*sig2,datamin2)
    if not max2:
        max2 = min(med2 + sigfac*sig2,datamax2)

    aspect = (max1-min1)/(max2-min2)
    X, Y = np.mgrid[min1:max1:100j, min2:max2:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([dist1, dist2])

    kernel = KDE(values,var_type='cc',bw=[sig1/sigsamp,sig2/sigsamp])
    Z = np.reshape(kernel.pdf(positions).T, X.shape)

#    kernel = gaussian_kde(values)
#    Z = np.reshape(kernel(positions).T, X.shape)
 
    ax = plt.subplot(position)
    ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,aspect=aspect,\
              extent=[min1, max1, min2, max2],origin='upper')
    clev = np.exp(np.log(np.max(Z))-0.5)
    cset = ax.contour(X,Y,Z,[clev],colors='w',linewidth=5,linestyles='dotted')
#    print "x range for column plot:"
#    print min1,max1
#    ax.plot(val1,val2, 'wx', markersize=3)
#    ax.set_xlim(min1, max1)
#    ax.set_ylim(min2, max2)
    ax.set_xticklabels(())
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(8) 
    ax.set_ylabel(ylabel,fontsize=12)

    return



def row_plot(dist1,dist2,position,val1=False,val2=False,sigfac=3.5,sigsamp=5.0,
             min1=False,max1=False,min2=False,max2=False,xlabel=None):
    from statsmodels.nonparametric.kernel_density import KDEMultivariate as KDE


#    len1 = np.size(dist1)
#    min1 = np.float(np.sort(dist1)[np.round(frac*len1)])
#    max1 = np.float(np.sort(dist1)[np.round((1.-frac)*len1)])
    med1 = np.median(dist1)
    sig1 = mad_std(dist1)
    datamin1 = np.min(dist1)
    datamax1 = np.max(dist1)
    if not min1:
        min1 = max(med1 - sigfac*sig1,datamin1)
    if not max1:
        max1 = min(med1 + sigfac*sig1,datamax1)

    
#    len2 = np.size(dist2)
#    min2 = np.float(np.sort(dist2)[np.round(frac*len2)])
#    max2 = np.float(np.sort(dist2)[np.round((1.-frac)*len2)])
    med2 = np.median(dist2)
    sig2 = mad_std(dist2)
    datamin2 = np.min(dist2)
    datamax2 = np.max(dist2)
    if not min2:
        min2 = max(med2 - sigfac*sig2,datamin2)
    if not max2:
        max2 = min(med2 + sigfac*sig2,datamax2)
    
    aspect = (max1-min1)/(max2-min2)
    X, Y = np.mgrid[min1:max1:100j, min2:max2:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([dist1, dist2])

    kernel = KDE(values,var_type='cc',bw=[sig1/sigsamp,sig2/sigsamp])
    Z = np.reshape(kernel.pdf(positions).T, X.shape)

#    kernel = gaussian_kde(values)
#    Z = np.reshape(kernel(positions).T, X.shape)

    ax = plt.subplot(position)
    ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,aspect=aspect,\
              extent=[min1, max1, min2, max2],origin='upper')
    clev = np.exp(np.log(np.max(Z))-0.5)
    cset = ax.contour(X,Y,Z,[clev],colors='w',linewidth=5,linestyles='dotted')
#    ax.plot(val1,val2, 'wx', markersize=3)
    ax.set_xlim(min1, max1)
    ax.set_ylim(min2, max2)
    ax.set_yticklabels(())
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8) 
        tick.label.set_rotation('vertical')
    ax.set_xlabel(xlabel,fontsize=12)
    return


def middle_plot(dist1,dist2,position,val1=False,val2=False,sigfac=3.5,sigsamp=5.0,
                min1=False,max1=False,min2=False,max2=False):
    from statsmodels.nonparametric.kernel_density import KDEMultivariate as KDE


#    len1 = np.size(dist1)
#    min1 = np.float(np.sort(dist1)[np.round(frac*len1)])
#    max1 = np.float(np.sort(dist1)[np.round((1.-frac)*len1)])
    med1 = np.median(dist1)
    sig1 = mad_std(dist1)
    datamin1 = np.min(dist1)
    datamax1 = np.max(dist1)
    if not min1:
        min1 = max(med1 - sigfac*sig1,datamin1)
    if not max1:
        max1 = min(med1 + sigfac*sig1,datamax1)

#    len2 = np.size(dist2)
#    min2 = np.float(np.sort(dist2)[np.round(frac*len2)])
#    max2 = np.float(np.sort(dist2)[np.round((1.-frac)*len2)])
    med2 = np.median(dist2)
    sig2 = mad_std(dist2)
    datamin2 = np.min(dist2)
    datamax2 = np.max(dist2)
    if not min2:
        min2 = max(med2 - sigfac*sig2,datamin2)
    if not max2:
        max2 = min(med2 + sigfac*sig2,datamax2)

    aspect = (max1-min1)/(max2-min2)
    X, Y = np.mgrid[min1:max1:100j, min2:max2:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([dist1, dist2])

    kernel = KDE(values,var_type='cc',bw=[sig1/sigsamp,sig2/sigsamp])
    Z = np.reshape(kernel.pdf(positions).T, X.shape)

#    kernel = gaussian_kde(values)
#    Z = np.reshape(kernel(positions).T, X.shape)
 

    ax = plt.subplot(position)
    ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,aspect=aspect,\
              extent=[min1, max1, min2, max2],origin='upper')
    clev = np.exp(np.log(np.max(Z))-0.5)
    cset = ax.contour(X,Y,Z,[clev],colors='w',linewidth=5,linestyles='dotted')
#    ax.plot(val1,val2, 'wx', markersize=3)
    ax.set_xlim(min1, max1)
    ax.set_ylim(min2, max2)
    ax.set_xticklabels(())
    ax.set_yticklabels(())
    return



def corner_plot(dist1,dist2,position,val1=False,val2=False,sigfac=3.5,sigsamp=5.0,
                min1=False,max1=False,min2=False,max2=False,
                xlabel=None,ylabel=None):

    from statsmodels.nonparametric.kernel_density import KDEMultivariate as KDE

    med1 = np.median(dist1)
    sig1 = mad_std(dist1)
    datamin1 = np.min(dist1)
    datamax1 = np.max(dist1)
    if not min1:
        min1 = max(med1 - sigfac*sig1,datamin1)
    if not max1:
        max1 = min(med1 + sigfac*sig1,datamax1)

    
#    len2 = np.size(dist2)
#    min2 = np.float(np.sort(dist2)[np.round(frac*len2)])
#    max2 = np.float(np.sort(dist2)[np.round((1.-frac)*len2)])
    med2 = np.median(dist2)
    sig2 = mad_std(dist2)
    datamin2 = np.min(dist2)
    datamax2 = np.max(dist2)
    if not min2:
        min2 = max(med2 - sigfac*sig2,datamin2)
    if not max2:
        max2 = min(med2 + sigfac*sig2,datamax2)

    aspect = (max1-min1)/(max2-min2)
    X, Y = np.mgrid[min1:max1:100j, min2:max2:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([dist1, dist2])

    kernel = KDE(values,var_type='cc',bw=[sig1/sigsamp,sig2/sigsamp])
    Z = np.reshape(kernel.pdf(positions).T, X.shape)

#    kernel = gaussian_kde(values)
#    Z = np.reshape(kernel(positions).T, X.shape) 

    ax = plt.subplot(position)
    ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,aspect=aspect,\
              extent=[min1, max1, min2, max2],origin='upper')
    clev = np.exp(np.log(np.max(Z))-0.5)
    cset = ax.contour(X,Y,Z,[clev],colors='w',linewidth=5,linestyles='dotted')
#    ax.plot(val1,val2, 'wx', markersize=3)
    ax.set_xlim(min1, max1)
    ax.set_ylim(min2, max2)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(8) 
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8) 
        tick.label.set_rotation('vertical')
    ax.set_xlabel(xlabel,fontsize=12)
    ax.set_ylabel(ylabel,fontsize=12)

    return

def get_limb_curve(ldc,limbmodel='quad'):


    """
    get_limb_curve:
    ---------------
    Function to compute limb darkening curve given models and parameters

    """

    gamma = np.linspace(0,np.pi/2.0,1000,endpoint=True)
    theta = gamma*180.0/np.pi
    mu = np.cos(gamma)
    
    if limbmodel == 'nlin':
        c1 = ldc[0]
        c2 = ldc[1]
        c3 = ldc[2]
        c4 = ldc[3]
        Imu = 1.0 - c1*(1.0 - mu**0.5) - c2*(1.0 - mu) - \
              c3*(1.0 - mu**1.5) - c4*(1.0 - mu**2.0)
    elif limbmodel == 'quad':
        c1 = ldc[0]
        c2 = ldc[1]
        Imu = 1.0 - c1*(1.0 - mu) - c2*(1.0 - mu)**2.0
    elif limbmodel == 'sqrt':
        c1 = ldc[0]
        c2 = ldc[1]
        Imu = 1.0 - c1*(1.0 - mu) - c2*(1.0 - mu**0.5)
    else: pass

    return theta, Imu



def plot_limb_curves(ldc=False,limbmodel='quad',write=False,network=None):

    """
    
    plot_limb_curves:
    -----------------
    

    """
    import constants as c
    plt.rcdefaults()

    if network == 'koi':
        net = None
    else:
        net = network


    Mstar = sdata[0,0]
    Rstar = sdata[0,1]
    Tstar = sdata[0,2]

    loggstar = np.log10( c.G * Mstar / Rstar**2. )
    
    a1,a2,a3,a4 = get_limb_coeff(Tstar,loggstar,limb='nlin',interp='nearest',network=net)
    a,b = get_limb_coeff(Tstar,loggstar,limb='quad',interp='nearest',network=net)
    c,d = get_limb_coeff(Tstar,loggstar,limb='sqrt',interp='nearest',network=net)

    thetaq,Imuq = get_limb_curve([a,b],limbmodel='quad')
    thetas,Imus = get_limb_curve([c,d],limbmodel='sqrt')
    thetan,Imun = get_limb_curve([a1,a2,a3,a4],limbmodel='nlin')
    if ldc:
        thetain,Iin = get_limb_curve(ldc,limbmodel=limbmodel)

    if write:
        plt.figure(1,figsize=(11,8.5))
    else:
        plt.ion()
        plt.figure()
        plt.plot(thetaq,Imuq,label='Quadratic LD Law')
        plt.plot(thetas,Imus,label='Root-Square LD Law')
        plt.plot(thetan,Imun,label='Non-Linear LD Law')

    if ldc:
        if limbmodel == 'nlin':
            label = '{0:0.2f}, {1:0.2f}, {2:0.2f}, {3:0.2f} ('+limbmodel+')'
            plt.plot(thetain,Iin,label=label.format((ldc[0],ldc[1],ldc[2],ldc[3])))
        else:
            label = '%.2f, ' % ldc[0] + '%.2f' % ldc[1]+' ('+limbmodel+')'
            plt.plot(thetain,Iin,label=label.format((ldc[0],ldc[1])))

    plt.ylim([0,1.0])
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xlabel(r"$\theta$ (degrees)",fontsize=18)
    plt.ylabel(r"$I(\theta)/I(0)$",fontsize=18)
    plt.title("KOI-"+str(koi)+" limb darkening",fontsize=20)
    plt.legend(loc=3)
    plt.annotate(r'$T_{\rm eff}$ = %.0f K' % sdata[0][2], [0.86,0.82],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
    plt.annotate(r'$\log(g)$ = %.2f' % loggstar, [0.86,0.77],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
    
    if write:
        plt.savefig(path+'MCMC/'+str(koi)+stag+ctag+rtag+lptag+'_'+limbmodel+'.png')
        plt.clf()
    else:
        plt.ioff()

    return




def get_limb_spread(q1s,q2s,sdata=None,factor=1,limbmodel='quad',
                    fontsz=18,write=False,plot=True,network=None):

    """
    
    get_limb_spread:
    -----------------

    To Do:
    ------
    Write out distributions of q1 values so that it does not have
    to be done every refit.

    """
    plt.rcdefaults()

    if sdata == None:
        print "Must supply stellar data!"
        return
    
    if limbmodel == 'quad':
        lname = 'Quadratic'
    if limbmodel == 'sqrt':
        lname = 'Root Square'
    if limbmodel == 'nlin':
        lname = '4 Parameter'

    Mstar  = sdata[0][0]
    eMstar = sdata[1][0]/c.Msun * factor

    Rstar  = sdata[0][1]
    eRstar = sdata[1][1]/c.Rsun * factor

    Tstar  = sdata[0][2]
    eTstar = sdata[1][2] * factor

    loggstar = np.log10( c.G * Mstar / Rstar**2. )


    if plot:
        if write:
            plt.figure(101)
        else:
            plt.ion()
            plt.figure(123)
            plt.clf()

    sz = len(q1s)
    for i in range(sz):
        u1,u2 = qtou(q1s[i],q2s[i],limb=limb)
        theta,Imu = get_limb_curve([u1,u2],limbmodel=limbmodel)
        plt.plot(theta,Imu,lw=0.1,color='blue')
        
    plt.ylim([0,1.4])
    plt.tick_params(axis='both', which='major', labelsize=fontsz-2)
    plt.xlabel(r"$\theta$ (degrees)",fontsize=fontsz)
    plt.ylabel(r"$I(\theta)/I(0)$",fontsize=fontsz)
    plt.title("KOI-"+str(koi)+" limb darkening prior distribution",fontsize=fontsz)
#    plt.legend(loc=3)

    plt.annotate(r'$\Delta T_{\rm eff}$ = %.0f K' % eTstar, [0.86,0.82],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
    plt.annotate(r'$\Delta M_\star$ = %.2f M$_\odot$' % eMstar, [0.86,0.77],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
    plt.annotate(r'$\Delta R_\star$ = %.2f R$_\odot$' % eRstar, [0.86,0.72],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
 
    plt.annotate(r'$T_{\rm eff}$ = %.0f K' % sdata[0][2], [0.16,0.82],horizontalalignment='left',
                 xycoords='figure fraction',fontsize='large')
    plt.annotate(r'$\log(g)$ = %.2f' % loggstar, [0.16,0.77],horizontalalignment='left',
                 xycoords='figure fraction',fontsize='large')
    plt.annotate(lname, [0.16,0.72],horizontalalignment='left',
                 xycoords='figure fraction',fontsize='large')
    

    if write:
        directory = path+'MCMC/'
        if not os.path.exists(directory):
            os.makedirs(directory)
        plt.savefig(directory+str(koi)+stag+ctag+rtag+lptag+'_'+limbmodel+'fit_LDspread.png')
        plt.clf()


    plt.hist(q1s[~np.isnan(q1s)],bins=sz/70,normed=True,label=r'$q_1$')
    plt.hist(q2s[~np.isnan(q2s)],bins=sz/70,normed=True,label=r'$q_2$')
    plt.tick_params(axis='both', which='major', labelsize=fontsz-2)
    plt.title('Distribution of Kipping $q$ values',fontsize=fontsz)
    plt.xlabel(r'$q$ value',fontsize=fontsz)
    plt.ylabel('Normalized Frequency',fontsize=fontsz)
    plt.legend(loc='upper right',prop={'size':fontsz-2},shadow=True)
    plt.xlim(0,1)

    if write:
        plt.savefig(directory+str(koi)+stag+ctag+rtag+lptag+'_'+limbmodel+'qdist.png',dpi=300)
        plt.clf()
    else:
        plt.ioff()

 
    return 



def thin_chains(koi,planet,thin=10,short=False,network=None,clip=False,limbmodel='quad',rprior=False):
    
    lc,pdata,sdata = get_koi_info(koi,planet,short=short,network=network,\
                                      clip=clip,limbmodel=limbmodel,rprior=rprior)
    t = time.time()
    print 'Importing MCMC chains'
    rprsdist = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_rchain.txt')
    ddist    = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_dchain.txt')*24.
    bdist    = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_bchain.txt')
    tdist    = (np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_t0chain.txt')) * 24.0 * 3600.0 + pdata[0,3]
    pdist    = (np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_pchain.txt')-pdata[0,4])*24.*3600.
    q1dist   = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_q1chain.txt')
    q2dist   = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_q2chain.txt')
    print done_in(t)

    rprsdist = rprsdist[0::thin]
    ddist    = ddist[0::thin]
    tdist    = tdist[0::thin]
    bdist    = bdist[0::thin]
    pdist    = pdist[0::thin]
    q1dist   = q1dist[0::thin]
    q2dist   = q2dist[0::thin]

    t = time.time()
    print 'Exporting thinned chains'
    np.savetxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_thin_rchain.txt',rdist)
    np.savetxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_thin_dchain.txt',ddist)
    np.savetxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_thin_bchain.txt',bdist)
    np.savetxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_thin_t0chain.txt',tdist)
    np.savetxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_thin_pchain.txt',pdist)
    np.savetxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_thin_q1chain.txt',q1dist)
    np.savetxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_thin_q2chain.txt',q2dist)
    print done_in(t)
    
    return

def get_rhostar(koi,planet,short=False,network=None,clip=False,limbmodel='quad',
                rprior=False,lprior=False,thin=100,bins=100):
    plt.rcdefaults()

    lc,pdata,sdata = get_koi_info(koi,planet,short=short,network=network,\
                                  clip=clip,limbmodel=limbmodel,rprior=rprior,lprior=lprior)
    t = time.time()

    print 'Importing MCMC chains'
    rprsdist = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_rchain.txt')
    print done_in(t)
    ddist    = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_dchain.txt')
    bdist    = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_bchain.txt')
    pdist    = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_pchain.txt')
    q1dist   = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_q1chain.txt')
    q2dist   = np.loadtxt(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_q2chain.txt')
    print done_in(t)

    if thin:
        rprsdist = rprsdist[0::thin]
        ddist = ddist[0::thin]
        pdist = pdist[0::thin]
        bdist = bdist[0::thin]
        q1dist = q1dist[0::thin]
        q2dist = q2dist[0::thin]


#    ttotdist,tfulldist = midtotot(ddist,rprsdist,bdist)

    ttotdist = ddist

#    t = time.time()
#    tfulldist =[]
#    for i in np.arange(len(rprsdist)):
#        tmodel,smoothmodel,tcontact = compute_trans(rprsdist[i],ddist[i],bdist[i],0.0,
#                                                 pdist[i],[q1dist[i],q2dist[i]],
#                                                 getcontact=True)
#        ttot=tcontact[0]
#        tfull=tcontact[1]
#
#        tfulldist = np.append(tfulldist,tfull)
#        ttotdist  = np.append(ttotdist,ttot)
#        if (i+1) % 500 == 0:
#            out = done_in(t)
#            pdone = (i+1.0)/len(rprsdist) * 100.0
#            print str(pdone)+"% "+out

    aorsdist = np.sqrt( ((1+rprsdist)**2 - bdist**2)/((np.sin(np.pi*ttotdist/pdist))**2) + bdist**2)
#    aorsdist = np.sqrt( ((1-rprsdist)**2 - bdist**2)/((np.sin(np.pi*tfulldist/pdist))**2) + bdist**2)

#    aorsdist = 2.0 * np.sqrt(rprsdist) * pdist / (np.pi*np.sqrt(ttotdist**2 - tfulldist**2))
    
    rhodist =  3.0*np.pi/( c.G * (pdist*24.*3600.)**2 ) * aorsdist**3

    plt.figure(999,figsize=(8.5,11))
    plt.clf()
    plt.subplot(2,1,1)
    
    plt.hist(aorsdist,bins=bins)
    nsamp = len(aorsdist)
    frac = 0.003
    aorsmin = np.float(np.sort(aorsdist)[np.round(frac*nsamp)])
    aorsmax = np.float(np.sort(aorsdist)[np.round((1-frac)*nsamp)])
    plt.xlim([aorsmin,aorsmax])
    plt.xlim([10.8,11.1])
#    plt.axvline(x=q1val,color='r',linestyle='--')
#    plt.axvline(x=q1lo,color='c',linestyle='--')
#    plt.axvline(x=q1med,color='c',linestyle='--')
#    plt.axvline(x=q1hi,color='c',linestyle='--')
    plt.xlabel(r'$a/R_{\star}$')
    plt.ylabel(r'$N$')

    plt.subplot(2,1,2)
    
    plt.hist(rhodist,bins=bins)
    nsamp = len(rhodist)
    frac = 0.003
    rhomin = np.float(np.sort(rhodist)[np.round(frac*nsamp)])
    rhomax = np.float(np.sort(rhodist)[np.round((1-frac)*nsamp)])
    plt.xlim([rhomin,rhomax])
    plt.xlim([3.95,4.3])
#    plt.axvline(x=q1val,color='r',linestyle='--')
#    plt.axvline(x=q1lo,color='c',linestyle='--')
#    plt.axvline(x=q1med,color='c',linestyle='--')
#    plt.axvline(x=q1hi,color='c',linestyle='--')
    plt.xlabel(r'$\rho_{\star}$')
    plt.ylabel(r'$N$')
    
    plt.savefig(directory+targname+stag+ctag+rtag+lptag+ltag+'fit_rhostar.png')

    pdb.set_trace()
    return rhodist

def do_fit(koi,planet,nwalkers=1000,burnsteps=1000,mcmcsteps=1000,clobber=False,
           network=None,thin=50,tthin=100,sigfac=3.5,clip=True,limbmodel='quad',doplots=False,
           short=False,rprior=True,lprior=False,notriangle=True,getsamp=True,
           errfac=3,rpmax=0.5,durmax=10,bindiv=75,sigsamp=5.0,pdf=False):
    
    import numpy as np
    import time
    import os
    import constants as c

    print ""
    print ""    
    print "Starting MCMC fitting for KOI-"+str(koi)+".0"+str(planet)

    check1 = isthere(koi,planet,short=short,network=network,clip=clip)
    if check1:
        lc,pdata,sdata = get_koi_info(koi,planet,short=short,network=network, \
                                      clip=clip,limbmodel=limbmodel,getsamp=getsamp, \
                                      rprior=rprior,lprior=lprior,errfac=errfac)
        chains,lp = fit_single(nwalkers=nwalkers,burnsteps=burnsteps, \
                                   mcmcsteps=mcmcsteps,clobber=clobber)
        if doplots != False:
            fit,rpdf = bestvals(chains=chains,lp=lp,thin=thin,rpmax=rpmax,
                                durmax=durmax,bindiv=bindiv,sigfac=sigfac,pdf=pdf)
            if notriangle != True:
                triangle_plot(thin=tthin,sigfac=sigfac,bindiv=bindiv,sigsamp=sigsamp,pdf=pdf)

    else: pass
 
    return
