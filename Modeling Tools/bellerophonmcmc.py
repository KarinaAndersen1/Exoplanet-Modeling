''' 
Bellerophon version of mcmc script, just run this then 
transit_mcmc(theta, t, f, ferr, param_info) and then 
bestvals() to get specific parameters

Running into issue using log_likelihood_standard; works 
but gives bad values with regular log_likelihood
'''


# Imports #
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
import corner



# Functions for generating the model given the parameters #
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

def compute_trans(theta, t, param_info):
    t0,P,rprs,aors,inc,ecc,omega,q1,q2 = theta
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
                            exp_time = param_info['exp_time'])    #initializes mode
    model = m.light_curve(params)
    return model



# Target info #
targname = 'TIC296780789.01'
obstype = 'Transit'
rastr = '09 29 23.33'
decstr = '-14 30 41.0'
date = '20200225'

#targname = 'TIC125442121.01'
#obstype = 'Transit'
#rastr = '23 26 19.13'
#decstr = '+38 33 10.7'
#date = '20191112'

coords = SkyCoord(rastr,decstr,unit=(u.hour,u.deg))
ra  = coords.ra.deg
dec = coords.dec.deg

ingress = 2458904.753
midtrans = 2458904.789
egress = 2458904.824

paths = pp.get_paths(targname=targname,obstype=obstype)


# Note: this takes a csv file with 3 columns: 'BJD', 'flux', and 'flux_error' which must be created and put in the right folder prior to running #
data = pd.read_csv(paths['output']+'/modelfit_measurements.csv')

t = data['BJD'].values
f = data['flux'].values
ferr = data['flux_error'].values

exp_time = 20./(86400.)

param_info = {'t0':2458532,'P':13.58 ,'rprs':np.sqrt(0.24) ,
              'aors':60.0 , 'inc':89.0 ,'ecc':0.3, 'omega':0.0,
              'q1':0.3 ,'q2':0.3 , 'exp_time':exp_time, 'supersample_factor':7,
              'limb_model':'quadratic'}


theta = [param_info['t0'],param_info['P'],param_info['rprs'],param_info['aors'],
         param_info['inc'],param_info['ecc'],param_info['omega'],param_info['q1'],param_info['q2']]

# Computing and plotting the preliminary model #
model = compute_trans(theta, t, param_info)

ffit = np.append(f[0:9],f[-20:])
fmean  = np.mean(ffit)
fnorm = f/fmean

plt.plot(t, fnorm, 'bo', zorder=1)
plt.plot(t, model, 'r-', zorder=2)
plt.show()


# Defining functions needed for fitting #
def log_prior(theta):

    t0,P,rprs,aors,inc,ecc,omega,q1,q2 = theta
 
    if rprs > 0.8 or rprs < 0:
        return -np.inf

    b = aors*np.cos(np.radians(inc))*(1-ecc**2)/(1+ecc*np.sin(np.radians(omega)))
    
    if b > (1 + rprs) or b < 0:
        return -np.inf

    if aors < 0.5:
        return -np.inf

    if ecc > 1 or ecc < 0:
        return -np.inf

    if q1<0 or q1>1 or q2<0 or q2>1:
        return -np.inf

    if P < 5 or P > 20:
        return -np.inf

    if inc > 90.0 or inc < 85.0:
        return -np.inf
    
    return 0.0

theta = [param_info['t0'],param_info['rprs'],param_info['aors'], param_info['inc'],param_info['q1'],param_info['q2']] 
         
def log_likelihood_standard(theta, t, f, ferr, param_info): # ValueError: too many values to unpack when called from transit_mcmc

    t0,rprs,aors,inc,q1,q2 = theta
    
    P = param_info['P']
    ecc = param_info['ecc']
    omega = param_info['omega']
    
    theta = t0,P,rprs,aors,inc,ecc,omega,q1,q2
    
    model = compute_trans(theta,t,param_info)
    
    sigma2 = ferr ** 2 

    loglike =  -0.5 * np.sum((f - model) ** 2 / sigma2 + np.log(2*np.pi*sigma2))

    logp = log_prior(theta)

    return loglike + logp


# The function that takes everything and does everything #
def transit_mcmc(t, f, ferr, param_info):
    initial = [param_info['t0'],param_info['rprs'],param_info['aors'], param_info['inc'],param_info['q1'],param_info['q2']] 
    nll = lambda *args: -log_likelihood_standard(*args)
    soln = minimize(nll, initial, args=(t, f, ferr, param_info))
    t0_ml, rprs_ml, aors_ml, inc_ml, q1, q2 = soln.x
    
    pos = soln.x + 1e-4 * np.random.randn(100, 6)
    nwalkers, ndim = pos.shape
    mcs = 2000

    directory = paths['output']+'/MCMC/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood_standard, args=(t, f, ferr, param_info))
    sampler.run_mcmc(pos, mcs, progress=True)
    
    fig, axes = plt.subplots(6, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    labels = ["t0", "rprs", "aors", "inc", "q1", "q2"]
    
    for i in range(ndim):   # Plotting the chains
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
    axes[-1].set_xlabel("step number")
    plt.savefig(directory+'final_chains.png',dpi=300)
    plt.show()

    tdist  = sampler.flatchain[:,0]
    rdist  = sampler.flatchain[:,1]
    adist  = sampler.flatchain[:,2]
    idist  = sampler.flatchain[:,3]
    q1dist  = sampler.flatchain[:,4]
    q2dist  = sampler.flatchain[:,5]

    np.savetxt(directory+targname+'_'+str(np.int(nwalkers))+'walkers_'+str(np.int(mcs))+'steps_t0.txt',tdist)
    np.savetxt(directory+targname+'_'+str(np.int(nwalkers))+'walkers_'+str(np.int(mcs))+'steps_rprs.txt',rdist)
    np.savetxt(directory+targname+'_'+str(np.int(nwalkers))+'walkers_'+str(np.int(mcs))+'steps_aors.txt',adist)
    np.savetxt(directory+targname+'_'+str(np.int(nwalkers))+'walkers_'+str(np.int(mcs))+'steps_inc.txt',idist)
    np.savetxt(directory+targname+'_'+str(np.int(nwalkers))+'walkers_'+str(np.int(mcs))+'steps_q1.txt',q1dist)
    np.savetxt(directory+targname+'_'+str(np.int(nwalkers))+'walkers_'+str(np.int(mcs))+'steps_q2.txt',q2dist)


    lp = sampler.lnprobability.flatten()
    np.savetxt(directory+targname+'_'+str(np.int(nwalkers))+'walkers_'+str(np.int(mcs))+'steps_logprob.txt',lp)

    flat_samples = sampler.get_chain(discard=1000, thin=15, flat=True)
    print(flat_samples.shape)

    fig = corner.corner(  # Plotting the corner plots
        flat_samples, labels=labels)
    plt.show()

    return lp

def bestvals():     # Parsed routine from thacher_transit to get the best values
    directory = paths['output']+'/MCMC/'
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
    chainfile = glob.glob(directory+targname+'*logprob.txt')    
    lp = np.loadtxt(chainfile[0])


#  Get maximum likelihood values
    imax = np.argmax(lp) 
    t0val = np.float(tdist[imax])
    rprsval = np.float(rdist[imax])
    aorsval = np.float(adist[imax])
    incval  = np.float(idist[imax])
    q1val   = np.float(q1dist[imax])
    q2val   = np.float(q2dist[imax])

    return t0val, rprsval, aorsval, incval, q1val, q2val
