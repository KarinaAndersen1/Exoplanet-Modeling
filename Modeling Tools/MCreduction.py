### Imports ###

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### Accessing Data ###

data = pd.read_csv('C:\Users\karina\Desktop\Linefit\Fit_Line\make_line_data.csv')
x = np.array(data[['x']])
y = np.array(data[['y']])
yerr = np.array(data[['yerr']]) # Arrays of x, y, and yerr values

x0 = data['x'].tolist()
y0 = data['y'].tolist()
yerr0 = data['yerr'].tolist() # Lists of x, y, and yerr values to be used with plt

### Maximum Likelihood ###

def log_likelihood(theta, x, y, yerr):
    m, b = theta
    model = m * x + b
    sigma2 = yerr ** 2
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2)) # LL equation

from scipy.optimize import minimize

nll = lambda *args: -log_likelihood(*args)
initial = np.array([1, 2])
soln = minimize(nll, initial, args=(x, y, yerr)) # Make sure nll and initial are both arrays!
m_ml, b_ml = soln.x

print("Maximum likelihood estimates:")
print("m = {0:.3f}".format(m_ml))
print("b = {0:.3f}".format(b_ml))

plt.errorbar(x0, y0, yerr=yerr0, fmt=".k", capsize=0)
plt.plot(x0, np.dot(np.vander(x0, 2), [m_ml, b_ml]), ":k", label="ML")
plt.legend(fontsize=14)
plt.xlim(0, 10)
plt.xlabel("x")
plt.ylabel("y")
plt.show() # Plots LL solution

### MCMC ###

def log_prior(theta):
    m, b = theta
    if -10.0 < m < 10.0 and -10.0 < b < 10.0:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

import emcee

pos = soln.x + 1e-4 * np.random.randn(32, 2)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, yerr))
sampler.run_mcmc(pos, 5000, progress=True) # Change the number of steps if wanted

fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["m", "b"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
    plt.show() # Did not work correctly, plots were plotted but they aren't what they should be







fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["m", "b"]
ax = axes[0]
ax.plot(samples[:, :, 0], "k", alpha=0.3)
ax.set_xlim(0, len(samples))
ax.set_ylabel(labels[0])
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.show()

fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["m", "b"]

ax = axes[1]
ax.plot(samples[:, :, 1], "k", alpha=0.3)
ax.set_xlim(0, len(samples))
ax.set_ylabel(labels[1])
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.show()
    
axes[-1].set_xlabel("step number")

tau = sampler.get_autocorr_time()
print(tau) # What is tau?

flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
print(flat_samples.shape)

### Corner Plot ###

import corner

fig = corner.corner(
    flat_samples, labels=labels)

inds = np.random.randint(len(flat_samples), size=100)
for ind in inds:
    sample = flat_samples[ind]
    plt.plot(x0, np.dot(np.vander(x0, 2), sample[:2]), "C1", alpha=0.1)
plt.errorbar(x0, y0, yerr=yerr0, fmt=".k", capsize=0)
plt.legend(fontsize=14)
plt.xlim(0, 10)
plt.xlabel("x")
plt.ylabel("y")
plt.show() # Plotted, but doesn't look like what it should be

from IPython.display import display, Math

for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    print mcmc[1]
    print q[0]
    print q[1]
    print labels[i]
    txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
    txt = txt.format(mcmc[1], q[0], q[1], labels[i])
    print Math(txt)
