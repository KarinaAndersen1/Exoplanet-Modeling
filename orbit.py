# Imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
import emcee
import corner

# Read data file
data = pd.read_csv('C:\\Users\\karina\\Astronomy\\mid_transit_times4.csv')
x = data['Epoch'].values
y = data['Deviation_From_Period(min)'].values
bjd = data['BJD'].values
yerr = data['Error'].values

# Starting conditions
P_0 = 3.0395803
dp_dt_0 = 0
t0 = 2455186.01879 # Remains constant

'''
Line fit 
'''
# Log likelihood for a constant period
def log_likelihood_line(theta, x, y, yerr):
    P = theta
    model_line = (bjd - (t0 + P * x)) * 1440
    sigma2_line = yerr ** 2 + model_line ** 2
    return -0.5 * np.sum((y - model_line) ** 2 / sigma2_line + np.log(sigma2_line))

np.random.seed(42)
nll_line = lambda *args: -log_likelihood_line(*args)
initial_line = np.array([P_0]) + 0.1 * np.random.randn(1)
soln_line = minimize(nll_line, initial_line, args=(x, y, yerr))
P_line = soln_line.x

print("Maximum likelihood estimates:")
print(f'P = {P_line}')

# Plotting the results
x0 = np.linspace(0, 1500, 1500)
t_tra_con = t0 + P_line * x0
y_line = (t_tra_con - (t0 + P_line * x0)) * 1440 # Weird conversion to Deviation from Period

plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
plt.plot(x0, y_line, ":k", label="Line")
plt.legend(fontsize=14)
plt.xlim(0, 1500)
plt.xlabel("x")
plt.ylabel("y")

def log_prior(theta):
    P, dp_dt = theta
    if 2 < P < 4:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_decay(theta, x, y, yerr)

pos = soln_decay.x + 1e-4 * np.random.randn(32, 1)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, yerr))
sampler.run_mcmc(pos, 5000, progress=True)

'''
dp_dt Fit
'''
# Log likelihood for decay model
def log_likelihood_decay(theta, x, y, yerr):
    P, dp_dt = theta
    model_decay = (bjd - (t0 + P * x + 0.5 * P * dp_dt * x**2)) * 1440
    sigma2_decay = yerr ** 2 + model_decay ** 2
    return -0.5 * np.sum((y - model_decay) ** 2 / sigma2_decay + np.log(sigma2_decay))

np.random.seed(44)
nll_decay = lambda *args: -log_likelihood_decay(*args)
initial_decay = np.array([P_0, dp_dt_0]) + 0.1 * np.random.randn(2)
soln_decay = minimize(nll_decay, initial_decay, args=(x, y, yerr))
P_ml, dp_dt_ml = soln_decay.x

print("Maximum likelihood estimates:")
print(f'P = {P_ml} \ndp_dt = {dp_dt_ml}')

# Plotting the results
x0 = np.linspace(0, 1500, 1500)
t_tra_vary = t0 + P_ml * x0 + 0.5 * P_ml * dp_dt_ml * (x0 ** 2)
y_decay = (t_tra_con - t_tra_vary) * 1440 # Not sure if this is exactly how to go about converting to deviation from period?

plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
plt.plot(x0, y_decay, ":k", label="Decay")
plt.legend(fontsize=14)
plt.xlim(0, 1500)
plt.xlabel("x")
plt.ylabel("y")

# Markov Chains for decay model
def log_prior(theta):
    P, dp_dt = theta
    if 2 < P < 4 and 0 < dp_dt < 1:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_decay(theta, x, y, yerr)

pos = soln_decay.x + 1e-4 * np.random.randn(32, 2)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, yerr))
sampler.run_mcmc(pos, 5000, progress=True)

# Plotting the chains
fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["P", "dp_dt"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number")

flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)

# Corner plot
fig = corner.corner(flat_samples, labels=labels, truths=[P_0, dp_dt_0])