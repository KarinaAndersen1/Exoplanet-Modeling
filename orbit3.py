# Imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize, leastsq
import emcee
import corner
plt.ion()

'''
Suggested updates:
1. Label plots correctly
2. Name variables in a way that can be easily relatable back to theory
3. You get period and dP/dt from your fit. You should be able to construct
   an appropriate model from the best fit outputs.
4. Revisit your ML fit until you get something reasonable. Make sure it is
   reasonable with a graph of the output
5. Use the ML fit as the input for the MCMC
'''

# Read data file
data = pd.read_csv('C:\\Users\\karina\\Astronomy\\mid_transit_times4.csv')
epoch = data['Epoch'].values.astype('float')
y = data['Deviation_From_Period(min)'].values.astype('float')
bjd = data['BJD'].values.astype('float')
yerr = data['Error'].values.astype('float')

# Starting conditions
P_0 = 3.0395803
dp_dt_0 = -0.00000001
t0 = 2455186.01879 # Remains constant

'''
Line fit 
'''
# Log likelihood for a constant period
def log_likelihood_line(theta, x, y, yerr):
    P = theta
    model_line = (bjd - (t0 + P * epoch)) * 1440
    sigma2_line = yerr ** 2 + model_line ** 2
    return -0.5 * np.sum((y - model_line) ** 2 / sigma2_line + np.log(sigma2_line))

np.random.seed(49)
nll_line = lambda *args: -log_likelihood_line(*args)
initial_line = np.array([P_0]) + 0.1 * np.random.randn(1)
soln_line = minimize(nll_line, initial_line, args=(epoch, y, yerr))
P_line = soln_line.x

print("Maximum likelihood estimates:")
print('P = %.5f'%P_line)

# Plotting the results
x0 = np.linspace(0, 1500, 1500)
t_tra_con = t0 + P_line * x0
y_line = (t_tra_con - (t0 + P_line * x0)) * 1440 # Weird conversion to Deviation from Period

plt.figure(1)
plt.clf()
plt.errorbar(epoch, y, yerr=yerr, fmt=".k", capsize=0)
plt.plot(x0, y_line, ":k", label="Line")
plt.legend(fontsize=14)
plt.xlim(0, 1500)
# See suggestion #1 above
plt.xlabel("Epoch")
plt.ylabel("Deviation From Period (min)")

def log_prior(theta):
    P = theta
    # Changed prior to be what is physically expected
    if 2 < P < 4:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    return log_prior(theta) + log_likelihood_line(theta, epoch, y, yerr)

rand_array = np.random.randn(100, 1)
rand_array[:, 0] = rand_array[:, 0] * 1e-6
pos = soln_line.x + rand_array
#pos = np.array([P_line,-1e-6]) + 1e-6 * np.random.randn(100, 3)
nwalkers, ndim = np.shape(pos)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(epoch, y, yerr))
sampler.run_mcmc(pos, 5000, progress=True)

samples = sampler.get_chain()
plt.plot(samples[:, :, 0], "k", alpha=0.3)
plt.xlim(0, len(samples))
plt.ylabel('P')
plt.xlabel("step number")

flat_samples = sampler.get_chain(discard=1000, thin=10, flat=True)
print(f'Mean: {np.mean(flat_samples)}')
print(f'Std: {np.std(flat_samples)}')

'''
dp_dt Fit
'''
# Log likelihood for decay model
# See suggestion #2 above
def log_likelihood_decay(theta, x, y, yerr):
    P, dp_dt = theta
    model_decay = -(0.5 * P * dp_dt * (epoch**2)) * 1440
    sigma2_decay = yerr ** 2
    return -0.5 * np.sum((y - model_decay) ** 2 / sigma2_decay + np.log(2*np.pi*sigma2_decay))

np.random.seed(60)
nll_decay = lambda *args: -log_likelihood_decay(*args)
initial_decay = np.array([P_0, dp_dt_0]) + 0.1 * np.random.randn(2)
soln_decay = minimize(nll_decay, initial_decay, args=(epoch, y, yerr))
P_ml, dp_dt_ml = soln_decay.x

print("Maximum likelihood estimates:")
print('P = %.10f'%P_ml)
print('dP/dt = %.10f'%dp_dt_ml)

# Plotting the results
x0 = np.linspace(0, 1500, 1500)
y_decay = -(0.5 * P_ml * dp_dt_ml * (x0 ** 2)) * 1440

plt.figure(3)
plt.clf()
plt.errorbar(epoch, y, yerr=yerr, fmt=".k", capsize=0)
plt.plot(x0, y_decay, "r--", label="Decay")
plt.plot(x0, y_line, "b--", label="Constant Period")
plt.legend(fontsize=14)
plt.xlim(0, 1500)
plt.xlabel("Epoch")
plt.ylabel("Deviation From Period (min)")
plt.title('HAT-P-27 b')

# Markov Chains for decay model
def log_prior(theta):
    P, dp_dt = theta
    # Changed prior to be what is physically expected
    if 2 < P < 4 and -0.01 < dp_dt < 0.01:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    return log_prior(theta) + log_likelihood_decay(theta, epoch, y, yerr)

rand_array = np.random.randn(100, 2)
rand_array[:, 0] = rand_array[:, 0] * 1e-6
rand_array[:, 1] = rand_array[:, 1] * 1e-12
pos = soln_decay.x + rand_array
pos = soln_decay.x + 1e-6 * np.random.randn(100, 2)
#pos = np.array([P_line,-1e-6]) + 1e-6 * np.random.randn(100, 3)
nwalkers, ndim = np.shape(pos)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(epoch, y, yerr))
sampler.run_mcmc(pos, 5000, progress=True)

# Plotting the chains
discard = 1000
fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["P", "dP/dt"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    #ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number")

flat_samples = sampler.get_chain(discard=discard, thin=10, flat=True)

# Corner plot
fig = corner.corner(flat_samples, labels=labels, truths=[P_0, dp_dt_0])

sum_1 = 0
for i in y:
    pred1 = 0
    sum_1 += ((i - pred1)**2)

sum_2 = 0
for i in y:
    pred2 = (0.5 * P * dp * (i ** 2)) * 1440
    print(pred2)
    sum_2 += ((i - pred2)**2)



#define F-test function
def f_test(sum_1, sum_2, ddof1, ddof2):
    f = ((sum_1 - sum_2) / (ddof1 - ddof2)) / (sum_2 / ddof2)
    return f

#perform F-test
f_test(sum_1, sum_2, 1, 2)