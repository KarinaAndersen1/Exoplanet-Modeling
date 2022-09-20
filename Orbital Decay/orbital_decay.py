'''
Calculate the stellar Q and error given rate of change of the system's period
* Ensure that parameters are of the same corresponding units before running
'''

# Imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
import emcee

filepath = 'Users\karina\Astronomy'

# Parameters of star/planet system
m_star = np.random.normal(0.945, 0.035, 100000) 
m_planet = np.random.normal(0.00063, 0.0000315, 100000)
r_star = np.random.normal(0.898, 0.039, 100000)
r_planet = np.random.normal(0.106, 0.00362, 100000)
a = np.random.normal(8.658, 0.107, 100000)
e = np.random.normal(0.078, 0.047, 100000)
period = np.random.normal(3.006, 0.3935, 100000)
dp_dt = np.random.normal(-1e-3, 5e-11, 100000)

q_star = np.random.normal(1e8, 1e2, 100000)
# Calculation of Q, eq. from Patra 2020/Goldreich+Soter 1966
dp_dt = (((a / r_star) ** 5) * (m_planet / m_star) * (-27 * np.pi)) / (2 * q_star)

mean = np.round(np.mean(q_star))
std = np.round(np.std(q_star), decimals=3)
# Plot histogram of Q, obtain mean and std
fig = plt.figure(figsize=(16, 16))
ax = fig.add_subplot(111)

plt.hist(q_star, bins=200)
plt.xlabel('Q_star')
plt.ylabel('Frequency')
plt.title('Distribution of Q_star for a Given dp_dt')

mean = np.round(np.mean(q_star))
std = np.round(np.std(q_star), decimals=3)

ax.text(1000000, 1000, f'Mean: {mean} \nStd: {std}', family = 'Calibri', bbox=dict(facecolor='none', edgecolor='blue', pad=10.0), fontsize=20)
plt.show()

data = pd.read_csv('C:\\Users\\karina\\Astronomy\\mid_transit_times4.csv')
x = data['Epoch'].values
bjd = data['BJD'].values
yerr = data['Error'].values

y0 = 1
fit, cov = np.polyfit(x, y, 1, w=1/yerr, cov=True)
m, b = fit
p_err = np.srqt(cov[0,0])

plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
plt.plot(x0, np.dot(np.vander(x0, 2), [m, b]), "--k", label="Polyfit" )
plt.legend(fontsize=14)
plt.xlim(0, 1500)
plt.xlabel("x")
plt.ylabel("y")


P_0 = 3.0395803
dp_dt_0 = 0.0001
t0 = 2455186.01879
fit, cov = np.polyfit(x, y, 1, cov=True)

# Fit line with errors
# Input is period, output is slope, wrap in minimizer
y = data['Deviation_From_Period(min)']
fit2, cov2 = np.polyfit(x, y, 2, w=1/yerr, cov=True)

# f-test, MCMC on dp_dt, compare slope vs. error on slope for linear 
def log_likelihood(theta, x, y, yerr):
    P = theta
    model = (bjd - (t0 + P * x)) * 1440
    sigma2 = yerr ** 2 + model ** 2
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

np.random.seed(42)
nll = lambda *args: -log_likelihood(*args)
initial = np.array([P_0]) + 0.1 * np.random.randn(1)
soln = minimize(nll, initial, args=(x, y, yerr))
P_line = soln.x

print("Maximum likelihood estimates:")
print(f'P = {P_line}')

x0 = np.linspace(0, 1500, 1500)
t_tra = t0 + P_line * x0
y_line = (t_tra - (t0 + P_line * x0)) * 1440

plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
plt.plot(x0, y_line, ":k", label="Line")
plt.legend(fontsize=14)
plt.xlim(0, 1500)
plt.xlabel("x")
plt.ylabel("y")

def log_prior(theta):
    P = theta
    if 2 < P < 4:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

pos = soln.x + 1e-4 * np.random.randn(32, 1)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, yerr))
sampler.run_mcmc(pos, 5000, progress=True)

samples = sampler.get_chain()
plt.plot(samples[:,:,0], "k", alpha=0.5)
plt.xlim(0, len(samples))
plt.ylabel("P")
plt.show()

fig, axes = plt.subplots(1, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["P"]
#for i in range(ndim):
   # ax = axes[i]
   # ax.plot(samples[:, :, i], "k", alpha=0.3)
   # ax.set_xlim(0, len(samples))
   # ax.set_ylabel(labels[i])
   # ax.yaxis.set_label_coords(-0.1, 0.5)

#axes[-1].set_xlabel("step number")

flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
print(flat_samples.shape)



'''

def log_likelihood(theta, x, y, yerr):
    P, dp_dt = theta
    model = (bjd - (t0 + P * x + 0.5 * P * dp_dt * x**2)) * 1440
    sigma2 = yerr ** 2 + model ** 2
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

from scipy.optimize import minimize

np.random.seed(42)
nll = lambda *args: -log_likelihood(*args)
initial = np.array([P_0, dp_dt_0]) + 0.1 * np.random.randn(2)
soln = minimize(nll, initial, args=(x, y, yerr))
P_ml, dp_dt_ml = soln.x

print("Maximum likelihood estimates:")
print("P = {0:.3f}".format(P_ml))
print("dp_dt = {0:.3f}".format(dp_dt_ml))

x0 = np.linspace(0, 1500, 1500)
y_ml = (bjd - (t0 + P_ml * x + 0.5 * P_ml * dp_dt_ml * x**2)) * 1440

plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
plt.plot(x, y_ml, ":k", label="ML")
plt.legend(fontsize=14)
plt.xlim(0, 15000)
plt.xlabel("x")
plt.ylabel("y")
'''