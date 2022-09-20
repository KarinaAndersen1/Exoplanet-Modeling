import numpy as np
import matplotlib.pyplot as plt

#q_star = np.random.normal(100, 5, 100000)

#k2 = np.random.normal(0.5, 0.01, 100000) # N at least 100,000

m_star = np.random.normal(0.945, 0.035, 100000) 

m_planet = np.random.normal(0.00063, 0.0000315, 100000)

r_star = np.random.normal(0.898, 0.039, 100000)

r_planet = np.random.normal(0.106, 0.00362, 100000)

a_cur = np.random.normal(8.658, 0.107, 100000)

e_cur = np.random.normal(0.078, 0.047, 100000)

period = np.random.normal(262569.6, 0.000001, 100000)


m_star = 0.945

m_planet = 0.00063

r_star = 0.898

r_planet = 0.106

a_cur = 8.658

e_cur = 0.078

period = 262569.6

G = 6.674e-11

Q_planet = 10 ** (6.5)

q_star = 10 ** (5.5)

#dp_dt = np.random.normal(-.00325, .00162, 100000)

#q_star = (((a / r_star) ** 5) * (m_planet / m_star) * (-27 * np.pi)) / (2 * dp_dt)

#dp_dt = ((-9 * np.pi * k2) / q_star) * (m_star / m_planet) * ((a / r_star) ** 5)

e = e_cur
a = a_cur

ecc_const = -((((63/4)*(G*(m_star**3))**0.5)*((r_planet**5)/(Q_planet*m_planet))+(171/16)*((G/m_star)**0.5)*((m_planet*(r_star)**5)/q_star)))
a_const_1 = ((63/2)*(G*(m_star**3))**0.5)*((r_planet**5)/(Q_planet*m_planet))
a_const_2 = (9/2)*((G/m_star)**0.5)*((m_planet*(r_star)**5)/q_star)

def evolution(e, a):
    a = axis(a)
    e = ecc(e)
    return np.mean(e)

def axis(a):
    a = a + (-a*(((63/2)*(G*(m_star**3))**0.5)*((r_planet**5)/(Q_planet*m_planet))*(e**2)+(9/2)*((G/m_star)**0.5)*((m_planet*(r_star)**5)/q_star))*(a**(-13/2)))
    return np.mean(a)

def ecc(e):

    e = e + (-e*(((63/4)*(G*(m_star**3))**0.5)*((r_planet**5)/(Q_planet*m_planet))+(171/16)*((G/m_star)**0.5)*((m_planet*(r_star)**5)/q_star))*(a**(-13/2)))
    return np.mean(e)

def period(period):
    period = period + (-(((7)*(G*(m_star**3))**0.5)*((r_planet**5)/(Q_planet*m_planet))*(e**2)+((G/m_star)**0.5)*((m_planet*(r_star)**5)/q_star))*((27*np.pi*(a**(-7/2)))/(period*m_star)))
    return np.mean(period)
'''
- model of a over time
- model of e over time
- model of P over time
- model of e vs. P
'''
from scipy.integrate import odeint

def f(y, t):
    Ei = y[0]
    Ai = y[1]

    f0 = ecc_const * (Ai ** (-13/2)) * Ei
    f1 = -1 * [a_const_1 * (Ei ** 2) + a_const_2] * (Ai ** (-11/2))
    
    return [f0, f1]

E0 = e_cur
A0 = a_cur
y0 = [E0, A0]
t = np.linspace(0, -15000000000, 15)

# solve the DEs
soln = odeint(f, y0, t)
E = soln[:, 0]
A = soln[:, 1]

plt.figure()
plt.plot(t, E, label='Ecc')
plt.plot(t, A, label='Axis')
plt.xlabel('Time')
plt.ylabel('Scale')
plt.title('Orbital Decay')
#plt.legend(loc=0)
plt.show()
