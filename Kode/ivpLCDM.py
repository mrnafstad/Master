import matplotlib.pylab as mpl
import numpy as np
import sys
from scipy.integrate import solve_ivp


def r(y, x):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda) 	

	return rr

def kollaps(y, x):
	return x[1]



kollaps.terminal = True

Omega_m0 = 0.26
Lambda = 0.74

delta_i = float(sys.argv[1])

y0 = np.log(1/(1 + float(sys.argv[2])))
N = 50000
y = np.linspace(y0, 1e-12, N)

r0 = np.exp(y0)
drdx0 = np.exp(y0)

radius = solve_ivp(r, [y0, 1e-12], [r0, drdx0], method = 'RK45', t_eval = y, events = kollaps)
print type(radius.t_events[0][0])
print "Collapse at z ={:2.5f}".format(np.exp(-radius.t_events[0][0]) - 1)

mpl.plot(radius.t, radius.y[1])
mpl.show()