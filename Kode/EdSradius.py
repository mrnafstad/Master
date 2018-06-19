import numpy as np
import matplotlib.pylab as mpl
import sys
from scipy.integrate import odeint

def general_virial(y, mix, delta_i, Kinetic, Potential, Overdensity, EdS):
	#beginning with the skin

	rdot = mix[:,1]
	r = mix[:,0]

	a = np.exp(y)
	s = 0
	p = 0
	rvir = False
	collapse = False
	k = 0

	s = np.argmax(r)
	rmax = r[s]
	amax = a[s]

	U = np.zeros(len(r))

	K = np.zeros(len(r))

	while s <= len(r)-1:

		if r[s] <= 0:
			r[s] = 0
			collapse = True

		K[s] = Kinetic(r[s], rdot[s], a[s], delta_i)
		U[s] = Potential(r[s], rdot[s], a[s], delta_i, EdS)

		if K[s] <= 4*U[s]:
			p = s


		s += 1

	if p > 1:
		rvir = r[p]
		avir = a[p]

		Overdensity(rmax, amax, rvir, avir, collapse, delta_i, EdS)

		while p + 1 <= len(r):
			r[p] = rvir
			p += 1

	return r, avir, U, K
"""
def potential():


def kinetic():


def overdensity():
"""

def EdSr(x, y, r_i, delta_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = -(1 + delta_i)*a**3/(2.0*r**2) + 3.0/2.0*drdy	

	return rr

N = 5e5

y0 = float(sys.argv[1])
y = np.linspace(y0, 1e-15, N)

file = open("Edsvalues", "w")

r0 = drdx0 = np.exp(y0)

delta_i = float(sys.argv[2])

radius = odeint(EdSr, [r0, drdx0], y, args = (0, delta_i))

mpl.plot(y, radius[:,0])
mpl.show()