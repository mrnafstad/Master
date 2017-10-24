import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint
import sys

"""
In this program we want to calculate the density contrast from the equations 
1) d'' + (3/a + E'/E)*d' - 4/3 * d'**2/(1+d) - 3/2 * omega/(a**5*E**2(a))*d(1+d) = 0, non linear
2) d'' + (3/a + E'/E)*d' - 3/2 * omega/(a**5*E**2(a))*d = 0
For d. We want to use the time factor ln(a), rather than a, and 
E = sqrt(omega_m,0/a**3 + Lambda/3) at least without radiation
"""

#We make the replacement x = ln(a)	

"""
with odeint
def d_nl():
	z = 0.0
	dz = 1e-5
	dzdt = - 3*np.exp(-x)*( 1 - Omega_m0/(2.*(Omega_m0*np.exp(-3*x) + Lambda/3.))) + 4/3.*ddelta**2/(1+d) + 4./3*Omega_m0/(Omega_m0*np.exp(-3*x) + Lambda/3)

def d_lin():

def Es(t, Lambda, Omega_m0):
	#t = variable x, keeping the symbol x open for use as array
	E = np.sqrt(Omega_m0*np.exp(-3*t) + Lambda/3.)
	dEoverE = -3./2 * Omega_m0*exp(-4*t)/(Omega_m0*np.exp(-3*t) + Lambda/3.)
	return [E,dEoverE]
"""
def function(x, t, Lambda, Omega_m0):
	#Es(t, Lambda, Omega_m0)
	#a = 3*(np.exp(-t) - Es[1]/2.)
	#b = 4/3.
	#c = 4/3. * Omega_m0*np.exp(-5*t)/Es[0]
	d = x[0]
	d1 = x[1]
	deltadd = [[],[]]
	deltadd[0] = d
	deltadd[1] = - 3*np.exp(-t)*(1 - Omega_m0/2./(Omega_m0*np.exp(-3*t) + Lambda/3.))*d1 + 4./3*d1**2/(1+d) + 4./3 * Omega_m0*np.exp(-5*t)/(np.sqrt(Omega_m0*np.exp(-3*t) + Lambda/3.))**2*d*(1+d)

	return deltadd

def nonlinear(x, a, Lambda, Omega_m0, Omega_K):
	#Es(t, Lambda, Omega_m0)
	#a = 3*(np.exp(-t) - Es[1]/2.)
	#b = 4/3.
	#c = 4/3. * Omega_m0*np.exp(-5*t)/Es[0]
	d = x[0]
	d1 = x[1]
	deltadd = [[],[]]
	deltadd[0] = d
	deltadd[1] = - 3./a*(1 - Omega_m0/2./a**3/(Omega_m0/a**3 - Omega_K/a**2 + Lambda/3.))*d1 + 4./3*d1**2/(1+d) + 4./3 * Omega_m0/a**5/(Omega_m0/a**3 - Omega_K/a**2 + Lambda/3.)*d*(1+d)

	return deltadd



def linear(x, a, Lambda, Omega_m0, Omega_K):

	d = x[0]
	d1 = x[1]
	deltadd = [[],[]]
	deltadd[0] = d
	deltadd[1] = - (3./a - 3/2. * Omega_m0/a**4/(Omega_m0/a**3 - Omega_K/a**2 + Lambda/3.))*d1 + 4./3 * Omega_m0/a**5/(Omega_m0/a**3 - Omega_K/a**2 + Lambda/3.)*d

	return deltadd

eps = 1e-4
N = 100000

a = np.linspace(eps, 1.1, N)
t = np.log(a)

Lambda = float(sys.argv[1])
Omega_m0 = float(sys.argv[2])
Omega_K = 1 - Omega_m0 - Lambda

delta_0 = 1.0
delta1_1 = 1e-6


#z2 = odeint(function, [0, 1e-5], t, args=(Lambda, Omega_m0)

deltanl = odeint(nonlinear, [delta_0, delta1_1], a, args=(Lambda, Omega_m0, Omega_K))

deltalin = odeint(linear, [delta_0, delta1_1], a, args=(Lambda, Omega_m0, Omega_K))

delta_c = np.array([1.686]*N)



mpl.plot(a, delta_c, "--r")

mpl.plot(a, deltanl[:, 0], "--b")
mpl.plot(a, deltalin[:, 0], "--g")
mpl.xlabel("a(t)")
mpl.ylabel(r" $ \delta $")
mpl.legend([r"$\delta _c$", r"$\delta _{non-linear}$", r"$\delta _{linear}$"])
mpl.xscale("log")
mpl.yscale("log")

mpl.show()

mpl.plot(a, delta_c, "--r")

mpl.plot(a, deltanl[:, 1], "--b")
mpl.plot(a, deltalin[:, 1], "--g")
mpl.xlabel("a(t)")
mpl.ylabel(r" $ \delta $")
mpl.legend([r"$\delta _c$", r"$\delta _{non-linear}$", r"$\delta _{linear}$"])
mpl.xscale("log")
mpl.yscale("log")

mpl.show()