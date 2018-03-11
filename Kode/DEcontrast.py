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


"""

H_0 = 1.	#km/s/Mpc

def expansionfunc(a, Omega_m0, Omega_K, Lambda):

	E = np.sqrt(Omega_m0/a**3 + Omega_K/a**2 + Lambda/3.)

	dE_E = - (3./2. *Omega_m0/a**4 + Omega_K/a**3)/(Omega_m0/a**3 + Omega_K/a**2 + Lambda/3.)

	return [E, dE_E]

	

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
"""


def nonlinear(x, a, Lambda, Omega_m0, Omega_K):
	#Es(t, Lambda, Omega_m0)
	#a = 3*(np.exp(-t) - Es[1]/2.)
	#b = 4/3.
	#c = 4/3. * Omega_m0*np.exp(-5*t)/Es[0]

	Es = expansionfunc(a, Omega_m0, Omega_K, Lambda)
	E = Es[0]
	dEdtE = Es[1]

	d = x[0]
	d1 = x[1]
	deltadd = [[],[]]
	deltadd[0] = d1
	deltadd[1] = - (3./a + dEdtE)*d1 + 4./3.*d1**2/(1+d) + 4./3. * float(Omega_m0)/a**5/E**2*d*(1+d)

	return deltadd



def linear(x, a, Lambda, Omega_m0, Omega_K):

	Es = expansionfunc(a, Omega_m0, Omega_K, Lambda)
	E = Es[0]
	dEdtE = Es[1]

	d = x[0]
	d1 = x[1]
	deltadd = [[],[]]
	deltadd[0] = d1
	deltadd[1] = - (3./a + dEdtE)*d1 + 4./3 * float(Omega_m0)/a**5/E**2*d

	return deltadd

eps = 1.63e-5
N = 100000

a = np.linspace(eps, 1.1, N)
t = np.log(a)

Lambda = float(sys.argv[1])
Omega_m0 = float(sys.argv[2])

if len(sys.argv) == 3:
	Omega_K = 1 - Omega_m0 - Lambda
else:
	Omega_K = float(sys.argv[3])

print Omega_K

delta_0 = 1e-4
delta1_1 = 5e-5


#z2 = odeint(function, [0, 1e-5], t, args=(Lambda, Omega_m0)

deltanl = odeint(nonlinear, [delta_0, delta1_1], a, args=(Lambda, Omega_m0, Omega_K))

deltalin = odeint(linear, [delta_0, delta1_1], a, args=(Lambda, Omega_m0, Omega_K))

print deltanl[-1,0], deltalin[-1,0]

delta_c = np.array([1.686]*N)



mpl.plot(a, delta_c, "--c", linewidth = 0.75, label=r"$\delta _c = 1.686$")

mpl.plot(a, deltanl[:, 0], "--b", linewidth = 0.75, label = r"$\delta _{non-linear}$")
mpl.plot(a, deltalin[:, 0], "-y", linewidth = 0.75, label = r"$\delta _{linear}$")
mpl.xlabel("a(t)")
mpl.ylabel(r" $ \delta $")
mpl.legend()#[r"$\delta _c = 1.686$", r"$\delta _{non-linear}$", r"$\delta _{linear}$"])
mpl.xscale("log")
mpl.yscale("log")

mpl.show()

"""
mpl.plot(a, delta_c, "--c", linewidth = 0.75)

mpl.plot(a, deltanl[:, 1], "--b", linewidth = 0.75)
mpl.plot(a, deltalin[:, 1], "--g", linewidth = 0.75)
mpl.xlabel("a(t)")
mpl.ylabel(r" $ \delta $")
mpl.legend([r"$\delta _c$", r"$\delta _{non-linear}$", r"$\delta _{linear}$"])
mpl.xscale("log")
mpl.yscale("log")

mpl.show()
"""

#Rewrite the following in terms of d/dx, rather than d/da!

#Er noe bug her!


def expansionfunc(a, Omega_m0, Omega_K, Lambda):

	E = np.sqrt(Omega_m0/a**3 + Omega_K/a**2 + Lambda/3.)

	dE_E = - (3./2. *Omega_m0/a**4 + Omega_K/a**3)/(Omega_m0/a**3 + Omega_K/a**2 + Lambda/3.)

	return [E, dE_E]

	

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
"""


def nonlinear(x, y, Lambda, Omega_m0, Omega_K):
	#Es(t, Lambda, Omega_m0)
	#a = 3*(np.exp(-t) - Es[1]/2.)
	#b = 4/3.
	#c = 4/3. * Omega_m0*np.exp(-5*t)/Es[0]

	Es = expansionfunc(x, Omega_m0, Omega_K, Lambda)
	E = Es[0]
	dEdtE = Es[1]

	a = np.exp(y)

	d = x[0]
	d1 = x[1]
	deltadd = [[],[]]
	deltadd[0] = d1
	deltadd[1] = - (2 + dEdtE)*d1 + 4./3.*d1**2/(1+d) + 4./3. * float(Omega_m0)/a**3/E**2*d*(1+d)

	return deltadd



def linear(x, a, Lambda, Omega_m0, Omega_K):

	Es = expansionfunc(a, Omega_m0, Omega_K, Lambda)
	E = Es[0]
	dEdtE = Es[1]

	d = x[0]
	d1 = x[1]
	deltadd = [[],[]]
	deltadd[0] = d1
	deltadd[1] = - (3./a + dEdtE)*d1 + 4./3 * float(Omega_m0)/a**5/E**2*d

	return deltadd

eps = 1.63e-5
N = 100000

x = np.linspace(-7, 0, N)


Lambda = float(sys.argv[1])
Omega_m0 = float(sys.argv[2])

if len(sys.argv) == 3:
	Omega_K = 1 - Omega_m0 - Lambda
else:
	Omega_K = float(sys.argv[3])

print Omega_K

delta_0 = 1e-4
delta1_1 = 5e-5


#z2 = odeint(function, [0, 1e-5], t, args=(Lambda, Omega_m0)

deltanl2 = odeint(nonlinear, [delta_0, delta1_1], x, args=(Lambda, Omega_m0, Omega_K))

deltalin2 = odeint(linear, [delta_0, delta1_1], x, args=(Lambda, Omega_m0, Omega_K))

print deltanl2[-1,0], deltalin2[-1,0]

delta_c = np.array([1.686]*N)



mpl.plot(x, delta_c, "--c", linewidth = 0.75, label=r"$\delta _c = 1.686$")

mpl.plot(x, deltanl2[:, 0], "--b", linewidth = 0.75, label = r"$\delta _{non-linear}$")
mpl.plot(x, deltalin2[:, 0], "-y", linewidth = 0.75, label = r"$\delta _{linear}$")
mpl.xlabel("a(t)")
mpl.ylabel(r" $ \delta $")
mpl.legend()#[r"$\delta _c = 1.686$", r"$\delta _{non-linear}$", r"$\delta _{linear}$"])
mpl.xscale("log")
mpl.yscale("log")

mpl.show()