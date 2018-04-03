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


def expansionfunc(a, Omega_m0, Omega_K, Lambda):

	E = np.sqrt(Omega_m0/a**3 + Omega_K/a**2 + Lambda/3.)

	dE_E = - (3./2. *Omega_m0/a**4 + Omega_K/a**3)/(Omega_m0/a**3 + Omega_K/a**2 + Lambda)

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

	# Skriv om til d/dx !! skriv om dEda til dE/dx

	a = np.exp(y)

	Es = expansionfunc(a, Omega_m0, Omega_K, Lambda)
	
	E = Es[0]
	dEdtE = Es[1]

	d = x[0]
	d1 = x[1]
	deltadd = [[],[]]
	deltadd[0] = d1
	deltadd[1] = - (2. + a*dEdtE)*d1 + 4./3.*d1**2/(1+d) + 4./3. * float(Omega_m0)/a**3/E**2*d*(1+d)

	return deltadd



def linear(x, y, Lambda, Omega_m0, Omega_K):

	a = np.exp(y)

	Es = expansionfunc(a, Omega_m0, Omega_K, Lambda)
	E = Es[0]
	dEdtE = Es[1]

	d = x[0]
	d1 = x[1]
	deltadd = [[],[]]
	deltadd[0] = d1
	deltadd[1] = - (2 + a*dEdtE)*d1 + 4./3 * float(Omega_m0)/a**3/E**2*d

	return deltadd

eps = float(sys.argv[4])
N = 100000

y = np.linspace(eps, 1e-15, N)


Lambda = 0
Omega_m0 = 1

Omega_K = 1



delta_0 = np.linspace( float(sys.argv[2]), float(sys.argv[1]), int(sys.argv[3]) )
delta1_1 = 5e-5


#z2 = odeint(function, [0, 1e-5], t, args=(Lambda, Omega_m0)
"""
deltanl = odeint(nonlinear, [delta_0[0], delta1_1], y, args=(Lambda, Omega_m0, Omega_K))

deltalin = odeint(linear, [delta_0[0], delta1_1], y, args=(Lambda, Omega_m0, Omega_K))

print deltanl[-1,0], deltalin[-1,0]

delta_c = np.array([1.686]*N)



mpl.plot(y, delta_c, "--c", linewidth = 0.75, label=r"$\delta _c = 1.686$")

mpl.plot(y, deltanl[:, 0], "--b", linewidth = 0.75, label = r"$\delta _{non-linear}$")
mpl.plot(y, deltalin[:, 0], "-y", linewidth = 0.75, label = r"$\delta _{linear}$")
mpl.xlabel("x = ln(a)")
mpl.ylabel(r" $ \delta $")
mpl.legend()#[r"$\delta _c = 1.686$", r"$\delta _{non-linear}$", r"$\delta _{linear}$"])


mpl.show()
"""
for i in range(len(delta_0)):

	deltanl = odeint(nonlinear, [delta_0[i], delta1_1], y, args=(Lambda, Omega_m0, Omega_K))

	mpl.plot(y, deltanl[:, 0], "--", linewidth = 0.75, label = r"$\delta _{nl, i} =$ %.2e" % delta_0[i])

for i in range(len(delta_0)):

	deltalin = odeint(linear, [delta_0[i], delta1_1], y, args=(Lambda, Omega_m0, Omega_K))

	mpl.plot(y, deltalin[:, 0], "-", linewidth = 0.75, label = r"$\delta _{i} =$ %.2e" % delta_0[i])

delta_c = np.array([1.686]*N)

mpl.plot(y, delta_c, "--c", linewidth = 0.75, label=r"$\delta _c = 1.686$")

mpl.xlabel("x = ln(a)")
mpl.ylabel(r"$\delta$")
mpl.legend( loc = 2 )
mpl.title(r"$\delta (x)$, $x_0 =$ %.2f" % eps)
mpl.grid( True )
mpl.show()