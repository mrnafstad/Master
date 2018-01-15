import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint
import sys

def s(x, z, omega_m, omega_k, Lambda):
	s = x[0]
	dsdz = np.sqrt(omega_m*(1 + z)**3 + omega_k*(1 + z)**2 + Lambda)/np.sqrt(omega_m*(1 + z)*(1./s) + omega_k + Lambda*s**2/(1 + z)**2)

	return dsdz


def s_a(x, a, omega_m, omega_k, Lambda, a0):
	s = x[0]
	dsda = np.sqrt(omega_m/a0**3/s + omega_k/a0**2 + Lambda/3.*s**2)/np.sqrt(a**2*(omega_m/a**3 + omega_k/a**2 + Lambda/3.))

	return dsda

N = 10000
a0 = 1.65e-5

z = np.linspace(20, 0, N)

a = np.linspace(a0, 1, N)

Lambda = float(sys.argv[1])

Omega_m0_min = 1.0
Omega_m0_max = float(sys.argv[2])
n=11

Omega_m0 = np.linspace(Omega_m0_min, Omega_m0_max, n)

if len(sys.argv) == 3:
	Omega_K = 1 - Omega_m0 - Lambda
else:
	Omega_K = np.zeros(n)
	for j in range(n):
		Omega_K[j] = float(sys.argv[3])

omega = Omega_K + Omega_m0 + Lambda

s0 = 0.5

#result = odeint(s, s0, z, args=(Lambda, Omega_m0, Omega_K))


for i in range(len(Omega_m0)):
	sofa = odeint(s_a, s0, a, args = (Omega_m0[i], Omega_K[i], Lambda, a0))

	mpl.plot(a, sofa, linewidth = 0.75, label = r"s, $\Omega_{m0} =$ %.1f" % Omega_m0[i])
	mpl.xlabel("a(t)")
	mpl.ylabel("s(a)")
	mpl.legend()
	mpl.yscale("log")
	mpl.xscale("log")


mpl.show()