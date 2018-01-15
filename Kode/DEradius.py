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
a0 = 1.63e-5

z = np.linspace(20, 0, N)

a = np.linspace(a0, 1, N)

Lambda = float(sys.argv[1])
Omega_m0 = float(sys.argv[2])

if len(sys.argv) == 3:
	Omega_K = 1 - Omega_m0 - Lambda
else:
	Omega_K = float(sys.argv[3])

omega = Omega_K + Omega_m0 + Lambda

s0 = 1.0

#result = odeint(s, s0, z, args=(Lambda, Omega_m0, Omega_K))

sofa = odeint(s_a, s0, a, args = (Omega_m0, Omega_K, Lambda, a0))

"""
mpl.plot(z, result, linewidth = 0.75, label = "s")
mpl.xlabel("z(t)")
mpl.ylabel("s(z)")
mpl.legend()
#mpl.yscale("log")
mpl.xscale("log")
mpl.show()
"""

mpl.plot(a, sofa, linewidth = 0.75, label = r"s, $\Omega =$ (%.1f, %.1f, %.1f) = %.1f " %(Omega_m0, Omega_K, Lambda, omega))
mpl.xlabel("a(t)")
mpl.ylabel("s(a)")
mpl.title(r"$s_0 = 1.0 $")
mpl.legend()
mpl.yscale("log")
mpl.xscale("log")
mpl.show()

s0 = 0.5

#result = odeint(s, s0, z, args=(Lambda, Omega_m0, Omega_K))

sofa = odeint(s_a, s0, a, args = (Omega_m0, Omega_K, Lambda, a0))

"""
mpl.plot(z, result, linewidth = 0.75, label = "s")
mpl.xlabel("z(t)")
mpl.ylabel("s(z)")
mpl.legend()
#mpl.yscale("log")
mpl.xscale("log")
mpl.show()
"""

mpl.plot(a, sofa, linewidth = 0.75, label = r"s, $\Omega =$ (%.1f, %.1f, %.1f) = %.1f " %(Omega_m0, Omega_K, Lambda, omega))
mpl.xlabel("a(t)")
mpl.ylabel("s(a)")
mpl.title(r"$s_0=0.5$")
mpl.legend()
mpl.yscale("log")
mpl.xscale("log")
mpl.show()

s0 = 2.0

#result = odeint(s, s0, z, args=(Lambda, Omega_m0, Omega_K))

sofa = odeint(s_a, s0, a, args = (Omega_m0, Omega_K, Lambda, a0))

"""
mpl.plot(z, result, linewidth = 0.75, label = "s")
mpl.xlabel("z(t)")
mpl.ylabel("s(z)")
mpl.legend()
#mpl.yscale("log")
mpl.xscale("log")
mpl.show()
"""

mpl.plot(a, sofa, linewidth = 0.75, label = r"s, $\Omega =$ (%.1f, %.1f, %.1f) = %.1f " %(Omega_m0, Omega_K, Lambda, omega))
mpl.xlabel("a(t)")
mpl.ylabel("s(a)")
mpl.title(r"$s_0=2.0$")
mpl.legend()
mpl.yscale("log")
mpl.xscale("log")
mpl.show()

