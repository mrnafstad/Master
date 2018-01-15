import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint
import sys

H0 = 2.2685455*10**12	#km/s/Mpc

def r(x, a, Omega_m0, Omega_K, Lambda):
	r = x[0]
	drda = x[1]

	rr = [[],[]]
	rr[0] = drda
	rr[1] = r/a**2*(-3.*Omega_m0/(8*np.pi) + Lambda/3. / (H0**2*(Omega_m0/a**3 + Omega_K/a**2 + Lambda/3.))) - \
	1./(Omega_m0/a**3 + Omega_K/a**2 + Lambda/3.)*(Lambda/3. - 2*Omega_m0/a**3)*drda

	return rr

eps = 1.63e-5
N = 10000

a = np.linspace(eps, 1.1, N)

Lambda = float(sys.argv[2])
Omega_m0 = float(sys.argv[3])

if len(sys.argv) == 4:
	Omega_K = 1 - Omega_m0 - Lambda
else:
	Omega_K = float(sys.argv[4])

print Omega_K

r0 = 1.0
drda0 = float(sys.argv[1])

radius = odeint(r, [r0, drda0], a, args=(Omega_m0, Omega_K, Lambda))

mpl.plot(a, radius[:, 0], "--c", linewidth = 0.75, label = "s(a)")
mpl.xlabel("a(t)")
mpl.ylabel(r"$s = \frac{r}{r_i}$")
mpl.legend()
mpl.yscale("log")
mpl.xscale("log")

mpl.show()