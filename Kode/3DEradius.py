import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint
import sys
from itertools import cycle

H0 = 2.2685455*10**12	#km/s/Mpc

def r(x, y, Omega_m0, Omega_K, Lambda, r_i, delta_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (Omega_m0/a**3 + Omega_K/a**2 + Lambda)**(-2) * (r*(-Omega_m0*(1+delta_i)/r**3 + Lambda) + drdy/2.*a*(Omega_m0/a**3 - 2*Lambda)) 	

	return rr

eps = 1.63e-5
N = 10000
 
y0 = float(sys.argv[3])					#-2 might be a decent number, roughly -7 corresponds to recombination

y = np.linspace( y0, -1e-10, N)

Lambda = 0.74
Omega_m0 = 0.26

Omega_K = 1 - Omega_m0 - Lambda

delta_i_max = float(sys.argv[1])

delta_int = int(sys.argv[2])

delta_i = np.linspace(0, delta_i_max, delta_int)

r0 = np.exp(y0)
drdx0 = np.exp(y0)

"""
radius = odeint(r, [r0, drdx0], y, args=(Omega_m0, Omega_K, Lambda, r0, delta_i[0]), mxstep=1000000000)#	, full_output=1)

mpl.plot(y, radius[:,1], "--c", linewidth = 0.75, label = r"Background ($\delta_i = 0$)") 		#TypeError: tuple indices must be integers, not tuple
mpl.xlabel("x = ln(a(t))")
mpl.ylabel(r"$s = \frac{r}{r_i}$")
mpl.legend()
#mpl.yscale("log")
#mpl.xscale("log")
"""

if len(sys.argv) == 5:
	#Run through different values of the initial overdensity and plot together
	for i in range(len(delta_i)):

		sofa = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_i[i]))

		if i == 0:

			mpl.plot(y, sofa[:,0], "-c", linewidth = 0.75, label = r"$\delta_{i} =$ %.1e" % delta_i[i])

		else:

			mpl.plot(y, sofa[:,0], "-.", linewidth = 0.75, label = r"$\delta_{i} =$ %.1e" % delta_i[i])

		print delta_i[i]




	mpl.xlabel("x = ln(a)")
	mpl.ylabel("r(a)")
	mpl.legend( loc=2)

	mpl.title(r"$x_0 =$ %0.2f, varying values of $\delta_i$" % y0)

	mpl.show()

else:
#Want to run through different values of y0 at



	y0_arr = np.linspace(-8, -1, 10)



	for i in range(len(y0_arr)):


		delta__i = 1.686

		y = np.linspace( y0_arr[i], -1e-10, N)
		sofa = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta__i))	

		mpl.plot(y, sofa[:,0], "--", linewidth = 0.75, label = r"$x_0 =$ %.1f, $\delta_i =$ %.2f" % (y0_arr[i], delta__i))



		delta__i = 0.0

		y = np.linspace( y0_arr[i], -1e-10, N)
		sofa = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta__i))

		mpl.plot(y, sofa[:,0], "-c", linewidth = 0.75)#, label = r"$x_0 =$ %.1f, $\delta_i =$ %.2f" % (y0_arr[i], delta__i))



		delta__i = 1.0

		y = np.linspace( y0_arr[i], -1e-10, N)
		sofa = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta__i))

		mpl.plot(y, sofa[:,0], ":", linewidth = 0.75, label = r"$x_0 =$ %.1f, $\delta_i =$ %.2f" % (y0_arr[i], delta__i))

		delta__i = 0.5

		y = np.linspace( y0_arr[i], -1e-10, N)
		sofa = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta__i))

		mpl.plot(y, sofa[:,0], "-.", linewidth = 0.75, label = r"$x_0 =$ %.1f, $\delta_i =$ %.2f" % (y0_arr[i], delta__i))


	mpl.xlabel("x = ln(a)")
	mpl.ylabel("r(a)")
	mpl.legend(loc=2)
	mpl.title(r"$\delta_i =$ %.2f, varying values of $x_0$" % delta__i)
	mpl.show()

"""
y0_arr = np.linspace(-8, -1, 10)

lines = ["--", "-."]



for i in range(len(delta_i)):

	for j in range(len(y0_arr)):

		for k in range(len(lines)):

			y = np.linspace( y0_arr[j], -1e-10, N)

			sofa = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_i[i]))

			mpl.plot(y, sofa[:,0], lines[k], linewidth = 0.75, label = r"$\delta_{i} =$ %.1f" % delta_i[i])


mpl.xlabel("x = ln(a)")
mpl.ylabel("r(a)")
mpl.legend()

mpl.title(r"$x_0 =$ %0.2f, varying values of $\delta_i$" % y0)

mpl.show()
"""