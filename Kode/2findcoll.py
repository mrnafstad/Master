import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint
import sys
from itertools import cycle

H0 = 2.2685455*10**12	#km/s/Mpc

def virialcheck(y, mix, Omega_m0, Omega_K, Lambda, delta_i, runer):
	rdot = mix[:,1]
	r = mix[:,0]

	a = np.exp(y)
	p = 0
	t = 0
	rvir = False
	collapse = False

	first = True
	values = [0, 0, 0, 0, collapse]

	for i in range(len(rdot)):
		
		if ( r[i] <= 0 and first ):
			collapse = True
			first = False
			t = i

		
		
		U = r[i]*r[i]/(Omega_m0/a[i]**3 + Omega_K/a[i]**2 + Lambda)*3./5.*(Omega_m0*(1+delta_i)/(2*r[i]**3) - Lambda) 
		T = rdot[i]*rdot[i]
		

		

		if T <= U:
			rvir = r[i]
			p = i


	rover = rvir/r[0]
	avir = a[p]
	
	if p >= 1:
		rvir = r[p]
	rover = rvir/r[0]
	avir = a[p]

	if rvir:
		odensity = (Omega_m0*(1+delta_i)/rvir**3 + Lambda)/(Omega_m0/avir**3 + Lambda)
		s = np.argmax(r)
		rmax = r[s]
		amax = a[s]

		odensitymax = (Omega_m0*(1+delta_i)/rmax**3 + Lambda)/(Omega_m0/amax**3 + Lambda)

		values = [odensitymax, odensity, rvir/rmax, avir/amax, collapse]

	if rvir:

		k = p

		while k + 1 <= len(r):
			r[k] = rvir
			k += 1

	return r, collapse, y[t], values


def r(x, y, Omega_m0, Omega_K, Lambda, r_i, delta_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (Omega_m0/a**3 + Lambda)**-1 * (r*(-Omega_m0*(1+delta_i)/(2.*r**3) + Lambda) + 3./2.*drdy*Omega_m0/a**3) 	

	return rr


def findcoll(delta_max, delta_min, y0):
	
	N = 500000

	eps = 1e-13

	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	outputvals = np.zeros(5)


	for i in range(15):

		delta_mid = (delta_max + delta_min)/2.0
		#print "dmin: {:.15e}, dmax: {:.15e}".format(delta_min, delta_max)

		#for deltamax
		radiusmax = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_max))

		radmax, collmax, colltime_max, currentvalue = virialcheck(y, radiusmax, Omega_m0, Omega_K, Lambda, delta_max, runer)

		#for deltamin
		radiusmin = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_min))

		radmin, collmin, colltime_min, currentmin = virialcheck(y, radiusmin, Omega_m0, Omega_K, Lambda, delta_min, runer)

		#for deltamid
		radiusmid = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_mid))

		radmid, collmid, colltime_mid, currentmid = virialcheck(y, radiusmid, Omega_m0, Omega_K, Lambda, delta_mid, runer)

		if ( collmax and collmid ):
			delta_max = delta_mid
			outputvals = currentvalue

		elif ( collmax and not collmid ):
			delta_min = delta_mid

		x_coll = colltime_max


	return delta_min, delta_max, x_coll, outputvals

file = open("fittingvalues.txt", "w")

file.write("      z0          |          ODMax         |          ODVir           |        rvir/rmax         |      avir/amax       |     bool collapse     |      collapse time\n")

Lambda = 0.74
Omega_m0 = 0.26

Omega_K = 1 - Omega_m0 - Lambda

runer = False

m = 50

y0_array = np.linspace(-11, -6.9, m)

dmax_array = np.zeros(m)

#print dmin, dmax, colltime

acceptance = float(sys.argv[1])
print acceptance


for i in range(len(y0_array)):

	dmin = 0.0
	dmax = 0.01

	colltime = 10		#just a dummy for the while loop

	while abs(colltime) > acceptance:
		print "---------------"*3
		print "Still working towards collapse today! {:3.2f}".format(y0_array[i]) 

		dmin, dmax, colltime, vals = findcoll(dmax, dmin, y0_array[i])
		diff = abs(abs(colltime) - acceptance)
		print "{:.2e}".format(diff)
		print "---------------"*3

	dmax_array[i] = dmax

	file.write("   {:<10.2f}     |       {:5.10f}     |     {:5.10e}     |     {:5.10e}     |     {:5.10f}     |         {}          |     {:5.10e}     |     {:5.10e}\n".format(np.exp(-y0_array[i]) -1, vals[0], vals[1], vals[2], vals[3], vals[4], np.exp(-colltime) - 1, dmax_array[i]))


	#print "A density between {:.6e} and {:.6e} collapse at redshift {:.2e} (for the highest value)".format(dmin, dmax, np.exp(-colltime) - 1)

file.close()

mpl.plot(y0_array, dmax_array, label = r"$\delta_i (y_0)$")
mpl.xlabel(r"$y_0$")
mpl.ylabel(r"$\delta_i$")
mpl.title("Overdensities collapsing today (or after z = {:.2f})".format(np.exp(-acceptance -1)))
mpl.legend()
mpl.show()

#varyinitial()