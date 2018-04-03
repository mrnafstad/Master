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
	rvir = False
	collapse = False
	for i in range(len(rdot)):
		
		if r[i] <= 0:
			collapse = True
		
		else:
			U = r[i]*r[i]/(Omega_m0/a[i]**3 + Omega_K/a[i]**2 + Lambda)*3./5.*(Omega_m0*(1+delta_i)/(2*r[i]**3) - Lambda) 
			T = rdot[i]*rdot[i]


		if T <= U:
			rvir = r[i]
			p = i

	rover = rvir/r[0]
	avir = a[p]
	
	"""
	if collapse :
		print "It collapses! delta_i = {:.5e}".format(delta_i)

	else:
		print "It does not collapse.. delta_i = {:.5e}".format(delta_i)
	"""
	if rvir:
		#print " %4.4e | %5.4e | %5.4e" % (rover, avir, delta_i)

		k = p
		print 
		while k + 1 <= len(r):
			r[k] = rvir
			k += 1
	elif runer==True:
		print "It does not collapse!", y[0]
	return r, collapse, y[p]


def r(x, y, Omega_m0, Omega_K, Lambda, r_i, delta_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (Omega_m0/a**3 + Lambda)**-1 * (r*(-Omega_m0*(1+delta_i)/(2.*r**3) + Lambda) + 3./2.*drdy*Omega_m0/a**3) 	

	return rr


def varyinitial():

	runer = True

	M = 200
	N = 500000

	Lambda = 0.74
	Omega_m0 = 0.26
	Omega_K = 1 - Omega_m0 - Lambda

	delta_i = 1.15e-3

	y0 = np.linspace(-8.5, -7, M)

	avir_arr = np.zeros(M)

	k = 0
	for i in range(len(y0)):
		yi = y0[i]
		k += 1
		if k == 10:
			percentage = float(i)/float(M)*100
			print "Still going strong at %.0f%%" % percentage
			k = 0
		
		r0 = np.exp(yi)
		drdx0 = np.exp(yi)
		y = np.linspace(yi, -1e-15, N)

		radius = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_i))

		rad, avir_arr[i] = virialcheck(y, radius, Omega_m0, Omega_K, Lambda, delta_i, runer)

	print len(avir_arr)

	mpl.plot(y0, avir_arr, "-.", linewidth = 0.75, label =  r"$a_{vir} (y_0)$")
	mpl.xlabel(r"$y_0$")
	mpl.ylabel(r"$a_{vir}$")
	mpl.legend()
	mpl.show()


def findcoll(delta_i_min, delta_i_max, delta_int, y0):
	
	N = 500000

	eps = 1e-13


	y = np.linspace(y0, -1e-15, N)


	#dummy numbers for tests
	delta_max = 1.0
	delta_min = 0.0

	delta_i = np.linspace(delta_i_min, delta_i_max, delta_int)

	r0 = np.exp(y0)
	drdx0 = np.exp(y0)
	x_coll = 0

	
	for i in range(len(delta_i)):
	
		sofa = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_i[i]))

		rad, collapse, x_coll1 = virialcheck(y, sofa, Omega_m0, Omega_K, Lambda, delta_i[i], runer)

		if (collapse and delta_i[i] < delta_max):
			delta_max = delta_i[i]
			delta_min = delta_i[i-1]
			x_coll = x_coll1
			print "---------------"*3
			print "Values stored, lets move on"
			break


		print "@ {:.1f} %".format(float(i)/float(delta_int)*100)

	if delta_max == 1.0:
		print "This aint working out.."

	return delta_min, delta_max, x_coll



Lambda = 0.74
Omega_m0 = 0.26

Omega_K = 1 - Omega_m0 - Lambda

runer = False

m = 50

y0_array = np.linspace(-11, -7.5, m)

dmax_array = np.zeros(m)

#print dmin, dmax, colltime

acceptance = float(sys.argv[5])
print acceptance


for i in range(len(y0_array)):

	dmin = float(sys.argv[1])
	dmax = float(sys.argv[2])
	delta_int = int(sys.argv[3])
	y0 = float(sys.argv[4])

	colltime = 10		#just a dummy for the while loop

	while abs(colltime) > acceptance:
		dmin, dmax, colltime = findcoll(dmin, dmax, delta_int, y0_array[i])
		diff = abs(abs(colltime) - acceptance)
		print "---------------"*3
		print "Still working towards collapse today! {:3.2f}".format(y0_array[i]) 
		print "{:.2e}".format(diff)
		print "---------------"*3

	dmax_array[i] = dmax

	print "A density between {:.6e} and {:.6} collapse at redshift {:.2e} (for the highest value)".format(dmin, dmax, np.exp(-colltime) - 1)




mpl.plot(y0_array, dmax_array, label = r"$\delta_i (y_0)$")
mpl.xlabel(r"$y_0$")
mpl.ylabel(r"$\delta_i$")
mpl.title("Overdensities collapsing today (or after z = {:.2f})".format(np.exp(-acceptance -1)))
mpl.legend()
mpl.show()

#varyinitial()