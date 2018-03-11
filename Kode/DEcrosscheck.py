import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint
import sys
from itertools import cycle



def virialcheck(y, mix, Omega_m0, Omega_K, Lambda, delta_i, runer):
	rdot = mix[:,1]
	r = mix[:,0]

	a = np.exp(y)
	p = 0
	rvir = False
	collapse = False
	this = True
	for i in range(len(rdot)):

		U = r[i]*r[i]/(Omega_m0/a[i]**3 + Omega_K/a[i]**2 + Lambda)*3./5.*(Omega_m0*(1+delta_i)/(2*r[i]**3) - Lambda) 
		T = rdot[i]*rdot[i]

		if r[i] <= 0:
			collapse = True

		if T <= U:
			rvir = r[i]
			p = i

	rover = rvir/r[0]
	avir = a[p]
	
	"""
	if collapse :
		print "It collapses!"

	else:
		print "It does not collapse.."
	"""
	if rvir:
		#print " %4.4e | %5.4e | %5.4e" % (rover, avir, delta_i)

		k = p
		while k + 1 <= len(r):
			r[k] = rvir
			k += 1
	"""
	if collapse:
		print "It does not collapse!", y[0]
		this = False
	"""
	return r, avir, collapse


def r(x, y, Omega_m0, Omega_K, Lambda, r_i, delta_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (Omega_m0/a**3 + Lambda)**-1 * (r*(-Omega_m0*(1+delta_i)/(2.*r**3) + Lambda) + 3./2.*drdy*Omega_m0/a**3) 	

	return rr


def crosschecker(y0max, y0end, delta_i_max, delta_i_end, y0number):


	M = y0number
	N = 500000

	Lambda = 0.74
	Omega_m0 = 0.26
	Omega_K = 1 - Omega_m0 - Lambda

	delta_i = 1.7e-3

	y0 = np.linspace(y0max, y0end, M)

	avir_arr = np.zeros(M)

	delta_i = np.linspace(delta_i_end, delta_i_max, M)

	intervall = y0number*0.1

	for j in range(len(delta_i)):		

		k = 0

		runer = True

		totalcoll = 0

		print delta_i[j]

		for i in range(len(y0)):

			yi = y0[i]

			k += 1
			
			if k >= intervall:
			
				percentage = float(i)/float(M)*100
				total = float(j)/float(M)*100

				print "delta_i = %.2e @ %.0f%%, total @ %.0f%%" % (delta_i[j], percentage, total)
				k = 0
			
			r0 = np.exp(yi)
			drdx0 = np.exp(yi)
			y = np.linspace(yi, -1e-15, N)

			radius = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_i[j]))

			rad, avir_arr[i], coll = virialcheck(y, radius, Omega_m0, Omega_K, Lambda, delta_i[j], runer)

			if not coll:
				totalcoll += 1

		if totalcoll >= 1:
			print "It did not collapse at any initial redshift!"

		mpl.plot(y0, avir_arr, "-.", linewidth = 0.75, label =  r"$\delta_i =$ %.1e" % delta_i[j])

	mpl.xlabel(r"$y_0$")
	mpl.ylabel(r"$a_{vir}$")
	mpl.legend()
	mpl.show()

y0max = float(sys.argv[1])
y0end = float(sys.argv[2])

y0number = int(sys.argv[3])

delta_i_max = float(sys.argv[4])
delta_i_end = float(sys.argv[5])

crosschecker(y0max, y0end, delta_i_max, delta_i_end, y0number)