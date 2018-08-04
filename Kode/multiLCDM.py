import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint, quad
import sys

from time import *

t0 = clock()

def virialcheck(y, mix, Omega_m0, Omega_K, Lambdavar, delta_i, model):
	rdot = mix[:,1]
	r = mix[:,0]

	a = np.exp(y)
	s = 0
	p = 0
	rvir = False
	collapse = False
	k = 0

	t = s = np.argmax(r)
	rmax = r[s]
	amax = a[s]

	U = np.zeros(len(r))

	T = np.zeros(len(r))

	controll = True

	First = True

	while s < len(r) - 1:



		if r[s] <= 0:
			collapse = True
			t = s
			break

		if abs(3./10.*(Omega_m0*(1+delta_i)*(1./r[t] - 1./(2*r[s])) + Lambdavar*(r[t]**2 - 2* r[s]**2))) <= 1e-4:
			controll = False
			p = s

		#elif abs(Omega_m0*(1 + delta_i)*(1./1/))			
		s += 1

	if p >= 1:
		rvir = r[p]
	rover = rvir/r[0]
	avir = a[p]

	if rvir:
		odensity = (Omega_m0*(1+delta_i)/rvir**3 + Lambdavar)/(Omega_m0/avir**3 + Lambdavar)


		odensitymax = (Omega_m0*(1+delta_i)/rmax**3 + Lambdavar)/(Omega_m0/amax**3 + Lambdavar)


	if rvir:
		
		k = p
		while k + 1 <= len(r):
			r[k] = rvir
			k += 1

	if ( collapse and rvir ):
		file.write("        {:5.10f} 	  | 		 {:5.10e}    	|     {:5.10e}  	|	{:5.10f}	|	{:5.10e} 	|	{:5.10f}    |     {} \n".format(odensity, odensitymax, rvir/rmax, avir/amax, np.exp(-y[t]) -1, delta_i, model))

	return r, collapse, y[t]




def r(x, y, Omega_m0, Omega_K, Lambda, r_i, delta_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda) 	

	return rr

def findcoll(delta_max, delta_min, y0, model):
	
	N = 500000

	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)


	if model == "EdS":
		Omega_m0 = 1.0
		Omega_K = 0.0
		Lambda = 0.0
	else:
		Omega_m0 = 0.26
		Omega_K = 0.0
		Lambda = 0.74

	for i in range(15):

		delta_mid = (delta_max + delta_min)/2.0
		#print "dmin: {:.15e}, dmax: {:.15e}".format(delta_min, delta_max)

		#for deltamax
		radiusmax = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_max))

		radmax, collmax, colltime_max = virialcheck(y, radiusmax, Omega_m0, Omega_K, Lambda, delta_max, model)

		#for deltamin
		radiusmin = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_min))

		radmin, collmin, colltime_min = virialcheck(y, radiusmin, Omega_m0, Omega_K, Lambda, delta_min, model)

		#for deltamid
		radiusmid = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_mid))

		radmid, collmid, colltime_mid = virialcheck(y, radiusmid, Omega_m0, Omega_K, Lambda, delta_mid, model)

		if ( collmax and collmid ):
			delta_max = delta_mid


		elif ( collmax and not collmid ):
			delta_min = delta_mid

		x_coll = colltime_max


	return delta_min, delta_max, x_coll

def findcollshell(y0, model):
	#y0 corresponds to the initial redshift of the perturbation
	#model indicates which model (EdS or LCDM) we want to fitt the collapse to
	acceptance = 0.0001
	print model

	fname = "Figures\LCDMfitting"  + str(int(np.exp(-y0) - 1)) + model + ".png"

	tykk = 0.75
	#tykk = 1.0

	dmin = 0.0	
	dmax = 0.01

	colltime = 10		#just a dummy for the while loop

	while abs(colltime) > acceptance:
		dmin, dmax, colltime = findcoll(dmax, dmin, y0, model)
		diff = abs(abs(colltime) - acceptance)

	mpl.plot(y, radback, "-c", linewidth = tykk, label = "Background")
	print dmax - dmin

	fitt = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, dmax))
	overvir, coll, ct = virialcheck(y, fitt, Omega_m0, Omega_K, Lambda, dmax, "LCDM")


	mpl.plot(y, overvir, "r-.", linewidth = tykk, label  = r"Fitted $\Lambda$CDM, $\delta_i = $ %.5e" % dmax)

	radius = odeint(r, [r0, drdx0], y, args = (1.0, 0.0, 0.0, r0, dmax))

	rad, coll, time = virialcheck(y, radius, 1.0, 0, 0, dmax, "EdS")


	mpl.plot(y, rad, ":b", linewidth = tykk, label = r"EdS ($\Lambda$CDM-fitted), $\delta_i =$ %.5e" % dmax )



	mpl.xlabel("x = ln(a)")
	mpl.ylabel("r(a)")
	mpl.legend( loc=2)

	mpl.title(r"$z =$ %0.2f" % (np.exp(-y0)-1))
	#mpl.show()
	mpl.savefig(fname)

	mpl.clf()



	if coll == False:
		file.write("EdS fitted perturbation does not collapse or virialize today in LCDM \n")


	return


N = 5000000



Lambda = 0.74
Omega_m0 = 0.26

Omega_K = 0








for j in range(len(sys.argv)-1):

	y0 = np.log(1/(1 + float(sys.argv[j+1])))		

	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	filename = "Numbers\LCDMvalues" + str(sys.argv[j+1]) + ".txt"

	y = np.linspace( y0, -1e-15, N)

	file = open(filename, "w")

	file.write("Overdensity at virialization     Overdensity at maximum radius		Radius ratio 		        Time ratio	           Z_coll			Initial Overdensity \n")
	backrad = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, 0))

	radback, coll, time = virialcheck(y, backrad, Omega_m0, Omega_K, Lambda, 0, "Background")
	findcollshell(y0, "LCDM")
	file.write("-----"*100)
	file.write("\n")
	findcollshell(y0, "EdS")
	print y0

	file.close()

t1 = clock()

dt = t1-t0
print dt