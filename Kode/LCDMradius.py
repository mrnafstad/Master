import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint, quad
import sys

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

		if abs(3*Omega_m0*(1+delta_i)/10*(1/(2*r[s])- 1/r[t]) +Lambdavar*(3*r[s]**2 - r[t]**2)) <= 1e-4:
			controll = False
			p = s			
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

	eps = 1e-13

	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	outputvals = np.zeros(5)

	if model == "EdS":
		Omega_m0 = 1.0
		Omega_K = 0.0
		Lambda = 0.0
	else:
		Omega_m0 = 0.26
		Omega_K = 0.0
		Lambda = 0.74

	for i in range(30):

		delta_mid = (delta_max + delta_min)/2.0
		#print "dmin: {:.15e}, dmax: {:.15e}".format(delta_min, delta_max)

		#for deltamax
		radiusmax = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_max), mxstep = 5000000)

		radmax, collmax, colltime_max = virialcheck(y, radiusmax, Omega_m0, Omega_K, Lambda, delta_max, model)

		#for deltamin
		radiusmin = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_min), mxstep = 5000000)

		radmin, collmin, colltime_min = virialcheck(y, radiusmin, Omega_m0, Omega_K, Lambda, delta_min, model)

		#for deltamid
		radiusmid = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_mid), mxstep = 5000000)

		radmid, collmid, colltime_mid = virialcheck(y, radiusmid, Omega_m0, Omega_K, Lambda, delta_mid, model)

		if ( collmax and collmid ):
			delta_max = delta_mid


		elif ( collmax and not collmid ):
			delta_min = delta_mid

		x_coll = colltime_max


	return delta_min, delta_max, x_coll



N = 500000
 
y0 = float(sys.argv[1])					

y = np.linspace( y0, -1e-15, N)

file = open("Numbers\LCDMvalues.txt", "w")

file.write("Overdensity at virialization     Overdensity at maximum radius		Radius ratio 		        Time ratio	           Z_coll			Initial Overdensity \n")

Lambda = 0.74
Omega_m0 = 0.26

Omega_K = 0



r0 = np.exp(y0)
drdx0 = np.exp(y0)

delta_EdS = 1.2e-3



#Background evolution
backrad = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, 0))

radback, coll, time = virialcheck(y, backrad, Omega_m0, Omega_K, Lambda, 0, "Background")

mpl.plot(y, radback, "-c", linewidth = 0.75, label = "Background")





acceptance = 0.0001



dmin = 0.0	
dmax = 0.01

colltime = 10		#just a dummy for the while loop

while abs(colltime) > acceptance:
	dmin, dmax, colltime = findcoll(dmax, dmin, y0, "LCDM")
	diff = abs(abs(colltime) - acceptance)
	
fitt = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, dmax))
overvir, coll, ct = virialcheck(y, fitt, Omega_m0, Omega_K, Lambda, dmax, "LCDM")

mpl.plot(y, overvir, "-.", linewidth = 0.75, label  = r"Fitted $\Lambda$CDM, $\delta_i = $ %.5e" % dmax)

radius = odeint(r, [r0, drdx0], y, args = (1.0, 0.0, 0.0, r0, dmax))

rad, coll, time = virialcheck(y, radius, 1.0, 0, 0, dmax, "EdS")


mpl.plot(y, rad, ":b", linewidth = 0.75, label = r"EdS ($\Lambda$CDM-fitted), $\delta_i =$ %.5e" % dmax )



mpl.xlabel("x = ln(a)")
mpl.ylabel("r(a)")
mpl.legend( loc=2)

mpl.title(r"$x_0 =$ %0.2f, varying values of $\delta_i$" % y0)
mpl.show()

dmin = 0.0	
dmax = 0.01

colltime = 10		#just a dummy for the while loop

while abs(colltime) > acceptance:
	dmin, dmax, colltime = findcoll(dmax, dmin, y0, "EdS")
	diff = abs(abs(colltime) - acceptance)

fitt = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, dmax), mxstep = 5000000)
overvir, coll, ct = virialcheck(y, fitt, Omega_m0, Omega_K, Lambda, dmax, "LCDM")

if coll == False:
	file.write("EdS fitted perturbation does not collapse or virialize today in LCDM \n")

mpl.plot(y, overvir, "-.", linewidth = 0.75, label  = r"Fitted $\Lambda$CDM, $\delta_i = $ %.5e" % dmax)

radius = odeint(r, [r0, drdx0], y, args = (1.0, 0.0, 0.0, r0, dmax), mxstep = 5000000)

rad, coll, time = virialcheck(y, radius, 1.0, 0, 0, dmax, "EdS")


mpl.plot(y, rad, ":b", linewidth = 0.75, label = r"EdS ($\Lambda$CDM-fitted), $\delta_i =$ %.5e" % dmax )



mpl.xlabel("x = ln(a)")
mpl.ylabel("r(a)")
mpl.legend( loc=2)

mpl.title(r"$x_0 =$ %0.2f, varying values of $\delta_i$" % y0)

mpl.show()

file.close()