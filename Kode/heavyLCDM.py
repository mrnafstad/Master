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

		if (abs(3*Omega_m0*(1+delta_i)/10*(1/(2*r[s])- 1/r[t]) +Lambdavar*(3*r[s]**2 - r[t]**2)) <= 1e-4 or abs(3*Omega_m0*(1+delta_i)/10*(1/(2*r[s])- 1/r[t])) <= abs(Lambdavar*(3*r[s]**2 - r[t]**2))):
			controll = False
			p = s
			break 			
		s += 1

	if p >= 1:
		rvir = r[p]
	rover = rvir/r[0]
	avir = a[p]

	if rvir:

		k = p
		while k + 1 <= len(r):
			r[k] = rvir
			k += 1

		odensity = (Omega_m0*(1+delta_i)/rvir**3 + Lambdavar)/(Omega_m0/avir**3 + Lambdavar)


		odensitymax = (Omega_m0*(1+delta_i)/rmax**3 + Lambdavar)/(Omega_m0/amax**3 + Lambdavar)

		vals = [[odensity, rvir, avir], [odensitymax, rmax, amax,]]

		return r, collapse, y[t], vals



		
		

	
	else:
		return r, collapse, y[t], 0




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

	outputvals = 0

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

		radmax, collmax, colltime_max, maxvals = virialcheck(y, radiusmax, Omega_m0, Omega_K, Lambda, delta_max, model)

		#for deltamin
		radiusmin = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_min))

		radmin, collmin, colltime_min, minvals = virialcheck(y, radiusmin, Omega_m0, Omega_K, Lambda, delta_min, model)

		#for deltamid
		radiusmid = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_mid))

		radmid, collmid, colltime_mid, midvals = virialcheck(y, radiusmid, Omega_m0, Omega_K, Lambda, delta_mid, model)

		if ( collmax and collmid ):
			delta_max = delta_mid
			outputvals = maxvals


		elif ( collmax and not collmid ):
			delta_min = delta_mid

		x_coll = colltime_max


	return delta_min, delta_max, x_coll, outputvals

def findcollshell(y0, model):
	#y0 corresponds to the initial redshift of the perturbation
	#model indicates which model (EdS or LCDM) we want to fitt the collapse to
	acceptance = 0.0001


	
	tykk = 0.75
	#tykk = 1.0

	dmin = 0.0	
	dmax = 0.01

	colltime = 10		#just a dummy for the while loop

	while abs(colltime) > acceptance:
		dmin, dmax, colltime, outvals = findcoll(dmax, dmin, y0, model)
		diff = abs(abs(colltime) - acceptance)

	#mpl.plot(y, radback, "-c", linewidth = tykk, label = "Background")

	fitt = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, dmax))
	overvir, coll, ct, vals2 = virialcheck(y, fitt, Omega_m0, Omega_K, Lambda, dmax, "LCDM")


	#mpl.plot(y, overvir, "r-.", linewidth = tykk, label  = r"Fitted $\Lambda$CDM, $\delta_i = $ %.5e" % dmax)

	radius = odeint(r, [r0, drdx0], y, args = (1.0, 0.0, 0.0, r0, dmax))

	rad, coll, time, vals1 = virialcheck(y, radius, 1.0, 0, 0, dmax, "EdS")


	#mpl.plot(y, rad, ":b", linewidth = tykk, label = r"EdS ($\Lambda$CDM-fitted), $\delta_i =$ %.5e" % dmax )



	#mpl.xlabel("x = ln(a)")
	#mpl.ylabel("r(a)")
	#mpl.legend( loc=2)

	#mpl.title(r"$z =$ %0.2f" % (np.exp(-y0)-1))
	#mpl.show()
	#mpl.savefig(fname)

	#mpl.clf()



	#if coll == False:
	#	file.write("EdS fitted perturbation does not collapse or virialize today in LCDM \n")


	return outvals


N = 5000000



Lambda = 0.74
Omega_m0 = 0.26

Omega_K = 0

M = 10
y0_list = np.log(1/(1 + np.linspace(2740, 1100, M)))

omax_list = np.zeros(M)
ovir_list = np.zeros(M)
rvir_list = np.zeros(M)
rmax_list = np.zeros(M)
avir_list = np.zeros(M)
amax_list = np.zeros(M)



for j in range(len(y0_list)):

	y0 = y0_list[j]		

	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	y = np.linspace( y0, -1e-15, N)


	backrad = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, 0))


	toplot = findcollshell(y0, "LCDM")


	ovir_list[j] = toplot[0][0]
	rvir_list[j] = toplot[0][1]
	avir_list[j] = toplot[0][2]

	omax_list[j] = toplot[1][0]
	rmax_list[j] = toplot[1][1]
	amax_list[j] = toplot[1][2]

mpl.subplot(231)
mpl.plot(y0_list, ovir_list, linewidth = 0.75, label = r"$\rho_{vir}/\rho_0$")
mpl.legend( loc=2)

mpl.subplot(232)
mpl.plot(y0_list, rvir_list, linewidth = 0.75, label = r"$r_{vir}$")
mpl.legend( loc=3)

mpl.subplot(233)
mpl.plot(y0_list, 1./avir_list - 1, linewidth = 0.75, label = r"$z_{vir}$")
mpl.legend( loc=2)

mpl.subplot(234)
mpl.plot(y0_list, omax_list, linewidth = 0.75, label = r"$\rho_{max}/\rho_0$")
mpl.legend( loc=2)

mpl.subplot(235)
mpl.plot(y0_list, rmax_list, linewidth = 0.75, label = r"$r_{max}$")
mpl.legend( loc=3)

mpl.subplot(236)
mpl.plot(y0_list, 1./amax_list - 1, linewidth = 0.75, label = r"$z_{max}$")
mpl.legend( loc=2)

mpl.show()