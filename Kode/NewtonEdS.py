import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint, RK45
import sys


def virialcheck(y, mix, delta_i):
	rdot = mix[:,1]
	r = mix[:,0]

	a = np.exp(y)
	s = 0
	p = 0
	rvir = False
	collapse = False
	k = t = 0

	t = s = np.argmax(r)
	rmax = r[s]
	amax = a[s]

	U = np.zeros(len(r))

	T = np.zeros(len(r))

	controll = True

	First = True
	count = 0

	while s < len(r) - 1:

		if r[s] <= 0:
			collapse = True
			t = s
			break
		
		U[s] = ( 1 + delta_i )*a[s]**3/r[s]
		T[s] = rdot[s]*rdot[s]


		if (abs(r[s] - 0.5*r[t]) <= 0.5e-5 ):
			
			controll = False
			p = s
			count +=1		


		s += 1


	if p >= 1:
		rvir = r[p]

	rover = rvir/r[0]
	avir = a[p]
	#print count


	if rvir:
		odensity = (1+delta_i)/rvir**3*avir**3


		odensitymax = (1+delta_i)/rmax**3*amax**3
			


	
		if collapse:
			file.write("        {:5.10f} 	  | 		 {:5.10e}    	|     {:5.10e}  	|	{:5.10f}	|	{:5.10e} 	|	{:5.10f} \n".format(odensity, odensitymax, rvir/rmax, avir/amax, np.exp(-y[t]) -1, delta_i))


	if rvir:

		k = p
		while k + 1 <= len(r):
			r[k] = rvir
			k += 1


	return r, collapse, y[t]





def findcoll(delta_max, delta_min, y0):
	
	N = 500000

	eps = 1e-13

	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	outputvals = np.zeros(5)



	for i in range(30):

		delta_mid = (delta_max + delta_min)/2.0
		#print "dmin: {:.15e}, dmax: {:.15e}".format(delta_min, delta_max)

		#for deltamax
		radiusmax = odeint(EdSr, [r0, drdx0], y, args = (1.0, delta_max))

		radmax, collmax, colltime_max = virialcheck(y, radiusmax, delta_max)

		#for deltamin
		radiusmin = odeint(EdSr, [r0, drdx0], y, args = (1.0, delta_min))

		radmin, collmin, colltime_min = virialcheck(y, radiusmin, delta_min)

		#for deltamid
		radiusmid = odeint(EdSr, [r0, drdx0], y, args = (1.0, delta_mid))

		radmid, collmid, colltime_mid = virialcheck(y, radiusmid, delta_mid)

		if ( collmax and collmid ):
			delta_max = delta_mid


		elif ( collmax and not collmid ):
			delta_min = delta_mid

		x_coll = colltime_max


	return delta_min, delta_max, x_coll







def EdSr(x, y, omega_mo, delta_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = -(1 + delta_i)*a**3/(2.0*r**2) + 3.0/2.0*drdy	

	return rr

file = open("EdSvalues.txt", "w")

file.write("Overdensity at virialization     Overdensity at maximum radius		Radius ratio 		        Time ratio	      Z_coll			Initial Overdensity \n")

delta_i = float(sys.argv[1])

y0 = float(sys.argv[2])
y = np.linspace( y0, -1e-10, 400000)


r0 = drdx0 = np.exp(y0)

background = odeint(EdSr, [r0, drdx0], y, args = (1.0, 0.0))

radius = odeint(EdSr, [r0, drdx0], y, args = (1.0, delta_i))
overvir, coll, ct = virialcheck(y, radius, delta_i)

mpl.plot(y, background[:,0], "c-", linewidth = 0.75, label = "Background")

mpl.plot(y, overvir, "-.", linewidth = 0.75, label  = "EdS, vir")
mpl.legend()
mpl.show()



#Implement rootfinding to find the exact density that collapses today. And consider this overdensity compared to the theoretical expectation.

acceptance = 0.0001



dmin = 0.0	
dmax = 0.01

colltime = 10		#just a dummy for the while loop

while abs(colltime) > acceptance:
	dmin, dmax, colltime = findcoll(dmax, dmin, y0)
	diff = abs(abs(colltime) - acceptance)
	
radius = RK45(EdSr, [r0, drdx0], y, args = (1.0, dmax))
overvir, coll, ct = virialcheck(y, radius, dmax)

mpl.plot(y, background[:,0], "c-", linewidth = 0.75, label = "Background")

mpl.plot(y, overvir, "-.", linewidth = 0.75, label  = "EdS, vir")
mpl.legend()
mpl.show()


file.close()