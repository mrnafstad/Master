import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint
import sys


def virialcheck(y, mix, delta_i):
	rdot = mix[:,1]
	r = mix[:,0]

	a = np.exp(y)
	s = 0
	p = 0
	rvir = False
	collapse = False
	k = 0

	s = np.argmax(r)
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
			k = s
			break
		
		U[s] = 3.0/5.0*( 1 + delta_i )*a[s]**3/r[s]
		T[s] = rdot[s]*rdot[s]



		if ( T[s] <= U[s] + 1e-4 and T[s] >= U[s] - 1e-4 ):
			
			controll = False
			p = s
			count +=1


		s += 1


	if p >= 1:
		rvir = r[p]

	rover = rvir/r[0]
	avir = a[p]
	print count


	if rvir:
		odensity = (1+delta_i)/rvir**3*avir**3


		odensitymax = (1+delta_i)/rmax**3*amax**3
			


	

		file.write("        {:5.10f} 	  | 		 {:5.10e}    	|     {:5.10e}  	|	{:5.10f}	|	{}\n".format(odensity, odensitymax, rvir/rmax, avir/amax, collapse))

	if rvir:
			#print " %4.4e | %5.4e | %5.4e" % (rover, avir, delta_i)
		
		print k, p
		k = p
		while k + 1 <= len(r):
			r[k] = rvir
			k += 1


	return r


def EdSr(x, y, omega_mo, delta_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = -(1 + delta_i)*a**3/(2.0*r**2) + 3.0/2.0*drdy	

	return rr

file = open("EdSvalues.txt", "w")

file.write("Overdensity at virialization     Overdensity at maximum radius		Radius ratio 		        Time ratio \n")

delta_i = float(sys.argv[1])

y0 = float(sys.argv[2])
y = np.linspace( y0, -1e-10, 400000)


r0 = drdx0 = np.exp(y0)
print "---"*10
background = odeint(EdSr, [r0, drdx0], y, args = (1.0, 0.0))
print "---"*10
radius = odeint(EdSr, [r0, drdx0], y, args = (1.0, delta_i))

mpl.plot(y, background[:,0], "--", linewidth = 0.75, label = "Background")

mpl.plot(y, radius[:,0], "-.", linewidth = 0.75, label  = "EdS")
mpl.legend()
mpl.show()

print "---"*10
radius = odeint(EdSr, [r0, drdx0], y, args = (1.0, delta_i))
print "---"*10
overvir = virialcheck(y, radius, delta_i)

mpl.plot(y, background[:,0], "--", linewidth = 0.75, label = "Background")

mpl.plot(y, overvir, "-.", linewidth = 0.75, label  = "EdS, vir")
mpl.legend()
mpl.show()

file.close()