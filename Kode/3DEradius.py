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
	elif runer==True:
		print "It does not collapse!", y[0]
	return r, avir


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

	M = 1000
	N = 500000

	Lambda = 0.74
	Omega_m0 = 0.26
	Omega_K = 1 - Omega_m0 - Lambda

	delta_i = 1.7e-3

	y0 = np.linspace(-10, -6, M)

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











eps = 1.63e-5
N = 500000
 
y0 = float(sys.argv[3])					#-2 might be a decent number, roughly -7 corresponds to recombination

y = np.linspace( y0, -1e-15, N)

Lambda = 0.74
Omega_m0 = 0.26

Omega_K = 1 - Omega_m0 - Lambda

delta_i_max = float(sys.argv[1])

delta_int = int(sys.argv[2])

delta_i = np.linspace(0, delta_i_max, delta_int)

r0 = np.exp(y0)
drdx0 = np.exp(y0)

virialization = False
if len(sys.argv) == 6:
	virialization = True
	print "You chose virialization!"

else:
	print "This system vil not virialize!"

print "%4s | %5s" % ("r_vir/r_i", "a_vir")

runer = False

radius = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, 0, r0, np.mean(delta_i)))

rad, avirr = virialcheck(y, radius, Omega_m0, Omega_K, 0, np.mean(delta_i), runer)


mpl.plot(y, rad, ":b", linewidth = 0.75, label = r"EdS, $\delta_i = 0$") 		#TypeError: tuple indices must be integers, not tuple






#mpl.yscale("log")
#mpl.xscale("log")



if len(sys.argv) >= 5:
	#Run through different values of the initial overdensity and plot together
	for i in range(len(delta_i)):
		
		sofa = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_i[i]))

		if virialization:
			rad, a_virr = virialcheck(y, sofa, Omega_m0, Omega_K, Lambda, delta_i[i], runer)

			if i == 0:

				mpl.plot(y, rad, "-c", linewidth = 0.75, label = r"$\delta_{i} =$ %.1e" % delta_i[i])

			else:

				mpl.plot(y, rad, "-.", linewidth = 0.75, label = r"$\delta_{i} =$ %.1e" % delta_i[i])

		else:
			if i == 0:

				mpl.plot(y, sofa, "-c", linewidth = 0.75, label = r"$\delta_{i} =$ %.1e" % delta_i[i])

			else:

				mpl.plot(y, sofa, "-.", linewidth = 0.75, label = r"$\delta_{i} =$ %.1e" % delta_i[i])


		#print delta_i[i]




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

varyinitial()