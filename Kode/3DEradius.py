import numpy as np
import matplotlib.pylab as mpl
from scipy.integrate import odeint, quad
import sys
from itertools import cycle

H0 = 2.2685455*10**12	#km/s/Mpc

def virialcheck(y, mix, Omega_m0, Omega_K, Lambdavar, delta_i, runer):
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

		#U[i] = (3.0*Omega_m0*(1+delta_i)/(10.0*r[i]) - 2*r[i]*r[i]*Lambdavar)/(Omega_m0/a[i]**3 + Lambdavar)
		#T[i] = rdot[i]*rdot[i]

		if r[s] <= 0:
			collapse = True
			k = s

		"""
		if T[i] <= U[i]:
			#Error seems to be in virialisation crit, but it did not entirely solve the problem. Still 1/10 off
			controll = False
			p = i
		"""
		if abs(3*Omega_m0*(1+delta_i)/10*(1/(2*r[s])- 1/r[t]) +Lambdavar*(3*r[s]**2 - r[t]**2)) <= 1e-4:
			controll = False
			p = s			
		s += 1

	if p >= 1:
		rvir = r[p]
	rover = rvir/r[0]
	avir = a[p]
	"""
	odensity = Omega_m0*(avir/rvir)**3*(1 + delta_i)
	s = np.argmax(r)
	rmax = r[s]
	amax = a[s]

	odensitymax = Omega_m0*(amax/rmax)**3*(1 + delta_i)
		
	print "rho_vir = {:.4f}, a_vir = {:.2e}, rvir = {:.2e}".format(odensity, avir, rvir)

	print "rho_max = {:.4f}, a_max = {:.2e}, rmax = {:.2e}".format(odensitymax, amax, rmax), s
	"""
	if rvir:
		odensity = (Omega_m0*(1+delta_i)/rvir**3 + Lambdavar)/(Omega_m0/avir**3 + Lambdavar)


		odensitymax = (Omega_m0*(1+delta_i)/rmax**3 + Lambdavar)/(Omega_m0/amax**3 + Lambdavar)
			


	

		file.write("        {:5.10f} 	  | 		 {:5.10e}    	|     {:5.10e}  	|	{:5.10f}	|	{}\n".format(odensity, odensitymax, rvir/rmax, avir/amax, collapse))

	if rvir:
			#print " %4.4e | %5.4e | %5.4e" % (rover, avir, delta_i)
		
		print k, p
		k = p
		while k + 1 <= len(r):
			r[k] = rvir
			k += 1
	elif runer==True:
		print "It does not collapse!", y[0]
	"""
	mpl.plot(y, U, "r--", linewidth = 0.75, label = "U")
	mpl.plot(y, T, "b--", linewidth = 0.75, label = "T")
	mpl.legend()
	mpl.title(delta_i)
	mpl.yscale("log")
	mpl.show()
	"""
	return r, avir, U, T


 




def LCDM_Kinetic(rmax, rdot, a, delta_i):

	Omega_m0 = 0.26

	density = Omega_m0*(1+delta_i)/rmax**3

	x2 = lambda r: density*rdot**2*r**2

	kin = quad(x2, 0, rmax)

	return kin[0]


def LCDM_Potential(rmax, rdot, a, delta_i, EdS):

	if EdS:
		Omega_m0 = 1.0
		Lambda = 0.0

	else:
		Omega_m0 = 0.26
		Lambda = 0.74
		
	density = Omega_m0*(1+delta_i)/rmax**3

	acc = (-Omega_m0*(1+delta_i)/(2.*rmax**2) + rmax*Lambda + 3./2.*rdot*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda)

	x3 = lambda r: 	( density*Omega_m0*r**3 * rdot/(Omega_m0/a**3 + Lambda)/a**3 - density*acc*r**3 )

	pot = quad(x3, 0, rmax)

	return pot[0]


def LCDM_overdensity(rmax, amax, rvir, avir, collapse, delta_i, EdS):

	if EdS:
		Omega_m0 = 1.0
		Lambda = 0.0

	else:
		Omega_m0 = 0.26
		Lambda = 0.74

	odensity = (Omega_m0*(1+delta_i)/rvir**3 + Lambda)/(Omega_m0/avir**3 + Lambda)

	print odensity


	odensitymax = (Omega_m0*(1+delta_i)/rmax**3 + Lambda)/(Omega_m0/amax**3 + Lambda)
			
	file.write("gen        {:5.10f} 	  | 		 {:5.10e}    	|     {:5.10e}  	|	{:5.10f}	|	{}\n".format(odensity, odensitymax, rvir/rmax, avir/amax, collapse))

	return

def general_virial(y, mix, delta_i, Kinetic, Potential, Overdensity, EdS):
	#beginning with the skin

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

	K = np.zeros(len(r))

	while s <= len(r)-1:

		if r[s] <= 0:
			r[s] = 0
			collapse = True

		K[s] = Kinetic(r[s], rdot[s], a[s], delta_i)
		U[s] = Potential(r[s], rdot[s], a[s], delta_i, EdS)

		if K[s] <= 4*U[s]:
			p = s


		s += 1

	if p > 1:
		rvir = r[p]
		avir = a[p]

		Overdensity(rmax, amax, rvir, avir, collapse, delta_i, EdS)

		while p + 1 <= len(r):
			r[p] = rvir
			p += 1

	return r, avir, U, K






def r(x, y, Omega_m0, Omega_K, Lambda, r_i, delta_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda) 	

	return rr


def EdSr(x, y, Omega_m0, r_i, delta_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = -(1 + delta_i)*a**3/(2.0*r**2) + 3.0/2.0*drdy	

	return rr




eps = 1.63e-5
N = 500000
 
y0 = float(sys.argv[4])					#-2 might be a decent number, roughly -7 corresponds to recombination

y = np.linspace( y0, -1e-15, N)

file = open("values.txt", "w")

file.write("Overdensity at virialization     Overdensity at maximum radius		Radius ratio 		        Time ratio \n")

Lambda = 0.74
Omega_m0 = 0.26

Omega_K = 0

delta_i_max = float(sys.argv[1])

delta_i_min = float(sys.argv[2])

delta_int = int(sys.argv[3])

delta_i = np.linspace(delta_i_min, delta_i_max, delta_int)

kinetic = [0, 0, 0, 0, 0, 0]
potential = [0, 0, 0, 0, 0, 0]

r0 = np.exp(y0)
drdx0 = np.exp(y0)

virialization = False
if len(sys.argv) == 7:
	virialization = True
	print "You chose virialization!"

else:
	print "This system vil not virialize!"

print "%4s | %5s" % ("r_vir/r_i", "a_vir")

runer = False

delta_EdS = 1.2e-3




print "EdS", "----"*5
radius, info = odeint(r, [r0, drdx0], y, args = (1.0, 0, 0, r0, delta_EdS), full_output = 1)

rad, avirr, potent, kinet = general_virial(y, radius, delta_EdS, LCDM_Kinetic, LCDM_Potential, LCDM_overdensity, True)

mpl.plot(y, rad, "-", linewidth = 0.75, label=r"EdS, $\delta_i =$ %.5e" % delta_EdS )
mpl.legend()
mpl.show()
print type(info)
"""
#not very important, but could give valuable info
for key,values in info:

	file.write("{} \n".format(values))
"""
rad, avirr, potential[0], kinetic[0] = virialcheck(y, radius, 1.0, 0, 0, delta_EdS, runer)


mpl.plot(y, rad, ":b", linewidth = 0.75, label = r"EdS, $\delta_i =$ %.5e" % delta_EdS )

radius = odeint(EdSr, [r0, drdx0], y, args = (1.0, r0, delta_EdS))

rad, avirr, potential[0], kinetic[0] = virialcheck(y, radius, 1.0, 0, 0, delta_EdS, runer)


mpl.plot(y, rad, ":b", linewidth = 0.75, label = r"EdS, $\delta_i =$ %.5e" % delta_EdS )		

print "Background", "----"*5
backrad = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, 0))

radback, acir, pot, kin = virialcheck(y, backrad, Omega_m0, Omega_K, Lambda, 0, runer)

mpl.plot(y, radback, "-c", linewidth = 0.75, label = "Background")

print "Loop", "----"*5




#mpl.yscale("log")
#mpl.xscale("log")



if len(sys.argv) >= 6:
	#Run through different values of the initial overdensity and plot together
	for i in range(len(delta_i)):
		
		sofa = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta_i[i]))

		if virialization:
			rad, a_virr, potential[i+1], kinetic[i+1] = virialcheck(y, sofa, Omega_m0, Omega_K, Lambda, delta_i[i], runer)

			mpl.plot(y, rad, "-.", linewidth = 0.75, label = r"$\delta_{i} =$ %.5e" % delta_i[i])

		else:
			
			mpl.plot(y, sofa, "-.", linewidth = 0.75, label = r"$\delta_{i} =$ %.5e" % delta_i[i])


		#print delta_i[i]


	file.close()

	mpl.xlabel("x = ln(a)")
	mpl.ylabel("r(a)")
	mpl.legend( loc=2)

	mpl.title(r"$x_0 =$ %0.2f, varying values of $\delta_i$" % y0)

	mpl.show()

	for i in range(len(potential)):
		if i == 0:
			mpl.plot(y, potential[i], "c-.", linewidth = 0.75, label = "U, EdS")
			mpl.plot(y, kinetic[i], "c-", linewidth = 0.75, label = "T, EdS")
			mpl.plot(y, potential[i] + kinetic[i], linewidth = 0.75, label = r"$E_{tot}$")

		else:
			mpl.plot(y, potential[i], "g-.", linewidth = 0.75, label = "U, {:.6e}".format(delta_i[i-1]))
			mpl.plot(y, kinetic[i], "g-", linewidth = 0.75, label = "T, {:.6e}".format(delta_i[i-1]))
			mpl.plot(y, potential[i] + kinetic[i], "b--", linewidth = 0.75, label = r"$E_{tot}$")	
		mpl.legend()
		mpl.yscale("log")
		mpl.show()

else:
#Want to run through different values of y0 at



	y0_arr = np.linspace(-8, -1, 10)



	for i in range(len(y0_arr)):


		delta__i = 1.686

		y = np.linspace( y0_arr[i], -1e-10, N)
		sofa = odeint(r, [r0, drdx0], y, args = (Omega_m0, Omega_K, Lambda, r0, delta__i))	

		mpl.plot(y, sofa[:,0], "--", linewidth = 0.75, label = r"$x_0 =$ %.1f, $\delta_i =$ %.4f" % (y0_arr[i], delta__i))



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

#varyinitial()