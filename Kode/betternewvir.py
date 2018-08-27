#import numpy for array handling
import numpy as np				
#import matplotlib for figures											
import matplotlib.pylab as mpl
#import odeint from scipy to solve DE
from scipy.integrate import odeint
#import sys to take redshifts from the command line
import sys


def vircheck(R, y, Omega_m0, Lambda, delta_i, acceleration):

	a = np.exp(y)

	rad = R[:, 0]
	drdy = R[:, 1]
	ddrddy = acceleration(Omega_m0, Lambda, delta_i, a, rad, drdy)
	j = np.argmax(rad)

	while abs(2*ddrddy[j] + drdy[j]) > 1e-4:
		j += 1
		if j-1 == len(rad):
			break

	i = j
	while i < len(rad):
		rad[j] = rad[i]
		i += 1

	return rad




def r(x, y, Omega_m0, Lambda, r_i, delta_i, acceleration):
	#input function to be used by odeint. Generic form for a seccond order ode
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = acceleration(Omega_m0, Lambda, delta_i, a, r, drdy)	

	return rr



def findcoll(tolerance, acceleration, model, y0):

	#Set initial conditions and time array
	N = 5000000

	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	tol = np.log(1/(1+tolerance))

	#Set density parameters according to argument "model"
	if model == "EdS":
		Omega_m0 = 1.0
		Lambda = 0.0
	else:
		Omega_m0 = 0.25
		Lambda = 0.75

	collmax = False 
	collmin = False 
	collmid = False

	colltime_max = 10

	delta_max = 0.01
	delta_min = 0.0000001
	j = 0

	while abs(colltime_max) >= abs(tol):

		#set bisection point
		delta_mid = (delta_max + delta_min)/2.0

		# solve for deltamax
		radiusmax = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_max, acceleration))


		#solve for deltamin
		radiusmin = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_min, acceleration))


		#solve for deltamid
		radiusmid = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_mid, acceleration))

		for i in range(len(radiusmax[:,0])):
			if radiusmax[i,0] <= 0:
				colltime_max = y[i]
				collmax = True
				break

		for i in range(len(radiusmid[:,0])):
			if radiusmid[i,0] <= 0:
				colltime_mid = y[i]
				collmid = True
				break

		for i in range(len(radiusmin[:,0])):
			if radiusmin[i,0] <= 0:
				colltime_min = y[i]
				collmin = True
				break

		if ( collmax and collmid ):
			#check wether deltamax and deltamid gives collapse
			#sets new delta_max to delta_mid for next iteration
			delta_max = delta_mid


		elif ( not collmid ):
			#checks wether deltamax gives collapse, but deltamid does not
			#sets new deltamin to deltamid for next iteration
			delta_min = delta_mid

		elif ( not collmax):
			print "well shiet"

		#set x_coll, which is returned from findcoll()
		x_coll = colltime_max

		collmax = False 
		collmin = False 
		collmid = False

		j += 1

		if j > 100:
			print "This may be infinite"
			break

	R = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_max, acceleration))

	mpl.plot(y, R[:,0], linewidth = 0.75, label = model)

	rvir = vircheck(R, y, Omega_m0, Lambda, delta_max, acceleration)
	mpl.plot(y, rvir, linewidth = 0.75, label = "vir")

	return 


def LCDMacc( Omega_m0, Lambda, delta_i, a, r, drdy):
	return (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda)


tolerance = float(sys.argv[1])
y0 = np.log(1/(1+float(sys.argv[2])))

findcoll(tolerance, LCDMacc, "LCDM",y0)
findcoll(tolerance, LCDMacc, "EdS", y0)

mpl.legend()
mpl.show()