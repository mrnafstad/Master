#import numpy for array handling
import numpy as np				
#import matplotlib for figures											
import matplotlib.pylab as mpl
#import odeint from scipy to solve DE
from scipy.integrate import odeint
#import sys to take redshifts from the command line
import sys


def vircheck(R, y, Omega_m0, Lambda, delta_i, acceleration, E, file):

	a = np.exp(y)

	rad = R[:, 0]
	drdy = R[:, 1]
	ddrddy = acceleration(Omega_m0, Lambda, delta_i, a, rad, drdy, E, 0, 0)
	s = j = np.argmax(rad)

	T = kinetic(drdy)
	W = potential(rad, drdy, ddrddy, a, E, Omega_m0, Lambda)

	W_ta = potential(rad[s], drdy[s], ddrddy[s], a[s], E, Omega_m0, Lambda)

	vir = False

	while ( j < len(rad) - 1 and rad[j] >= 0):

		if abs( 2*T[j] + W[j] ) <= 1e-4:

			i = k = j
			vir = True


		j += 1
	
	if not vir:
		i = j
		print "It did not work"

	while i < len(rad):
		rad[i] = rad[k]
		i += 1

	odensity = (Omega_m0*(1+delta_i)/rad[k]**3)/(Omega_m0/a[k]**3)
	odensitymax = (Omega_m0*(1+delta_i)/rad[s]**3)/(Omega_m0/a[s]**3)

	file.write("{:1.5f} \n".format(odensity))
	file.write("{:1.5f} \n".format(odensitymax))

	rviroverrta = rad[k]/rad[s]
	aviroverata = a[k]/a[s]

	file.write("{:1.10e} \n".format(rviroverrta))
	file.write("{:1.5f} \n".format(aviroverata))

	return rad, T, W


def kinetic(drdy):
	return 3./10.*drdy**2

def potential(radius, dotR, ddotR, a, E, Omega_m0, Lambda):
	W = 3./5.*radius*(-3*Omega_m0/(2*a**3*E(Omega_m0, Lambda, a, 0, 0))*dotR + ddotR)
	return W


def r(x, y, Omega_m0, Lambda, r_i, delta_i, acceleration, E):
	#input function to be used by odeint. Generic form for a seccond order ode
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = acceleration(Omega_m0, Lambda, delta_i, a, r, drdy, E, 0, 0)	

	return rr



def findcoll(tolerance, acceleration, model, E, y0, file):

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
		radiusmax = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_max, acceleration, E))


		#solve for deltamin
		radiusmin = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_min, acceleration, E))


		#solve for deltamid
		radiusmid = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_mid, acceleration, E))

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

	ct = np.exp(-x_coll) -1
	file.write("{:1.7f} \n".format(delta_max))
	file.write("{:1.7e} \n".format(ct))

	R = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_max, acceleration, E))

	#mpl.plot(y, R[:,0], linewidth = 0.75, label = model)

	rvir, T, W = vircheck(R, y, Omega_m0, Lambda, delta_max, acceleration, E, file)
	mpl.plot(y, rvir, "-.", linewidth = 0.75, label = model)

	return T, W, y


def LCDMacc( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta):
	return (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/E(Omega_m0, Lambda, a, 0, 0)

def ELCDMnorad(Omega_m0, Lambda, a, gamma, beta):
	return Omega_m0/a**3 + Lambda

def chameleonacc( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta):
	return  ( Lambda*r - gamma()*Omega_m0*(1 + delta_i)/(2*r**2) + 3./2.*gamma()*Omega_m0*drdy/a**3)/E(Omega_m0, Lambda, a, gamma)

def Echamnorad(Omega_m0, Lambda, a, gamma):
	return gamma()*Omega_m0/a**3 + Lambda

def gammaHuSawicki():
	return 


tolerance = float(sys.argv[1])
y0 = np.log(1/(1+float(sys.argv[2])))
"""
findcoll(tolerance, LCDMacc, "LCDM", ELCDMnorad, y0)
findcoll(tolerance, LCDMacc, "EdS", ELCDMnorad, y0)


"""

fileLCDM = open("Numbers\LCDMviracc.txt", "w")
fileEdS = open("Numbers\EdSvirracc.txt", "w")

T1, W1, y = findcoll(tolerance, LCDMacc, "LCDM", ELCDMnorad, y0, fileLCDM)
T2, W2, y = findcoll(tolerance, LCDMacc, "EdS", ELCDMnorad, y0, fileEdS)

fileLCDM.close()
fileEdS.close()

mpl.legend()
mpl.savefig("Figures\Evolution.png")
mpl.clf()

mpl.plot(y, T1, "c-", linewidth = 0.75, label = r"$T_{\Lambda CDM}$")
mpl.plot(y, -W1, "c--", linewidth = 0.75, label = r"$W_{\Lambda CDM}$")

mpl.plot(y, 0.5*T2, "r:", linewidth = 0.75, label = r"$T_{EdS}$")
mpl.plot(y, -W2, "r-.", linewidth = 0.75, label = r"$W_{EdS}$")

mpl.yscale("log")

mpl.legend()
mpl.savefig("Figures\Energies.png")
mpl.clf()

relLCDM = -W1/T1
relEdS = -W2/T2

mpl.plot(y, relLCDM, "b-", linewidth = 0.75, label = r"$\frac{W_{\Lambda CDM}}{T_{\Lambda CDM}}$")
mpl.plot(y, relEdS, "r:", linewidth = 0.75, label = r"$\frac{W_{EdS}}{T_{EdS}}$")
mpl.yscale("log")
mpl.legend()
mpl.savefig("Figures\RelativeEnergies.png")
mpl.clf()