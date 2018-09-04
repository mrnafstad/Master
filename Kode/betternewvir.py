#import numpy for array handling
import numpy as np				
#import matplotlib for figures											
import matplotlib.pylab as mpl
#import odeint from scipy to solve DE
from scipy.integrate import odeint
#import sys to take redshifts from the command line
import sys


H_0 = 70 #km/s/Mpc
G = 6.67408e-11 		#M^3/kg/s^2
M_pl = 1./np.sqrt(8*np.pi*G)

def vircheck(R, y, Omega_m0, Lambda, delta_i, acceleration, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n, file):

	a = np.exp(y)

	rad = R[:, 0]
	drdy = R[:, 1]
	ddrddy = acceleration( Omega_m0, Lambda, delta_i, a, rad, drdy, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n )
	s = j = np.argmax(rad)

	T = kinetic(drdy)
	W = potential(rad, drdy, ddrddy, a, E, Omega_m0, Lambda, gamma, beta)

	W_ta = potential(rad[s], drdy[s], ddrddy[s], a[s], E, Omega_m0, Lambda, gamma, beta)

	vir = False

	while ( j < len(rad) - 1 and rad[j] >= 0):

		if abs( 2*T[j] + W[j] ) <= 1e-4:

			i = k = j
			vir = True


		j += 1
	
	if not vir:
		i = j
		print "It did not work"

	if vir:
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
	else:
		file.write("These parameters do not lead to collapse.")
	return rad, T, W


def kinetic( drdy ):
	return 3./10.*drdy**2

def potential( radius, dotR, ddotR, a, E, Omega_m0, Lambda, gamma, beta ):
	W = 3./5.*radius*(-3*Omega_m0/(2*a**3*E(Omega_m0, Lambda, a, gamma, beta))*dotR + ddotR)
	return W


def r( x, y, acceleration, Omega_m0, Lambda, delta_i, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n ):
	#input function to be used by odeint. Generic form for a seccond order ode
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = acceleration( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n )	

	return rr



def findcoll(tolerance, acceleration, model, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n, y0, file):

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
	delta_min = 0.00000001
	j = 0

	while abs(colltime_max) >= abs(tol):

		#set bisection point
		delta_mid = (delta_max + delta_min)/2.0

		# solve for deltamax
		radiusmax = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_max, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n) )


		#solve for deltamin
		radiusmin = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_min, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n) )


		#solve for deltamid
		radiusmid = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_mid, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n) )

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

		if ( collmin ):
			if delta_min == 0:
				print "Well damn.."
				delta_max = 0
				break
			delta_min = 0
			

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
	print delta_max

	R = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_max, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n ) )

	#mpl.plot(y, R[:,0], linewidth = 0.75, label = model)
	

	rvir, T, W = vircheck( R, y, Omega_m0, Lambda, delta_max, acceleration, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n, file )
	if model == "Chameleon":
		mpl.plot(y, rvir, "-.", linewidth = 0.75, label = r"Chameleon, M = %1.1e, $\frac{r_i}{a_i}$ = %1.0f, $f_{R0}$ = %1.5e" % (M, r_iovera_i, f_R0))
	elif model == "LCDM":
		mpl.plot(y, rvir, "--", linewidth = 0.75, label = r"$\Lambda$CDM, $\delta_i$ = %1.5e" % delta_max)
	elif model == "EdS":
		mpl.plot(y, rvir, ":", linewidth = 0.75, label = model)

	return T, W, y

def d_rho(Omega_m0, delta_i, R, r_iovera_i):
	return Omega_m0*(1 + delta_i)*r_iovera_i**3/R**3


def LCDMacc( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n ):
	return (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/E(Omega_m0, Lambda, a, 0, 0)

def ELCDMnorad(Omega_m0, Lambda, a, gamma, beta):
	return Omega_m0/a**3 + Lambda

def chameleonacc( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta, M, f_R0, r_iovera_i, delta_rho, n ):
	return  ( Lambda*r - gamma(Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r_iovera_i, delta_rho, r, a, True)*Omega_m0*(1 + delta_i)/(2*r**2) + 3./2.*gamma(Omega_m0, Lambda, 0, beta, 0, 0, 0, 0, 0, 0, 0, False)*Omega_m0*drdy/a**3)/E(Omega_m0, Lambda, a, gamma, beta)

def Echamnorad( Omega_m0, Lambda, a, gamma, beta ):
	return gamma(Omega_m0, Lambda, 0, beta, 0, 0, 0, 0, 0, 0, 0, False)*Omega_m0/a**3 + Lambda

def gammaHuSawicki1( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r_iovera_i, delta_rho, r, a, perturbation ):
	if perturbation:
		return 1 - 2*beta**2*( 1- abs(1 - f_R0)*r_iovera_i*r/(12*beta*G*M)*( (Omega_m0 + 4*Lambda)/(delta_rho(Omega_m0, delta_i, r, r_iovera_i) - Omega_m0/a**3) )**(n/(n + 1)) )

	else:
		return 1 - 2*beta**2

def gammaHuSawicki( Omega_m0, Lambda, delta_i, a, M_c, R_c, f_R0, beta, R_i, a_i, phi_inf, n ):
	R_c = R_c(Omega_m0, Lambda, delta_max, phi_inf, n, f_R0)
	return 1+ 2*beta**2*( 1- abs(1 - f_R0)*r*R_i/(12*beta**2*(G*M_c*a_i))*( (Omega_m0 + 4*Lambda)/(3*M_c/(4*rho_c0*np.pi*R_c**3) - Omega_m0/a**3) ))

def R_c(Omega_m0, Lambda, delta_i, phi_inf, n, f_R0):
	R_0 = 3*H_0*(Omega_m0 + 4*Lambda)
	return np.sqrt(delta_i)*phi_inf**2*2*(n + 1)/(M_pl**2*(f_R0 - 1)*R_0)*(M_pl*(1 - f_R0)/(2*beta*phi_inf))**(n/(n + 1))


tolerance = float(sys.argv[1])
y0 = np.log(1/(1+float(sys.argv[2])))
a_i = np.exp(y0)

#findcoll(tolerance, LCDMacc, "LCDM", ELCDMnorad, y0)
#findcoll(tolerance, LCDMacc, "EdS", ELCDMnorad, y0)




"""

fileLCDM = open("Numbers\LCDMviracc.txt", "w")
fileEdS = open("Numbers\EdSvirracc.txt", "w")
filecham = open("Numbers\Chameleon10.txt", "w")

print "Working on chameleon"
Tcham, Wcham, y = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), 1e10, 1.00001, 100, d_rho, 0.1, y0, filecham)
print "Working on LCDM"
T1, W1, y = findcoll(tolerance, LCDMacc, "LCDM", ELCDMnorad, 0, 0, 0, 0, 0, 0, 0, y0, fileLCDM)
print "Working on EdS"
T2, W2, y = findcoll(tolerance, LCDMacc, "EdS", ELCDMnorad, 0, 0, 0, 0, 0, 0, 0, y0, fileEdS)

fileLCDM.close()
fileEdS.close()
filecham.close()

mpl.xlabel("ln(a)")
mpl.ylabel(r"$\tilde{R}$", rotation = 0)
mpl.legend()
mpl.savefig("Figures\Evolution.png")
mpl.clf()

mpl.plot(y, T1, "c-", linewidth = 0.75, label = r"$T_{\Lambda CDM}$")
mpl.plot(y, -W1, "c--", linewidth = 0.75, label = r"$W_{\Lambda CDM}$")

mpl.plot(y, T2, "r:", linewidth = 0.75, label = r"$T_{EdS}$")
mpl.plot(y, -W2, "r-.", linewidth = 0.75, label = r"$W_{EdS}$")

mpl.plot(y, Tcham, "b-", linewidth = 0.75, label = r"$T_{c10}$")
mpl.plot(y, Wcham, "b--", linewidth = 0.75, label = r"$W_{c10}$")

mpl.yscale("log")
mpl.xlabel("ln(a)")
mpl.ylabel("T / W", rotation = 0)
mpl.legend()
mpl.savefig("Figures\Energies.png")
mpl.clf()

relLCDM = -W1/T1
relEdS = -W2/T2
relCham = - Wcham/Tcham

mpl.plot(y, relLCDM, "c-", linewidth = 0.75, label = r"$\frac{W_{\Lambda CDM}}{T_{\Lambda CDM}}$")
mpl.plot(y, relEdS, "r:", linewidth = 0.75, label = r"$\frac{W_{EdS}}{T_{EdS}}$")
mpl.plot(y, relCham, "b-.", linewidth = 0.75, label = r"$\frac{W_{c10}}{T_{c10}}$")
mpl.yscale("log")
mpl.xlabel("ln(a)")
mpl.ylabel(r"$-\frac{W}{T}$", rotation = 0)
mpl.legend()
mpl.savefig("Figures\RelativeEnergies.png")
mpl.clf()
"""
r_iovera_i_list = [1, 10, 25, 50, 100]

print "On loop over r_i/a_i"
riloop = open("Numbers\Riloop.txt", "w")
for rioai in r_iovera_i_list:
	riloop.write("r_i/a_i = {:1.0f}".format(rioai))
	Tcham, Wcham, y = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), 1e10, 1.00001, rioai, d_rho, 0.1, y0, riloop)
riloop.close()
mpl.xlabel("ln(a)")
mpl.ylabel(r"$\tilde{R}$", rotation = 0)
mpl.legend()
mpl.savefig("Figures\Riloop.png")
mpl.clf()


print "Loop over f_R0"
eps = 5e-5
f_R0_arr = np.linspace(1 - eps, 1 + eps, 10)
f_R0_loop = open("Numbers\Loopf_R0.txt", "w")
for f_R0 in f_R0_arr:
	f_R0_loop.write("f_R0 = {:1.5e}".format(f_R0))
	Tcham, Wcham, y = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), 1e10, f_R0, 100, d_rho, 0.1, y0, f_R0_loop)

f_R0_loop.close()
mpl.xlabel("ln(a)")
mpl.ylabel(r"$\tilde{R}$", rotation = 0)
mpl.legend()
mpl.savefig("Figures\Loopf_R0.png")
mpl.clf()
