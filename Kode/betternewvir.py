#import numpy for array handling
import numpy as np				
#import matplotlib for figures											
import matplotlib.pylab as mpl
#import odeint from scipy to solve DE
from scipy.integrate import odeint
#import sys to take redshifts from the command line
import sys

"""
H_0 = 3.44e-35				#eV
G = 6.71e-57				#eV^-2
M_pl = 2.44e27				#eV
rho_c0 = 3*H_0**2/(8*np.pi*G)
Msun = 2e30*5.61e26
GtimeM_sun = G*Msun
"""
#in natural (Planck) units:
H_0 = 3.44e-35										#eV 
G = 1.0												#1
rho_c0 = 3.*H_0**2/(8.*np.pi*G)
Msun = 5.291e12 #2.357e14							#eV-1
M_pl = 1/np.sqrt(8*np.pi)
GtimeM_sun = G*Msun

"""
G = 4.73e-58										#eV^-2, revisit the value of G in NU
Msun = 1.99e30*1.78e36								#eV, should be 1.12e66
M_pl = 1/np.sqrt(8*np.pi*G)
H_0 = 3.44e-35										#eV
rho_c0 = 3.*H_0**2/(8.*np.pi*G)
GtimeM_sun = G*Msun
"""
print H_0, G, rho_c0, Msun, M_pl

def vircheck(R, y, Omega_m0, Lambda, delta_i, acceleration, E, gamma, beta, M, f_R0, delta_rho, n, file):

	a = np.exp(y)

	rad = R[:, 0]
	drdy = R[:, 1]
	ddrddy = np.zeros(len(a))
	for i in range(len(a)):
		ddrddy[i] = acceleration( Omega_m0, Lambda, delta_i, a[i], rad[i], drdy[i], E, gamma, beta, M, f_R0, delta_rho, n )
	s = j = np.argmax(rad)
	W_ta = potential(rad[s], drdy[s], ddrddy[s], a[s], E, Omega_m0, Lambda, gamma, beta, delta_rho, M, r, delta_i, n )
	T = kinetic(drdy)
	W = potential(rad, drdy, ddrddy, a, E, Omega_m0, Lambda, gammaBack, beta, delta_rho, M, r, delta_i, n)

	

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
		file.write("{:1.5e} \n".format(a[k]))
	else:
		file.write("These parameters do not lead to collapse. \n")
	return rad, T, W


def kinetic( drdy ):
	return 3./10.*drdy**2

def potential( radius, dotR, ddotR, a, E, Omega_m0, Lambda, gamma, beta, delta_rho, M, r, delta_i, n ):
	
	W = 3./5.*radius*(-3*gammaBack(beta)*Omega_m0/(2*a**3*E(Omega_m0, Lambda, a, gamma, beta))*dotR + ddotR)
	return W


def r( x, y, acceleration, Omega_m0, Lambda, delta_i, E, gamma, beta, M, f_R0, delta_rho, n ):
	#input function to be used by odeint. Generic form for a seccond order ode
	rad = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = acceleration( Omega_m0, Lambda, delta_i, a, rad, drdy, E, gamma, beta, M, f_R0, delta_rho, n )	

	return rr



def findcoll(tolerance, acceleration, model, E, gamma, beta, M, f_R0, delta_rho, n, y0, file):

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

	delta_max = 0.001
	delta_min = 0.00000001
	j = 0
	c = 0

	while ( abs(colltime_max) >= abs(tol) ):

		#set bisection point
		delta_mid = (delta_max + delta_min)/2.0

		# solve for deltamax
		radiusmax = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_max, E, gamma, beta, M, f_R0, delta_rho, n) )


		#solve for deltamin
		radiusmin = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_min, E, gamma, beta, M, f_R0, delta_rho, n) )


		#solve for deltamid
		radiusmid = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_mid, E, gamma, beta, M, f_R0, delta_rho, n) )

		for i in range(len(radiusmax[:,0])):
			if radiusmax[i,0] <= 0:
				colltime_max = y[i]
				collmax = True
				#assert np.isnan(radiusmax[i,0]), "radius is %.3e, delta_max = %.3e" % (radiusmax[i,0], delta_max)
				#print "max"
				break
			else:
				collmax = False

		for i in range(len(radiusmid[:,0])):
			if radiusmid[i,0] <= 0:
				colltime_mid = y[i]
				collmid = True
				#print "mid"
				break
			else:
				collmid = False

		for i in range(len(radiusmin[:,0])):
			if radiusmin[i,0] <= 0:
				colltime_min = y[i]
				collmin = True
				#print "min"
				break
			else:
				collmin = False


		if ( collmax and collmid ):
			#check wether deltamax and deltamid gives collapse
			#sets new delta_max to delta_mid for next iteration
			delta_max = delta_mid


		elif ( not collmid ):
			#checks wether deltamax gives collapse, but deltamid does not
			#sets new deltamin to deltamid for next iteration
			delta_min = delta_mid

		elif ( not collmax):
			print "Increase delta_max"
			break

		if ( collmin ):
		
			print "Well damn.."
			#delta_max = 1e-15
			delta_min *= 1e-7
			c += 1
			if c > 10:
				break
			
			

		#set x_coll, which is returned from findcoll()
		x_coll = colltime_max
		"""
		collmax = False 
		collmin = False 
		collmid = False
		"""
		

		j += 1

		if j > 100:
			print "This may be infinite"
			break


	ct = np.exp(-x_coll) -1
	file.write("{:1.7f} \n".format(delta_max))
	file.write("{:1.7e} \n".format(ct))
	print delta_max

	R = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_max, E, gamma, beta, M, f_R0, delta_rho, n ) )
	k = np.argmin(R[:,0])
	print k, R[k,0]


	#mpl.plot(y, R[:,0], linewidth = 0.75, label = model)
	

	rvir, T, W = vircheck( R, y, Omega_m0, Lambda, delta_max, acceleration, E, gamma, beta, M, f_R0, delta_rho, n, file )
	if model == "Chameleon":
		mpl.plot(y, rvir, "-.", linewidth = 0.75, label = r"Chameleon, M = %1.0e, $f_{R0} - 1$ = %1.0e, $\delta_i$ = %1.5e" % (M, f_R0 - 1, delta_max))
		#mpl.plot(y, R[:,0], "-.", linewidth = 0.75, label = "Chameleon no vir")
		#Rmin = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_min, E, gamma, beta, M, f_R0, delta_rho, n ) )
		#rvirmin, Tmin, Wmin = vircheck( Rmin, y, Omega_m0, Lambda, delta_min, acceleration, E, gamma, beta, M, f_R0, delta_rho, n, file )
		#mpl.plot(y, rvirmin, "-.", linewidth = 0.75, label = r"Chameleon, M = %1.0e, $f_{R0} - 1$ = %1.0e, $\delta_i$ = %1.5e" % (M, f_R0 - 1, delta_min))
	elif model == "LCDM":
		mpl.plot(y, rvir, "--", linewidth = 0.75, label = r"$\Lambda$CDM, $\delta_i$ = %1.5e" % delta_max)
	elif model == "EdS":
		mpl.plot(y, rvir, ":", linewidth = 0.75, label = model)

	return T, W, y, delta_max

def gammaBack(beta):
	return 1 + 2*beta**2

def d_rho(Omega_m0, Lambda, delta_i, M, R, beta, n, f_R0):
	return Omega_m0*(1 + delta_i)/R**3 								#This!
	

def d_rho_complex(Omega_m0, Lambda, delta_i, M, a, beta, n, f_R0):
	R0 = 3.*H_0**2*(Omega_m0 + 4*Lambda)
	R_c = np.sqrt(delta_i)*(n + 1)*abs(1 - f_R0)/(2*beta**2*R0)*( (Omega_m0 + 4*Lambda)/(Omega_m0/a**3 + 4*Lambda) )**(n + 2)
	#R_c = 1e-4
	return 3.*M/(4.*rho_c0*np.pi*R_c**3*(2*M/(H_0**2*Omega_m0*(1+delta_i))))

def LCDMacc( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta, M, f_R0, delta_rho, n ):
	return (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/E(Omega_m0, Lambda, a, 0, 0)

def ELCDMnorad(Omega_m0, Lambda, a, gamma, beta):
	return Omega_m0/a**3 + Lambda

def chameleonacc( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta, M, f_R0, delta_rho, n ):
	return  ( Lambda*r - gamma(Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho, r, a, True)*Omega_m0*(1 + delta_i)/(2*r**2) + \
		3./2.*gammaBack(beta)*Omega_m0*drdy/a**3)/E(Omega_m0, Lambda, a, gammaBack, beta)

def Echamnorad( Omega_m0, Lambda, a, gamma, beta,  ):
	return gammaBack(beta)*Omega_m0/a**3 + Lambda

def Phi_N( r, M, Omega_m0, delta_i ):
	#assert r < 0, "Unphysical radius in Phi_N for %.2e" % delta_i
	Phi_N = ( (GtimeM_sun*M*H_0)**(2./3.)*(Omega_m0*(1 + delta_i)/2.)**(1./3.)/r )
	#Phi_N = (M/(np.sqrt(3)*np.pi))**(2./3.)*(2*rho_c0*Omega_m0*(1 + delta_i))**(1./3.)/(M_pl**2*r)
	return Phi_N

"""
def DRoverR( Omega_m0, Lambda, delta_i, beta, ):
	return abs(1-f_R0)/(12*beta**2*Phi)*( (Omega_m0 + 4*Lambda)/(delta - Omega_m0/a**3) )**(n + 1)
"""

first = True
values = open("Numbers\Vals.txt", "w")
def gammaHuSawicki1( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho, r, a, perturbation ):
	
	#print Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho(M), r, a
	Phi = Phi_N(r, M, Omega_m0, delta_i)
	delta = delta_rho(Omega_m0, Lambda, delta_i, M, r, beta, n, f_R0)
	eps = 1e-40


	if not isinstance(a, float):
		for i in range(len(a)):
			assert np.isnan(r), "r is NaN"
			if ( beta**2 <= 0 or Phi[i] < eps or (delta[i] - Omega_m0/a**3) < 0 ):
				assert (delta - Omega_m0/a**3 < 0), "Delta_rho is probably NaN"
				deltRoverR = 0
			else:
				deltRoverR = abs(1-f_R0)/(12*beta**2*Phi[i])*( (Omega_m0 + 4*Lambda)/(delta[i] - Omega_m0/a[i]**3) )**(n + 1)

	else:
		if ( beta**2 <= 0 or Phi < eps or delta - Omega_m0/a**3 < 0 ):
			deltRoverR = 0
		elif r < 0:
			deltRoverR = 0
		else:
			#assert np.isnan(r), "r is %.2e" % r
			deltRoverR = abs(1-f_R0)/(12*beta**2*Phi)*( (Omega_m0 + 4*Lambda)/(delta - Omega_m0/a**3) )**(n + 1)
	#print type(Phi_N), type(delta), type(deltRoverR), type(delta_i)
	#values.write("{:1.3e}        {:2.3e}         {:2.3e}         {:2.3e} \n".format(float(Phi_N), float(delta), float(deltRoverR), float(delta_i)))


	return 1 + 2*beta**2*( 1 - deltRoverR**3 )
	#return 1 - 2*beta**2*( 1- abs(1 - f_R0)*r_iovera_i*r/(12*beta*G*M)*( (Omega_m0 + 4*Lambda)/(delta_rho(Omega_m0, delta_i, r, r_iovera_i) - Omega_m0/a**3) )**(n/(n + 1)) )


	

def R_c(Omega_m0, Lambda, delta_i, phi_inf, n, f_R0):
	R_0 = 3*H_0*(Omega_m0 + 4*Lambda)
	return np.sqrt(delta_i)*phi_inf**2*2*(n + 1)/(M_pl**2*(f_R0 - 1)*R_0)*(M_pl*(1 - f_R0)/(2*beta*phi_inf))**(n/(n + 1))

def controll( model1, model2, E1, E2, Gamma, beta, M, f_R0, n, delta_rho, delta_i ):

	N = 5000000
	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	f1 = odeint( r, [r0, drdx0], y, args = ( model1, Omega_m0, Lambda, delta_i, E1, Gamma, beta, M, f_R0, delta_rho, n ) )
	f2 = odeint( r, [r0, drdx0], y, args = ( model2, Omega_m0, Lambda, delta_i, E2, Gamma, beta, M, f_R0, delta_rho, n ) )

	diff = f1[:,0] - f2[:,0]

	mpl.plot(y, diff, "-.", linewidth = 0.75, label = r"$R_{\Lambda CDM} - R_{Cham}$ with $\beta$ = %.2f" % beta)
	mpl.xlabel("ln(a)")
	mpl.ylabel("Diff")
	return

tolerance = 0.001
y0 = np.log(1e-4)
a_i = np.exp(y0)

#findcoll(tolerance, LCDMacc, "LCDM", ELCDMnorad, y0)
#findcoll(tolerance, LCDMacc, "EdS", ELCDMnorad, y0)


M = 1e14

beta = 1/np.sqrt(6)
f_R0 = 1.00001
n = 0.02
Omega_m0 = 0.25
Lambda = 0.75

fileLCDM = open("Numbers\LCDMviracc.txt", "w")

fileEdS = open("Numbers\EdSvirracc.txt", "w")
filecham = open("Numbers\Chameleon10.txt", "w")
print "Working on chameleon"
Tcham, Wcham, y, d_cham = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), M, 1.00001, d_rho, 0.02, y0, filecham)
a = np.exp(y)



print "Working on LCDM"
T1, W1, y, d_LCDM = findcoll(tolerance, LCDMacc, "LCDM", ELCDMnorad, 0, 0, 0, 0, 0, 0, y0, fileLCDM)
print "Working on EdS"
T2, W2, y, d_EdS = findcoll(tolerance, LCDMacc, "EdS", ELCDMnorad, 0, 0, 0, 0, 0, 0, y0, fileEdS)

fileLCDM.close()

fileEdS.close()

filecham.close()

mpl.xlabel("ln(a)")
mpl.ylabel(r"$\tilde{R}$", rotation = 0)
mpl.legend()
mpl.savefig("Figures\Evolution.png", dpi = 1000)
mpl.clf()

"""
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
mpl.savefig("Figures\Energies.pdf")
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
mpl.savefig("Figures\RelativeEnergies.pdf")
mpl.clf()

mpl.clf()

"""
N = 5000000
y = np.linspace(y0, -1e-15, N)
r0 = np.exp(y0)
drdx0 = np.exp(y0)


masses = [ 1e10, 1e12, 1e14, 1e16, 1e20]
#want to loop over masses, but not within the rootfinding. I want a figure showing how Delta=rho_p/rho_b varies with the mass and curvature
print "On loop over M"
Mloop = open("Numbers\Massloop.txt", "w")
for Mass in masses:
	Mloop.write("M = {:1.1e} \n".format(Mass))
	print  GtimeM_sun*Mass, H_0**2
	rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gammaHuSawicki1, beta, Mass, f_R0, d_rho, n) )
	rvir, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gammaHuSawicki1, beta, Mass, f_R0, d_rho, n, Mloop )
	mpl.plot(y, rvir, linewidth = 0.75, label = r"$M_c$ = %1.1e" % Mass)
Mloop.close()
mpl.xlabel("ln(a)")
mpl.ylabel(r"$\tilde{R}$", rotation = 0)
mpl.legend()
mpl.savefig("Figures\Massloop.png")
mpl.clf()


print "Loop over f_R0"
M = 1e14
f_R0_list = np.array([1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]) + 1
f_R0_loop = open("Numbers\Loop_f_R0.txt", "w")
for f_R0_i in f_R0_list:
	f_R0_loop.write("f_R0 = {:1.8f} \n".format(f_R0_i))
	rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gammaHuSawicki1, beta, M, f_R0_i, d_rho, n) )
	rvir, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gammaHuSawicki1, beta, M, f_R0_i, d_rho, n, f_R0_loop )
	mpl.plot(y, rvir, linewidth = 0.75, label = r"$f_{R0}-1$ = %1.1e" % (f_R0_i-1))
f_R0_loop.close()
mpl.xlabel("ln(a)")
mpl.ylabel(r"$\tilde{R}$", rotation = 0)
mpl.legend()
mpl.savefig("Figures\Loop_f_R0.pdf")
mpl.clf()

print "Loop over n"
M = 1e14
f_R0 = 1 + 1e-5
n_loop = open("Numbers\Loop_n.txt", "w")
n_list = np.linspace(1e-4, 1, 5)
for ns in n_list:
	n_loop.write("n = {:1.8f} \n".format(ns))
	rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gammaHuSawicki1, beta, M, f_R0, d_rho, ns) )
	rvir, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gammaHuSawicki1, beta, M, f_R0, d_rho, ns, n_loop )
	mpl.plot(y, rvir, linewidth = 0.75, label = r"$n$ = %1.1e" % (ns))
n_loop.close()
mpl.xlabel("ln(a)")
mpl.ylabel(r"$\tilde{R}$", rotation = 0)
mpl.legend()
mpl.savefig("Figures\Loop_n.pdf")
mpl.clf()

file = open("Numbers\Random.txt", "w")
rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gammaHuSawicki1, beta, M, f_R0, d_rho, n) )
rvir, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gammaHuSawicki1, beta, M, f_R0, d_rho, n, file )
file.close
a = np.exp(y)
geff = np.zeros(len(a))
for i in range(len(a)):
	geff[i] = gammaHuSawicki1(Omega_m0, Lambda, d_LCDM, beta, n, M, f_R0, d_rho, rvir[i], a[i], True) - 1
values.close()
mpl.plot(a, geff, "-.", linewidth = 0.75)
mpl.ylabel(r"$\gamma - 1$")
mpl.xlabel("a")
mpl.show()

Delta = d_rho( Omega_m0, Lambda, d_cham, M, rvir, beta, n, f_R0 )

mpl.plot(a, Delta, "-.", linewidth = 0.75)
mpl.ylabel(r"$\Delta_{\rho}$")
mpl.xlabel("a")
mpl.show()

phi = Phi_N( rvir, M, Omega_m0, d_cham )

mpl.plot(a, phi, "-.", linewidth = 0.75)
mpl.ylabel(r"$\Phi_N$")
mpl.xlabel("a")
mpl.show()
print type(Delta), type(phi), len(a), len(Delta), len(phi)


print "----"
controll( LCDMacc, chameleonacc, ELCDMnorad, Echamnorad, gammaHuSawicki1, 0, M, f_R0, n, d_rho, d_cham)

print "----"
controll( LCDMacc, chameleonacc, ELCDMnorad, Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), M, f_R0, n, d_rho, d_cham)

mpl.legend()
mpl.show()