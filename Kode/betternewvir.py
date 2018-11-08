"""
cd c:\Users\Halvor\Documents\Master\Kode\
python betternewvir.py

"""
#import numpy for array handling
import numpy as np				
#import matplotlib for figures											
import matplotlib.pylab as mpl
#import odeint from scipy to solve DE
from scipy.integrate import odeint
#import sys to take redshifts from the command line
import sys
import time

"""
H_0 = 3.44e-35				#eV
G = 6.71e-57				#eV^-2
M_pl = 2.44e27				#eV
rho_c0 = 3*H_0**2/(8*np.pi*G)
Msun = 2e30*5.61e26
GtimeM_sun = G*Msun
"""
#in natural (Planck) units:
H_0 = 9.369e-32										#eV 
G = 1.0												#1
rho_c0 = 3.*H_0**2/(8.*np.pi*G)
Msun = 1.989e30*5.977e-22 #2.357e14							#eV-1
M_pl = 4.341e-9*5.977e-22 #2.176e-8*5.977e-22 
GtimeM_sun = G*Msun

"""
hbar = 6.58e-16
c = 3e8
G = 6.71e-39										#GeV^-2, revisit the value of G in NU
Msun = 1.99e30*5.61e26								#GeV, should be 1.12e66
M_pl = 2.44e18										#GeV
H_0 = 70/3.09e19*1.519e24							#GeV
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
	W_ta = potential(rad[s], drdy[s], ddrddy[s], a[s], E, Omega_m0, Lambda, gamma, beta, delta_rho, M, delta_i, n )
	T = kinetic(drdy)
	W = potential(rad, drdy, ddrddy, a, E, Omega_m0, Lambda, gamma, beta, delta_rho, M, delta_i, n)

	

	vir = False

	while ( j < len(rad) - 1 and rad[j] >= 0):

		if abs( 2*T[j] + W[j] ) <= 1e-4:

			i = k = j
			vir = True

		elif ( 2*T[j] < abs(W[j]) ):

			i = k = j
			vir = True

		j += 1
	
	if not vir:
		i = j
		print "It did not work"

	t = 0
	while rad[t] > 0:
		t += 1

	if vir:
		while i < len(rad):
			rad[i] = rad[k]
			i += 1


		odensity = (Omega_m0*(1+delta_i)/rad[k]**3)/(Omega_m0/a[k]**3)
		odensitymax = (Omega_m0*(1+delta_i)/rad[s]**3)/(Omega_m0/a[s]**3)
		odensitycoll = (Omega_m0*(1+delta_i)/rad[k]**3)/(Omega_m0/a[t]**3)

		file.write("{:1.15f} \n".format(odensity))
		file.write("{:1.15f} \n".format(odensitymax))
		file.write("{:1.15f} \n".format(odensitycoll))

		rviroverrta = rad[k]/rad[s]
		aviroverata = a[k]/a[s]

		file.write("{:1.15f} \n".format(rviroverrta))
		file.write("{:1.15f} \n".format(aviroverata))
		file.write("{:1.15f} \n".format(a[k]))
		file.write("{:1.15f} \n".format(a[t]))
	else:
		file.write("These parameters do not lead to collapse. \n")
	return rad, T, W


def kinetic( drdy ):
	return 3./10.*drdy**2

def potential( radius, dotR, ddotR, a, E, Omega_m0, Lambda, gamma, beta, delta_rho, M, delta_i, n ):
	if gamma == 0:
		W = 3./5.*radius*(-3*Omega_m0/(2*a**3*E(Omega_m0, Lambda, a, gamma, beta, n, M, f_R0, delta_rho, radius))*dotR + ddotR)

	else:
		W = 3./5.*radius*(-3*gammaBack( Omega_m0, Lambda, 0, beta, n, M, f_R0, delta_rho, a, a, True )*Omega_m0/(2*a**3*E(Omega_m0, Lambda, a, gamma, beta, n, M, f_R0, delta_rho, radius))*dotR + ddotR)
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



def findcoll(tolerance, acceleration, model, E, gamma, beta, M, f_R0, delta_rho, n, y0, plot, file):

	#Set initial conditions and time array
	N = 5000000

	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	a = np.exp(y)

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
	delta_min = 0.00000000001
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
			delta_next = delta_min*1e-7
			nextmin = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_next, E, gamma, beta, M, f_R0, delta_rho, n) )
			for i in range(len(radiusmin[:,0])):
				if radiusmin[i,0] <= 0:
					nextmintime = y[i]
					nextcoll = True
					#print "min"
					break
				else:
					nextcoll = False

			if abs(nextmintime) > abs(colltime_min):
				delta_min *= 2

			elif ( nextcoll == False or abs(nextmintime) < colltime_min): 
				delta_min = delta_next
				c += 1
				if c > 10:
					print colltime_min
					break
			
			

		#set x_coll, which is returned from findcoll()
		x_coll = colltime_max
		"""
		collmax = False 
		collmin = False 
		collmid = False
		"""
		

		j += 1

		if j > 200:
			print "This may be infinite"
			break


	ct = np.exp(-x_coll) -1
	
	#file.write("{:1.7e} \n".format(ct))
	print delta_max

	R = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_max, E, gamma, beta, M, f_R0, delta_rho, n ) )
	k = np.argmin(R[:,0])
	#print k, R[k,0]


	#mpl.plot(y, R[:,0], linewidth = 0.75, label = model)
	

	rvir, T, W = vircheck( R, y, Omega_m0, Lambda, delta_max, acceleration, E, gamma, beta, M, f_R0, delta_rho, n, file )
	if plot:
		if model == "Chameleon":
			mpl.plot(a, rvir, "-.", linewidth = 0.75, label = "Chameleon, $M = 10^{%.0f} M_{\odot} h^{-1}$, \n$f_{R0} - 1 = 10^{%.0f}$, $\delta_i$ = %1.5e" % (np.log10(M/Msun*0.7), np.log10(abs(f_R0 - 1)), delta_max))
			#mpl.plot(y, R[:,0], "-.", linewidth = 0.75, label = "Chameleon no vir")
			#Rmin = odeint(r, [r0, drdx0], y, args = ( acceleration, Omega_m0, Lambda, delta_min, E, gamma, beta, M, f_R0, delta_rho, n ) )
			#rvirmin, Tmin, Wmin = vircheck( Rmin, y, Omega_m0, Lambda, delta_min, acceleration, E, gamma, beta, M, f_R0, delta_rho, n, file )
			#mpl.plot(y, rvirmin, "-.", linewidth = 0.75, label = r"Chameleon, M = %1.0e, $f_{R0} - 1$ = %1.0e, $\delta_i$ = %1.5e" % (M, f_R0 - 1, delta_min))
		if model == "Symmetron":
			mpl.plot(a, rvir, "-", linewidth = 0.75, label = "Symmetron, $M = 10^{%.0f} M_{\odot} h^{-1}$, \n$z_{ssb} = %.1f, L = %.0f Mpc h^{-1}, \delta_i$ = %1.5e" % (np.log10(M/Msun*0.7), f_R0, n/3.09e22*1.97e-7*0.7, delta_max))
		elif model == "LCDM":
			mpl.plot(a, rvir, "--", linewidth = 0.75, label = r"$\Lambda$CDM, $\delta_i$ = %1.5e" % delta_max)
		elif model == "EdS":
			mpl.plot(a, rvir, ":", linewidth = 0.75, label = r"EdS$\delta_i$ = %1.5e" % delta_max)
	file.write("{:1.7f} \n".format(delta_max))
	#return T, W, y, delta_max

	return rvir, delta_max

def gammaBack(Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho, r, a, perturbation):
	return 1 

def d_rho(Omega_m0, Lambda, delta_i, R):
	return Omega_m0*(1 + delta_i)/R**3 								#This!
	

def LCDMacc( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta, M, f_R0, delta_rho, n ):
	return (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/E(Omega_m0, Lambda, a, 0, 0, 0, 0, 0, 0, 0 )

def ELCDMnorad(Omega_m0, Lambda, a, gamma, beta, s, t, u, v, w):
	return Omega_m0/a**3 + Lambda

def chameleonacc( Omega_m0, Lambda, delta_i, a, r, drdy, E, gamma, beta, M, f_R0, delta_rho, n ):
	return  ( Lambda*r - gamma( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho, r, a, True )*Omega_m0*(1 + delta_i)/(2*r**2) + \
		3./2.*gamma( Omega_m0, Lambda, 0, beta, n, M, f_R0, delta_rho, a, a, True )*Omega_m0*drdy/a**3)/E(Omega_m0, Lambda, a, gamma, beta, n, M, f_R0, delta_rho, r)

def Echamnorad( Omega_m0, Lambda, a, gamma, beta, n, M, f_R0, delta_rho, r ):
	return gamma( Omega_m0, Lambda, 0, beta, n, M, f_R0, delta_rho, a, a, True )*Omega_m0/a**3 + Lambda

def Phi_N( r, M, Omega_m0, delta_i ):
	#assert r < 0, "Unphysical radius in Phi_N for %.2e" % delta_i
	#Phi_N = ( (GtimeM_sun*M*H_0)**(2./3.)*(Omega_m0*(1 + delta_i)/2.)**(1./3.)/r )
	#Phi_N = (M/(np.sqrt(3)*np.pi))**(2./3.)*(2*rho_c0*Omega_m0*(1 + delta_i))**(1./3.)/(M_pl**2*r)
	#Phi_N = ( 2*M*np.sqrt(Omega_m0*(1 + delta_i))/H_0**2 )**(2./3.)*rho_c0/(6*M_pl**2*r)
	Phi_N = ((G*M*H_0)**2*Omega_m0*(1 + delta_i)/2)**(1./3.)/r
	return Phi_N

def phi_cham( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r, a, delta_rho):
	if delta_i == 0:
		delta = delta_rho(Omega_m0, Lambda, 0, a)
	else:
		delta = delta_rho(Omega_m0, Lambda, delta_i, r)

	phi = M_pl*abs(1 - f_R0)/(2*beta)*( (Omega_m0 + 4*Lambda)/(delta + 4*Lambda) )**(n + 1)

	#phi = M_pl*abs(1-f_R0)/(2*beta)*( (n + 1)**2 *(Omega_m0 + 4*Lambda)/(delta + 4*Lambda) )**(1/(n**2 + n + 1))

	return phi

def DROR( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r, a, delta_rho, perturbation ):

	if beta**2 <= 0:
		return 0

	else:
		if perturbation:			
			Phi = Phi_N(r, M, Omega_m0, delta_i)
			phi_c = phi_cham(Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r, a, delta_rho)
			phi_inf = phi_cham(Omega_m0, Lambda, 0, beta, n, M, f_R0, r, a, delta_rho)
			d = abs( phi_inf - phi_c )/(6*beta*M_pl*Phi)


		else:

			phi_c = 0
			Phi = Phi_N(a, M, Omega_m0, 0)
			phi_inf = phi_cham(Omega_m0, Lambda, 0, beta, n, M, f_R0, a, a, delta_rho)
			d = abs( phi_inf - phi_c )/(6*beta*M_pl*Phi)
			
			#d = 1
		
		try:
			for i in range(len(a)):
				if d[i] < 0:
					d[i] = 0
				elif d[i] > 1:
					d[i] = 1
		except:
			if d < 0:
				d = 0
			elif d > 1:
				d = 1

		
		return d

def phi_symm(rho, rho_ssb, L, M, beta):
	"""
	try:			
		for i in range(len(rho)):
			if rho[i] > rho_ssb:
				return 0
			else:
				return np.sqrt(np.sqrt(2)*L*beta/M_pl)*np.sqrt(rho_ssb - rho[i])
	except:
	"""
	if rho > rho_ssb:
		return 0
	else:
		return np.sqrt(np.sqrt(2)*L*beta/M_pl)*np.sqrt(rho_ssb - rho)

def DRORsymm( Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, r, a, delta_rho ):

	if beta**2 <= 0:
		return 0

	else:
		
		mu = 1/(np.sqrt(2)*L)
		a_ssb = 1./float(1. + z_ssb)
		
		M_phisquare = 2*Omega_m0*rho_c0*L**2/a_ssb**3
		l = M_pl**2 * a_ssb**3/(4*Omega_m0*rho_c0*beta**2*L**2)
		
		delta = delta_rho(Omega_m0, Lambda, delta_i, r)
		rho_ssb = rho_c0*Omega_m0/a_ssb**3
		rho_p = rho_c0*delta

		Phi = Phi_N(r, M, Omega_m0, delta_i)
		try:
			d = np.zeros(len(phi))
		except:
			d = 0	
		try:
						
			
			for i in range(len(r)):
				if rho_p[i] < rho_ssb:
					phi_inf = np.sqrt(np.sqrt(2)*L*beta/M_pl)*np.sqrt(rho_ssb - rho_p[i])
				else:
					phi_inf = 0
				#phi_c = phi_symm(rho_p, rho_ssb, L, M, beta)
				d[i] = abs(phi_inf)/(6*beta*Phi*M_pl)

				if d[i] < 0:
					d[i] = 0
				elif d[i] > 1:
					d[i] = 1
		
		
		except ZeroDivisionError:
			d = 1
		except TypeError:
			try:
				if rho_p < rho_ssb:
					phi_inf = np.sqrt(np.sqrt(2)*L*beta/M_pl)*np.sqrt(rho_ssb - rho_p)
				else:
					phi_inf = 0

				#phi_c = phi_symm(rho_p, rho_ssb, L, M, beta)
				d = abs(phi_inf)/(6*beta*Phi*M_pl)

				if d < 0:
					d = 0
				elif d > 1:
					d = 1
			except ValueError:
				d = np.zeros(len(rho_p))
				for i in range(len(rho_p)):
					if rho_p[i] < rho_ssb:
						phi_inf = np.sqrt(np.sqrt(2)*L*beta/M_pl)*np.sqrt(rho_ssb - rho_p[i])
					else:
						phi_inf = 0
					#phi_c = phi_symm(rho_p, rho_ssb, L, M, beta)
					d[i] = abs(phi_inf)/(6*beta*Phi[i]*M_pl)

					if d[i] < 0:
						d[i] = 0
					elif d[i] > 1:
						d[i] = 1
		
		return d

			#print d
		


first = True
values = open("Numbers\Vals.txt", "w")
def gammaHuSawicki1( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, delta_rho, r, a, perturbation ):
	
	deltRoverR = DROR( Omega_m0, Lambda, delta_i, beta, n, M, f_R0, r, a, delta_rho, perturbation )
	"""
	if 3*deltRoverR <= 1:
	
		return 1 + 2*beta**2*( 3*deltRoverR )

	else:
		return 1 + 2*beta**2
	"""
	Rs = 1 - deltRoverR
	return 1 + 2*beta**2*( 1 - Rs**3)


def gamma_symm( Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, delta_rho, r, a, perturbation ):
	#In comparison to chameleon, g -> beta, M_phi -> f_R0 notationally

	deltaR_overR = DRORsymm( Omega_m0, Lambda, delta_i, beta, L, M, z_ssb, r, a, delta_rho )
	"""
	if 3*deltaR_overR <= 1:

		return 1 + 2*g**2*( 3*deltaR_overR )

	else:
		return 1 + 2*g**2
	"""
	
	Rs = 1 - deltaR_overR
	return 1 + 2*beta**2*( 1 - Rs**3)

def controll( model1, model2, E1, E2, Gamma, beta, M, f_R0, n, delta_rho, delta_i ):

	N = 5000000
	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)
	a = np.exp(y)

	f1 = odeint( r, [r0, drdx0], y, args = ( model1, Omega_m0, Lambda, delta_i, E1, Gamma, beta, M, f_R0, delta_rho, n ) )
	f2 = odeint( r, [r0, drdx0], y, args = ( model2, Omega_m0, Lambda, delta_i, E2, Gamma, beta, M, f_R0, delta_rho, n ) )

	diff = f1[:,0] - f2[:,0]

	mpl.plot(a, diff, "-.", linewidth = 0.75, label = r"$R_{\Lambda CDM} - R_{Cham}$ with $\beta$ = %.2f" % beta)
	mpl.xlabel("a")
	mpl.ylabel("Diff")
	mpl.xscale("log")
	return

def placement( filename, label1, seccondvar, rootfinding ):

	file = open(filename, "r")
	d_vir = []
	d_ta = []
	d_c = []
	rvirrta = []
	avirata = []
	variable = []
	avir = []
	acoll = []
	if rootfinding:
		delta_c = []
	i = 0

	for j in file.readlines():
		try:
			value = float(j)
		except ValueError:
			break
		if i == 0:
			variable.append(value)
			i += 1
		elif i == 1:
			d_vir.append(value)
			i += 1
		elif i == 2:
			d_ta.append(value)
			i += 1
		elif i == 3:
			d_c.append(value)
			i += 1

		elif i == 4:
			rvirrta.append(value)
			i += 1
		elif i == 5:
			avirata.append(value)
			i += 1
		elif i == 6:
			avir.append(value)
			i += 1
		elif ( i == 7 and not rootfinding ):
			acoll.append(value)
			i = 0

		elif ( i == 7 and rootfinding):
			acoll.append(value)
			i += 1

		elif i == 8:
			delta_c.append(value)
			i = 0
		

	file.close()
	
	mpl.figure(1)
	mpl.plot(variable, d_vir, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$\Delta_{vir}$", rotation = 0)
	mpl.xscale("log")
	
	mpl.figure(2)
	mpl.plot(variable, d_ta, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$\Delta_{ta}$", rotation = 0)
	mpl.xscale("log")

	mpl.figure(3)
	mpl.plot(variable, rvirrta, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$\frac{\tilde{R}_{vir}}{\tilde{R}_{ta}}$", rotation = 0)
	mpl.xscale("log")

	mpl.figure(4)
	mpl.plot(variable, avirata, linewidth = 0.75, label = seccondvar)
	
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$\frac{a_{vir}}{a_{ta}}$", rotation = 0)
	mpl.xscale("log")

	mpl.figure(5)
	mpl.plot(variable, avir, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	
	mpl.ylabel(r"$a_{vir}$", rotation = 0)
	mpl.xscale("log")

	mpl.figure(6)
	mpl.plot(variable, d_c, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$\Delta_{c}$", rotation = 0)
	mpl.xscale("log")

	mpl.figure(7)
	mpl.plot(variable, acoll, linewidth = 0.75, label = seccondvar)
	mpl.legend()
	mpl.xlabel(label1)
	mpl.ylabel(r"$a_c$", rotation = 0)
	mpl.xscale("log")

	if rootfinding:
		mpl.figure(8)
		mpl.plot(variable, delta_c, linewidth = 0.75, label = seccondvar)
		mpl.legend()
		mpl.xlabel(label1)
		mpl.ylabel(r"$\delta_i$", rotation = 0)
		mpl.xscale("log")
		mpl.yscale("log")

	return 



tolerance = 0.001
y0 = np.log(1e-4)
a_i = np.exp(y0)
N = 5000000
y = np.linspace(y0, -1e-15, N)
r0 = np.exp(y0)
drdx0 = np.exp(y0)



#findcoll(tolerance, LCDMacc, "LCDM", ELCDMnorad, y0)
#findcoll(tolerance, LCDMacc, "EdS", ELCDMnorad, y0)


M = 1e14*Msun/0.7

beta = 1/np.sqrt(6)
f_R0 = 1.00001
n = 0.02
Omega_m0 = 0.25
Lambda = 0.75
"""
beta = .5
z_ssb = 1.
L = 1.
filesymm = open("Numbers\Symmetronacc.txt", "w")
Tsymm, Wsymm, y, d_symm = findcoll(tolerance, chameleonacc, "Symmetron", Echamnorad, gamma_symm, beta, M, z_ssb, d_rho, L*3.09e22/1.97e-7/0.7, y0, True, filesymm)


fileLCDM = open("Numbers\LCDMviracc.txt", "w")

fileEdS = open("Numbers\EdSvirracc.txt", "w")

filecham = open("Numbers\Chameleon10.txt", "w")
print "Working on chameleon"
Tcham, Wcham, y, d_cham = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), M, f_R0, d_rho, 1.0, y0, True, filecham)




print "Working on LCDM"
T1, W1, y, d_LCDM = findcoll(tolerance, LCDMacc, "LCDM", ELCDMnorad, 0, 0, 0, 0, 0, 0, y0, True, fileLCDM)
print "Working on EdS"
T2, W2, y, d_EdS = findcoll(tolerance, LCDMacc, "EdS", ELCDMnorad, 0, 0, 0, 0, 0, 0, y0, True, fileEdS)

filesymm.close()
filecham.close()
""
fileLCDM.close()
fileEdS.close()


mpl.xlabel("a")
mpl.ylabel(r"$\tilde{R}$", rotation = 0)
mpl.legend()
mpl.savefig("Figures\Evolution_all.pdf")#, addition_artists = mpl.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.), bbox_inches = "tight")
mpl.clf()
a = np.exp(y)

mpl.plot(a, T1, "c-", linewidth = 0.75, label = r"$T_{\Lambda CDM}$")
mpl.plot(a, -W1, "c--", linewidth = 0.75, label = r"$W_{\Lambda CDM}$")

mpl.plot(a, T2, "r:", linewidth = 0.75, label = r"$T_{EdS}$")
mpl.plot(a, -W2, "r-.", linewidth = 0.75, label = r"$W_{EdS}$")

#mpl.plot(y, Tcham, "b-", linewidth = 0.75, label = r"$T_{c10}$")
#mpl.plot(y, Wcham, "b--", linewidth = 0.75, label = r"$W_{c10}$")

mpl.yscale("log")
mpl.xscale("log")
mpl.xlabel("a")
mpl.ylabel("T / W", rotation = 0)
mpl.legend()
mpl.savefig("Figures\Energies.pdf")
mpl.clf()

relLCDM = -W1/T1
relEdS = -W2/T2
#relCham = - Wcham/Tcham

mpl.plot(a, relLCDM, "c-", linewidth = 0.75, label = r"$\frac{W_{\Lambda CDM}}{T_{\Lambda CDM}}$")
mpl.plot(a, relEdS, "r:", linewidth = 0.75, label = r"$\frac{W_{EdS}}{T_{EdS}}$")
#mpl.plot(y, relCham, "b-.", linewidth = 0.75, label = r"$\frac{W_{c10}}{T_{c10}}$")
mpl.yscale("log")
mpl.xscale("log")
mpl.xlabel("a")
mpl.ylabel(r"$-\frac{W}{T}$", rotation = 0)
mpl.legend()
mpl.savefig("Figures\RelativeEnergies.pdf")
mpl.clf()

mpl.clf()



d_cham = 5.9083e-4
masses = np.logspace(10, 18, 30)
#want to loop over masses, but not within the rootfinding. I want a figure showing how Delta=rho_p/rho_b varies with the mass and curvature
print "On loop over M"
f_R0_forM = np.array([1e-9, 1e-7, 1e-5, 1e-3]) + 1
t01 = time.time()
for f in f_R0_forM:
	t11 = time.time()
	filename = "Numbers\Gback\Massloopf_R0%.0f.txt" % np.log10(abs(1-f))
	

	Mloop = open(filename, "w")
	
	for Mass in masses:
		t21 = time.time()
		print "1 - f_R0 = 10^ %.0f" %  np.log10(abs(1 - f))
		Mloop.write("{:1.1f} \n".format(Mass))
		Tcham, Wcham, y, d_cham = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), Mass*Msun/0.7, f, d_rho, 1.0, y0, False, Mloop)
		#rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gammaHuSawicki1, beta, Mass*Msun/0.7, f, d_rho, n) )
		#rvir, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gammaHuSawicki1, beta, Mass*Msun/0.7, f, d_rho, n, Mloop )
		#mpl.plot(y, rvir, linewidth = 0.75, label = r"$M_c$ = %1.1e" % Mass)
		t22 = time.time()
		print "Time spent on individual mass: %.5e" % (t22 - t21)
	Mloop.close()
	
	t12 = time.time()
	
	print "Time spent on whole massloop: %.5e" % (t12 - t11)
	placement(filename, r"$M_{\odot}h^{-1}$", r"$|1 - f_{R0}| = 10^{%.0f}$ " %  np.log10(abs(1 - f)), True)

t02 = time.time()

print "time spent on nested massloop: %.5e" % (t02 - t01)

mpl.figure(1)
mpl.axhline(y = 294.605, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Massmaindvir.pdf")
mpl.clf()

mpl.figure(2)
mpl.axhline(y = 7.09, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Massmaindta.pdf")
mpl.clf()

mpl.figure(3)
mpl.axhline(y = 0.48181, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Massmainrvirrta.pdf")
mpl.clf()

mpl.figure(4)
mpl.axhline(y = 1.66881, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Massmainavirata.pdf")
mpl.clf()

mpl.figure(5)
mpl.axhline(y = 0.916595229171737, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Massmainavir.pdf")
mpl.clf()

mpl.figure(6)
mpl.axhline(y = 382.562, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Massmaindelta_c.pdf")
mpl.clf()

mpl.figure(7)
mpl.axhline(y = 1, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Massmaina_c.pdf")
mpl.clf()

mpl.figure(8)
mpl.axhline(y = 0.00059037, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Massmaindelta_i.pdf")
mpl.clf()


#mpl.xlabel("ln(a)")
#mpl.ylabel(r"$\tilde{R}$", rotation = 0)
#mpl.legend()
#mpl.show()
#mpl.close("all")
#mpl.savefig("Figures\Massloop.png")
#mpl.clf()




d_cham = 5.9083e-4
print "Loop over f_R0"
Mforf = [1e11, 1e13, 1e15, 1e17]
f_R0_list = np.logspace(-12, -3, 30) + 1
for M in Mforf:
	filename = "Numbers\Gback\Loop_f_R0M" + str(np.log10(M)) + "rootfinding.txt"
	f_R0_loop = open(filename, "w")
	for f_R0_i in f_R0_list:
		print "$M = %.2e M_sun/h" % M
		f_R0_loop.write("{:1.8f} \n".format(f_R0_i - 1))
		Tcham, Wcham, y, d_cham = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), M*Msun/0.7, f_R0_i, d_rho, 1.0, y0, False, f_R0_loop)
		#rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gammaHuSawicki1, beta, M*Msun/0.7, f_R0_i, d_rho, n) )
		#rvir, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gammaHuSawicki1, beta, M*Msun/0.7, f_R0_i, d_rho, n, f_R0_loop )
		#mpl.plot(y, rvir, linewidth = 0.75, label = r"$f_{R0}-1$ = %1.1e" % (f_R0_i-1))
	f_R0_loop.close()
	
	placement(filename, r"$|f_{R0} - 1|$", r"$M = 10^{%.0f} M_{\odot} h^{-1}$" % np.log10(M), True)

mpl.figure(1)
mpl.axhline(y = 294.605, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Fmaindvir.pdf")
mpl.clf()

mpl.figure(2)
mpl.axhline(y = 7.09, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Fmaindta.pdf")
mpl.clf()

mpl.figure(3)
mpl.axhline(y = 0.48181, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Fmainrvirrta.pdf")
mpl.clf()

mpl.figure(4)
mpl.axhline(y = 1.66881, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Fmainavirata.pdf")
mpl.clf()

mpl.figure(5)
mpl.axhline(y = 0.916595229171737, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Fmainavir.pdf")
mpl.clf()

mpl.figure(6)
mpl.axhline(y = 382.562, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Fmaindelta_c.pdf")
mpl.clf()

mpl.figure(7)
mpl.axhline(y = 1, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Fmainac.pdf")
mpl.clf()

mpl.figure(8)
mpl.axhline(y = 0.00059037, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\Gback\Fmaindelta_i.pdf")
mpl.clf()

mpl.close("all")
#mpl.xlabel("ln(a)")
#mpl.ylabel(r"$\tilde{R}$", rotation = 0)
#mpl.legend()
#mpl.show()
#mpl.close("all")
#mpl.savefig("Figures\Loop_f_R0.pdf")
#mpl.clf()


















masses = np.logspace(10, 18, 30)
#want to loop over masses, but not within the rootfinding. I want a figure showing how Delta=rho_p/rho_b varies with the mass and curvature
print "On loop over M"
z_ssbs = np.linspace(0.5, 4, 4)
t01 = time.time()
for z_ssb in z_ssbs:
	t11 = time.time()
	filename = "Numbers\SymmetronMassloopf_R0%.1f.txt" % z_ssb
	

	Mloop = open(filename, "w")
	
	for Mass in masses:
		t21 = time.time()
		Mloop.write("{:1.1f} \n".format(Mass))
		rvirsymm = findcoll(tolerance, chameleonacc, "Symmetron", Echamnorad, gamma_symm, 1, Mass*Msun/0.7, z_ssb, d_rho, 3.09e22/1.97e-7/0.7, y0, False, Mloop)
		#rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gammaHuSawicki1, beta, Mass*Msun/0.7, f, d_rho, n) )
		#rvir, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gammaHuSawicki1, beta, Mass*Msun/0.7, f, d_rho, n, Mloop )
		#mpl.plot(y, rvir, linewidth = 0.75, label = r"$M_c$ = %1.1e" % Mass)
		t22 = time.time()
		print "Time spent on individual mass: %.5e" % (t22 - t21), z_ssb, M
	Mloop.close()
	
	t12 = time.time()
	
	print "Time spent on whole massloop: %.5e" % (t12 - t11)
	placement(filename, r"$M_{\odot}h^{-1}$", r"$z_ssb = %.1f$ " %  z_ssb, True)

t02 = time.time()

print "time spent on nested massloop: %.5e" % (t02 - t01)

mpl.figure(1)
mpl.axhline(y = 294.605, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.title("Symmetron")
mpl.savefig("Figures\SymmetronMassmaindvir.pdf")
mpl.clf()

mpl.figure(2)
mpl.title("Symmetron")
mpl.axhline(y = 7.09, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\SymmetronMassmaindta.pdf")
mpl.clf()

mpl.figure(3)
mpl.title("Symmetron")
mpl.axhline(y = 0.48181, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\SymmetronMassmainrvirrta.pdf")
mpl.clf()

mpl.figure(4)
mpl.title("Symmetron")
mpl.axhline(y = 1.66881, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\SymmetronMassmainavirata.pdf")
mpl.clf()

mpl.figure(5)
mpl.axhline(y = 0.916595229171737, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.title("Symmetron")
mpl.legend()
mpl.savefig("Figures\SymmetronMassmainavir.pdf")
mpl.clf()

mpl.figure(6)
mpl.axhline(y = 382.562, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.title("Symmetron")
mpl.legend()
mpl.savefig("Figures\SymmetronMassmaindelta_c.pdf")
mpl.clf()

mpl.figure(7)
mpl.axhline(y = 1, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.savefig("Figures\SymmetronMassmaina_c.pdf")
mpl.clf()

mpl.figure(8)
mpl.axhline(y = 0.00059037, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.title("Symmetron")
mpl.legend()
mpl.savefig("Figures\SymmetronMassmaindelta_i.pdf")
mpl.clf()

Mforf = [1e12, 1e13, 1e15, 1e17]
z_ssb_list = np.linspace(0.5, 4, 30)
for M in Mforf:
	filename = "Numbers\SymmetronLoop_f_R0M" + str(np.log10(M)) + "rootfinding.txt"
	f_R0_loop = open(filename, "w")
	for z_ssb in z_ssb_list:
		print M, z_ssb
		f_R0_loop.write("{:1.8f} \n".format(f_R0_i - 1))
		rvirsymm = findcoll(tolerance, chameleonacc, "Symmetron", Echamnorad, gamma_symm, 1, M*Msun/0.7, z_ssb, d_rho, 3.09e22/1.97e-7/0.7, y0, False, f_R0_loop)
		#rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gammaHuSawicki1, beta, M*Msun/0.7, f_R0_i, d_rho, n) )
		#rvir, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gammaHuSawicki1, beta, M*Msun/0.7, f_R0_i, d_rho, n, f_R0_loop )
		#mpl.plot(y, rvir, linewidth = 0.75, label = r"$f_{R0}-1$ = %1.1e" % (f_R0_i-1))
	f_R0_loop.close()
	
	placement(filename, r"$z_ssb$", r"$M = 10^{%.0f} M_{\odot} h^{-1}$" % np.log10(M), True)

mpl.figure(1)
mpl.axhline(y = 294.605, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.title("Symmetron")
mpl.savefig("Figures\SymmetronFmaindvir.pdf")
mpl.clf()

mpl.figure(2)
mpl.axhline(y = 7.09, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.title("Symmetron")
mpl.savefig("Figures\SymmetronFmaindta.pdf")
mpl.clf()

mpl.figure(3)
mpl.axhline(y = 0.48181, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.title("Symmetron")
mpl.savefig("Figures\SymmetronFmainrvirrta.pdf")
mpl.clf()

mpl.figure(4)
mpl.axhline(y = 1.66881, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.title("Symmetron")
mpl.savefig("Figures\SymmetronFmainavirata.pdf")
mpl.clf()

mpl.figure(5)
mpl.axhline(y = 0.916595229171737, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.title("Symmetron")
mpl.savefig("Figures\SymmetronFmainavir.pdf")
mpl.clf()

mpl.figure(6)
mpl.axhline(y = 382.562, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.title("Symmetron")
mpl.savefig("Figures\SymmetronFmaindelta_c.pdf")
mpl.clf()

mpl.figure(7)
mpl.axhline(y = 1, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.title("Symmetron")
mpl.savefig("Figures\SymmetronFmainac.pdf")
mpl.clf()

mpl.figure(8)
mpl.axhline(y = 0.00059037, color = "red", linestyle = "--", linewidth = 0.75, label = r"$\Lambda$CDM")
mpl.legend()
mpl.title("Symmetron")
mpl.savefig("Figures\SymmetronFmaindelta_i.pdf")
mpl.clf()







mpl.close("all")















#placement("Numbers\Massloop.txt", r"$h^2 M_{\odot}$")
#placement("Numbers\Loop_f_R0.txt", r"$f_{R0} - 1$")




mpl.close("all")
#d_cham = 5.73246e-4
d_cham = 5.9083e-4

M_S = 1e-4*M_pl

file = open("Numbers\Random.txt", "w")
linestyles = ["-.", "--", ":"]
colors = ["grey", "black", "lightgrey"]


masses = np.logspace(12, 16, 3)
fs = np.logspace(-8, -4, 3) + 1
a = np.exp(y)
#rlcdm = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 0, 0, 0, d_rho, 1.0, y0, False, file)
mpl.figure(1)
#mpl.plot(a, rlcdm, "c-", label = r"$\Lambda$CDM")
beta = 1
z_ssb = 2
L = 1
for Mass in masses:
		rsymm = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gamma_symm, beta, Mass*Msun/0.7, z_ssb, d_rho, L, y0, False, file)

#Need to implement a whole new variable to everything!

		mpl.plot(a, rsymm, "r-.", linewidth = 0.75, label = r"$M = 10^{%.0f} M_{\odot} h^{-1}$" % Mass)

for i in range(len(masses)):
	rvircham = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), masses[2-i]*Msun/0.7, fs[i], d_rho, 1.0, y0, False, file)
	mpl.plot(a, rvircham, linestyle = linestyles[i], label = r"$|f_{R0} - 1| = 10^{%.0f}, M = 10^{%.0f} M_{\odot} h^{-1}$" %(np.log10(abs(fs[i] - 1)), np.log10(masses[i])))

mpl.legend()
mpl.xlabel("a")
mpl.ylabel(r"$\tilde{R}$", rotation = 0)
mpl.show()




k = 0
for M in masses:
	maxes = []
	mins = []
	j = 0
	for f_R0 in fs:
		#rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gammaHuSawicki1, beta, M, f_R0, d_rho, n) )
		#rvircham, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gammaHuSawicki1, beta, M, f_R0, d_rho, n, file )

		rvircham = findcoll(tolerance, chameleonacc, "Chameleon", Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), M*Msun/0.7, f_R0, d_rho, 1.0, y0, False, file)

		geff_cham = np.zeros(len(a))
		for i in range(len(a)):
			geff_cham[i] = gammaHuSawicki1(Omega_m0, Lambda, d_cham, beta, n, M, f_R0, d_rho, rvircham[i], a[i], True) - 1

		Delta_cham = d_rho( Omega_m0, Lambda, d_cham, rvircham )
		Phi_cham_newt = Phi_N( rvircham, M, Omega_m0, d_cham )
		phi_c = phi_cham( Omega_m0, Lambda, d_cham, beta, n, M, f_R0, rvircham, a, d_rho)
		if j == 0:
			maxes.append(geff_cham)
			maxes.append(Delta_cham)
			maxes.append(Phi_cham_newt)
			maxes.append(phi_c)


		elif j == 2:
			mins.append(geff_cham)
			mins.append(Delta_cham)
			mins.append(Phi_cham_newt)
			mins.append(phi_c)			


		mpl.figure(1)
		mpl.plot(a, geff_cham, linestyles[k], linewidth = 0.75, label = r"1 - $f_{R0} = 10^{%.0f}$" % np.log10(abs(1 - f_R0)))

		mpl.figure(2)
		mpl.plot(a, Delta_cham, linestyles[k], linewidth = 0.75, label = r"1 - $f_{R0} = 10^{%.0f}$" % np.log10(abs(1 - f_R0)))
	
		mpl.figure(3)
		mpl.plot(a, Phi_cham_newt, linestyles[k], linewidth = 0.75, label = r"1 - $f_{R0} = 10^{%.0f}$" % np.log10(abs(1 - f_R0)))

		mpl.figure(4)
		mpl.plot(a, phi_c, linestyles[k], linewidth = 0.75, label = r"1 - $f_{R0} = 10^{%.0f}$" % np.log10(abs(1 - f_R0)))
		j += 1

	

	mpl.figure(1)
	mpl.legend()
	mpl.title(r"M =$ 10^{%.0f} M_{\odot}h^{-1}$" % (np.log10(M)))
	mpl.xlabel("a")
	mpl.ylabel(r"$\gamma - 1$", rotation = 0)
	mpl.xscale("log")
	mpl.savefig("Figures\Gback\GammaMrootfind" + str(np.log10(M)) + ".pdf")
	mpl.clf()

	mpl.figure(2)
	mpl.legend()
	mpl.title(r"M =$ 10^{%.0f} M_{\odot}h^{-1}$" % (np.log10(M)))
	mpl.xlabel("a")
	mpl.ylabel(r"$\Delta_{\rho}$", rotation = 0)
	mpl.xscale("log")
	mpl.yscale("log")
	mpl.savefig("Figures\Gback\DeltaMrootfind" + str(np.log10(M)) + ".pdf")
	mpl.clf()

	mpl.figure(3)
	mpl.legend()
	mpl.title(r"M =$ 10^{%.0f} M_{\odot}h^{-1}$" % (np.log10(M)))
	mpl.xlabel("a")
	mpl.ylabel(r"$\Phi_N$", rotation = 0)
	mpl.xscale("log")
	mpl.yscale("log")
	mpl.savefig("Figures\Gback\Phi_NMrootfind" + str(np.log10(M)) + ".pdf")
	mpl.clf()

	mpl.figure(4)
	mpl.legend()
	mpl.title(r"M =$ 10^{%.0f} M_{\odot}h^{-1}$" % (np.log10(M)))
	mpl.xlabel("a")
	mpl.ylabel(r"$\phi_{cham}$", rotation = 0)
	mpl.xscale("log")
	mpl.yscale("log")
	mpl.savefig("Figures\Gback\Phi_cMrootfind" + str(np.log10(M)) + ".pdf")
	mpl.clf()

	mpl.figure(5)
	mpl.plot(a, maxes[0], linestyles[k], color = colors[k], linewidth = 0.75, label = r"$|1-f_{R0}| = \{10^{-8}, 10^{-4}\}$")
	mpl.plot(a, mins[0], linestyles[k], color = colors[k], linewidth = 0.75)
	mpl.fill_between(a, maxes[0], mins[0], facecolor = colors[k], alpha = 0.5, label = r"$M = 10^{%.0f}$" % np.log10(M), rasterized = True)

	mpl.figure(6)
	mpl.plot(a, maxes[1], linestyles[k], color = colors[k], linewidth = 0.75, label = r"$|1-f_{R0}| = \{10^{-8}, 10^{-4}\}$")
	mpl.plot(a, mins[1], linestyles[k], color = colors[k], linewidth = 0.75)
	mpl.fill_between(a, maxes[1], mins[1], facecolor = colors[k], alpha = 0.5, label = r"$M = 10^{%.0f}$" % np.log10(M), rasterized = True)

	mpl.figure(7)
	mpl.plot(a, maxes[2], linestyles[k], color = colors[k], linewidth = 0.75, label = r"$|1-f_{R0}| = \{10^{-8}, 10^{-4}\}$")
	mpl.plot(a, mins[2], linestyles[k], color = colors[k], linewidth = 0.75)
	mpl.fill_between(a, maxes[2], mins[2], facecolor = colors[k], alpha = 0.5, label = r"$M = 10^{%.0f}$" % np.log10(M), rasterized = True)

	mpl.figure(8)
	mpl.plot(a, maxes[3], linestyles[k], color = colors[k], linewidth = 0.75, label = r"$|1-f_{R0}| = \{10^{-8}, 10^{-4}\}$")
	mpl.plot(a, mins[3], linestyles[k], color = colors[k], linewidth = 0.75)
	mpl.fill_between(a, maxes[3], mins[3], facecolor = colors[k], alpha = 0.5, label = r"$M = 10^{%.0f}$" % np.log10(M), rasterized = True)

	k += 1

	#mpl.show()

mpl.figure(5)
mpl.legend()
mpl.xlabel("a")
mpl.ylabel(r"$\gamma - 1$", rotation = 0)
mpl.xscale("log")
mpl.savefig("Figures\Gback\GammaMrootfind.pdf")
mpl.clf()

mpl.figure(6)
mpl.legend()
mpl.xlabel("a")
mpl.ylabel(r"$\Delta_{\rho}$", rotation = 0)
mpl.xscale("log")
mpl.yscale("log")
mpl.savefig("Figures\Gback\DeltaMrootfind.pdf")
mpl.clf()

mpl.figure(7)
mpl.legend()
mpl.xlabel("a")
mpl.ylabel(r"$\Phi_N$", rotation = 0)
mpl.xscale("log")
mpl.yscale("log")
mpl.savefig("Figures\Gback\Phi_NMrootfind.pdf")
mpl.clf()

mpl.figure(8)
mpl.legend()
mpl.xlabel("a")
mpl.ylabel(r"$\phi_{c}$", rotation = 0)
mpl.xscale("log")
mpl.yscale("log")
mpl.savefig("Figures\Gback\Phi_cMrootfind.pdf")
mpl.clf()

mpl.close("all")

file.close()


rad = odeint(r, [r0, drdx0], y, args = ( LCDMacc, Omega_m0, Lambda, d_cham, ELCDMnorad, gammaHuSawicki1, 0, M, f_R0, d_rho, n) )
rvirLCDM, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, LCDMacc, ELCDMnorad, gammaHuSawicki1, 0, M, f_R0, d_rho, n, file )

rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gammaHuSawicki1, beta, M, f_R0, d_rho, n) )
rvircham, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gammaHuSawicki1, beta, M, f_R0, d_rho, n, file )

rad = odeint(r, [r0, drdx0], y, args = ( chameleonacc, Omega_m0, Lambda, d_cham, Echamnorad, gamma_symm, beta, M, M_S, d_rho, n) )
rvirsymm, T, W = vircheck( rad, y, Omega_m0, Lambda, d_cham, chameleonacc, Echamnorad, gamma_symm, beta, M, M_S, d_rho, n, file )

file.close()
a = np.exp(y)
mpl.plot(a, rvirLCDM, "--", linewidth = 0.75, label = r"$\Lambda$ CDM")
mpl.plot(a, rvircham, "-.", linewidth = 0.75, label = "Chameleon")
mpl.plot(a, rvirsymm, ":", linewidth = 0.75, label = "Symmetron")
mpl.legend()
mpl.ylabel(r"$\tilde{R}$", rotation = 0)
mpl.xlabel("a")
mpl.xscale("log")
mpl.show()



mpl.close("all")
"""
a = np.exp(y)
file = open("Numbers\Random.txt", "w")

beta = 0.5
z_ssb = 1.
L = 1.
rvirsymm, dsymm = findcoll(tolerance, chameleonacc, "Symmetron", Echamnorad, gamma_symm, beta, M, z_ssb, d_rho, L*3.09e22/1.97e-7/0.7, y0, False, file)
geff_symm = np.zeros(len(a))
for i in range(len(a)):
	geff_symm[i] = gamma_symm(Omega_m0, Lambda, dsymm, beta, L*3.09e22/1.97e-7/0.7, M, z_ssb, d_rho, rvirsymm[i], a[i], True) - 1
mpl.plot(a, geff_symm, "-.", label = "\gamma - 1")
mpl.legend()
mpl.show()
file.close()
"""
geff_cham = np.zeros(len(a))
DRoR_cham = np.zeros(len(a))
geff_symm = np.zeros(len(a))
DRoR_symm = np.zeros(len(a))
for i in range(len(a)):
	geff_cham[i] = gammaHuSawicki1(Omega_m0, Lambda, d_cham, beta, n, M, f_R0, d_rho, rvircham[i], a[i], True) - 1
	geff_symm[i] = gamma_symm(Omega_m0, Lambda, d_cham, beta, n, M, M_S, d_rho, rvirsymm[i], a[i], True) - 1
	DRoR_cham[i] = DROR(Omega_m0, Lambda, d_cham, beta, n, M, f_R0, rvircham[i], a[i], d_rho, True) 
	DRoR_symm[i] = DRORsymm(Omega_m0, Lambda, d_cham, beta, M, M_S, rvirsymm[i], d_rho ) 
values.close()
mpl.subplot(2, 2, 1)
mpl.plot(a, geff_cham, "-.", linewidth = 0.75, label = "Chameleon")
mpl.plot(a, geff_symm, "--", linewidth = 0.75, label = "Symmetron")
mpl.legend()
mpl.ylabel(r"$\gamma - 1$", rotation = 0)
mpl.xscale("log")
mpl.xlabel("a")



mpl.subplot(2, 2, 2)
mpl.plot(a, DRoR_cham, "-.", linewidth = 0.75, label = "Chameleon")
mpl.plot(a, DRoR_symm, "--", linewidth = 0.75, label = "Symmetron")
mpl.legend()
mpl.ylabel(r"$\frac{\Delta R}{R}$")
mpl.xlabel("a")
mpl.xscale("log")



Delta_cham = d_rho( Omega_m0, Lambda, d_cham, rvircham )
Delta_symm = d_rho( Omega_m0, Lambda, d_cham, rvirsymm )

mpl.subplot(2, 2, 3)
mpl.plot(a, Delta_cham, "-.", linewidth = 0.75, label = "Chameleon")
mpl.plot(a, Delta_symm, "--", linewidth = 0.75, label = "Symmetron")
mpl.legend()
mpl.ylabel(r"$\Delta_{\rho}$")
mpl.xscale("log")
mpl.yscale("log")
mpl.xlabel("a")


phi_cham = Phi_N( rvircham, M, Omega_m0, d_cham )
Phi_symm = Phi_N( rvirsymm, M, Omega_m0, d_cham)

mpl.subplot(2, 2, 4)
mpl.plot(a, phi_cham, "-.", linewidth = 0.75, label = "Chameleon")
mpl.plot(a, Phi_symm, "-", linewidth = 0.75, label = "Symmetron")
mpl.ylabel(r"$\Phi_N$")
mpl.legend()
mpl.xscale("log")
mpl.yscale("log")
mpl.xlabel("a")
mpl.show()
#print type(Delta), type(phi), len(a), len(Delta), len(phi)

	
d_cham = 5.9083e-4
print "----"
controll( LCDMacc, chameleonacc, ELCDMnorad, Echamnorad, gammaHuSawicki1, 0, M, f_R0, n, d_rho, d_cham)

print "----"
controll( LCDMacc, chameleonacc, ELCDMnorad, Echamnorad, gammaHuSawicki1, 1/np.sqrt(6), M, f_R0, n, d_rho, d_cham)

print "----"
controll( LCDMacc, chameleonacc, ELCDMnorad, Echamnorad, gammaHuSawicki1, 0, M, f_R0, n, d_rho, 0)

mpl.legend()
mpl.show()
"""