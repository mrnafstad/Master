#import numpy for array handling
import numpy as np				
#import matplotlib for figures											
import matplotlib.pylab as mpl
#import odeint from scipy to solve DE
from scipy.integrate import odeint
#import sys to take redshifts from the command line
import sys

def virialcheck(y, mix, Omega_m0, Lambdavar, delta_i, model, full):
	#This function checks wether the system virializes and collapses
	#Fetch radial variables from mix (returned ndarray from odeint)
	#returns a new radial array, a boolean "collapse" and the "redshift" at virialization
	rdot = mix[:,1]
	r = mix[:,0]

	#Make scalefactor array
	a = np.exp(y)

	#set dummy indices to zero
	s = 0
	p = 0
	k = 0
	
	#set dummy variables to false, they will be used in if-tests
	rvir = False
	collapse = False

	#Find the index of turnaround and set the ta radius and scalefactor
	t = s = np.argmax(r)
	rmax = r[s]
	amax = a[s]


	while s < len(r) - 1:
		#while loop that starts at ta, i.e. after contraction starts


		if r[s] <= 0:
			#finding collapse and breaking the loop when collapse occurs. t is the collapse index
			collapse = True
			t = s
			break

		if abs(3./10.*(Omega_m0*(1+delta_i)*(1./r[t] - 1./(2*r[s])) + Lambdavar*(r[t]**2 - 2* r[s]**2))) <= 1e-4:
			#finding virialization. p is the virialization index
			controll = False
			p = s
			if collapse:
				break

		elif abs(Omega_m0*(1+delta_i)*(1./r[t] - 1./(2*r[s]))) <= abs(Lambdavar*(r[t]**2 - 2* r[s]**2)):
			#safety elif test in case the difference set in the if test is too small, 
			#i.e. there is no index where the statement is true, but virialization will still occur
			p = s
			if collapse:
				break		
		s += 1

	#set virialization radius and scalefactor 
	if p > 1:
		#if test is just to make sure virilization actually did occur
		rvir = r[p]
		avir = a[p]

	if rvir:
		#calculate the overdensity of the perturbation at TA (odensitymax) and virialization (odensity)
		odensity = (Omega_m0*(1+delta_i)/rvir**3)/(Omega_m0/avir**3)


		odensitymax = (Omega_m0*(1+delta_i)/rmax**3)/(Omega_m0/amax**3)
		
		#set all elements in the radial array after virialization to r_vir
		k = p
		while k + 1 <= len(r):
			r[k] = rvir
			k += 1

	if ( collapse and rvir and not full):
		#write values to file
		file.write("        {:5.10f} 	  | 		 {:5.10e}    	|     {:5.10e}  	|	{:5.10f}	|	{:5.10e} 	|	{:5.10f}    |     {} \n".format(odensity, odensitymax, rvir/rmax, avir/amax, np.exp(-y[t]) -1, delta_i, model))


	return r, collapse, y[t]




def r(x, y, Omega_m0, Lambda, r_i, delta_i):
	#input function to be used by odeint. Generic form for a seccond order ode
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda) 	

	return rr

def RsoverRspecial(r, r_iovera_i, a):
	return 1 - (1.6e20/r_iovera_i)**2*r*(a**(3./(0.02 + 1)) - r**(3./(0.02 +1)))

def RsoverRgeneral(r, r_iovera_i, a, beta, n, M_pl, Omega_m0, delta_i, L, rho_c):
	return 1 - 2./(beta*Omega_m0*(1+delta_i))*(M_pl*L/(rho_c*r_iovera_i**2))*(n*M_pl*L**3/(3*beta*Omega_m0*rho_c))**(1./(n+1)) \
	*r*(a**(3./(n+1)) - r**(3./(n+1))*(1 + delta_i)**(-1./(n+1)))


def chamscalarRsoR(fact, beta, Omega_phi, Omega_phi, 	):
	return abs(fact)/12

def chameleon(x, y, Omega_m0, Omega_phi, delta_i, beta, r_iovera_i):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	#Rs = 
	g_p = 1 + 2.*beta**2*(1 - RsoverRspecial(r, r_iovera_i, a))
	g_b = 1 + 2.*beta**2
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (Omega_phi*r - g_p*Omega_m0*(1+delta_i)/(2.*r**2) + 3./2.*drdy*g_b*Omega_m0/a**3)/(g_b*Omega_m0/a**3 + Omega_phi)

	return rr


def findcoll(delta_max, delta_min, y0, model, full):
	#This is a part of the rootfinding algorithm, it is used inside findcollshell()
	#The function employs the bisection method

	#The argument "model" is used to set the density parameters
	

	#Set initial conditions and time array
	N = 5000000

	y = np.linspace(y0, -1e-15, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	#Set density parameters according to argument "model"
	if model == "EdS":
		Omega_m0 = 1.0
		Lambda = 0.0
	else:
		Omega_m0 = 0.25
		Lambda = 0.75

	for i in range(15):

		#set bisection point
		delta_mid = (delta_max + delta_min)/2.0

		# solve for deltamax
		radiusmax = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_max))
		#check virialization and collapse for deltamax, 
		#this is the reason for returning the boolean collapse in virialcheck()
		radmax, collmax, colltime_max = virialcheck(y, radiusmax, Omega_m0, Lambda, delta_max, model, full)

		#solve for deltamin
		radiusmin = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_min))
		#check virialization for deltamin
		radmin, collmin, colltime_min = virialcheck(y, radiusmin, Omega_m0, Lambda, delta_min, model, full)

		#solve for deltamid
		radiusmid = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, delta_mid))
		#check virialization for deltamid
		radmid, collmid, colltime_mid = virialcheck(y, radiusmid, Omega_m0, Lambda, delta_mid, model, full)
		

		if ( collmax and collmid ):
			#check weth410er deltamax and deltamid gives collapse
			#sets new delta_max to delta_mid for next iteration
			delta_max = delta_mid


		elif ( collmax and not collmid ):
			#checks wether deltamax gives collapse, but deltamid does not
			#sets new deltamin to deltamid for next iteration
			delta_min = delta_mid

		#set x_coll, which is returned from findcoll()
		x_coll = colltime_max



	return delta_min, delta_max, x_coll, radmax

def findcollshell(y0, model, full):
	#shell for rootfinding algorithm, may be combined into one function
	#y0 corresponds to the initial redshift of the perturbation
	#model indicates which model (EdS or LCDM) we want to fitt the collapse to

	#set the acceptance of the rootfinding, i.e. how big redshift we can accept
	#for collapse today
	acceptance = 0.0001
	#print np.exp(acceptance)-1

	#generate filename for plot based on the model and initial redshift
	fname = "Figures\LCDMproject"  + str(int(np.exp(-y0) - 1)) + model + ".png"


	#set initial minimum and maximum for rootfinding
	dmin = 0.0000001	
	dmax = 0.1

	colltime = 10		#just an initial dummy for the while loop
	diff1 = 0
	diff = 1
	while abs(colltime) > acceptance:
		diff1 = diff
		#loop through rootfind until collapse happens ~today
		dmin, dmax, colltime, fitt = findcoll(dmax, dmin, y0, model, full)
		diff = abs(abs(colltime) - acceptance)
		if diff1 == diff:
			#break if the rootfinding does not converge
			print "We may have an infinite loop on our hands"
			break
	if not full:
		#plot background
		#mpl.xlim(3000, 0)
		z = np.exp(y)
		mpl.plot(z, radback, "-c", linewidth = 0.75, label = "Background")

		if model == "EdS":
			#solve for the fitted overdensity in LCDM, and check virialization
			fitt = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, dmax))
			fitt, coll, ct = virialcheck(y, fitt, Omega_m0, Lambda, dmax, "LCDM", full)

		#mpl.xlim(2800, 0.0000001)

		#plot LCDM evolution
		mpl.plot(z, fitt, "r-.", linewidth = 0.75, label  = r"Fitted $\Lambda$CDM, $\delta_i = $ %.5e" % dmax)

		#solve for the fitted overdensity in EdS and check virialization 
		radius = odeint(r, [r0, drdx0], y, args = (1.0, 0.0, r0, dmax))

		rad, coll, time = virialcheck(y, radius, 1.0, 0, dmax, "EdS", full)


		#plot EdS evolution
		mpl.plot(z, rad, ":b",linewidth = 0.75, label = r"EdS ($\Lambda$CDM-fitted), $\delta_i =$ %.5e" % dmax )


		#create and save figure
		mpl.xlabel("x = a", labelpad = 10)
		#mpl.xscale("log")
		mpl.ylabel(r"$\~R$", rotation = 0, labelpad = 10)
		mpl.legend( loc=2)

		mpl.title(r"$z_i =$ %0.2f" % (np.exp(-y0)-1))
		mpl.savefig(fname)

		mpl.clf()



		if coll == False:
			#write a message to the file if we fitt the overdensity to EdS and generate the evolution
			#for LCDM, it will not collapse or virialize
			file.write("EdS fitted perturbation does not collapse or virialize today in LCDM \n")


	return dmax


N = 5000000


#set density parameters for the background evolution
Lambda = 0.75
Omega_m0 = 0.25






#The following exception test determines wether the program finds collapse for values given on the commandline or fitts collapse
#from matter-radiation equality to recombination.
#if there last argument on the command line is a string, fitting from mre to rec happens. Else fitting to the values given on the commandline happens
try:
	y = float(sys.argv[-1])
	print "only given values"
	Bigloop = False
	
except ValueError:
	print "only from mre to rec"
	Bigloop = True

if Bigloop == False:
	for j in range(len(sys.argv)-1):
		#loop through redshifts given on the command line

		#convert from redshift to variable x (in text)
		y0 = np.log(1/(1 + float(sys.argv[j+1])))		

		#set initial conditions
		r0 = np.exp(y0)
		drdx0 = np.exp(y0)

		#generate filename based on initial redshift
		filename = "Numbers\LCDMvaluesproject" + str(sys.argv[j+1]) + ".txt"

		#time array
		y = np.linspace( y0, -1e-15, N)

		#open, and start writing to file
		file = open(filename, "w")

		file.write("Overdensity at virialization     Overdensity at maximum radius		Radius ratio 		        Time ratio	           Z_coll			Initial Overdensity \n")
		#Solve for the background
		backrad = odeint(r, [r0, drdx0], y, args = (Omega_m0, Lambda, r0, 0))
		radback, coll, time = virialcheck(y, backrad, Omega_m0, Lambda, 0, "Background", False)
		
		#find which delta collapses today in LCDM
		findcollshell(y0, "LCDM", False)
		file.write("-----"*100)
		file.write("\n")
		#find which delta collapses today in EdS
		findcollshell(y0, "EdS", False)

		file.close()

		delta_i = 1e-4
		for i in [0.001, 0.01, 0.025, 0.05, 0.1]:
			rcham = odeint(chameleon, [r0, drdx0], y, args = (Omega_m0, Lambda, delta_i, 0.01, i))
			mpl.plot(y, rcham[:,0], "-.", linewidth = 0.75, label = r"$\frac{r_i}{a_i}$ = %.1f Mpc" % i)

		mpl.title(r"Chameleon, $\delta_i$ = %.1e" % delta_i)	
		mpl.legend()
		mpl.show()




if Bigloop:
	n = 10
	y0_arr = np.linspace(np.log(1./2701.), np.log(1./1101.), n)
	z_arr = np.exp(-y0_arr) - 1
	delta_i_LCDM = np.zeros(n)
	delta_i_EdS = np.zeros(n)

	for i in range(len(y0_arr)):
		delta_i_LCDM[i] = findcollshell(y0_arr[i], "LCDM", True)
		delta_i_EdS[i] = findcollshell(y0_arr[i], "EdS", True)

	Z = z_arr[::-1]
	delta_LCDM = delta_i_LCDM[::-1] 
	delta_EdS = delta_i_EdS[::-1]

	mpl.xlim(2750, 1000)

	mpl.plot(Z, delta_LCDM, "r--", linewidth = 0.75, label = r"$\Lambda$CDM")
	mpl.plot(Z, delta_EdS, "b-.", linewidth = 0.75, label = "EdS")
	mpl.ylabel(r"$\delta _i$", rotation = 0, labelpad = 10)
	mpl.xlabel(r"$z_i$", labelpad = 10)
	mpl.legend(loc = 2)
	mpl.title(r"Fitted values of $\delta_i$")

	mpl.savefig("Figures\Fittedforproject.png", bbox_inches = "tight")