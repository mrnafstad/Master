#import matplotlib for figures, numpy for array handling, sys for command line arguments and scipy.integrat.solve_ivp for ODE solution
import matplotlib.pylab as mpl
import numpy as np
import sys
from scipy.integrate import solve_ivp

def virialcheck(y, mix, Omega_m0, Lambdavar, delta_i, model):
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
		odensity = (Omega_m0*(1+delta_i)/rvir**3 + Lambdavar)/(Omega_m0/avir**3 + Lambdavar)


		odensitymax = (Omega_m0*(1+delta_i)/rmax**3 + Lambdavar)/(Omega_m0/amax**3 + Lambdavar)
		
		#set all elements in the radial array after virialization to r_vir
		k = p
		while k + 1 <= len(r):
			r[k] = rvir
			k += 1

	if ( collapse and rvir ):
		#write values to file
		file.write("        {:5.10f} 	  | 		 {:5.10e}    	|     {:5.10e}  	|	{:5.10f}	|	{:5.10e} 	|	{:5.10f}    |     {} \n".format(odensity, odensitymax, rvir/rmax, avir/amax, np.exp(-y[t]) -1, delta_i, model))


	return r

def r(y, x):
	#define the rhs of the ode. r and drdy is the radius and rate of change of the perturbation
	#a is the scalefactor
	#rr is  a list storing the return values of drdy and r
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda) 	

	return rr

def kollaps(y, x):
	#This is an event function. Defines r = 0 as an event
	return x[1]

def finddelta(y0):

	N = 50000
	y = np.linspace(y0, 1e-12, N)
	r0 = np.exp(y0)
	drdx0 = np.exp(y0)

	delta_max = 0.1
	delta_min = 0.000001
	xcoll = np.zeros(N)
	ycoll = np.zeros(N)
	zcoll = np.zeros(N)
	koll = [0, 0, 0]
	colltime = 1
	count = 0
	while colltime > 0.001:
		count += 1
		delta_mid = 0.5*(delta_max + delta_min)
		y = np.linspace(y0, 1e-12, N)

		delta_i = delta_max
		radiusmax = solve_ivp(rmax, [y0, 1e-12], [r0, drdx0], method = 'RK45', t_eval = y, events = kollaps)
		print radiusmax
		xcoll[count - 1] = radiusmax.t_events[0][0]
		collmax = np.exp(-xcoll[count - 1]) - 1
		koll[0] = radiusmax.status

		mpl.plot(radiusmax.t, radiusmax.y[1], "-.", linewidth = 0.5)
		y = np.linspace(y0, 1e-12, N)
		delta_i = delta_min
		radiusmin = solve_ivp(rmin, [y0, 1e-12], [r0, drdx0], method = 'RK45', t_eval = y, events = kollaps)
		ycoll[count - 1] = radiusmin.t_events[0][0]
		collmin = np.exp(-ycoll[count - 1]) - 1
		koll[1] = radiusmin.status

		mpl.plot(radiusmin.t, radiusmin.y[1], ".", linewidth = 0.5)
		y = np.linspace(y0, 1e-12, N)
		delta_i = delta_mid
		radiusmid = solve_ivp(rmid, [y0, 1e-12], [r0, drdx0], method = 'RK45', t_eval = y, events = kollaps)
		zcoll[count - 1] = radiusmid.t_events[0][0]
		collmid = np.exp(-zcoll[count - 1]) - 1
		koll[2] = radiusmid.status

		mpl.plot(radiusmid.t, radiusmid.y[1], "--", linewidth = 0.5)



		if ( koll[0] == 1 and koll[2] == 1):
			delta_max = delta_mid
			print "yah", delta_max, collmax

		elif ( koll[0] == 1 and koll[2] != 1):
			delta_min = delta_mid
			"print nah"



		colltime = collmax
		#print len(radiusmax.t_events[0])
		if abs(delta_max - delta_min) < 0.00001:
			print "We may be stuck in a loop at", colltime, count
			print xcoll
			mpl.show()

			break

		#for some reason radiusmax.t_events[0][0] does not update, find out why
		#print "{:.5e} {:.5e} {:.5e}".format(delta_max, delta_min, delta_mid)

	return delta_max

#define that if r = 0 kollaps is a terminal event for solve_ivp, i.e. stopping integration if r = 0
kollaps.terminal = True

#define density parameters for the background
Omega_m0 = 0.26
Lambda = 0.74

#import initial overdensity of the perturbation as a cmd argument
delta_i = float(sys.argv[1])

#define the initial time variable from a cmd argument, given as the redshift z
y0 = np.log(1/(1 + float(sys.argv[2])))
#Define time-array to be returned from solve_ivp
N = 500000
y = np.linspace(y0, 1e-12, N)

#defining initial conditions
r0 = np.exp(y0)
drdx0 = np.exp(y0)

#open a file to write results into
file = open("Numbers\ivpvalues.txt", "w")
file.write("Overdensity at virialization     Overdensity at maximum radius		Radius ratio 		        Time ratio	           Z_coll			Initial Overdensity \n")

#Solve the ODE using solve_ivp
radiuslcdm = solve_ivp(r, [y0, 1e-12], [r0, drdx0], method = 'RK45', t_eval = y, events = kollaps)
print radiuslcdm


mpl.plot(radiuslcdm.t, radiuslcdm.y[1], "r--", linewidth = 0.75, label = "RK45")


radiuslcdm = solve_ivp(r, [y0, 1e-12], [r0, drdx0], method = 'LSODA', t_eval = y, events = kollaps)
print radiuslcdm


mpl.plot(radiuslcdm.t, radiuslcdm.y[1], "b-", linewidth = 0.75, label = "LSODA")

radiuslcdm = solve_ivp(r, [y0, 1e-12], [r0, drdx0], method = 'RK23', t_eval = y, events = kollaps)

print radiuslcdm


mpl.plot(radiuslcdm.t, radiuslcdm.y[1], "c:", linewidth = 0.75, label = "RK23")

radiuslcdm = solve_ivp(r, [y0, 1e-12], [r0, drdx0], method = 'BDF', t_eval = y, events = kollaps)
print radiuslcdm


mpl.plot(radiuslcdm.t, radiuslcdm.y[1], "g--", linewidth = 0.75, label = "BDF")

radiuslcdm = solve_ivp(r, [y0, 1e-12], [r0, drdx0], method = 'Radau', t_eval = y, events = kollaps)
print radiuslcdm


mpl.plot(radiuslcdm.t, radiuslcdm.y[1], "m-.", linewidth = 0.75, label = "Radau")


"""

delta_i = 0.002
radius = solve_ivp(r, [y0, 1e-12], [r0, drdx0], method = 'RK45', t_eval = y, events = kollaps)
print radius.t_events[0][0]
#print "Collapse at z = {:2.5f}".format(np.exp(-radius.t_events[0][0]) - 1)
x = radius.t_events[0][0]
print x
mpl.plot(radius.t, radius.y[1])


delta_i = 0.0009
radius = solve_ivp(r, [y0, 1e-12], [r0, drdx0], method = 'RK45', t_eval = y, events = kollaps)

file.close()
print radius.t_events[0][0]
#print "Collapse at z = {:2.5f}".format(np.exp(-radius.t_events[0][0]) - 1)
x = radius.t_events[0][0]
print x
mpl.plot(radius.t, radius.y[1])

#delta_i = finddelta(y0)

#radius = solve_ivp(r, [y0, 1e-12], [r0, drdx0], method = 'RK45', t_eval = y, events = kollaps)

#mpl.plot(radius.t, radius.y[1], "r-.", label = "fitted")
"""
file.close()
mpl.legend()
mpl.xlabel("x(a)", labelpad = 10)
mpl.ylabel(r"$\~R$  ", rotation = 0, labelpad = 10)
mpl.show()


