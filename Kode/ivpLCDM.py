import matplotlib.pylab as mpl
import numpy as np
import sys
from scipy.integrate import solve_ivp

def virialcheck(y, mix, Omega_m0, Omega_K, Lambdavar, delta_i, model):
	r = mix

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


		if abs(3/10*(-Omega_m0*(1+delta_i)*(1/(2*r[s])- 1/r[t]) + Lambdavar*r[t]**2)) <= 1e-2:
			controll = False
			p = s			
		s += 1

	if p >= 1:
		rvir = r[p]
	rover = rvir/r[0]
	avir = a[p]

	if rvir:
		odensity = (Omega_m0*(1+delta_i)/rvir**3 + Lambdavar)/(Omega_m0/avir**3 + Lambdavar)


		odensitymax = (Omega_m0*(1+delta_i)/rmax**3 + Lambdavar)/(Omega_m0/amax**3 + Lambdavar)


	if rvir:
		
		k = p
		while k + 1 <= len(r):
			r[k] = rvir
			k += 1

	if ( collapse and rvir ):
		file.write("        {:5.10f} 	  | 		 {:5.10e}    	|     {:5.10e}  	|	{:5.10f}	|	{:5.10e} 	|	{:5.10f}    |     {} \n".format(odensity, odensitymax, rvir/rmax, avir/amax, np.exp(-y[t]) -1, delta_i, model))

	return r


def r(y, x):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (-Omega_m0*(1+delta_i)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda) 	

	return rr

def rmax(y, x):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (-Omega_m0*(1+delta_imax)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda) 	

	return rr

def rmin(y, x):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (-Omega_m0*(1+delta_imin)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda) 	

	return rr

def rmid(y, x):
	r = x[0]
	drdy = x[1]
	a = np.exp(y)
	rr = [[],[]]
	rr[0] = drdy
	rr[1] = (-Omega_m0*(1+delta_imid)/(2.*r**2) + r*Lambda + 3./2.*drdy*Omega_m0/a**3)/(Omega_m0/a**3 + Lambda) 	

	return rr

def kollaps(y, x):
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

kollaps.terminal = True

Omega_m0 = 0.26
Lambda = 0.74

delta_i = float(sys.argv[1])

y0 = np.log(1/(1 + float(sys.argv[2])))
N = 50000
y = np.linspace(y0, 1e-12, N)

r0 = np.exp(y0)
drdx0 = np.exp(y0)

file = open("ivpvalues.txt", "w")

radius = solve_ivp(r, [y0, 1e-12], [r0, drdx0], method = 'RK45', t_eval = y, events = kollaps)
newrad = virialcheck(radius.t, radius.y[1], 0.26, 0.00, 0.74, delta_i, "LCDM")
print radius.t_events[0][0]
x = radius.t_events[0][0]
print x
mpl.plot(radius.t, newrad, "-.", label = "virialcheck")

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
mpl.legend()
mpl.show()
