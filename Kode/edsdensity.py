# coding=utf-8
import matplotlib.pylab as mpl
import numpy as np

N = 1000
eps = 1e-1

c = 1										#not true
a_0 = 1										#not true
H_0 = 1										#not true
G = 1										#not true

omega_m0 = 1 + c**2/(a_0*H_0)**2
A = omega_m0/(omega_m0 - 1)
B = omega_m0/(2*H_0*(omega_m0 - 1)**(3/2))

theta = np.linspace(0+9*eps, 2*np.pi-eps, N)
t = B*(theta - np.sin(theta))
R_p = A*(1 - np.cos(theta))


#Det er nummerisk feil i delta, men det er også nodt til å vaere teoretisk feil, da verdiene ikke stemmer i det hele tatt. Trenden, at delta forst synker, og saa oker virker fornuftig for meg i det minste!

"""
delta = omega_m0*3/2.*H_0*t/(R_p)**3
count = 0
for i in range(len(delta)):
	count += 1
	if count == 50:
		print i, delta[i]
		count = 0
"""
rho = 1/(6*np.pi*G*t**2)

mpl.plot(t, R_p)
mpl.xlabel("Time [t]")
mpl.ylabel("$R_p$")
mpl.show()

mpl.plot(t, rho)
mpl.xlabel("Time [t]")
mpl.ylabel(r"Density $ \rho $")

mpl.show()