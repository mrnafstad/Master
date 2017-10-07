import matplotlib.pylab as mpl
import numpy as np

N = 1000
eps = 1e-5

c = 1										#not true
a_0 = 1										#not true
H_0 = 1										#not true

omega_m0 = 1+c**2/(a_0*H_0)**2
A = omega_m0/(omega_m0 - 1)
B = omega_m0/(2*H_0*(omega_m0 - 1)**(3/2))

theta = np.linspace(0+eps, 2*np.pi-eps, N)
t = B*(theta - np.sin(theta))
R_p = A*(1 - np.cos(theta))


delta = omega_m0*(3/2*H_0*t/(R_p))**3

mpl.plot(t, R_p)
mpl.show()

mpl.plot(t, delta)
mpl.show()